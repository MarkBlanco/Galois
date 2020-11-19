/*
 * This file belongs to the Galois project, a C++ library for exploiting parallelism.
 * The code is being released under the terms of the 3-Clause BSD License (a
 * copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2018, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 */

#include "Lonestar/BoilerPlate.h"
#include "pi_constants.h"
#include "galois/Galois.h"
#include "galois/LargeArray.h"
#include "galois/Timer.h"
#include "galois/graphs/LCGraph.h"
#include "galois/graphs/TypeTraits.h"
#include "galois/gstl.h"

unsigned long long rdtsc()
{
 unsigned a, d;

 __asm__ volatile("rdtsc" : "=a" (a), "=d" (d));

 return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}


// #define DEBUG 1
const char* desc =
    "Computes protein interaction ranks a la Page and Brin. This is a pull-style algorithm.";

enum Algo { Topo = 0, Residual };
const char* const ALGO_NAMES[] = {"Topological", "Residual"};

static cll::opt<Algo> algo("algo", cll::desc("Choose an algorithm:"),
                           cll::values(clEnumVal(Topo, "Topological"),
                                       clEnumVal(Residual, "Residual"),
                                       clEnumValEnd),
                           cll::init(Residual));

constexpr static const unsigned CHUNK_SIZE = 32;

struct LNode {
  PRTy value;
  PRTy restart_init;
  uint32_t nout;
  PRTy normalized_out_weight;
};

typedef galois::graphs::LC_CSR_Graph<LNode, int>::with_no_lockable<
    true>::type ::with_numa_alloc<true>::type Graph;
typedef typename Graph::GraphNode GNode;

using DeltaArray    = galois::LargeArray<PRTy>;
using ResidualArray = galois::LargeArray<PRTy>;

//! [example of no_stats]
// void initNodeDataTopological(Graph& g) {
//   galois::do_all(galois::iterate(g),
//                  [&](const GNode& n) {
//                    auto& sdata = g.getData(n, galois::MethodFlag::UNPROTECTED);
//                    sdata.value = INIT_RESIDUAL;
//                    sdata.nout  = 0;
//                  },
//                  galois::no_stats(), galois::loopname("initNodeData"));
// }
//! [example of no_stats]

void initNodeDataResidual(Graph& g) {
  // std::cout<<"Init\n";
  galois::do_all(galois::iterate(g),
                 [&](const GNode& n) {
                   auto& sdata = g.getData(n, galois::MethodFlag::UNPROTECTED);
                   sdata.value = (1.0/g.size());
                   sdata.restart_init = INIT_RESIDUAL*(1.0/g.size());
                   sdata.nout  = 0;
                   sdata.normalized_out_weight = 1.0;
                   // std::cout<<sdata.value<<std::endl;
                 },
                 galois::no_stats(), galois::loopname("initNodeData-HI"));
}

void computeoutdegsums(Graph& graph, galois::LargeArray<std::atomic<size_t>>& vec) {
  //initialize all out degrees to 0
  galois::do_all(galois::iterate(graph),
                 [&](const GNode& src) { vec.constructAt(src, 0ul); },
                 galois::no_stats(), galois::loopname("InitDegVec-IGNORE"));

  //compute the sum of outgoing weights of each node
  galois::do_all(galois::iterate(graph),
                 [&](const GNode& src) {
                   // std::cout<<src<< std::endl;
                   for (auto nbr : graph.edges(src)) {
                     GNode dst = graph.getEdgeDst(nbr);
                     // auto weight = graph.getEdgeData(nbr);
                     // std::cout<<"Edge Weight:"<<weight<<" "<<graph.getEdgeData(nbr) <<"\t";
                     //TODO:use real edgeweight
                     //std::cout<<"ew"<<weight<<std::endl;
                     //assume ew=1
                     vec[dst].fetch_add(graph.getEdgeData(nbr));
                   };
                   //std::cout<<std::endl;
                 },
                 galois::steal(), galois::chunk_size<CHUNK_SIZE>(),
                 galois::no_stats(), galois::loopname("computeOutWeight-IGNORE"));
/*
	for (uint32_t i = 0; i < 10; i++){
		std::cout << vec[i] << std::endl;
	}*/
}

// Computing outdegrees in the tranpose graph is equivalent to computing the
// indegrees in the original graph
void computeoutdeg(Graph& graph, galois::LargeArray<std::atomic<size_t>>& vec) {
  galois::StatTimer outDegreeTimer("computeOutDegFunc-HI");
  outDegreeTimer.start();

  //store the reciprocal of weights
  galois::do_all(galois::iterate(graph),
                 [&](const GNode& src) {
                   auto& srcData =
                       graph.getData(src, galois::MethodFlag::UNPROTECTED);
                  if(vec[src]>0){
                    // std::cout<<src<<" : "<<vec[src]<< " "<< 1.0/vec[src]<<std::endl;
                     srcData.normalized_out_weight = 1.0/vec[src];
                   }
                 },
                 galois::no_stats(), galois::loopname("CopyNorm-HI"));

  outDegreeTimer.stop();
}

//! [scalarreduction]
void computePRResidual(Graph& graph) {
  unsigned int iterations = 0;
  galois::GAccumulator<PRTy> accum;


  while (true) {
    galois::do_all(galois::iterate(graph),
                   [&](const GNode& src) {
                     auto& sdata = graph.getData(src);
                     PRTy update = 0.0;
                     //pull updates from source

                     // std::cout<<src<<" = " << ALPHA <<"( ";
                     for (auto nbr : graph.edges(src)) {
                       GNode dst = graph.getEdgeDst(nbr);
                       auto ew = graph.getEdgeData(nbr);
                       // std::cout << ew <<"\t";
                       auto& ddata = graph.getData(dst);
                       // std::cout<< ew <<"("<< dst <<")" <<"*" <<ddata.value << "*" <<ddata.normalized_out_weight<< " + ";
                       update += ew*ddata.value*ddata.normalized_out_weight;
                     }

                     update *= ALPHA;
                     update += sdata.restart_init;
                     accum += fabs(update-sdata.value);
                     // std::cout<< "\b\b) + "<<sdata.restart_init<<"\n = "<<update<<"" <<" change: " <<update-sdata.value << "\n\n\n";
                     sdata.value = update;

                   },
                   galois::steal(), galois::chunk_size<CHUNK_SIZE>(),
                   galois::no_stats(), galois::loopname("PageRank-HI"));

#if DEBUG
    std::cout << "iteration: " << iterations << "\n";
    std::cout<<"cummulative error"<<accum.reduce()<<std::endl;
#endif



    iterations++;

    if (iterations >= maxIterations || accum.reduce()<tolerance) {
      break;
    }
    accum.reset();
  } // end while(true)

  if (iterations >= maxIterations) {
    std::cerr << "ERROR: failed to converge in " << iterations << " iterations"
              << std::endl;
  }
}
//! [scalarreduction]

void prResidual(Graph& graph, galois::LargeArray<std::atomic<size_t>>& in_weights) {
  //DeltaArray delta;
  //delta.allocateInterleaved(graph.size());
  // ResidualArray residual;
  // residual.allocateInterleaved(graph.size());

  //galois::StatTimer prTimer("INSIDE PR");
  //prTimer.start();
  initNodeDataResidual(graph);
  computeoutdeg(graph, in_weights);
  computePRResidual(graph);
  //prTimer.stop();
  //std::cout<<"Time: "<<prTimer.get()<<std::endl;
}

int main(int argc, char** argv) {
  std::cout<<"This is pull\n";

  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url);

  galois::StatTimer overheadTime("OverheadTime");
  overheadTime.start();

  Graph transposeGraph;
  std::cout << "WARNING: pull style algorithms work on the transpose of the "
               "actual graph\n"
            << "WARNING: this program assumes that " << filename
            << " contains transposed representation\n\n"
            << "Reading graph: " << filename << std::endl;

  galois::graphs::readGraph(transposeGraph, filename);
  std::cout << "Read " << transposeGraph.size() << " nodes, "
            << transposeGraph.sizeEdges() << " edges\n";

  galois::preAlloc(2 * numThreads + (3 * transposeGraph.size() *
                                     sizeof(typename Graph::node_data_type)) /
                                        galois::runtime::pagePoolSize());
  galois::reportPageAlloc("MeminfoPre");

  // switch (algo) {
  // case Topo: {
  //   std::cout << "Running Pull Topological version, tolerance:" << tolerance
  //             << ", maxIterations:" << maxIterations << "\n";
  //   prTopological(transposeGraph);
  //   break;
  // }
  // case Residual: {
#define ITERS 16
	galois::StatTimer t_main("LOOK AT ME");
	uint64_t acc = 0;
	galois::LargeArray<std::atomic<size_t>>  in_weights;
	in_weights.allocateInterleaved(transposeGraph.size());
	computeoutdegsums(transposeGraph, in_weights);
	for (uint32_t i = 0; i < ITERS; i++){
		std::cout << "Running Pull Residual version, tolerance:" << tolerance
			<< ", maxIterations:" << maxIterations << "\n";
		uint64_t st = rdtsc();//omp_get_wtime();
		//t_main.start();	
		prResidual(transposeGraph, in_weights);
		uint64_t nd = rdtsc();//omp_get_wtime();
		//t_main.stop();
		acc += nd-st; //t_main.get();
	}
	double time = acc / ((double)ITERS );
	printf("Average time for %d trials: %f CYCLES CONVERT THIS\n", ITERS, time);
	
  //   break;
  // }
  // default: { std::abort(); }
  // }

  //galois::reportPageAlloc("MeminfoPost");

  /*
	 * if (!skipVerify) {
    printTop(transposeGraph);
  }*/

// #if DEBUG
//   printPageRank(transposeGraph);
// #endif

  overheadTime.stop();
  return 0;
}
