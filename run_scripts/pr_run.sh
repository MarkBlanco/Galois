#! /bin/bash

# Get DATADIR, SOURCESDIR, and GKCBASE, NUM_THREADS
source ./setup.sh

EXECS="pagerank-pull" 
ALGO="PageRank"
OUTDIR="${GKCBASE}/${ALGO}_outputs/"
mkdir ${OUTDIR}

for EXEC in ${EXECS} 
do
	echo ${EXEC}
	# for symmetric graphs
	for GRAPH in road kron urand
	do
		name=${GRAPH}
		echo $name
		OUTPUT="${OUTDIR}${GRAPH}_${EXEC}.dat"
		echo $DATE >> ${OUTPUT}
		hostname   >> ${OUTPUT}
		echo ${EXEC} >> ${OUTPUT}
		CMD="./${EXEC} \
			-t=${NUM_THREADS} \
			${DATADIR}${name}.sgr" 
		echo $CMD
		./${EXEC} \
			"-t=${NUM_THREADS}" \
			"${DATADIR}${name}.sgr" >> ${OUTPUT} 2>&1
	done

	# For non-symmetric graphs
	for GRAPH in web twitter
	do
		name=${GRAPH}
		echo $name
		OUTPUT="${OUTDIR}${GRAPH}_${EXEC}.dat"
		echo $DATE >> ${OUTPUT}
		hostname   >> ${OUTPUT}
		echo ${EXEC} >> ${OUTPUT}
		CMD="./${EXEC} \
			-t=${NUM_THREADS} \
			${DATADIR}${name}.tgr" 
		echo $CMD
		./${EXEC} \
			"-t=${NUM_THREADS}" \
			"${DATADIR}${name}.tgr" >> ${OUTPUT} 2>&1
	done
done
