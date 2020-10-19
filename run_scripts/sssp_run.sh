#! /bin/bash

# Get DATADIR, SOURCESDIR, and GKCBASE, NUM_THREADS
source ./setup.sh

EXECS="sssp-pull" 
ALGO="SSSP"
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
			${DATADIR}${name}.sgr < ${SOURCESDIR}${GRAPH}.sources" 
		echo $CMD
		./${EXEC} \
			"-t=${NUM_THREADS}" \
			"${DATADIR}${name}.sgr" < "${SOURCESDIR}${GRAPH}.sources" >> ${OUTPUT} 2>&1
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
			${DATADIR}${name}.gr < ${SOURCESDIR}${GRAPH}.sources" 
		echo $CMD
		./${EXEC} \
			"-t=${NUM_THREADS}" \
			"${DATADIR}${name}.gr" < "${SOURCESDIR}${GRAPH}.sources" >> ${OUTPUT} 2>&1
	done
done
