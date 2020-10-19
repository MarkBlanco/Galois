#! /bin/bash

OUTPUT=run_log.dat
echo $DATE >> ${OUTPUT}
hostname   >> ${OUTPUT}

./pi_run.sh 	| tee -a ${OUTPUT}
./pr_run.sh 	| tee -a  ${OUTPUT}
#./sssp_run.sh | tee -a ${OUTPUT}
#./sir_run.sh 	| tee -a  ${OUTPUT}

echo "All done."
