#!/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${sregmi}
export HOME=/sphenix/user/${sregmi}

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
export MYINSTALL=/sphenix/user/sregmi/install/
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

printenv 

echo ' '
echo 'START: '`date`
echo ' '


cd /sphenix/user/sregmi/WORKING_AREA/truth_par_ana/pi0ClusterAna/macros/

z=$(($1+11))

echo starting job for loop

 root.exe -b <<EOF
	 .L run_loop_truth.C 
	  run_loop_truth($z);
 .q

EOF

echo ' '
echo ' '
echo 'END: '`date`
echo ' '
