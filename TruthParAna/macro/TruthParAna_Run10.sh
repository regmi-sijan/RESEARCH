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

cd /sphenix/user/sregmi/WORKING_AREA/truth_par_ana/study_truth_particle/macro/

z=$(($1+101))

echo starting processing to make ntuple from truth particle

 root.exe -b <<EOF
	 .L run_loop_truth.C
	 run_loop_truth($z)
 .q

EOF

echo ' '
echo ' '
echo 'END: '`date`
echo ' '
