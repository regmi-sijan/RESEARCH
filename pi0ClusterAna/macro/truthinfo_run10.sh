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

z=$(($1+40000))

echo starting processing to make ntuple from truth particle

 root.exe -b <<EOF
		.L Fun4All_RunPi0ClusterAna.C
		Fun4All_RunPi0ClusterAna($z)
 .q

EOF

echo ' '
echo ' '
echo 'END: '`date`
echo ' '
