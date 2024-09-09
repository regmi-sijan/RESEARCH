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


echo job to create background and foreground

cd /sphenix/user/sregmi/WORKING_AREA/ETA_MESON/macro/

 
z=$(($1+101))

echo starting job for loop


 root.exe -b <<EOF
	 .L RunEtaLoop.C 
	  RunEtaLoop($z);
 .q

EOF

echo ' '
echo ' '
echo 'END: '`date`
echo ' '
