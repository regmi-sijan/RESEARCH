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
 

echo to make a plot of start time of waveform


cd /sphenix/user/sregmi/WORKING_AREA/PI_0/wvfm_pi0_tbt/macros/

z=$(($1+1)) 

# Define the path to the filelist
file_path="/sphenix/user/sregmi/WORKING_AREA/PI_0/wvfm_pi0_tbt/macros/file_list/filelist_$z"

# Read and process all lines
while IFS= read -r line; do

    # Split the line into two variables based on space
    n_run="$(echo "$line" | awk '{print $1}')"
    seg="$(echo "$line" | awk '{print $2}')"

		for((i = 0; i < seg; i++))
		do
			echo "$n_run, $i"
			
			root.exe -b << EOF
				.L Fun4All_pi0_wvfm_tbt.C
				Fun4All_pi0_wvfm_tbt($n_run, $i);
			.q

EOF

		done
    
done < "$file_path"

