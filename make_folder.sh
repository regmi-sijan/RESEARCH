#!/bin/bash

# Define the path to the filelist
file_path="/sphenix/user/sregmi/WORKING_AREA/PI_0/wvfm_pi0_tbt/macros/file_list/filelist"

# Read and process all lines
while IFS= read -r line; do

    # Split the line into two variables based on space
    n_run="$(echo "$line" | awk '{print $1}')"
    seg="$(echo "$line" | awk '{print $2}')"

		mkdir "run_$n_run"
		hadd run_"$n_run"_CEMC-result_combined.root run_"$n_run"_*.root
		mv run_"$n_run"_CEMC-result* run_"$n_run" 
    
done < "$file_path"

