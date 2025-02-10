#!/bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/user/${LOGNAME} 

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
export MYINSTALL=/sphenix/user/sregmi/install/
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

 
this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`

echo rsyncing from $this_dir
echo running: $this_script $*

# input file list to pass into the macro
inputfile=$1

if [[ -z "$inputfile" ]]; then
   echo "Error: No input file provided."
   echo "Usage: $0 <inputfile>"
   exit 1
fi

#source /opt/sphenix/core/bin/sphenix_setup.sh -n
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]
then
  cd $_CONDOR_SCRATCH_DIR
  rsync -av $this_dir/*.C .
  rsync -av $this_dir/*.sh .
  rsync -av $this_dir/new_dst_truth.list .
  rsync -av $this_dir/*.job .
  getinputfiles.pl $inputfile
else
  echo condor scratch NOT set
  exit -1
fi

#cd /sphenix/user/sregmi/WORKING_AREA/truth_par_ana/EtaDecayAna/macros/
 
# this is how you run your Fun4All_G4_sPHENIX.C macro which was rsynced from your initial dir: 
root.exe -b <<EOF
	.L Fun4All_Eta_SIMPLE_EMBED.C
	Fun4All_Eta_SIMPLE_EMBED("$inputfile")
.q

EOF 

echo ' '
echo ' '
echo 'END: '`date`
echo ' '
