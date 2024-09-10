#include <pi0clusterana/pi0ClusterAna.h>
//#include "GetTChainMacro.C"
#include <GetTChainMacro.C>

// check the input file, output file and their path before running this macro
void run_loop_truth( int n_file = 11 )
{
	int nevents = -1;
	 
	const std::string i_file_name = "truth_run10_ntuple_";
	const std::string macro_loc = "/sphenix/user/sregmi/WORKING_AREA/truth_par_ana/pi0ClusterAna/macros";

	// check this and down than this to confirm input filelist
	const std::string ifile= macro_loc + "/file_list/" + i_file_name + std::to_string(n_file); // input file-list
	const std::string ofile= macro_loc + "/RESULTS/truth_loop/results_" + i_file_name + std::to_string(n_file) + ".root";  // output file
	
	R__LOAD_LIBRARY(libpi0ClusterAna.so)	

	pi0ClusterAna truth_obj("pi0ClusterAna", ofile.c_str());
  truth_obj.InitRun(0);
  TTree * intree1 =  GetTChainMacro(ifile.c_str());
  truth_obj.Loop(nevents, ifile.c_str(), intree1);
  truth_obj.End(0);
	gSystem->Exit(0);

}
