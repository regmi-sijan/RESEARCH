#include <truthparticleana/TruthParticleAna.h>
//#include "GetTChainMacro.C"
#include <GetTChainMacro.C>

// check the input file, output file and their path before running this macro
void run_loop_truth(int n_file = 101)
{
	int nevents = -1;

	const std::string i_file_name = "run10_truth_multi-par_embedding_";
	//const std::string i_file_name = "run10_no_sec_eta_";
	//const std::string i_file_name = "run10_truth_par_pri_";
	//const std::string i_file_name = "final_data_list_";
	
	const std::string macro_loc = "/sphenix/user/sregmi/WORKING_AREA/truth_par_ana/study_truth_particle/macro";
	
	// check this and down than this to confirm input filelist
	const std::string ifile= macro_loc + "/file_list/" + i_file_name + std::to_string(n_file); // input file-list
	const std::string ofile= macro_loc + "/RESULTS/JulyAugust/delR_pT_loop_" + i_file_name + std::to_string(n_file) + ".root";  // output file
	
	const std::string ref_file = macro_loc + "/RESULTS/JulyAugust/ref_file_delR_weights.root"; // ref file taken for weight
	
	R__LOAD_LIBRARY(libTruthParticleAna.so)	

	TruthParticleAna truth_obj("TruthParticleAna", ofile.c_str());
	//TruthParticleAna truth_obj("TruthParticleAna", ofile.c_str());
  truth_obj.InitRun(0);
  TTree * intree1 =  GetTChainMacro(ifile.c_str());

  //truth_obj.Loop(nevents, ifile.c_str(), intree1);
  //truth_obj.Loop_ETAtoPIONS(nevents, ifile.c_str(), intree1);
	//truth_obj.check_primary_secondary(nevents, ifile.c_str(), intree1);
 
	//truth_obj.Loop_background_position_swappping(nevents, ifile.c_str(), intree1); 
	
	truth_obj.Loop_background_event_mixing(nevents, ifile.c_str(), ref_file.c_str(), intree1);
	
  truth_obj.End(0);
	std::cout << "All Done" << std::endl;
	gSystem->Exit(0);

}
