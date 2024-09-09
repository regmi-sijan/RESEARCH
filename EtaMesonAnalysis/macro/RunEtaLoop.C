#include <calib_emc_eta/CaloCalibEmc_eta.h>
#include "GetTChainMacro.C"
//#include <GetTChainMacro.C>

void RunEtaLoop(const int n_file = 101)
{
  int nevents = -1;
	
  //const std::string in_fname = "run24_pp_ntuple_";
  const std::string in_fname = "run10_hijing_ntuple_";

  const std::string macro_loc = "/sphenix/user/sregmi/WORKING_AREA/ETA_MESON/macro";
	
  // check this and down than this to confirm input filelist
  const std::string ifile= macro_loc + "/file_list/" + in_fname + std::to_string(n_file); // input file-list
	
  const std::string ofile= macro_loc + "/RESULTS/hijing/cent_run10_loop_" + std::to_string(n_file) + ".root";  
  //const std::string ofile= macro_loc + "/RESULTS/RUN_24/new_run24_pp_loop_" + std::to_string(n_file) + ".root";
	
  std::cout << " input file list " << ifile << std::endl;
  std::cout << " output file list " << ofile << std::endl;
	
  R__LOAD_LIBRARY(libcalibCaloEmc_eta.so)	

  CaloCalibEmc_eta calo_obj("CaloCalibEmc_eta", ofile.c_str());
  calo_obj.InitRun(0);
  TTree * intree1 =  GetTChainMacro(ifile.c_str());
  
  //calo_obj.Loop(nevents, ifile.c_str(), intree1); // creates foreground
  calo_obj.Loop_background_event_mixing(nevents, ifile.c_str(), intree1); // for creating background (using event mixing technique)
  //calo_obj.Loop_background_position_swapping(nevents, ifile.c_str(), intree1); // for creating background (using event mixing technique)
  
  calo_obj.End(0);
  std::cout << "All Done" << std::endl;
  gSystem->Exit(0);
}

