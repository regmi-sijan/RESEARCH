// these include guards are not really needed, but if we ever include this
// file somewhere they would be missed and we will have to refurbish all macros
#ifndef MACRO_FUN4ALLTRUTHPARANA_C
#define MACRO_FUN4ALLTRUTHPARANA_C

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <../src/TruthParticleAna.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libTruthParticleAna.so)

// this macro can be used to read DST file and save information from there
void Fun4All_TruthParAna(
			const int file_num = 0
			)

{	
	const int nEvents = -1;

	// input and output file	
	//string in_truthFile = "DST_TRUTH_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000010-";

	string in_truthFile = "DST_TRUTH_pythia8_pp_mb_3MHz-0000000011-";
  string outputFile = "/sphenix/user/sregmi/WORKING_AREA/saved_ntuples/test_run11_primary_par_ntuple_";
	
	//"/sphenix/user/sregmi/WORKING_AREA/saved_ntuples/run10_truth_par_pri/run10_primary_par_ntuple_";
	
	int ynum_int = 100000 + file_num;
  TString yn_tstr = "";
  yn_tstr += ynum_int;
  yn_tstr.Remove(0,1);

  in_truthFile += yn_tstr.Data();
  outputFile += yn_tstr.Data();

  in_truthFile += ".root";
  outputFile += ".root";
	
	// this convenience library knows all our i/o objects so you don't
  // have to figure out what is in each dst type
  gSystem->Load("libg4dst.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);  // set it to 1 if you want event printouts

	/*
  Fun4AllInputManager *inCluster = new Fun4AllDstInputManager("DSTClusters");
  std::cout << "Adding file list " << clusterFile << std::endl;
  inCluster -> AddListFile(clusterFile,1);
  se->registerInputManager(inCluster);
	*/

  Fun4AllInputManager *in_file = new Fun4AllDstInputManager("DSTCaloTruth");
  std::cout << "Adding truth file  " << in_truthFile  << std::endl;
	in_file->AddFile(in_truthFile);
  se -> registerInputManager(in_file);

  TruthParticleAna *teval = new TruthParticleAna("Truth_Ana", outputFile);
  se->registerSubsystem(teval);
  
  se->run(nEvents);
  se->End();
  
  delete se;
  cout << "Analysis Completed" << endl;
  
  gSystem->Exit(0);
}

#endif  //MACRO_FUN4ALLG4SLOPECAL_C
