#include "TruthParticleAna.h"

//Fun4all stuff
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <ffaobjects/EventHeader.h>

//ROOT stuff
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <vector>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>

//for emc clusters
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>

//for vetex information
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

//tracking
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
//#include <trackbase_historic/SvtxVertex.h>
//#include <trackbase_historic/SvtxVertexMap.h>


//truth information
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/IteratorRange.h>
#include <HepMC/SimpleVector.h> 
#include <HepMC/GenParticle.h>
#pragma GCC diagnostic pop
#include <CLHEP/Geometry/Point3D.h>


#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>

#include <algorithm>
#include <vector>
#include <TLorentzVector.h>

#include "TRandom3.h"
#include <random>

#include <THn.h>
#include <TFile.h>
#include <TRandom.h>
#include <iostream>
#include <THnSparse.h>

#include <map>
#include <cmath>
#include <iomanip>  // for std::setprecision
#include <sstream>  // for std::ostringstream

//____________________________________________________________________________..
TruthParticleAna::TruthParticleAna(const std::string &name, const std::string &outName = "TruthParticleAnaOut"):
SubsysReco(name)
, clusters_Towers(nullptr)
  , truth_photon(nullptr)
  , truth_pi0(nullptr)
  , alphaCut(-1.0)
  , n_event(0)
  , f_temp(0)
, Outfile(outName)

{
  std::cout << "TruthParticleAna::TruthParticleAna(const std::string &name) Calling constructor" << std::endl;
}

//____________________________________________________________________________..
TruthParticleAna::~TruthParticleAna()
{
  std::cout << "TruthParticleAna::~TruthParticleAna() Calling destructor" << std::endl;
}

//____________________________________________________________________________..
int TruthParticleAna::Init(PHCompositeNode *topNode)
{
  n_event = 0; // initialization of event number

  std::cout << "TruthParticleAna::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TruthParticleAna::InitRun(PHCompositeNode *topNode)
{
  std::cout << "TruthParticleAna::InitRun(PHCompositeNode *topNode) Initializing" << std::endl;
	
  out = new TFile(Outfile.c_str(),"RECREATE");


  pairInvMassPtdelR = new TH3F("pairInvMassPtdelR", "Pair_Mass_Histo_PtdelR", 300, 0.0, 3.0, 200, 0.0, 20.0, 150, 0.0, 1.5);
  pairInvMassPtdelRBkgd = new TH3F("pairInvMassPtdelRBkgd", "Pair_Mass_Histo_PtdelR_Bkgd", 300, 0.0, 3.0, 200, 0.0, 20.0, 150, 0.0, 1.5);

  pairInvMassTotal = new TH1F("pairInvMassTotal", "invariant mass histogram", 500, 0.0, 4.0);
  pairInvMassTotalBkgd = new TH1F("pairInvMassTotalBkgd", "Pair_Mass_Histo_Bkgd", 500, 0.0, 4.0);
	
  pairInvMassPtEta = new TH3F("pairInvMassPtEta", "invariant mass histogram", 500, 0.0, 4.0, 150, 0.0, 15.0, 250, -1.5, 1.5);
  pairInvMassPtEtaBkgd = new TH3F("pairInvMassPtEtaBkgd", "Pair Mass Histo Bkgd", 500, 0.0, 4.0, 150, 0.0, 15.0, 250, -1.5, 1.5);
	
  pairpTDelPhiDelEta = new TH3F("pairpTDelPhiDelEta", "", 150, 0.0, 15.0, 130, -6.5, 6.5, 250, -4.0, 4.0);
  pairpTDelPhiDelEtaBkgd = new TH3F("pairpTDelPhiDelEtaBkgd", "", 150, 0.0, 15.0, 130, -6.5, 6.5, 250, -4.0, 4.0);
	
  delphi_f = new TH1F("del_phi_f", "del_phi_hist", 130, -6.5, 6.5);  
  delphi_b = new TH1F("del_phi_b", "del_phi_hist", 130, -6.5, 6.5);
	  
  deleta_f = new TH1F("del_eta_f", "del_eta_hist", 250, -4.0, 4.0);
  deleta_b = new TH1F("del_eta_b", "del_eta_hist", 250, -4.0, 4.0);

  Pair_PT_f = new TH1F("Pair_PT_f", "pair pT Dist", 150, 0.0, 15.0);
  Pair_PT_b = new TH1F("Pair_PT_b", "pair pT Dist", 150, 0.0, 15.0);
	
  delR_f = new TH1F("delR_f", "delR_hist_foreground", 200, 0.0, 2.0);
  delR_b = new TH1F("delR_b", "delR_hist_background", 200, 0.0, 2.0);
	
  DelPhi_DelEta_f = new TH2F("DelPhiDelEta_f", "DelPhi_DelEta", 130, -6.5, 6.5, 250, -4.0, 4.0);
  DelPhi_DelEta_b = new TH2F("DelPhiDelEta_b", "DealPhi_DelEta", 130, -6.5, 6.5, 250, -4.0, 4.0);	

  DelPhi_PairpT_f = new TH2F("DelPhi_PairpT_f", "DelPhi_pairpT", 130, -6.5, 6.5, 150, 0.0, 15.0);
  DelPhi_PairpT_b = new TH2F("DelPhi_PairpT_b", "DelPhi_pairpT", 130, -6.5, 6.5, 150, 0.0, 15.0);
	
  DelEta_PairpT_f = new TH2F("DelEta_PairpT_f", "DelEta_pairpT", 250, -4.0, 4.0, 150, 0.0, 15.0);
  DelEta_PairpT_b = new TH2F("DelEta_PairpT_b", "DelEta_pairpT", 250, -4.0, 4.0, 150, 0.0, 15.0);
	
  Pri1_Pri2_pairpT = new TH3F("Pri1_Pri2_pairpT", "primary particles", 1000, -3400, 3400, 1000, -3400, 3400, 140, 0.0, 14.0);

  DelR_pairpT_f = new TH2F("DelR_pairpT_f", "DelR_pairpT", 200, 0.0, 2.0, 150, 0.0, 15.0);
  DelR_pairpT_b = new TH2F("DelR_pairpT_b", "DelR_pairpT", 200, 0.0, 2.0, 150, 0.0, 15.0);


	
  if (topNode != 0)
    {
      // TTree declare
      _eventTree = new TTree("_eventTree", "An event level tree info for truth particle");
      // TTree branches
      _eventTree->Branch("_eventNumber", &_eventNumber, "_eventNumber/I");
      _eventTree->Branch("_nFourVector", &_nFourVector, "_nFourVector/I");
      _eventTree->Branch("_fv_energy", _fv_energy, "_fv_energy[_nFourVector]/F");
      _eventTree->Branch("_fv_px", _fv_px, "_fv_px[_nFourVector]/F");
      _eventTree->Branch("_fv_py", _fv_py, "_fv_py[_nFourVector]/F");
      _eventTree->Branch("_fv_pz", _fv_pz, "_fv_pz[_nFourVector]/F");
      _eventTree->Branch("_fv_Eta", _fv_Eta, "_fv_Eta[_nFourVector]/F");
      _eventTree->Branch("_primary_id", _primary_id, "_primary_id[_nFourVector]/F"); // this is same for both primary and secondary particles
      _eventTree->Branch("_parent_id", _parent_id, "_parent_id[_nFourVector]/F"); // this information is useful (used) for secondary particles only
      _eventTree->Branch("_pid_primary", _pid_primary, "_pid_primary[_nFourVector]/F"); // this is pid of primary particles
      _eventTree->Branch("_pid_secondary", _pid_secondary, "_pid_secondary[_nFourVector]/F"); // this is pid for secondar particles
      _eventTree->Branch("_embedding", _embedding, "_embedding[_nFourVector]/F"); // indicatior (flag) to differentiate embedded particles with others

    }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TruthParticleAna::process_event(PHCompositeNode *topNode)
{
  _eventNumber = n_event; // updating event number to the Ttree

  if (n_event%10 == 0) { std::cout << " processing event number " << n_event << std::endl;}

  //truth particle information
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if(!truthinfo)
    {
      std::cout << PHWHERE << "pi0ClusterAna::process_event Could not find node G4TruthInfo"  << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
	
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

  // vectors to store the information for ntuple
  std::vector<float> _PX;
  std::vector<float> _PY;
  std::vector<float> _PZ;
  std::vector<float> _E;
  std::vector<float> _Eta;
  std::vector<int> _PID_PRI;
  std::vector<int> _PID_SEC;
  std::vector<int> _primary_id;

  // preparing all vectors 
  _PX.clear();
  _PY.clear();
  _PZ.clear();
  _E.clear();
  _Eta.clear();
  _PID_PRI.clear();
  _PID_SEC.clear();
  _primary_id.clear();

  /*
  // identify photon and electron (looking for secondary particle)
  PHG4TruthInfoContainer::Range truthRange_sec = truthinfo->GetSecondaryParticleRange();
  PHG4TruthInfoContainer::ConstIterator truthIter_sec;

  PHG4Particle *truth_sec;

  for(truthIter_sec = truthRange_sec.first; truthIter_sec != truthRange_sec.second; truthIter_sec++)
  {
  truth_sec = truthIter_sec->second;

  // do not accept if pseudorapidity falls outside the range of detector (sPHENX) acceptance
  if (fabs(getEta(truth_sec)) > 1.1)  continue;
						
  if ((truth_sec->get_parent_id()) ==0) 
  {
  std::cout << "Primary Particle Spotted" << std::endl;
  continue;
  }
					 
  int pid_sec = truth_sec->get_pid();
  //if ((truth_sec->get_parent_id()) == 221) {std::cout << "Eta meson spotted" << std::endl;}

  if ((pid_sec != 11) && (pid_sec != -11) && (pid_sec != 22)) continue; // ignore if we do not have electron or photon cluster as secondary particle 

  //if ((truth_sec->get_e()) < 0.5) continue;

  _primary_id.push_back(truth_sec->get_primary_id());
  _PID_SEC.push_back(truth_sec->get_pid()); // for secondary particle
  _Eta.push_back(getEta(truth_sec)); // updating pseudorapidity of secondary particle
  _PX.push_back(truth_sec->get_px());
  _PY.push_back(truth_sec->get_py());
  _PZ.push_back(truth_sec->get_pz());
  _E.push_back(truth_sec->get_e());
  */

  /*
    auto iter_ = std::find(primary_id.begin(), primary_id.end(), primary_id_sec);
			
    if (iter_ != primary_id.end()) // if we find a match
    {
    int pos_ = std::distance(primary_id.begin(), iter_); // this gives the index (position) where we find a match 
			
    // updating everything in vector (before ntuple)
    _PID_PRI.push_back(pid[pos_]); // for primary particle
    _PID_SEC.push_back(truth_sec->get_pid()); // for secondary particle
    _Eta.push_back(getEta(truth_sec)); // updating pseudorapidity of secondary particle
    _PX.push_back(truth_sec->get_px());
    _PY.push_back(truth_sec->get_py());
    _PZ.push_back(truth_sec->get_pz());
    _E.push_back(truth_sec->get_e());

    }
  */
  /*
    }

    _PID_PRI.resize(_E.size(), 0); // resize to the number of 4 vectors and initialize with zero

    if (_PID_PRI.size() != _primary_id.size()) {std::cout << "there is mismatch in the size" << std::endl;}
	
  */
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // get primary particle
  PHG4TruthInfoContainer::Range truthRange = truthinfo->GetPrimaryParticleRange();
  PHG4TruthInfoContainer::ConstIterator truthIter;

  PHG4Particle *truth_pri;
	
  /*
    std::vector<int> primary_id;
    std::vector<int> pid;
	
    primary_id.clear();
    pid.clear();
  */

  for(truthIter = truthRange.first; truthIter != truthRange.second; truthIter++)
    {
      truth_pri = truthIter->second;

      // do not accept if pseudorapidity falls outside the range of detector (sPHENX) acceptance
      if (fabs(getEta(truth_pri)) > 1.1)  continue;
			
      if ((truth_pri->get_parent_id()) != 0)
	{
	  std::cout << "Secondary Particle Spotted" << std::endl;
	  continue; 
	}

      /*
	int primary_id_pri = truth_pri->get_primary_id(); // get primary id of primary particle

	// Iterate through _primary_id to find matches between primary and secondary particles and update _PID_PRI accordingly
	for (size_t i = 0; i < _primary_id.size(); ++i) 
	{
        if (_primary_id[i] == primary_id_pri) {_PID_PRI[i] = truth_pri->get_pid();} // update "_PID_PRI"
	}
      */
      /*			 
				 primary_id.push_back(truth_pri->get_primary_id());
				 pid.push_back(truth_pri->get_pid());
      */

      _primary_id.push_back(truth_pri->get_primary_id());
      _PID_PRI.push_back(truth_pri->get_pid()); // for primary particle
      _PID_SEC.push_back(0);
      _Eta.push_back(getEta(truth_pri)); // updating pseudorapidity of primary particle
      _PX.push_back(truth_pri->get_px());
      _PY.push_back(truth_pri->get_py());
      _PZ.push_back(truth_pri->get_pz());
      _E.push_back(truth_pri->get_e());


    }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // updating in TTree
  _eventNumber = n_event;
  _nFourVector = _E.size();

  // if the size of vector is more than 600K (our limitation of the vector) we ignore that event
  if (_E.size() > 10000) {std::cout << " The size of vector is = " << _E.size() << " . " << std::endl;}
  if (_E.size() < 1) {std::cout << "Missing data in event number = " << n_event << std::endl;}
	
  // updating Ttree
  if ((_E.size() > 0) && (_E.size() < 10000))
    {
      for (int kk = 0; kk < _nFourVector; kk++)
	{
	  // update all the branches of the Ttree
	  // if size of four vector is zero then, we expect there is nothing inside the branches of the Ttree
	  _fv_energy[kk] = _E[kk];
	  _fv_px[kk] = _PX[kk];
	  _fv_py[kk] = _PY[kk];
	  _fv_pz[kk] = _PZ[kk];
	  _fv_Eta[kk] = _Eta[kk];
	  _pid_primary[kk] = _PID_PRI[kk];
	  _pid_secondary[kk] = _PID_SEC[kk];

	  //std::cout << _PA_ID[kk] << " , " << _PZ[kk] << std::endl;
	  //std::cout << _Parent_ID[kk] << " , " << _fv_pz[kk] << std::endl;
	}
    }
		  
  n_event++; // updating event number
  _eventTree->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int TruthParticleAna::ResetEvent(PHCompositeNode *topNode)
{
  //std::cout << "TruthParticleAna::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TruthParticleAna::EndRun(const int runnumber)
{
  std::cout << "TruthParticleAna::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int TruthParticleAna::End(PHCompositeNode *topNode)
{
  std::cout << "TruthParticleAna::End(PHCompositeNode *topNode) This is the End" << std::endl;

  if (topNode == 0 && f_temp)
    {
      out->Close();
      f_temp->Close();
      delete f_temp;
      delete out;
      return Fun4AllReturnCodes::EVENT_OK;
    }

  out->cd();
  //_data->Write(); 
  out->Write();
  out->Close();
  delete out; 

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TruthParticleAna::Reset(PHCompositeNode *topNode)
{
  std::cout << "TruthParticleAna::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void TruthParticleAna::Print(const std::string &what) const
{
  std::cout << "TruthParticleAna::Print(const std::string &what) const Printing info for " << what << std::endl;
}
//____________________________________________________________________________.. 
float TruthParticleAna::getEta(PHG4Particle *particle)
{
  float px = particle -> get_px();
  float py = particle -> get_py();
  float pz = particle -> get_pz();
  float p = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));

  return 0.5*log((p+pz)/(p-pz));
}


//____________________________________________________________________________.. 
float TruthParticleAna::FourVecEta(TLorentzVector& four_vector)
{
  float px = four_vector.Px();
  float py = four_vector.Py();
  float pz = four_vector.Pz();
  float p = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));

  return 0.5*log((p+pz)/(p-pz));
}


//______________________________________________________________________________..
// this function can be used to track among primary and primary particles and analyze them
void TruthParticleAna::check_primary_secondary(int nevts, TString _filename, TTree * intree)
{

  std::cout << "starting loop to read ntuples and make the list of all particles among primary and secondary particle range" << std::endl;

  TTree * t1 = intree;
  if (!intree)
    {
      TFile *f = new TFile(_filename);
      t1 = (TTree *) f->Get("_eventTree");
    }
  
  // Set Branches
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nFourVector", &_nFourVector);
  t1->SetBranchAddress("_fv_energy", _fv_energy);
  t1->SetBranchAddress("_fv_px", _fv_px);
  t1->SetBranchAddress("_fv_py", _fv_py);
  t1->SetBranchAddress("_fv_pz", _fv_pz);
  t1->SetBranchAddress("_fv_Eta", _fv_Eta);
  t1->SetBranchAddress("_primary_id", _primary_id);
  t1->SetBranchAddress("_pid_primary", _pid_primary);
  t1->SetBranchAddress("_pid_secondary", _pid_secondary);
  t1->SetBranchAddress("_embedding", _embedding);



  int nEntries = (int) t1->GetEntries();
  int nevts2 = nevts;

  if (nevts < 0 || nEntries < nevts)
    nevts2 = nEntries;

  int counter = 0;
	
  //for (int i = 0; i < 1; i++)
  for (int i = 0; i < nevts2; i++)
    {
      if (i%10000 == 0) {std::cout << "Event Number = " << i << std::endl;}
      // load the ith instance of the TTree
      t1->GetEntry(i);

      //std::cout << "----------------------------------------------------------------------------------------------------------" << std::endl;
      //std::cout << "----------------------------------------------------------------------------------------------------------" << std::endl;
      //std::cout << "running event number = " << i << std::endl;
		
      for (int j = 0; j < _nFourVector; j++)
	{
	  if (_embedding[j] != 0) continue;
			 
	  if (abs(_pid_primary[j]) == 11) {counter += 1;}
			
	  //if ((abs(_pid_primary[j]) == 411) || (abs(_pid_primary[j]) == 421) || (abs(_pid_primary[j]) == 431)) {counter += 1;}
	  //if (_primary_id[j] != _parent_id[j]) continue;

	  //if (_primary_id[j] != 7398) continue;
	  //std::cout << " primary id , parent id , secondary pid = " << _primary_id[j] << " , " << _parent_id[j] << " , " << _pid_secondary[j] << " . " << std::endl;
	}
    } // reading every event

  std::cout << "Number of particle of our concern as primary particle = " << counter << std::endl;
	
}

//______________________________________________________________________________..
void TruthParticleAna::Loop(int nevts, TString _filename, TTree * intree)
{
  // set of parameters used for different selection criteria (cuts)
  int ncluster_cut = 2000; // cut to limit the maximum number of clusters
  /*
    float first_pt_cut = 1.0; // first photon pt cut
    float second_pt_cut = 0.6; // second photon pt cut
    float alpha_cut = 0.6; // energy asymetry cut
    float delR_cut = 1.1; // cone cut
  */
  int min_FV_cut = 0; // some peripheral events we are ignoring


  std::cout << "starting loop to make four vector and hence reconstructing invariant mass" << std::endl;

  TTree * t1 = intree;
  if (!intree)
    {
      TFile *f = new TFile(_filename);
      t1 = (TTree *) f->Get("_eventTree");
    }
  
  // Set Branches
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nFourVector", &_nFourVector);
  t1->SetBranchAddress("_fv_energy", _fv_energy);
  t1->SetBranchAddress("_fv_px", _fv_px);
  t1->SetBranchAddress("_fv_py", _fv_py);
  t1->SetBranchAddress("_fv_pz", _fv_pz);
  t1->SetBranchAddress("_fv_Eta", _fv_Eta);
  t1->SetBranchAddress("_primary_id", _primary_id);
  t1->SetBranchAddress("_pid_primary", _pid_primary);
  t1->SetBranchAddress("_pid_secondary", _pid_secondary);
  t1->SetBranchAddress("_embedding", _embedding);

  // make TLorentz vector to save 4 vector data
  std::vector<TLorentzVector> saveFV;

  //  int nEntries = (int) t1->GetEntriesFast();
  int nEntries = (int) t1->GetEntries();
  int nevts2 = nevts;

  if (nevts < 0 || nEntries < nevts)
    nevts2 = nEntries;

  for (int i = 0; i < nevts2; i++)
    {
      // load the ith instance of the TTree
      t1->GetEntry(i);

      if (i % 10000 == 0) {std::cout << "event number = " << i << std::endl;}

      int nFourVector = _nFourVector;

      if (nFourVector < min_FV_cut) continue; // if the count of 4-vector is not 2 (min) then we can not make invariant mass plot.

      if (nFourVector > ncluster_cut) continue;
		
      saveFV.clear(); // clearing fourvector for new start

      for (int j = 0; j < nFourVector; j++)
	{		 
	  TLorentzVector LV_d;

	  LV_d.SetPxPyPzE(_fv_px[j], _fv_py[j], _fv_pz[j], _fv_energy[j]);
	  saveFV.push_back(LV_d);
	}
    
				
      for (int jCs = 0; jCs < nFourVector; jCs++) // outer loop
	{
	  TLorentzVector& pho1 = saveFV[jCs];
			
	  // combined use of embedding and pid_primary will control what we are looking in primary particle range
	  if (_embedding[jCs] != 0) continue; // use either of embedded or non embedded particle (0 is hijing only, 1 is embedding only)
	  if (_pid_primary[jCs] == 111) continue; 
	  //if (_pid_secondary[jCs] != 22) continue; // use only photons

	  int pri1 = _pid_primary[jCs];
			
	  //if (fabs(pho1.Pt()) < first_pt_cut)	continue;

	  // etabin for leading cluster
	  //float eta1 =  _fv_Eta[jCs];
			
	  // another loop to go into the saved cluster
	  // we have removed the situation for double counting
	  for (int kCs = jCs+1; kCs < nFourVector; kCs++)
	    {
	      TLorentzVector& pho2 = saveFV[kCs];
				
	      // combined use of embedding and pid_primary will control what we are looking in primary particle range
	      if (_embedding[kCs] != 0) continue; // use either of embedded or non embedded particle (0 is hijing only, 1 is embedding only)
	      //if (fabs(_pid_primary[kCs]) > 500) continue; // only use eta meson in primary particle 
	      //if (_pid_secondary[kCs] != 22) continue; // use only photons
				
	      //if (_pid_primary[kCs] != 22) continue;

	      if (_pid_primary[kCs] == 111) continue; 
	      //if (abs(_pid_primary[kCs]) != 321) continue; 

	      int pri2 = _pid_primary[kCs];
			
	      //if (fabs(pho2.Pt()) < second_pt_cut) continue;

	      // etabin for leading cluster
	      //float eta2 =  _fv_Eta[kCs];
	      float _deleta = pho1.Eta() - pho2.Eta();
	      float _delphi = pho1.Phi() - pho2.Phi();
				
	      if (_delphi > M_PI) {_delphi -= 2 * M_PI;}       
	      if (_delphi < -1 * M_PI) {_delphi += 2 * M_PI;}  
				
	      float _delR = TMath::Sqrt((_deleta * _deleta) + (_delphi * _delphi));
				
			
	      //alphaCut = fabs((pho1.E() - pho2.E())/(pho1.E()+ pho2.E()));
	      //if (alphaCut > alpha_cut) continue;

	      //if (pho1.DeltaR(pho2) > delR_cut) continue;

	      TLorentzVector pi0lv;
	      pi0lv = pho1 + pho2;
				
	      //if (pi0lv.M() < 0.02) continue;

	      Pri1_Pri2_pairpT->Fill(pri1, pri2, pi0lv.Pt());
	  		 
	      pairInvMassTotal->Fill(pi0lv.M());
	      pairInvMassPtEta->Fill(pi0lv.M(), pi0lv.Pt(), pi0lv.Eta());


	      delR_f->Fill(pho1.DeltaR(pho2));
	      delR_b->Fill(_delR);

	      delphi_f->Fill(_delphi);
	      deleta_f->Fill(_deleta);
				
	      //Pair_PX_dist->Fill(pi0lv.X());
	      //Pair_PY_dist->Fill(pi0lv.Y());
	      //PXPYPT->Fill(pi0lv.X(), pi0lv.Y(), pi0lv.Pt());
				
	      /*
		pairInvMassPt1Pt2->Fill(pi0lv.M(), pho1.Pt(), pho2.Pt());
		pairInvMassPX1PX2->Fill(pi0lv.M(), pho1.Px(), pho2.Px());
		pairInvMassPY1PY2->Fill(pi0lv.M(), pho1.Py(), pho2.Py());
		pairInvMassPZ1PZ2->Fill(pi0lv.M(), pho1.Pz(), pho2.Pz());
		pairInvMassE1E2->Fill(pi0lv.M(), pho1.E(), pho2.E());
		pairInvMassEta1Eta2->Fill(pi0lv.M(), eta1, eta2);
	      */
	  		 
	    } // inner loop for clusters (truth particle)
	} // outer loop for clusters (truth particle)
    } // reading every event
}

//______________________________________________________________________________..
void TruthParticleAna::Loop_background_position_swappping(int nevts, TString _filename, TTree * intree)
{
  // set of parameters used for different selection criteria (cuts)
  int ncluster_cut = 700; // cut to limit the maximum number of clusters
  float first_pt_cut = 1.0; // first photon pt cut
  float second_pt_cut = 0.6; // second photon pt cut
  float alpha_cut = 0.6; // energy asymetry cut
  float delR_cut = 1.1; // cone cut
  int min_FV_cut = 50; // min number of four vector we are keeping

  std::cout << "starting loop to make four vector and hence reconstructing invariant mass" << std::endl;

  TTree * t1 = intree;
  if (!intree)
    {
      TFile *f = new TFile(_filename);
      t1 = (TTree *) f->Get("_eventTree");
    }
  
  // Set Branches
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nFourVector", &_nFourVector);
  t1->SetBranchAddress("_fv_energy", _fv_energy);
  t1->SetBranchAddress("_fv_px", _fv_px);
  t1->SetBranchAddress("_fv_py", _fv_py);
  t1->SetBranchAddress("_fv_pz", _fv_pz);
  t1->SetBranchAddress("_fv_Eta", _fv_Eta);
  t1->SetBranchAddress("_primary_id", _primary_id);
  t1->SetBranchAddress("_pid_primary", _pid_primary);
  t1->SetBranchAddress("_pid_secondary", _pid_secondary);
  t1->SetBranchAddress("_embedding", _embedding);


  //  int nEntries = (int) t1->GetEntriesFast();
  int nEntries = (int) t1->GetEntries();
  int nevts2 = nevts;

  if (nevts < 0 || nEntries < nevts)
    nevts2 = nEntries;
	
  for (int i = 0; i < nevts2; i++)
    {
      // load the ith instance of the TTree
      t1->GetEntry(i);

      if (i % 1000 == 0) {std::cout << "event number = " << i << std::endl;}

      int nFourVector = _nFourVector;

      if (nFourVector < min_FV_cut) continue; // if the count of 4-vector is not "min_FV_cut" then we ignore that event

      if (nFourVector > ncluster_cut) continue; // if 4 vector are more than this threshold we will ignore the event

      int nloop_bck = 25 + (nFourVector - min_FV_cut)/10; 
	
      // make TLorentz vector to save 4 vector data
      std::vector<TLorentzVector> saveFV;
      std::vector<int> secondary_pid;

      saveFV.clear(); // clearing fourvector for new start
      secondary_pid.clear();

      for (int j = 0; j < nFourVector; j++)
	{		 
	  TLorentzVector LV_d;

	  LV_d.SetPxPyPzE(_fv_px[j], _fv_py[j], _fv_pz[j], _fv_energy[j]);
	  saveFV.push_back(LV_d);

	  secondary_pid.push_back(_pid_secondary[j]);
	}
    
      // outer loop to iterate over cluster (get first photon)		
      for (int jCs = 0; jCs < nFourVector; jCs++) // outer loop
	{
	  TLorentzVector& pho1 = saveFV[jCs];
			
	  // combined use of embedding and pid_primary will control what we are looking in primary particle range
	  if (_embedding[jCs] != 0) continue; // use either of embedded or non embedded particle (0 is hijing only, 1 is embedding only)
	  //if (fabs(_pid_primary[jCs]) > 500) continue; // only use eta meson in primary particle 
	  //if (_pid_secondary[jCs] != 22) continue; // use only photons
	  //if ((_pid_primary[jCs] != 310) && (fabs(_pid_primary[jCs]) != 321) && (_pid_primary[jCs] != 130)) continue; 

	  if (fabs(pho1.Pt()) < first_pt_cut)	continue;

			
	  // we will have two inner loop to ensure one will make foreground and another will make background			
	  // this inner loop will make foreground
	  for (int kCs = jCs+1; kCs < nFourVector; kCs++)
	    {
	      TLorentzVector& pho2 = saveFV[kCs];
				
	      // combined use of embedding and pid_primary will control what we are looking in primary particle range
	      if (_embedding[kCs] != 0) continue; // use either of embedded or non embedded particle (0 is hijing only, 1 is embedding only)
	      //if (fabs(_pid_primary[kCs]) > 500) continue; // only use eta meson in primary particle 
	      //if (_pid_secondary[kCs] != 22) continue; // use only photons
	      //if ((_pid_primary[kCs] != 310) && (fabs(_pid_primary[kCs]) != 321) && (_pid_primary[kCs] != 130)) continue; 
				
	      if (fabs(pho2.Pt()) < second_pt_cut) continue;
				
	      alphaCut = fabs((pho1.E() - pho2.E())/(pho1.E()+ pho2.E()));
	      if (alphaCut > alpha_cut) continue;
			
	      if (pho1.DeltaR(pho2) > delR_cut) continue;

	      float _delphi = pho1.Phi() - pho2.Phi();
	      float _deleta = pho1.Eta() - pho2.Eta();
       
	      // getting the 4-vector for photon pair
	      TLorentzVector mesonlv= pho1 + pho2;
					
	      if((mesonlv.M() < 0.02) || (mesonlv.Pt() < 1.0)) continue;
						
	      pairInvMassTotal->Fill(mesonlv.M());
	      pairInvMassPtEtaBkgd->Fill(mesonlv.M(), mesonlv.Pt(), mesonlv.Eta());
				
	      Pair_PT_f->Fill(mesonlv.Pt());
	      delphi_f->Fill(_delphi);
	      deleta_f->Fill(_deleta);
	    }
			
	  // Now we are doing inner loop for making background via position swapping						
	  TLorentzVector dummy_LV; // this is created for 4-vector in swapping
	  dummy_LV.SetPtEtaPhiE(0, 0, 0, 0);

	  // second inner loop for background (background by position swapping)
	  for (int lCs = jCs+1; lCs < nFourVector; lCs++)
	    {
	      TLorentzVector& pho3 = saveFV[lCs];
				
	      // combined use of embedding and pid_primary will control what we are looking in primary particle range
	      if (_embedding[lCs] != 0) continue; // use either of embedded or non embedded particle (0 is hijing only, 1 is embedding only)
	      //if (fabs(_pid_primary[lCs]) > 500) continue; // only use eta meson in primary particle 
	      //if (_pid_secondary[lCs] != 22) continue; // use only photons
	      //if ((_pid_primary[lCs] != 310) && (fabs(_pid_primary[lCs]) != 321) && (_pid_primary[lCs] != 130)) continue; 

				
	      if (fabs(pho3.Pt()) < second_pt_cut) continue;
				
	      alphaCut = fabs((pho1.E() - pho3.E())/(pho1.E()+ pho3.E()));
	      if (alphaCut > alpha_cut) continue;
							
				
	      if (pho1.DeltaR(pho3) > delR_cut) continue; // we might be able to use this cone cut only for high pT clusters not for lower one
				
       
	      // we will first randomly select the third cluster and extract the "z co-ordinate" from third cluster and swap with second cluster
	      std::vector<int> third_clus = {jCs, lCs};
				
	      for (int nloop = 0; nloop < nloop_bck; nloop++) // starting the loop that contributes to background
		{
		  int mCs = lCs;
		  do{
		    TRandom3 randgen(0);
		    mCs = randgen.Integer(nFourVector);
		  } while (std::find(third_clus.begin(), third_clus.end(), mCs) != third_clus.end());

		  third_clus.push_back(mCs); //push back the randomly generated cluster number

		  if (_embedding[mCs] != 0) continue; // use either of embedded or non embedded particle (0 is hijing only, 1 is embedding only)
		  //if (secondary_pid[mCs] != 22) continue; // test of which (either of photon/lepton) we are using
		  //if (primary_pid[mCs] != 221) continue; // ignore if primary particle not eta-meosn						
				
		  TLorentzVector &pho4 = saveFV[mCs];
					
					
		  if (fabs(pho4.Pt()) < second_pt_cut) continue;
					
		  float alphaCut_4 = fabs((pho1.E() - pho4.E())/(pho1.E()+ pho4.E()));
		  if (alphaCut_4 > alpha_cut) continue;
				

		  if (pho1.DeltaR(pho4) > delR_cut) continue;

		  //if (abs(pho3.Eta() - pho4.Eta()) > 0.6) continue;

		  //if (abs((pho3.M()- pho4.M())/(pho3.M()+ pho4.M())) > 0.1) continue; // not taking anything less than 10 percent difference
		  //if (abs((pho3.Pt()- pho4.Pt())/(pho3.Pt()+ pho4.Pt())) > 0.1) continue; // not taking anything less than 10 percent difference

		  TLorentzVector &pho_bck = dummy_LV;
					
		  // "pho3" and "pho4" are photons for background
		  //pho_bck.SetPtEtaPhiE(pho4.Pt(), pho3.Eta(), pho3.Phi(), pho4.E()); // this is the updated second cluster
		  pho_bck.SetXYZM(pho4.X(), pho4.Y(), pho3.Z(), pho3.M()); // this is the updated second cluster
					
		  //pho_bck.SetXYZM(pho4.X(), pho4.Y(), pho3.Z(), pho3.M()); // this is the updated second cluster					
					
		  //pho_bck.SetPtEtaPhiE(pho3.Pt(), pho4.Eta(), pho4.Phi(), pho3.E()); // this is the updated second cluster					
		  //pho_bck.SetPtEtaPhiE(pho4.Pt(), pho3.Eta(), pho3.Phi(), pho4.E()); // this is the updated second cluster					
					
					
		  if (fabs(pho_bck.Pt()) < second_pt_cut) continue;
					
		  float alphaCut_bck = fabs((pho1.E() - pho_bck.E())/(pho1.E()+ pho_bck.E()));
		  if (alphaCut_bck > alpha_cut) continue;
					
		  if (pho1.DeltaR(pho_bck) > delR_cut) continue;
					

		  // getting the 4-vector for eta-meson
		  TLorentzVector mesonlv_bck= pho1 + pho_bck;
					
		  float _delphi_b = pho1.Phi() - pho4.Phi();
		  float _deleta_b = pho1.Eta() - pho4.Eta();
					
		  if((mesonlv_bck.M() < 0.02) || (mesonlv_bck.Pt() < 1.0)) continue;
						
		  pairInvMassTotalBkgd->Fill(mesonlv_bck.M()); // fill the tower by tower histograms with invariant mass
		  pairInvMassPtEtaBkgd->Fill(mesonlv_bck.M(), mesonlv_bck.Pt(), mesonlv_bck.Eta());
				
		  Pair_PT_b->Fill(mesonlv_bck.Pt());
					
		  delphi_b->Fill(_delphi_b);
		  deleta_b->Fill(_deleta_b);
															
		} // loop for background (event mixing)
			 
	    } // inner loop for clusters/ truth particles (for background)
	} // outer loop for clusters (or truth particle) covering all (both foreground and background)

    } // reading every event

}


//______________________________________________________________________________..
void TruthParticleAna::Loop_background_event_mixing(int nevts, TString _filename, TString ref_file, TTree * intree)
{
  // set of parameters used for different selection criteria (cuts)
  int ncluster_cut = 2000; // cut to limit the maximum number of clusters
  float first_pt_cut = 1.0; // initialize first photon pt cut
  float second_pt_cut = 0.6; // initialize second photon pt cut
  float alpha_cut = 0.6; // energy asymetry cut
  float delR_cut = 1.1; // cone cut
  int bckgnd_evnts = 50; // number of background events we are interested in

  std::cout << "starting loop to construct foreground and background using event mixing method" << std::endl;
	
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Open ref root file that has histogram for calculating weight
  TFile *file = TFile::Open(ref_file, "UPDATE");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: cannot open file " << ref_file << std::endl;
    return;}

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  // Retrieve the 3D histogram for pairpT_DelEta_DelPhi
  TH3F *fg_w = dynamic_cast<TH3F*>(file->Get("pairpTDelPhiDelEta"));
  TH3F *bg_w = dynamic_cast<TH3F*>(file->Get("pairpTDelPhiDelEtaBkgd"));

  if ((!fg_w) || (!bg_w)) {
    std::cerr << "Error: cannot find one or more 3D histogram in file." << std::endl;
    file->Close();
    return;}
	
  fg_w->Divide(bg_w);

	
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Retrieve two different kind of 2D histogram for "DelEta_pairpT" and "DelPhi_pairpT"
  // first we do delphi and pairpT
  TH2F *fg_phi_pT = dynamic_cast<TH2F*>(file->Get("DelPhi_PairpT_f"));
  TH2F *bg_phi_pT = dynamic_cast<TH2F*>(file->Get("DelPhi_PairpT_b"));

  if ((!fg_phi_pT) || (!bg_phi_pT)) {
    std::cerr << "Error: cannot find one or more 2D histogram in file." << std::endl;
    file->Close();
    return;}
	
  fg_phi_pT->Divide(bg_phi_pT);

  // second we do deleta and pairpT
  TH2F *fg_eta_pT = dynamic_cast<TH2F*>(file->Get("DelEta_PairpT_f"));
  TH2F *bg_eta_pT = dynamic_cast<TH2F*>(file->Get("DelEta_PairpT_b"));

  if ((!fg_eta_pT) || (!bg_eta_pT)) {
    std::cerr << "Error: cannot find one or more 2D histogram in file." << std::endl;
    file->Close();
    return;}
	
  fg_eta_pT->Divide(bg_eta_pT);

  // third we do deleta and delphi
  TH2F *fg_phi_eta = dynamic_cast<TH2F*>(file->Get("DelPhiDelEta_f"));
  TH2F *bg_phi_eta = dynamic_cast<TH2F*>(file->Get("DelPhiDelEta_b"));

  if ((!fg_phi_eta) || (!bg_phi_eta)) {
    std::cerr << "Error: cannot find one or more 2D histogram in file." << std::endl;
    file->Close();
    return;}
	
  fg_phi_eta->Divide(bg_phi_eta);
	

  // fourth we do delR and pairpT
  TH2F *fg_delR_pt = dynamic_cast<TH2F*>(file->Get("DelR_pairpT_f"));
  TH2F *bg_delR_pt = dynamic_cast<TH2F*>(file->Get("DelR_pairpT_b"));

  if ((!fg_delR_pt) || (!bg_delR_pt)) {
    std::cerr << "Error: cannot find one or more 2D histogram in file." << std::endl;
    file->Close();
    return;}
	
  fg_delR_pt->Divide(bg_delR_pt);
	
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  // Retrieve the 1D histogram for pairpT only
  TH1F *fg_pT = dynamic_cast<TH1F*>(file->Get("Pair_PT_f"));
  TH1F *bg_pT = dynamic_cast<TH1F*>(file->Get("Pair_PT_b"));

  if ((!fg_pT) || (!bg_pT)) {
    std::cerr << "Error: cannot find one or more 1D histogram in file." << std::endl;
    file->Close();
    return;}
	
  fg_pT->Divide(bg_pT);

  TH1F *fg_delphi = dynamic_cast<TH1F*>(file->Get("del_phi_f"));
  TH1F *bg_delphi = dynamic_cast<TH1F*>(file->Get("del_phi_b"));

  if ((!fg_delphi) || (!bg_delphi)) {
    std::cerr << "Error: cannot find one or more 1D histogram in file." << std::endl;
    file->Close();
    return;}
	
  fg_delphi->Divide(bg_delphi);

  TH1F *fg_deleta = dynamic_cast<TH1F*>(file->Get("del_eta_f"));
  TH1F *bg_deleta = dynamic_cast<TH1F*>(file->Get("del_eta_b"));

  if ((!fg_deleta) || (!bg_deleta)) {
    std::cerr << "Error: cannot find one or more 1D histogram in file." << std::endl;
    file->Close();
    return;}
	
  fg_deleta->Divide(bg_deleta);


  TH1F *fg_delR = dynamic_cast<TH1F*>(file->Get("delR_f"));
  TH1F *bg_delR = dynamic_cast<TH1F*>(file->Get("delR_b"));

  if ((!fg_delR) || (!bg_delR)) {
    std::cerr << "Error: cannot find one or more 1D histogram in file." << std::endl;
    file->Close();
    return;}
	
  fg_delR->Divide(bg_delR);
	
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  TTree * t1 = intree;
  if (!intree)
    {
      TFile *f = new TFile(_filename);
      t1 = (TTree *) f->Get("_eventTree");
    }  
	
  // Set Branches
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nFourVector", &_nFourVector);
  t1->SetBranchAddress("_fv_energy", _fv_energy);
  t1->SetBranchAddress("_fv_px", _fv_px);
  t1->SetBranchAddress("_fv_py", _fv_py);
  t1->SetBranchAddress("_fv_pz", _fv_pz);
  t1->SetBranchAddress("_fv_Eta", _fv_Eta);
  t1->SetBranchAddress("_primary_id", _primary_id);
  t1->SetBranchAddress("_pid_primary", _pid_primary);
  t1->SetBranchAddress("_pid_secondary", _pid_secondary);
  t1->SetBranchAddress("_embedding", _embedding);


  //  int nEntries = (int) t1->GetEntriesFast();
  int nEntries = (int) t1->GetEntries();
  int nevts2 = nevts;

  if (nevts < 0 || nEntries < nevts) {nevts2 = nEntries;}
	

  for (int i = 0; i < nevts2; i++) // starting to load Ttree ntuple
    {
      // load the ith instance of the TTree
      t1->GetEntry(i);

      if (i % 1000 == 0) {std::cout << "event number = " << i << std::endl;}

      int nFourVector1 = _nFourVector;

      if (nFourVector1 < 10) continue; 

      if (nFourVector1 > ncluster_cut) continue;

      // make TLorentz vector to save 4 vector data
      std::vector<TLorentzVector> saveFV1;
      std::vector<int> primary_pid1;
      std::vector<int> secondary_pid1;
      std::vector<int> embedding1;

      saveFV1.clear(); // clearing fourvector for new start
      primary_pid1.clear();
      secondary_pid1.clear();
      embedding1.clear();

      // total energy of all four-vector
      float TE_1 = 0;

      for (int j = 0; j < nFourVector1; j++)
	{		 
	  TLorentzVector LV_d1;

	  LV_d1.SetPxPyPzE(_fv_px[j], _fv_py[j], _fv_pz[j], _fv_energy[j]);
	  saveFV1.push_back(LV_d1);

	  // now push back pid of primary and secondary particles
	  primary_pid1.push_back(_pid_primary[j]);
	  secondary_pid1.push_back(_pid_secondary[j]);
	  embedding1.push_back(_embedding[j]);
		 
	  TE_1 +=  _fv_energy[j]; // adding up total energy (for first event)
	}
    
      // getting into another event, if two events are same then it will make foreground
      // but if two events are different then it will make background

      for (int i2 = i - bckgnd_evnts; i2 < i + bckgnd_evnts; i2++)
	{
	  if ((i2 < 0) || (i2 > nevts2)) continue; //discard negative number events and more than it exists

	  t1->GetEntry(i2);

	  int nFourVector2 = _nFourVector;

	  if (nFourVector2 < 10) continue; // if the count of 4-vector is not 2 (min) then we can not make invariant mass plot.

	  if (nFourVector2 > ncluster_cut) continue;

	  if (abs(nFourVector1 - nFourVector2) > 50) continue; // ignore events which are extremely different

	  // make TLorentz vector to save 4 vector data
	  std::vector<TLorentzVector> saveFV2;
	  std::vector<int> primary_pid2;
	  std::vector<int> secondary_pid2;
	  std::vector<int> embedding2;

	  saveFV2.clear(); // clearing fourvector for new start
	  primary_pid2.clear();
	  secondary_pid2.clear();
	  embedding2.clear();

	  // total energy of all four-vector
	  float TE_2 = 0;

	  for (int k = 0; k < nFourVector2; k++)
	    {		 
	      TLorentzVector LV_d2;

	      LV_d2.SetPxPyPzE(_fv_px[k], _fv_py[k], _fv_pz[k], _fv_energy[k]);
	      saveFV2.push_back(LV_d2);

	      // now push back pid of primary and secondary particles
	      primary_pid2.push_back(_pid_primary[k]);
	      secondary_pid2.push_back(_pid_secondary[k]);
	      embedding2.push_back(_embedding[k]);

	      TE_2 +=  _fv_energy[k]; // adding up total energy (for first event)
	    }

	  // trying to avoid two different events with similar kinematics
	  // this is event level cuts 
	  // we do not have z_vertex cut but it's nice to test with omitting two extremely close events
	  if ((i != i2) && ((((TE_1 - TE_2)/ (TE_1 + TE_2)) < 0.01) || ((nFourVector1 - nFourVector2) < 2))) continue;

	  for (int ev1 = 0; ev1 < nFourVector1; ev1++) // outer loop for events
	    {
	      TLorentzVector& pho1 = saveFV1[ev1];

	      if(fabs(pho1.Pt()) < first_pt_cut) continue;

	      if (primary_pid1[ev1] == 22) continue; // check which primary particle we want to use
	      if (embedding1[ev1] != 0) continue; // use either of embedded or non embedded particle (0 is hijing only, 1 is embedding only)
	      //if (secondary_pid1[ev1] != 22) continue; // test of which (either of photon/lepton) we are using

				 
	      for (int ev2 = ev1; ev2 < nFourVector2; ev2++) // inner loop for events
		{

		  TLorentzVector& pho2 = saveFV2[ev2];

		  if(fabs(pho2.Pt()) < second_pt_cut) continue;
		  //if(fabs(pho2.Pt()) < (second_pt_cut - 0.1*(nFourVector2 - 600)/2000)) continue;
						
		  if (primary_pid2[ev2] == 22) continue; // check which primary particle we want to use 						
		  if (embedding2[ev2] != 0) continue; // use either of embedded or non embedded particle (0 is hijing only, 1 is embedding only)
		  //if (secondary_pid2[ev2] != 22) continue; // test of which (either of photon/lepton) we are using

		  alphaCut = fabs((pho1.E() - pho2.E())/(pho1.E()+ pho2.E()));
		  if (alphaCut > alpha_cut) continue;	
								
		  if (pho1.DeltaR(pho2) > delR_cut) continue;

		  TLorentzVector mesonlv;
		  mesonlv = pho1 + pho2;
							
		  if ((mesonlv.M() < 0.02) || (mesonlv.Pt() < 1.0)) continue;
					
		  double _delphi = pho1.Phi() - pho2.Phi();
		  double _deleta = pho1.Eta() - pho2.Eta();
		  double _pairpT = mesonlv.Pt();
		  double _delR = pho1.DeltaR(pho2);

		  if (i == i2) // histogram for foreground
		    {
						 
		      pairInvMassTotal->Fill(mesonlv.M());
		      pairInvMassPtEta->Fill(mesonlv.M(), _pairpT, mesonlv.Eta());
						 
		      pairInvMassPtdelR->Fill(mesonlv.M(), _pairpT, _delR);
						 
		      delphi_f->Fill(_delphi);
		      deleta_f->Fill(_deleta);
		      delR_f->Fill(_delR);

		      Pair_PT_f->Fill(_pairpT);

		      pairpTDelPhiDelEta->Fill(_pairpT, _delphi, _deleta); // for foreground

		      DelPhi_PairpT_f->Fill(_delphi, _pairpT);
		      DelEta_PairpT_f->Fill(_deleta, _pairpT);						 
		      DelPhi_DelEta_f->Fill(_delphi, _deleta);
						 
		      DelR_pairpT_f->Fill(_delR, _pairpT);

		    }
						 						 
		  else if (i != i2) // histogram for background
		    {
						 
		      ////////////////////////////////////////////////////////////////////////////////////////////////////
						 
		      // While Using 3D Hist
		      // Find the bin corresponding to the phi and eta values
		      int BinPt = fg_w->GetXaxis()->FindBin(_pairpT);
		      int BinPhi = fg_w->GetYaxis()->FindBin(_delphi);
		      int BinEta = fg_w->GetZaxis()->FindBin(_deleta);

		      // Get the bin content (which will be weight)
		      double weight_3D = fg_w->GetBinContent(BinPt, BinPhi, BinEta);
						 
		      if (abs(weight_3D) > 100) {std::cout << "check weights" << std::endl;}						 

		      ////////////////////////////////////////////////////////////////////////////////////////////////////
		      // While Using 2D Hist
						 
		      // Find the bin corresponding to the pT and delphi values
		      int BinPhi1 = fg_phi_pT->GetXaxis()->FindBin(_delphi);
		      int BinpT1 = fg_phi_pT->GetYaxis()->FindBin(_pairpT);						 
		      // Get the bin content (which will be weight)
		      double weight1 = fg_phi_pT->GetBinContent(BinPhi1, BinpT1);

		      // Find the bin corresponding to the pT and deleta values
		      int BinEta2 = fg_eta_pT->GetXaxis()->FindBin(_deleta);
		      int BinpT2 = fg_eta_pT->GetYaxis()->FindBin(_pairpT);						 
		      // Get the bin content (which will be weight)
		      double weight2 = fg_eta_pT->GetBinContent(BinEta2, BinpT2);

		      // Find the bin corresponding to the delphi and deleta values
		      int BinPhi3 = fg_phi_eta->GetXaxis()->FindBin(_delphi);
		      int BinEta3 = fg_phi_eta->GetYaxis()->FindBin(_deleta);						 
		      // Get the bin content (which will be weight)
		      double weight3 = fg_phi_eta->GetBinContent(BinPhi3, BinEta3);
		      if ((abs(weight1) > 100) || (abs(weight2) > 100) || (abs(weight3) > 100)) {std::cout << "check weights" << std::endl;}
						
						
		      // Find the bin corresponding to the delR and pt values
		      int BindelR = fg_delR_pt->GetXaxis()->FindBin(_delR);
		      int BinpT4 = fg_delR_pt->GetYaxis()->FindBin(_pairpT);	

		      double weight4 = fg_delR_pt->GetBinContent(BindelR, BinpT4);
						

		      if(abs(weight4 * weight2 * weight3 * weight1 ) > 1e6){std::cout << "check weights" << std::endl;}
		      // Get the final weight from weight1 anf weight2
						 
		      ////////////////////////////////////////////////////////////////////////////////////////////////////
						 
		      // while using 1D Hist
		      int BinPairpT = fg_pT->GetXaxis()->FindBin(_pairpT);					 
		      double weight_pairpT = fg_pT->GetBinContent(BinPairpT);
						 
		      int BinDelPhi = fg_delphi->GetXaxis()->FindBin(_delphi);					 
		      double weight_delphi = fg_delphi->GetBinContent(BinDelPhi);

		      int BinDelEta = fg_deleta->GetXaxis()->FindBin(_deleta);					 
		      double weight_deleta = fg_deleta->GetBinContent(BinDelEta);
							
		      if ((abs(weight_pairpT) > 100) || (abs(weight_delphi) > 100) || (abs(weight_deleta) > 100)) {std::cout << "check weights" << std::endl;}
						 
		      //double weight = weight3 * weight11;
						 
		      int BinDelR = fg_delR->GetXaxis()->FindBin(_delR);					 
		      double weight_delR = fg_delR->GetBinContent(BinDelR);

		      if (abs(weight_delR) > 100) {std::cout << "check" << std::endl;}

		      //double weight = weight_pairpT * weight_delR;
						 
		      ////////////////////////////////////////////////////////////////////////////////////////////////////
					 	 
		      //double weight = weight_3D;
		      double weight;
						 
		      // Check for valid weight if we do not find valid one then assign it to be unity
		      if (std::isnan(weight) || std::isinf(weight)) {
			std::cout << "Nan or Inf weight found" << std::endl;
			weight = 1.0;
		      }

		      weight = 1.0;

		      pairInvMassTotalBkgd->Fill(mesonlv.M(), weight);
		      pairInvMassPtEtaBkgd->Fill(mesonlv.M(), _pairpT, mesonlv.Eta(), weight);

		      pairInvMassPtdelRBkgd->Fill(mesonlv.M(), _pairpT, _delR, weight);

		      delphi_b->Fill(_delphi, weight);
		      deleta_b->Fill(_deleta, weight);
		      delR_b->Fill(_delR, weight);
						 
		      Pair_PT_b->Fill(_pairpT, weight);
						 
		      pairpTDelPhiDelEtaBkgd->Fill(_pairpT, _delphi, _deleta, weight);	// for background	
						 
		      DelPhi_PairpT_b->Fill(_delphi, _pairpT, weight);
		      DelEta_PairpT_b->Fill(_deleta, _pairpT, weight);
		      DelPhi_DelEta_b->Fill(_delphi, _deleta, weight);
						 
		      DelR_pairpT_b->Fill(_delR, _pairpT, weight);
		    }
						
		} // inner loop closure 
	    } // outer loop closure

	} // event2 loop (outer-inner loop)

    } // event1 loop (outer-outer loop)
  file->Close();

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
// this version of position swapping can be used by exchanging PX and PY among two foreground clusters
//______________________________________________________________________________..
void TruthParticleAna::Loop_background_position_swappping_v2(int nevts, TString _filename, TTree * intree)
{
// set of parameters used for different selection criteria (cuts)
int ncluster_cut = 6000; // cut to limit the maximum number of clusters
float first_pt_cut = 2.0; // initialize first photon pt cut
float second_pt_cut = 1.0; // initialize second photon pt cut
float alpha_cut = 0.5; // energy asymetry cut
float delR_cut = 0.8; // cone cut

std::cout << "starting loop to construct foreground and background using position swapping method" << std::endl;

TTree * t1 = intree;
if (!intree)
{
TFile *f = new TFile(_filename);
t1 = (TTree *) f->Get("_eventTree");
}
  
// Set Branches
t1->SetBranchAddress("_eventNumber", &_eventNumber);
t1->SetBranchAddress("_nFourVector", &_nFourVector);
t1->SetBranchAddress("_fv_energy", _fv_energy);
t1->SetBranchAddress("_fv_px", _fv_px);
t1->SetBranchAddress("_fv_py", _fv_py);
t1->SetBranchAddress("_fv_pz", _fv_pz);
t1->SetBranchAddress("_fv_Eta", _fv_Eta);
//t1->SetBranchAddress("_pid_primary", _pid_primary);
t1->SetBranchAddress("_pid_secondary", _pid_secondary);


//  int nEntries = (int) t1->GetEntriesFast();
int nEntries = (int) t1->GetEntries();
int nevts2 = nevts;

if (nevts < 0 || nEntries < nevts)
nevts2 = nEntries;

// pre-loop to save all the clusters LorentzVector
//TLorentzVector *saveFV[10000];
//std::vector<std::unique_ptr<TLorentzVector>> saveFV;

for (int i = 0; i < nevts2; i++)
{
		
TLorentzVector *saveFV[10000];
// load the ith instance of the TTree
t1->GetEntry(i);

if (i % 1000 == 0) {std::cout << "event number = " << i << std::endl;}

int nFourVector = _nFourVector;
//std::cout << "four vector number = " << nFourVector << std::endl;

if (nFourVector < 5) continue; // if the count of 4-vector is not 2 (min) then we can not make invariant mass plot.

if (nFourVector > ncluster_cut) continue;

// save parent ID in a separate list
//std::vector<int> primary_pid(10000);
std::vector<int> secondary_pid(10000);

//primary_pid.clear();
secondary_pid.clear();

for (int j = 0; j < nFourVector; j++)
{
float px, py, pz, En;
En = _fv_energy[j];
px = _fv_px[j];
py = _fv_py[j];
pz = _fv_pz[j];

saveFV[j] = new TLorentzVector();
saveFV[j]->SetPxPyPzE(px, py, pz, En);
//saveFV.push_back(std::make_unique<TLorentzVector>(px, py, pz, En));

//primary_pid[j] = _pid_primary[j]; // save parent ID in separate vector
secondary_pid[j] = _pid_secondary[j]; // save p-ID in separate vector
}

//std::cout << p_id.size() << std::endl;

TLorentzVector *pho1, *pho2;

for (int jCs = 0; jCs < nFourVector; jCs++) // outer loop
{

pho1 = saveFV[jCs];
			
//if ((primary_pid[jCs] != 111) && (primary_pid[jCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
//if (primary_pid[jCs] != 221) continue; // ignore if primary particle is not eta-meosn
//if (secondary_pid[jCs] != 22) continue; // test of which (either of photon/lepton) we are using

//if (fabs(pho1->Pt()) < (first_pt_cut + (nFourVector - 600)/2000))	continue;
if (fabs(pho1->Pt()) < first_pt_cut)	continue;
			
// another loop to go into the saved cluster
// we have removed the situation for double counting
TLorentzVector dummy_LV;
dummy_LV.SetPtEtaPhiE(0, 0, 0, 0);
				
for (int kCs = jCs+1; kCs < nFourVector; kCs++)
{
//if (jCs == kCs) continue;

pho2 = saveFV[kCs];

//if ((primary_pid[kCs] != 111) && (primary_pid[kCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
//if (primary_pid[kCs] != 221) continue; // ignore if primary particle not eta-meosn
//if (secondary_pid[kCs] != 22) continue; // test of which (either of photon/lepton) we are using

//if (fabs(pho2->Pt()) < (second_pt_cut + (nFourVector - 600)/2000)) continue;
if (fabs(pho2->Pt()) < second_pt_cut ) continue;
				
alphaCut = fabs((pho1->E() - pho2->E())/(pho1->E()+ pho2->E()));
if (alphaCut > alpha_cut) continue;

if (pho1->DeltaR(*pho2) > delR_cut) continue;

//if (secondary_pid[lCs] != 22) continue; // test of which (either of photon/lepton) we are using
//if (primary_pid[lCs] != 221) continue; // ignore if primary particle not eta-meosn
						
TLorentzVector &pho1_c = dummy_LV;				

TLorentzVector &pho2_c = dummy_LV;
pho1_c.SetXYZM(pho2->X(), pho2->Y(), pho1->Z(), pho1->M()); // this is the updated first cluster
pho2_c.SetXYZM(pho1->X(), pho1->Y(), pho2->Z(), pho2->M()); // this is the updated second cluster		
				
if (fabs(pho1_c.Pt()) < first_pt_cut) continue; //second P_T cut
if (fabs(pho2_c.Pt()) < second_pt_cut) continue; //second P_T cut
	
float alphaCut_c1 = fabs((pho1->E() - pho1_c.E())/(pho1->E()+ pho1_c.E()));
if (alphaCut_c1 > alpha_cut) continue;
				
float alphaCut_c2 = fabs((pho2->E() - pho2_c.E())/(pho2->E()+ pho2_c.E()));
if (alphaCut_c2 > alpha_cut) continue;

if (pho1->DeltaR(pho1_c) > delR_cut) continue;
if (pho2->DeltaR(pho2_c) > delR_cut) continue;
				
// getting the 4-vector for eta-meson
TLorentzVector mesonlv_1 = *pho1 + pho2_c;
TLorentzVector mesonlv_2 = *pho2 + pho1_c;

TLorentzVector pi0lv;
pi0lv = *pho1 + *pho2;

//  updating forground and background
if ((pi0lv.M() > 0.02) && (mesonlv_1.M() > 0.02) && (mesonlv_2.M() > 0.02))
{
// fill the tower by tower histograms with invariant mass
pairInvMassTotalBkgd->Fill(mesonlv_1.M());
pairInvMassTotalBkgd->Fill(mesonlv_2.M());
pairInvMassTotal->Fill(pi0lv.M());
}
				
} // inner loop for clusters
} // outer loop for clusters	
for (int j = 0; j < nFourVector; j++) {delete saveFV[j];}
//delete [] saveFV;
} // reading every event
//delete [] saveFV;
}
*/



/*
//______________________________________________________________________________..
void TruthParticleAna::Loop_background_position_swapppingF(int nevts, TString _filename, TTree * intree)
{
// set of parameters used for different selection criteria (cuts)
int ncluster_cut = 6000; // cut to limit the maximum number of clusters
float first_pt_cut = 3.0; // initialize first photon pt cut
float second_pt_cut = 2.0; // initialize second photon pt cut
float alpha_cut = 0.5; // energy asymetry cut
float delR_cut = 0.8; // cone cut

std::cout << "starting loop to construct foreground and background using position swapping method" << std::endl;

TTree * t1 = intree;
if (!intree)
{
TFile *f = new TFile(_filename);
t1 = (TTree *) f->Get("_eventTree");
}
  
// Set Branches
t1->SetBranchAddress("_eventNumber", &_eventNumber);
t1->SetBranchAddress("_nFourVector", &_nFourVector);
t1->SetBranchAddress("_fv_energy", _fv_energy);
t1->SetBranchAddress("_fv_px", _fv_px);
t1->SetBranchAddress("_fv_py", _fv_py);
t1->SetBranchAddress("_fv_pz", _fv_pz);
t1->SetBranchAddress("_fv_Eta", _fv_Eta);
//t1->SetBranchAddress("_pid_primary", _pid_primary);
t1->SetBranchAddress("_pid_secondary", _pid_secondary);


//  int nEntries = (int) t1->GetEntriesFast();
int nEntries = (int) t1->GetEntries();
int nevts2 = nevts;

if (nevts < 0 || nEntries < nevts)
nevts2 = nEntries;

// pre-loop to save all the clusters LorentzVector
//TLorentzVector *saveFV[10000];
//std::vector<std::unique_ptr<TLorentzVector>> saveFV;

for (int i = 0; i < nevts2; i++)
{
		
TLorentzVector *saveFV[10000];
// load the ith instance of the TTree
t1->GetEntry(i);

if (i % 1000 == 0) {std::cout << "event number = " << i << std::endl;}

int nFourVector = _nFourVector;
//std::cout << "four vector number = " << nFourVector << std::endl;

if (nFourVector < 10) continue; // if the count of 4-vector is not 2 (min) then we can not make invariant mass plot.

if (nFourVector > ncluster_cut) continue;

// save parent ID in a separate list
//std::vector<int> primary_pid(10000);
std::vector<int> secondary_pid(10000);

//primary_pid.clear();
secondary_pid.clear();

for (int j = 0; j < nFourVector; j++)
{
float px, py, pz, En;
En = _fv_energy[j];
px = _fv_px[j];
py = _fv_py[j];
pz = _fv_pz[j];

saveFV[j] = new TLorentzVector();
saveFV[j]->SetPxPyPzE(px, py, pz, En);
//saveFV.push_back(std::make_unique<TLorentzVector>(px, py, pz, En));

//primary_pid[j] = _pid_primary[j]; // save parent ID in separate vector
secondary_pid[j] = _pid_secondary[j]; // save p-ID in separate vector
}

//std::cout << p_id.size() << std::endl;

TLorentzVector *pho1, *pho2;

for (int jCs = 0; jCs < nFourVector; jCs++) // outer loop
{

pho1 = saveFV[jCs];
			
//if ((primary_pid[jCs] != 111) && (primary_pid[jCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
//if (primary_pid[jCs] != 221) continue; // ignore if primary particle is not eta-meosn
if (secondary_pid[jCs] != 22) continue; // test of which (either of photon/lepton) we are using

//if (fabs(pho1->Pt()) < (first_pt_cut + (nFourVector - 600)/2000))	continue;
if (fabs(pho1->Pt()) < first_pt_cut)	continue;
			
// another loop to go into the saved cluster
// we have removed the situation for double counting
TLorentzVector dummy_LV;
dummy_LV.SetPtEtaPhiE(0, 0, 0, 0);
				
for (int kCs = jCs+1; kCs < nFourVector; kCs++)
{
//if (jCs == kCs) continue;

pho2 = saveFV[kCs];

//if ((primary_pid[kCs] != 111) && (primary_pid[kCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
//if (primary_pid[kCs] != 221) continue; // ignore if primary particle not eta-meosn
if (secondary_pid[kCs] != 22) continue; // test of which (either of photon/lepton) we are using

//if (fabs(pho2->Pt()) < (second_pt_cut + (nFourVector - 600)/2000)) continue;
if (fabs(pho2->Pt()) < second_pt_cut ) continue;
			
alphaCut = fabs((pho1->E() - pho2->E())/(pho1->E()+ pho2->E()));
if (alphaCut > alpha_cut) continue;

if (pho1->DeltaR(*pho2) > delR_cut) continue;
				
TLorentzVector pi0lv;
pi0lv = *pho1 + *pho2;

//  updating forground and background
if ((pho1->E()  > 1.0) && (pho2->E() > 0.6) && (fabs(pi0lv.Pt()) > 1.0))
//if ((fabs(pi0lv.Pt()) > 1.0) && (fabs(pi0lv.Pt()) < 4.0))
{
// fill the tower by tower histograms with invariant mass
//pairInvMassTotalBkgd->Fill(mesonlv_c.M());
pairInvMassTotal->Fill(pi0lv.M());
}				
} // inner loop for clusters
} // outer loop for clusters	
for (int j = 0; j < nFourVector; j++) {delete saveFV[j];}
//delete [] saveFV;
} // reading every event
//delete [] saveFV;
}
*/


/*
//______________________________________________________________________________..
void TruthParticleAna::Loop_background_position_swapppingB(int nevts, TString _filename, TTree * intree)
{
// set of parameters used for different selection criteria (cuts)
int ncluster_cut = 6000; // cut to limit the maximum number of clusters
float first_pt_cut = 3.0; // initialize first photon pt cut
float second_pt_cut = 2.0; // initialize second photon pt cut
float alpha_cut = 0.5; // energy asymetry cut
float delR_cut = 0.8; // cone cut

std::cout << "starting loop to construct foreground and background using position swapping method" << std::endl;

TTree * t1 = intree;
if (!intree)
{
TFile *f = new TFile(_filename);
t1 = (TTree *) f->Get("_eventTree");
}
  
// Set Branches
t1->SetBranchAddress("_eventNumber", &_eventNumber);
t1->SetBranchAddress("_nFourVector", &_nFourVector);
t1->SetBranchAddress("_fv_energy", _fv_energy);
t1->SetBranchAddress("_fv_px", _fv_px);
t1->SetBranchAddress("_fv_py", _fv_py);
t1->SetBranchAddress("_fv_pz", _fv_pz);
t1->SetBranchAddress("_fv_Eta", _fv_Eta);
//t1->SetBranchAddress("_pid_primary", _pid_primary);
t1->SetBranchAddress("_pid_secondary", _pid_secondary);


//  int nEntries = (int) t1->GetEntriesFast();
int nEntries = (int) t1->GetEntries();
int nevts2 = nevts;

if (nevts < 0 || nEntries < nevts)
nevts2 = nEntries;

// pre-loop to save all the clusters LorentzVector
//TLorentzVector *saveFV[10000];
//std::vector<std::unique_ptr<TLorentzVector>> saveFV;

for (int i = 0; i < nevts2; i++)
{
		
TLorentzVector *saveFV[10000];
// load the ith instance of the TTree
t1->GetEntry(i);

if (i % 1000 == 0) {std::cout << "event number = " << i << std::endl;}

int nFourVector = _nFourVector;
//std::cout << "four vector number = " << nFourVector << std::endl;

if (nFourVector < 10) continue; // if the count of 4-vector is not 2 (min) then we can not make invariant mass plot.

if (nFourVector > ncluster_cut) continue;

// save parent ID in a separate list
//std::vector<int> primary_pid(10000);
std::vector<int> secondary_pid(10000);

//primary_pid.clear();
secondary_pid.clear();

for (int j = 0; j < nFourVector; j++)
{
float px, py, pz, En;
En = _fv_energy[j];
px = _fv_px[j];
py = _fv_py[j];
pz = _fv_pz[j];

saveFV[j] = new TLorentzVector();
saveFV[j]->SetPxPyPzE(px, py, pz, En);
//saveFV.push_back(std::make_unique<TLorentzVector>(px, py, pz, En));

//primary_pid[j] = _pid_primary[j]; // save parent ID in separate vector
secondary_pid[j] = _pid_secondary[j]; // save p-ID in separate vector
}

//std::cout << p_id.size() << std::endl;

TLorentzVector *pho1, *pho2;

for (int jCs = 0; jCs < nFourVector; jCs++) // outer loop
{

pho1 = saveFV[jCs];
			
//if ((primary_pid[jCs] != 111) && (primary_pid[jCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
//if (primary_pid[jCs] != 221) continue; // ignore if primary particle is not eta-meosn
//if (secondary_pid[jCs] != 22) continue; // test of which (either of photon/lepton) we are using

//if (fabs(pho1->Pt()) < (first_pt_cut + (nFourVector - 600)/2000))	continue;
if (fabs(pho1->Pt()) < first_pt_cut)	continue;
			
// another loop to go into the saved cluster
// we have removed the situation for double counting
TLorentzVector dummy_LV;
dummy_LV.SetPtEtaPhiE(0, 0, 0, 0);
				
for (int kCs = jCs+1; kCs < nFourVector; kCs++)
{
//if (jCs == kCs) continue;

pho2 = saveFV[kCs];

//if ((primary_pid[kCs] != 111) && (primary_pid[kCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
//if (primary_pid[kCs] != 221) continue; // ignore if primary particle not eta-meosn
//if (secondary_pid[kCs] != 22) continue; // test of which (either of photon/lepton) we are using

// we will first randomly select the third cluster and extract the "z co-ordinate" from third cluster and swap with second cluster
std::vector<int> third_clus = {jCs, kCs};

int max_loop = 5;

if (max_loop > int (nFourVector/2)) {max_loop = int (nFourVector/2);}
				
for (int nloop = 0; nloop < max_loop; nloop++) // starting the loop that contributes to background
{
int lCs = kCs;
do{
TRandom3 randgen(0);
lCs = randgen.Integer(nFourVector);
} while (std::find(third_clus.begin(), third_clus.end(), lCs) != third_clus.end());

third_clus.push_back(lCs); //push back the randomly generated cluster number

//if (secondary_pid[lCs] != 22) continue; // test of which (either of photon/lepton) we are using
//if (primary_pid[lCs] != 221) continue; // ignore if primary particle not eta-meosn
						
TLorentzVector &pho3 = *(saveFV[lCs]);				

TLorentzVector &pho2_c = dummy_LV;
//pho2_c.SetXYZM(pho2->X(), pho2->Y(), pho3.Z(), pho3.M()); // this is the updated second cluster
pho2_c.SetXYZM(pho3.X(), pho3.Y(), pho2->Z(), pho2->M()); // this is the updated second cluster
				
					
if (fabs(pho2_c.Pt()) < second_pt_cut) continue; //second P_T cut
//if ((fabs(pho2_c.Pt()) < second_pt_cut) || (fabs(pho2_c.Pt()) > (2 * second_pt_cut)))  continue; //second P_T cut
//if (fabs(pho2_c.Pt()) < first_pt_cut) continue; //first P_T cut
	
float alphaCut_c = fabs((pho1->E() - pho2_c.E())/(pho1->E()+ pho2_c.E()));
if (alphaCut_c > alpha_cut) continue;

//float alpha_pt_b = fabs((pho1->Pt() - pho2_c.Pt())/(pho1->Pt()+ pho2_c.Pt()));
//if (alpha_pt_b > 0.95) continue;
//if ((pho2_c.Pt() > (1.5 * second_pt_cut)) && (alpha_pt_b < 0.5)) continue;
//if ((pho2_c.Pt() > (1.5 * second_pt_cut)) || (pho1->Pt() > (1.5 * first_pt_cut))) continue;

if (pho1->DeltaR(pho2_c) > delR_cut) continue;					

//if (fabs(pho2.Z() - pho2_c.Z()) < 2) continue;
				
// getting the 4-vector for eta-meson
TLorentzVector mesonlv_c= *pho1 + pho2_c;

//TLorentzVector pi0lv;
//pi0lv = *pho1 + *pho2;

//  updating forground and background
if ((pho1->E()  > 1.0) && (pho2_c.E() > 0.6) && (fabs(mesonlv_c.M()) > 0.02))
//if ((fabs(mesonlv_c.Pt()) > 1.0) && (fabs(mesonlv_c.Pt()) < 4.0))
{
// fill the tower by tower histograms with invariant mass
pairInvMassTotalBkgd->Fill(mesonlv_c.M());
//pairInvMassTotal->Fill(pi0lv.M());
}				
}// loop for multiple background
third_clus.clear();
} // inner loop for clusters
} // outer loop for clusters	
for (int j = 0; j < nFourVector; j++) {delete saveFV[j];}
} // reading every event
}
*/


/*
//______________________________________________________________________________..
void TruthParticleAna::Loop_background_90deg_rotation(int nevts, TString _filename, TTree * intree)
{
// set of parameters used for different selection criteria (cuts)
int ncluster_cut = 600; // cut to limit the maximum number of clusters
float first_pt_cut = 1.2; // initialize first photon pt cut
float second_pt_cut = 0.8; // initialize second photon pt cut
float alpha_cut = 0.7; // energy asymetry cut
float delR_cut = 0.8; // cone cut

std::cout << "starting loop to construct foreground and background using 90 degree rotation method" << std::endl;

TTree * t1 = intree;
if (!intree)
{
TFile *f = new TFile(_filename);
t1 = (TTree *) f->Get("_eventTree");
}
  
// Set Branches
t1->SetBranchAddress("_eventNumber", &_eventNumber);
t1->SetBranchAddress("_nFourVector", &_nFourVector);
t1->SetBranchAddress("_fv_energy", _fv_energy);
t1->SetBranchAddress("_fv_px", _fv_px);
t1->SetBranchAddress("_fv_py", _fv_py);
t1->SetBranchAddress("_fv_pz", _fv_pz);
t1->SetBranchAddress("_fv_Eta", _fv_Eta);
//t1->SetBranchAddress("_pid_primary", _pid_primary);
t1->SetBranchAddress("_pid_secondary", _pid_secondary);


//  int nEntries = (int) t1->GetEntriesFast();
int nEntries = (int) t1->GetEntries();
int nevts2 = nevts;

if (nevts < 0 || nEntries < nevts)
nevts2 = nEntries;

// pre-loop to save all the clusters LorentzVector
//TLorentzVector *saveFV[10000];
//std::vector<std::unique_ptr<TLorentzVector>> saveFV;

for (int i = 0; i < nevts2; i++)
{
		
TLorentzVector *saveFV[10000];
// load the ith instance of the TTree
t1->GetEntry(i);

if (i % 1000 == 0) {std::cout << "event number = " << i << std::endl;}

int nFourVector = _nFourVector;
//std::cout << "four vector number = " << nFourVector << std::endl;

if (nFourVector < 10) continue; // if the count of 4-vector is not 2 (min) then we can not make invariant mass plot.

if (nFourVector > ncluster_cut) continue;

// save parent ID in a separate list
//std::vector<int> primary_pid(10000);
std::vector<int> secondary_pid(10000);

//primary_pid.clear();
secondary_pid.clear();

for (int j = 0; j < nFourVector; j++)
{
float px, py, pz, En;
En = _fv_energy[j];
px = _fv_px[j];
py = _fv_py[j];
pz = _fv_pz[j];

saveFV[j] = new TLorentzVector();
saveFV[j]->SetPxPyPzE(px, py, pz, En);
//saveFV.push_back(std::make_unique<TLorentzVector>(px, py, pz, En));

//primary_pid[j] = _pid_primary[j]; // save parent ID in separate vector
secondary_pid[j] = _pid_secondary[j]; // save p-ID in separate vector
}

//std::cout << p_id.size() << std::endl;

TLorentzVector *pho1, *pho2;

for (int jCs = 0; jCs < nFourVector; jCs++) // outer loop
{

pho1 = saveFV[jCs];
			
//if ((primary_pid[jCs] != 111) && (primary_pid[jCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
//if (primary_pid[jCs] != 221) continue; // ignore if primary particle is not eta-meosn
//if (secondary_pid[jCs] != 22) continue; // test of which (either of photon/lepton) we are using

//if (fabs(pho1->Pt()) < (first_pt_cut + (nFourVector - 600)/2000))	continue;
if (fabs(pho1->Pt()) < first_pt_cut)	continue;
			
// another loop to go into the saved cluster
// we have removed the situation for double counting
//TLorentzVector dummy_LV1, dummy_LV2;
//dummy_LV1.SetPtEtaPhiE(0, 0, 0, 0);
//dummy_LV2.SetPtEtaPhiE(0, 0, 0, 0);
				
for (int kCs = jCs+1; kCs < nFourVector; kCs++)
{
//if (jCs == kCs) continue;

pho2 = saveFV[kCs];

//if ((primary_pid[kCs] != 111) && (primary_pid[kCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
//if (primary_pid[kCs] != 221) continue; // ignore if primary particle not eta-meosn
//if (secondary_pid[kCs] != 22) continue; // test of which (either of photon/lepton) we are using

//TLorentzVector &pho2_1 = dummy_LV1;
//TLorentzVector &pho1_1 = dummy_LV2;

//float new_phi1 = pho1->Phi()+1.571;
//float new_phi2 = pho2->Phi()+1.571;

//if (new_phi1 > 3.14159) {new_phi1 = new_phi1 - 6.283;}
//if (new_phi2 < 3.14159) {new_phi2 = new_phi2 - 6.283;}

pho1_1.SetPtEtaPhiE(pho1->Pt(), pho1->Eta(), new_phi1, pho1->E());  // this is for background
//pho2_1.SetPtEtaPhiE(pho2->Pt(), pho2->Eta(), new_phi2, pho2->E());  // this is for background
				
//if (fabs(pho2->Pt()) < (second_pt_cut + (nFourVector - 600)/2000)) continue;
if (fabs(pho2->Pt()) < second_pt_cut ) continue;
//if (fabs(pho2_1.Pt()) < second_pt_cut ) continue;
if (fabs(pho1_1.Pt()) < first_pt_cut ) continue;
				
alphaCut = fabs((pho1->E() - pho2->E())/(pho1->E()+ pho2->E()));
if (alphaCut > alpha_cut) continue;


//float alphaCut_b = fabs((pho1_1.E() - pho2_1.E())/(pho1_1.E()+ pho2_1.E()));
//if (alphaCut_b > alpha_cut) continue;

if (pho1->DeltaR(*pho2) > delR_cut) continue;
//if (pho1_1.DeltaR(pho2_1) > delR_cut) continue;

// getting the 4-vector for eta-meson

TLorentzVector pi0lv;
pi0lv = *pho1 + *pho2;

pho1.RotateZ(TMath::Pi() / 2);
pho2.RotateZ(TMath::Pi() / 2);

TLorentzVector pi0lv_b;
pi0lv_b = *pho1 + *pho2;


//  updating forground and background
if (pho1->E()  > 1.0 && pho2->E() > 0.6 && fabs(pi0lv.Pt()) > 1.0)
{
// fill the tower by tower histograms with invariant mass
pairInvMassTotal->Fill(pi0lv.M());

pairInvMassTotalBkgd->Fill(pi0lv_b.M());
}
} // inner loop for clusters
} // outer loop for clusters	
for (int j = 0; j < nFourVector; j++) {delete saveFV[j];}
	 
} // reading every event

}
*/



/*
//______________________________________________________________________________..
void TruthParticleAna::Loop_background_position_swapppingB(int nevts, TString _filename, TTree * intree)
{
// set of parameters used for different selection criteria (cuts)
int ncluster_cut = 6000; // cut to limit the maximum number of clusters
float first_pt_cut = 2.0; // initialize first photon pt cut
float second_pt_cut = 1.0; // initialize second photon pt cut
float alpha_cut = 0.3; // energy asymetry cut
float delR_cut = 0.5; // cone cut

std::cout << "starting loop to construct foreground and background using position swapping method" << std::endl;

TTree * t1 = intree;
if (!intree)
{
TFile *f = new TFile(_filename);
t1 = (TTree *) f->Get("_eventTree");
}
  
// Set Branches
t1->SetBranchAddress("_eventNumber", &_eventNumber);
t1->SetBranchAddress("_nFourVector", &_nFourVector);
t1->SetBranchAddress("_fv_energy", _fv_energy);
t1->SetBranchAddress("_fv_px", _fv_px);
t1->SetBranchAddress("_fv_py", _fv_py);
t1->SetBranchAddress("_fv_pz", _fv_pz);
t1->SetBranchAddress("_fv_Eta", _fv_Eta);
//t1->SetBranchAddress("_pid_primary", _pid_primary);
t1->SetBranchAddress("_pid_secondary", _pid_secondary);


//  int nEntries = (int) t1->GetEntriesFast();
int nEntries = (int) t1->GetEntries();
int nevts2 = nevts;

if (nevts < 0 || nEntries < nevts)
nevts2 = nEntries;

// pre-loop to save all the clusters LorentzVector
//TLorentzVector *saveFV[10000];
//std::vector<std::unique_ptr<TLorentzVector>> saveFV;

for (int i = 0; i < nevts2; i++)
{
		
TLorentzVector *saveFV[10000];
// load the ith instance of the TTree
t1->GetEntry(i);

if (i % 1000 == 0) {std::cout << "event number = " << i << std::endl;}

int nFourVector = _nFourVector;
//std::cout << "four vector number = " << nFourVector << std::endl;

if (nFourVector < 10) continue; // if the count of 4-vector is not 2 (min) then we can not make invariant mass plot.

if (nFourVector > ncluster_cut) continue;

// save parent ID in a separate list
//std::vector<int> primary_pid(10000);
std::vector<int> secondary_pid(10000);

//primary_pid.clear();
secondary_pid.clear();

for (int j = 0; j < nFourVector; j++)
{
float px, py, pz, En;
En = _fv_energy[j];
px = _fv_px[j];
py = _fv_py[j];
pz = _fv_pz[j];

saveFV[j] = new TLorentzVector();
saveFV[j]->SetPxPyPzE(px, py, pz, En);
//saveFV.push_back(std::make_unique<TLorentzVector>(px, py, pz, En));

//primary_pid[j] = _pid_primary[j]; // save parent ID in separate vector
secondary_pid[j] = _pid_secondary[j]; // save p-ID in separate vector
}

//std::cout << p_id.size() << std::endl;

TLorentzVector *pho1, *pho2;

for (int jCs = 0; jCs < nFourVector; jCs++) // outer loop
{

pho1 = saveFV[jCs];
			
//if ((primary_pid[jCs] != 111) && (primary_pid[jCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
//if (primary_pid[jCs] != 221) continue; // ignore if primary particle is not eta-meosn
//if (secondary_pid[jCs] != 22) continue; // test of which (either of photon/lepton) we are using

//if (fabs(pho1->Pt()) < (first_pt_cut + (nFourVector - 600)/2000))	continue;
if (fabs(pho1->Pt()) < first_pt_cut)	continue;
			
// another loop to go into the saved cluster
// we have removed the situation for double counting
TLorentzVector dummy_LV;
dummy_LV.SetPtEtaPhiE(0, 0, 0, 0);
				
for (int kCs = jCs+1; kCs < nFourVector; kCs++)
{
//if (jCs == kCs) continue;

pho2 = saveFV[kCs];

//if ((primary_pid[kCs] != 111) && (primary_pid[kCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
//if (primary_pid[kCs] != 221) continue; // ignore if primary particle not eta-meosn
//if (secondary_pid[kCs] != 22) continue; // test of which (either of photon/lepton) we are using

//if (fabs(pho2->Pt()) < (second_pt_cut + (nFourVector - 600)/2000)) continue;
if (fabs(pho2->Pt()) < second_pt_cut ) continue;
			
alphaCut = fabs((pho1->E() - pho2->E())/(pho1->E()+ pho2->E()));
if (alphaCut > alpha_cut) continue;

//float alpha_pt = fabs((pho1->Pt() - pho2->Pt())/(pho1->Pt()+ pho2->Pt()));
//if ((pho2->Pt() > (1.5 * second_pt_cut)) && (alpha_pt < 0.5)) continue;

//if ((pho2->Pt() > (1.5 * second_pt_cut)) || (pho1->Pt() > (1.5 * first_pt_cut))) continue;

if (pho1->DeltaR(*pho2) > delR_cut) continue;
				
TLorentzVector pi0lv;
pi0lv = *pho1 + *pho2;

//  updating forground and background
//if ((pho1->E()  > 1.0) && (pho2->E() > 0.6) && ((fabs(pi0lv.Pt()) > 3.0) || (fabs(pi0lv.Pt()) < 4.0)))
if ((fabs(pi0lv.Pt()) > 1.0) && (fabs(pi0lv.Pt()) < 4.0))
{
// fill the tower by tower histograms with invariant mass
//pairInvMassTotalBkgd->Fill(mesonlv_c.M());
pairInvMassTotal->Fill(pi0lv.M());
}

// we will first randomly select the third cluster and extract the "z co-ordinate" from third cluster and swap with second cluster
std::vector<int> third_clus = {jCs, kCs};

int max_loop = 50;

if (max_loop > int (nFourVector/2)) {max_loop = int (nFourVector/2);}
				
for (int nloop = 0; nloop < max_loop; nloop++) // starting the loop that contributes to background
{
int lCs = kCs;
do{
TRandom3 randgen(0);
lCs = randgen.Integer(nFourVector);
} while (std::find(third_clus.begin(), third_clus.end(), lCs) != third_clus.end());

third_clus.push_back(lCs); //push back the randomly generated cluster number

//if (secondary_pid[lCs] != 22) continue; // test of which (either of photon/lepton) we are using
//if (primary_pid[lCs] != 221) continue; // ignore if primary particle not eta-meosn
						
TLorentzVector &pho3 = *(saveFV[lCs]);				

TLorentzVector &pho2_c = dummy_LV;
pho2_c.SetXYZM(pho3.X(), pho3.Y(), pho2->Z(), pho2->M()); // this is the updated second cluster
//pho2_c.SetPtEtaPhiE(pho2->Pt(), pho2->Eta(), pho3.Phi(), pho2->E());  // this is the updated second cluster
				
					
if (fabs(pho2_c.Pt()) < second_pt_cut) continue; //second P_T cut
//if ((fabs(pho2_c.Pt()) < second_pt_cut) || (fabs(pho2_c.Pt()) > (2 * second_pt_cut)))  continue; //second P_T cut
//if (fabs(pho2_c.Pt()) < first_pt_cut) continue; //first P_T cut
	
float alphaCut_c = fabs((pho1->E() - pho2_c.E())/(pho1->E()+ pho2_c.E()));
if (alphaCut_c > alpha_cut) continue;

//float alpha_pt_b = fabs((pho1->Pt() - pho2_c.Pt())/(pho1->Pt()+ pho2_c.Pt()));
//if (alpha_pt_b > 0.95) continue;
//if ((pho2_c.Pt() > (1.5 * second_pt_cut)) && (alpha_pt_b < 0.5)) continue;
//if ((pho2_c.Pt() > (1.5 * second_pt_cut)) || (pho1->Pt() > (1.5 * first_pt_cut))) continue;

if (pho1->DeltaR(pho2_c) > delR_cut) continue;					

//if (fabs(pho2.Z() - pho2_c.Z()) < 2) continue;
				
// getting the 4-vector for eta-meson
TLorentzVector mesonlv_c= *pho1 + pho2_c;

//TLorentzVector pi0lv;
//pi0lv = *pho1 + *pho2;

//  updating forground and background
//if ((pho1->E()  > 1.0) && (pho2_c.E() > 0.6) && ((fabs(mesonlv_c.Pt()) > 3.0) || (fabs(mesonlv_c.Pt()) < 4.0)) )
if ((fabs(mesonlv_c.Pt()) > 1.0) && (fabs(mesonlv_c.Pt()) < 4.0))
{
// fill the tower by tower histograms with invariant mass
pairInvMassTotalBkgd->Fill(mesonlv_c.M());
//pairInvMassTotal->Fill(pi0lv.M());
}
//delete &(pho3);
				
}// loop for multiple background
third_clus.clear();
} // inner loop for clusters
} // outer loop for clusters	
for (int j = 0; j < nFourVector; j++) {delete saveFV[j];}
//delete [] saveFV;
} // reading every event
//delete [] saveFV;
}
*/



//______________________________________________________________________________..
void TruthParticleAna::Loop_back_up(int nevts, TString _filename, TTree * intree)
{
  // set of parameters used for different selection criteria (cuts)
  int ncluster_cut = 2000; // cut to limit the maximum number of clusters
  float first_pt_cut = 3.0; // first photon pt cut
  float second_pt_cut = 2.0; // second photon pt cut
  float alpha_cut = 0.5; // energy asymetry cut
  float delR_cut = 0.5; // cone cut

  /*
  // setting the vector for pdg values of the data
  std::vector<int> ref_pdg = {-3334, -3322, -3312, -3222, -3212, -3122, -3112, -2212, -2112, -431, -421, -411, -321, -211, -13, -11, 
  11, 22, 111, 130, 211, 310, 321, 411, 421, 431, 2112, 2212, 3112, 3122, 3212, 3222, 3312, 3322, 3334};
  */

  std::cout << "starting loop to make four vector and hence reconstructing invariant mass" << std::endl;

  TTree * t1 = intree;
  if (!intree)
    {
      TFile *f = new TFile(_filename);
      t1 = (TTree *) f->Get("_eventTree");
    }
  
  // Set Branches
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nFourVector", &_nFourVector);
  t1->SetBranchAddress("_fv_energy", _fv_energy);
  t1->SetBranchAddress("_fv_px", _fv_px);
  t1->SetBranchAddress("_fv_py", _fv_py);
  t1->SetBranchAddress("_fv_pz", _fv_pz);
  t1->SetBranchAddress("_fv_Eta", _fv_Eta);
  //t1->SetBranchAddress("_primary_id", _primary_id);
  //t1->SetBranchAddress("_parent_id", _parent_id);
  t1->SetBranchAddress("_pid_primary", _pid_primary);
  t1->SetBranchAddress("_pid_secondary", _pid_secondary);
  //t1->SetBranchAddress("_embedding", _embedding);

  // pre-loop to save all the clusters LorentzVector

  //TLorentzVector *saveFV[10000];
  // make TLorentz vector to save 4 vector data
  std::vector<TLorentzVector> saveFV;

  //  int nEntries = (int) t1->GetEntriesFast();
  int nEntries = (int) t1->GetEntries();
  int nevts2 = nevts;

  if (nevts < 0 || nEntries < nevts)
    nevts2 = nEntries;
	
  //int eta_meson = 0; // initiate the count of eta meson
  //int pion_meson = 0; // initiate the count of pi0

  // check put all the pid particles (primary)
  //std::vector<int> list_pri_pid;
  //list_pri_pid.clear();
	
  for (int i = 0; i < nevts2; i++)
    {
      // load the ith instance of the TTree
      t1->GetEntry(i);

      if (i % 10000 == 0) {std::cout << "event number = " << i << std::endl;}

      int nFourVector = _nFourVector;

      if (nFourVector < 2) continue; // if the count of 4-vector is not 2 (min) then we can not make invariant mass plot.

      if (nFourVector > ncluster_cut) continue;
		
      /*
      // save parent ID in a separate list
      std::vector<int> primary_pid(10000);
      std::vector<int> secondary_pid(10000);

      primary_pid.clear();
      secondary_pid.clear();
      */

      saveFV.clear(); // clearing fourvector for new start

      for (int j = 0; j < nFourVector; j++)
	{
		 
	  /*
	    float px, py, pz, En;
	    En = _fv_energy[j];
	    px = _fv_px[j];
	    py = _fv_py[j];
	    pz = _fv_pz[j];
	  */

	  //saveFV[j] = new TLorentzVector();
	  //saveFV[j]->SetPxPyPzE(px, py, pz, En);
		
		 
	  TLorentzVector LV_d;

	  LV_d.SetPxPyPzE(_fv_px[j], _fv_py[j], _fv_pz[j], _fv_energy[j]);
	  saveFV.push_back(LV_d);
		 
	  //primary_pid[j] = _pid_primary[j]; // save parent ID in separate vector
	  //secondary_pid[j] = _pid_secondary[j]; // save p-ID in separate vector
		 

	  //if ((_pid_secondary[j] < 100) && (_pid_secondary[j] > 0)){std::cout << "primary particles = " << _pid_secondary[j] << std::endl;}
		 
	  //if (_pid_primary[j] == 221){eta_meson += 1;}
	  //if (_pid_primary[j] == 111){pion_meson += 1;}

	  //if (std::find(list_pri_pid.begin(), list_pri_pid.end(), _pid_primary[j]) == list_pri_pid.end())
	  //	{list_pri_pid.push_back(_pid_primary[j]);} //push back the new primary particle not previously listed in the vector

	}
    

      //std::cout << p_id.size() << std::endl;
		
		
      for (int jCs = 0; jCs < nFourVector; jCs++) // outer loop
	{
	  TLorentzVector& pho1 = saveFV[jCs];
	  //pho1 = saveFV[jCs];

	  //int pri1 = primary_pid[jCs];
	  //int sec1 = secondary_pid[jCs];
			
	  //if (primary_pid[jCs] != 1) continue; // eta embedding in hijing (jan27 version) has -1 parimary_pid for eta and 1 for hijing
	  //if (secondary_pid[jCs] == 22) continue; // test of which (either of photon/lepton) we are using

	  if (fabs(pho1.Pt()) < first_pt_cut)	continue;
			
	  // another loop to go into the saved cluster
	  // we have removed the situation for double counting
	  for (int kCs = jCs+1; kCs < nFourVector; kCs++)
	    {
	      //if (jCs == kCs) continue;

	      //pho2 = saveFV[kCs];

	      TLorentzVector& pho2 = saveFV[kCs];

	      //int pri2 = primary_pid[kCs];
	      //int sec2 = secondary_pid[kCs];

	      //if ((primary_pid[kCs] != 111) && (primary_pid[kCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
	      //if (primary_pid[kCs] != 1) continue;  // eta embedding in hijing (jan27 version) has -1 parimary_pid for eta and 1 for hijing
	      //if (secondary_pid[kCs] == 22) continue; // test of which (either of photon/lepton) we are using

	      if (fabs(pho2.Pt()) < second_pt_cut) continue;

	      if (fabs(fabs(pho1.Pt()) - fabs(pho2.Pt())) > 1.0) continue;
				
	      alphaCut = fabs((pho1.E() - pho2.E())/(pho1.E()+ pho2.E()));
	      if (alphaCut > alpha_cut) continue;

	      if (pho1.DeltaR(pho2) > delR_cut) continue;

	      TLorentzVector pi0lv;
	      pi0lv = pho1 + pho2;
	      if ((pho1.E()  > 1.0) && (pho2.E() > 0.6) && (fabs(pi0lv.Pt()) > 1.0) && (fabs(pi0lv.M()) > 0.02))
		{
		  //if ((abs(pri1) > 3500) || (abs(pri2) > 3500)) {std::cout << "primary particle = " << pri1 << " , " << pri2 << " , " << pi0lv.M() << std::endl;}
		  // fill the tower by tower histograms with invariant mass
		  //if ((sec1 == -11) || (sec2 == -11)) {std::cout << "pi0 mass = " << pi0lv.M() << std::endl;}
					
		  pairInvMassTotal->Fill(pi0lv.M());

		  /*
		    if ((pi0lv.M() > 0.2) && (pi0lv.M() < 0.4))
		    {
		    auto find_pri1 = std::find(ref_pdg.begin(), ref_pdg.end(), pri1);
		    auto x = std::distance(ref_pdg.begin(), find_pri1); 

		    auto find_pri2 = std::find(ref_pdg.begin(), ref_pdg.end(), pri2);
		    auto y = std::distance(ref_pdg.begin(), find_pri2); 

		    check_PriP1_PriP2->Fill(x, y);
		    }
		  */
					 
		  //pairInvMass_SecP1_SecP2->Fill(pi0lv.M(), sec1, sec2);
		  //pairInvMass_PriP1_PriP2->Fill(pi0lv.M(), pri1, pri2);
		}
	    } // inner loop for clusters
	} // outer loop for clusters

    } // reading every event

  // now reading all the contents of the vector (containing all the particles (primary) species)

  /*
    std::sort(list_pri_pid.begin(), list_pri_pid.end());

    std::cout << "The total number of different species of particles are = " << list_pri_pid.size() << std::endl;

    for (unsigned int kk = 0; kk < list_pri_pid.size(); kk++)
    {std::cout << "Primary particles are = " << list_pri_pid[kk] << std::endl;}
  */


  /*
    std::cout << "The number of eta meson is = " << eta_meson << " . " << std::endl;
    std::cout << "The number of pion is = " << pion_meson << " . " << std::endl;
  */

}

//______________________________________________________________________________..
void TruthParticleAna::Loop_ETAtoPIONS(int nevts, TString _filename, TTree * intree)
{
  // set of parameters used for different selection criteria (cuts)
  int ncluster_cut = 200; // cut to limit the maximum number of clusters
  float first_pt_cut = 0.2; // first photon pt cut
  float second_pt_cut = 0.1; // second photon pt cut, second pT cut is always smaller than first pT cut
  float alpha_cut = 0.8; // energy asymetry cut
  float delR_cut = 2.0; // cone cut

  std::cout << "starting loop to make four vector, identify pi0 and eventually get eta meson from 3 pi0s" << std::endl;

  TTree * t1 = intree;
  if (!intree)
    {
      TFile *f = new TFile(_filename);
      t1 = (TTree *) f->Get("_eventTree");
    }
  
  // Set Branches
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nFourVector", &_nFourVector);
  t1->SetBranchAddress("_fv_energy", _fv_energy);
  t1->SetBranchAddress("_fv_px", _fv_px);
  t1->SetBranchAddress("_fv_py", _fv_py);
  t1->SetBranchAddress("_fv_pz", _fv_pz);
  t1->SetBranchAddress("_fv_Eta", _fv_Eta);
  //t1->SetBranchAddress("_primary_id", _primary_id);
  //t1->SetBranchAddress("_parent_id", _parent_id);
  t1->SetBranchAddress("_pid_primary", _pid_primary);
  t1->SetBranchAddress("_pid_secondary", _pid_secondary);
  t1->SetBranchAddress("_embedding", _embedding);

  //TLorentzVector *saveFV[10000];
  // make TLorentz vector to save 4 vector data
  std::vector<TLorentzVector> saveFV;

  // make TLorentz vector to save pion candidate (possible eta meson decay products)
  std::vector<TLorentzVector> pion;

  //  int nEntries = (int) t1->GetEntriesFast();
  int nEntries = (int) t1->GetEntries();
  int nevts2 = nevts;

  if (nevts < 0 || nEntries < nevts)
    nevts2 = nEntries;
		
  for (int i = 0; i < nevts2; i++)
    {
      // load the ith instance of the TTree
      t1->GetEntry(i);

      if (i % 100 == 0) {std::cout << "event number = " << i << std::endl;}

      int nFourVector = _nFourVector;

      if (nFourVector < 5) continue; 

      if (nFourVector > ncluster_cut) continue;
		
      saveFV.clear(); // clearing fourvector for new start
		
      // clean up vector that are designed to keep eta meson possible decay products
      pion.clear();

      for (int j = 0; j < nFourVector; j++)
	{
	  // more efficient method is not to save clusters that are not required (something we will cut off later)	
	  if (_embedding[j] != 0) continue; // "0" means only hijing and "1" means only embedded particles
	  //if (_pid_secondary[jCs] == 22) continue;
		 
	  TLorentzVector LV_d;
	  LV_d.SetPxPyPzE(_fv_px[j], _fv_py[j], _fv_pz[j], _fv_energy[j]);

	  if (LV_d.Pt() < second_pt_cut) continue; // we do not save fourvector if their pT is less than lower bound in pT
		 
	  saveFV.push_back(LV_d);
	}

      for (size_t jCs = 0; jCs < saveFV.size(); jCs++) // outer loop
	{
	  TLorentzVector& pho1 = saveFV[jCs];
      
	  if (fabs(pho1.Pt()) < first_pt_cut)	continue; // this must exist here as it is higher than second pT cut
			
	  // another loop to go into the saved cluster
	  for (size_t kCs = jCs+1; kCs < saveFV.size(); kCs++)
	    {
	      TLorentzVector& pho2 = saveFV[kCs];

	      if (fabs(pho2.Pt()) < second_pt_cut) continue;

	      alphaCut = fabs((pho1.E() - pho2.E())/(pho1.E()+ pho2.E()));
	      if (alphaCut > alpha_cut) continue;

	      if (pho1.DeltaR(pho2) > delR_cut) continue;

	      TLorentzVector pi0lv;
	      pi0lv = pho1 + pho2;

	      // check if "pi0lv" falls under category of pi0 mass range
	      // if yes then append pi0lv into new vector collection of pi0
	      if ((pi0lv.M() > 0.134) || (pi0lv.M() < 0.136)) {pion.push_back(pi0lv);}
	    } // innner loop for cluster
	} // outer loop for cluster
	
      if (pion.size() < 5) continue;
	
      // once we save all possible pi0 candidate we will loop over them and combine in the group of 3 to get possible eta meson candidate
      for (size_t C1 = 0; C1 < pion.size(); C1++) {
		
	TLorentzVector& pi0_1 = pion[C1];
	//if (pi0_1.Pt() < 3.0) continue;
		
	for (size_t C2 = C1+1; C2 < pion.size(); C2++) {

	  TLorentzVector& pi0_2 = pion[C2];
	  //if ((pi0_2.Pt() > 2.5) || (pi0_2.Pt() < 1.0))  continue;
			
	  for (size_t C3 = C2+1; C3 < pion.size(); C3++) {
				
	    TLorentzVector& pi0_3 = pion[C3];
				
	    //if (pi0_3.Pt() > 1.0) continue;

	    TLorentzVector etalv;
	    etalv = pi0_1 + pi0_2 + pi0_3;
				
	    pairInvMassTotal->Fill(etalv.M());
	  }
	}
      }
    } // reading every event
}

