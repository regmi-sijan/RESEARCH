#include "pi0ClusterAna.h"

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

//____________________________________________________________________________..
pi0ClusterAna::pi0ClusterAna(const std::string &name, const std::string &outName = "pi0ClusterAnaOut"):
SubsysReco(name)
  , clusters_Towers(nullptr)
  , truth_photon(nullptr)
  , truth_pi0(nullptr)
//  , caloevalstack(NULL)
  , m_eta_center()
  , m_phi_center()
  , m_tower_energy()
	, m_cluster_eta()
  , m_cluster_phi()
	, m_cluster_e()
  , m_cluster_chi2()
  , m_cluster_prob()
	, m_cluster_nTowers()
	, alphaCut(-1.0)
//  , m_asym()
//  , m_deltaR()
//  , m_lead_E()
//  , m_sublead_E()
//  , m_lead_phi()
//  , m_lead_eta()
//  , m_sublead_phi()
//  , m_sublead_eta()
//  , m_pi0_E()
//  , m_pi0_eta()
//  , m_pi0_phi()
, n_event(0)
, f_temp(0)
, Outfile(outName)

{
  std::cout << "pi0ClusterAna::pi0ClusterAna(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
pi0ClusterAna::~pi0ClusterAna()
{
  std::cout << "pi0ClusterAna::~pi0ClusterAna() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int pi0ClusterAna::Init(PHCompositeNode *topNode)
{
	n_event = 0; // initialization of event number

  
  std::cout << "pi0ClusterAna::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int pi0ClusterAna::InitRun(PHCompositeNode *topNode)
{
  std::cout << "pi0ClusterAna::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
	
  out = new TFile(Outfile.c_str(),"RECREATE");
	
	/*
	Pi0_pt_eta_mass = new TH3F("Pi0_pt_eta_mass", "3D_Pi0_properties", 300, 0.0, 30, 240, -1.2, 1.2, 400, 0.0, 40);
	Eta_pt_eta_mass = new TH3F("Eta_pt_eta_mass", "3D_Eta_properties", 300, 0.0, 30, 240, -1.2, 1.2, 400, 0.0, 40);
	Pi0_pt = new TH1F("Pi0_pt", "Pi0_Pt", 300, 0.0, 30);
	Eta_pt = new TH1F("Eta_pt", "Eta_Pt", 300, 0.0, 30);
	*/

	//pairInvMassTotal = new TH1F("pairInvMassTotal", "invariant mass histogram", 240, 0.0, 1.2);
	
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
    _eventTree->Branch("_pid_primary", _pid_primary, "_pid_primary[_nFourVector]/F"); // this is pid of primary particles
    _eventTree->Branch("_pid_secondary", _pid_secondary, "_pid_secondary[_nFourVector]/F"); // this is pid for secondar particles
    _eventTree->Branch("_embedding", _embedding, "_embedding[_nFourVector]/F"); // indicatior (flag) to differentiate embedded particles with others
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int pi0ClusterAna::process_event(PHCompositeNode *topNode)
{
	_eventNumber = n_event; // updating event number to the Ttree

	if (n_event%10 == 0) {std::cout << " processing event number " << n_event << std::endl;}

  //truth particle information
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if(!truthinfo)
    {
      std::cout << PHWHERE << "pi0ClusterAna::process_event Could not find node G4TruthInfo"  << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
   
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// SCANNING PRIMARY PARTICLE RANGE
	 
  // get primary particle
  PHG4TruthInfoContainer::Range truthRange = truthinfo->GetPrimaryParticleRange();
  PHG4TruthInfoContainer::ConstIterator truthIter;

  PHG4Particle *truth_pri;
	
	std::vector<int> primary_id; // primary id and primary track id should be same but lets keep both
	std::vector<int> primary_track_id; // this will be one of the crucial to track final particles in list of all secondary particles
	std::vector<int> pid;
	std::vector<int> embbed_flag;
	
	primary_id.clear();
	primary_track_id.clear();
	pid.clear();
	embbed_flag.clear();


	int eta_flag = 0;	
  for(truthIter = truthRange.first; truthIter != truthRange.second; truthIter++)
    {
      truth_pri = truthIter->second;
			
			if ((truth_pri->get_parent_id()) != 0) {std::cout << "Secondary Particle Spotted" << std::endl; continue;}
			
			// update pid, primary id and track id for all primary particles (hijing and embedded as well)
			pid.push_back(truth_pri->get_pid()); // update pid of primary particle (eta-meson in this case)
		  primary_id.push_back(truth_pri->get_primary_id()); // append primary id of the primary particle (we will use it to match with secondary particle)
		  primary_track_id.push_back(truth_pri->get_track_id()); // append track id of the primary particle (useful to identify the first match in secondary particle)
	
			// since we do not have eta meson as primary particle (so once eta meson is found it indicates the start of embedding particles)
			if (((truth_pri->get_pid()) == 221)||(eta_flag == 1))
			 {
				eta_flag = 1 ; // once we enterd into this loop we just keep it as this is the start of embeding
				embbed_flag.push_back(1); // indicating embedding
			 }
			else {embbed_flag.push_back(0);} // for hijng events								
	  }
	

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// SCANNING SECONDARY PARTICLE RANGE
	
	// we will first scan over all secondary particles and save few things like pid, parent id and track id of secondary particles
  PHG4TruthInfoContainer::Range truthRange_sec = truthinfo->GetSecondaryParticleRange();
  PHG4TruthInfoContainer::ConstIterator truthIter_sec;

  PHG4Particle *truth_sec;
	
	// vectors to store the information for ntuple
	std::vector<int> _pid_sec;
	std::vector<int> _parentid_sec;
	std::vector<int> _trackid_sec;
	std::vector<int> _primaryid_sec; // primary id for secondary particle must be same as primary id for primary particles (but keep both of them)
	
	
	_pid_sec.clear();
	_parentid_sec.clear();
	_trackid_sec.clear();
	_primaryid_sec.clear();
	
  for(truthIter_sec = truthRange_sec.first; truthIter_sec != truthRange_sec.second; truthIter_sec++)
   {
     truth_sec = truthIter_sec->second;
						
		 if ((truth_sec->get_primary_id()) < 1) {std::cout << "non positive primary id spotted" << std::endl;}

		 // updating everything in vector (this is to prepare for selecting useful secondary particles before final update in TTree Ntuples)
		 _pid_sec.push_back(truth_sec->get_pid()); // for secondary particle
		 _parentid_sec.push_back(truth_sec->get_parent_id()); // update the value of parent id of secondary particles
		 _trackid_sec.push_back(truth_sec->get_track_id()); //update the value of track id of secondary particles
		 _primaryid_sec.push_back(truth_sec->get_primary_id()); //update the value of primary id of secondary particles
		 
	 }
	
	// now we will finalize which one to keep and which one not-to in Ttree (we will stop when we encounter at photon, electron, positron)	 
		
  // using unordered_map to store elements of _parentid_sec and their indices
  std::unordered_map<int, std::vector<size_t>> map_parentid_sec;

  for (size_t i = 0; i < _parentid_sec.size(); ++i) {map_parentid_sec[_parentid_sec[i]].push_back(i);} // updating element in unordered map
		
	// first two vectors are to store track id in temporary manner and the third one is to store final version of track id that we need
	std::vector<int> tmp_trackid_sec;
	std::vector<int> temp_trackid_sec;
	std::vector<int> final_trackid_sec; // this saves track id of all tracks we need

  // Iterate over each element in primary_track_id and check if it exists in map_parentid_sec (secondary particles)
  for (size_t i = 0; i < primary_track_id.size(); ++i) 
	 {
		tmp_trackid_sec.clear();
		tmp_trackid_sec.push_back(primary_track_id[i]); // update "primary track id" we will be looking at for this iteration  

		while (tmp_trackid_sec.size() > 0) // we will use until we are done with finding all the track id particles
		 {
			temp_trackid_sec.clear();

			for (size_t j = 0; j < tmp_trackid_sec.size(); ++j)
			 { 
				if (map_parentid_sec.find(tmp_trackid_sec[j]) != map_parentid_sec.end()) // trying to find a match
				 {
					// If there's a match, we save the trackid
					for (auto index : map_parentid_sec[tmp_trackid_sec[j]]) // index is the location saved in parent_id map which can map to any other vector
					 {
						 
						 	if (_pid_sec[index] == 221) std::cout << "pid sec, track id sec, primaryid sec, parent id sec " << _pid_sec[index] << " , " << _trackid_sec[index] << " , " << _primaryid_sec[index] << " , " << _parentid_sec[index] << std::endl;


						//if (_pid_sec[index] == 221) continue; // immediately skip if there is eta meson
						
						// update as a final match if it is photon, electron or positron
						if ((_pid_sec[index] == 22) || (_pid_sec[index] == 11) || (_pid_sec[index] == -11)) {final_trackid_sec.push_back(_trackid_sec[index]);}
						else {temp_trackid_sec.push_back(_trackid_sec[index]);} // we have to further dig down
					 }
				 }
			 }
			// now updating all temp_trackid_sec to tmp_trackid_sec
			tmp_trackid_sec.clear();
			for (size_t jj = 0; jj < temp_trackid_sec.size(); ++jj) {tmp_trackid_sec.push_back(temp_trackid_sec[jj]);}
		 }
   }

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// EXTRACTING THE INFORMATION FROM TRACK_ID (SECONDARY PARTICLE)

	// now with the help of trackid we will extract all the information of secondary particles and save it in TTree histogram	
	
	// vectors to store the information for ntuple
	std::vector<float> _PX;
	std::vector<float> _PY;
	std::vector<float> _PZ;
	std::vector<float> _E;
	std::vector<float> _Eta;
	std::vector<int> _PID_PRI;
	std::vector<int> _PID_SEC;
	std::vector<int> _PRIMARY_ID;
	std::vector<int> _EMBEDDING_ID;

	// preparing all vectors 
	_PX.clear();
	_PY.clear();
	_PZ.clear();
	_E.clear();
	_Eta.clear();
	_PID_PRI.clear();
	_PID_SEC.clear();
	_PRIMARY_ID.clear();
	_EMBEDDING_ID.clear();


	// looping over "final_trackid_sec" trackid and extracting all the essentail information
	for (unsigned int kk = 0; kk < final_trackid_sec.size(); kk++)
	 {
		PHG4Particle *final_par = truthinfo->GetParticle(final_trackid_sec[kk]);
			
		if (final_par == nullptr) {std::cout << " null pointer" << std::endl;}

		if ((final_par->get_e() < 0.4) || (fabs(getEta(final_par)) > 1.1)) continue; // apply some conditions (to help to reduce the size of output data)
		//if (fabs(getEta(final_par)) > 1.1) continue; // apply some conditions (to help to reduce the size of output data)
		
		auto iter_ = std::find(primary_id.begin(), primary_id.end(), final_par->get_primary_id());			
		if (iter_ != primary_id.end()) // if we find a match
		 {
			int pos_ = std::distance(primary_id.begin(), iter_); // this gives the index (position) where we find a match 
	 		
			_PID_PRI.push_back(pid[pos_]); // update idea of pid of primary particle
			_EMBEDDING_ID.push_back(embbed_flag[pos_]); // update if this came from embedding particles or not 

			_PID_SEC.push_back(final_par->get_pid()); // pid for secondary particle
			_PRIMARY_ID.push_back(final_par->get_primary_id()); // update the value of primary id of primary/secondary particles (same primary id means they are connected)
		 	_Eta.push_back(getEta(final_par)); // updating pseudorapidity of secondary particle
			_PX.push_back(final_par->get_px());
			_PY.push_back(final_par->get_py());
		  _PZ.push_back(final_par->get_pz());
			_E.push_back(final_par->get_e());
		 }
		else {std::cout << "There is no match in primary id (while doing GetParticle with track-id." << std::endl;}		
	 }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// SAVING IN TTREE
	
	// updating in TTree
	_eventNumber = n_event;
	_nFourVector = _E.size();

	// if the size of vector is more than 600K (our limitation of the vector) we ignore that event
	if (_E.size() > 50000) {std::cout << " The size of vector is = " << _E.size() << " . " << std::endl;}
	if (_E.size() < 1) {std::cout << "Missing data in event number = " << n_event << std::endl;}
	
	// updating Ttree
	if ((_E.size() > 0) && (_E.size() < 99000))
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
			_primary_id[kk] = _PRIMARY_ID[kk];
			_embedding[kk] = _EMBEDDING_ID[kk];
		}
	}	
		  
  n_event++; // updating event number
	_eventTree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int pi0ClusterAna::ResetEvent(PHCompositeNode *topNode)
{
  //std::cout << "pi0ClusterAna::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;

  m_eta_center.clear();
  m_phi_center.clear();
  m_tower_energy.clear();
  m_cluster_eta.clear();
  m_cluster_phi.clear();
  m_cluster_e.clear();
  m_cluster_chi2.clear();
  m_cluster_prob.clear();
  m_cluster_nTowers.clear();
  m_asym.clear();
  m_deltaR.clear();
  m_lead_E.clear();
  m_sublead_E.clear();
  m_lead_phi.clear();
  m_lead_eta.clear();
  m_sublead_phi.clear();
  m_sublead_eta.clear();
  //m_pi0_E.clear();
  //m_pi0_eta.clear();
  //m_pi0_phi.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int pi0ClusterAna::EndRun(const int runnumber)
{
  std::cout << "pi0ClusterAna::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int pi0ClusterAna::End(PHCompositeNode *topNode)
{
  std::cout << "pi0ClusterAna::End(PHCompositeNode *topNode) This is the End" << std::endl;

	if (topNode == 0 && f_temp)
	{
		out->Close();
		f_temp->Close();
		delete f_temp;
		delete out;
		return Fun4AllReturnCodes::EVENT_OK;
	}
  out->cd();
  
  //truth_pi0 -> Write();
  //truth_photon -> Write();
  //clusters_Towers -> Write();
	
  out->Write();
  out->Close();
  delete out; 

  //hm -> dumpHistos(Outfile.c_str(),"UPDATE");


  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int pi0ClusterAna::Reset(PHCompositeNode *topNode)
{
  std::cout << "pi0ClusterAna::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void pi0ClusterAna::Print(const std::string &what) const
{
  std::cout << "pi0ClusterAna::Print(const std::string &what) const Printing info for " << what << std::endl;
}
//____________________________________________________________________________.. 
float pi0ClusterAna::getEta(PHG4Particle *particle)
{
  float px = particle -> get_px();
  float py = particle -> get_py();
  float pz = particle -> get_pz();
  float p = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));

  return 0.5*log((p+pz)/(p-pz));
}


//______________________________________________________________________________..
void pi0ClusterAna::Loop(int nevts, TString _filename, TTree * intree)
{
	// set of parameters used for different selection criteria (cuts)
	int ncluster_cut = 6000; // cut to limit the maximum number of clusters
	float first_pt_cut = 3.0; // first photon pt cut
	float second_pt_cut = 2.0; // second photon pt cut
	float alpha_cut = 0.5; // energy asymetry cut
	float delR_cut = 0.8; // cone cut

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
  t1->SetBranchAddress("_pid_primary", _pid_primary);
  t1->SetBranchAddress("_pid_secondary", _pid_secondary);

  // pre-loop to save all the clusters LorentzVector

  TLorentzVector *saveFV[10000];

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

		if (nFourVector < 2) continue; // if the count of 4-vector is not 2 (min) then we can not make invariant mass plot.

		if (nFourVector > ncluster_cut) continue;

		// save parent ID in a separate list
		std::vector<int> primary_pid(10000);
		std::vector<int> secondary_pid(10000);

		primary_pid.clear();
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

		 primary_pid[j] = _pid_primary[j]; // save parent ID in separate vector
		 secondary_pid[j] = _pid_secondary[j]; // save p-ID in separate vector
    }

		//std::cout << p_id.size() << std::endl;

    TLorentzVector *pho1, *pho2;

    for (int jCs = 0; jCs < nFourVector; jCs++) // outer loop
    {
      pho1 = saveFV[jCs];
			
			//if ((primary_pid[jCs] != 111) && (primary_pid[jCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
			if (primary_pid[jCs] != 221) continue; // ignore if primary particle is neither of pi0 or eta-meosn
			if (secondary_pid[jCs] != 22) continue; // test of which (either of photon/lepton) we are using

      if (fabs(pho1->Pt()) < first_pt_cut)	continue;
			
      // another loop to go into the saved cluster
			// we have removed the situation for double counting
      for (int kCs = jCs+1; kCs < nFourVector; kCs++)
      {
        //if (jCs == kCs) continue;

        pho2 = saveFV[kCs];

				//if ((primary_pid[kCs] != 111) && (primary_pid[kCs] != 221)) continue; // ignore if primary particle is neither of pi0 or eta-meosn
				if (primary_pid[kCs] != 221) continue; // ignore if primary particle is neither of pi0 or eta-meosn
				if (secondary_pid[kCs] != 22) continue; // test of which (either of photon/lepton) we are using

				if (fabs(pho2->Pt()) < second_pt_cut) continue;
				
				alphaCut = fabs((pho1->E() - pho2->E())/(pho1->E()+ pho2->E()));
				if (alphaCut > alpha_cut) continue;

        if (pho1->DeltaR(*pho2) > delR_cut) continue;

				TLorentzVector pi0lv;
        pi0lv = *pho1 + *pho2;
				if (pho1->E()  > 1.0 && pho2->E() > 0.6 && fabs(pi0lv.Pt()) > 1.0)
	  		{
					// fill the tower by tower histograms with invariant mass
					pairInvMassTotal->Fill(pi0lv.M());
	  		}
       } // inner loop for clusters
      } // outer loop for clusters
  } // reading every event
}


