#include "CaloCalibEmc_eta.h"

// Calobase headers
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>

// Global Vertex headers
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

// Fun4All headers
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

// Phool headers
#include <phool/getClass.h>
#include <phool/phool.h>

// ROOT headers
#include <TF1.h>
#include <TVectorF.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TList.h>
#include <TArray.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TMath.h>
#include "TVector3.h"
#include <deque>
#include <sys/time.h>
#include <TRandom.h>

// CLHEP header
#include <CLHEP/Vector/ThreeVector.h>

// Standard Library headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <string>
#include <random>
#include <fstream>
#include <ctime>

// Centrality headers
#include <centrality/CentralityInfo.h>
#include <centrality/CentralityInfov1.h>
#include <ffarawobjects/Gl1Packet.h>

//using namespace std;

//____________________________________________________________________________..
CaloCalibEmc_eta::CaloCalibEmc_eta(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_ievent(0)
  , m_Filename(filename)
  , cal_output(0)
  , _caloname("CEMC")
    //, fitp1_eta_phi2d(0)
    //, pairInvMassTotal(0)
  ,_eventTree(0)
  ,_eventNumber(-1)
  ,_nClusters(-1)
  , maxTowerEta(-1)
  , maxTowerPhi(-1)
  , alphaCut(-1.0)   
  , f_temp(0)  
{
  // different type of fore-ground we are making
  pairInvMassTotal = NULL;
  pairInvMassTotalBkgd = NULL;

  // different types of background we are making
  pairInvMassPtEta = NULL;
  pairInvMassPtEtaBkgd = NULL;
	
  //m_cent_nclus_cut = 1000;

}

//____________________________________________________________________________..
int CaloCalibEmc_eta::InitRun(PHCompositeNode *topNode)
{
  std::cout << "CaloCalibEmc_eta::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  m_ievent = 0;

  cal_output = new TFile(m_Filename.c_str(), "RECREATE");

  pairInvMassTotal = new TH1F("pairInvMassTotal", "Pair_Mass_Histo", 300, 0.0, 3.0);
  pairInvMassTotal->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  pairInvMassTotal->SetYTitle("Counts");

  pairInvMassTotalBkgd = new TH1F("pairInvMassTotalBkgd", "Pair_Mass_Histo_Bkgd", 300, 0.0, 3.0);

  pairInvMassPtEta = new TH3F("pairInvMassPtEta", "Pair_Mass_Histo_PtEta", 300, 0.0, 3.0, 200, 0.0, 20.0, 50, -2.5, 2.5);
  pairInvMassPtEtaBkgd = new TH3F("pairInvMassPtEtaBkgd", "Pair_Mass_Histo_PtEta_Bkgd", 300, 0.0, 3.0, 200, 0.0, 20.0, 50, -2.5, 2.5);

  pairInvMassPtdelR = new TH3F("pairInvMassPtdelR", "Pair_Mass_Histo_PtdelR", 300, 0.0, 3.0, 200, 0.0, 20.0, 150, 0.0, 1.5);
  pairInvMassPtdelRBkgd = new TH3F("pairInvMassPtdelRBkgd", "Pair_Mass_Histo_PtdelR_Bkgd", 300, 0.0, 3.0, 200, 0.0, 20.0, 150, 0.0, 1.5);

  pairpTDelPhiDelEta = new TH3F("pairpTDelPhiDelEta", "pairpT_delphi_deleta", 200, 0.0, 20.0, 100, -3.2, 3.2, 250, -2.5, 2.5);
  pairpTDelPhiDelEtaBkgd = new TH3F("pairpTDelPhiDelEtaBkgd", "pairpT_delphi_deleta", 200, 0.0, 20.0, 100, -3.2, 3.2, 250, -2.5, 2.5);
		
  DelR_pairpT_f = new TH2F("DelR_pairpT_f", "DelR_pairpT", 300, 0.0, 1.5, 200, 0.0, 20.0);
  DelR_pairpT_b = new TH2F("DelR_pairpT_b", "DelR_pairpT", 300, 0.0, 1.5, 200, 0.0, 20.0);
  
  h_etaphi_clus = new TH2F("h_etaphi_clus", "etaphi_clus", 250, -2.5, 2.5, 64, -1 * M_PI, M_PI);  
  
  if (topNode != 0)
    {
      // TTree declare
      _eventTree = new TTree("_eventTree", "An event level info Tree");
      // TTree branches
      _eventTree->Branch("_eventNumber", &_eventNumber, "_eventNumber/I");
      _eventTree->Branch("_nClusters", &_nClusters, "_nClusters/I");
      _eventTree->Branch("_triggerVector", _triggerVector, "_triggerVector[64]/I");
      //_eventTree->Branch("_nCentrality", &_nCentrality, "_nCentrality/F");
      _eventTree->Branch("_vertex", _vertex, "_vertex[3]/F"); // save co-ordinates of global-vertex (x, y and z co-ordinates)
      //_eventTree->Branch("_nTowers", &_nTowers, "_nTowers[_nClusters]/I");
      //_eventTree->Branch("_clusterIDs", _clusterIDs, "_clusterIDs[_nClusters]/I");
      _eventTree->Branch("_clusterEnergies", _clusterEnergies, "_clusterEnergies[_nClusters]/F");
      _eventTree->Branch("_clusterPts", _clusterPts, "_clusterPts[_nClusters]/F");
      _eventTree->Branch("_clusterEtas", _clusterEtas, "_clusterEtas[_nClusters]/F");
      _eventTree->Branch("_clusterPhis", _clusterPhis, "_clusterPhis[_nClusters]/F");
      _eventTree->Branch("_maxTowerEtas", _maxTowerEtas, "_maxTowerEtas[_nClusters]/I");
      _eventTree->Branch("_maxTowerPhis", _maxTowerPhis, "_maxTowerPhis[_nClusters]/I");
      _eventTree->Branch("_clusterChi2", _clusterChi2, "_clusterChi2[_nClusters]/F");
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
//____________________________________________________________________________..

int CaloCalibEmc_eta::process_event(PHCompositeNode *topNode)
{
  // list of all cuts we are going to use within process event
  float pt1_cut = 0.0;
  float pt2_cut = 0.0;
  int nCluster_cut = 10;
  float alpha_cut = 0.6;
  float delR_cut = 1.1;

  if (m_ievent % 500 == 0)
    {
      std::cout << std::endl;
      std::cout << "Beginning of the event " << m_ievent << std::endl;
      std::cout << "====================================" << std::endl;
    }

  _eventNumber = m_ievent++;

  //--------------------------- trigger and GL1-------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------------------------------//
  /*
  //check if our event is/not min bias
  bool isMinBias = false;
  Gl1Packet* gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  if (!gl1PacketInfo)
  {
  std::cout << PHWHERE << "CaloValid::process_event: GL1Packet node is missing" << std::endl;
  }

  if (gl1PacketInfo)
  {
  uint64_t triggervec = gl1PacketInfo->getScaledVector();
  if ((triggervec >> 10U) & 0x1U  || (triggervec >> 11U) & 0x1U|| (triggervec >> 12U) & 0x1U )
	
  {
  isMinBias = true;
  }

  }

  if (isMinBias != true)
  {
  return Fun4AllReturnCodes::EVENT_OK;
  }

  */	
  // saving all the trigger-bit information in data	
  Gl1Packet *gl1_packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  if (gl1_packet)
    {
      uint64_t gl1_scaledtriggervector = gl1_packet->lValue(0, "ScaledVector");

      for (int i = 0; i < 64; i++)
	{
	  int triggerbit = ((gl1_scaledtriggervector & 0x1U) == 0x1U);
	  gl1_scaledtriggervector = (gl1_scaledtriggervector >> 1U) & 0xffffffffU;

	  // updating in ntuples
	  _triggerVector[i] =  triggerbit;
	}
    }

  //----------------------------------get vertex-----------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------------------------------//
  GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
    {
      // std::cout << PHWHERE << " Fatal Error - GlobalVertexMap node is missing"<< std::endl;
      std::cout << "pi0EtaByEta GlobalVertexMap node is missing" << std::endl;
      // return Fun4AllReturnCodes::ABORTRUN;
    }
  
  // Get Vertex
  float vx = 0;
  float vy = 0;
  float vz = 0;

  bool found_vertex = false;
  if (vertexmap && !vertexmap->empty()) 
    {
      GlobalVertex *vtx = vertexmap->begin()->second;
      if (vtx)
	{
	  if (m_use_vertextype) 
	    {
	      auto typeStartIter = vtx->find_vertexes(m_vertex_type);
	      auto typeEndIter = vtx->end_vertexes();
	      for (auto iter = typeStartIter; iter != typeEndIter; ++iter)
		{
		  const auto &[type, vertexVec] = *iter;
		  if (type != m_vertex_type) { continue; }
		  for (const auto *vertex : vertexVec)
		    {
		      if (!vertex) { continue; }
		      //vx = vertex->get_x();
		      //vy = vertex->get_y();
		      vz = vertex->get_z();
		      found_vertex = true;
		    }
		}
	    } 
	  else 
	    {
	      //vx = vtx->get_x();
	      //vy = vtx->get_y();
	      vz = vtx->get_z();
	      found_vertex = true;
	    }
	}
    }
  
  if (!found_vertex && reqVertex) 
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }

	
  // saving vertex in ntuple
  _vertex[0] = vx;
  _vertex[1] = vy;
  _vertex[2] = vz;

  //-------------------------------------------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------------------------------//

  // create a cluster object
  //std::string clusnodename = "CLUSTERINFO_POS_COR_" + _caloname; // two options (pos-corr node and raw cluster
  std::string clusnodename = "CLUSTERINFO_" + _caloname; // two options (pos-corr node and raw cluster

  RawClusterContainer *recal_clusters = findNode::getClass<RawClusterContainer>(topNode,clusnodename.c_str());

  if (!recal_clusters)
    {
      std::cout << PHWHERE << clusnodename << " node is missing" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
	
  //-------------------------------------------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------------------------------//
  
  // get towers and tower-level information
  TowerInfoContainer *_towerinfos = NULL;
		
  std::string towernode = "TOWERINFO_CALIB_" + _caloname;
  _towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towernode.c_str());
			
  if (!_towerinfos) 
    {
      std::cout << PHWHERE << " ERROR: Can't find " << towernode << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  //-------------------------------------------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------------------------------//
  // create a tower geometry object
  // geometry object is same in tower info and raw tower
  std::string towergeomnode = "TOWERGEOM_" + _caloname;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnode.c_str());

  if (!towergeom)
    {
      std::cout << PHWHERE << ": Could not find node " << towergeomnode << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
	
  //-------------------------------------------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------------------------------//
  /*
  //get centrality node
  CentralityInfo *cent = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!cent)
  {
  std::cout << " Centrality Node is missing." << std::endl;
  }

  //get centrality of event
  float cent_val = cent->get_centile(CentralityInfo::PROP::mbd_NS);

  _nCentrality = cent_val;
  */
  
  //-------------------------------------------------------------------------------------------------------------------//
  //-------------------------------------------------------------------------------------------------------------------//
  
  CLHEP::Hep3Vector vertex(vx, vy, vz); // creating 3-vector

  // loop over the clusters
  RawClusterContainer::ConstRange t_rbegin_end = recal_clusters->getClusters();
  RawClusterContainer::ConstIterator t_rclusiter;

  RawCluster *savCs[10000];  // savingClusters that has 1 GeV or more
  int iCs = 0;

  // saving the clusters
  for (t_rclusiter = t_rbegin_end.first; t_rclusiter != t_rbegin_end.second; ++t_rclusiter)
    {
      RawCluster *t_recalcluster = t_rclusiter->second;

      if ((t_recalcluster->get_ecore() > 0.6) &&  (t_recalcluster->get_chi2() < 6.0))  savCs[iCs++] = t_recalcluster; 
    }

  _nClusters = iCs;

  // looping on the saved clusters savCs
  // outer loop (we want to do pair of the loops to get invariant mass distribution)
  for (int jCs = 0; jCs < _nClusters; jCs++)
    {
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*savCs[jCs], vertex);

      //_clusterIDs[jCs] = savCs[jCs]->get_id();

      _clusterChi2[jCs] = savCs[jCs]->get_chi2();

      //vector to hold all the towers etas, phis, and energy in this cluster
      std::vector<int> toweretas;
      std::vector<int> towerphis;
      std::vector<float> towerenergies;
	
      TowerInfo * towinfo;
		
      // loop over the towers from the outer loop cluster and find the max tower location
      RawCluster::TowerConstRange towers = savCs[jCs]->get_towers();
      RawCluster::TowerConstIterator toweriter;
    
      for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
	{
	  //RawTowerGeom *tower_geom_instance = geom->get_tower_geometry(toweriter->first);
	  //unsigned int towerindex = _towerinfos->decode_key(toweriter->first);
	  int iphi = RawTowerDefs::decode_index2(toweriter->first);  // index2 is phi in CYL
	  int ieta = RawTowerDefs::decode_index1(toweriter->first);  // index1 is eta in CYL	    
	  unsigned int towerkey = iphi + (ieta << 16U);
	  unsigned int towerindex =  _towerinfos->decode_key(towerkey);
	    
	  towinfo = _towerinfos->get_tower_at_channel(towerindex);
	    
	  double towerenergy = towinfo->get_energy();

	  // put the eta, phi, energy into corresponding vectors
	  toweretas.push_back(ieta);
	  towerphis.push_back(iphi);
	  towerenergies.push_back(towerenergy);
	}
    
      int maxTowerIndex = max_element(towerenergies.begin(), towerenergies.end()) - towerenergies.begin();
      maxTowerEta = toweretas[maxTowerIndex];
      maxTowerPhi = towerphis[maxTowerIndex];
		
      _maxTowerEtas[jCs] = maxTowerEta;
      _maxTowerPhis[jCs] = maxTowerPhi;
   
      // extracting parameters and saving to Ttree
      float tt_clus_energy = E_vec_cluster.mag();
      float tt_clus_eta = E_vec_cluster.pseudoRapidity();
      float tt_clus_pt = E_vec_cluster.perp();
      float tt_clus_phi = E_vec_cluster.getPhi();

      _clusterEnergies[jCs] = tt_clus_energy;
      _clusterPts[jCs] = tt_clus_pt;
      _clusterEtas[jCs] = tt_clus_eta;
      _clusterPhis[jCs] = tt_clus_phi;

    
      if ((tt_clus_pt < pt1_cut) || (_nClusters < nCluster_cut)) continue; //first photon cut and centrality cut (based on number of clusters)
				
      // another loop to go into the saved cluster
      for (int kCs = 0; kCs < _nClusters; kCs++)
	{
	  if (jCs == kCs) continue;

	  CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetECoreVec(*savCs[kCs], vertex);
	      
	  float tt2_clus_energy = E_vec_cluster2.mag();
	  float tt2_clus_pt = E_vec_cluster2.perp();
							
	  if (tt2_clus_pt < pt2_cut) continue;  // second photon cut
				
	  // applying alpha cut (energy assymetry cut)
	  alphaCut = abs(tt_clus_energy - tt2_clus_energy) / (tt_clus_energy + tt_clus_energy);
				
	  if (alphaCut > alpha_cut) continue;
						
	  float tt2_clus_eta = E_vec_cluster2.pseudoRapidity();
	  float tt2_clus_phi = E_vec_cluster2.getPhi();
		      
	  TLorentzVector pho1, pho2, mesonlv;
		      
	  pho1.SetPtEtaPhiE(tt_clus_pt, tt_clus_eta, tt_clus_phi, tt_clus_energy);
	  pho2.SetPtEtaPhiE(tt2_clus_pt, tt2_clus_eta, tt2_clus_phi, tt2_clus_energy);
		      
	  if (pho1.DeltaR(pho2) > delR_cut) continue;
											
	  mesonlv = pho1 + pho2;
	  float pairInvMass = mesonlv.M();

	  //std::cout << pairInvMass << std::endl;
		      

	  if (fabs(mesonlv.M()) > 0.02)
	    { 
	      pairInvMassTotal->Fill(pairInvMass);
	    }

	} // closing for second loop cluster
    } // closing for first loop cluster
							
							
  _eventTree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloCalibEmc_eta::End(PHCompositeNode * topNode)
{
  if (topNode == 0 && f_temp)
    {
      cal_output->Close();
      f_temp->Close();
      delete f_temp;
      delete cal_output;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  
  cal_output->cd();
  //_eventTree->Write();
  cal_output->Write();
  cal_output->Close();
  delete cal_output;

  return Fun4AllReturnCodes::EVENT_OK;
}


//______________________________________________________________________________..
// get the inv mass histogram for pi0/eta meson after collecting information from Ttree
// for time being we can consider this as signal+background

void CaloCalibEmc_eta::Loop(std::vector<int> desired_triggers, std::vector<double> CutParms,
			    int nevts, TString _filename, TTree * intree)

{
  std::cout << "Running in loop mode to get invariant mass plot" << std::endl;
  
  // unpack cuts in the exact order as they were packed
  // cluster-level cuts
  float pt1_cut = static_cast<float>(CutParms[0]);
  float pt2_cut = static_cast<float>(CutParms[1]);
  float chi2_cut = static_cast<float>(CutParms[2]);
  
  // pair-wise cuts
  float alpha_cut = static_cast<float>(CutParms[3]);
  float delR_cut  = static_cast<float>(CutParms[4]);
  
  // event level cuts
  int max_nCluster_cut = static_cast<int>(CutParms[5]);
  int min_nCluster_cut = static_cast<int>(CutParms[6]);
  float zver_cut = static_cast<float>(CutParms[7]);


  // modify to use valid triggers only (get from input of loop)
  std::vector<int> valid_triggers;
  for (int trg : desired_triggers) {
    if (trg >= 0 && trg <= 63) {
      valid_triggers.push_back(trg);
    } 
  }

  // exit if there are no valid triggers
  if (valid_triggers.empty()) {
    std::cerr << "No valid triggers to check. Exiting..." << std::endl;
    return;
  }

  // Get TTree from the provided pointer or from file
  TTree *t1 = intree;
  TFile* file = nullptr;
  if (!t1)
    {
      file = TFile::Open(_filename.Data()); // open input-file
      if (!file || file->IsZombie())
	{
	  std::cerr << "Error opening file " << _filename << std::endl;
	  return;
	}
      // read event Tree
      t1 = dynamic_cast<TTree*>(file->Get("_eventTree"));
      if (!t1)
	{
	  std::cerr << "Error: TTree '_eventTree' not found in file " << _filename << std::endl;
	  file->Close();
	  return;
	}
    }

  // Set branch addresses
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nClusters", &_nClusters);
  t1->SetBranchAddress("_triggerVector", _triggerVector); // 64 sized array
  t1->SetBranchAddress("_vertex", _vertex);
  t1->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  t1->SetBranchAddress("_clusterPts", _clusterPts);
  t1->SetBranchAddress("_clusterEtas", _clusterEtas);
  t1->SetBranchAddress("_clusterPhis", _clusterPhis);
  t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  t1->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);
  t1->SetBranchAddress("_clusterChi2", _clusterChi2);

  Long64_t nEntries = t1->GetEntries();
  int nEvents = (nevts < 0 || nEntries < nevts) ? nEntries : nevts;

  // store photon-candidates that pass photon-level cuts.
  std::vector<TLorentzVector> _Photons;

  // Reserve memory based on max_number of cluster (we will revise this later)
  _Photons.reserve(max_nCluster_cut);

  // Loop over events
  for (Long64_t i = 0; i < nEvents; i++)
    {
      t1->GetEntry(i); // get i_th Event from Ntuples

      if (i % 1000 == 0) {std::cout << "Processing event " << i << std::endl;}

      // Event loop check for trigger
      bool triggered = false;
      for (int trg : valid_triggers) {
	if (_triggerVector[trg] == 1) {
	  triggered = true;
	  break;
	}
      }
			
      // if we do not find any trigger then, skip the rest of the loop
      if (!triggered) continue; 

      // event level cluster cut
      if ((_nClusters > max_nCluster_cut) || (_nClusters < min_nCluster_cut)) continue;
    
      float vertex_z = _vertex[2];
      if (std::abs(vertex_z) > zver_cut) continue;

      // modify the enough space for photon-vector
      _Photons.clear();
      _Photons.reserve(_nClusters);

      // Loop over all clusters
      for (int j = 0; j < _nClusters; j++)
	{
	  // pT cut (lower of pt1 and pt2)
	  if (_clusterPts[j] < std::min(pt1_cut, pt2_cut)) continue;
      
	  if (_clusterChi2[j] > chi2_cut) continue;
      
	  TLorentzVector candidate;
	  candidate.SetPtEtaPhiE(_clusterPts[j], _clusterEtas[j], _clusterPhis[j], _clusterEnergies[j]);
	  _Photons.push_back(candidate);
	}

      // Loop over all unique pairs of photon candidates (vectors)
      const size_t nCandidates = _Photons.size();
      for (size_t ii = 0; ii < nCandidates; ii++)
	{
	  const TLorentzVector &pho1 = _Photons[ii];
      
	  if (pho1.Pt() < pt1_cut) continue;
      
	  // single-photon histogram for eta and phi
	  h_etaphi_clus->Fill(pho1.Eta(), pho1.Phi());
      
	  for (size_t ij = ii + 1; ij < nCandidates; ij++)
	    {
	      const TLorentzVector &pho2 = _Photons[ij];
        
	      if (pho2.Pt() < pt2_cut) continue;

	      float alpha = std::abs(pho1.E() - pho2.E()) / (pho1.E() + pho2.E());
	      if (alpha > alpha_cut) continue;

	      if (pho1.DeltaR(pho2) > delR_cut) continue;

	      TLorentzVector mesonCandidate = pho1 + pho2;

	      // Fill the histograms with the invariant mass and other variables
	      pairInvMassTotal->Fill(mesonCandidate.M());
	      pairInvMassPtEta->Fill(mesonCandidate.M(), mesonCandidate.Pt(), mesonCandidate.Eta());

	      // Calculate del_phi (adjust to -pi to pi)
	      float deltaEta = pho1.Eta() - pho2.Eta();
	      float deltaPhi = pho1.Phi() - pho2.Phi();
        
	      if (deltaPhi > M_PI)  deltaPhi -= 2 * M_PI;
	      if (deltaPhi < -M_PI) deltaPhi += 2 * M_PI;
        
	      pairpTDelPhiDelEta->Fill(mesonCandidate.Pt(), deltaPhi, deltaEta);
	      DelR_pairpT_f->Fill(pho1.DeltaR(pho2), mesonCandidate.Pt());
      
	    } // end inner loop over photon candidates
	} // end outer loop over photon candidates
    } // end event loop
  
  // to clean (close) input file
  if (file)
    {
      file->Close();
      delete file;
    }
}



//______________________________________________________________________________..
// getting background (using inv mass histogram) for pi0/eta meson using Ttree
// this code is inspired (mostly based on) from Loop function that creates fore-ground/signal
void CaloCalibEmc_eta::Loop_Event_Mixing(std::vector<double> evt_mix_parms,
					 std::vector<int> desired_triggers, std::vector<double> CutParms,
					 int nevts, TString _filename, TTree * intree)
{
  std::cout << "Creating Combinatorial Background and  Foreground (signal)" << std::endl;
  
  // unpack cuts in the exact order as they were packed
  // cluster-level cuts
  float pt1_cut = static_cast<float>(CutParms[0]);
  float pt2_cut = static_cast<float>(CutParms[1]);
  float chi2_cut = static_cast<float>(CutParms[2]);
  
  // pair-wise cuts
  float alpha_cut = static_cast<float>(CutParms[3]);
  float delR_cut  = static_cast<float>(CutParms[4]);
  
  // event level cuts
  int max_nCluster_cut = static_cast<int>(CutParms[5]);
  int min_nCluster_cut = static_cast<int>(CutParms[6]);
  float zver_cut = static_cast<float>(CutParms[7]);

  // cuts for event mixing
  float diff_zver_cut = static_cast<float>(evt_mix_parms[0]);
  int window_size = static_cast<int>(evt_mix_parms[1]);
  int del_nClus_cut = static_cast<int>(evt_mix_parms[2]);

  // modify to use valid triggers only
  std::vector<int> valid_triggers;
  for (int trg : desired_triggers) {
    if (trg >= 0 && trg <= 63) {
      valid_triggers.push_back(trg);
    } 
  }

  // exit if there are no valid triggers
  if (valid_triggers.empty()) {
    std::cerr << "No valid triggers to check." << std::endl;
    return;
  }

  // Open TTree from input file.
  TTree *t1 = intree;
  TFile *file = nullptr;
  if (!t1)
    {
      file = TFile::Open(_filename.Data()); // load the input root file
      if (!file || file->IsZombie())
	{
	  std::cerr << "Error opening file " << _filename << std::endl;
	  return;
	}
      
      t1 = dynamic_cast<TTree*>(file->Get("_eventTree")); // get Ttree Ntuple
      if (!t1)
	{
	  std::cerr << "Error: TTree '_eventTree' not found in file " << _filename << std::endl;
	  file->Close();
	  return;
	}
    }

  // Set branch addresses for NTuples
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nClusters", &_nClusters);
  t1->SetBranchAddress("_triggerVector", _triggerVector); // 64 sized array
  t1->SetBranchAddress("_vertex", _vertex);
  t1->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  t1->SetBranchAddress("_clusterPts", _clusterPts);
  t1->SetBranchAddress("_clusterEtas", _clusterEtas);
  t1->SetBranchAddress("_clusterPhis", _clusterPhis);
  t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  t1->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);
  t1->SetBranchAddress("_clusterChi2", _clusterChi2);

  // get total event number
  Long64_t nEntries = t1->GetEntries();
  int nevtsToProcess = (nevts < 0 || nEntries < nevts) ? nEntries : nevts;

  // Structure to hold the event’s data (z-vertex and 4-vector)
  struct EventData {
    int originalIndex;
    float vertex_z;
    int nClus;
    std::vector<TLorentzVector> clus_candidate;
  };

  // function to select event and build the window
  // we will have all trigger based cuts and event level cuts in this function
  auto SelectAndBuild = [&](int entry, EventData& out)-> bool {    
    // get event related to entry from NTuple 
    // we will return false if anything is not satisfied
    t1->GetEntry(entry);
    
    // logic for trigger
    bool hasTrig = false;
    for(int trg : valid_triggers) {
      if(_triggerVector[trg] == 1) { hasTrig = true; break; }
    }   
    if(!hasTrig) return false;
    
    // new we will have event-level cuts
    // first with nCluster related cuts
    if((_nClusters > max_nCluster_cut) || (_nClusters < min_nCluster_cut)) return false;
    	
    // z-vertex cuts        
    float vertex_z = _vertex[2];
    if(std::abs(vertex_z) > zver_cut) return false;
            
    // build cluster candidates
    std::vector<TLorentzVector> list;
    list.reserve(_nClusters);

    // getting into all clusters for particular event        
    for(int j = 0; j < _nClusters; ++j) {
                
      if(_clusterPts[j] < std::min(pt1_cut, pt2_cut)) continue;
      if(_clusterChi2[j] > chi2_cut) continue;
                
      TLorentzVector v;
      v.SetPtEtaPhiE(_clusterPts[j], _clusterEtas[j], _clusterPhis[j], _clusterEnergies[j]);
      list.push_back(v);
    }
    if(list.empty()) return false;

    // now update everything to the structure EventData
    out.originalIndex   = entry;
    out.vertex_z    = vertex_z;
    out.nClus = _nClusters;
    out.clus_candidate = std::move(list);
    return true;
  };

  // defining sliding window to hold as many events as required
  // as we move along each event (as first loop is done) we keep adding newer event and remove older

  // lets start to build sliding window (this will be the starting point)
  std::deque<EventData> window;
  int nextEntry = 0; // initiate the event entry

  // either we fill complete window or we are not left with any event
  while(((int)window.size() < window_size) && (nextEntry < nevtsToProcess)) {

    EventData ed;
    // update to the window
    if(SelectAndBuild(nextEntry, ed)) {window.push_back(std::move(ed));}
		
    // update to next entry
    nextEntry++;
  }

  // make signal/foreground and background using a moving window
  
  // get until we are done with all events
  while(!window.empty()) { 
            
    // this is the primary entry (or first entry)
    const EventData prim = window.front(); 
           
    // primary with primary gives signal and primary with others gets background
    // "prim" is the reference event and "mix" is moving event
    for(const auto& mix : window) {
                
      bool isSignal = (mix.originalIndex == prim.originalIndex); // boolean for signal
      if(!isSignal) {
	// background-specific cuts
        if(std::abs(prim.vertex_z - mix.vertex_z) > diff_zver_cut) continue;
        if(std::abs(prim.nClus - mix.nClus) > del_nClus_cut) continue;
      }

      // cluster-pair loops with in and with other events
      for (std::size_t i1 = 0; i1 < prim.clus_candidate.size(); ++i1) {

	const auto& pho1 = prim.clus_candidate[i1];
				
	if(pho1.Pt() < pt1_cut) continue; // first pT cut(s)

	// decide where our inner loop begins
	// for signal we do not want to double count but for background there is no possibility of double counting
	std::size_t jStart = 0;
	if (isSignal) {jStart = i1 + 1;} // setting correct place for starting location

	// starting inner loop (from correct location)
	for (std::size_t i2 = jStart; i2 < mix.clus_candidate.size(); ++i2) {

	  const auto& pho2 = mix.clus_candidate[i2];
					
	  if(pho2.Pt() < pt2_cut) continue; // second pT cut(s)
					
	  // energy assymetry cut
          float alpha = std::abs(pho1.E() - pho2.E()) / (pho1.E() + pho2.E());
          if(alpha > alpha_cut) continue;
          
	  // delR cut(s)
	  if (pho1.DeltaR(pho2) > delR_cut) continue;

          // combine two photon pairs to make meson candidate
	  TLorentzVector MesonCandidate = pho1 + pho2;
		
	  double deltaPhi = pho1.Phi() - pho2.Phi();
	  double deltaEta = pho1.Eta() - pho2.Eta();
          
	  if (deltaPhi > M_PI)  deltaPhi -= 2 * M_PI;
	  if (deltaPhi < -M_PI) deltaPhi += 2 * M_PI;
          
	  double pairPt = MesonCandidate.Pt();
	  double deltaR = pho1.DeltaR(pho2);
	  if(MesonCandidate.M() < 0.02) continue;
          
	  // fill all histogram for all signal/foreground
	  if(isSignal) {
	    pairInvMassTotal->Fill(MesonCandidate.M());
	    pairInvMassPtEta->Fill(MesonCandidate.M(), MesonCandidate.Pt(), MesonCandidate.Eta());
	    pairInvMassPtdelR->Fill(MesonCandidate.M(), MesonCandidate.Pt(), deltaR);
	    pairpTDelPhiDelEta->Fill(MesonCandidate.Pt(), deltaPhi, deltaEta);
	    DelR_pairpT_f->Fill(deltaR, pairPt);
          }

	  // fill all histogram for event mixing background
	  else {
	    pairInvMassTotalBkgd->Fill(MesonCandidate.M());
	    pairInvMassPtEtaBkgd->Fill(MesonCandidate.M(), MesonCandidate.Pt(), MesonCandidate.Eta());
	    pairInvMassPtdelRBkgd->Fill(MesonCandidate.M(), MesonCandidate.Pt(), deltaR);
	    pairpTDelPhiDelEtaBkgd->Fill(pairPt, deltaPhi, deltaEta);
	    DelR_pairpT_b->Fill(deltaR, pairPt);
          }
        
	} // second photon loop
      
      } // first photon loop
            
    } // slide window loop
    
    // remove primary first event in slide-window
    window.pop_front();
    
    // Refill one new passing event if available
    while(((int)window.size() < window_size) && (nextEntry < nevtsToProcess)) {
      EventData ed;
        
      if(SelectAndBuild(nextEntry, ed)) {window.push_back(std::move(ed));} // update next eligible event
      nextEntry++;
				
      if(nextEntry%10000 == 0){std::cout << nextEntry << " th event running." << std::endl;}
    }
	
  } // end of event
	
  // clean/remove file	
  if(file) { 
    file->Close(); 
    delete file; 
  }

}


//______________________________________________________________________________..
// getting background (using inv mass histogram) for pi0/eta meson using Ttree
// this code is inspired (mostly based on) from Loop function that creates fore-ground
void CaloCalibEmc_eta::Loop_Event_Mixing_Old(std::vector<double> evt_mix_parms,
					     std::vector<int> desired_triggers, std::vector<double> CutParms,
					     int nevts, TString _filename, TTree * intree)
{
  std::cout << "Creating Combinatorial Background based on Event Mixing. This will create Signal/Foreground as well." << std::endl;

  // unpack cuts in the exact order as they were packed
  // cluster-level cuts
  float pt1_cut = static_cast<float>(CutParms[0]);
  float pt2_cut = static_cast<float>(CutParms[1]);
  float chi2_cut = static_cast<float>(CutParms[2]);
  
  // pair-wise cuts
  float alpha_cut = static_cast<float>(CutParms[3]);
  float delR_cut  = static_cast<float>(CutParms[4]);
  
  // event level cuts
  int max_nCluster_cut = static_cast<int>(CutParms[5]);
  int min_nCluster_cut = static_cast<int>(CutParms[6]);
  float zver_cut = static_cast<float>(CutParms[7]);

  // cuts for event mixing
  float diff_zver_cut = static_cast<float>(evt_mix_parms[0]);
  int bckgnd_evnts = static_cast<int>(evt_mix_parms[1]);
  int del_nClus_cut = static_cast<int>(evt_mix_parms[2]);

  // modify to use valid triggers only
  std::vector<int> valid_triggers;
  for (int trg : desired_triggers) {
    if (trg >= 0 && trg <= 63) {
      valid_triggers.push_back(trg);
    } 
  }

  // exit if there are no valid triggers
  if (valid_triggers.empty()) {
    std::cerr << "No valid triggers to check." << std::endl;
    return;
  }

  // making vectors that has TLorentz vectors
  std::vector<TLorentzVector> savClusLV1;
  std::vector<TLorentzVector> savClusLV2;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TTree * t1 = intree;
  if (!intree)
    {
      TFile *f = new TFile(_filename);
      t1 = (TTree *) f->Get("_eventTree");
    }

  // Set Branches
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nClusters", &_nClusters);
  t1->SetBranchAddress("_triggerVector", _triggerVector); // 64 sized array
  //t1->SetBranchAddress("_nCentrality", &_nCentrality);
  //t1->SetBranchAddress("_clusterIDs", _clusterIDs);
  t1->SetBranchAddress("_vertex", _vertex); // contains co-ordinates of global-vertex	
  t1->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  t1->SetBranchAddress("_clusterPts", _clusterPts);
  t1->SetBranchAddress("_clusterEtas", _clusterEtas);
  t1->SetBranchAddress("_clusterPhis", _clusterPhis);
  t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  t1->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);
  t1->SetBranchAddress("_clusterChi2", _clusterChi2);

  // getting total number of events
  int nEntries = (int) t1->GetEntries();
  int nevts2 = nevts;

  if (nevts < 0 || nEntries < nevts) nevts2 = nEntries;
	

  for (int i1 = 0; i1 < nevts2; i1++) // iterating for each/every event(s)
    {
      // load the i1_th instance of the TTree
      t1->GetEntry(i1);

      if (i1 % 1000 == 0) {std::cout << "event number = " << i1 << std::endl;}

      // Event loop check for trigger
      bool triggered = false;
      for (int trg : valid_triggers) {
			
	if (_triggerVector[trg] == 1) {
	  triggered = true;
	  break;
	}
      }
			
      // if we do not find any trigger then, skip the rest of the loop
      if (!triggered) continue;

      // ncluster cuts (this can always be treated similar to centrality cut)
      int nClusters1 = _nClusters;
      if ((nClusters1 > max_nCluster_cut) || (nClusters1 < min_nCluster_cut)) continue; // will adjust this to look either peripheral or central event

      float vertex_z1 = _vertex[2]; // co-ordinates are 0, 1, 2 for "x" "y" and "z" respectively
      if (abs(vertex_z1) > zver_cut) continue;

      /*
	float Cent1 = _nCentrality; // centraity of the event
	if ((Cent1 < min_cent_cut) || (Cent1 > max_cent_cut)) continue; // centrality cut
      */

      // saving all cluster in vector of TLorentz Vector
      savClusLV1.clear();
      for (int i = 0; i < nClusters1; i++) // for each cluster in event
	{
	  if (_clusterPts[i] < std::min(pt1_cut, pt2_cut)) continue;
	  if (_clusterChi2[i] > chi2_cut) continue; 

	  TLorentzVector LV_Clus;

	  LV_Clus.SetPtEtaPhiE(_clusterPts[i], _clusterEtas[i], _clusterPhis[i], _clusterEnergies[i]);
	  savClusLV1.push_back(LV_Clus);
	}

      for (int i2 = i1 - bckgnd_evnts; i2 < i1 + bckgnd_evnts; i2++) // we take some events up and down each wrt to event of our concern
	{
	  if (i2 < 0) continue; // discard negative event number (which does not exist at all)
	  if (i2 > nevts2) continue; // when max event number is reached, we exit our iteration
	
	  t1->GetEntry(i2);

	  // Event loop check for trigger
	  bool triggered2 = false;
	  for (int trg : valid_triggers) {
			
	    if (_triggerVector[trg] == 1) {
	      triggered2 = true;
	      break;
	    }
	  }
			
	  // if we do not find any trigger then, skip the rest of the loop
	  if (!triggered2) continue;

	  int nClusters2 = _nClusters;
	  if ((nClusters2 > max_nCluster_cut) || (nClusters2 < min_nCluster_cut)) continue;

	  float vertex_z2 = _vertex[2]; // co-ordinates are 0, 1, 2
	  if (abs(vertex_z2) > zver_cut) continue;

	  if (abs(vertex_z1 - vertex_z2) > diff_zver_cut) continue; // ignore those events if their z-vertex is far from 5 units
				
	  // looking for similar cluster (not significantly relavant for p+p collision)
	  if (abs(nClusters2 - nClusters1) > del_nClus_cut) continue;

	  /*
	    float Cent2 = _nCentrality; // centraity of the event
	    if ((Cent2 < min_cent_cut) || (Cent2 > max_cent_cut)) continue; // centrality cut
				
	    if (abs(Cent1-Cent2) > cent_range) continue; // we use two events with almost same centrality
	  */

	  // saving all cluster in vector of TLorentz Vector
	  savClusLV2.clear();
	  for (int j = 0; j < nClusters2; j++) // for each cluster in event
	    {
	      if (_clusterPts[j] < std::min(pt1_cut, pt2_cut)) continue;
	      if (_clusterChi2[j] > chi2_cut) continue; 

	      TLorentzVector LV_Clus2;

	      LV_Clus2.SetPtEtaPhiE(_clusterPts[j], _clusterEtas[j], _clusterPhis[j], _clusterEnergies[j]);
	      savClusLV2.push_back(LV_Clus2);
	    }

			
	  for (size_t iCs = 0; iCs < savClusLV1.size(); iCs++) // outer loop to iterate clusters
	    {
	      TLorentzVector& pho1 = savClusLV1[iCs];

	      if (fabs(pho1.Pt()) < pt1_cut)	continue; // first P_T cut
						
	      for (size_t jCs = iCs+1; jCs < savClusLV2.size(); jCs++) // inner loop to iterate clusters
		{
		  TLorentzVector& pho2 = savClusLV2[jCs];

		  if (fabs(pho2.Pt()) < pt2_cut) continue; //second P_T cut
								
		  alphaCut = fabs((pho1.E() - pho2.E())/(pho1.E()+ pho2.E()));
		  if (alphaCut > alpha_cut) continue;	
								
		  //if ((pho1.DeltaR(pho2) > delR_cut) || (pho1.DeltaR(pho2) < 0.05))  continue;
		  if (pho1.DeltaR(pho2) > delR_cut)  continue;

		  TLorentzVector mesonlv = pho1 + pho2;

		  if ((mesonlv.M() < 0.02) || (mesonlv.Pt() < 1.0)) continue;

		  double _delphi = pho1.Phi() - pho2.Phi();
		  double _deleta = pho1.Eta() - pho2.Eta();
		  double _pairpT = mesonlv.Pt();
		  double _delR = pho1.DeltaR(pho2);

		  if (_delphi > M_PI) {_delphi -= 2 * M_PI;}       
		  if (_delphi < -1 * M_PI) {_delphi += 2 * M_PI;}
								
		  if (i1 == i2) // histogram for foreground
		    {

		      pairInvMassTotal->Fill(mesonlv.M());
		      pairInvMassPtEta->Fill(mesonlv.M(), mesonlv.Pt(), mesonlv.Eta());
		      pairInvMassPtdelR->Fill(mesonlv.M(), mesonlv.Pt(), _delR);

		      pairpTDelPhiDelEta->Fill(mesonlv.Pt(), _delphi, _deleta); // for foreground
		      DelR_pairpT_f->Fill(_delR, _pairpT);

		    }

		  else if (i1 != i2) // histogram for background
		    {

		      DelR_pairpT_b->Fill(_delR, _pairpT);

		      pairInvMassTotalBkgd->Fill(mesonlv.M());
		      pairInvMassPtEtaBkgd->Fill(mesonlv.M(), mesonlv.Pt(), mesonlv.Eta());
		      pairInvMassPtdelRBkgd->Fill(mesonlv.M(), mesonlv.Pt(), _delR);
		      pairpTDelPhiDelEtaBkgd->Fill(_pairpT, _delphi, _deleta);	// for background

		    }
		  else {std::cout << "Something is not good for event mixing function." << std::endl;} 
									
		} // inner loop (cluster)
	    } // outer loop (cluster)
	} // second-set of events (i2) (when i1 == i2, it makes foreground)
    }// iterating over all events (first set :: i1)
}


// this function is useful for position swapping method of background reconstruction
// we have a flag that will make either signal (foreground) or swapped background, so this function will do either of that

void CaloCalibEmc_eta::Loop_Position_Swapping(std::vector<double> pos_swap_parms,
					      std::vector<int> desired_triggers, std::vector<double> CutParms,
					      int nevts, TString _filename, TTree * intree)
{
  std::cout << "Creating Position Swapping and/or Foreground (signal)" << std::endl;
  
  // unpack cuts in the exact order as they were packed
  // cluster-level cuts
  float pt1_cut = static_cast<float>(CutParms[0]);
  float pt2_cut = static_cast<float>(CutParms[1]);
  float chi2_cut = static_cast<float>(CutParms[2]);
  
  // pair-wise cuts
  float alpha_cut = static_cast<float>(CutParms[3]);
  float delR_cut  = static_cast<float>(CutParms[4]);
  
  // event level cuts
  int max_nCluster_cut = static_cast<int>(CutParms[5]);
  int min_nCluster_cut = static_cast<int>(CutParms[6]);
  float zver_cut = static_cast<float>(CutParms[7]);

  // parameters for position swapping
  bool doSwapBkg = static_cast<bool>(pos_swap_parms[0]);
  int n_Swap = static_cast<int>(pos_swap_parms[1]);
  float delR_bkg = static_cast<float>(pos_swap_parms[2]);
  float delE_bkg = static_cast<float>(pos_swap_parms[3]);

  
  // modify to use valid triggers only
  std::vector<int> valid_triggers;
  for (int trg : desired_triggers) {
    if (trg >= 0 && trg <= 63) {
      valid_triggers.push_back(trg);
    } 
  }

  // exit if there are no valid triggers
  if (valid_triggers.empty()) {
    std::cerr << "No valid triggers to check. Exiting..." << std::endl;
    return;
  }

  // use random number generator to randomly selecting photons
  struct timeval tv;
  gettimeofday(&tv, nullptr);

  // use only the microseconds field:
  UInt_t seed = static_cast<UInt_t>(tv.tv_usec);
  TRandom3 rnd(seed);


  // Open TTree from input file.
  TTree *t1 = intree;
  TFile *file = nullptr;
  if (!t1)
    {
      file = TFile::Open(_filename.Data()); // load the input root file
      if (!file || file->IsZombie())
	{
	  std::cerr << "Error opening file " << _filename << std::endl;
	  return;
	}
      t1 = dynamic_cast<TTree*>(file->Get("_eventTree")); // get Ttree Ntuple
      if (!t1)
	{
	  std::cerr << "Error: TTree '_eventTree' not found in file " << _filename << std::endl;
	  file->Close();
	  return;
	}
    }

  // Set branch addresses for NTuples
  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nClusters", &_nClusters);
  t1->SetBranchAddress("_triggerVector", _triggerVector); // 64 sized array
  t1->SetBranchAddress("_vertex", _vertex);
  t1->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  t1->SetBranchAddress("_clusterPts", _clusterPts);
  t1->SetBranchAddress("_clusterEtas", _clusterEtas);
  t1->SetBranchAddress("_clusterPhis", _clusterPhis);
  t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  t1->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);
  t1->SetBranchAddress("_clusterChi2", _clusterChi2);
  //t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  //t1->SetBranchAdderss("_maxTowerPhis", _maxTowerPhis);

  // getting total number of events
  Long64_t nEntries = t1->GetEntries();
  int nEvents = (nevts < 0 || nEntries < nevts) ? nEntries : nevts;

  // store photon-candidates that pass photon-level cuts.
  std::vector<TLorentzVector> _Photons;

  // Reserve memory based on max_number of cluster (we will revise this later)
  _Photons.reserve(max_nCluster_cut);

  // Loop over events
  for (int i = 0; i < nEvents; i++)
    {
      t1->GetEntry(i); // get i_th Event from Ntuples

      if (i % 1000 == 0) {std::cout << "Processing event " << i << std::endl;}

      // Event loop check for trigger
      bool triggered = false;
      for (int trg : valid_triggers) {
	if (_triggerVector[trg] == 1) {
	  triggered = true;
	  break;
	}
      }
			
      // if we do not find any trigger then, skip the rest of the loop
      if (!triggered) continue; 

      // event level cluster cut
      if ((_nClusters > max_nCluster_cut) || (_nClusters < min_nCluster_cut)) continue;
    
      float vertex_z = _vertex[2];
      if (std::abs(vertex_z) > zver_cut) continue;

      // modify the enough space for photon-vector
      _Photons.clear();
      _Photons.reserve(_nClusters);

      // Loop over all clusters
      for (int j = 0; j < _nClusters; j++)
	{
	  // pT cut (lower of pt1 and pt2)
	  if (_clusterPts[j] < std::min(pt1_cut, pt2_cut)) continue;
      
	  if (_clusterChi2[j] > chi2_cut) continue;
      
	  TLorentzVector candidate;
	  candidate.SetPtEtaPhiE(_clusterPts[j], _clusterEtas[j], _clusterPhis[j], _clusterEnergies[j]);
	  _Photons.push_back(candidate);
			
	}
			
      // we will ignore if we can not make any pair-wise cluster
      const size_t nPho = _Photons.size();
      if (nPho < 2) continue;


      // now we will use flag either to make signal or background
      if (!doSwapBkg) {

	// make signal (foreground)
	for (size_t i1 = 0; i1 < nPho; ++i1) {

	  const auto &pho1 = _Photons[i1];
	  if (pho1.Pt() < pt1_cut) continue; // first pT cut
        
	  for (size_t j1 = i1+1; j1 < nPho; ++j1) {
          
	    const auto &pho2 = _Photons[j1];
	    if (pho2.Pt() < pt2_cut) continue; // second pT cut
          
	    // energy assymetry cut
	    float alpha = std::abs(pho1.E() - pho2.E()) / (pho1.E() + pho2.E());
	    if (alpha > alpha_cut) continue;
          
	    // delR cuts
	    if (pho1.DeltaR(pho2) > delR_cut) continue;
          
	    TLorentzVector mesonCandidate = pho1 + pho2;

	    // Fill the histograms with the invariant mass and other variables
	    pairInvMassTotal->Fill(mesonCandidate.M());
	    pairInvMassPtEta->Fill(mesonCandidate.M(), mesonCandidate.Pt(), mesonCandidate.Eta());

	    // Calculate del_phi (adjust to -pi to pi)
	    float deltaEta = pho1.Eta() - pho2.Eta();
	    float deltaPhi = pho1.Phi() - pho2.Phi();
        
	    if (deltaPhi > M_PI)  deltaPhi -= 2 * M_PI;
	    if (deltaPhi < -M_PI) deltaPhi += 2 * M_PI;
        
	    pairpTDelPhiDelEta->Fill(mesonCandidate.Pt(), deltaPhi, deltaEta);
	    DelR_pairpT_f->Fill(pho1.DeltaR(pho2), mesonCandidate.Pt());
				
	  } // inner loop for photons
	} // outer loop for photons

      } // Swap flag 

      else {

	// we will ignore if we can not make any pair-wise cluster
	size_t nPho2 = _Photons.size();
	if (nPho2 < static_cast<size_t>(n_Swap)) continue;

	// to make position swapped background
	for (size_t is = 0; is < nPho; ++is) {
        
	  const auto &phob1 = _Photons[is];
	  if (phob1.Pt() < pt1_cut) continue;

	  for (size_t js = is+1; js < nPho; ++js) {
					 
	    // perform n_Swap independent swaps
	    // we do not make any cuts or similar things before we swap
	    for (int swapIter = 0; swapIter < n_Swap; ++swapIter) {
              
	      auto phob2 = _Photons[js];  // fresh copy each iteration
              
	      // Pick random k != i, j
	      size_t swap_k = rnd.Integer(nPho);
	      while (swap_k == is || swap_k == js) {
		swap_k = rnd.Integer(nPho);
	      }

	      // delR and delE cuts between j and k before making swapping
	      float dR_jk = _Photons[js].DeltaR(_Photons[swap_k]);
	      float dE_jk = std::abs(_Photons[js].E() - _Photons[swap_k].E()) / (_Photons[js].E() + _Photons[swap_k].E());

	      if ((dR_jk <= delR_bkg) || (dE_jk <= delE_bkg)) continue;

	      // Perform position swap: swap E and keep everything else same on "j"
	      double newE = _Photons[swap_k].E();
	      double newEta  = _Photons[js].Eta();
	      double newPhi  = _Photons[js].Phi();
	      double newpT = newE / std::cosh(newEta);

	      phob2.SetPtEtaPhiE(newpT, newEta, newPhi, newE);

	      // second pT cuts
	      if (phob2.Pt() < pt2_cut) continue;
					
	      // energy assymtery cuts
	      float alpha2 = std::abs(phob1.E() - phob2.E()) / (phob1.E() + phob2.E());
	      if (alpha2 > alpha_cut) continue;
					
	      // delR cuts
	      if (phob1.DeltaR(phob2) > delR_cut) continue;

          
	      TLorentzVector mesonCandidate = phob1 + phob2;

	      // Fill the histograms with the invariant mass and other variables
	      pairInvMassTotalBkgd->Fill(mesonCandidate.M());
	      pairInvMassPtEtaBkgd->Fill(mesonCandidate.M(), mesonCandidate.Pt(), mesonCandidate.Eta());
	      pairInvMassPtdelRBkgd->Fill(mesonCandidate.M(), mesonCandidate.Pt(), phob1.DeltaR(phob2)); 
					
	      // Calculate del_phi (adjust to -pi to pi)
	      float deltaEta = phob1.Eta() - phob2.Eta();
	      float deltaPhi = phob1.Phi() - phob2.Phi();
        
	      if (deltaPhi > M_PI)  deltaPhi -= 2 * M_PI;
	      if (deltaPhi < -M_PI) deltaPhi += 2 * M_PI;
        
	      pairpTDelPhiDelEtaBkgd->Fill(mesonCandidate.Pt(), deltaPhi, deltaEta);
	      DelR_pairpT_b->Fill(phob1.DeltaR(phob2), mesonCandidate.Pt());
					 
	    } // loop for swapping background

	  } // inner loop
      
	} // outer loop

      } // swap flag for background
  
    } // end event loop


  // Clean up
  if (file) {
    file->Close();
    delete file;
  }
}


