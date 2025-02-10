#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic pop

//Fun4all stuff
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <ffaobjects/EventHeader.h>

//truth information
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

//ROOT stuff
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <vector>
#include <TTree.h>
#include <TMath.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <unordered_map>

#include <CLHEP/Geometry/Point3D.h>

#include "pi0ClusterAna.h"

//____________________________________________________________________________..
pi0ClusterAna::pi0ClusterAna(const std::string &name, const std::string &outName = "pi0ClusterAnaOut"):
  SubsysReco(name)
  , n_event(0)
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
      _eventTree->Branch("_parent_id", _parent_id, "_parent_id[_nFourVector]/F"); // parent id (0 for primary particle)and track id will help to track particle(s) together
      _eventTree->Branch("_primary_id", _primary_id, "_primary_id[_nFourVector]/F"); // this is same for both primary and secondary particles
      _eventTree->Branch("_track_id", _track_id, "_track_id[_nFourVector]/F"); // save track id to check if two particle have same immediate parent
      _eventTree->Branch("_pid_particle", _pid_particle, "_pid_particle[_nFourVector]/F"); // this is pid of particle for both primary and secondary particles
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

  // we can embed multiple events and any number of them but make sure to embed eta meson at first
  // since, we do not have eta meson as primary particle we can track start of embedding tracking first eta meson
  // in macro we have to keep eta meson at the very last, that makes sure eta meson comes first here
  // primary id is same for primary and all secondary particle(s) means they belong to same parent (primary particle)
  // to see parent particle and daughter particle connection in secondary particle, we match "track id" of parent particle with "parent id" of daughter particle
	

  // save the information so that we can re-trace it later
  // we can save only track id and flag to check if it is embedded or not
  // this will be for all particles, wether primary or secondary
  // we want to have final listing of track id
  std::vector<int> ref_track_id; // this saves track id of all tracks (primary + secondary)
  std::vector<int> _embbed_flag; // flag to indicate if it is embeded or not
  std::vector<int> _primary_id_PRI; // number to save primary id of primary particle
	
  // clear vector to save information later
  ref_track_id.clear();
  _embbed_flag.clear();
  _primary_id_PRI.clear();

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we will save track id in two places (one to use to scan secondary particles and another to scan all particles at the end)
  std::vector<int> primary_track_id_temp;	
  primary_track_id_temp.clear();
	
  // SCANNING PRIMARY PARTICLE RANGE
	 
  // get primary particle
  PHG4TruthInfoContainer::Range truthRange = truthinfo->GetPrimaryParticleRange();
  PHG4TruthInfoContainer::ConstIterator truthIter;

  PHG4Particle *truth_pri;
	
  bool eta_flag = false; // flag to identify first eta meson (start of embedded particles)
  for(truthIter = truthRange.first; truthIter != truthRange.second; truthIter++)
    {
      truth_pri = truthIter->second;
			
      if ((truth_pri->get_parent_id()) != 0) {std::cout << "Secondary Particle Spotted" << std::endl; continue;}
			
      // update track id for primary particles
      ref_track_id.push_back(truth_pri->get_track_id());
      primary_track_id_temp.push_back(truth_pri->get_track_id()); // to use while working with secondary particle
			
      // save primary id to be used later (have same dimension as embedd flag)
      _primary_id_PRI.push_back(truth_pri->get_primary_id()); 

      // since we do not have eta meson as primary particle (so once eta meson is found it indicates the start of embedding particles)
      // we have to make sure we embed eta meson at first (just to simplify things)
      if (((truth_pri->get_pid()) == 221) || (eta_flag == true))
	{
	  eta_flag = true ; // once we enterd into this loop we just keep it as this is the start of embeding
	  _embbed_flag.push_back(1); // indicating embedding
	}
      else {_embbed_flag.push_back(0);} // for hijng events								
    }
	
  std::cout << "Size of primary particle array is = " << _primary_id_PRI.size() << std::endl;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  // SCANNING SECONDARY PARTICLE RANGE
	
  // we will first scan over all secondary particles and save few things like pid, parent id and track id of secondary particles
  PHG4TruthInfoContainer::Range truthRange_sec = truthinfo->GetSecondaryParticleRange();
  PHG4TruthInfoContainer::ConstIterator truthIter_sec;

  PHG4Particle *truth_sec;
	
  // vectors to store the information for all secondary particles
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

  // Using unordered_map to store elements of _parentid_sec and their indices
  std::unordered_map<int, std::vector<size_t>> map_parentid_sec;

  // Populate map with parent IDs and their corresponding indices
  for (size_t i = 0; i < _parentid_sec.size(); ++i) 
    {
      map_parentid_sec[_parentid_sec[i]].push_back(i);
    }

  // Vectors to store track IDs (temporary and final)
  std::vector<int> main_trackid_sec;       // Track IDs currently being processed
  std::vector<int> temporary_trackid_sec;  // Temporary storage for new track IDs

  // Iterate over each primary track ID
  for (int _Primary_Track_ID : primary_track_id_temp) 
    {
      main_trackid_sec.clear();
      main_trackid_sec.push_back(_Primary_Track_ID); // Start with the primary particle

      while (!main_trackid_sec.empty())	// Continue until no new track IDs are found
	{  				
	  temporary_trackid_sec.clear();

	  // Process each track ID in the current main list
	  for (int track_id : main_trackid_sec) 
	    {
	      auto it = map_parentid_sec.find(track_id);
	      if (it != map_parentid_sec.end()) 
		{  
		  // If track_id has secondary particles
		  for (size_t index : it->second) // to get index (multiple indices exists)
		    {
		      // Always save the secondary particle
		      // if we wish to save only first electron photon and positron then we can keep these inside if statement
		      // right now we are saving track id at final vector every time until we get 22, 11 or -11
		      if(AllMode){ref_track_id.push_back(_trackid_sec[index]);} // save all generation of secondary particle
									
		      // If particle is photon (22), electron (11), or positron (-11), do NOT track further
		      // Skip tracking deep further, but continue processing other particles
		      if (_pid_sec[index] == 22 || _pid_sec[index] == 11 || _pid_sec[index] == -11)
			{
			  if (!AllMode){ref_track_id.push_back(_trackid_sec[index]);} // save only final 11, 22 and -11 of secondary particle
			  continue;	
			}

		      // continue tracking this particle if it not photon
		      temporary_trackid_sec.push_back(_trackid_sec[index]);
		    }
		}
	    }
				
	  // clear the main_trackid_sec vector before you push temp vector
	  main_trackid_sec.clear();
					
	  // Move temporary_trackid_sec to main_trackid_sec for the next iteration
	  main_trackid_sec = std::move(temporary_trackid_sec);
	}
    }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  std::cout << "Size of total particles in an event is = " << ref_track_id.size() << std::endl;
	
  // do not save if there are sgnificantly large number of particles
  if ((ref_track_id.size() < 1) || (ref_track_id.size() > 9900)){

    std::cout << "The final size of the primary and secondary particles is/are = " << ref_track_id.size() << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  // EXTRACTING THE INFORMATION FROM TRACK_ID (SECONDARY PARTICLE)

  // now with the help of trackid we will extract all the information of primary and secondary particles and save it in TTree histogram	

  // vectors to store the information for ntuple
  std::vector<float> _PX;
  std::vector<float> _PY;
  std::vector<float> _PZ;
  std::vector<float> _E;
  std::vector<float> _Eta;
  std::vector<int> _PARENT_ID;
  std::vector<int> _PRIMARY_ID;
  std::vector<int> _TRACK_ID;
  std::vector<int> _PAR_PID;
  std::vector<int> _EMBEDDING_ID;
  
  // preparing all vectors 
  _PX.clear();
  _PY.clear();
  _PZ.clear();
  _E.clear();
  _Eta.clear();
  _PARENT_ID.clear();
  _PRIMARY_ID.clear();
  _TRACK_ID.clear();
  _PAR_PID.clear();
  _EMBEDDING_ID.clear();

  // if two size matches we will continue
  if (_embbed_flag.size() != _primary_id_PRI.size()) {std::cout << "primary id and embedded id vectors mismatch" << std::endl;}

  // Create a mapping before the loop
  std::unordered_map<int, int> primaryIdToEmbedFlag;

  for (size_t i = 0; i < _primary_id_PRI.size(); ++i) {
    primaryIdToEmbedFlag[_primary_id_PRI[i]] = _embbed_flag[i];
  }
	
  // looping over "ref_track_id" trackid and extracting all the essentail information
  for (unsigned int kk = 0; kk < ref_track_id.size(); kk++)
    {
      PHG4Particle *final_par = truthinfo->GetParticle(ref_track_id[kk]);
			
      if (final_par == nullptr) {std::cout << " null pointer" << std::endl;}
		
      // some conditions to reduce the size of particle
      if (fabs(getEta(final_par)) > 1.1) continue; // condition of pseudorapidity
		
      float _truth_pT = std::sqrt((final_par->get_px()) * (final_par->get_px()) + (final_par->get_py()) *(final_par->get_py()));
      if (_truth_pT < 0.3) continue; // condition of pT cut
		
      // update all the vectors that are essentials
      _PX.push_back(final_par->get_px());
      _PY.push_back(final_par->get_py());
      _PZ.push_back(final_par->get_pz());
      _E.push_back(final_par->get_e());
      _Eta.push_back(getEta(final_par)); // updating pseudorapidity of secondary particle		
      _PARENT_ID.push_back(final_par->get_parent_id());
      _PRIMARY_ID.push_back(final_par->get_primary_id());
      _TRACK_ID.push_back(ref_track_id[kk]);
      _PAR_PID.push_back(final_par->get_pid()); // pid for secondary particle
	
      // find the location of primary id to eventually get embedded flag
      auto it = primaryIdToEmbedFlag.find(final_par->get_primary_id());
		
      if (it != primaryIdToEmbedFlag.end()){
	_EMBEDDING_ID.push_back(it->second);
      }
      else{
	_EMBEDDING_ID.push_back(-1); // this indicates some issue with embedding flag
      } 
		
    }
	
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  // SAVING IN TTREE
  // updating in TTree
  _eventNumber = n_event;
  _nFourVector = _E.size();	
		
  for (int kl = 0; kl < _nFourVector; kl++) // update all vectors
    {
      // update all the branches of the Ttree
      // if size of four vector is zero then, we expect there is nothing inside the branches of the Ttree
      _fv_energy[kl] = _E[kl];
      _fv_px[kl] = _PX[kl];
      _fv_py[kl] = _PY[kl];
      _fv_pz[kl] = _PZ[kl];
      _fv_Eta[kl] = _Eta[kl];
      _parent_id[kl] = _PARENT_ID[kl];
      _primary_id[kl] = _PRIMARY_ID[kl];
      _track_id[kl] = _TRACK_ID[kl];
      _pid_particle[kl] = _PAR_PID[kl];
      _embedding[kl] = _EMBEDDING_ID[kl];
		
      //std::cout << "primary id, track id, parent id , pid = " << _PRIMARY_ID[kl] << " , " << _TRACK_ID[kl] << " , " << _PARENT_ID[kl] << " , " << _PAR_PID[kl] << std::endl;
			 
    }

  n_event++; // updating event number
  _eventTree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int pi0ClusterAna::ResetEvent(PHCompositeNode *topNode)
{
  std::cout << "pi0ClusterAna::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;

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

  if (topNode == 0)
    {
      out->Close();
      delete out;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  out->cd();
  
  out->Write();
  out->Close();
  delete out; 

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

