#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TSystem.h"
#include "TVector3.h"

#include <iostream>
#include <cmath>

// Convert pseudorapidity (eta) to polar angle (theta)
double Eta2Theta(double etaVal) { 
  return 2.0 * atan(exp(-etaVal)); 
}

void super_fastmc_theta() {
  
  const char* outputName = "test.root";
  int numPions = 1e6;
  float maxPt = 20.0;
  int NumBins = 500;
  
  // Open a ROOT file in UPDATE mode.
  TFile* out = new TFile(outputName, "UPDATE");

  // Define simple flat functions for pt, rapidity, and phi.
  TF1* flatPt = new TF1("flatpt", "1", 0.02, maxPt);
  float MYPI = 3.14159;
  TF1* flat_rapidity = new TF1("flatRap", "1", -1, 1);
  TF1* flat_phi = new TF1("flatPhi", "1", 0, 2 * MYPI);
  TF1* real_theta = new TF1("flatTheta", "sin(x)", 0, MYPI);
  
  // Woods-Saxon type function (currently defined but not used in this macro)
  TF1 woodsSax("woodsSax", "1.0 - 1.0/(1.0 + exp((0.5*x - 0.25)/0.02))", 0.0, maxPt);
  woodsSax.SetNpx(20000);

  // Create histograms for pT spectra, invariant mass, and other correlations.
  TH1F* picorr = new TH1F("picorr", "corrected pi0 spect, according to fit values", NumBins, 0.0, maxPt);
  TH1F* gamcorr = new TH1F("gamcorr", "correct gamma spectrum according to pi0 fit values", NumBins, 0.0, maxPt);
  TH2F* invMass = new TH2F("invMassPt", "", 160, 0.0, 0.8, 4, 2.0, 6.0);
  invMass->Sumw2();
  
  // Clone histogram for pT distributions.
  TH1F* pipt = new TH1F(*picorr);
  pipt->Reset();
  
  TH3F* ptp_pt1_dang = new TH3F("ptp_pt1_dang", "", 50, 1.0, 20.0, 50, 1.0, 20.0, 600, 0, 0.600);
  
  TH1F* ratiopt = new TH1F(*pipt);
  ratiopt->SetName("ratiopt");
  std::cout << "processing " << ratiopt->GetName() << std::endl;
    
  TH1F* gampt = new TH1F(*pipt);
  TH1F* gampipt = new TH1F(*pipt);
  TH1F* ratiogampipt = new TH1F(*pipt);

  // Eta histograms (currently created but not further used)
  TH1F* hetapt = new TH1F(*pipt);
  TH1F* gametapt = new TH1F(*pipt);
  TH1F* ratiogametapt = new TH1F(*pipt);

  // Set names and titles for clarity.
  pipt->SetName("pipt");
  gampt->SetName("gampt");
  gampipt->SetName("gampipt");
  ratiogampipt->SetName("ratiogampipt");

  pipt->SetTitle("Pi0 Pt Dist");
  gampt->SetTitle("All Gamma Pt Dist");
  ratiopt->SetTitle("All Gamma / Pi0 Ratio (Pt)");
  gampipt->SetTitle("Pi0 Gam. Pt Dist");
  ratiogampipt->SetTitle("Gamma_Pi0 / Pi0 Ratio (Pt)");

  TH1F* thetaDist = new TH1F("thetaDist", "distribution of opening angles of the two gamma", 3140, 0, 3.14);
	  		
  // Define the pi0 and eta pT spectral functions.
  TF1* pi0pt = new TF1("pi0pt", "[0]*(x^0.5)/([1]+x)^[2]", 0.001, 10);
  pi0pt->SetParameters(1.0e11, 0.5, 8.0);
  
  TF1* etapt = new TF1("etapt", "[0]*x/([1]+x)^[2]", 0.001, 10);
  etapt->SetParameters(1.0e8, 1.1, 5.0);  // Changed 10^8 to 1.0e8 for proper exponentiation.

  // Define a resolution function and a Gaussian to simulate detector smearing.
  TF1* f1 = new TF1("f1reso", "sqrt([0]*[0]/x+[1]*[1])", 0.5, 20.);
  f1->SetParameter(0, 0.165);
  f1->SetParameter(1, 0.26);
	  
  TF1* fgaus = new TF1("mygaus", "gaus", -1.5, 1.5);
  fgaus->SetParameter(0, 1e6);
  fgaus->SetParameter(1, 0);

  TH2F* h2esm = new TH2F("h2esm", "", 80, 0, 20, 200, 0, 20);
  TH2F* h2resoe = new TH2F("h2resoe", "", 80, 0, 20, 200, -1, 1);

  // Main Monte Carlo loop over pions (and etas)
  for (int i = 0; i < numPions; i++) { 
    for (int ketapi = 0; ketapi < 2; ketapi++) {
      float massPi = 0.13497;
      if (ketapi == 1) { 
        massPi = 0.547;
      }
	      
      float thetaOrig = Eta2Theta(flat_rapidity->GetRandom());
      float ptPi = flatPt->GetRandom(); 
      float pPi = ptPi / sin(thetaOrig);
      float ePi = sqrt(pPi * pPi + massPi * massPi);

      if (ePi < 100) {
		
        // Build the pion 3-vector and 4-vector.
        TVector3 pion3v(2.0, 0.0, 0.0);
        pion3v.SetTheta(thetaOrig);	  
        pion3v.SetPhi(flat_phi->GetRandom());
        pion3v.SetMag(pPi);
		
        TLorentzVector pcylorig;
        pcylorig.SetVect(pion3v);
        pcylorig.SetE(ePi);	  
		
        // Simulate the two-photon decay in the pion’s rest frame.
        TLorentzVector photon1;
        TLorentzVector photon2;
        TLorentzVector* g1 = &photon1;
        TLorentzVector* g2 = &photon2;

        TVector3 photon3v(2.0, 0.0, 0.0);
        photon3v.SetMag(massPi / 2);
        photon3v.SetTheta(real_theta->GetRandom());
        photon3v.SetPhi(flat_phi->GetRandom());
        g1->SetVect(photon3v);
        g1->SetE(massPi / 2);
		
        g2->SetVect(-photon3v);
        g2->SetE(massPi / 2);
		
        // Boost photons from the rest frame to the lab frame.
        g1->Boost(pcylorig.BoostVector());
        g2->Boost(pcylorig.BoostVector());
		
        // Apply smearing to mimic detector resolution.
        TLorentzVector pho1, pho2, pi0lv;

        float tt_clus_energy = g1->E();
        float reso1 = f1->Eval(tt_clus_energy);
        fgaus->SetParameter(2, reso1);
        float smear1 = fgaus->GetRandom();
        tt_clus_energy += tt_clus_energy * smear1;
        h2esm->Fill(g1->E(), tt_clus_energy);
        h2resoe->Fill(g1->E(), tt_clus_energy / g1->E() - 1.0);
        float tt_clus_pt = g1->Pt();
        tt_clus_pt += tt_clus_pt * smear1;
        float tt_clus_phi = g1->Phi();
        float tt_clus_eta = g1->Eta();

        float tt_clus_energy2 = g2->E();
        float reso2 = f1->Eval(tt_clus_energy2);
        fgaus->SetParameter(2, reso2);
        float smear2 = fgaus->GetRandom();
        tt_clus_energy2 += tt_clus_energy2 * smear2;
        h2esm->Fill(g2->E(), tt_clus_energy2);
        h2resoe->Fill(g2->E(), tt_clus_energy2 / g2->E() - 1.0);
        float tt_clus_pt2 = g2->Pt();
        tt_clus_pt2 += tt_clus_pt2 * smear2;
        float tt_clus_phi2 = g2->Phi();
        float tt_clus_eta2 = g2->Eta();

        pho1.SetPtEtaPhiE(tt_clus_pt, tt_clus_eta, tt_clus_phi, tt_clus_energy);
        pho2.SetPtEtaPhiE(tt_clus_pt2, tt_clus_eta2, tt_clus_phi2, tt_clus_energy2);
		
        // Reconstruct the invariant mass.
        pi0lv = pho1 + pho2;
        float pairInvMass = pi0lv.M();

        float weight = 1.0;
	      
        if (ketapi == 0) {
          weight = pi0pt->Eval(ptPi);
          pipt->Fill(ptPi, weight);
        } else {
          weight = etapt->Eval(ptPi);
          hetapt->Fill(ptPi, weight);
        }
		
        float alph = fabs(pho1.E() - pho2.E()) / (pho1.E() + pho2.E());
        if (((pho1.E() > 1.0 && pho2.E() > 0.60) ||
             (pho2.E() > 1.0 && pho1.E() > 0.60)) &&
            fabs(pi0lv.Pt()) > 2.00 && alph < 0.6) {
          invMass->Fill(pairInvMass, pi0lv.Pt(), weight);
        }

        if (pho2.E() > 0.60 && pho1.E() > 0.6 && alph < 1.0) {
          ptp_pt1_dang->Fill(pi0lv.Pt(), pho1.Pt(), pho1.DeltaR(pho2));
        }
		    
        gampt->Fill(g1->Pt(), weight);
        gampt->Fill(g2->Pt(), weight);
        gamcorr->Fill(g1->Pt(), weight);
        gamcorr->Fill(g2->Pt(), weight);
		
        thetaDist->Fill(g1->Vect().Angle(g2->Vect()));
        if (ketapi == 0) {
          gampipt->Fill(g1->Pt(), weight);
          gampipt->Fill(g2->Pt(), weight);
        } else {
          gametapt->Fill(g1->Pt(), weight);
          gametapt->Fill(g2->Pt(), weight);
        }
      } // if (ePi < 100)
    } // loop over ketapi (pi0/eta)
  } // loop over numPions
	  
  // Write out histograms and functions.
  invMass->Write();
  ptp_pt1_dang->Write();
  pi0pt->Write();
  etapt->Write();

  // Overwrite picorr bins using the pi0pt function.
  for (int k = 1; k <= NumBins; k++) {
    float binCent = picorr->GetBinCenter(k);
    float ptval = pi0pt->Eval(binCent);
    picorr->SetBinContent(k, ptval);
  }
  
  picorr->Scale(1e-5);
  gamcorr->Scale(1e-15);

  ratiogampipt->Divide(gampipt, pipt);
  ratiopt->Divide(gampt, pipt, 1.02, 1.0);
	  
  out->Write(0, TObject::kOverwrite);
  std::cout << "Done" << std::endl;
  
  // Close the file. Remove the following exit call if you want to remain in an interactive session.
  out->Close();
  // gSystem->Exit(0);
}

