#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TFile.h"



double Eta2Theta(double etaVal) {return 2.0 * atan(exp(-1.0* etaVal));}
/*
  double getTransFunVal(double ptVal)
  {
  if (ptVal < 2.0)
  return 0.0;
  
  if (ptVal >= 2.0 && ptVal < 4.0)
  return 0.5 * ptVal - 1.0;

  if (ptVal >=4.0)
  return 1.0;
  }

  double getPi0FunVal(TF1 * fit1, TF1 * fit2, float ptval)
  {

  float val1 = fit1->Eval(ptval);
  float val2 = fit2->Eval(ptval);
  
  return val1 * (1.0 - getTransFunVal(ptval)) +
  va12 * getTransFunVal();

  }
*/

void super_fastmc_theta(){
	
  const char * outputName = "test.root";
  int  numPions = 1e6;

  //now takes input pi0 file, with fits already done, uses those fits
  //  gSystem->Load("libdirPhotons.so");
  //gSystem->Load("libPhysics.so");

  float maxPt = 20.0;
  int NumBins = 500;

  /*
    strstream cpCom;
    cpCom << "cp -f " << inputFile << "  " << outputName << ends;
    cout << "updating file with command: " << cpCom.str() << endl;
    gSystem->Exec(cpCom.str());
    TFile * out = new TFile(outputName, "UPDATE");
  */

  TFile * out = new TFile(outputName, "UPDATE");


  // PIMASS^2 = 0.018225
  // ETAMASS^2 = 0.299209
  // at high pt, pi0/eta = 0.55

  TF1 * flatPt = new TF1("flatpt" , "1", 0.02, maxPt);
  
  float MYPI = 3.14159;
  // flat in rapidity means ~flat in eta means d(eta)/d(theta) in theta
  //  TF1 *f_pi0theta = new TF1("pi0Theta", "0.5*tan(0.5*x)*(1/cos(0.5*x))*(1/cos(0.5*x))",6*MYPI/18,12*MYPI/18); 
  
  TF1 *flat_rapidity = new TF1("flatRap", "1", -1, 1);
  TF1 *flat_phi = new TF1("flatPhi", "1", 0, 2*MYPI);
  //  TF1 *flat_theta = new TF1("flatTheta","1",0,MYPI);
  TF1 *real_theta = new TF1("flatTheta","sin(x)",0,MYPI);
  
  //    TF1 powlaw1("powlaw1","[0] / (1.72 + x)^[1]", 0, 6.0); 
  //    TF1 powlaw2("powlaw2","[0] / (x)^[1]", 4.0, 15.0); 

  //    TF1 etapowlaw1("powlaw1","[0] / (1.72 + sqrt(x^2 + 0.299209 - 0.018225) )^[1]", 0, 6.0); 
  //    TF1 etapowlaw2("powlaw2","[0] / (sqrt(x^2 + 0.299209 - 0.018225))^[1]", 4.0, 15.0); 

  TF1 woodsSax("woodsSax","1.0- 1.0/(1.0 + exp((0.5*x-0.25)/0.02))", 0.0, maxPt);
  woodsSax.SetNpx(20000);

  TH1F * picorr = new TH1F("picorr",
			   "corrected pi0 spect, according to fit values", 
			   NumBins,0.0,maxPt);
  TH1F * gamcorr = new TH1F("gamcorr",
			    "correct gamma spectrum according to pi0 fit values",
			    NumBins,0.0,maxPt);
  
  TH2F * invMass = new TH2F("invMassPt","",160,0.0,0.8,4,2.0,6.0);
  invMass->Sumw2();
  
  TH1F * pipt = new TH1F(*picorr); // gets same binning
  pipt->Reset();
  

  TH3F * ptp_pt1_dang = new TH3F("ptp_pt1_dang","",50,1.0,20.0,50,1.0,20.0,600,0,0.600);
  
  TH1F * ratiopt = new TH1F(*pipt);
  ratiopt->SetName("ratiopt");
  //ratiopt->SetName(mcgpiTa->ptHist[jj][icut][idet]->GetName());
  //ratiopt->SetDirectory(centDir);
  cout << "processing " << ratiopt->GetName() << endl; 
  //delete mcgpiTa->ptHist[jj][icut][idet];
  //mcgpiTa->ptHist[jj][icut][idet] = ratiopt;
    
  
  TH1F * gampt = new TH1F(*pipt);
  TH1F * gampipt = new TH1F(*pipt);
  TH1F * ratiogampipt = new TH1F(*pipt);

  //eta stuff
  TH1F * hetapt = new TH1F(*pipt);
  TH1F * gametapt = new TH1F(*pipt);	  
  TH1F * ratiogametapt = new TH1F(*pipt);


  pipt->SetName("pipt");
  gampt->SetName("gampt");
  //ratiopt->SetName("ratiopt");
  gampipt->SetName("gampipt");
  ratiogampipt->SetName("ratiogampipt");
  //eta stuff
  /*
    hetapt->SetName("hetapt");
    gametapt->SetName("gametapt");
    ratiogametapt->SetName("ratiogametapt");
  */
	  
  pipt->SetTitle("Pi0 Pt Dist");
  gampt->SetTitle("All Gamma Pt Dist");
  ratiopt->SetTitle("All Gamma / Pi0 Ratio (Pt)");
  gampipt->SetTitle("Pi0 Gam. Pt Dist");
  ratiogampipt->SetTitle("Gamma_Pi0 / Pi0 Ratio (Pt)");
	
  /*
  //eta stuff
  hetapt->SetTitle("Eta Pt Dist"); 
  gametapt->SetTitle("Eta Gam. Pt Dist");
  ratiogametapt->SetTitle("Gamma_Eta / Pi0 Ratio (Pt)");
  */
	  
  TH1F * thetaDist = new TH1F("thetaDist", "distribution of opening angles of the two gamma", 3140, 0,3.14);
	  		
  /* strstream fit1name, fit2name; */
  /* fit1name << "fit1" << ptH->GetName() << ends; */
  /* fit2name << "fit2" << ptH->GetName() << ends; */
	  
  /* TF1 * fit1 = (TF1 *) centDir->Get(fit1name.str()); */
  /* TF1 * fit2 = (TF1 *) centDir->Get(fit2name.str()); */

  /*
    TF1 * etapt = new TF1("etapt","10^15/x^[0]", 0.1, maxPt);
    TF1 * pi0pt = new TF1("pi0pt","10^15/x^[0]", 0.1, maxPt);

    pi0pt->SetParameter(0,6.5);
    etapt->SetParameter(0,6.5);
  */

  TF1 * pi0pt = new TF1("pi0pt","[0]*(x^0.5)/([1]+x)^[2]",0.001,10);
  //pi0pt->SetParameters(10^8,1.5,11.0);
  //
  pi0pt->SetParameters(1.0e11,0.5,8.0);


	  
  TF1 * etapt = new TF1("etapt","[0]*x/([1]+x)^[2]",0.001,10);
  etapt->SetParameters(10^8,1.1,5.0);
	  

  //	  TF1 * pi0pt = new TF1("pi0pt","(1.0-1.0/(1.0 + exp((0.5*x-0.2)/0.02)))*[0]*x/([1]+x)^[2]",0.001,10);



  /*

    TF1 * etapt = new TF1("etapt","10^15*([0]*x/((sqrt(x^2 + 0.299209 - 0.018225)+[1])^[2])) * woodsSax  + 10^15 * (1.0 - woodsSax) * ([3] * x / (sqrt(x^2 + 0.299209 - 0.018225))^[4])" ,0.0, maxPt);
    TF1 * pi0pt = new TF1("pi0pt","10^15*([0]*x/((x+[1])^[2])) * woodsSax  + (1.0 - woodsSax) * 10^15 * ([3]*x / x^[4])" ,0.0, maxPt);
	       
	  
    pi0pt->SetParameter(0, fit1->GetParameter(0));
    pi0pt->SetParameter(1, 1.72);
    pi0pt->SetParameter(2, fit1->GetParameter(1));
    pi0pt->SetParameter(3, fit2->GetParameter(0));
    pi0pt->SetParameter(4, fit2->GetParameter(1));
	  
    etapt->SetParameter(0, fit1->GetParameter(0));
    etapt->SetParameter(1, 1.72);
    etapt->SetParameter(2, fit1->GetParameter(1));
    etapt->SetParameter(3, fit2->GetParameter(0));
    etapt->SetParameter(4, fit2->GetParameter(1));
  */


  TF1 * f1 = new TF1("f1reso","sqrt([0]*[0]/x+[1]*[1])",0.5,20.);
  f1->SetParameter(0, 0.165);
  f1->SetParameter(1, 0.26);
	  
  TF1 * fgaus = new TF1("mygaus","gaus", -1.5,1.5);
  fgaus->SetParameter(0,10e5);
  fgaus->SetParameter(1,0); //mean

  TH2F * h2esm = new TH2F("h2esm","",80,0,20,200, 0, 20);
  TH2F * h2resoe = new TH2F("h2resoe","",80, 0, 20,200,-1,1);

  // decay stuff
  for (int i = 0; i < numPions; i++ ) { 
    for (int ketapi = 0; ketapi < 2; ketapi++) {
	      
      //	      if (((jj > 2 && jj < 6) || jj == 9 || jj == 8) && numPions > 3000)
      //	continue;
		
      //	      if (numPions % 1000 == 0) 
      //	cout << numPions << endl;

      float massPi = 0.13497;
      if (ketapi == 1) 
	{
	  massPi = 0.547;
	  //eta's should only decay 39 % of the time to 2 gammas
	  //combined with the eta/pi == .55, we should only decay
	  // eta's 22/100 times 
	  //		  if ((i % 100) < 78) continue;
	}
	      
      //if (i % 100 == 0) cout << i << " pions processed" << endl;
      float thetaOrig =  Eta2Theta(flat_rapidity->GetRandom());
      //f_theta->GetRandom();      
	      
      float ptPi = flatPt->GetRandom(); 
      float pPi = ptPi / sin(thetaOrig);
	      
      //  	      float pPi = flatPt->GetRandom(); 
      //  	      float ptPi = pPi * sin(thetaOrig);


      float ePi = sqrt(pPi * pPi + massPi * massPi);

      if (ePi < 100) {
		
	TVector3 pion3v(2.0,0.0,0.0);
	pion3v.SetTheta(thetaOrig);	  
	pion3v.SetPhi(flat_phi->GetRandom());
	pion3v.SetMag(pPi);
		
	TLorentzVector pcylorig;
	pcylorig.SetVect(pion3v);
	pcylorig.SetE(ePi);	  
		
	// make photons in CM and transform them back
	TLorentzVector photon1;
	TLorentzVector photon2;
		
	TLorentzVector * g1 = &photon1;
	TLorentzVector * g2 = &photon2;

	TVector3 photon3v(2.0,0.0,0.0);
	photon3v.SetMag(massPi/2);
	photon3v.SetTheta( real_theta->GetRandom() );
	photon3v.SetPhi(flat_phi->GetRandom() );
	g1->SetVect(1 * photon3v);
	g1->SetE(massPi/2);
		
	g2->SetVect( -1 * photon3v);
	g2->SetE(massPi/2);
	g1->Boost(pcylorig.BoostVector());
	g2->Boost(pcylorig.BoostVector());
	// cout << "E was " << g1->E() << " and " << g2->E()<< endl;
	// check that backwards boost yeilds at rest pion.
	/*    
	      TLorentzVector test = pcylorig;
	      test.Boost(-1.0 * pcylorig.BoostVector());
	      cout << "X: (Px)" << test.X() << "(" << test.Px() << ")" << endl;
	      cout << "Y,Z,T " << test.Y()<< endl << test.Z() << endl << test.T() endl;
	*/

		
	// reco'd variables
	TLorentzVector pho1, pho2, pi0lv;

	float tt_clus_energy = g1->E();
	float reso1 =  f1->Eval(tt_clus_energy);
	fgaus->SetParameter(2,reso1);
	float smear1 = fgaus->GetRandom();
	tt_clus_energy = tt_clus_energy + tt_clus_energy*smear1;
	h2esm->Fill(g1->E(), tt_clus_energy);
	h2resoe->Fill(g1->E(),tt_clus_energy/g1->E()-1.0);
	float tt_clus_pt = g1->Pt();
	tt_clus_pt = tt_clus_pt + tt_clus_pt*smear1;
	float tt_clus_phi = g1->Phi();
	float tt_clus_eta = g1->Eta();


	float tt_clus_energy2 = g2->E();
	float reso2 =  f1->Eval(tt_clus_energy2);
	fgaus->SetParameter(2,reso2);
	float smear2 = fgaus->GetRandom();
	tt_clus_energy2 = tt_clus_energy2 + tt_clus_energy2*smear2;
	h2esm->Fill(g2->E(), tt_clus_energy2);
	h2resoe->Fill(g2->E(),tt_clus_energy2/g2->E()-1.0);
	float tt_clus_pt2 = g2->Pt();
	tt_clus_pt2 = tt_clus_pt2 + tt_clus_pt2*smear2;
	float tt_clus_phi2 = g2->Phi();
	float tt_clus_eta2 = g2->Eta();

	pho1.SetPtEtaPhiE(tt_clus_pt, tt_clus_eta, tt_clus_phi, tt_clus_energy);
	pho2.SetPtEtaPhiE(tt_clus_pt2, tt_clus_eta2, tt_clus_phi2, tt_clus_energy2);
		
	/*
	  if (pho1.DeltaR(pho2) > 0.8) continue;
	  if (pho1.Eta()/pho2.Eta() < 0) continue;
	*/
		
	pi0lv = pho1 + pho2;
	float pairInvMass = pi0lv.M();

	float weight = 1.0;
		
	if (ketapi == 0)
	  {
	    weight = pi0pt->Eval(ptPi);
	    //weight = 1.0;
	    pipt->Fill(ptPi, weight);
	  }
	else
	  {
	    weight = etapt->Eval(ptPi);
	    hetapt->Fill(ptPi,weight);
	  }
		
	float alph = fabs(pho1.E() - pho2.E())/(pho1.E() + pho2.E());
	if (  ((pho1.E()  > 1.0 && pho2.E() > 0.60) ||
	       (pho2.E()  > 1.0 && pho1.E() > 0.60) ) &&
	      fabs(pi0lv.Pt()) > 2.00 &&
	      alph < 0.6)
	  {
	    invMass->Fill(pairInvMass, pi0lv.Pt(), weight);
	  }

	if (  pho2.E() > 0.60 && pho1.E() > 0.6 && 
	      alph < 1.0)
	  {
	    //		    ptp_pt1_dang->Fill(pi0lv.Pt(),pho1.Pt(),pho1.Angle(pho2.Vect()));
	    //		    ptp_pt1_dang->Fill(pi0lv.E(),pho1.E(),pho1.Angle(pho2.Vect()));
	    ptp_pt1_dang->Fill(pi0lv.Pt(),pho1.Pt(),pho1.DeltaR(pho2));
	    //ptp_pt1_dang->Fill(ePi,g1->E(),g1->Angle(g2->Vect()));
	    //		    ptp_pt1_dang->Fill(ptPi,g1->Pt(),g1->DeltaR(*g2));
	  }

		    
	gampt->Fill(g1->Pt(), weight);
	gampt->Fill(g2->Pt(), weight);
	gamcorr->Fill(g1->Pt(), weight);
	gamcorr->Fill(g2->Pt(), weight);
		
	thetaDist->Fill(g1->Vect().Angle(g2->Vect()));
	if (ketapi == 0)
	  {
	    gampipt->Fill(g1->Pt(), weight);
	    gampipt->Fill(g2->Pt(), weight);
	  }
	else 
	  {
	    gametapt->Fill(g1->Pt(), weight);
	    gametapt->Fill(g2->Pt(), weight);
	  }
		
      }// if ePi < 100
    } // ketapi (eta or pi)
  }// for numPions
	  
	  
	  
  invMass->Write();
  ptp_pt1_dang->Write();
  /* invMass2->Write(); */
  /* invMass3->Write(); */
  /* invMass4->Write(); */
	  
  /* fitResults->Write(); */

  // fit gampt
  pi0pt->Write();
  etapt->Write();

  // fix me....fit the gam fitting
  // so that this can be used as an input into
  // the eff calculation
	  
  for ( int k = 1; k < NumBins+1; k++)
    {
      float binCent = picorr->GetBinCenter(k);
      float ptval = pi0pt->Eval(binCent);
      float gptval =  0.0; 
      picorr->SetBinContent(k,ptval);
	      
    }
  //  gampt->GetFunction("pi0pt")->Eval(binCent);
  //pi0pt->SetParLimits(1, 1.72,1.72);	  T
  //  	  gampt->Fit("pi0pt");
	    

  picorr->Scale(1e-5);
  gamcorr->Scale(1e-15);

  ratiogampipt->Divide(gampipt, pipt);
  //	  ratiogametapt->Divide(gametapt,pipt);
	  
  ratiopt->Divide(gampt,pipt, 1.02, 1.0);
  // scale by up 2% to account for other decays
	  
  //    TCanvas * c1 = new TCanvas();
  //    ratiopt->Draw();
  //    ratiogampipt->SetLineColor(2);
  //    ratiogampipt->Draw("same");
  //    ratiogametapt->SetLineColor(7);
  //    ratiogametapt->Draw("same");
  
  //    TCanvas * c2 = new TCanvas();
  //    picorr->Draw();
  
  //    TCanvas * c3 = new TCanvas();
  //    gamcorr->Draw();

  
  out->Write(0, TObject::kOverwrite);

	std::cout << "Done" << std::endl;
	gSystem->Exit(0);  
  //  out->Close();

  //  exit(0);
}



