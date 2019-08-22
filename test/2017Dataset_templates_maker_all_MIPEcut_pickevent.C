//7/23, 2019 
//This is the SinglePhoton analysis template histograms maker for 2017 whole dataset (B,C,D,E,F)

#define temp_cxx
#include "temp.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include <fstream>

using namespace std;

/*
1. out-of-time photons only
2. in-time photons only
3. out-of-time photons and in-time photons but we only need out-of-time photons


1.)  If there is no OUT of time photon, and there is an IN time photon that can satisfy the trigger, take the IN time photon that satisfies the trigger.
2.)  If there is no IN time photon, and there is an OUT of time photon that can satisfy the trigger, then take the OUT of time photon that satisfies the trigger.
3.)  If there is BOTH an IN time photon AND an OUT of time photon, we must select the proper one to use:
a.)  If the IN time photon can satisfy the trigger and the OUT of time photon cannot, take the IN time photon.
b.)  If there is an OUT of time photon that can satisfy the trigger and there is not an IN time photon that can satisfy the trigger, take the OUT of time photon.
c.)  If there are BOTH an IN time photon and an OUT of time photon that can satisfy the trigger, take the OUT of time photon.

That handles the photon disambiguation.

The MET will require the following:
1.)  If you selected an IN time photon, you don't need to do anything.
2.)  If you selected an OUT of time photon:
  a.)  If there is NO IN time photon, then you need to correct the MET for the ET of the out of time photon.
  b.)  If there is both an IN time photon and an OUT of time photon, then you need to remove the IN time photon from the MET and then correct the MET for the OUT of time photon.

*/


Float_t MIPChi2_over_MIPNhitCone(Float_t MIPChi2, Float_t MIPNhitCone)
{
	Float_t phoMIPChi2 = MIPChi2;
	Float_t phoMIPNhitCone = MIPNhitCone;
	Float_t MIPChi2_over_MIPNhitCone = MIPChi2 / MIPNhitCone;

	return MIPChi2_over_MIPNhitCone;
}

Float_t newmMET(Float_t met, Float_t metphi, Float_t phi, Float_t et)
{
	Float_t phoEt;
	Float_t phoPhi;
	Float_t pfMET = met;
	Float_t pfMETPhi = metphi;
	Float_t ophoEt = et;
	Float_t ophoPhi = phi;

	Float_t pfMETX = met * cos(pfMETPhi);
	Float_t pfMETY = met * sin(pfMETPhi);
	Float_t ITEtX = phoEt * cos(phoPhi);
	Float_t ITEtY = phoEt * sin(phoPhi);

	pfMETX += ITEtX;
	pfMETY += ITEtY;

	Float_t OOTEtX = et * cos(ophoPhi);
	Float_t OOTEtY = et * sin(ophoPhi);

	pfMETX -= OOTEtX;
	pfMETY -= OOTEtY;

	Float_t newmMET = sqrt(pfMETX * pfMETX + pfMETY * pfMETY);

	return newmMET;
}

Float_t newuMET(Float_t met, Float_t metphi, Float_t phi, Float_t et)
{
	Float_t ophoEt = et;
	Float_t ophoPhi = phi;
	Float_t pfMET = met;
	Float_t pfMETPhi = metphi;

	Float_t pfMETX = met * cos(pfMETPhi);
	Float_t pfMETY = met * sin(pfMETPhi);
	Float_t OOTEtX = et * cos(ophoPhi);
	Float_t OOTEtY = et * sin(ophoPhi);

	pfMETX -= OOTEtX;
	pfMETY -= OOTEtY;

	Float_t newuMET = sqrt(pfMETX * pfMETX + pfMETY * pfMETY);

	return newuMET;
}



//calulate invariant mass Z

Float_t InvariMass(Float_t Et1, Float_t Et2, Float_t Phi1, Float_t Phi2, Float_t Eta1, Float_t Eta2)
{
	Float_t Theta1 = 2 * atan(exp(-1.*Eta1));
	Float_t Theta2 = 2 * atan(exp(-1.*Eta2));
	Float_t phoPhi1 = Phi1;
	Float_t phoPhi2 = Phi2;
	Float_t Etot1 = Et1 / sin(Theta1); //E_tot1
	Float_t Etot2 = Et2 / sin(Theta2); //E_tot2
	
	//reconstruct the vectors for x, y and z

	
	Float_t phoX1 = Etot1 * cos(Phi1) * sin(Theta1); 
	Float_t phoY1 = Etot1 * sin(Phi1) * sin(Theta1);
	Float_t phoZ1 = Etot1 * cos(Theta1);

	Float_t phoX2 = Etot2 * cos(Phi2) * sin(Theta2);
	Float_t phoY2 = Etot2 * sin(Phi2) * sin(Theta2);
	Float_t phoZ2 = Etot2 * cos(Theta2);

	Float_t EX1 = Et1*sin(Phi1);
	Float_t EY1 = Et1*cos(Phi1);
	Float_t EZ1 = Etot1*cos(Theta1);

	Float_t EX2 = Et2*sin(Phi2);
	Float_t EY2 = Et2*cos(Phi2);
	Float_t EZ2 = Etot2*cos(Theta2);

	Float_t E1 = sqrt(Etot1*Etot1);
	Float_t E2 = sqrt(Etot2*Etot2);

	//Float_t InvariMassSq = 2 * E1*E2 - phoMag1 * phoMag1 - phoMag2 * phoMag2 - 2 * DotProd12;

	//1	
	Float_t InvariMassSq = (Etot1 + Etot2)*(Etot1 + Etot2) - (EX1 + EX2)*(EX1 + EX2) - (EY1 + EY2)*(EY1 + EY2) - (EZ1 + EZ2)*(EZ1 + EZ2);

	Float_t InvMass = sqrt(InvariMassSq);

	
	
	return InvMass;
	
}


//calulate transverse mass W (we cannot calcuate the invariant mass of W, since we don't know the mass of neutrino)

Float_t TransMassW(Float_t Et1, Float_t Phi1, Float_t Eta1)
{
	Float_t Theta1 = 2 * atan(exp(-1.*Eta1));
	Float_t phoPhi1 = Phi1;
	Float_t Etot1 = Et1 / sin(Theta1); //E_tot1

	//reconstruct the vectors for x, y and z

	Float_t phoX1 = Etot1 * cos(Phi1) * sin(Theta1);
	Float_t phoY1 = Etot1 * sin(Phi1) * sin(Theta1);
	Float_t phoZ1 = Etot1 * cos(Theta1);

	Float_t EX1 = Et1 * sin(Phi1);
	Float_t EY1 = Et1 * cos(Phi1);
	Float_t EZ1 = Etot1 * cos(Theta1);

	Float_t TransMassSqW = (Etot1)*(Etot1) - (EZ1)*(EZ1);

	Float_t TranMassW = sqrt(TransMassSqW);



	return TranMassW;

}


//for checking matching

Float_t DeltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2)
{

	Float_t ophoPhi = phi1;
	Float_t phoPhi = phi2;
	Float_t dphi = fabs(phoPhi - ophoPhi);
	Float_t tdphi = dphi;
	if (dphi > TMath::Pi()) tdphi = TMath::Pi()*2.0 - dphi;
	dphi = tdphi;

	Float_t ophoEta = eta1;
	Float_t phoEta = eta2;
	Float_t deta = fabs(phoEta - ophoEta);

	Float_t deltaR = sqrt(deta*deta + dphi * dphi);
	return deltaR;
}


void temp::Loop()
{
	//   In a ROOT session, you can do:
	//      root> .L temp.C
	//      root> temp t
	//      root> t.GetEntry(12); // Fill t data members with entry number 12
	//      root> t.Show();       // Show values of entry 12
	//      root> t.Show(16);     // Read and show values of entry 16
	//      root> t.Loop();       // Loop on all entries
	//

	//     This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch
	

	TFile *monophoseedtime = new TFile("2017Dataset_SinglePhoton_analysis_all_MIPEcut_pickevent_metfilterplot.root", "RECREATE");

	//------------timing distribution for IT and OOT-------------
	
	TH1F *ITST = new TH1F("ITST", "ITSeedTime", 200, -25, 25);
	ITST->GetXaxis()->SetTitle("ITSeedTime (ns)");
	ITST->GetYaxis()->SetTitle("Entries");

	TH1F *OOTST = new TH1F("OOTST", "OOTSeedTime", 200, -25, 25);
	OOTST->GetXaxis()->SetTitle("OOTSeedTime (ns)");
	OOTST->GetYaxis()->SetTitle("Entries");

	TH1F *MITST = new TH1F("MITST", "ITSeedTime with MET cut", 200, -25, 25);
	MITST->GetXaxis()->SetTitle("ITSeedTime (ns)");
	MITST->GetYaxis()->SetTitle("Entries");

	TH1F *MOOTST = new TH1F("MOOTST", "OOTSeedTime with MET cut", 200, -25, 25);
	MOOTST->GetXaxis()->SetTitle("OOTSeedTime (ns)");
	MOOTST->GetYaxis()->SetTitle("Entries");

	//--------------------------------------------------------------


	//-----------------full timing distribution---------------------

	TH1F *FullTime = new TH1F("FullTime", "Timing Distribution", 200, -25, 25);
	FullTime->GetXaxis()->SetTitle("SeedTime (ns)");
	FullTime->GetYaxis()->SetTitle("Entries");

	//MET cut
	TH1F *MFullTime = new TH1F("MFullTime", "Timing Distribution_MET", 200, -25, 25);
	MFullTime->GetXaxis()->SetTitle("SeedTime (ns)");
	MFullTime->GetYaxis()->SetTitle("Entries");

	//--------------------------------------------------------------


	//----------------------Templates-----------------------

	//Candidate events
	TH1F *Candidate_Events = new TH1F("Candidate_Events", "Candidate Events", 200, -25, 25);
	Candidate_Events->GetXaxis()->SetTitle("SeedTime (ns)");
	Candidate_Events->GetYaxis()->SetTitle("Entries / 0.25 ns");

	//W Prompt Template
	TH1F *promptWTemp = new TH1F("promptWTemp", "prompt Template (W)", 200, -25, 25);
	promptWTemp->GetXaxis()->SetTitle("SeedTime (ns)");
	promptWTemp->GetYaxis()->SetTitle("Entries");

	//MIP total Energy from W Prompt
	TH1F *MIPW = new TH1F("MIPW", "MIP totE from W", 2000, 20, 20);
	MIPW->GetXaxis()->SetTitle("MIP totE (GeV)");
	MIPW->GetYaxis()->SetTitle("Entries");
	
	//Transverse Mass of W
	TH1F *TransMass_W = new TH1F("TransMass_W", "Transverse mass W", 2000, 20, 20);
	TransMass_W->GetXaxis()->SetTitle("Transverse mass (GeV)");
	TransMass_W->GetYaxis()->SetTitle("Entries");
	
	//Z Prompt Template
	TH1F *promptZTemp = new TH1F("promptZTemp", "prompt Template (Z)", 200, -25, 25);
	promptZTemp->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZTemp->GetYaxis()->SetTitle("Entries");

	//MIP total Energy from Z Prompt
	TH1F *MIPZ = new TH1F("MIPZ", "MIP totE from Z", 2000, 20, 20);
	MIPZ->GetXaxis()->SetTitle("MIP totE (GeV)");
	MIPZ->GetYaxis()->SetTitle("Entries");

	//Invariant Mass of Z
	TH1F *InvMass_Z = new TH1F("InvMass_Z", "invariant mass Z", 2000, 20, 20);
	InvMass_Z->GetXaxis()->SetTitle("Invariant mass (GeV)");
	InvMass_Z->GetYaxis()->SetTitle("Entries");

	//Spike Template
	TH1F *SpikeTemp = new TH1F("SpikeTemp", "Spike Template", 200, -25, 25);
	SpikeTemp->GetXaxis()->SetTitle("SeedTime (ns)");
	SpikeTemp->GetYaxis()->SetTitle("Entries");

	//Beam Halo Template
	TH1F *BeamHaloTemp = new TH1F("BeamHaloTemp", "Beam Halo Template", 200, -25, 25);
	BeamHaloTemp->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHaloTemp->GetYaxis()->SetTitle("Entries / 0.25 ns");

	//Check shower shape vs. time
	TH2F *BeamHaloSieie_vs_time = new TH2F("BeamHaloSieie_vs_time", "SigmaIEtaIEta vs Time", 200, 20, 20, 200, 20.0, 20.0);
	BeamHaloSieie_vs_time->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHaloSieie_vs_time->GetYaxis()->SetTitle("SigmaIEtaIEta");

	//pick event
	TH2F *BeamHalo_prompt_region = new TH2F("BeamHalo_prompt_region", "SigmaIEtaIEta vs Time prompt region", 200, 20, 20, 200, 20.0, 20.0);
	BeamHalo_prompt_region->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHalo_prompt_region->GetYaxis()->SetTitle("SigmaIEtaIEta");

	//Check shower shape
	TH1F *BeamHaloSieie = new TH1F("BeamHaloSieie", "SigmaIEtaIEta from Beam Halo", 2000, 20, 20);
	BeamHaloSieie->GetXaxis()->SetTitle("SigmaIEtaIEta");
	BeamHaloSieie->GetYaxis()->SetTitle("Entries");

	//METfilter for beamhalo
	TH1F *metfilter_beamhalo = new TH1F("metfilter_beamhalo", "MET Filter Beam Halo", 200, 20, 20);
	metfilter_beamhalo->GetXaxis()->SetTitle("MET Filter");
	metfilter_beamhalo->GetYaxis()->SetTitle("Entries");

	//METfilter for candidate event
	TH1F *metfilter_candidate = new TH1F("metfilter_candidate", "MET Filter Candidate", 200, 20, 20);
	metfilter_candidate->GetXaxis()->SetTitle("MET Filter");
	metfilter_candidate->GetYaxis()->SetTitle("Entries");
	

	//--------------------------------------------------------

	//----------------------All MIP variables for W prompt with MIP totE < 4.9 GeV-------------------

	TH1F *MIPTotalE_Wprompt = new TH1F("MIPTotalE_Wprompt", "MIP totE from W prompt", 2000, 20, 20);
	MIPTotalE_Wprompt->GetXaxis()->SetTitle("MIP totE (GeV)");
	MIPTotalE_Wprompt->GetYaxis()->SetTitle("Entries");

	TH1F *MIPChi2_Wprompt = new TH1F("MIPChi2_Wprompt", "MIP Chi2 from W prompt", 20000, 20, 20);
	MIPChi2_Wprompt->GetXaxis()->SetTitle("MIP Chi2");
	MIPChi2_Wprompt->GetYaxis()->SetTitle("Entries");

	TH1F *MIPSlope_Wprompt = new TH1F("MIPSlope_Wprompt", "MIP Slope from W prompt", 2000, 20, 20);
	MIPSlope_Wprompt->GetXaxis()->SetTitle("MIP Slope");
	MIPSlope_Wprompt->GetYaxis()->SetTitle("Entries");

	TH1F *MIPIntercept_Wprompt = new TH1F("MIPIntercept_Wprompt", "MIP Intercept from W prompt", 2000, 20, 20);
	MIPIntercept_Wprompt->GetXaxis()->SetTitle("MIP Intercept");
	MIPIntercept_Wprompt->GetYaxis()->SetTitle("Entries");

	TH1F *MIPNhitCone_Wprompt = new TH1F("MIPNhitCone_Wprompt", "MIP NhitCone from W prompt", 18000, 0, 180);
	MIPNhitCone_Wprompt->GetXaxis()->SetTitle("MIP NhitCone");
	MIPNhitCone_Wprompt->GetYaxis()->SetTitle("Entries");

	TH1F *MIPisHalo_Wprompt = new TH1F("MIPisHalo_Wprompt", "MIP isHalo from W prompt", 2000, 20, 20);
	MIPisHalo_Wprompt->GetXaxis()->SetTitle("MIP isHalo");
	MIPisHalo_Wprompt->GetYaxis()->SetTitle("Entries");

	TH1F *MIPCHi2_over_MIPNhitCone_Wprompt = new TH1F("MIPCHi2_over_MIPNhitCone_Wprompt", "MIPChi2/MIPNhitCone from Wprompt", 2000, 20, 20);
	MIPCHi2_over_MIPNhitCone_Wprompt->GetXaxis()->SetTitle("MIPChi2/MIPNhitCone");
	MIPCHi2_over_MIPNhitCone_Wprompt->GetYaxis()->SetTitle("Entries");

	//------------------------------------------------------------------------------

	//----------------------All MIP variables for Z prompt with MIP totE < 4.9 GeV--------------------

	TH1F *MIPTotalE_Zprompt = new TH1F("MIPTotalE_Zprompt", "MIP totE from Z prompt", 2000, 20, 20);
	MIPTotalE_Zprompt->GetXaxis()->SetTitle("MIP totE (GeV)");
	MIPTotalE_Zprompt->GetYaxis()->SetTitle("Entries");

	TH1F *MIPChi2_Zprompt = new TH1F("MIPChi2_Zprompt", "MIP Chi2 from Z prompt", 20000, 20, 20);
	MIPChi2_Zprompt->GetXaxis()->SetTitle("MIP Chi2");
	MIPChi2_Zprompt->GetYaxis()->SetTitle("Entries");

	TH1F *MIPSlope_Zprompt = new TH1F("MIPSlope_Zprompt", "MIP Slope from Z prompt", 2000, 20, 20);
	MIPSlope_Zprompt->GetXaxis()->SetTitle("MIP Slope");
	MIPSlope_Zprompt->GetYaxis()->SetTitle("Entries");

	TH1F *MIPIntercept_Zprompt = new TH1F("MIPIntercept_Zprompt", "MIP Intercept from Z prompt", 2000, 20, 20);
	MIPIntercept_Zprompt->GetXaxis()->SetTitle("MIP Intercept");
	MIPIntercept_Zprompt->GetYaxis()->SetTitle("Entries");

	TH1F *MIPNhitCone_Zprompt = new TH1F("MIPNhitCone_Zprompt", "MIP NhitCone from Z prompt", 9000, 0, 90);
	MIPNhitCone_Zprompt->GetXaxis()->SetTitle("MIP NhitCone");
	MIPNhitCone_Zprompt->GetYaxis()->SetTitle("Entries");

	TH1F *MIPisHalo_Zprompt = new TH1F("MIPisHalo_Zprompt", "MIP isHalo from Z prompt", 2000, 20, 20);
	MIPisHalo_Zprompt->GetXaxis()->SetTitle("MIP isHalo");
	MIPisHalo_Zprompt->GetYaxis()->SetTitle("Entries");

	//------------------------------------------------------------------------------


	//----------------------All MIP variables for Beam Halo WITHOUT MIP totE > 4.9 GeV--------------------

	TH1F *woMIPTotalE_BeamHalo = new TH1F("woMIPTotalE_BeamHalo", "MIP totE from Beam Halo", 2000, 20, 20);
	woMIPTotalE_BeamHalo->GetXaxis()->SetTitle("MIP totE (GeV)");
	woMIPTotalE_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *woMIPChi2_BeamHalo = new TH1F("woMIPChi2_BeamHalo", "MIP Chi2 from Beam Halo", 20000, 20, 20);
	woMIPChi2_BeamHalo->GetXaxis()->SetTitle("MIP Chi2");
	woMIPChi2_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *woMIPSlope_BeamHalo = new TH1F("woMIPSlope_BeamHalo", "MIP Slope from Beam Halo", 2000, 20, 20);
	woMIPSlope_BeamHalo->GetXaxis()->SetTitle("MIP Slope");
	woMIPSlope_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *woMIPIntercept_BeamHalo = new TH1F("woMIPIntercept_BeamHalo", "MIP Intercept from Beam Halo", 2000, 20, 20);
	woMIPIntercept_BeamHalo->GetXaxis()->SetTitle("MIP Intercept");
	woMIPIntercept_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *woMIPNhitCone_BeamHalo = new TH1F("woMIPNhitCone_BeamHalo", "MIP NhitCone from Beam Halo", 15000, 20, 20);
	woMIPNhitCone_BeamHalo->GetXaxis()->SetTitle("MIP NhitCone");
	woMIPNhitCone_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *woMIPisHalo_BeamHalo = new TH1F("woMIPisHalo_BeamHalo", "MIP isHalo from Beam Halo", 2000, 20, 20);
	woMIPisHalo_BeamHalo->GetXaxis()->SetTitle("MIP isHalo");
	woMIPisHalo_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *wophi_BeamHalo = new TH1F("wophi_BeamHalo", "Phi from Beam Halo", 2000, 20, 20);
	wophi_BeamHalo->GetXaxis()->SetTitle("Phi");
	wophi_BeamHalo->GetYaxis()->SetTitle("Entries");

	//------------------------------------------------------------------------------


	//----------------------All MIP variables for Beam Halo WITH MIP totE > 4.9 GeV--------------------

	TH1F *MIPTotalE_BeamHalo = new TH1F("MIPTotalE_BeamHalo", "MIP totE from Beam Halo", 2000, 20, 20);
	MIPTotalE_BeamHalo->GetXaxis()->SetTitle("MIP totE (GeV)");
	MIPTotalE_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *MIPChi2_BeamHalo = new TH1F("MIPChi2_BeamHalo", "MIP Chi2 from Beam Halo", 20000, 20, 20);
	MIPChi2_BeamHalo->GetXaxis()->SetTitle("MIP Chi2");
	MIPChi2_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *MIPSlope_BeamHalo = new TH1F("MIPSlope_BeamHalo", "MIP Slope from Beam Halo", 2000, 20, 20);
	MIPSlope_BeamHalo->GetXaxis()->SetTitle("MIP Slope");
	MIPSlope_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *MIPIntercept_BeamHalo = new TH1F("MIPIntercept_BeamHalo", "MIP Intercept from Beam Halo", 2000, 20, 20);
	MIPIntercept_BeamHalo->GetXaxis()->SetTitle("MIP Intercept");
	MIPIntercept_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *MIPNhitCone_BeamHalo = new TH1F("MIPNhitCone_BeamHalo", "MIP NhitCone from Beam Halo", 15000, 20, 20);
	MIPNhitCone_BeamHalo->GetXaxis()->SetTitle("MIP NhitCone");
	MIPNhitCone_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *MIPisHalo_BeamHalo = new TH1F("MIPisHalo_BeamHalo", "MIP isHalo from Beam Halo", 2000, 20, 20);
	MIPisHalo_BeamHalo->GetXaxis()->SetTitle("MIP isHalo");
	MIPisHalo_BeamHalo->GetYaxis()->SetTitle("Entries");

	//------------------------------------------------------------------------------

	//----------------------look into the most populated area from Beam Halo with MIP totE > 4.9 GeV----------------

	//get the phi distribution of beam halo
	TH1F *phi_BeamHalo = new TH1F("phi_BeamHalo", "Phi from Beam Halo", 2000, 20, 20);
	phi_BeamHalo->GetXaxis()->SetTitle("Phi");
	phi_BeamHalo->GetYaxis()->SetTitle("Entries");

	//from the most populated area
	TH1F *Phi_MIPTotalE_BeamHalo = new TH1F("Phi_MIPTotalE_BeamHalo", "MIP totE from Beam Halo", 2000, 20, 20);
	Phi_MIPTotalE_BeamHalo->GetXaxis()->SetTitle("MIP totE (GeV)");
	Phi_MIPTotalE_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *Phi_MIPChi2_BeamHalo = new TH1F("Phi_MIPChi2_BeamHalo", "MIP Chi2 from Beam Halo", 20000, 20, 20);
	Phi_MIPChi2_BeamHalo->GetXaxis()->SetTitle("MIP Chi2");
	Phi_MIPChi2_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *Phi_MIPSlope_BeamHalo = new TH1F("Phi_MIPSlope_BeamHalo", "MIP Slope from Beam Halo", 2000, 20, 20);
	Phi_MIPSlope_BeamHalo->GetXaxis()->SetTitle("MIP Slope");
	Phi_MIPSlope_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *Phi_MIPIntercept_BeamHalo = new TH1F("Phi_MIPIntercept_BeamHalo", "MIP Intercept from Beam Halo", 2000, 20, 20);
	Phi_MIPIntercept_BeamHalo->GetXaxis()->SetTitle("MIP Intercept");
	Phi_MIPIntercept_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *Phi_MIPNhitCone_BeamHalo = new TH1F("Phi_MIPNhitCone_BeamHalo", "MIP NhitCone from Beam Halo", 15000, 20, 20);
	Phi_MIPNhitCone_BeamHalo->GetXaxis()->SetTitle("MIP NhitCone");
	Phi_MIPNhitCone_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *Phi_MIPisHalo_BeamHalo = new TH1F("Phi_MIPisHalo_BeamHalo", "MIP isHalo from Beam Halo", 2000, 20, 20);
	Phi_MIPisHalo_BeamHalo->GetXaxis()->SetTitle("MIP isHalo");
	Phi_MIPisHalo_BeamHalo->GetYaxis()->SetTitle("Entries");

	TH1F *Phi_MIPCHi2_over_MIPNhitCone_BeamHalo = new TH1F("Phi_MIPCHi2_over_MIPNhitCone_BeamHalo", "MIPChi2/MIPNhitCone from Beam Halo", 2000, 20, 20);
	Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->GetXaxis()->SetTitle("MIPChi2/MIPNhitCone");
	Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->GetYaxis()->SetTitle("Entries");

	//----------------------------------------------------------------------------------------------
	

	/*
	//check run number for prompt and beamhalo(in prompt window)
	TH1F *Run_promptZTemp = new TH1F("Run_promptZTemp", "Run from prompt Template (Z)", 2000, 25, 25);
	Run_promptZTemp->GetXaxis()->SetTitle("Run");
	Run_promptZTemp->GetYaxis()->SetTitle("Entries");
	
	TH1F *Run_promptWTemp = new TH1F("Run_promptWTemp", "Run from prompt Template (W)", 2000, 25, 25);
	Run_promptWTemp->GetXaxis()->SetTitle("Run");
	Run_promptWTemp->GetYaxis()->SetTitle("Entries");

	TH1F *Run_BeamHaloprompt= new TH1F("Run_BeamHaloprompt", "Run from BeamHalo in promptregion", 2000, 25, 25);
	Run_BeamHaloprompt->GetXaxis()->SetTitle("Run");
	Run_BeamHaloprompt->GetYaxis()->SetTitle("Entries");

	ofstream RunZ, RunW, RunBeamhalo;
	RunZ.open("Run_Z.dat", ios::out);
	RunW.open("Run_W.dat", ios::out);
	RunBeamhalo.open("RunBeamhalo.dat", ios::out);
	*/


	cout << "Start analysis...." << endl;

	auto timenow_start = chrono::system_clock::to_time_t(chrono::system_clock::now());

	cout << "start time: " << ctime(&timenow_start) << " (CT)" << endl;

	ofstream file1, file2;
	file1.open("Analysis_log.txt", ios::out);
	file2.open("Beamhalo_issue_pick.txt, ios::out");
	file1 << "Analysis started at " << ctime(&timenow_start) << " (CT)" << endl;
	file1 << "Analysis in process....." << endl;

	//===========================no need to modify===================================
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	   // loop event begin
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	//for (Long64_t jentry = 0; jentry < 100000; jentry++)
	{
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		if ((jentry % 10000) == 0)
			// to print the number of processed entries
			std::cout << "Processed: " << jentry << std::endl;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		//Event Total for B: 14820000
		//Event Total for C: 36920000
		//Event Total for D: 9350000
		//Event Total for E: 18000000
		//Event Total for F: 26720000
		//Event Total for all: 105930000

		//Error in <TNetXNGFile::Open> : [ERROR] Server responded with an error : [3011] No servers are available to read the file.
		//Processed: 12250000
		//Processed: 12790000
		



		//===============================================================================


		//-----------------Analysis starts here--------------------

		//---------------full timing distribution------------------
		/*
		Prompt Z Template:
		basic photon ID selection
		without MET cut
		has pixel seed -> electron
		sieie or sipip > 0.01
		PhoSep > 0.2
		Require 2nd photon ET > 10 GeV
		Invariant mass ~ 91 GeV (85~100 window)
		*/

		//In-time photon
		vector<Int_t> itp;
		for (int iPho = 0; iPho < nPho; iPho++)
		{
			//in-time photon ID selection (basic)
			if ((*phoEt)[iPho] > 230.0 && fabs((*phoEta)[iPho]) < 1.442 && (*phoR9)[iPho] > 0.8)
			{
				itp.push_back(iPho);
				ITST->Fill((*phoSeedTime)[iPho]);
				FullTime->Fill((*phoSeedTime)[iPho]);

				//prompt Z template
				//if ((*phoMIPTotEnergy)[iPho] < 4.9 && (*phohasPixelSeed)[iPho] == 1 && ((*phoSigmaIEtaIEtaFull5x5)[iPho] > 0.001 || (*phoSigmaIPhiIPhiFull5x5)[iPho] > 0.001))
				if ((*phohasPixelSeed)[iPho] == 1 && (*phoSigmaIEtaIEtaFull5x5)[iPho] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[iPho] > 0.001)
				{
					float PhoSep = 999;
					for (int iiPho = iPho + 1; iiPho < nPho; iiPho++)
					{
						Float_t deltaR = DeltaR((*phoEta)[iPho], (*phoPhi)[iPho], (*phoEta)[iiPho], (*phoPhi)[iiPho]);
						if (deltaR < PhoSep)
						{
							PhoSep = deltaR;
						}


						if ((*phoMIPTotEnergy)[iPho] > 4.9)
						{
							//the second photon
							if ((*phoMIPTotEnergy)[iiPho] > 4.9 && (*phoEt)[iiPho] > 10.0 && fabs((*phoEta)[iiPho]) < 1.442 && (*phoR9)[iiPho] > 0.8 && (*phohasPixelSeed)[iiPho] == 1 && (*phoSigmaIEtaIEtaFull5x5)[iiPho] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[iiPho] > 0.001 && PhoSep > 0.2)
							{
								Float_t InvM = InvariMass((*phoEt)[iPho], (*phoEt)[iiPho], (*phoPhi)[iPho], (*phoPhi)[iiPho], (*phoEta)[iPho], (*phoEta)[iiPho]);
								InvMass_Z->Fill(InvM);

								if (InvM > 85 && InvM < 100)
								{
									promptZTemp->Fill((*phoSeedTime)[iPho]);
									MIPZ->Fill((*phoMIPTotEnergy)[iPho]);

									//check MIP variables
									MIPTotalE_Zprompt->Fill((*phoMIPTotEnergy)[iPho]);

									MIPChi2_Zprompt->Fill((*phoMIPChi2)[iPho]);

									MIPSlope_Zprompt->Fill((*phoMIPSlope)[iPho]);

									MIPIntercept_Zprompt->Fill((*phoMIPIntercept)[iPho]);

									MIPNhitCone_Zprompt->Fill((*phoMIPNhitCone)[iPho]);

									MIPisHalo_Zprompt->Fill((*phoMIPIsHalo)[iPho]);
								}
							}
						}						
					}
				}
			}
		}

		//out-of-time photon
		vector<Int_t> otp;
		for (int oPho = 0; oPho < onPho; oPho++)
		{
			//out-of-time photon ID selection (basic)
			if ((*ophoEt)[oPho] > 230.0 && fabs((*ophoEta)[oPho]) < 1.442 && (*ophoR9)[oPho] > 0.8)
			{
				otp.push_back(oPho);
				OOTST->Fill((*ophoSeedTime)[oPho]);
				FullTime->Fill((*ophoSeedTime)[oPho]);

				//prompt Z template
				//if ((*ophoMIPTotEnergy)[oPho] < 4.9 && (*ophohasPixelSeed)[oPho] == 1 && ((*ophoSigmaIEtaIEtaFull5x5)[oPho] > 0.001 || (*ophoSigmaIPhiIPhiFull5x5)[oPho] > 0.001))
				if ((*ophohasPixelSeed)[oPho] == 1 && (*ophoSigmaIEtaIEtaFull5x5)[oPho] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[oPho] > 0.001)
				{

					float PhoSep = 999;
					for (int ooPho = oPho + 1; ooPho < onPho; ooPho++)
					{
						Float_t deltaR = DeltaR((*ophoEta)[oPho], (*ophoPhi)[oPho], (*ophoEta)[ooPho], (*ophoPhi)[ooPho]);
						if (deltaR < PhoSep)
						{
							PhoSep = deltaR;
						}

						if ((*ophoMIPTotEnergy)[oPho] > 4.9)
						{
							//the second photon
							if ((*ophoMIPTotEnergy)[ooPho] > 4.9 && (*ophoEt)[ooPho] > 10.0 && fabs((*ophoEta)[ooPho]) < 1.442 && (*ophoR9)[ooPho] > 0.8 && (*ophohasPixelSeed)[ooPho] == 1 && (*ophoSigmaIEtaIEtaFull5x5)[ooPho] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[ooPho] > 0.001 && PhoSep > 0.2)
							{
								Float_t InvM = InvariMass((*ophoEt)[oPho], (*ophoEt)[ooPho], (*ophoPhi)[oPho], (*ophoPhi)[ooPho], (*ophoEta)[oPho], (*ophoEta)[ooPho]);
								InvMass_Z->Fill(InvM);

								if (InvM > 85 && InvM < 100)
								{
									promptZTemp->Fill((*ophoSeedTime)[oPho]);
									MIPZ->Fill((*ophoMIPTotEnergy)[oPho]);

									//check MIP variables
									MIPTotalE_Zprompt->Fill((*ophoMIPTotEnergy)[oPho]);

									MIPChi2_Zprompt->Fill((*ophoMIPChi2)[oPho]);

									MIPSlope_Zprompt->Fill((*ophoMIPSlope)[oPho]);

									MIPIntercept_Zprompt->Fill((*ophoMIPIntercept)[oPho]);

									MIPNhitCone_Zprompt->Fill((*ophoMIPNhitCone)[oPho]);

									MIPisHalo_Zprompt->Fill((*ophoMIPIsHalo)[oPho]);
								}
							}
						}
						
					}
				}
			}
		}



		if (itp.size() == 0 && otp.size() == 0) continue;



		Bool_t isOOT = kFALSE; //if it's OOT photon, select the OOT, can be true or false but false makes more sense.
		Int_t IDXSEL = -1; //set index for photon, initial value could be something non-physical, negative something.

		if (itp.size() > 0 && otp.size() == 0)
		{
			IDXSEL = itp[0];
			FullTime->Fill((*phoSeedTime)[IDXSEL]);

			//prompt Z template
			//if ((*phoMIPTotEnergy)[iPho] < 4.9 && (*phohasPixelSeed)[iPho] == 1 && ((*phoSigmaIEtaIEtaFull5x5)[iPho] > 0.001 || (*phoSigmaIPhiIPhiFull5x5)[iPho] > 0.001))
			if ((*phohasPixelSeed)[IDXSEL] == 1 && (*phoSigmaIEtaIEtaFull5x5)[IDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[IDXSEL] > 0.001)
			{
				float PhoSep = 999;
				for (int iiPho = IDXSEL + 1; iiPho < itp.size(); iiPho++)
				{
					Float_t deltaR = DeltaR((*phoEta)[IDXSEL], (*phoPhi)[IDXSEL], (*phoEta)[iiPho], (*phoPhi)[iiPho]);
					if (deltaR < PhoSep)
					{
						PhoSep = deltaR;
					}

					if ((*phoMIPTotEnergy)[IDXSEL] > 4.9)
					{
						//the second photon
						if ((*phoMIPTotEnergy)[iiPho] > 4.9 && (*phoEt)[iiPho] > 10.0 && fabs((*phoEta)[iiPho]) < 1.442 && (*phoR9)[iiPho] > 0.8 && (*phohasPixelSeed)[iiPho] == 1 && (*phoSigmaIEtaIEtaFull5x5)[iiPho] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[iiPho] > 0.001 && PhoSep > 0.2)
						{
							Float_t InvM = InvariMass((*phoEt)[IDXSEL], (*phoEt)[iiPho], (*phoPhi)[IDXSEL], (*phoPhi)[iiPho], (*phoEta)[IDXSEL], (*phoEta)[iiPho]);
							InvMass_Z->Fill(InvM);

							if (InvM > 85 && InvM < 100)
							{
								promptZTemp->Fill((*phoSeedTime)[IDXSEL]);
								MIPZ->Fill((*phoMIPTotEnergy)[IDXSEL]);

								//check MIP variables
								MIPTotalE_Zprompt->Fill((*phoMIPTotEnergy)[IDXSEL]);

								MIPChi2_Zprompt->Fill((*phoMIPChi2)[IDXSEL]);

								MIPSlope_Zprompt->Fill((*phoMIPSlope)[IDXSEL]);

								MIPIntercept_Zprompt->Fill((*phoMIPIntercept)[IDXSEL]);

								MIPNhitCone_Zprompt->Fill((*phoMIPNhitCone)[IDXSEL]);

								MIPisHalo_Zprompt->Fill((*phoMIPIsHalo)[IDXSEL]);
							}
						}
					}
				}
			}
		}

		if (itp.size() == 0 && otp.size() > 0)
		{
			IDXSEL = otp[0];
			isOOT = kTRUE;
			FullTime->Fill((*ophoSeedTime)[IDXSEL]);

			//prompt Z template
			//if ((*ophoMIPTotEnergy)[oPho] < 4.9 && (*ophohasPixelSeed)[oPho] == 1 && ((*ophoSigmaIEtaIEtaFull5x5)[oPho] > 0.001 || (*ophoSigmaIPhiIPhiFull5x5)[oPho] > 0.001))
			if ((*ophohasPixelSeed)[IDXSEL] == 1 && (*ophoSigmaIEtaIEtaFull5x5)[IDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[IDXSEL] > 0.001)
			{
				float PhoSep = 999;
				for (int ooPho = IDXSEL + 1; ooPho < otp.size(); ooPho++)
				{
					Float_t deltaR = DeltaR((*ophoEta)[IDXSEL], (*ophoPhi)[IDXSEL], (*ophoEta)[ooPho], (*ophoPhi)[ooPho]);
					if (deltaR < PhoSep)
					{
						PhoSep = deltaR;
					}

					if ((*ophoMIPTotEnergy)[IDXSEL] > 4.9)
					{
						//the second photo
						if ((*ophoMIPTotEnergy)[ooPho] > 4.9 && (*ophoEt)[ooPho] > 10.0 && fabs((*ophoEta)[ooPho]) < 1.442 && (*ophoR9)[ooPho] > 0.8 && (*ophohasPixelSeed)[ooPho] == 1 && (*ophoSigmaIEtaIEtaFull5x5)[ooPho] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[ooPho] > 0.001 && PhoSep > 0.2)
						{
							Float_t InvM = InvariMass((*ophoEt)[IDXSEL], (*ophoEt)[ooPho], (*ophoPhi)[IDXSEL], (*ophoPhi)[ooPho], (*ophoEta)[IDXSEL], (*ophoEta)[ooPho]);
							InvMass_Z->Fill(InvM);

							if (InvM > 85 && InvM < 100)
							{
								promptZTemp->Fill((*ophoSeedTime)[IDXSEL]);
								MIPZ->Fill((*ophoMIPTotEnergy)[IDXSEL]);
								
								//check MIP variables
								MIPTotalE_Zprompt->Fill((*ophoMIPTotEnergy)[IDXSEL]);

								MIPChi2_Zprompt->Fill((*ophoMIPChi2)[IDXSEL]);

								MIPSlope_Zprompt->Fill((*ophoMIPSlope)[IDXSEL]);

								MIPIntercept_Zprompt->Fill((*ophoMIPIntercept)[IDXSEL]);

								MIPNhitCone_Zprompt->Fill((*ophoMIPNhitCone)[IDXSEL]);

								MIPisHalo_Zprompt->Fill((*ophoMIPIsHalo)[IDXSEL]);
							}
						}
					}
				}
			}
		}

		if (itp.size() > 0 && otp.size() > 0)
		{
			IDXSEL = otp[0];
			isOOT = kTRUE;
			FullTime->Fill((*ophoSeedTime)[IDXSEL]);

			//prompt Z template
			//if ((*ophoMIPTotEnergy)[oPho] < 4.9 && (*ophohasPixelSeed)[oPho] == 1 && ((*ophoSigmaIEtaIEtaFull5x5)[oPho] > 0.001 || (*ophoSigmaIPhiIPhiFull5x5)[oPho] > 0.001))
			if ((*ophohasPixelSeed)[IDXSEL] == 1 && (*ophoSigmaIEtaIEtaFull5x5)[IDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[IDXSEL] > 0.001)
			{
				float PhoSep = 999;
				for (int ooPho = IDXSEL + 1; ooPho < otp.size(); ooPho++)
				{
					Float_t deltaR = DeltaR((*ophoEta)[IDXSEL], (*ophoPhi)[IDXSEL], (*ophoEta)[ooPho], (*ophoPhi)[ooPho]);
					if (deltaR < PhoSep)
					{
						PhoSep = deltaR;
					}

					if ((*ophoMIPTotEnergy)[IDXSEL] > 4.9)
					{
						//the second photon
						 
						if ((*ophoMIPTotEnergy)[ooPho] > 4.9 && (*ophoEt)[ooPho] > 10.0 && fabs((*ophoEta)[ooPho]) < 1.442 && (*ophoR9)[ooPho] > 0.8 && (*ophohasPixelSeed)[ooPho] == 1 && (*ophoSigmaIEtaIEtaFull5x5)[ooPho] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[ooPho] > 0.001 && PhoSep > 0.2)
						{
							Float_t InvM = InvariMass((*ophoEt)[IDXSEL], (*ophoEt)[ooPho], (*ophoPhi)[IDXSEL], (*ophoPhi)[ooPho], (*ophoEta)[IDXSEL], (*ophoEta)[ooPho]);
							InvMass_Z->Fill(InvM);

							if (InvM > 85 && InvM < 100)
							{
								promptZTemp->Fill((*ophoSeedTime)[IDXSEL]);
								MIPZ->Fill((*ophoMIPTotEnergy)[IDXSEL]);

								//check MIP variables
								MIPTotalE_Zprompt->Fill((*ophoMIPTotEnergy)[IDXSEL]);

								MIPChi2_Zprompt->Fill((*ophoMIPChi2)[IDXSEL]);

								MIPSlope_Zprompt->Fill((*ophoMIPSlope)[IDXSEL]);

								MIPIntercept_Zprompt->Fill((*ophoMIPIntercept)[IDXSEL]);

								MIPNhitCone_Zprompt->Fill((*ophoMIPNhitCone)[IDXSEL]);

								MIPisHalo_Zprompt->Fill((*ophoMIPIsHalo)[IDXSEL]);
							}
						}
					}
				}
			}
		}
		//-----------------------------Done timing distribution before MET & Z prompt template--------------------------------------------------


		//---------------full timing distribution with MET cut------------------
		/*
		Candidate Events:
		basic photon ID selection
		MET cut > 210 GeV
		has no pixel seed -> photon
		sieie or sipip > 0.01 for spike cleaning
		sieie < 0.01015 for discriminating jet
		MIP TotE < 4.9 GeV for Beam Halo reducing
		*/

		/*
		Prompt W Template:
		basic photon ID selection
		MET cut > 210 GeV
		has pixel seed -> electron
		sieie or sipip > 0.01
		PhoSep > 0.2
		Require 2nd photon ET > 10 GeV
		Invariant mass ~ 80 GeV (70~90 window)
		*/

		/*
		Spike Template:
		basic photon ID selection
		MET cut > 210 GeV
		has no pixel seed -> photon
		sieie or sipip < 0.01
		MIP TotE < 4.9 GeV for Beam Halo reducing
		*/

		/*
		Beam Halo Template:
		basic photon ID selection
		MET cut > 210 GeV
		has no pixel seed -> photon
		sieie or sipip > 0.01
		MIP TotE > 4.9 GeV 
		*/

		
		//In-time photon with MET
		vector<Int_t> mitp;
		for (int iPho = 0; iPho < nPho; iPho++)
		{
			//in-time photon ID selection (basic)
			if ((*phoEt)[iPho] > 230.0 && fabs((*phoEta)[iPho]) < 1.442 && (*phoR9)[iPho] > 0.8 && pfMET > 210)
			{
				mitp.push_back(iPho);
				MITST->Fill((*phoSeedTime)[iPho]);
				MFullTime->Fill((*phoSeedTime)[iPho]);

				//Candidate Events
				if ((*phohasPixelSeed)[iPho] == 0 && (*phoMIPTotEnergy)[iPho] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[iPho] < 0.01015 && (*phoSigmaIEtaIEtaFull5x5)[iPho] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[iPho] > 0.001)
				{
					Candidate_Events->Fill((*phoSeedTime)[iPho]);
					metfilter_candidate->Fill(metFilters);
				}

				//Spike Template
				if ((*phohasPixelSeed)[iPho] == 0 && (*phoMIPTotEnergy)[iPho] < 4.9 && ((*phoSigmaIEtaIEtaFull5x5)[iPho] < 0.001 || (*phoSigmaIPhiIPhiFull5x5)[iPho] < 0.001))
				{
					SpikeTemp->Fill((*phoSeedTime)[iPho]);
				}

				//Beam Halo Template
				if ((*phohasPixelSeed)[iPho] == 0 && (*phoSigmaIEtaIEtaFull5x5)[iPho] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[iPho] > 0.001)
				{
					wophi_BeamHalo->Fill((*phoPhi)[iPho]);

					woMIPTotalE_BeamHalo->Fill((*phoMIPTotEnergy)[iPho]);

					woMIPChi2_BeamHalo->Fill((*phoMIPChi2)[iPho]);

					woMIPSlope_BeamHalo->Fill((*phoMIPSlope)[iPho]);

					woMIPIntercept_BeamHalo->Fill((*phoMIPIntercept)[iPho]);
					
					woMIPNhitCone_BeamHalo->Fill((*phoMIPNhitCone)[iPho]);

					woMIPisHalo_BeamHalo->Fill((*phoMIPIsHalo)[iPho]);

					Float_t MIPChi2_over_MIPNhitCone = ((*phoMIPChi2)[iPho] / (*phoMIPNhitCone)[iPho]);
					if ((*phoMIPTotEnergy)[iPho] > 4.9 && MIPChi2_over_MIPNhitCone < 0.2)
					{
						BeamHaloTemp->Fill((*phoSeedTime)[iPho]);

						metfilter_beamhalo->Fill(metFilters);

						BeamHaloSieie_vs_time->Fill((*phoSeedTime)[iPho], (*phoSigmaIEtaIEtaFull5x5)[iPho]);

						if ((*phoSeedTime)[iPho] > -2.0 && (*phoSeedTime)[iPho] < 2.0 && (*phoSigmaIEtaIEtaFull5x5)[iPho] > 0.0 && (*phoSigmaIEtaIEtaFull5x5)[iPho] < 0.015)
						{
							BeamHalo_prompt_region->Fill((*phoSeedTime)[iPho], (*phoSigmaIEtaIEtaFull5x5)[iPho]);
							file2 << run << " " << event << " " << lumis << endl;
						}

						BeamHaloSieie->Fill((*phoSigmaIEtaIEtaFull5x5)[iPho]);

						phi_BeamHalo->Fill((*phoPhi)[iPho]);

						MIPTotalE_BeamHalo->Fill((*phoMIPTotEnergy)[iPho]);

						MIPChi2_BeamHalo->Fill((*phoMIPChi2)[iPho]);

						MIPSlope_BeamHalo->Fill((*phoMIPSlope)[iPho]);

						MIPIntercept_BeamHalo->Fill((*phoMIPIntercept)[iPho]);

						MIPNhitCone_BeamHalo->Fill((*phoMIPNhitCone)[iPho]);

						MIPisHalo_BeamHalo->Fill((*phoMIPIsHalo)[iPho]);

						//check MIP variables
						if ((*phoPhi)[iPho] > -3.2 && (*phoPhi)[iPho] < -2.8)
						{
							Phi_MIPTotalE_BeamHalo->Fill((*phoMIPTotEnergy)[iPho]);

							Phi_MIPChi2_BeamHalo->Fill((*phoMIPChi2)[iPho]);

							Phi_MIPSlope_BeamHalo->Fill((*phoMIPSlope)[iPho]);

							Phi_MIPIntercept_BeamHalo->Fill((*phoMIPIntercept)[iPho]);

							Phi_MIPNhitCone_BeamHalo->Fill((*phoMIPNhitCone)[iPho]);

							Phi_MIPisHalo_BeamHalo->Fill((*phoMIPIsHalo)[iPho]);

							Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*phoMIPChi2)[iPho] / (*phoMIPNhitCone)[iPho]);

							
						}

						if ((*phoPhi)[iPho] > -0.5 && (*phoPhi)[iPho] < 0.5)
						{
							Phi_MIPTotalE_BeamHalo->Fill((*phoMIPTotEnergy)[iPho]);

							Phi_MIPChi2_BeamHalo->Fill((*phoMIPChi2)[iPho]);

							Phi_MIPSlope_BeamHalo->Fill((*phoMIPSlope)[iPho]);

							Phi_MIPIntercept_BeamHalo->Fill((*phoMIPIntercept)[iPho]);

							Phi_MIPNhitCone_BeamHalo->Fill((*phoMIPNhitCone)[iPho]);

							Phi_MIPisHalo_BeamHalo->Fill((*phoMIPIsHalo)[iPho]);

							Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*phoMIPChi2)[iPho] / (*phoMIPNhitCone)[iPho]);
						}

						if ((*phoPhi)[iPho] > 2.8 && (*phoPhi)[iPho] < 3.2)
						{
							Phi_MIPTotalE_BeamHalo->Fill((*phoMIPTotEnergy)[iPho]);

							Phi_MIPChi2_BeamHalo->Fill((*phoMIPChi2)[iPho]);

							Phi_MIPSlope_BeamHalo->Fill((*phoMIPSlope)[iPho]);

							Phi_MIPIntercept_BeamHalo->Fill((*phoMIPIntercept)[iPho]);

							Phi_MIPNhitCone_BeamHalo->Fill((*phoMIPNhitCone)[iPho]);

							Phi_MIPisHalo_BeamHalo->Fill((*phoMIPIsHalo)[iPho]);

							Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*phoMIPChi2)[iPho] / (*phoMIPNhitCone)[iPho]);
						}

					}
					
				}

			
				//prompt W template
				if ((*phohasPixelSeed)[iPho] == 1 && (*phoSigmaIEtaIEtaFull5x5)[iPho] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[iPho] > 0.001)
				{
					if ((*phoMIPTotEnergy)[iPho] > 4.9)
					{
						promptWTemp->Fill((*phoSeedTime)[iPho]);
						Float_t TransMW = TransMassW((*phoEt)[iPho], (*phoPhi)[iPho], (*phoEta)[iPho]);
						TransMass_W->Fill(TransMW);
						MIPW->Fill((*phoMIPTotEnergy)[iPho]);

						//check MIP variables
						MIPTotalE_Wprompt->Fill((*phoMIPTotEnergy)[iPho]);

						MIPChi2_Wprompt->Fill((*phoMIPChi2)[iPho]);

						MIPSlope_Wprompt->Fill((*phoMIPSlope)[iPho]);

						MIPIntercept_Wprompt->Fill((*phoMIPIntercept)[iPho]);

						MIPNhitCone_Wprompt->Fill((*phoMIPNhitCone)[iPho]);

						MIPisHalo_Wprompt->Fill((*phoMIPIsHalo)[iPho]);

						MIPCHi2_over_MIPNhitCone_Wprompt->Fill((*phoMIPChi2)[iPho] / (*phoMIPNhitCone)[iPho]);
					}
				}
			}
		}

		//Out-of-time photon with MET 
		vector<Int_t> motp;
		for (int oPho = 0; oPho < onPho; oPho++) 
		{
			Float_t NewMET = newuMET(pfMET, pfMETPhi, (*ophoPhi)[oPho], (*ophoEt)[oPho]);
			if (NewMET > 210 && (*ophoEt)[oPho] > 230.0 && fabs((*ophoEta)[oPho]) < 1.442 && (*ophoR9)[oPho] > 0.8)
			{
				motp.push_back(oPho);
				MOOTST->Fill((*ophoSeedTime)[oPho]);
				MFullTime->Fill((*ophoSeedTime)[oPho]);

				//Candidate Events
				if ((*ophohasPixelSeed)[oPho] == 0 && (*ophoMIPTotEnergy)[oPho] < 4.9 && (*ophoSigmaIEtaIEtaFull5x5)[oPho] < 0.01015 && (*ophoSigmaIEtaIEtaFull5x5)[oPho] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[oPho] > 0.001)
				{
					Candidate_Events->Fill((*ophoSeedTime)[oPho]);
					metfilter_candidate->Fill(metFilters);
				}

				//Spike Template
				if ((*ophohasPixelSeed)[oPho] == 0 && (*ophoMIPTotEnergy)[oPho] < 4.9 && ((*ophoSigmaIEtaIEtaFull5x5)[oPho] < 0.001 || (*ophoSigmaIPhiIPhiFull5x5)[oPho] < 0.001))
				{
					SpikeTemp->Fill((*ophoSeedTime)[oPho]);
				}

				//Beam Halo Template
				if ((*ophohasPixelSeed)[oPho] == 0 && (*ophoSigmaIEtaIEtaFull5x5)[oPho] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[oPho] > 0.001)
				{
					woMIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[oPho]);

					woMIPChi2_BeamHalo->Fill((*ophoMIPChi2)[oPho]);

					woMIPSlope_BeamHalo->Fill((*ophoMIPSlope)[oPho]);

					woMIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[oPho]);

					woMIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[oPho]);

					woMIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[oPho]);

					wophi_BeamHalo->Fill((*ophoPhi)[oPho]);

					Float_t MIPChi2_over_MIPNhitCone = ((*ophoMIPChi2)[oPho] / (*ophoMIPNhitCone)[oPho]);
					if ((*ophoMIPTotEnergy)[oPho] > 4.9 && MIPChi2_over_MIPNhitCone < 0.2)
					{
						BeamHaloTemp->Fill((*ophoSeedTime)[oPho]);

						metfilter_beamhalo->Fill(metFilters);

						BeamHaloSieie_vs_time->Fill((*ophoSeedTime)[oPho], (*ophoSigmaIEtaIEtaFull5x5)[oPho]);


						if ((*ophoSeedTime)[oPho] > -2.0 && (*ophoSeedTime)[oPho] < 2.0 && (*ophoSigmaIEtaIEtaFull5x5)[oPho] > 0.0 && (*ophoSigmaIEtaIEtaFull5x5)[oPho] < 0.015)
						{
							BeamHalo_prompt_region->Fill((*ophoSeedTime)[oPho], (*ophoSigmaIEtaIEtaFull5x5)[oPho]);
							file2 << run << " " << event << " " << lumis << endl;
						}

						BeamHaloSieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[oPho]);

						phi_BeamHalo->Fill((*ophoPhi)[oPho]);

						MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[oPho]);

						MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[oPho]);

						MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[oPho]);

						MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[oPho]);

						MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[oPho]);

						MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[oPho]);

						phi_BeamHalo->Fill((*ophoPhi)[oPho]);

						//check MIP variables
						if ((*ophoPhi)[oPho] > -3.2 && (*ophoPhi)[oPho] < -2.8)
						{
							Phi_MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[oPho]);

							Phi_MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[oPho]);

							Phi_MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[oPho]);

							Phi_MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[oPho]);

							Phi_MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[oPho]);

							Phi_MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[oPho]);

							Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*ophoMIPChi2)[oPho] / (*ophoMIPNhitCone)[oPho]);
						}

						if ((*ophoPhi)[oPho] > -0.5 && (*ophoPhi)[oPho] < 0.5)
						{
							Phi_MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[oPho]);

							Phi_MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[oPho]);

							Phi_MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[oPho]);

							Phi_MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[oPho]);

							Phi_MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[oPho]);

							Phi_MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[oPho]);

							Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*ophoMIPChi2)[oPho] / (*ophoMIPNhitCone)[oPho]);
						}

						if ((*ophoPhi)[oPho] > 2.8 && (*ophoPhi)[oPho] < 3.2)
						{
							Phi_MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[oPho]);

							Phi_MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[oPho]);

							Phi_MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[oPho]);

							Phi_MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[oPho]);

							Phi_MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[oPho]);

							Phi_MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[oPho]);

							Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*ophoMIPChi2)[oPho] / (*ophoMIPNhitCone)[oPho]);
						}

					}
				}

				//prompt W template
				if ((*ophohasPixelSeed)[oPho] == 1 && (*ophoSigmaIEtaIEtaFull5x5)[oPho] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[oPho] > 0.001)
				{
					if ((*ophoMIPTotEnergy)[oPho] > 4.9)
					{
						promptWTemp->Fill((*ophoSeedTime)[oPho]);
						MIPW->Fill((*ophoMIPTotEnergy)[oPho]);
						Float_t TransMW = TransMassW((*ophoEt)[oPho], (*ophoPhi)[oPho], (*ophoEta)[oPho]);
						TransMass_W->Fill(TransMW);

						//check MIP variables
						MIPTotalE_Wprompt->Fill((*ophoMIPTotEnergy)[oPho]);

						MIPChi2_Wprompt->Fill((*ophoMIPChi2)[oPho]);

						MIPSlope_Wprompt->Fill((*ophoMIPSlope)[oPho]);

						MIPIntercept_Wprompt->Fill((*ophoMIPIntercept)[oPho]);

						MIPNhitCone_Wprompt->Fill((*ophoMIPNhitCone)[oPho]);

						MIPisHalo_Wprompt->Fill((*ophoMIPIsHalo)[oPho]);

						MIPCHi2_over_MIPNhitCone_Wprompt->Fill((*ophoMIPChi2)[oPho] / (*ophoMIPNhitCone)[oPho]);
					}
				}
			}
		}

		Bool_t misOOT = kFALSE; //if it's OOT photon, select the OOT, can be true or false but false makes more sense.
		Int_t mIDXSEL = -1; //set index for photon, initial value could be something non-physical, negative something.
		if (motp.size() == 0 && mitp.size() == 0) continue;

		if (motp.size() == 0 && mitp.size() > 0)
		{
			mIDXSEL = mitp[0];
			MFullTime->Fill((*phoSeedTime)[mIDXSEL]);

			//Candidate Events
			if ((*phohasPixelSeed)[mIDXSEL] == 0 && (*phoMIPTotEnergy)[mIDXSEL] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.01015 && (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
			{
				Candidate_Events->Fill((*phoSeedTime)[mIDXSEL]);
				metfilter_candidate->Fill(metFilters);
			}

			//Spike Template
			if ((*phohasPixelSeed)[mIDXSEL] == 0 && (*phoMIPTotEnergy)[mIDXSEL] < 4.9 && ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.001 || (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] < 0.001))
			{
				SpikeTemp->Fill((*phoSeedTime)[mIDXSEL]);
			}

			//Beam Halo Template
			if ((*phohasPixelSeed)[mIDXSEL] == 0 && (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
			{
				wophi_BeamHalo->Fill((*phoPhi)[mIDXSEL]);

				woMIPTotalE_BeamHalo->Fill((*phoMIPTotEnergy)[mIDXSEL]);

				woMIPChi2_BeamHalo->Fill((*phoMIPChi2)[mIDXSEL]);

				woMIPSlope_BeamHalo->Fill((*phoMIPSlope)[mIDXSEL]);

				woMIPIntercept_BeamHalo->Fill((*phoMIPIntercept)[mIDXSEL]);

				woMIPNhitCone_BeamHalo->Fill((*phoMIPNhitCone)[mIDXSEL]);

				woMIPisHalo_BeamHalo->Fill((*phoMIPIsHalo)[mIDXSEL]);


				Float_t MIPChi2_over_MIPNhitCone = ((*phoMIPChi2)[mIDXSEL] / (*phoMIPNhitCone)[mIDXSEL]);
				if ((*phoMIPTotEnergy)[mIDXSEL] > 4.9 && MIPChi2_over_MIPNhitCone < 0.2)
				{
					BeamHaloTemp->Fill((*phoSeedTime)[mIDXSEL]);

					metfilter_beamhalo->Fill(metFilters);

					BeamHaloSieie_vs_time->Fill((*phoSeedTime)[mIDXSEL], (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]);

					if ((*phoSeedTime)[mIDXSEL] > -2.0 && (*phoSeedTime)[mIDXSEL] < 2.0 && (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.0 && (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.015)
					{
						BeamHalo_prompt_region->Fill((*phoSeedTime)[mIDXSEL], (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]);
						file2 << run << " " << event << " " << lumis << endl;
					}

					BeamHaloSieie->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]);
									
					phi_BeamHalo->Fill((*phoPhi)[mIDXSEL]);

					MIPTotalE_BeamHalo->Fill((*phoMIPTotEnergy)[mIDXSEL]);

					MIPChi2_BeamHalo->Fill((*phoMIPChi2)[mIDXSEL]);

					MIPSlope_BeamHalo->Fill((*phoMIPSlope)[mIDXSEL]);

					MIPIntercept_BeamHalo->Fill((*phoMIPIntercept)[mIDXSEL]);

					MIPNhitCone_BeamHalo->Fill((*phoMIPNhitCone)[mIDXSEL]);

					MIPisHalo_BeamHalo->Fill((*phoMIPIsHalo)[mIDXSEL]);

					//check the MIP variables
					if ((*phoPhi)[mIDXSEL] > -3.2 && (*phoPhi)[mIDXSEL] < -2.8)
					{
						Phi_MIPTotalE_BeamHalo->Fill((*phoMIPTotEnergy)[mIDXSEL]);

						Phi_MIPChi2_BeamHalo->Fill((*phoMIPChi2)[mIDXSEL]);

						Phi_MIPSlope_BeamHalo->Fill((*phoMIPSlope)[mIDXSEL]);

						Phi_MIPIntercept_BeamHalo->Fill((*phoMIPIntercept)[mIDXSEL]);

						Phi_MIPNhitCone_BeamHalo->Fill((*phoMIPNhitCone)[mIDXSEL]);

						Phi_MIPisHalo_BeamHalo->Fill((*phoMIPIsHalo)[mIDXSEL]);

						Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*phoMIPChi2)[mIDXSEL] / (*phoMIPNhitCone)[mIDXSEL]);
					}

					if ((*phoPhi)[mIDXSEL] > -0.5 && (*phoPhi)[mIDXSEL] < 0.5)
					{
						Phi_MIPTotalE_BeamHalo->Fill((*phoMIPTotEnergy)[mIDXSEL]);

						Phi_MIPChi2_BeamHalo->Fill((*phoMIPChi2)[mIDXSEL]);

						Phi_MIPSlope_BeamHalo->Fill((*phoMIPSlope)[mIDXSEL]);

						Phi_MIPIntercept_BeamHalo->Fill((*phoMIPIntercept)[mIDXSEL]);

						Phi_MIPNhitCone_BeamHalo->Fill((*phoMIPNhitCone)[mIDXSEL]);

						Phi_MIPisHalo_BeamHalo->Fill((*phoMIPIsHalo)[mIDXSEL]);

						Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*phoMIPChi2)[mIDXSEL] / (*phoMIPNhitCone)[mIDXSEL]);
					}

					if ((*phoPhi)[mIDXSEL] > 2.8 && (*phoPhi)[mIDXSEL] < 3.2)
					{
						Phi_MIPTotalE_BeamHalo->Fill((*phoMIPTotEnergy)[mIDXSEL]);

						Phi_MIPChi2_BeamHalo->Fill((*phoMIPChi2)[mIDXSEL]);

						Phi_MIPSlope_BeamHalo->Fill((*phoMIPSlope)[mIDXSEL]);

						Phi_MIPIntercept_BeamHalo->Fill((*phoMIPIntercept)[mIDXSEL]);

						Phi_MIPNhitCone_BeamHalo->Fill((*phoMIPNhitCone)[mIDXSEL]);

						Phi_MIPisHalo_BeamHalo->Fill((*phoMIPIsHalo)[mIDXSEL]);

						Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*phoMIPChi2)[mIDXSEL] / (*phoMIPNhitCone)[mIDXSEL]);
					}
				}			
			}

			//prompt W template
			if ((*phohasPixelSeed)[mIDXSEL] == 1 && (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
			{
				

				if ((*phoMIPTotEnergy)[mIDXSEL] > 4.9)
				{
					promptWTemp->Fill((*phoSeedTime)[mIDXSEL]);
					MIPW->Fill((*phoMIPTotEnergy)[mIDXSEL]);
					Float_t TransMW = TransMassW((*phoEt)[mIDXSEL], (*phoPhi)[mIDXSEL], (*phoEta)[mIDXSEL]);
					TransMass_W->Fill(TransMW);

					//check MIP variables
					MIPTotalE_Wprompt->Fill((*phoMIPTotEnergy)[mIDXSEL]);

					MIPChi2_Wprompt->Fill((*phoMIPChi2)[mIDXSEL]);

					MIPSlope_Wprompt->Fill((*phoMIPSlope)[mIDXSEL]);

					MIPIntercept_Wprompt->Fill((*phoMIPIntercept)[mIDXSEL]);

					MIPNhitCone_Wprompt->Fill((*phoMIPNhitCone)[mIDXSEL]);

					MIPisHalo_Wprompt->Fill((*phoMIPIsHalo)[mIDXSEL]);

					MIPCHi2_over_MIPNhitCone_Wprompt->Fill((*phoMIPChi2)[mIDXSEL] / (*phoMIPNhitCone)[mIDXSEL]);
				}
			}
		}

		if (motp.size() > 0 && mitp.size() == 0)
		{
			mIDXSEL = motp[0];
			misOOT = kTRUE;
			MFullTime->Fill((*ophoSeedTime)[mIDXSEL]);

			//Candidate Events
			if ((*ophohasPixelSeed)[mIDXSEL] == 0 && (*ophoMIPTotEnergy)[mIDXSEL] < 4.9 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.01015 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
			{
				Candidate_Events->Fill((*ophoSeedTime)[mIDXSEL]);

				metfilter_candidate->Fill(metFilters);
			}

			//Spike Template
			if ((*ophohasPixelSeed)[mIDXSEL] == 0 && (*ophoMIPTotEnergy)[mIDXSEL] < 4.9 && ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.001 || (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] < 0.001))
			{
				SpikeTemp->Fill((*ophoSeedTime)[mIDXSEL]);
			}

			//Beam Halo Template
			if ((*ophohasPixelSeed)[mIDXSEL] == 0 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
			{
				woMIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

				woMIPChi2_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL]);

				woMIPSlope_BeamHalo->Fill((*ophoMIPSlope)[mIDXSEL]);

				woMIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[mIDXSEL]);

				woMIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[mIDXSEL]);

				woMIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[mIDXSEL]);

				wophi_BeamHalo->Fill((*ophoPhi)[mIDXSEL]);

				Float_t MIPChi2_over_MIPNhitCone = ((*ophoMIPChi2)[mIDXSEL] / (*ophoMIPNhitCone)[mIDXSEL]);
				if ((*ophoMIPTotEnergy)[mIDXSEL] > 4.9 && MIPChi2_over_MIPNhitCone < 0.2)
				{
					BeamHaloTemp->Fill((*ophoSeedTime)[mIDXSEL]);

					metfilter_beamhalo->Fill(metFilters);

					BeamHaloSieie_vs_time->Fill((*ophoSeedTime)[mIDXSEL], (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]);

					if ((*ophoSeedTime)[mIDXSEL] > -2.0 && (*ophoSeedTime)[mIDXSEL] < 2.0 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.0 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.015)
					{
						BeamHalo_prompt_region->Fill((*ophoSeedTime)[mIDXSEL], (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]);
						file2 << run << " " << event << " " << lumis << endl;
					}

					BeamHaloSieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]);

					MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

					MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL]);

					MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[mIDXSEL]);

					MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[mIDXSEL]);

					MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[mIDXSEL]);

					MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[mIDXSEL]);

					phi_BeamHalo->Fill((*ophoPhi)[mIDXSEL]);

					//check MIP variables
					if ((*ophoPhi)[mIDXSEL] > -3.2 && (*ophoPhi)[mIDXSEL] < -2.8)
					{
						Phi_MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

						Phi_MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL]);

						Phi_MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[mIDXSEL]);

						Phi_MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[mIDXSEL]);

						Phi_MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[mIDXSEL]);

						Phi_MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[mIDXSEL]);

						Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL] / (*ophoMIPNhitCone)[mIDXSEL]);
					}

					if ((*ophoPhi)[mIDXSEL] > -0.5 && (*ophoPhi)[mIDXSEL] < 0.5)
					{
						Phi_MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

						Phi_MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL]);

						Phi_MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[mIDXSEL]);

						Phi_MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[mIDXSEL]);

						Phi_MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[mIDXSEL]);

						Phi_MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[mIDXSEL]);
				
						Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL] / (*ophoMIPNhitCone)[mIDXSEL]);
					}

					if ((*ophoPhi)[mIDXSEL] > 2.8 && (*ophoPhi)[mIDXSEL] < 3.2)
					{
						Phi_MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

						Phi_MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL]);

						Phi_MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[mIDXSEL]);

						Phi_MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[mIDXSEL]);

						Phi_MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[mIDXSEL]);

						Phi_MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[mIDXSEL]);

						Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL] / (*ophoMIPNhitCone)[mIDXSEL]);
					}
				}

				
			}

			//prompt W template
			if ((*ophohasPixelSeed)[mIDXSEL] == 1 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
			{
				

				if ((*ophoMIPTotEnergy)[mIDXSEL] > 4.9)
				{
					promptWTemp->Fill((*ophoSeedTime)[mIDXSEL]);
					MIPW->Fill((*ophoMIPTotEnergy)[mIDXSEL]);
					Float_t TransMW = TransMassW((*ophoEt)[mIDXSEL], (*ophoPhi)[mIDXSEL], (*ophoEta)[mIDXSEL]);
					TransMass_W->Fill(TransMW);

					//check MIP variables
					MIPTotalE_Wprompt->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

					MIPChi2_Wprompt->Fill((*ophoMIPChi2)[mIDXSEL]);

					MIPSlope_Wprompt->Fill((*ophoMIPSlope)[mIDXSEL]);

					MIPIntercept_Wprompt->Fill((*ophoMIPIntercept)[mIDXSEL]);

					MIPNhitCone_Wprompt->Fill((*ophoMIPNhitCone)[mIDXSEL]);

					MIPisHalo_Wprompt->Fill((*ophoMIPIsHalo)[mIDXSEL]);

					MIPCHi2_over_MIPNhitCone_Wprompt->Fill((*ophoMIPChi2)[mIDXSEL] / (*ophoMIPNhitCone)[mIDXSEL]);
				}
			}
		}

		if (motp.size() > 0 && mitp.size() > 0)
		{
			mIDXSEL = motp[0];
			misOOT = kTRUE;
			MFullTime->Fill((*ophoSeedTime)[mIDXSEL]);

			//Candidate Events
			if ((*ophohasPixelSeed)[mIDXSEL] == 0 && (*ophoMIPTotEnergy)[mIDXSEL] < 4.9 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.01015 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
			{
				Candidate_Events->Fill((*ophoSeedTime)[mIDXSEL]);

				metfilter_candidate->Fill(metFilters);
			}

			//Spike Template
			if ((*ophohasPixelSeed)[mIDXSEL] == 0 && (*ophoMIPTotEnergy)[mIDXSEL] < 4.9 && ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.001 || (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] < 0.001))
			{
				SpikeTemp->Fill((*ophoSeedTime)[mIDXSEL]);
			}

			//Beam Halo Template
			if ((*ophohasPixelSeed)[mIDXSEL] == 0 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
			{
				woMIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

				woMIPChi2_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL]);

				woMIPSlope_BeamHalo->Fill((*ophoMIPSlope)[mIDXSEL]);

				woMIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[mIDXSEL]);
				
				woMIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[mIDXSEL]);

				woMIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[mIDXSEL]);

				wophi_BeamHalo->Fill((*ophoPhi)[mIDXSEL]);

				Float_t MIPChi2_over_MIPNhitCone = ((*ophoMIPChi2)[mIDXSEL] / (*ophoMIPNhitCone)[mIDXSEL]);
				if ((*ophoMIPTotEnergy)[mIDXSEL] > 4.9 && MIPChi2_over_MIPNhitCone < 0.2)
				{
					BeamHaloTemp->Fill((*ophoSeedTime)[mIDXSEL]);

					metfilter_beamhalo->Fill(metFilters);

					BeamHaloSieie_vs_time->Fill((*ophoSeedTime)[mIDXSEL], (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]);

					if ((*ophoSeedTime)[mIDXSEL] > -2.0 && (*ophoSeedTime)[mIDXSEL] < 2.0 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.0 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.015)
					{
						BeamHalo_prompt_region->Fill((*ophoSeedTime)[mIDXSEL], (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]);
						file2 << run << " " << event << " " << lumis << endl;
					}

					BeamHaloSieie->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]);

					MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

					MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL]);

					MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[mIDXSEL]);

					MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[mIDXSEL]);

					MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[mIDXSEL]);

					MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[mIDXSEL]);

					phi_BeamHalo->Fill((*ophoPhi)[mIDXSEL]);

					//check MIP variables
					if ((*ophoPhi)[mIDXSEL] > -3.2 && (*ophoPhi)[mIDXSEL] < -2.8)
					{
						Phi_MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

						Phi_MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL]);

						Phi_MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[mIDXSEL]);

						Phi_MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[mIDXSEL]);

						Phi_MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[mIDXSEL]);

						Phi_MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[mIDXSEL]);

						Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL] / (*ophoMIPNhitCone)[mIDXSEL]);
					}

					if ((*ophoPhi)[mIDXSEL] > -0.5 && (*ophoPhi)[mIDXSEL] < 0.5)
					{
						Phi_MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

						Phi_MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL]);

						Phi_MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[mIDXSEL]);

						Phi_MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[mIDXSEL]);

						Phi_MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[mIDXSEL]);

						Phi_MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[mIDXSEL]);

						Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL] / (*ophoMIPNhitCone)[mIDXSEL]);
					}

					if ((*ophoPhi)[mIDXSEL] > 2.8 && (*ophoPhi)[mIDXSEL] < 3.2)
					{
						Phi_MIPTotalE_BeamHalo->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

						Phi_MIPChi2_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL]);

						Phi_MIPSlope_BeamHalo->Fill((*ophoMIPSlope)[mIDXSEL]);

						Phi_MIPIntercept_BeamHalo->Fill((*ophoMIPIntercept)[mIDXSEL]);

						Phi_MIPNhitCone_BeamHalo->Fill((*ophoMIPNhitCone)[mIDXSEL]);

						Phi_MIPisHalo_BeamHalo->Fill((*ophoMIPIsHalo)[mIDXSEL]);

						Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Fill((*ophoMIPChi2)[mIDXSEL] / (*ophoMIPNhitCone)[mIDXSEL]);
					}
				}
				
			}

			//prompt W template
			if ((*ophohasPixelSeed)[mIDXSEL] == 1 && (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 && (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001)
			{


				if ((*ophoMIPTotEnergy)[mIDXSEL] > 4.9)
				{
					promptWTemp->Fill((*ophoSeedTime)[mIDXSEL]);
					MIPW->Fill((*ophoMIPTotEnergy)[mIDXSEL]);
					Float_t TransMW = TransMassW((*ophoEt)[mIDXSEL], (*ophoPhi)[mIDXSEL], (*ophoEta)[mIDXSEL]);
					TransMass_W->Fill(TransMW);

					//check MIP variables
					MIPTotalE_Wprompt->Fill((*ophoMIPTotEnergy)[mIDXSEL]);

					MIPChi2_Wprompt->Fill((*ophoMIPChi2)[mIDXSEL]);

					MIPSlope_Wprompt->Fill((*ophoMIPSlope)[mIDXSEL]);

					MIPIntercept_Wprompt->Fill((*ophoMIPIntercept)[mIDXSEL]);

					MIPNhitCone_Wprompt->Fill((*ophoMIPNhitCone)[mIDXSEL]);

					MIPisHalo_Wprompt->Fill((*ophoMIPIsHalo)[mIDXSEL]);

					MIPCHi2_over_MIPNhitCone_Wprompt->Fill((*ophoMIPChi2)[mIDXSEL] / (*ophoMIPNhitCone)[mIDXSEL]);
				}
			}
		}
		//-------------Done timing distribution with MET, Candidate Events & W prompt, spike, beam halo template--------------------
	}

	//-----------------full timing distribution---------------------

	FullTime->SetLineColor(1);
	FullTime->Write();

	MFullTime->SetLineColor(1);
	MFullTime->Write();

	//--------------------------------------------------------------


	//----------------------Templates-----------------------

	//Candidate events
	Candidate_Events->SetLineColor(1);
	Candidate_Events->Write();

	//W Prompt Template
	promptWTemp->SetLineColor(1);
	promptWTemp->Write();

	//MIP total Energy from W Prompt
	MIPW->SetLineColor(1);
	MIPW->Write();

	//Invariant Mass of W
	TransMass_W->SetLineColor(1);
	TransMass_W->Write();

	//Z Prompt Template
	promptZTemp->SetLineColor(1);
	promptZTemp->Write();

	//MIP total Energy from Z Prompt
	MIPZ->SetLineColor(1);
	MIPZ->Write();

	//Invariant Mass of Z
	InvMass_Z->SetLineColor(1);
	InvMass_Z->Write();

	//Spike Template
	SpikeTemp->SetLineColor(1);
	SpikeTemp->Write();

	//Beam Halo Template
	BeamHaloTemp->SetLineColor(1);
	BeamHaloTemp->Write();

	//Check shower shape vs. time
	BeamHaloSieie_vs_time->SetOption("COLZ");
	BeamHaloSieie_vs_time->Write();

	//Pick the event
	BeamHalo_prompt_region->SetOption("COLZ");
	BeamHalo_prompt_region->Write();

	//Check shower shape
	BeamHaloSieie->SetLineColor(1);
	BeamHaloSieie->Write();

	//MIP variables
	MIPTotalE_BeamHalo->SetLineColor(1);
	MIPTotalE_BeamHalo->Write();

	MIPChi2_BeamHalo->SetLineColor(1);
	MIPChi2_BeamHalo->Write();

	MIPSlope_BeamHalo->SetLineColor(1);
	MIPSlope_BeamHalo->Write();

	MIPIntercept_BeamHalo->SetLineColor(1);
	MIPIntercept_BeamHalo->Write();

	MIPNhitCone_BeamHalo->SetLineColor(1);
	MIPNhitCone_BeamHalo->Write();

	MIPisHalo_BeamHalo->SetLineColor(1);
	MIPisHalo_BeamHalo->Write();

	metfilter_candidate->SetLineColor(1);
	metfilter_candidate->Write();

	metfilter_beamhalo->SetLineColor(1);
	metfilter_beamhalo->Write();
	//--------------------------------------------------------

	//MIP variable w/o totE cut
	woMIPTotalE_BeamHalo->SetLineColor(1);
	woMIPTotalE_BeamHalo->Write();

	woMIPChi2_BeamHalo->SetLineColor(1);
	woMIPChi2_BeamHalo->Write();

	woMIPSlope_BeamHalo->SetLineColor(1);
	woMIPSlope_BeamHalo->Write();

	woMIPIntercept_BeamHalo->SetLineColor(1);
	woMIPIntercept_BeamHalo->Write();

	woMIPNhitCone_BeamHalo->SetLineColor(1);
	woMIPNhitCone_BeamHalo->Write();

	woMIPisHalo_BeamHalo->SetLineColor(1);
	woMIPisHalo_BeamHalo->Write();

	wophi_BeamHalo->SetLineColor(1);
	wophi_BeamHalo->Write();

	//------------MIP variable in the poputlated region Beam halo-----------

	phi_BeamHalo->SetLineColor(1);
	phi_BeamHalo->Write();

	Phi_MIPTotalE_BeamHalo->SetLineColor(1);
	Phi_MIPTotalE_BeamHalo->Write();

	Phi_MIPChi2_BeamHalo->SetLineColor(1);
	Phi_MIPChi2_BeamHalo->Write();

	Phi_MIPSlope_BeamHalo->SetLineColor(1);
	Phi_MIPSlope_BeamHalo->Write();

	Phi_MIPIntercept_BeamHalo->SetLineColor(1);
	Phi_MIPIntercept_BeamHalo->Write();

	Phi_MIPNhitCone_BeamHalo->SetLineColor(1);
	Phi_MIPNhitCone_BeamHalo->Write();

	Phi_MIPisHalo_BeamHalo->SetLineColor(1);
	Phi_MIPisHalo_BeamHalo->Write();

	Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->SetLineColor(1);
	Phi_MIPCHi2_over_MIPNhitCone_BeamHalo->Write();
	//----------------------------------------------------------

	//-------------MIP variables from W prompt---------------
	MIPTotalE_Wprompt->SetLineColor(1);
	MIPTotalE_Wprompt->Write();

	MIPChi2_Wprompt->SetLineColor(1);
	MIPChi2_Wprompt->Write();

	MIPSlope_Wprompt->SetLineColor(1);
	MIPSlope_Wprompt->Write();

	MIPIntercept_Wprompt->SetLineColor(1);
	MIPIntercept_Wprompt->Write();

	MIPNhitCone_Wprompt->SetLineColor(1);
	MIPNhitCone_Wprompt->Write();

	MIPisHalo_Wprompt->SetLineColor(1);
	MIPisHalo_Wprompt->Write();

	MIPCHi2_over_MIPNhitCone_Wprompt->SetLineColor(1);
	MIPCHi2_over_MIPNhitCone_Wprompt->Write();
	//--------------------------------------------------------

	//-------------MIP variables from Z prompt---------------
	MIPTotalE_Zprompt->SetLineColor(1);
	MIPTotalE_Zprompt->Write();

	MIPChi2_Zprompt->SetLineColor(1);
	MIPChi2_Zprompt->Write();

	MIPSlope_Zprompt->SetLineColor(1);
	MIPSlope_Zprompt->Write();

	MIPIntercept_Zprompt->SetLineColor(1);
	MIPIntercept_Zprompt->Write();

	MIPNhitCone_Zprompt->SetLineColor(1);
	MIPNhitCone_Zprompt->Write();

	MIPisHalo_Zprompt->SetLineColor(1);
	MIPisHalo_Zprompt->Write();
	//--------------------------------------------------------






	cout << "Done filling histograms" << endl;

	/*
	RunZ.close();
	RunW.close();
	RunBeamhalo.close();
	*/

	cout << "Analysis completed" << endl;
	auto timenow_end = chrono::system_clock::to_time_t(chrono::system_clock::now());

	cout << "start time: " << ctime(&timenow_start) << " (CT)" << endl;

	cout << "End time: " << ctime(&timenow_end) << " (CT)" << endl;

	file1 << "Analysis completed!" << endl;
	file1 << "Completed time: " << ctime(&timenow_end) << endl;

	file1.close();
	file2.close();
}



