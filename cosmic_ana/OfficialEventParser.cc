#define __OFFICIALEVENTPARSER_CXX

#include "OfficialEventParser.h"
#include "Event.h"
#include <string>
#include <math.h>
#include <cmath>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <map>
#include <TMath.h>
//#include "pdg_codes.h"

#define __FC_RUN_SELECTION__  true 

using namespace std;
using namespace global;

OfficialEventParser::OfficialEventParser(){
  E = new Event();
  set_version(-1); // initialize
  set_nrgSepType(0);
  set_maxE(999999);
  set_osc(false); // not oscillate
  set_datatype(fcData); // fcdata, fcpc, etc. written in global.cc
  
  for (int ireco=0; ireco<nRecos; ireco++)  {
    for (int imc=0; imc<nSample; imc++)  {
        nentries_norm[imc][ireco] = 0;
	efile[imc][ireco] = 0;
    }
  }
}

void OfficialEventParser::InitEvent(){
  set_from_wall_flag(false);
}

TH1D* OfficialEventParser::MakeNewHisto(const char* name, const char* title, Int_t nbinsx, const Double_t *xbins, vector<TH1 *>& h_vect) {
  
  TH1D *h = new TH1D(name, title, nbinsx, xbins);
  
  h->Sumw2();
  h_vect.push_back(h);
  
  return h;
}

TH1D* OfficialEventParser::MakeNewHisto(const char* name, const char* title, Int_t nbinsx, const Double_t xlow, const Double_t xhigh, vector<TH1 *>& h_vect) {
  
  TH1D *h = new TH1D(name, title, nbinsx, xlow, xhigh);
  
  h->Sumw2();
  h_vect.push_back(h);
  
  return h;
}

TH2D* OfficialEventParser::MakeNewHisto(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, vector<TH1 *>&h_vect) {
  
  TH2D *h = new TH2D(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);
  
  h->Sumw2();
  h_vect.push_back(h);
  
  return h;
}

TH3D* OfficialEventParser::MakeNewHisto(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup, vector<TH1 *>&h_vect) {

  TH3D *h = new TH3D(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup);
  h->Sumw2();
  h_vect.push_back(h);

  return h;

}




void OfficialEventParser::InitHistograms(int ireco, int imc){
  // Fill histograms in a root file.
  string outputnameroot="";
  
  //for (int ireco=0; ireco<nRecos; ireco++) 
  {
	
    char *h_name[2] = {"h","hmc"};
    //for (int imc=0; imc<nSample; imc++) 
    {

      string outputname=outputnameroot+Form("%s_reco%d_nrg%d_Elt%d_file%d.root",h_name[imc],ireco,_nrgSepType,_maxE,_fileNum);
      efile[imc][ireco] = new TFile(outputname.c_str(),"recreate");

      TString histName;
      
      for (int ise=0; ise<nSE; ise++) {
	
	histName = Form("%s_cosmic_ise%d_mom_binning_reco%d",h_name[imc],ise,ireco);
	
	if (ise==ise_muon)
	  h_cosmic_mom_binning[imc][ireco][ise] = MakeNewHisto(histName.Data(), Form("%s; Momentum (MeV/c)",histName.Data()), nMuonMomBinsLike, muonMomBinsLike, h_vect_hists[imc][ireco]);
	else if (ise==ise_dcye)
	  h_cosmic_mom_binning[imc][ireco][ise] = MakeNewHisto(histName.Data(), Form("%s; Momentum (MeV/c)",histName.Data()), nDcyeMomBinsLike, dcyeMomBinsLike, h_vect_hists[imc][ireco]);
      }

      //h_cosmic_emom_vs_year[imc][ireco] = new TH2D(Form("%s_cosmic_emom_vs_year_reco%d",h_name[imc],ireco), Form("%s_cosmic_emom_vs_year_reco%d; Year; Momentum (MeV/c)",h_name[imc],ireco), nYearBins, yearRange[0], yearRange[1], nMomBins[ise_dcye], momLow[ise_dcye], momHigh[ise_dcye]);
      //h_vect_hists[imc][ireco].push_back(h_cosmic_emom_vs_year[imc][ireco]);
      //h_cosmic_emom_vs_year[imc][ireco] ->Sumw2();
      

      histName = Form("%s_cosmic_emom_vs_year_reco%d",h_name[imc],ireco);
      h_cosmic_emom_vs_year[imc][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Momentum (MeV/c)",histName.Data()), nYearBins, yearRange[0], yearRange[1], nMomBins[ise_dcye], momLow[ise_dcye], momHigh[ise_dcye], h_vect_hists[imc][ireco]);

      histName = Form("%s_cosmic_mumom_vs_year_reco%d",h_name[imc],ireco);
      h_cosmic_mumom_vs_year[imc][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Momentum (MeV/c)",histName.Data()), nYearBins, yearRange[0], yearRange[1], nMomBins[ise_muon], momLow[ise_muon], momHigh[ise_muon], h_vect_hists[imc][ireco]);
      
      histName = Form("%s_cosmic_lifetime_vs_year_reco%d",h_name[imc],ireco);
      h_cosmic_lifetime_vs_year[imc][ireco] = MakeNewHisto(histName.Data(),Form("%s; Year; Decay-e Time, #Deltat (#mus)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 250, 0, 25, h_vect_hists[imc][ireco]); 
      
      histName = Form("%s_cosmic_emom_vs_lifetime_vs_year_reco%d",h_name[imc],ireco);
      h_cosmic_emom_vs_lifetime_vs_year[imc][ireco]     = MakeNewHisto(histName.Data(), Form("%s; Year; Decay-e Time, #Deltat (#mus); Momentum (MeV/c)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 250, 0, 25, nMomBins[ise_dcye], momLow[ise_dcye], momHigh[ise_dcye], h_vect_hists[imc][ireco]);

      histName = Form("%s_cosmic_emom_vs_lifetime_vs_year_reco%d; Year; Decay-e Time, #Deltat (#mus); Momentum (MeV/c)",h_name[imc],ireco);
      h_cosmic_range_vs_mumom_vs_year[imc][ireco]     = MakeNewHisto(histName.Data(), Form("%s;Year; Momentum (GeV/c);Range (m)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 60, 0, 12, 100, 0, 50, h_vect_hists[imc][ireco]);
      
      
      for (int ise=0; ise<nSE; ise++) {
	
	for (int icut=0; icut<nCuts; icut++) {
	  histName = Form("%s_cosmic_cut%d_ise%d_z_vs_year_reco%d",h_name[imc],icut,ise,ireco);
	  h_cosmic_z_vs_year[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s;Year; z (cm)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, -2000, 2000, h_vect_hists[imc][ireco]);
	  
	  histName = Form("%s_cosmic_cut%d_ise%d_r2_vs_year_reco%d",h_name[imc],icut,ise,ireco);
	  h_cosmic_r2_vs_year[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(),  Form("%s;Year; R^{2} (cm^{2})",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0000000, 4000000, h_vect_hists[imc][ireco]);
	  //h_cosmic_r2_vs_year[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(),  Form("%s;Year; R^{2} (cm^{2})",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, 100000, h_vect_hists[imc][ireco]);

	  histName = Form("%s_cosmic_cut%d_ise%d_zzoom_vs_year_reco%d",h_name[imc],icut,ise,ireco);
	  h_cosmic_zzoom_vs_year[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s;Year; z (cm)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 1400, 2200, h_vect_hists[imc][ireco]);
	  
	  histName = Form("%s_cosmic_cut%d_ise%d_r2zoom_vs_year_reco%d",h_name[imc],icut,ise,ireco);
	  h_cosmic_r2zoom_vs_year[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(),  Form("%s;Year; R^{2} (cm^{2})",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 1400000, 4200000, h_vect_hists[imc][ireco]);

	  histName = Form("%s_cosmic_cut%d_ise%d_z_nearwall_reco%d",h_name[imc],icut,ise,ireco);
	  h_cosmic_z_nearwall[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s;z (cm)",histName.Data()), 100, 1650, 1900, h_vect_hists[imc][ireco]);
	  histName = Form("%s_cosmic_cut%d_ise%d_r2_nearwall_reco%d",h_name[imc],icut,ise,ireco);
	  h_cosmic_r2_nearwall[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s;R^{2} (cm^{2})",histName.Data()), 100, 2400000, 3400000, h_vect_hists[imc][ireco]);
	  
	}

	histName = Form("%s_cosmic_ise%d_wall_vs_year_reco%d",h_name[imc],ise,ireco);
	if (ise==0) h_cosmic_wall_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(),  Form("%s;Year; Distance from Nearest Wall (cm)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 120, -600, 600, h_vect_hists[imc][ireco]);
	else h_cosmic_wall_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(),  Form("%s;Year; Distance from Nearest Wall (cm)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 120, -600, 1800, h_vect_hists[imc][ireco]);
	
	histName = Form("%s_cosmic_ise%d_vtx_prfvtx_diff_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_vtx_prfvtx_diff_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; |x_{full}-x_{prefit}| (cm)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, 4000, h_vect_hists[imc][ireco]);
	
	histName = Form("%s_cosmic_ise%d_pcflg_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_pcflg_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; PC Flag",histName.Data()), nYearBins, yearRange[0], yearRange[1], 5, -2.5, 2.5, h_vect_hists[imc][ireco]);
	
	
      }

      histName = Form("%s_cosmic_nmue_vs_year_reco%d",h_name[imc],ireco); 
      h_cosmic_nmue_vs_year[imc][ireco] = MakeNewHisto(histName.Data(),  Form("%s;Year; Number of Decay-e",histName.Data()), nYearBins, yearRange[0], yearRange[1], 7, -0.5, 6.5, h_vect_hists[imc][ireco]);

      histName = Form("%s_cosmic_ingate_vs_year_reco%d",h_name[imc],ireco); 
      h_cosmic_ingate_vs_year[imc][ireco] = MakeNewHisto(histName.Data(),  Form("%s;Year; In-gate",histName.Data()), nYearBins, yearRange[0], yearRange[1], 5, -0.5, 4.5, h_vect_hists[imc][ireco]);
      
      for (int icut=0; icut<nRangeCuts; icut++) {
	for (int irange=0; irange<nRanges; irange++) {
	
	  double rangeLowEdge = (rangeLow+irange*rangeStep);
	  double rangeHighEdge = (rangeLow+(irange+1)*rangeStep);
	
	  double MomOverRangeHigh = 3.5;
	  if (irange==0) MomOverRangeHigh = 5;

	  histName = Form("%s_cosmic_mom_over_range_vs_year_reco%d_rangesep%03d_%03d_icut%d",h_name[imc],ireco,(int)rangeLowEdge,(int)rangeHighEdge,icut);
	  h_cosmic_mom_over_range_rangesep_vs_year[imc][ireco][irange][icut] = MakeNewHisto(histName.Data(), Form("%s; Year; p/Range (MeV/c/cm)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 40, 1.5, MomOverRangeHigh, h_vect_hists[imc][ireco]);
	}
      }



      histName = Form("%s_cosmic_n50_vs_year_reco%d",h_name[imc],ireco);
      h_cosmic_n50_vs_year[imc][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year;  n50 (# of Hits)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, 500, h_vect_hists[imc][ireco]);

      histName = Form("%s_cosmic_q50_vs_year_reco%d",h_name[imc],ireco);
      h_cosmic_q50_vs_year[imc][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; q50 (p.e.)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, 800, h_vect_hists[imc][ireco]);
      
      double minLength = 400;  // ns
      double maxHits[nSE] = {12000, 600};
      double maxCharge[nSE] = {3e5, 1200};

      for (int ise=0; ise<nSE; ise++) {
	
	histName = Form("%s_cosmic_ise%d_clnhits_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_clnhits_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Hits in spliTChan Cluster",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, maxHits[ise], h_vect_hists[imc][ireco]);

	histName = Form("%s_cosmic_ise%d_cltotq_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_cltotq_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Charge in spliTChan Cluster",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, maxCharge[ise], h_vect_hists[imc][ireco]);


	histName = Form("%s_cosmic_ise%d_cltstart_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_cltstart_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Cluster Start Time (#mus)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, 30, h_vect_hists[imc][ireco]);


	histName = Form("%s_cosmic_ise%d_cltend_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_cltend_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Cluster End Time (#mus)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, 30, h_vect_hists[imc][ireco]);

	histName = Form("%s_cosmic_ise%d_cllength_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_cllength_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Cluster Length (#mus)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, 2, h_vect_hists[imc][ireco]);

	histName = Form("%s_cosmic_ise%d_cllength_vs_cltstart_reco%d",h_name[imc],ise,ireco);
	h_cosmic_cllength_vs_cltstart[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Cluster Start (#mus); Cluster Length (#mus)",histName.Data()), 100, 0, 30, 100, 0, 2, h_vect_hists[imc][ireco]);


	histName = Form("%s_cosmic_ise%d_cltotqnorm_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_cltotqnorm_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; spliTChan Cluster Charge Rate (pe / ns)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, maxCharge[ise]/minLength, h_vect_hists[imc][ireco]);


	histName = Form("%s_cosmic_ise%d_clnhitsnorm_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_clnhitsnorm_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; spliTChan Cluster Hit Rate (Hits / ns)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, maxHits[ise]/minLength, h_vect_hists[imc][ireco]);


	histName = Form("%s_cosmic_ise%d_fqtotqnorm_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_fqtotqnorm_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Charge Rate (pe / ns)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, maxCharge[ise]/minLength, h_vect_hists[imc][ireco]);


	histName = Form("%s_cosmic_ise%d_fqnhitsnorm_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_fqnhitsnorm_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Hit Rate (Hits / ns)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, maxHits[ise]/minLength, h_vect_hists[imc][ireco]);


	histName = Form("%s_cosmic_ise%d_fqnhits_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_fqnhits_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Hits in fiTQun Sub-Event",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, maxHits[ise], h_vect_hists[imc][ireco]);

	histName = Form("%s_cosmic_ise%d_fqtotq_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_fqtotq_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Charge in fiTQun Sub-Event",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, maxCharge[ise], h_vect_hists[imc][ireco]);

	histName = Form("%s_cosmic_ise%d_fq1rtotmu_vs_year_reco%d",h_name[imc],ise,ireco);
	h_cosmic_fq1rtotmu_vs_year[imc][ise][ireco] = MakeNewHisto(histName.Data(), Form("%s; Year; Predicted Charge (p.e.)",histName.Data()), nYearBins, yearRange[0], yearRange[1], 100, 0, maxCharge[ise], h_vect_hists[imc][ireco]);

	float lnLrat1R_low[nRecos][nSE] = {
	  -1, 0, //apfit
	  //-25000, -500,
	  -12000, -500
	};
  
	float lnLrat1R_high[nRecos][nSE] = {
	  1, 1, //apfit
	  //2000, 1000,
	  2000, 700
	};
	
	histName = Form("%s_cosmic_ise%d_lnLrat1R_vs_mom_reco%d",h_name[imc],ise,ireco);
	if (ise==0) {
	  h_cosmic_lnLrat1R_vs_mom[imc][ise][ireco]     = MakeNewHisto(histName.Data(), Form("%s;p_{e} (MeV/c); ln(L_{e}/L_{#mu})",histName.Data()), 250, 0, min(_maxE,5000),  (int)((lnLrat1R_high[ireco][ise]-lnLrat1R_low[ireco][ise])/20), lnLrat1R_low[ireco][ise], lnLrat1R_high[ireco][ise], h_vect_hists[imc][ireco]);
	} else {
	  h_cosmic_lnLrat1R_vs_mom[imc][ise][ireco]     = MakeNewHisto(histName.Data(), Form("%s;p_{e} (MeV/c); ln(L_{e}/L_{#mu})",histName.Data()), 100, 0, 100, 150, lnLrat1R_low[ireco][ise], lnLrat1R_high[ireco][ise], h_vect_hists[imc][ireco]);
	}
		
	histName = Form("%s_cosmic_ise%d_lnLrat2R1R_vs_mom_reco%d",h_name[imc],ise,ireco);
	if (ise==0) {
	  h_cosmic_lnLrat2R1R_vs_mom[imc][ise][ireco]     = MakeNewHisto(histName.Data(), Form("%s;p_{e} (MeV/c); ln(L_{2R}/L_{1R})",histName.Data()), 250, 0, 10000,  200, -20000, 20000, h_vect_hists[imc][ireco]);
	} else {
	  h_cosmic_lnLrat2R1R_vs_mom[imc][ise][ireco]     = MakeNewHisto(histName.Data(), Form("%s;p_{e} (MeV/c); ln(L_{2R}/L_{1R})",histName.Data()), 100, 0, 2000,  100, -3000, 3000, h_vect_hists[imc][ireco]);
	}

      }

      
      //double yearBins[nYearBins];
      //for (int idate=0; idate<nYearBins; idate++)
      //yearBins[idate] = yearRange[0] + idate*(yearRange[1]-yearRange[0])/nYearBins;
      
      // Muon sub-event likelihoods
      for (int ilike=0; ilike<nMuonLikes; ilike++) { 

	double slopeTop = (muonLikeTopRight[ilike]-muonLikeTopLeft[ilike])/(muonMomBinsLike[nMuonMomBinsLike]-muonMomBinsLike[0]);
	double slopeBot = (muonLikeBotRight[ilike]-muonLikeBotLeft[ilike])/(muonMomBinsLike[nMuonMomBinsLike]-muonMomBinsLike[0]);
	
	for (int imom=0; imom<nMuonMomBinsLike; imom++) {
	  
	  double topLike = slopeTop*muonMomBinsLike[imom] + muonLikeTopLeft[ilike];
	  double botLike = slopeBot*muonMomBinsLike[imom] + muonLikeBotLeft[ilike];
	  
	  // Put 0 on a bin edge
	  if (botLike<0 && topLike>0) {
	    int nBinsLeftSide = fabs(botLike)/(topLike-botLike)*nLikeBins;
	    double binWidth = fabs(botLike)/nBinsLeftSide;
	    topLike = botLike+binWidth*nLikeBins;
	  }

	  //double likeBins[nLikeBins];
	  //for (int ilike=0; ilike<nLikeBins; ilike++) 
	  //likeBins[ilike] = botLike + ilike*(topLike-botLike)/nLikeBins;
	  
	  histName = Form("%s_cosmic_ise0_%s_vs_year_reco%d_momsep%d_%d",h_name[imc],muonLikeString[ilike].Data(),ireco,(int)muonMomBinsLike[imom],(int)muonMomBinsLike[imom+1]);
	  h_cosmic_Lmuon_vs_year[imc][ireco][ilike][imom] = MakeNewHisto(histName.Data(), Form("%s; Year; %s",histName.Data(), muonLikeTitle[ilike].Data()), nYearBins, yearRange[0], yearRange[1], nLikeBins, botLike, topLike, h_vect_hists[imc][ireco]);	  


	  histName = Form("%s_cosmic_ise0_%s_vs_zenith_reco%d_momsep%d_%d",h_name[imc],muonLikeString[ilike].Data(),ireco,(int)muonMomBinsLike[imom],(int)muonMomBinsLike[imom+1]);
	  h_cosmic_Lmuon_vs_zenith[imc][ireco][ilike][imom] = MakeNewHisto(histName.Data(), Form("%s; cos(#theta_{zenith}); %s",histName.Data(), muonLikeTitle[ilike].Data()), 100, -1, 1, nLikeBins, botLike, topLike, h_vect_hists[imc][ireco]);	  
	  
	}
	
      }
      
      // Electron sub-event likelihoods
      for (int ilike=0; ilike<nDcyeLikes; ilike++) { 

	double slopeTop = (dcyeLikeTopRight[ilike]-dcyeLikeTopLeft[ilike])/(dcyeMomBinsLike[nDcyeMomBinsLike]-dcyeMomBinsLike[0]);
	double slopeBot = (dcyeLikeBotRight[ilike]-dcyeLikeBotLeft[ilike])/(dcyeMomBinsLike[nDcyeMomBinsLike]-dcyeMomBinsLike[0]);
	
	for (int imom=0; imom<nDcyeMomBinsLike; imom++) {
	  
	  double topLike = slopeTop*dcyeMomBinsLike[imom] + dcyeLikeTopLeft[ilike];
	  double botLike = slopeBot*dcyeMomBinsLike[imom] + dcyeLikeBotLeft[ilike];

	  // Put 0 on a bin edge
	  if (botLike<0 && topLike>0) {
	    int nBinsLeftSide = fabs(botLike)/(topLike-botLike)*nLikeBins;
	    double binWidth = fabs(botLike)/nBinsLeftSide;
	    topLike = botLike+binWidth*nLikeBins;
	  }

	  //double likeBins[nLikeBins];
	  //for (int ilike=0; ilike<nLikeBins; ilike++) 
	  //likeBins[ilike] = botLike + ilike*(topLike-botLike)/nLikeBins;
	    
	  histName = Form("%s_cosmic_ise1_%s_vs_year_reco%d_momsep%d_%d",h_name[imc],dcyeLikeString[ilike].Data(),ireco,(int)dcyeMomBinsLike[imom],(int)dcyeMomBinsLike[imom+1]);
	  h_cosmic_Ldcye_vs_year[imc][ireco][ilike][imom] = MakeNewHisto(histName.Data(), Form("%s; Year; %s",histName.Data(), dcyeLikeTitle[ilike].Data()), nYearBins, yearRange[0], yearRange[1], nLikeBins, botLike, topLike, h_vect_hists[imc][ireco]);	  

	  histName = Form("%s_cosmic_ise1_%s_vs_zenith_reco%d_momsep%d_%d",h_name[imc],dcyeLikeString[ilike].Data(),ireco,(int)dcyeMomBinsLike[imom],(int)dcyeMomBinsLike[imom+1]);
	  h_cosmic_Ldcye_vs_zenith[imc][ireco][ilike][imom] = MakeNewHisto(histName.Data(), Form("%s; cos(#theta_{zenith}); %s",histName.Data(), dcyeLikeTitle[ilike].Data()), 100, -1, 1, nLikeBins, botLike, topLike, h_vect_hists[imc][ireco]);	  
	  
	}
	
      }

      /*
      for (int imom=0; imom<nMuonMomBinsLike; imom++) {
	  
	histName = Form("%s_cosmic_ise0_muonmisid_vs_year_reco%d_momsep%d_%d",h_name[imc],ireco,(int)muonMomBinsLike[imom],(int)muonMomBinsLike[imom+1]);
	h_cosmic_Lmuon_vs_year[imc][ireco][ilike][imom] = MakeNewHisto(histName.Data(), Form("%s; Year; %s",histName.Data(), muonLikeTitle[ilike].Data()), nYearBins, yearRange[0], yearRange[1], 50, 0, 5, h_vect_hists[imc][ireco]);	  	
      }

      for (int imom=0; imom<nDcyeMomBinsLike; imom++) {
	  	    
	histName = Form("%s_cosmic_ise1_dcyemisid_vs_year_reco%d_momsep%d_%d",h_name[imc],ireco,(int)dcyeMomBinsLike[imom],(int)dcyeMomBinsLike[imom+1]);
	h_cosmic_Ldcye_vs_year[imc][ireco][ilike][imom] = MakeNewHisto(histName.Data(), Form("%s; Year; %s",histName.Data(), dcyeLikeTitle[ilike].Data()), nYearBins, yearRange[0], yearRange[1], 50, 0, 5, h_vect_hists[imc][ireco]);	  
      }
      /**/

      double dirResMax[nSE] = {15, 60};
      double resLow[nSE] = {-30, -100};
      double resHigh[nSE] = {30, 100};

      for (int icut=0; icut<nCuts; icut++) {
	for (int ise=0; ise<nSE; ise++) {
	  
	  histName = Form("%s_cosmic_ise%d_mom_res_reco%d_cut%d",h_name[imc],ise,ireco,icut);
	  h_cosmic_mom_res[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s; (p_{Rec.}-p_{True})/p_{True} (%%)",histName.Data()), 200, resLow[ise], resHigh[ise], h_vect_hists[imc][ireco]);

	  histName = Form("%s_cosmic_ise%d_mom_res_vs_trumom_reco%d_cut%d",h_name[imc],ise,ireco,icut);
	  h_cosmic_mom_res_vs_trumom[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s; p_{True} (MeV/c); (p_{Rec.}-p_{True})/p_{True} (%%)",histName.Data()), nMomBins[ise], momLow[ise], momHigh[ise], 200, resLow[ise], resHigh[ise], h_vect_hists[imc][ireco]);
	  
	  histName = Form("%s_cosmic_ise%d_vtx_res_reco%d_cut%d",h_name[imc],ise,ireco,icut);
	  h_cosmic_vtx_res[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s; |Vtx_{Rec.} - Vtx_{True}| (cm)",histName.Data()), 200, 0, 200, h_vect_hists[imc][ireco]);

	  histName = Form("%s_cosmic_ise%d_tvtx_res_reco%d_cut%d",h_name[imc],ise,ireco,icut);
	  h_cosmic_tvtx_res[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s; |Vtx_{#perp, Rec.} - Vtx_{#perp, True}| (cm)",histName.Data()), 150, 0, 150, h_vect_hists[imc][ireco]);

	  histName = Form("%s_cosmic_ise%d_lvtx_res_reco%d_cut%d",h_name[imc],ise,ireco,icut);
	  h_cosmic_lvtx_res[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s; Vtx_{#parallel, Rec.} - Vtx_{#parallel, True} (cm)",histName.Data()), 300, -150, 150, h_vect_hists[imc][ireco]);
	  
	  histName = Form("%s_cosmic_ise%d_dir_res_reco%d_cut%d",h_name[imc],ise,ireco,icut);
	  h_cosmic_dir_res[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s; |#theta_{Rec.} - #theta_{True}| (#circ)",histName.Data()), 100, 0, dirResMax[ise], h_vect_hists[imc][ireco]);
	  
	  histName = Form("%s_cosmic_ise%d_time_res_reco%d_cut%d",h_name[imc],ise,ireco,icut);
	  h_cosmic_time_res[imc][ise][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s; #Delta{t}_{Rec.} - #Delta{t}_{True} (ns)",histName.Data()), 100, -500, 500, h_vect_hists[imc][ireco]);

	}
	
	histName = Form("%s_cosmic_ddir_res_reco%d_cut%d",h_name[imc],ireco,icut);
	h_cosmic_ddir_res[imc][ireco][icut] = MakeNewHisto(histName.Data(), Form("%s;Vertex-based #mu Direction Resolution (#circ)",histName.Data()),  150*10, 0, 30, h_vect_hists[imc][ireco]);
	
	histName = Form("%s_cosmic_range_res_reco%d_cut%d",h_name[imc],ireco,icut);
	h_cosmic_range_res[imc][ireco][icut]     = MakeNewHisto(histName.Data(), Form("%s; Range_{Rec.} - Range_{True} (cm)",histName.Data()), 300, -200, 200, h_vect_hists[imc][ireco]);
	  
  
      }
    }
  }
}

void OfficialEventParser::Parse(int ireco){
  // This is called event by event, and filling histograms.

  // initialize flag event by event if need
  InitEvent();
	
  if (get_datatype() == fcData || get_datatype() == fcMC){

    if ( get_datatype() == fcData && ! __FC_RUN_SELECTION__ ) {
    } else {


      if (isCosmic(ireco)) FillCosmic(ireco);
      

    }
  }
}


bool OfficialEventParser::isCosmic(int ireco){

  bool isTrue = true;
  
  return isTrue;

}


void OfficialEventParser::FillCosmic(int ireco){

  int imc=-1;
  if (get_datatype() == fcMC) imc = 1; // mc
  else if (get_datatype() == fcData) imc = 0;

  // Count events for normalization
  nentries_norm[imc][ireco]++;

  // Get Truth Info
  for (int ise=0; ise<nSE; ise++) 
    found_true_particle[ise] = get_truth_info(imc, ireco, ise, pid_tru[ise], mom_tru[ise], v_dir_tru[ise], v_vtx_tru[ise], etimev);

  // Get Reconstructed info

  // Reset
  for (int imrhyp=0; imrhyp<nMRhyps; imrhyp++) {
    MRnll[imrhyp]=0;
    //maxMom[imrhyp]=0; minMom[imrhyp]=0;
	  
    for (int imrhyp2=0; imrhyp2<nMRhyps; imrhyp2++) 
      MrPID[imrhyp][imrhyp2]=0;
  }
	
  MRnll[0] = getMin1Rnll(ireco, MrPID[0][0]);
  MRnll[1] = getMin2Rnll(ireco, MrPID[1][0], MrPID[1][1]);
  MRnll[2] = getMin3Rnll(ireco, MrPID[2][0], MrPID[2][1], MrPID[2][2]);
  MRnll[3] = getMin4Rnll(ireco, MrPID[3][0], MrPID[3][1], MrPID[3][2], MrPID[3][3]);

  //lnLrat1R2R_muon = getLnLrat1R2R(ireco,fillPID); 

  etime = get_dcye_time(ireco, ise_dcye);
      
  for (int ise=0; ise<nSE; ise++) {
    
    // Momentum
    mom_rec[ise] = get_1Rmom(ireco, ise);

    imom_rec_idx[ise] = h_cosmic_mom_binning[imc][ireco][ise]->FindBin(mom_rec[ise])-1;
    
    //likelihood_mom[ise] = mom_rec[ise];
    //if ( ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) && ise==ise_dcye) likelihood_mom[ise_dcye]=E->fq1rmom[ise_dcye][ie];
    
    // Vertex
    TVector3 vertex(get_vertex(ireco, ise));
    v_vtx_rec[ise] = vertex;
    
    float vertex_xyz[3];
    v_vtx_rec[ise].GetXYZ(vertex_xyz);
    wall[ise] = get_wall(vertex_xyz);  
    
    // Direction
    TVector3 dir(get_1Rdir(ireco, ise));
    v_dir_rec[ise] = dir;

    // Likelihoods
    lnLrat1R[ise] = getLnLrat1R(ireco, ise);
    
    r2[ise] = v_vtx_rec[ise][0]*v_vtx_rec[ise][0]+v_vtx_rec[ise][1]*v_vtx_rec[ise][1]; 
    z[ise] = v_vtx_rec[ise][2];

    // Prefit vertex
    TVector3 Tprftvtx(E->fqtwnd_prftvtx[ise]);
    TVector3 Tvtx_diff = vertex - Tprftvtx;
    abs_vtx_prfvtx_diff[ise] = Tvtx_diff.Mag();  
  }
  
  // Resolutions
  if (imc) {
    for (int ise=0; ise<nSE; ise++) {
      if (found_true_particle[ise]>0) {
	mom_res[ise] = (mom_rec[ise]-mom_tru[ise])/mom_tru[ise]*100;
	dir_res[ise] = v_dir_rec[ise].Angle(v_dir_tru[ise])*180/TMath::Pi();
	vtx_res[ise] = (v_vtx_rec[ise]-v_vtx_tru[ise]).Mag();
	lvtx_res[ise] = (v_vtx_rec[ise]-v_vtx_tru[ise]).Dot(v_dir_tru[ise]);
	tvtx_res[ise] = ((v_vtx_rec[ise]-v_vtx_tru[ise]).Cross(v_dir_tru[ise])).Mag();
	if (ise==1) time_res = (etime-etimev)*1000;

	int icut=0;
	h_cosmic_mom_res_vs_trumom[imc][ise][ireco][icut]->Fill(mom_tru[ise],mom_res[ise]);
	//FillResolutionHistos(imc, ireco, ise, 0);

	//h_cosmic_momv[imc][ise][ireco][icut]->Fill(mom_tru[ise], weight());
      }
    }
  }
  ddir_res = -1; // Needs 2 sub-events below
  range_res = -999999; // Needs 2 sub-events below

  double year = get_year(ireco);
  double run = get_run(ireco);
      
  //cout << "year run " << year << " " << run << endl;

  // Direction of 2 rings
  //float dir2R[2][3] = {{0}};
  //get_2Rdir(ireco, dir2R[0], dir2R[1]);
  
  // Get total ID charge
  double potot = get_potot(ireco);
  //h_cosmic_potot[imc][ireco]->Fill(potot);


  int    itwnd[nSE];
  int    icluster[nSE];

  double cltstart[nSE];
  double cltend[nSE];
  double cllength[nSE];

  int    clnhits[nSE];
  float  cltotq[nSE];
  double clnhitsnorm[nSE];
  double cltotqnorm[nSE];


  int    fqnhits[nSE];
  float  fqtotq[nSE];
  double fqnhitsnorm[nSE];
  double fqtotqnorm[nSE];

  float  fq1rtotmu[nSE];

  for (int ise=0; ise<nSE; ise++) {
    itwnd[ise] = E->fqitwnd[ise];		      
    icluster[ise] = E->fqtwnd_iclstr[itwnd[ise]];	      
    
    cltstart[ise] = E->cluster_tstart[icluster[ise]]/1000;    
    cltend[ise] = E->cluster_tend[icluster[ise]]/1000;	      
    cllength[ise] = cltend[ise]-cltstart[ise];		      
                                                     
    clnhits[ise] = E->cluster_nhits[icluster[ise]];	      
    cltotq[ise] = E->cluster_totq[icluster[ise]];	      
    clnhitsnorm[ise] = clnhits[ise]/cllength[ise]/1000;	      
    cltotqnorm[ise] = cltotq[ise]/cllength[ise]/1000;	      
    
    fqnhits[ise] = E->fqnhitpmt[ise];		      
    fqtotq[ise] = E->fqtotq[ise];		      
    fqnhitsnorm[ise] = fqnhits[ise]/cllength[ise]/1000;	      
    fqtotqnorm[ise] = fqtotq[ise]/cllength[ise]/1000;	      
    
    fq1rtotmu[ise] = E->fq1rtotmu[ise][get1rType(ireco,ise)];     
  }


  // Cut: potot for sufficient light in event
  //if (potot > 1000) 
  h_cosmic_fqtotqnorm_vs_year[imc][ise_muon][ireco]->Fill(year, fqtotqnorm[ise_muon]);
  h_cosmic_cllength_vs_year[imc][ise_muon][ireco]->Fill(year, cllength[ise_muon]);
  if ( ( ( (ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) && E->stpotot<110e3 ) ||
       //( ( (ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun ) && E->potot<110e3 ) //fqtotqnorm[ise_muon]<75 )
       ( ( (ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun ) && fqtotqnorm[ise_muon]<75 )
       ) {
  

      h_cosmic_mom_binning[imc][ireco][ise_muon]->Fill(mom_rec[ise_muon]);

      
    // Cut: PID must be muon-like cut (fiTQun), or Muon-goodness (stmufit) 
    //if ( (ireco==apfit && get1rType(ireco,ise_muon) == 3 ) ||
    if ( ( ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) && lnLrat1R[ise_muon] >= -0.9 ) ||
	 ( (ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun ) ) { // && get1rType(ireco,ise_muon) == imu ) ) {

      // Single ring likelihood
      
      // Cut: 1-ring
      //if (is1R(ireco)) 
      { 
	
	// Sub-GeV or Multi-GeV (currently based on momentum)
	
	bool bIsSubGeV = mom_rec[ise_muon]<1330; //isSubGeV(ireco);
	if ( (_nrgSepType == 0 && mom_rec[ise_muon]<_maxE) || 
	     (_nrgSepType == 1 && bIsSubGeV) ||
	     (_nrgSepType == 2 && !bIsSubGeV && mom_rec[ise_muon]<_maxE) ) {
	  
	  int nmue = muedcy(ireco);
	  h_cosmic_nmue_vs_year[imc][ireco]->Fill(year, nmue);

	  // Cut: 1 decay-e
	  //if (nmue==1 && cllength>0.7) {
	  if (nmue==1) {

	    TVector3 vertex_diff = v_vtx_rec[ise_dcye]-v_vtx_rec[ise_muon];
	    TVector3 vertex_diff_mc = v_vtx_tru[ise_dcye]-v_vtx_tru[ise_muon];
	    double range = vertex_diff.Mag();
	    double range_mc = vertex_diff_mc.Mag();
	    float MomOverRange = mom_rec[ise_muon]/range;

	    ddir_res = vertex_diff.Angle(v_dir_rec[ise_muon])*180/TMath::Pi();
	    range_res = range - range_mc;

	    int n50 = get_n50(ireco,ise_dcye);
	    float q50 = get_q50(ireco,ise_dcye);
	
	    int inGate = isInGate(ireco,ise_dcye);
	    h_cosmic_ingate_vs_year[imc][ireco]->Fill(year, inGate);	    
	    if (!inGate) {
	      
	      h_cosmic_lifetime_vs_year[imc][ireco]->Fill(year, etime);
	      if (1.2<etime) {
		
		// Cut: Fit goodness (stmufit) and pcflg (fiTQun)
		for (int ise=0; ise<2; ise++) h_cosmic_pcflg_vs_year[imc][ise][ireco]->Fill(year, E->fq1rpcflg[ise][get1rType(ireco,ise)]);
		if ( ( ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) && lnLrat1R[ise_dcye]>0.5) ||
		     ( ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) && E->fq1rpcflg[ise_dcye][ie]==0 && E->fq1rpcflg[ise_muon][imu]==0 
		       //&& getLnLrat1R(ireco,ise_dcye)<0 
		       ) ) { // && get1rType(ireco,ise_dcye) == ie ) ) {    
		  
		  for (int ise=0; ise<2; ise++) if (etime<10) { h_cosmic_vtx_prfvtx_diff_vs_year[imc][ise][ireco]->Fill(year, abs_vtx_prfvtx_diff[ise]); }
	      	  
		  // Vertex at top or side
		  double r = sqrt(r2[ise_muon]);
		  Bool_t topEnter = (z[ise_muon] > 1750.0 && z[ise_muon]<1850.0) && r < 1600.0;	      
		  Bool_t sideEnter = (z[ise_muon] < 1750.0) && (r> 1600.0 && r < 1750.0);

		  for (int ise=0; ise<2; ise++) {
		    if (etime<10) {
		      
		      if ((ise==ise_dcye && (topEnter || sideEnter)) || ise==ise_muon) {
			h_cosmic_wall_vs_year[imc][ise][ireco]->Fill(year, wall[ise]);
			h_cosmic_z_vs_year[imc][ise][ireco][0]->Fill(year, v_vtx_rec[ise][2]); 
			h_cosmic_r2_vs_year[imc][ise][ireco][0]->Fill(year, r2[ise]);
			h_cosmic_zzoom_vs_year[imc][ise][ireco][0]->Fill(year, v_vtx_rec[ise][2]); 
			h_cosmic_r2zoom_vs_year[imc][ise][ireco][0]->Fill(year, r2[ise]);
			
			if (r<1400) h_cosmic_z_nearwall[imc][ise][ireco][0]->Fill(v_vtx_rec[ise][2]); 
			else if (fabs(z[ise_muon])<1500) h_cosmic_r2_nearwall[imc][ise][ireco][0]->Fill(r2[ise]);
		      }
		    }
		  }


		  // Dcye-e vertex in ID
		  if (wall[ise_dcye]>100) {
		    
		    if (etime<10) {

		      if (topEnter || sideEnter) {
			
			for (int ilike=0; ilike<nMuonLikes; ilike++) {
		    
			  double likelihood;
			  int fillPID;
			  if (ilike==0) likelihood = -E->fq1rnll[ise_muon][ie]+E->fq1rnll[ise_muon][imu]; //getLnLrat1R(ireco,ise_muon);
			  else if (ilike==1) likelihood = -E->fq1rnll[ise_muon][ie]+E->fqpi0nll[fq_pi0_mode];
			  else if (ilike==2) likelihood = -E->fq1rnll[ise_muon][imu]+E->fqpi0nll[fq_pi0_mode];
			  else if (ilike==3) likelihood = getLnLrat1R2R(ireco, fillPID);
			  else if (ilike==4) likelihood = getRingCountingLikelihood(ireco);
		    
			  int mom_idx = imom_rec_idx[ise_muon];
			  //mom_idx = min(mom_idx, nMuonMomBinsLike-1);
			  //mom_idx = max(0, mom_idx);
		    
			  if (mom_idx<0 || mom_idx>=nMuonMomBinsLike) continue;
			  if (likelihood!=likelihood) continue;
		    
			  h_cosmic_Lmuon_vs_year[imc][ireco][ilike][mom_idx]->Fill(year, likelihood);
			  h_cosmic_Lmuon_vs_zenith[imc][ireco][ilike][mom_idx]->Fill(-v_dir_rec[ise_muon][2], likelihood);		  
			}
			
			float mom_like = get_1Rmom(ireco, ise_muon);
			if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun)
			  mom_like = E->fq1rmom[ise_muon][ie];
			h_cosmic_lnLrat1R_vs_mom[imc][ise_muon][ireco]->Fill(mom_like,-E->fq1rnll[ise_muon][ie]+E->fq1rnll[ise_muon][imu]);
			h_cosmic_lnLrat2R1R_vs_mom[imc][ise_muon][ireco]->Fill(mom_like,getLnLrat1R2R(ireco, fillPID));
			h_cosmic_lnLrat2R1R_vs_mom[imc][1][ireco]->Fill(mom_like,getLnLrat1R2R(ireco, fillPID));
		      }
		    }
				  
		    h_cosmic_range_vs_mumom_vs_year[imc][ireco]->Fill(year, mom_rec[ise_muon]/1000, range/100);

		    if (range<rangeHigh*100) {
		      int range_idx = (int)((range/100-rangeLow)/rangeStep);
		    
		      // Official E-Scale
		      if (etime>2 && v_dir_rec[ise_muon][2]<-0.94 && hasRvtxInFid(ireco)) 
			if ( (ireco==0 && fabs(z[ise_muon]-1810)<0.1) || (ireco==1 && topEnter) )
			  h_cosmic_mom_over_range_rangesep_vs_year[imc][ireco][range_idx][0]->Fill(year, MomOverRange);
		      if (etime<10) {
			if (topEnter || sideEnter)
			  h_cosmic_mom_over_range_rangesep_vs_year[imc][ireco][range_idx][1]->Fill(year, MomOverRange);
			if (topEnter)
			  h_cosmic_mom_over_range_rangesep_vs_year[imc][ireco][range_idx][2]->Fill(year, MomOverRange);
			else if (sideEnter)
			  h_cosmic_mom_over_range_rangesep_vs_year[imc][ireco][range_idx][3]->Fill(year, MomOverRange);
		      }
		      
		    }
		  }

		  // Cut: Decay-e vertex in FV
                  if (wall[ise_dcye]>100 && etime<10 && (topEnter || sideEnter)) {
		
		    for (int ise=0; ise<nSE; ise++) {

		      int icut=0;
		      //if (wall[ise_dcye]<200)
		      //if (8.5 < mom_tru[ise_dcye] && mom_tru[ise_dcye] < 11.5)
		      FillResolutionHistos(imc, ireco, ise, icut);
		    
		      h_cosmic_cltstart_vs_year[imc][ise][ireco]->Fill(year, cltstart[ise]);
		      h_cosmic_cltend_vs_year[imc][ise][ireco]->Fill(year, cltend[ise]);
		      if (ise==ise_dcye) h_cosmic_cllength_vs_year[imc][ise][ireco]->Fill(year, cllength[ise]);
		      h_cosmic_cllength_vs_cltstart[imc][ise][ireco]->Fill(cltstart[ise], cllength[ise]);
		    		    
		      h_cosmic_clnhits_vs_year[imc][ise][ireco]->Fill(year, clnhits[ise]);
		      h_cosmic_cltotq_vs_year[imc][ise][ireco]->Fill(year, cltotq[ise]);
		      h_cosmic_cltotqnorm_vs_year[imc][ise][ireco]->Fill(year, cltotqnorm[ise]);
		      h_cosmic_clnhitsnorm_vs_year[imc][ise][ireco]->Fill(year, clnhitsnorm[ise]);
		    
		      h_cosmic_fqnhits_vs_year[imc][ise][ireco]->Fill(year, fqnhits[ise]);
		      h_cosmic_fqtotq_vs_year[imc][ise][ireco]->Fill(year, fqtotq[ise]);
		      h_cosmic_fqnhitsnorm_vs_year[imc][ise][ireco]->Fill(year, fqnhitsnorm[ise]);
		      if (ise==ise_dcye) h_cosmic_fqtotqnorm_vs_year[imc][ise][ireco]->Fill(year, fqtotqnorm[ise]);
		    }

		    h_cosmic_n50_vs_year[imc][ireco]->Fill(year, n50);
		    //if (n50>60) 
		    {			
		    
		      //h_cosmic_lifetime_vs_year[imc][ireco]->Fill(year, etime);
		      h_cosmic_emom_vs_lifetime_vs_year[imc][ireco]->Fill(year, etime, mom_rec[ise_dcye]);

		      // Final Decay-e spectrum
		      //if (1.2<etime && etime<10) 
		      {

			for (int ise=0; ise<2; ise++)
			  h_cosmic_fq1rtotmu_vs_year[imc][ise][ireco]->Fill(year, fq1rtotmu[ise]);

			h_cosmic_q50_vs_year[imc][ireco]->Fill(year, q50);
			
			h_cosmic_emom_vs_year[imc][ireco]->Fill(year, mom_rec[ise_dcye]);
			h_cosmic_mumom_vs_year[imc][ireco]->Fill(year, mom_rec[ise_muon]);

			h_cosmic_mom_binning[imc][ireco][ise_dcye]->Fill(mom_rec[ise_dcye]);
			
			//for (int ise=0; ise<2; ise++) h_cosmic_vtx_prfvtx_diff_vs_year[imc][ise][ireco]->Fill(year, abs_vtx_prfvtx_diff[ise]);

			for (int ilike=0; ilike<nDcyeLikes; ilike++) {
			
			  //if (20<mom_rec[ise_dcye] && mom_rec[ise_dcye]<60) 
			  //if (r2[ise_dcye]<10000) 
			  {

			    double likelihood;
			    if (ilike==0) likelihood = -E->fq1rnll[ise_dcye][ie]+E->fq1rnll[ise_dcye][imu]; //getLnLrat1R(ireco,ise_dcye);
			    //else if (ilike==1) likelihood = getLnLrat1R2R(ireco, fillPID);
			    //else if (ilike==1) likelihood = -E->fq1rnll[ise_dcye][ie]+E->fqpi0nll[fq_pi0_mode];
			    //else if (ilike==2) likelihood = -E->fq1rnll[ise_dcye][imu]+E->fqpi0nll[fq_pi0_mode];
			
			    int mom_idx = imom_rec_idx[ise_dcye];
			    //mom_idx = min(mom_idx, nDcyeMomBinsLike-1);
			    //mom_idx = max(0, mom_idx);
		  
			    if (mom_idx<0 || mom_idx>=nDcyeMomBinsLike) continue;
			    if (likelihood!=likelihood) continue;
		  
			    h_cosmic_Ldcye_vs_year[imc][ireco][ilike][mom_idx]->Fill(year, likelihood);
			    h_cosmic_Ldcye_vs_zenith[imc][ireco][ilike][mom_idx]->Fill(-v_dir_rec[ise_dcye][2], likelihood);		  
			  }
			}

			float mom_like = get_1Rmom(ireco, ise_dcye);
			if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun)
			  mom_like = E->fq1rmom[ise_dcye][ie];
			h_cosmic_lnLrat1R_vs_mom[imc][ise_dcye][ireco]->Fill(mom_like,-E->fq1rnll[ise_dcye][ie]+E->fq1rnll[ise_dcye][imu]);      
		      }
		    }
		  }
		}
		}
	      }
	  }
	}
      }
    }  
  }
}


  
void OfficialEventParser::FillResolutionHistos(int imc, int ireco, int ise, int icut) {

  if (ise==0) h_cosmic_ddir_res[imc][ireco][icut]->Fill(ddir_res);
  
  if (!imc) return;
  	
  if (found_true_particle[ise]<1) return;

  h_cosmic_mom_res[imc][ise][ireco][icut]->Fill(mom_res[ise]);
  //h_cosmic_mom_res_vs_trumom[imc][ise][ireco][icut]->Fill(mom_tru[ise],mom_res[ise]);
  h_cosmic_dir_res[imc][ise][ireco][icut]->Fill(dir_res[ise]);
  h_cosmic_vtx_res[imc][ise][ireco][icut]->Fill(vtx_res[ise]);			    
  h_cosmic_lvtx_res[imc][ise][ireco][icut]->Fill(lvtx_res[ise]);
  h_cosmic_tvtx_res[imc][ise][ireco][icut]->Fill(tvtx_res[ise]);
  
  if (ise==1) h_cosmic_time_res[imc][ise][ireco][icut]->Fill(time_res);
  else if (ise==0) h_cosmic_range_res[imc][ireco][icut]->Fill(range_res);

}

//void OfficialEventParser::FillLikelihoodHistos(int imc, int ireco, int ise, int icut) {
//  
//  h_cosmic_lnLrat1R_vs_mom_year[imc][ise][ireco]->Fill(likelihood_mom[ise],lnLrat1R[ise]);
//
//  if (ise==ise_muon) {
//    h_cosmic_lnLrat1R2R_vs_mom_year[imc][ireco]->Fill(mom_rec[ise],lnLrat1R2R_muon); 
//    h_cosmic_lnLrat1R2R_vs_lnLrat1R[imc][ireco]->Fill(lnLrat1R[ise],lnLrat1R2R_muon); 
//  }
//}

double OfficialEventParser::get_run(int ireco){

  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    return E->stnrun;
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    return E->nrun;
  }

}

double OfficialEventParser::get_year(int ireco){

  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    return E->stnday[0]+1900 + (E->stnday[1]-0.5)/12.;
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    return E->date[0]+1900 + (E->date[1]-0.5)/12.;
  }

}



int OfficialEventParser::get_truth_info(int imc, int ireco, int ise, int &pid_tru, float &mom_tru, TVector3 &v_dir_tru, TVector3 &v_vtx_tru, float &etimev){

  if (!imc) return 0;
  
  int ise_muon = 0;
  int ise_dcye = 1;

  int found_true_particle = 0;
  
  // Muon information
  if (ise==ise_muon) {

    if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
      
      if (E->ipv[0]==5) pid_tru = -13;
      else if (E->ipv[0]==6) pid_tru = 13;
      else {
	cout << "Error: Unknown PID: " << E->ipv[0] << endl;
	cout << "  Setting as most likely mu-" << endl;
	pid_tru = 13;
	//exit (-1);
      }

      mom_tru = E->pmomv[0]; //Abspvc[0];
      for (int ix=0; ix<3; ix++) {
	v_dir_tru[ix] = E->dirv[0][ix]; //Pvc[0][ix]/E->Abspvc[0];
	v_vtx_tru[ix] = E->posv[ix];
      }
    }
    
    else if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {

      if (E->stdcyepidv==-11) pid_tru = -13;
      else if (E->stdcyepidv==11) pid_tru = 13;
      else if (E->stdcyepidv==22) pid_tru = 13;
      else {
	cout << "Error: Unknown PID: " << E->ipv[0] << endl;
	cout << "  Setting as most likely mu-" << endl;
	pid_tru = 13;
	//exit (-1);
      }

      mom_tru = E->stamomv;
      for (int ix=0; ix<3; ix++) {
	v_dir_tru[ix] = E->stdirvmu[ix];
	v_vtx_tru[ix] = E->stposvmu[ix];
      }
      etimev = E->stdecaytv/1000;
    }
    
    found_true_particle = 1;
    
  }
  
  // Electron information
  else if (ise==ise_dcye) {
  
    if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
      // Find true dcy-e info in particle stack
      int ip=0;
      for (ip=0; ip<E->nscndprt; ip++) {
			      
	for (int ix=0; ix<3; ix++) {
	  v_dir_tru[ix] = E->pscnd[ip][ix];
	  v_vtx_tru[ix] = E->vtxscnd[ip][ix];
	}
	mom_tru = v_dir_tru.Mag();
	v_dir_tru *= 1/mom_tru;
	
	etimev = E->tscnd[ip]/1000;
	
	pid_tru = E->iprtscnd[ip];

	//if ((abs(E->iprtscnd[ip])==11 || abs(E->iprtscnd[ip])==22) && 
	if (abs(E->iprtscnd[ip])==11 && 
	    E->lmecscnd[ip]==5) {
	  found_true_particle++;
	  break;
	}
      }
      if (found_true_particle>1) {
	cout << "Error expected 1 true decay-e but found " << found_true_particle << " " << E->tscnd[ip-1] << " " << E->iprtscnd[ip-1] << endl;
	//exit (1);
      }
    }
    
    
    else if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {

      pid_tru = E->stdcyepidv;

      mom_tru = E->stamomve;
      for (int ix=0; ix<3; ix++) {
	v_dir_tru[ix] = E->stdirve[ix];
	v_vtx_tru[ix] = E->stposve[ix];
      }
      
      found_true_particle = 1;
    }

  }
  
  return found_true_particle;
  
}



float OfficialEventParser::get_nhits_over_time(int ireco){

  float hits_time = 0;
  
  for (int icluster=0; icluster<E->cluster_ncand; icluster++) {
    
    if (E->cluster_goodflag[icluster]==1) {
      hits_time = E->cluster_nhits[icluster]/(E->cluster_tend[icluster]-E->cluster_tstart[icluster]);
      break;
    }
  }
  
  return hits_time;
  
}

int OfficialEventParser::get_n50(int ireco, int ise){
  
  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    if (ise<=0) return 0;
    else return E->stn50;
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    return E->fqn50[ise];
  }
  
}


int OfficialEventParser::get_q50(int ireco, int ise){
  
  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    if (ise<=0) return 0;
    else return E->stq50;
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    return E->fqq50[ise];
  }
  
}

float OfficialEventParser::get_dcye_time(int ireco, int ise){
  
  float time = -1;
  
  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    if (ise==0) time = 0;
    //else if (ise>0) time = E->etime[ise-1];
    else time = E->stt;
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    if (ise>=0) {
      int i1RPID0 = get1rType(ireco, 0);
      int i1RPID1 = get1rType(ireco, ise);

      time  = (E->fq1rt0[ise][i1RPID1] - E->fq1rt0[0][i1RPID0])/1000;
  
      
    }
  }
  
  return time;
}


int OfficialEventParser::get1rType(int ireco, int ise){
  
  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    return E->ip[0];
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {

    //double lnLrat1R=getLnLrat1R(ireco,ise);
    //if (lnLrat1R > E->fq1rmom[ise][ie]*0.2) //PID'd as electron
    //  return ie;
    //else //muon
    //  return imu;
    
    //if (getLnLrat1R(ireco, ise) < 0) return imu;
    //else return ie;

    if (ise==0) return imu;
    else return ie;

  }
}

double OfficialEventParser::getLnLrat1R(int ireco, int ise){
  /*
  double val = 0;

  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    if (ise==0) val = E->stbgood;
    else val = E->stegood;
    // val = E->egood[ise-1];
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    
    val = -E->fq1rnll[ise][ie]+E->fq1rnll[ise][imu];
  }

  //if (val!=val) cout << "getLnLrat1R() Warning: val = " << val << endl;


  return val;
  /**/

  double likelihood=-999999;

  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    if (ise==0) likelihood = E->stbgood;
    else likelihood = E->stegood;
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    
    // For sub-GeV cut
    double slope = 0.2;
    
    double lnLrat1R=TMath::Sign(1.,slope)*(-E->fq1rnll[ise][ie]+E->fq1rnll[ise][imu]);
    
    double x = 999999;
    TVector2 v(-x,-x*slope);
    TVector2 w(x,x*slope);
    TVector2 p(E->fq1rmom[ise][ie],lnLrat1R);
    double distance = minimum_distance(v,w,p);
    
    if (lnLrat1R < E->fq1rmom[ise][ie]*slope) likelihood = -distance;
    else likelihood = distance;      
    
  }
  
  //if (likelihood<PID_Likelihood_multi_range[ireco][0]) likelihood=PID_Likelihood_multi_range[ireco][0]+.1;
  //else if (likelihood>PID_Likelihood_multi_range[ireco][1]) likelihood=PID_Likelihood_multi_range[ireco][1]-.1;

  return likelihood;

  
}



double OfficialEventParser::getLnLrat1R2R(int ireco, int &fillPID){

  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) return 0;
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    //int ise = 0;
    //int i1RPID = -1; // = get1rType(ireco,ise);
    //if (getLnLrat1R(ireco, ise) < 0) i1RPID = imu;
    //else i1RPID = ie;
    //
    //return -getMin2Rnll(ireco,fillPID) + E->fq1rnll[ise][i1RPID];

    return -MRnll[1] + MRnll[0];
  }

}


double OfficialEventParser::getMin1Rnll(int ireco, int &PID1){
  
  double minNll = 9999999;

  // For minimum (e.g. for generic AtmNu analysis)
  /*
  for (int iPID1=1; iPID1<=3; iPID1++) {
    if (iPID1==imu) continue; 
    if (E->fq1rnll[0][iPID1]<minNll) { 
      minNll = E->fq1rnll[0][iPID1];
      PID1 = iPID1;
    }
  }
  /**/

  // For T2K 1R analysis
  int ise=0;
  //PID1 = get1rType(ireco, ise);
  if (getLnLrat1R(ireco, ise) < 0) PID1 = imu;
  else PID1 = ie;

  minNll = E->fq1rnll[ise][PID1];

  //for (int ix=0; ix<3; ix++) 
  //  MrDir[0][0][ix] = E->fq1rdir[0][PID1][ix];

  
  return minNll;

}

double OfficialEventParser::getMin2Rnll(int ireco, int &PID1, int &PID2){
  
  double minNll = 9e9;
  for (int iPID1=0; iPID1<2; iPID1++) 
    for (int iPID2=0; iPID2<2; iPID2++) 
      if (E->fq2rnll[iPID1][iPID2]<minNll) {
	minNll = E->fq2rnll[iPID1][iPID2];
	PID1 = iPID1;
	PID2 = iPID2;

	//for (int ix=0; ix<3; ix++) {
	//  MrDir[1][0][ix] = E->fq2rdir1[iPID1][iPID2][ix];
	//  MrDir[1][1][ix] = E->fq2rdir2[iPID1][iPID2][ix];
	//}
      }

  return minNll;

}

double OfficialEventParser::getMin3Rnll(int ireco, int &PID1, int &PID2, int &PID3){
  
  double minNll = 9999999;
  for (int iPID1=0; iPID1<2; iPID1++) 
    for (int iPID2=0; iPID2<2; iPID2++) 
      for (int iPID3=0; iPID3<2; iPID3++) 
	if (E->fq3rnll[iPID1][iPID2][iPID3]<minNll) {
	  minNll = E->fq3rnll[iPID1][iPID2][iPID3];
	  PID1 = iPID1;
	  PID2 = iPID2;
	  PID3 = iPID3;

	  //for (int ix=0; ix<3; ix++) {
	  //  MrDir[2][0][ix] = E->fq2rdir1[iPID1][iPID2][ix];
	  //  MrDir[2][1][ix] = E->fq2rdir2[iPID1][iPID2][ix];
	  //  MrDir[2][2][ix] = E->fq3rdir3[iPID1][iPID2][iPID3][ix];
	  //}

    	}

  return minNll;

}


double OfficialEventParser::getMin4Rnll(int ireco, int &PID1, int &PID2, int &PID3, int &PID4){
  
  double minNll = 9999999;
  for (int iPID1=0; iPID1<2; iPID1++) 
    for (int iPID2=0; iPID2<2; iPID2++) 
      for (int iPID3=0; iPID3<2; iPID3++) 
	for (int iPID4=0; iPID4<2; iPID4++) 
	  if (E->fq4rnll[iPID1][iPID2][iPID3][iPID4]<minNll) {
	    minNll = E->fq4rnll[iPID1][iPID2][iPID3][iPID4];
	    PID1 = iPID1;
	    PID2 = iPID2;
	    PID3 = iPID3;
	    PID4 = iPID4;
	    
	    //for (int ix=0; ix<3; ix++) {
	    //  MrDir[3][0][ix] = E->fq2rdir1[iPID1][iPID2][ix];
	    //  MrDir[3][1][ix] = E->fq2rdir2[iPID1][iPID2][ix];
	    //  MrDir[3][2][ix] = E->fq3rdir3[iPID1][iPID2][iPID3][ix];
	    //  MrDir[3][3][ix] = E->fq4rdir4[iPID1][iPID2][iPID3][iPID4][ix];
	    //}
	  }

  return minNll;

}


double OfficialEventParser::getRingCountingLikelihood(int ireco){

  double likelihood=-999999;

  if (ireco==apfit) return 0; //likelihood = E->Dlfct;
  else if (ireco==fitqun) {

    likelihood = MRnll[0]-MRnll[1];

    // Muon
    if (MrPID[0][0]==imu) {
      
      int PID1;
      int ise = 0;
      if (getLnLrat1R(ireco, ise) < 0) PID1 = imu;
      else PID1 = ie;
      
      // The following was for when ipip was used for 1R PID and caused
      // the swoosh in the 2R/1R likelihood (actually imu causes swoosh in v3 now too)
      double f = E->fq1rmom[ise][ie]; //mom_rec[0]; //maxMom[0];
      double g = likelihood;

      // y = pow( (x+a)/b , 2 ) + c      
      
      // Prior to v3r1
      //double a = 300;
      //double b = 75;
      //double c = 100;

      // v3r1
      double a = 0;
      double b = 60;
      double c = 150;

      double b2 = pow(b,2);
      double b4 = pow(b,4);

      double x = pow( 108*a*b4 + sqrt( pow(108*b4*f + 108*a*b4 , 2) + 4*pow( 6*b4 + 12*b2*c - 12*b2*g , 3) ) + 108*b4*f , 1./3.) / (6*pow(2,1./3.)) - (6*b4 + 12*b2*c - 12*b2*g) / ( 3*pow(2,2./3.)*pow(108*a*b4 + sqrt( pow(108*b4*f + 108*a*b4 , 2) + 4*pow(6*b4 + 12*b2*c - 12*b2*g , 3) ) + 108*b4*f, 1./3.) ) - a;
      
      double d = sqrt( pow(f-x,2) + pow( pow((x+a)/b ,2) + c - g , 2 ) );
      if (g < pow((f+a)/b,2)+c) d *= -1;

      likelihood = d;
      /**/

      // For v2r1 when using imu
      //likelihood -= 100; // For generic AtmNu analysis
      //likelihood -= 150;  // For T2K 1R analyis
    }
  
    // Electron
    else if (MrPID[0][0]==ie) {

      //likelihood -= 100; // For generic AtmNu analysis
      likelihood -= 150; // For T2K 1R analyis
      
    }
  }

  return likelihood;
}

bool OfficialEventParser::is1R(int ireco){
  // 1 ring event
  bool isTrue = false;
  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    //if (E->nring == 1) isTrue = true;
    isTrue = true;
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
	  
    int fillPID;
    double lnLrat1R2R=getLnLrat1R2R(ireco,fillPID); //likelihood ratio: best 2R vs. best 1R
    if (lnLrat1R2R<150.) isTrue = true;
	  	    
  }
	
  return isTrue;
}

// Currently based on momentum of best PID
int OfficialEventParser::isInGate(int ireco, int ise){

  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    //if (!E->etype[ise-1]) return true;
    if (E->stmuetype!=1) return true;
    else return false;
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    return E->fqipeak[ise];
  }

}


bool OfficialEventParser::isSubGeV(int ireco){

  bool isTrue = false;
  float mom = get_1Rmom(ireco, 0);
  if (mom < 1330) isTrue = true;
  return isTrue;

}

// Warning: slightly different from E->wall
//float OfficialEventParser::get_wall(float *vertex) {
//      
//  double r = vertex[0]*vertex[0] + vertex[1]*vertex[1]; 
//  r = sqrt( r ) ;
//  float wall = RINTK - r ; 
//  double z1 = ZPINTK - vertex[2] ;
//  double z2 = vertex[2] - ZMINTK ;
//  if ( z1<wall ) wall = z1 ; 
//  if ( z2<wall ) wall = z2 ; 
//  return wall ; 
//  
//}

float OfficialEventParser::get_wall(float *vertex) {
  
  double x = vertex[0];
  double y = vertex[1];
  double z = vertex[2];

  double r = sqrt(x*x+y*y);
  double wallc;
  double capwal;
  double cylwal;
  
  // First check if we're inside ID:
  if(fabs(z)>ZPINTK || r>RINTK){
    // we are outside ID.  Are we above top endcap?
    if (z>ZPINTK) {
      // Above top. If inside (outside) cylinder projection, closest point is endcap (corner):
      if (r<RINTK) {
	wallc=ZPINTK-z;
      }      
      else {
	wallc=-sqrt((ZPINTK-z)*(ZPINTK-z)+(r-RINTK)*(r-RINTK));
      }
    }
    else if (z>-ZPINTK) {
      // Vertically between caps.  Closest point is cylinder:
      wallc = RINTK-r; 
    } 
    else {
      // Below bottom.
      if (r<RINTK) {
	wallc=ZPINTK+z;
      } 
      else {
	wallc=-sqrt((ZPINTK+z)*(ZPINTK+z)+(r-RINTK)*(r-RINTK));
      }
    } 
  }
  else { // We are inside ID
    cylwal = RINTK-r;
    capwal = ZPINTK-fabs(z);
    if (cylwal<capwal) {wallc=cylwal;}
    else {wallc=capwal;}
  }

  return wallc;

}


bool OfficialEventParser::hasRvtxInFid(int ireco){
  // vertex is in the Fid in r-direction :  r=sqrt(x*x + y*y)  is less than 1490.
  bool isTrue = false;

  float *vertex = get_vertex(ireco);

  isTrue = RvtxInFid(vertex);

  return isTrue;
  // R (sqrt(pos[1]**2+pos(2)**2)<1490.)
}

bool OfficialEventParser::RvtxInFid(float *vertex){

  bool isTrue = false;

  double r = sqrt(vertex[0] * vertex[0] + vertex[1] * vertex[1]);
  if (r < RINTK-200) isTrue = true;

  return isTrue;
}

bool OfficialEventParser::hasZvtxInFid(int ireco){
  // vertex is in the Fid in z-direction : z=abs(pos(3))<1610.
  bool isTrue = false;

  float *vertex = get_vertex(ireco);
	
  isTrue = ZvtxInFid(vertex);
	
  return isTrue;
}

bool OfficialEventParser::ZvtxInFid(float *vertex){

  bool isTrue = false;
	
  double z = TMath::Abs(vertex[2]);
  if (z < ZPINTK-200) isTrue = true;
	
  return isTrue;
}

float OfficialEventParser::cut_nhitac(){
  if (get_version() == 1) return 10.;      // nhitac<=9 for sk1
  else if (get_version() == 2) return 16.; // nhitac<=15 for sk2
  else if (get_version() == 3) return 16.; // nhitac<=15 for sk3
  else if (get_version() == 4) return 16.; // nhitac<=15 for sk4
  else { cerr << " Need to set nhitac cut in cut_nhitac()" << endl; exit(-1);}

}

float OfficialEventParser::cut_potot(){
  if (get_version() == 1) return 3000.;      //  for sk1
  else if (get_version() == 2) return 1500.; //  for sk2
  else if (get_version() == 3) return 1500.; //  for sk3
  else if (get_version() == 4) return 1500.; //  for sk4
  else { cerr << " Need to set nhitac cut in cut_potot()" << endl; exit(-1);}

}

int OfficialEventParser::cut_ehit( int i ){
  //cut value of ehit(i) in muedcy 
  // i==1 or 2 means just first appeared ehit(i) or second in code

  if ( i==1) { 
    if (get_version() == 1) return 60;      //  for sk1
    else if (get_version() == 2) return 30; //  for sk2
    else if (get_version() == 3) return 60; //  for sk3
    else if (get_version() == 4) return 60; //  for sk4
    else { cerr << " Need to set ehit cut in cut_ehit()" << endl; exit(-1);}
  } else if (i==2) {
    if (get_version() == 1) return 40;      //  for sk1
    else if (get_version() == 2) return 20; //  for sk2
    else if (get_version() == 3) return 40; //  for sk3
    else if (get_version() == 4) return 40; //  for sk4
    else { cerr << " Need to set ehit cut in cut_ehit()" << endl; exit(-1);}
  }
}

int OfficialEventParser::muedcy(int ireco){
  int muedcy;
  int nmuemax;

  /* temporary because SK-IV is not ready  */
  if ( get_version() == 4) {
    if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit)
      //return E->nmue;
      return E->stnmue;
    else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
      
      return E->fqnse - 1;
      //return E->nse - 1;

      //int nUnmatchedPeaks = 0; // Including parent peak
      //for (int ipeak=0; ipeak<E->muechk_ncand[imethod]; ipeak++) {
      //	if (E->muechk_icluster[imethod]<=0) nUnmatchedPeaks++; 
      //} 
      //return E->nse - nUnmatchedPeaks;
      
      //return E->muechk_ncand[0] - 1;

    }
  }
	
  if ( E->nmue <= 0 ) return 0;

  nmuemax = E->nmue;
  muedcy = 0;
  if ( E->nmue >= 10 ) nmuemax = 10;
	
  for ( int i=0; i< nmuemax ; i++){
    if ( E->evis < 1330.  && E->etime[i] < 0.1 ) continue;
    if ( E->evis >= 1330. && E->etime[i] < 1.2 ) continue;
    if ( E->etime[i] > 0.8 && E->etime[i] < 1.2 ) continue;
    if ( E->etype[i] == 1 && E->ehit[i] >= cut_ehit(1) && E->egood[i] > 0.5  )  { muedcy++; continue; }
    if ( (E->etype[i] == 2 || E->etype[i] == 3 || E->etype[i]== 4 ) && E->ehit[i]>=cut_ehit(2) ) { muedcy++; continue; }
  }

  return muedcy;

}

int OfficialEventParser::get_num_ring(int ireco){
  
  int rings = 0;
  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) rings = E->nring;
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    if (is1R(ireco)) rings = 1;
    else rings = 2;
  }
  return rings;
}

float OfficialEventParser::get_evis(int ireco){
  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) return E->evis;
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) return E->fq1rmom[0][ie];
}

float OfficialEventParser::get_potot(int ireco){

  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) return E->stpotot;
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) return E->potot;
}

float *OfficialEventParser::get_vertex(int ireco, int ise){

  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    if (ise==0) return E->stposmu; //E->pos;
    else return E->stpose; //E->epos[ise-1];
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    int i1rPID = get1rType(ireco,ise);
    return E->fq1rvtx[ise][i1rPID];
  }
  
}

void *OfficialEventParser::get_eststop_vertex(int ireco, int ise, float *stop_vertex){
  
  float *vertex = get_vertex(ireco,ise);
  float mom = get_1Rmom(ireco,ise);
  float *dir = get_1Rdir(ireco,ise);
  
  for (int i=0; i<3; i++) stop_vertex[i] = vertex[i] + (mom/2.3)*dir[i];

}

float *OfficialEventParser::get_1Rdir(int ireco, int ise){

  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    //if (ise==0) return E->dir[0];
    //else return E->edir[ise-1];    
    if (ise==0) return E->stdirmu;
    else return E->stdire;
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    int i1rPID = get1rType(ireco,ise);
    return E->fq1rdir[ise][i1rPID];
  }

}

void OfficialEventParser::get_2Rdir(int ireco, float *dir0, float *dir1){

  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    return;
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {

    for (int i=0; i<3; i++) {
      dir0[i] = E->fq2rdir1[MrPID[1][0]][MrPID[1][1]][i];
      dir1[i] = E->fq2rdir2[MrPID[1][0]][MrPID[1][0]][i];
    }
  }
  
  return;

}

float OfficialEventParser::get_1Rmom(int ireco, int ise){

  float mom = -1;

  if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    if (ise==0) {
    //  if (E->ip[0]==2) mom = E->amome[0];
      mom = E->stamom;
    } else {
    //  else if (E->ip[0]==3) mom = E->amomm[0];
      mom = E->stamomdcye;
    }
    
  }
  else if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
    
    int i1rPID = get1rType(ireco, ise);
    mom = E->fq1rmom[ise][i1rPID];	
  }

  if (mom!=mom) {
    //cout << "get_1Rmom() Warning: momentum = " << mom << ", setting to -1" << endl;
    mom = -1;
  }

  return mom;

}

double OfficialEventParser::get_fromwall( float *vertex, float *dir ){
  double from_wall;

  set_from_wall_flag(true);
  // inverse the ring direction, calculate "towall" instead
  from_wall = get_towall(vertex, dir);
  set_from_wall_flag(false);// reset the flag

  return from_wall;
}

double OfficialEventParser::get_towall(float *vertex, float *dir){
  /*
    Calculate towall for nth ring in the event.
    Original :: $ATMPD_ROOT/src/analysis/loe/lelib/towall.F
  */
  double x1, y1, z1;
  double x, y, z;
  double dx, dy, dz;
  double a0, b0, c0;
  double cz;
  double towall;
  double  rcir;
  double  hpos[3]; // 
  bool   might_go_to_side;
  double hanbe;
  double cz1;
  
  // get_towall_dir() can obtain where the ring direction goes.
  //  0 : out of tank,  1 : photon goes to barrel,   2 : photon goes to top,  3 : photon goes to bottom

  // initialize
  set_towall_dir(out_of_tank); // vertex is out of the tank
  might_go_to_side = false;
  cz = 0.;
  cz1 = 0.;
  towall = -999.; // initialize
  rcir = 0.;
  hanbe = 0.;

  // get vertex
  x1=vertex[0];	
  y1=vertex[1];	
  z1=vertex[2];	

  // get direction of the nth ring
  if ( get_from_wall_flag()) { // calclulate "from wall". So, inverse the ring direction,
    dx = -dir[0];
    dy = -dir[1];
    dz = -dir[2];
  } else { // default
    dx = dir[0];
    dy = dir[1];
    dz = dir[2];
  }

  //      IF(DZ) 12,20,11   --> negative : 12,   zero : 20,  plus : 11
  // the ring is down-going
  if (dz < 0.) {
    if (z1 < (double)ZMINTK) { // z vertex is at underneath of ID barrel.
      set_towall_dir(out_of_tank);
      return towall;
    } else {
      cz = (ZMINTK - z1) / (double)dz;
      x  = cz * dx + x1;
      y  = cz * dy + y1;
      z = (double)ZMINTK;
      rcir = TMath::Sqrt(x*x+y*y);

      if (rcir > (double)RINTK) {
	might_go_to_side = true;	
      } else {
	// photon goes to lower plate (bottom)
	set_towall_dir(to_bottom);
	towall = calc_towall(x,y,z, vertex);
      }
    }
  }  
  
  // the ring is upper-going
  if (dz > 0.){
    if ( z1 >= (double)ZPINTK){ // z vertex is higher than the top of the tank
      set_towall_dir(out_of_tank);
      return towall;
    } else {
      cz = ((double)ZPINTK - z1) / (double)dz;
      x  = cz * dx + x1;
      y  = cz * dy + y1;
      z  = (double)ZPINTK;
      rcir = TMath::Sqrt(x*x+y*y);

      if (rcir >= (double)RINTK){
	might_go_to_side = true;
      } else {
	// photon goes to upper plate (top)
	set_towall_dir(to_top);
	towall = calc_towall(x,y,z, vertex);
      }
    }
  }

  // the ring goes horizontally
  if (dz == 0. || might_go_to_side == true ) {
    double _rintk = (double)RINTK;
    a0 = dx*dx + dy*dy;
    b0 = dx*x1 + dy*y1;
    c0 = x1*x1 + y1*y1 - _rintk*_rintk; // RINTK*RINTK;
		
    hanbe = b0*b0 - c0*a0;

    if (hanbe < 0.0) {
      set_towall_dir(out_of_tank);
      return towall;
    }

    cz1 = (double)( -b0 + TMath::Sqrt(hanbe)) / (double)a0;

    if ( cz1 <= 0.0 ) {
      set_towall_dir(out_of_tank);
      return towall;
    }
		
    x = cz1 * dx + x1;
    y = cz1 * dy + y1;
    z = cz1 * dz + z1;

    if ( z >= ZPINTK || z <= ZMINTK) {
      set_towall_dir(out_of_tank);
      return towall;
    }

    set_towall_dir(to_side);
    towall = calc_towall(x,y,z,vertex);
  } 

  return towall;
}

double OfficialEventParser::calc_towall(double x, double y, double z, float *vertex){
  double towall;
  double hpos[3];
  
  hpos[0] = x;
  hpos[1] = y;
  hpos[2] = z;
  
  towall = TMath::Sqrt( (vertex[0] - hpos[0]) * (vertex[0] - hpos[0])
			+ (vertex[1] - hpos[1]) * (vertex[1] - hpos[1])
			+ (vertex[2] - hpos[2]) * (vertex[2] - hpos[2]));
  return towall;
}

void OfficialEventParser::rotateVector(float *vector, double theta, double phi) {
  
  float tmpvec[3];
  
  double COSTH = cos(theta);
  double SINTH = sin(theta);
  double COSPHI = cos(phi);
  double SINPHI = sin(phi);

  tmpvec[2] = COSTH*vector[2]  -  SINTH*vector[0];
  tmpvec[0] = SINTH*vector[2]  +  COSTH*vector[0];
  tmpvec[1] = vector[1];
  
  vector[0] = COSPHI*tmpvec[0] - SINPHI*vector[1];
  vector[1] = SINPHI*tmpvec[0] + COSPHI*vector[1];
  vector[2] = tmpvec[2];
}

double OfficialEventParser::weight(){
  return 1;
}

void OfficialEventParser::Print(){
}

void OfficialEventParser::Terminate(){
  
  /*
  int imc=0;
  for (int ireco=0; ireco<nRecos; ireco++) {
    
    int apr2009bin = h_cosmic_nmue_vs_year[imc][ireco]->GetXaxis()->FindBin(2009+3./12);
    double dataval = h_cosmic_nmue_vs_year[imc][ireco]->GetBinContent(apr2009bin, 2); // 1 decay-e at April 2009
    double mcval = h_cosmic_nmue_vs_year[1][ireco]->GetBinContent(apr2009bin, 2); // 1 decay-e at April 2009
    
    double scale = dataval/mcval;

    cout << ireco << " " << dataval << " " << mcval << " " << scale << endl;

    int nHists = h_vect_hists[imc][ireco].size();
    for (int ihist=0; ihist<nHists; ihist++) {
      h_vect_hists[imc][ireco][ihist]->Scale(scale);
    }
  }
  */

  for (int imc=0; imc<nSample; imc++) {
    for (int ireco=0; ireco<nRecos; ireco++) {
      if (efile[imc][ireco]) {
	efile[imc][ireco]->Write();
	efile[imc][ireco]->Close();
      }
    }
  }


}
