#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <math.h>


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TText.h>
#include <TLine.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TColor.h>
#include <TStyle.h>

using namespace std;

void set2Dcolor(int i);

bool drawNormalized = 0;
void drawTH1D(TH1 *h, TString args) {
  
  string hname = h->GetName();
  if (hname.find("mc")!=string::npos) args += "hist";

  if (!drawNormalized)
    h->Draw(args);
  else
    h->DrawNormalized(args);
}

bool doLifeTimeEff = 0;

int main(int argc, char *argv[]){

  gROOT->ProcessLine(".x ~/.rootlogon.C");  
  
  bool doPrint = 1;
  bool doLog = 0;
  int nrgSepType = 0;
  int maxE = 999999;
  string doHistoName = "";
  if (argc>1)
    doPrint = atoi(argv[1]);
  if (argc>2)
    doLog = atoi(argv[2]);
  if (argc>3)
    drawNormalized = atoi(argv[3]);
  if (argc>4)
    nrgSepType = atoi(argv[4]);
  if (argc>5)
    maxE = atoi(argv[5]);
  if (argc>6)
    doHistoName = argv[6];
    
    
  //cout << doPrint << " " << doLog << " " << drawNormalized << " " << nrgSepType << " " << doHistoName.c_str() << endl;

  string outputdir = Form("images_Elt%d/norm%d_log%d/nrg%d",maxE,drawNormalized,doLog,nrgSepType);

  gROOT->ProcessLine(Form(".! mkdir -p %s",outputdir.c_str())); 


  using namespace std;
  
  const int nMCs = 2;
  const int nRecos = 2;
  const int nSE = 2;
  
  TFile *file[nMCs][nRecos];
  
  char *h_name[nMCs] = {"h","hmc"};
  
  ofstream outtext(Form("%s/outtext.txt",outputdir.c_str()));

  double nEntriesForLifetimeNorm[nMCs][nRecos];

  for (int ireco=0; ireco<nRecos; ireco++) {
    for (int i=0; i<2; i++) {
      string filename = Form("my_cosmics_reco%d_nrg%d_Elt%d_%s.root",ireco, nrgSepType, maxE, h_name[i]);
      cout << filename.c_str() << endl;
      file[i][ireco] = new TFile(filename.c_str());

      nEntriesForLifetimeNorm[i][ireco] = ((TH1D*)file[i][ireco]->Get(Form("%s_cosmic_nmue_reco%d",h_name[i],ireco)))->Integral();
    }
  }


  const int nPlots = 4;
  string plotsToFind[nPlots] = {
    "lnLrat1R2R_vs_mom",
    "lnLrat1R2R_vs_lnLrat1R",
    "lnLrat1R_vs_mom_ise0",
    "lnLrat1R_vs_mom_ise1"
  };

  const int nCuts = 3;
  
  TH1 *tempH1[nMCs][nRecos][nPlots][nCuts] = {{{0}}};
    
  TIter next(file[0][0]->GetListOfKeys());
  TKey *key;
  while ((key=(TKey*)next())) {
    
    string strNameOrig  = key->GetName();
    string strClass = key->GetClassName();
    
    int iplot, icut;
    for (iplot=0; iplot<nPlots; iplot++) {
      
      if (strNameOrig.find(plotsToFind[iplot].c_str())==string::npos) continue;
      
      for (icut=0; icut<nCuts; icut++) {				
	if (strNameOrig.find(Form("cut%d",icut))!=string::npos)
	  break;
      }
      
      if (icut<nCuts) break;
      
    }
    
    cout << strNameOrig.c_str() << " " << iplot << " " << icut << endl;

    if (iplot>=nPlots) continue;
    
    int ise = -1;
    if (strNameOrig.find("_ise0_")!=string::npos) ise=0;
    else if (strNameOrig.find("_ise1_")!=string::npos) ise=1;

    for (int ireco=0; ireco<nRecos; ireco++) {
      
      string strName = strNameOrig;
      
      strName.replace(strName.find("reco0"), std::string("reco0").length(), Form("reco%d",ireco));
      
      for (int imc=0; imc<2; imc++) {
	
	strName.replace(strName.find("h_"), std::string("h_").length(), Form("%s_",h_name[imc]));
	
	tempH1[imc][ireco][iplot][icut] = (TH1*)(((TH1*)file[imc][ireco]->Get(strName.c_str()))->Clone());

	tempH1[imc][ireco][iplot][icut]->SetTitle("");
	 
	if (drawNormalized) 
	  tempH1[imc][ireco][iplot][icut]->Scale(1/tempH1[imc][ireco][iplot][icut]->Integral());
	
	cout << tempH1[imc][ireco][iplot][icut]->GetName() << endl;
      }
    }
  }
   
  for (int ireco=0; ireco<nRecos; ireco++) {

    for (int icut=0; icut<nCuts; icut++) {	
     
      for (int imc=0; imc<2; imc++) {
      	TCanvas *c_reco = new TCanvas(1);
	c_reco->Divide(2,2);
	
	for (int iplot=0; iplot<nPlots; iplot++) {
	  
	  if (doLog) c_reco->cd(iplot+1)->SetLogz(1);
	  c_reco->cd(iplot+1);
	  
	  tempH1[imc][ireco][iplot][icut]->Draw("COLZ");
	  
	  
	}
	
	if (doPrint) {
	  c_reco->Print(Form("%s/all_likelihoods_reco%d_imc%d_icut%d.png",outputdir.c_str(),ireco,imc,icut));
	  //c_reco->Print(Form("%s/all_likelihoods_reco%d_imc%d_icut%d.pdf",outputdir.c_str(),ireco,imc,icut));
	}
      }
      
      // 2D ratio
      for (int imc=0; imc<2; imc++) {
	for (int iplot=0; iplot<nPlots; iplot++) {
	  ((TH2D*)tempH1[imc][ireco][iplot][icut])->Rebin2D(2,2);
	  
	  if (iplot<=2) 
	    ((TH2D*)tempH1[imc][ireco][iplot][icut])->Rebin2D(8,8);
	  
	}
      }
      
      // Set binary color scale
      set2Dcolor(0);
      
      TCanvas *c_reco_rat = new TCanvas(1);
      c_reco_rat->Divide(2,2);
      
      for (int iplot=0; iplot<nPlots; iplot++) {
	c_reco_rat->cd(iplot+1);

	int imc=0;
	TH2D *tempH1_rat = (TH2D*)(tempH1[imc][ireco][iplot][icut])->Clone();
	tempH1_rat->Reset();
	tempH1_rat->Divide(tempH1[1][ireco][iplot][icut],tempH1[0][ireco][iplot][icut]);
	
	tempH1_rat->SetMinimum(0);
	tempH1_rat->SetMaximum(2);
	tempH1_rat->Draw("COLZ");
      
      
	if (doPrint) {
	  c_reco_rat->Print(Form("%s/all_likelihoods_reco%d_cut%d_rat.png",outputdir.c_str(),ireco,icut));
	  //c_reco_rat->Print(Form("%s/all_likelihoods_reco%d_cut%d_rat.pdf",outputdir.c_str(),ireco,icut));
	}
      }
      // Reset color scale
      set2Dcolor(1);
	
    }
      
  }
}


void set2Dcolor(int i) {

  if (i) {
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  else {
    const Int_t NRGBs = 3;
    const Int_t NCont = 255;
    
    Double_t stops[NRGBs] = { 0.00, 0.5, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 1.0, 1.00 };
    Double_t green[NRGBs] = { 0.00, 1.0, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 1.0, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
}
