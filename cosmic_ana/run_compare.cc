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
#include <TMath.h>
#include <TFitResult.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TText.h>
#include <TLatex.h>
#include <TLine.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TColor.h>
#include <TStyle.h>
#include <TArrow.h>

#include "geotank.h"

const int nMonths = 12;
TString monthNames[nMonths] = {"Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

using namespace std;

double get_res1sigma(TH1 * h1tmp);
void set2Dcolor(int i);

bool drawNormalized = 0;
void drawTH1D(TH1 *h, TString args) {
  
  string hname = h->GetName();
  if (hname.find("mc")!=string::npos) args += "hist";

  //if (!drawNormalized)
    h->Draw(args);
    //else
    //h->DrawNormalized(args);
}

bool doLifeTimeEff = 0;

int font=132;

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
  const int nRecos = 1;
  const int nSE = 2;

  const int nPars = 5;
  const TString parNames[nPars] = {"gausmean","sigma", "mean", "rms", "mispid"};
  TString parTitles[nPars] = {"Gaus. Mean", "#sigma", "Mean", "RMS","Mis-PID Rate (%)"};
  enum parEnum {igausmean, isigma, imean, irms, imispid};

  TFile *file[nMCs][nRecos];
  
  char *h_name[nMCs] = {"h","hmc"};
  
  ofstream outtext(Form("%s/outtext.txt",outputdir.c_str()));

  double nEntriesForLifetimeNorm[nMCs][nRecos];

  double normToData[nMCs][nRecos];
  double dataScale[nRecos];

  TFile *outfile = 0;
  if (!doPrint) outfile = new TFile("pars_vs_mom.root","RECREATE");

  for (int ireco=0; ireco<nRecos; ireco++) {
    for (int i=0; i<2; i++) {
      string filename = Form("%s_reco%d_nrg%d_Elt%d_file-1.root", h_name[i], ireco, nrgSepType, maxE);
      cout << filename.c_str() << endl;
      file[i][ireco] = new TFile(filename.c_str());

      //nEntriesForLifetimeNorm[i][ireco] = ((TH1D*)file[i][ireco]->Get(Form("%s_cosmic_nmue_reco%d",h_name[i],ireco)))->Integral();
      
      //TH2D *htmp = (TH2D*)file[i][ireco]->Get(Form("%s_cosmic_nmue_vs_year_reco%d",h_name[i],ireco));
      //int apr2009bin = htmp->GetXaxis()->FindBin(2009.+3/12.);
      //if (i) normToData[i][ireco] = htmp->ProjectionX("_pxfornorm",2,2)->Integral();
      //else normToData[i][ireco] = htmp->GetBinContent(apr2009bin,2);

      TH1D *htmp = (TH1D*)file[i][ireco]->Get(Form("%s_cosmic_nmue_reco%d",h_name[i],ireco));
      normToData[i][ireco] = htmp->Integral();
      cout << ireco << " " << i << " " << normToData[i][ireco] << endl;
    }

    dataScale[ireco] = normToData[0][ireco]/normToData[1][ireco];
    cout << dataScale[ireco] << endl;
  }

  TH1 *tempH1[nMCs][nRecos] = {{0}};


  TLatex *text;
  TLine *line;

  TIter next(file[0][0]->GetListOfKeys());
  TKey *key;
  int nPlots=0;
  while ((key=(TKey*)next())) {
    
    string strNameOrig  = key->GetName();
    string strClass = key->GetClassName();

    //if (strNameOrig.find("momsep")!=string::npos)
    //if (strNameOrig.find("lifetime")==string::npos)
    //if (strNameOrig.find("zenith")==string::npos)
    //if (strNameOrig.find("_res_")==string::npos)
    //if (strNameOrig.find("mom_over_range")==string::npos)
    if (strNameOrig.find(doHistoName.c_str())==string::npos)
      continue;
    
    //cout << strNameOrig.c_str() << endl;

    int ise = -1;
    if (strNameOrig.find("_ise0_")!=string::npos) ise=0;
    else if (strNameOrig.find("_ise1_")!=string::npos) ise=1;

    for (int ireco=0; ireco<nRecos; ireco++) {
      
      string strName = strNameOrig;
      
      strName.replace(strName.find("reco0"), std::string("reco0").length(), Form("reco%d",ireco));
      
      for (int i=0; i<2; i++) {
	
	strName.replace(strName.find("h_"), std::string("h_").length(), Form("%s_",h_name[i]));
	
	//cout << strName.c_str() << endl;
	tempH1[i][ireco] = (TH1*)(((TH1*)file[i][ireco]->Get(strName.c_str()))->Clone());

	tempH1[i][ireco]->SetTitle("");

	if (drawNormalized || 
	    strNameOrig.find("_res_")!=string::npos
	    ) {
	  tempH1[i][ireco]->Scale(1/tempH1[i][ireco]->Integral());
	  
	}
	else if (i==1) {
	  tempH1[1][ireco]->Scale(dataScale[ireco]);
	}
      }

    }

    TCanvas *c = new TCanvas(1);
    if (doLog) c->SetLogy(1);

    if ( !strClass.compare("TH1D") )  {
      
      
      //TCanvas *c = new TCanvas(Form("c%d",nPlots),Form("c%d",nPlots),0,0,800,600);

      double maxY = 0;

      TText textVal[nMCs][nRecos];
      double mean[nMCs][nRecos], meanerr[nMCs][nRecos];
      double err[nMCs][nRecos], errerr[nMCs][nRecos];
      
      TLine vertLine;
      vertLine.SetLineColor(kAzure-8);
      vertLine.SetLineWidth(3);

      for (int ireco=0; ireco<nRecos; ireco++) {
		
	for (int i=0; i<2; i++) {

	  textVal[i][ireco].SetNDC(true);
	  textVal[i][ireco].SetTextSize(0.05);
	  textVal[i][ireco].SetTextFont(font);
	  
	  if (!doLog) {
	    
	    if ( (strNameOrig.find("zenith_ise1_reco0")==string::npos) )
	      tempH1[i][ireco]->SetMinimum(0);
	    
	  }
	  else
	    tempH1[i][ireco]->SetMinimum(1);

	  if (drawNormalized) tempH1[i][ireco]->GetYaxis()->SetTitle("Area Normalized");
	  else tempH1[i][ireco]->GetYaxis()->SetTitle("Number of Events");
	  
	  // Control binning
	  if ( (strNameOrig.find("z_ise0")!=string::npos) ) 
	  tempH1[i][ireco]->Rebin(2);

	  else if ( (strNameOrig.find("_wall_")!=string::npos) ) {
	    //tempH1[i][ireco]->Rebin(4);
	  }
	  else if ( (strNameOrig.find("zenith_ise1")!=string::npos) ) {
	    tempH1[i][ireco]->Rebin(4);
	  }


	  // Control x-axis
	  if ( (strNameOrig.find("r2")!=string::npos) ) {
	    //tempH1[i][ireco]->GetXaxis()->SetRangeUser(2200000, 3500000);
	    tempH1[i][ireco]->GetXaxis()->SetRangeUser(0, 3500000);
	  } 
	  else if ( (strNameOrig.find("z_ise0")!=string::npos) ) {
	    tempH1[i][ireco]->GetXaxis()->SetRangeUser(1650, 2000);
	  } 
	  else if ( (strNameOrig.find("cosmic_mom_ise0_reco")!=string::npos) ) {
	    tempH1[i][ireco]->GetXaxis()->SetRangeUser(100, 10000);
	  }
	  else if (strNameOrig.find("_ddir_")!=string::npos) {
	    if (!doLog) tempH1[i][ireco]->GetXaxis()->SetRangeUser(0, 15);
	  }
	  else if ( (strNameOrig.find("_ise0_zzoom_")!=string::npos) ) {
	    tempH1[i][ireco]->GetXaxis()->SetRangeUser(1650, 1950);
	  } 
	  else if ( (strNameOrig.find("_ise0_r2zoom_")!=string::npos) ) {
	    tempH1[i][ireco]->GetXaxis()->SetRangeUser(2506100, 3206100);
	  } 
	  else if ( (strNameOrig.find("_ise0_r2zoom_")!=string::npos) ) {
	    tempH1[i][ireco]->GetXaxis()->SetRangeUser(2506100, 3206100);
	  } 
	  else if ( (strNameOrig.find("_ise0_wall_")!=string::npos) ) {
	    tempH1[i][ireco]->GetXaxis()->SetRangeUser(-150, 150);
	  } 
	  else if ( (strNameOrig.find("_ise0_Le_Lmu_")!=string::npos) ) {
	    if ( (strNameOrig.find("_momsep0_500")!=string::npos) || 
		 (strNameOrig.find("_momsep500_1000")!=string::npos) ) { 
	      tempH1[i][ireco]->GetXaxis()->SetRangeUser(-1400, 500);
	    }
	  } 	      
	  else if ( (strNameOrig.find("_ise1_Le_Lmu_")!=string::npos) ) {
	    if ( (strNameOrig.find("_momsep40_45")!=string::npos) )
	      tempH1[i][ireco]->GetXaxis()->SetRangeUser(-150, 500);
	  } 

	  // Resolution (1-sigma from 0)
	  if (strNameOrig.find("_res_")!=string::npos) {
	    
	    tempH1[i][ireco]->GetYaxis()->SetTitle("Area Normalized");
	    
	    mean[i][ireco] = get_res1sigma(tempH1[i][ireco]);
	    err[i][ireco] = 0;
	    textVal[i][ireco].SetText(0.6, 0.8-i*0.1, Form("%.2f", mean[i][ireco]));
	    
	    if (strNameOrig.find("_dir_")!=string::npos) 
	      tempH1[i][ireco]->SetTitle(Form("<1#sigma: %.1f#circ", mean[i][ireco]));
	    else
	      tempH1[i][ireco]->SetTitle(Form("<1#sigma: %.1f cm", mean[i][ireco]));
	    
	    if (strNameOrig.find("_ddir_")!=string::npos) {
	      if (i==0) tempH1[i][ireco]->SetTitle(Form("<1#sigma: Data = %.3f#circ", mean[i][ireco]));
	      else {
		tempH1[i][ireco]->SetTitle(Form("%s, MC = %.3f#circ",tempH1[0][ireco]->GetTitle(), mean[i][ireco]));
		tempH1[0][ireco]->SetTitle(tempH1[i][ireco]->GetTitle());
	      }
	      tempH1[0][ireco]->SetTitle("");
	      tempH1[1][ireco]->SetTitle("");
	      tempH1[i][ireco]->Rebin(10);
	      if (!doLog) tempH1[i][ireco]->GetXaxis()->SetRangeUser(0, 12);
	    }

	    outtext << tempH1[i][ireco]->GetName() << " res = " << mean[i][ireco] << endl;

	  }

	  // (Mean and RMS values)
	  if  ( (i && 
		 (strNameOrig.find("_lvtx_res_")!=string::npos || 
		  strNameOrig.find("_range_res_")!=string::npos) 
		 ) ) {
	    
	    mean[i][ireco] = tempH1[i][ireco]->GetMean();
	    err[i][ireco] = tempH1[i][ireco]->GetRMS();

	    meanerr[i][ireco] = tempH1[i][ireco]->GetMeanError();
	    errerr[i][ireco] = tempH1[i][ireco]->GetRMSError();
	    
	    outtext << tempH1[i][ireco]->GetName() << " mean/rms = " << mean[i][ireco] << " " << err[i][ireco] << endl;
	    
	    textVal[i][ireco].SetText(0.6, 0.8-i*0.1, Form("%.2f %.2f", mean[i][ireco], err[i][ireco]));
	    tempH1[i][ireco]->SetTitle(Form("Mean: %.2f #pm %.2f cm, RMS: %.2f #pm %.2f cm", mean[i][ireco], meanerr[i][ireco], err[i][ireco], errerr[i][ireco]));
	    
	  } 
	  
	  else if (strNameOrig.find("_mom_ise1_")!=string::npos) { 
	    
	    mean[i][ireco] = tempH1[i][ireco]->GetMean();
	    if (sqrt(tempH1[i][ireco]->GetEntries()))
		err[i][ireco]  = tempH1[i][ireco]->GetRMS()/sqrt(tempH1[i][ireco]->GetEntries());

	    textVal[i][ireco].SetText(0.6, 0.8-i*0.1, Form("%.3f %.3f", mean[i][ireco], err[i][ireco]));
	  }

	  // Likelihood Ratios
	  if (i && strNameOrig.find("_lnLrat1R")!=string::npos)
	    tempH1[i][ireco]->Scale(tempH1[0][ireco]->Integral()/tempH1[i][ireco]->Integral());


	  maxY = tempH1[i][ireco]->GetMaximum() > maxY ? tempH1[i][ireco]->GetMaximum() : maxY;
	  
	  
	  //tempH1[i][ireco]->SetLineColor(ireco+1);
	  //tempH1[i][ireco]->SetMarkerColor(ireco+1);
	  tempH1[i][ireco]->SetLineColor(i+1);
	  tempH1[i][ireco]->SetMarkerColor(i+1);
	  if (i==0 && ireco==1) {
	    tempH1[i][ireco]->SetLineColor(kBlack);
	    tempH1[i][ireco]->SetMarkerColor(kBlack);
	    tempH1[i][ireco]->SetMarkerStyle(20);
	  }
	  tempH1[i][ireco]->SetMarkerSize(1.2);

	  // Likelihood projection
	  if (strNameOrig.find("_lnLrat1R")!=string::npos) {
	    tempH1[0][ireco]->SetLineColor(1);
	    tempH1[0][ireco]->SetMarkerColor(1);

	    tempH1[1][ireco]->SetLineColor(2);
	    tempH1[1][ireco]->SetMarkerColor(2);
	  }
		

	  if (i==1)
	    tempH1[i][ireco]->SetLineWidth(2);
	  
	}
      }

      double maxScale = 1.1;
      if (doLog) maxScale = 2;
      
      TF1 *f[nMCs][nRecos];
      for (int ireco=0; ireco<nRecos; ireco++) {
	TCanvas *c_reco = new TCanvas(1);
	if (doLog) c_reco->SetLogy(1);

	if ( (strNameOrig.find("mom_ise0")!=string::npos) ) {
	  c_reco->SetLogx(0);
	} 
	
	TCanvas *c_mcdata[2] = {0};

	bool isDrawn = 0;
	for (int i=0; i<2; i++) {
	//for (int i=1; i>=0; i--) {
	  	  
	  if (!tempH1[i][ireco]->GetEntries()) continue;
	 
	  //tempH1[i][ireco]->SetMaximum(TMath::Max( TMath::Max(tempH1[0][0]->GetMaximum(),tempH1[1][0]->GetMaximum()),TMath::Max(tempH1[0][1]->GetMaximum(),tempH1[1][1]->GetMaximum()) )*maxScale );

	  // Lifetime
	  if (strNameOrig.find("lifetime")!=string::npos) {

	    if (!doLog) {
	      tempH1[i][ireco]->GetXaxis()->SetRangeUser(0,3);
	    }
	    else 
	      tempH1[i][ireco]->Rebin(10);
	    if (!doLifeTimeEff) {
	      
	      //f[i][ireco] = new TF1(Form("f_reco%d_i%d",ireco,i),"[0]*exp(-[1]*x)",1.5,8);
	      //tempH1[i][ireco]->Fit(f[i][ireco],"NRQ");
	      
	      //if (f[i][ireco]->GetParameter(1) && f[i][ireco]->GetParError(1)==f[i][ireco]->GetParError(1)) {
	      //double lifetime = 1/f[i][ireco]->GetParameter(1);
	      //double lifetime_err = (1/f[i][ireco]->GetParameter(1)) * (f[i][ireco]->GetParError(1)/f[i][ireco]->GetParameter(1));
	      //outtext << tempH1[i][ireco]->GetName() << " " << lifetime << " ± " << lifetime_err << endl;
	      
	      //textVal[i][ireco].SetText(0.6, 0.8-i*0.2, Form("%.3f %.3f",lifetime, lifetime_err));
	      
	      if (doPrint) {
		if (!isDrawn) {
		  drawTH1D(tempH1[i][ireco],"");
		  isDrawn = 1;
		}
		else {
		  drawTH1D(tempH1[i][ireco],"SAME");
		}
	      }
	      //f[i][ireco]->Draw("SAME");
	      //if (doPrint) textVal[i][ireco].Draw("SAME");
	      //}
	    }
	    
	    
	    else {
	      tempH1[i][ireco]->Rebin(2);

	      double maxT = 25;
	      f[i][ireco] = new TF1(Form("f_reco%d_i%d",ireco,i),"[0]*(0.578*exp(-x/2.197) + 0.422*exp(-x/1.795))",0,maxT);
	      f[i][ireco]->SetParameter(0,1);
	      double funcIntegral = f[i][ireco]->Integral(0,maxT);
	      f[i][ireco]->SetParameter(0,nEntriesForLifetimeNorm[i][ireco]/funcIntegral*tempH1[i][ireco]->GetBinWidth(1));	
	    
	      double eff = 100*tempH1[i][ireco]->Integral()/f[i][ireco]->Integral(0,maxT)*tempH1[i][ireco]->GetBinWidth(1);
	      double eff_err = eff/sqrt(tempH1[i][ireco]->GetEntries());
	      textVal[i][ireco].SetText(0.6, 0.8, Form("%.3f %.3f",eff, eff_err));
	      
	      outtext << tempH1[i][ireco]->GetName() << " " << eff << " ± " << eff_err << endl;
	      	      
	      c_mcdata[i] = new TCanvas(1);
	      if (doLog) c_mcdata[i]->SetLogy(1);
	    	    
	      tempH1[i][ireco]->SetMaximum(f[i][ireco]->Eval(0)*1.05);
	      tempH1[i][ireco]->GetXaxis()->SetRangeUser(0,10);
	    
	      tempH1[i][ireco]->SetFillColor(2*((i+1)%2)+2);
	      tempH1[i][ireco]->SetLineColor(2*((i+1)%2)+2);
	      
	      if (doPrint) {
		tempH1[i][ireco]->Draw("hist");
		f[i][ireco]->Draw("SAME");
		textVal[i][ireco].Draw("SAME");
	      }

	    }
	    
	  } 

	  else {
	    
	    if (doPrint) {
	      if (!isDrawn) {
		drawTH1D(tempH1[i][ireco],"");
		isDrawn = 1;
		
		if (strNameOrig.find("_res_")!=string::npos) tempH1[i][ireco]->SetMinimum(1e-5);

	      }
	      else {
		drawTH1D(tempH1[i][ireco],"SAME");
		
	      }
	    }
	   
	    
	    // Control cut line
	    //vertLine.Set(1330,0,1330,tempH1[0][ireco]->GetYaxis()->GetXmax());
	    if ( (strNameOrig.find("mom_ise0")!=string::npos) ) {
	      vertLine.DrawLine(1330,0,1330,tempH1[i][ireco]->GetYaxis()->GetXmax());
	    }
	  
	    
	  }
	  
	  
	  // Things that need to be fit by a Gaussian
	  int gausDist = -1;
	  double gausRange[2];
	  if (i && strNameOrig.find("_mom_res_")!=string::npos) gausDist = 0;
	  else if (strNameOrig.find("mom_over_range")!=string::npos) gausDist = 1;
	  else if (strNameOrig.find("_nearwall_")!=string::npos) gausDist = 2;
	  else if (strNameOrig.find("_emom_")!=string::npos) gausDist = 3;

	  if (gausDist>=0) {
	    double gausMean = tempH1[i][ireco]->GetMean();
	    double gausRMS = tempH1[i][ireco]->GetRMS();
	    gausRange[0]=gausMean-gausRMS*(1+gausDist); 
	    gausRange[1]=gausMean+gausRMS*(1+gausDist/4);

	    //cout << "range = " << i << " " << ireco << " " << gausRange[0] << " " << gausRange[1] << endl;
	    
	    f[i][ireco] = new TF1(Form("f_reco%d_i%d",ireco,i),"gaus",gausRange[0],gausRange[1]);
	    tempH1[i][ireco]->Fit(f[i][ireco],"NRQ");
	    
	    mean[i][ireco] = f[i][ireco]->GetParameter(1);
	    meanerr[i][ireco] = f[i][ireco]->GetParError(1);
	    if (gausDist==0) { // momentum resolution
	      err[i][ireco] = f[i][ireco]->GetParameter(2); 
	      errerr[i][ireco] = f[i][ireco]->GetParError(2); 
	    }
	    //else if (gausDist==1) err[i][ireco] = f[i][ireco]->GetParError(1);   // momentum/range
	    else if (gausDist>=1) {
	      mean[i][ireco] = tempH1[i][ireco]->GetMean();
	      meanerr[i][ireco] = tempH1[i][ireco]->GetMeanError();
	      err[i][ireco] = tempH1[i][ireco]->GetRMS();
	      errerr[i][ireco] = tempH1[i][ireco]->GetRMSError();

	      if (strNameOrig.find("_r2_")!=string::npos) {
		mean[i][ireco]    = sqrt(mean[i][ireco]   ); 
		meanerr[i][ireco] = sqrt(meanerr[i][ireco]);
		err[i][ireco]     = sqrt(err[i][ireco]    ); 
		errerr[i][ireco]  = sqrt(errerr[i][ireco] ); 
	      }
	    }

	    outtext << tempH1[i][ireco]->GetName() << endl;
	    outtext << "Mean: " << mean[i][ireco] << " +/- " <<  meanerr[i][ireco]  << endl; 
	    outtext << "RMS(or sigma): " << err[i][ireco] << " +/- " <<  errerr[i][ireco]  << endl; 
	      
	    //if (f[i][ireco]->GetParameter(1) && f[i][ireco]->GetParError(1)==f[i][ireco]->GetParError(1)) {
	    //cout << tempH1[i][ireco]->GetName() << " " << mean[i][ireco] << " ± " << f[i][ireco]->GetParameter(2) << " " << err[i][ireco]/f[i][ireco]->GetParameter(1)*100 << endl;
	    
	    //cout << tempH1[i][ireco]->GetName() << " " << mean[i][ireco] << " ± " << err[i][ireco] << endl; 
	    
	    textVal[i][ireco].SetText(0.6, 0.8-i*0.1, Form("Mean: %.3f, #sigma: %.3f", mean[i][ireco], err[i][ireco]));
	    if (gausDist==0) tempH1[i][ireco]->SetTitle(Form("Mean: %.2f #pm %.2f%%, #sigma: %.2f #pm %.2f%%", mean[i][ireco], meanerr[i][ireco], err[i][ireco], errerr[i][ireco]));
	    
	    if (doPrint) 
	      f[i][ireco]->Draw("SAME");
	    if (doPrint) textVal[i][ireco].Draw("SAME");
	    
	  }

	  
	  // Resolution
	  //else if ( (i && strNameOrig.find("_res_")!=string::npos) ||
	  //	    (strNameOrig.find("_ddir_")!=string::npos) ||
	  //	    (strNameOrig.find("_mom_ise1_")!=string::npos) ) {
	  //  //if (doPrint) textVal[i][ireco].Draw("SAME");
	  //}
	  
	  // Zenith 
	  //else if (strNameOrig.find("_zenith_")!=string::npos) {
	  //  
	  //  f[i][ireco] = new TF1(Form("f_reco%d_i%d",ireco,i),"[0]*x+[1]",-1,.9);
	  //  
	  //  tempH1[i][ireco]->Fit(f[i][ireco],"NR");
	  //  f[i][ireco]->Draw("SAME");
	  //
	  //  cout << f[i][ireco]->Eval(-1) << " " << f[i][ireco]->Eval(1) << " " << f[i][ireco]->Eval(-1)/f[i][ireco]->Eval(1) << endl; 
	  //}
	}
	outtext << endl;
	
	// Calculate differences
	TText textValDiff;
	textValDiff.SetNDC(true);
	textValDiff.SetTextSize(0.05);
	textValDiff.SetTextFont(font);
	int i=0;
	if ( (strNameOrig.find("mom_over_range")!=string::npos) ||
	     (strNameOrig.find("_mom_ise1_")!=string::npos) ||
	     (strNameOrig.find("_ddir_")!=string::npos) ) {
	
	  //double diff = 100*(mean[1][ireco]-mean[0][ireco])/mean[0][ireco];
	  //double errdiff = 100*sqrt(pow(err[0][ireco],2)+pow(err[1][ireco],2))/mean[0][ireco];
	  //textValDiff.SetText(0.6, 0.6, Form("%.3f %.3f", diff, errdiff));
	  ////if (doPrint) textValDiff.Draw("SAME");
	  //
	  //outtext << tempH1[i][ireco]->GetName() << " " << mean[1][ireco] << " ± " << err[1][ireco] << " " <<  mean[0][ireco] << " ± " << err[0][ireco] << " " << diff << " ± " << errdiff << endl; 
		      
	  double diff = 100*(mean[1][ireco]-mean[0][ireco])/mean[0][ireco];
	  double errdiff = 100*sqrt(pow(meanerr[0][ireco],2)+pow(meanerr[1][ireco],2))/mean[0][ireco];
	  textValDiff.SetText(0.6, 0.6, Form("%.3f %.3f", diff, errdiff));
	  //if (doPrint) textValDiff.Draw("SAME");
	  
	  outtext << tempH1[i][ireco]->GetName() << " " << mean[1][ireco] << " ± " << meanerr[1][ireco] << " " <<  mean[0][ireco] << " ± " << meanerr[0][ireco] << " " << diff << " ± " << errdiff << endl; 

	} 	
	
	
	
	if (doPrint) {
	  string strNamePrint = strNameOrig;
	  strNamePrint.replace(strNamePrint.find("h_"), 2, "");
	  strNamePrint.replace(strNamePrint.find("reco0"), 5, Form("reco%d",ireco));

	  //if (doLog) strNamePrint.append("_log");
	  
	  if (strNameOrig.find("lifetime")!=string::npos && doLifeTimeEff) {
	    for (int imc=0; imc<nMCs; imc++) {
	      if (c_mcdata[imc])
		c_mcdata[imc]->Print(Form("%s/%s_imc%d.png",outputdir.c_str(),strNamePrint.c_str(),imc));
	    }
	  }
	  else {
	    //c_reco->Print(Form("%s/%s.png",outputdir.c_str(),strNamePrint.c_str()));
	    c_reco->Print(Form("%s/%s.pdf",outputdir.c_str(),strNamePrint.c_str()));
	  }
	}
      }

      // Overlay all data/MC
      c->cd();
      bool isDrawn = 0;
      for (int ireco=0; ireco<nRecos; ireco++) {

	for (int i=0; i<2; i++) {
 	  tempH1[i][ireco]->SetMaximum(maxY*maxScale);
	  
	  if (!tempH1[i][ireco]->GetEntries()) continue;
	  
	  if (doPrint) {
	    if (!isDrawn) {
	      drawTH1D(tempH1[i][ireco],"");
	      isDrawn = 1;
	    }
	    else {
	      drawTH1D(tempH1[i][ireco],"SAME");
	    }
	  }

	}
	
	
  
	if (doPrint) {
	  string strNameAll = strNameOrig;
	  strNameAll.replace(strNameAll.find("h_"), 2, "");   
	  strNameAll.replace(strNameAll.find("reco0"), 5,"");
	  //if (doLog) strNameAll.append("_log");
	  //c->Print(Form("%s/%sall.png",outputdir.c_str(),strNameAll.c_str()));
	  c->Print(Form("%s/%sall.pdf",outputdir.c_str(),strNameAll.c_str()));
	}
      }

      
      isDrawn = 0;
      TCanvas *c_reco_rat = new TCanvas(1);
      c_reco_rat->cd();
      TH1D *tempH1_rat[2];
      for (int ireco=0; ireco<nRecos; ireco++) {
	
	// Ratio
	tempH1_rat[ireco] = (TH1D*)(tempH1[0][ireco])->Clone();
	tempH1_rat[ireco]->Reset();
	tempH1_rat[ireco]->Divide(tempH1[1][ireco],tempH1[0][ireco],1/tempH1[1][ireco]->Integral(),1/tempH1[0][ireco]->Integral());
	
	tempH1_rat[ireco]->SetMinimum(0.5);
	tempH1_rat[ireco]->SetMaximum(1.5);
	//tempH1_rat[ireco]->SetMaximum(2);
	if (doPrint) {
	  if (!isDrawn) {
	    tempH1_rat[ireco]->Draw("hist");
	    isDrawn = 1;
	  } else {
	    tempH1_rat[ireco]->Draw("hist same");
	  }
	}
      }
      
      if (doPrint) {
	string strNameAll = strNameOrig;
	strNameAll.replace(strNameAll.find("h_"), 2, "");   
	strNameAll.replace(strNameAll.find("reco0"), 5,"");
	//c_reco_rat->Print(Form("%s/%s_rat.png",outputdir.c_str(),strNameAll.c_str()));
	c_reco_rat->Print(Form("%s/%s_rat.pdf",outputdir.c_str(),strNameAll.c_str()));
      }
	
    }
      
    if ( !strClass.compare("TH2D") )  {
      
      // Cut lines
      TF1 *f_cut = 0;
      if (strNameOrig.find("_lnLrat1R_vs_mom_")!=string::npos) {
	f_cut = new TF1("f_cut","x*0.2",0,tempH1[0][0]->GetBinLowEdge(tempH1[0][0]->GetNbinsX()+1));
      }
      else if (strNameOrig.find("_lnLrat2R1R_")!=string::npos) {
      	//f_cut = new TF1("f_cut","150",0,1000);
	f_cut = new TF1("f_cut","pow(x/60.,2)+150",0,10000);
      }
      if (f_cut) f_cut->SetLineColor(kGray+1);
      
      for (int ireco=0; ireco<nRecos; ireco++) {
	
	for (int i=0; i<2; i++) {
	  
	  if (drawNormalized) {
	    
	    if (strNameOrig.find("z_vs_zenith")!=string::npos) {
	      ((TH2D*)tempH1[i][ireco])->Rebin2D(2,2);
	      
	      for (int jbin=1; jbin<=tempH1[i][ireco]->GetNbinsY(); jbin++) {
		double row_integral = 0;
		
		for (int ibin=1; ibin<=tempH1[i][ireco]->GetNbinsX(); ibin++) {
		  row_integral += tempH1[i][ireco]->GetBinContent(ibin,jbin);
		}
		
		// Row normalize
		for (int ibin=1; ibin<=tempH1[i][ireco]->GetNbinsX(); ibin++) 
		  if (row_integral)
		    tempH1[i][ireco]->SetBinContent(ibin,jbin, tempH1[i][ireco]->GetBinContent(ibin,jbin)/row_integral);
	      }
	    }
	  } 
	  
	  //else if (strNameOrig.find("cost_vs_mom")!=string::npos) {
	  else if (strNameOrig.find("_mom_res_vs_trumom_")!=string::npos ||
		   strNameOrig.find("_lnLrat1R_vs_mom_")!=string::npos) {

	    for (int ibin=1; ibin<=tempH1[i][ireco]->GetNbinsX(); ibin++) {
	  
	      double col_integral = 0;      
	      for (int jbin=1; jbin<=tempH1[i][ireco]->GetNbinsY(); jbin++) {
	  	col_integral += tempH1[i][ireco]->GetBinContent(ibin,jbin);
	      }
	      
	      // Column normalize
	      for (int jbin=1; jbin<=tempH1[i][ireco]->GetNbinsY(); jbin++) {
	  	if (col_integral)
	  	  tempH1[i][ireco]->SetBinContent(ibin,jbin, tempH1[i][ireco]->GetBinContent(ibin,jbin)/col_integral);
	      }
	    }
	    
	    if (strNameOrig.find("_ise1_")!=string::npos) 
	      tempH1[i][ireco]->GetXaxis()->SetRangeUser(0,65);
	    else 
	      ((TH2D*)tempH1[i][ireco])->Rebin2D(2,5);

	    tempH1[i][ireco]->GetZaxis()->SetTitle("Column Normalized");
	      //tempH1[i][ireco]->SetMaximum(0.065);


	  }      	  

	  gStyle->SetPadRightMargin(0.2);
	  
	  TCanvas *c_reco = new TCanvas(1);
	  if (doLog) c_reco->SetLogz(1);

	  if (doPrint) {
	    tempH1[i][ireco]->Draw("COLZ");
	    
	    if (f_cut) f_cut->Draw("SAME");
	  }
	  
	  // Time variation fits
	  if (strNameOrig.find("_vs_year_")!=string::npos) {
	    
	    gStyle->SetPadRightMargin(0.08);
	    gStyle->SetPadLeftMargin(0.16);

	    gStyle->SetStatX(0.91);
	    gStyle->SetStatY(0.88);
	    gStyle->SetStatW(0.3);
	    gStyle->SetStatH(0.20);
	    
	    if (i==1) continue;
	    if (ireco==0 && 
		(strNameOrig.find("_Le_")!=string::npos  || 
		 strNameOrig.find("_Lmu_")!=string::npos ||
		 strNameOrig.find("_Lpi0_")!=string::npos||
		 strNameOrig.find("_fq")!=string::npos   ||
		 strNameOrig.find("_cl")!=string::npos ) ) 
	      continue;

	    string histname = tempH1[i][ireco]->GetName();
	    //cout << "WTF 1 " << histname.c_str() << endl;

	    TH1D *hmc_par[nPars], *h_par[nPars];

	    for (int ipar=0; ipar<nPars; ipar++) {
	    
	      TString histTitle = Form("%s_%s_vs_year",tempH1[i][ireco]->GetName(),parNames[ipar].Data());
	      h_par[ipar] = ((TH2D*)tempH1[i][ireco])->ProjectionX(histTitle);
	      h_par[ipar]->Reset();
	      h_par[ipar]->SetTitle(histTitle);
	      if (ipar==imispid) {
		if (strNameOrig.find("lifetime")!=string::npos) 
		  h_par[ipar]->GetYaxis()->SetTitle("Lifetime from Single-Exp. (#mus)");
		else 
		  h_par[ipar]->GetYaxis()->SetTitle(Form("%s of %s",parTitles[ipar].Data(),tempH1[i][ireco]->GetYaxis()->GetTitle()));
	      }		
	      
	      h_par[ipar]->SetMarkerColor(kBlue);
	      h_par[ipar]->SetLineColor(kBlue);
	    }
	    
	    TCanvas *c_proj_year = new TCanvas(1);
	    if (doPrint) 
	      c_proj_year->Print(Form("%s/%s_projyear.pdf[",outputdir.c_str(),tempH1[i][ireco]->GetName()));

	    if (doLog) c_proj_year->SetLogy(1);

	    // MC
	    TH1D *hmc_projy = 0;
	    
	    double MCparVal[nPars] = {0};
	    double MCparErr[nPars] = {0};

	    //cout << "WTF 2 " << histname.c_str() << endl;

	    if (i==0) {
	      hmc_projy = ((TH2D*)tempH1[1][ireco])->ProjectionY(Form("%s_py",tempH1[1][ireco]->GetName()));

	      if (hmc_projy->GetEntries()) {
		if (drawNormalized) hmc_projy->Scale(1/hmc_projy->Integral());
	      
		for (int ipar=0; ipar<nPars; ipar++) {
		  TString histTitle = Form("%s_%s_vs_year",tempH1[1][ireco]->GetName(),parNames[ipar].Data());
		  hmc_par[ipar] = ((TH2D*)tempH1[1][ireco])->ProjectionX(histTitle);
		  hmc_par[ipar]->Reset();
		  hmc_par[ipar]->Rebin(hmc_par[ipar]->GetNbinsX());
		  hmc_par[ipar]->SetTitle(histTitle);
		  hmc_par[ipar]->GetYaxis()->SetTitle(Form("%s of %s",parTitles[ipar].Data(),tempH1[1][ireco]->GetYaxis()->GetTitle()));
		  hmc_par[ipar]->SetLineColor(kRed);
		  hmc_par[ipar]->SetMarkerColor(kRed);
		}
	      
		int maxBin = hmc_projy->GetMaximumBin();
		float maxX = hmc_projy->GetBinCenter(maxBin);
		float rms = hmc_projy->GetRMS();
	      
		float fitLow = maxX-rms;
		float fitHigh = maxX+rms;
		if (rms<=0) {
		  fitLow = maxX - maxX/3;
		  fitHigh = maxX + maxX/3;
		}

		MCparVal[imean] = hmc_projy->GetMean();
		MCparVal[irms] = rms;
		MCparErr[imean] = hmc_projy->GetMeanError();;
		MCparErr[irms] = hmc_projy->GetRMSError();

	    

		double leftIntegralErr, rightIntegralErr, fullIntegralErr;
		int zeroBin = hmc_projy->FindBin(1e-4);
		double leftIntegral = hmc_projy->IntegralAndError(0, zeroBin-1, leftIntegralErr);
		double rightIntegral = hmc_projy->IntegralAndError(zeroBin, hmc_projy->GetNbinsX()+1, rightIntegralErr);
		double fullIntegral = hmc_projy->IntegralAndError(0, hmc_projy->GetNbinsX()+1, fullIntegralErr);
	      
		if (strNameOrig.find("lifetime")!=string::npos) {
		  TF1 *f = new TF1(Form("f_reco%d_i%d",ireco,i),"[0]*exp(-x/[1])",1.2,10);
		  f->SetParameter(1,2);

		  hmc_projy->Fit(f,"RQ");
		
		  MCparVal[imispid] = f->GetParameter(1);
		  MCparErr[imispid] = f->GetParError(1);

		} else {
		  if (ise==1) {
		    MCparVal[imispid] = 100*(leftIntegral/fullIntegral);
		    MCparErr[imispid] = MCparVal[imispid]*sqrt(pow(leftIntegralErr/leftIntegral,2)+pow(fullIntegralErr/fullIntegral,2));
		  } else {
		    MCparVal[imispid] = 100*(rightIntegral/fullIntegral);
		    MCparErr[imispid] = MCparVal[imispid]*sqrt(pow(rightIntegralErr/rightIntegral,2)+pow(fullIntegralErr/fullIntegral,2));
		  }
		}
	      
		TFitResultPtr r = hmc_projy->Fit("gaus","NQRS","",fitLow, fitHigh);
		//r->Print("V");
	      
		//cout << "WTF 3 " << histname.c_str() << " " << r << " " << r->IsValid() << endl;
		if (r->IsValid()) {
		  MCparVal[igausmean] = r->Parameter(1);
		  MCparVal[isigma] = r->Parameter(2);
		  MCparErr[igausmean] = r->ParError(1);
		  MCparErr[isigma] = r->ParError(2);
		}
		      
		for (int ipar=0; ipar<nPars; ipar++) {
		  hmc_par[ipar]->SetBinContent(1,MCparVal[ipar]);
		  hmc_par[ipar]->SetBinError(1,MCparErr[ipar]);
		}

		gStyle->SetOptStat(2220);
		gStyle->SetOptFit(111);
		hmc_projy->SetMarkerColor(kRed);
		hmc_projy->SetLineColor(kRed);
		hmc_projy->SetLineWidth(3);
		if (doPrint) {
		  hmc_projy->Draw();
		  c_proj_year->Print(Form("%s/%s_projyear.pdf",outputdir.c_str(),tempH1[i][ireco]->GetName()));
		}
	      }
	    }
	    //cout << "WTF 3 " << histname.c_str() << endl;
	    // Find year range in bins
	    double binRange[2] = {0};

	    int apr2009bin = -1;

	    for (int ibin=1; ibin<=tempH1[i][ireco]->GetNbinsX(); ibin++) {
	      
	      TH1D *h_projy = ((TH2D*)tempH1[i][ireco])->ProjectionY(Form("%s_py",tempH1[i][ireco]->GetName()),ibin,ibin);
	      
	      int nentries = h_projy->GetEntries();

	      double date = h_par[igausmean]->GetBinLowEdge(ibin);
	      int year = TMath::FloorNint(date);
	      int imonth = TMath::Nint((date-year)*12);
	      TString histTitle = Form("%s %d",monthNames[imonth].Data(),year);
	      	      
	      if (binRange[0] && !binRange[1] && !nentries) binRange[1] = date;

	      if (!nentries) continue;
	      
	      if (!binRange[0]) binRange[0] = date;	      

	      if (drawNormalized) h_projy->Scale(1/h_projy->Integral());
	      
	      if (drawNormalized) {
		h_projy->GetYaxis()->SetTitle("Area Normalized");
		hmc_projy->GetYaxis()->SetTitle("Area Normalized");
	      }
	      else { 
		h_projy->GetYaxis()->SetTitle("Number of Events");
		hmc_projy->GetYaxis()->SetTitle("Number of Events");
	      }

	      h_projy->SetTitle(Form("%s (%s)", tempH1[i][ireco]->GetName(), histTitle.Data()));
	      
	      int maxBin = h_projy->GetMaximumBin();
	      float maxX = h_projy->GetBinCenter(maxBin);
	      float rms = h_projy->GetRMS();

	      float fitLow = maxX-rms;
	      float fitHigh = maxX+rms;
	      if (rms<=0) {
		fitLow = maxX - maxX/3;
		fitHigh = maxX + maxX/3;
	      }

	      double parVal[nPars] = {0};
	      double parErr[nPars] = {0};

	      parVal[imean] = h_projy->GetMean();
	      parVal[irms] = rms;
	      parErr[imean] = h_projy->GetMeanError();;
	      parErr[irms] = h_projy->GetRMSError();


	      double leftIntegralErr, rightIntegralErr, fullIntegralErr;
	      int zeroBin = h_projy->FindBin(1e-4);
	      double leftIntegral = h_projy->IntegralAndError(0, zeroBin-1, leftIntegralErr);
	      double rightIntegral = h_projy->IntegralAndError(zeroBin, h_projy->GetNbinsX()+1, rightIntegralErr);
	      double fullIntegral = h_projy->IntegralAndError(0, h_projy->GetNbinsX()+1, fullIntegralErr);
	      
	      //cout << "WTF 3." << ibin << " " << histname.c_str() << " " << nentries << " " << fitLow << " " << fitHigh << endl;

	      if (strNameOrig.find("lifetime")!=string::npos) {
		TF1 *f = new TF1(Form("f_reco%d_i%d",ireco,i),"[0]*exp(-x/[1])",1.2,10);
		f->SetParameter(1, 2);
		h_projy->Fit(f,"NRQ");
		
		parVal[imispid] = f->GetParameter(1);
		parErr[imispid] = f->GetParError(1);
		
	      } else {
		if (ise==1) {
		  parVal[imispid] = 100*(leftIntegral/fullIntegral);
		  parErr[imispid] = parVal[imispid]*sqrt(pow(leftIntegralErr/leftIntegral,2)+pow(fullIntegralErr/fullIntegral,2));
		} else {
		  parVal[imispid] = 100*(rightIntegral/fullIntegral);
		  parErr[imispid] = parVal[imispid]*sqrt(pow(rightIntegralErr/rightIntegral,2)+pow(fullIntegralErr/fullIntegral,2));
		}
	      }
	      
	      TFitResultPtr r = h_projy->Fit("gaus","NQRS","",fitLow, fitHigh);
	      //r->Print("V");
	      
	      //cout << "WTF 3." << ibin << " " << histname.c_str() << " " << r << " " << r->IsValid() << endl;
	      if (r->IsValid()) {
		parVal[igausmean] = r->Parameter(1);
		parVal[isigma] = r->Parameter(2);
		parErr[igausmean] = r->ParError(1);
		parErr[isigma] = r->ParError(2);
	      }

	      for (int ipar=0; ipar<nPars; ipar++) {
		h_par[ipar]->SetBinContent(ibin,parVal[ipar]);
		h_par[ipar]->SetBinError(ibin,parErr[ipar]);
	      }
	      
	      //gStyle->SetOptStat(2220);
	      gStyle->SetOptStat(0);
	      //gStyle->SetOptFit(111);
	      gStyle->SetOptFit(0);

	      h_projy->SetLineColor(kBlack);
	      h_projy->SetMarkerColor(kBlack);
	      h_projy->SetMarkerStyle(20);
	      h_projy->SetMarkerSize(1.2);
	      

	      if (doPrint) {
		if (doLog && hmc_projy) hmc_projy->Draw();
		h_projy->Draw();
	      }

	      // CUT LINES
	      TArrow *arw = new TArrow();
	      arw->SetLineWidth(3);
	      
	      double ymax = h_projy->GetMaximum();
	      
	      if (hmc_projy) {
		if (doPrint) 
		  if (ymax<hmc_projy->GetMaximum()) hmc_projy->Draw("HIST");
		
		if (doLog) {
		  if (!drawNormalized) {
		    h_projy->SetMinimum(0.05); 
		    hmc_projy->SetMinimum(0.05); 
		  } else {
		    h_projy->SetMinimum(1e-5); 
		    hmc_projy->SetMinimum(1e-5); 
		  }
		}
	      }

	      if (drawNormalized) {
		ymax *= 1.2;
		if (doLog) ymax*=10;
	      }
	      else {
		ymax *= 1.13;
		if (doLog) ymax*=1.5;
	      }
	      
	      
	      
	      double ybend=ymax*0.95;
	      if (doLog) ybend*=0.5;
	      double xlen = (h_projy->GetXaxis()->GetXmax()-h_projy->GetXaxis()->GetXmin())*0.06;

	      double arwlen = 0.05;

	      int ncuts=0;
	      double cut0, cut1, cut2;

	      // Tank boundaries
	      if (strNameOrig.find("_r2")!=string::npos) {
		arw->SetLineColor(kGray+2);
		arw->DrawLine(RINTK*RINTK,0.,RINTK*RINTK,ybend);
	      } 
	      else if (strNameOrig.find("_z")!=string::npos) {
		arw->SetLineColor(kGray+2);
		arw->DrawLine(ZPINTK,0.,ZPINTK,ybend);
		arw->DrawLine(-ZPINTK,0.,-ZPINTK,ybend);
	      }

	      // DEFINE CUTS
	      if (strNameOrig.find("_nmue_")!=string::npos) {
		ncuts = 2; cut0 = 0.5; cut1 = 1.5;
	      }
	      else if (strNameOrig.find("_ingate_")!=string::npos) {
		ncuts = 2; cut0 = -0.5; cut1 = 0.5;
	      }
	      else if (strNameOrig.find("_pcflg_")!=string::npos) {
		ncuts = 2; cut0 = -0.5; cut1 = 0.5;
	      }
	      else if (strNameOrig.find("_lifetime_")!=string::npos) {
		ncuts = 2; cut0 = 1.2; cut1 = 10;
	      }
	      else if (strNameOrig.find("_ise1_wall_")!=string::npos) {
		//ncuts = 3; cut0 = 100; cut1 = 200;
		ncuts = 1; cut0 = 100; 
	      }
	      else if (strNameOrig.find("_ise0_r2_")!=string::npos) {
		ncuts = 5; cut0 = 1600*1600;
	      }
	      else if (strNameOrig.find("_ise0_z_")!=string::npos) {
		ncuts = 5; cut0 = 1750;
	      }
	      //else if (strNameOrig.find("_ise0_r2zoom_")!=string::npos) {
	      //	ncuts = 2; cut0 = 1600*1600; cut1 = 1750*1750;
	      //}
	      //else if (strNameOrig.find("_ise0_zzoom_")!=string::npos) {
	      //	ncuts = 2; cut0 = 1750; cut1 = 1850; 
	      //}
	      //else if (strNameOrig.find("_n50_")!=string::npos) {
	      //ncuts = 1; cut0 = 60;
	      //}
	      else if (strNameOrig.find("_ise0_fqtotqnorm_")!=string::npos) {
		ncuts = 5; cut0 = 75;
	      }

	      arw->SetLineColor(kGreen-2);

	      if (ncuts>0) {
		h_projy->SetMaximum(ymax);

		if (ncuts==2 || ncuts==4) {
		  arw->DrawLine(cut0,0.,cut0,ybend);
		  arw->DrawArrow(cut0,ybend,cut0+xlen,ybend,arwlen,">");
		  arw->DrawLine(cut1,0.,cut1,ybend);
		  arw->DrawArrow(cut1,ybend,cut1-xlen,ybend,arwlen,">");

		  if (ncuts==4) {
		    arw->SetLineColor(kMagenta);
		    arw->DrawLine(cut2,0.,cut2,ybend);
		    arw->DrawArrow(cut2,ybend,cut2-xlen,ybend,arwlen,">");
		  }

		}
		else if (ncuts==1) {
		  arw->DrawLine(cut0,0.,cut0,ybend);
		  arw->DrawArrow(cut0,ybend,cut0+xlen,ybend,arwlen,">");
		}
		else if (ncuts==5) {
		  //arw->SetLineColor(kMagenta);
		  arw->DrawLine(cut0,0.,cut0,ybend);
		  arw->DrawArrow(cut0,ybend,cut0-xlen,ybend,arwlen,">");
		}
		else if (ncuts==3) {
		  arw->SetLineColor(kMagenta);
		  arw->DrawLine(cut0,0.,cut0,ybend);
		  arw->DrawArrow(cut0,ybend,cut0+xlen,ybend,arwlen,">");
		  arw->SetLineColor(kGreen-2);
		  arw->DrawLine(cut1,0.,cut1,ybend);
		  arw->DrawArrow(cut1,ybend,cut1+xlen,ybend,arwlen,">");
		}
	      }
	      // END CUT LINES
	      

	      // Control x-axis range
	      double norm = 1;
	      if ( (strNameOrig.find("_ise0_zzoom_")!=string::npos) ) {
		h_projy->GetXaxis()->SetRangeUser(1650, 1950);
		hmc_projy->GetXaxis()->SetRangeUser(1650, 1950);

		//norm = h_projy->Integral(h_projy->FindBin(1750), h_projy->FindBin(1850));
		//norm /= hmc_projy->Integral(hmc_projy->FindBin(1750), hmc_projy->FindBin(1850));
	      } 
	      else if ( (strNameOrig.find("_ise0_r2zoom_")!=string::npos) ) {
		h_projy->GetXaxis()->SetRangeUser(2506100, 3206100);
		hmc_projy->GetXaxis()->SetRangeUser(2506100, 3206100);

		//norm = h_projy->Integral(h_projy->FindBin(2506100), h_projy->FindBin(3206100));
		//norm /= hmc_projy->Integral(hmc_projy->FindBin(2506100), hmc_projy->FindBin(3206100));
		
	      } 
	      else if ( (strNameOrig.find("_ise0_r2zoom_")!=string::npos) ) {
		h_projy->GetXaxis()->SetRangeUser(2506100, 3206100);
		hmc_projy->GetXaxis()->SetRangeUser(2506100, 3206100);
	      } 
	      else if ( (strNameOrig.find("_ise0_wall_")!=string::npos) ) {
		h_projy->GetXaxis()->SetRangeUser(-150, 150);
		hmc_projy->GetXaxis()->SetRangeUser(-150, 150);
	      } 
	      else if ( (strNameOrig.find("_ise0_Le_Lmu_")!=string::npos) ) {
		if ( (strNameOrig.find("_momsep0_500")!=string::npos) || 
		     (strNameOrig.find("_momsep500_1000")!=string::npos) ) { 
		  h_projy->GetXaxis()->SetRangeUser(-1400, 500);
		  hmc_projy->GetXaxis()->SetRangeUser(-1400, 500);
		}
	      } 	      
	      else if ( (strNameOrig.find("_ise1_Le_Lmu_")!=string::npos) ) {
		if ( (strNameOrig.find("_momsep40_45")!=string::npos) )
		  h_projy->GetXaxis()->SetRangeUser(-150, 500);
		  hmc_projy->GetXaxis()->SetRangeUser(-150, 500);
	      } 
	      //hmc_projy->Scale(norm);


	      if (doPrint) {
		if (hmc_projy) {
		  hmc_projy->Draw("SAME HIST");
		  //h_projy->SetMinimum(TMath::Min(hmc_projy->GetMinimum(), h_projy->GetMinimum()));
		}
		h_projy->Draw("E1 SAME");
	      }

	      if (strNameOrig.find("_lifetime_")!=string::npos) {
		hmc_projy->Draw("SAME HIST");
	      }

	      if (doPrint) 
		c_proj_year->Print(Form("%s/%s_projyear.pdf",outputdir.c_str(),tempH1[i][ireco]->GetName()));

	      if (year==2009 && imonth==3) {
		apr2009bin = ibin;
		h_projy->SetTitle("");
		if (doPrint) {
		  c_proj_year->Print(Form("%s/%s_projyear_apr2009.png",outputdir.c_str(),tempH1[i][ireco]->GetName()));
		  c_proj_year->Print(Form("%s/%s_projyear_apr2009.pdf",outputdir.c_str(),tempH1[i][ireco]->GetName()));
		}
		
		if (strNameOrig.find("_lifetime_")!=string::npos) {
		  outtext << endl << h_par[imispid]->GetName() << " Lifetime" << endl;
		  outtext << "Data: " << parVal[imispid] << " +/- " << parErr[imispid]  << endl;
		  outtext << "MC: " << MCparVal[imispid] << " +/- " << MCparErr[imispid]  << endl;
		  outtext << "Diff: " << 100*(parVal[imispid]-MCparVal[imispid])/MCparVal[imispid] << " +/- " << 100*parVal[imispid]/MCparVal[imispid]*sqrt(pow(parErr[imispid]/parVal[imispid],2) + pow(MCparErr[imispid]/MCparVal[imispid],2))  << endl;
		}
		
	      }

	    }

	    //cout << "WTF 4 " << histname.c_str() << endl;
	    
	    gStyle->SetOptStat(0);
	    //gStyle->SetOptFit(111);
	    gStyle->SetOptFit(0);

	    gStyle->SetStatX(0.91);
	    gStyle->SetStatY(0.88);
	    gStyle->SetStatW(0.2);
	    gStyle->SetStatH(0.20);
	    
	    TString drawOptions = "";
	    //if (strNameOrig.find("_emom_")!=string::npos) drawOptions += "P";

	    c_proj_year->SetLogy(0);
	    if (doLog) c_proj_year->SetLogz(1);
	    
	    if (doPrint) {
	      tempH1[i][ireco]->Draw("COLZ");
	      h_par[imean]->Draw(drawOptions+"sames");
	      c_proj_year->Print(Form("%s/%s_projyear.pdf",outputdir.c_str(),tempH1[i][ireco]->GetName()));
	      //c_proj_year->Print(Form("%s/%s_projyear.png",outputdir.c_str(),tempH1[i][ireco]->GetName()));
	    }

	    double refYear = 2009+3.5/12;

	    TString fitOptions = "QR";
	    //if (strNameOrig.find("_emom_")!=string::npos) fitOptions += "w";

	    if (outfile) outfile->cd();

	 
	    for (int ipar=0; ipar<nPars; ipar++) {

	      double min;
	      double max;
	      
	      min=MCparVal[ipar];
	      max=MCparVal[ipar];
	      for (int ibin=1; ibin<=h_par[ipar]->GetNbinsX(); ibin++) {
		double val = h_par[ipar]->GetBinContent(ibin);
		if (val!=0) {
		  min = TMath::Min(val, min);
		  max = TMath::Max(val, max);
		}
	      }
	      double range = max-min;
	      h_par[ipar]->SetMinimum(min-range/2);
	      h_par[ipar]->SetMaximum(max+range/2);
	      
	      h_par[ipar]->SetMarkerColor(kBlue);
	      h_par[ipar]->SetLineColor(kBlue);
	      if (doPrint) h_par[ipar]->Draw(drawOptions);
	      
	      //TF1 *f_pol1 = new TF1(Form("f_pol1_%s",h_par[ipar]->GetName()),Form("[0]+(x-%f)*[1]",refYear),binRange[0],binRange[1]);
	      TF1 *f_pol1 = new TF1(Form("f_pol1_%s",h_par[ipar]->GetName()),Form("[0]+(x-%f)*[1]",refYear),2008,2014);
	      //TF1 *f_pol1 = new TF1(Form("f_pol1_%s",h_par[ipar]->GetName()),Form("[0]+(x-%f)*[1]+[2]*(x-%f)*x*x",refYear,refYear),2008,2014);
	      h_par[ipar]->Fit(f_pol1,fitOptions,"",binRange[0],binRange[1]);

	      outtext << endl << h_par[ipar]->GetName() << endl;
	      double f_pol1_p[2], f_pol1_p_err[2];
	      for (int ipar2=0; ipar2<2; ipar2++) {
		f_pol1_p[ipar2] = f_pol1->GetParameter(ipar2);
		f_pol1_p_err[ipar2] = f_pol1->GetParError(ipar2);
		
		outtext << "p" << ipar2 << ": " << f_pol1_p[ipar2] << " +/- " << f_pol1_p_err[ipar2]  << endl;
	      }

	      h_par[ipar]->SetTitle("");

	      if (ipar==imean) {
		
		if (strNameOrig.find("_fqtotqnorm_")!=string::npos) 
		  h_par[ipar]->GetYaxis()->SetRangeUser(0.65,0.8);

		else if (strNameOrig.find("_fqnhitsnorm_")!=string::npos) 
		  h_par[ipar]->GetYaxis()->SetRangeUser(0.42,0.46);

		else if (strNameOrig.find("_emom_")!=string::npos) 
		  h_par[ipar]->GetYaxis()->SetRangeUser(33,39);

		//else if (strNameOrig.find("_mom_over_range_")!=string::npos) 
		//h_par[ipar]->GetYaxis()->SetRangeUser(2.2,2.8);

	      }

	      if (doPrint) h_par[ipar]->Draw(drawOptions);
	      //f_pol1->Draw("SAME");

	      //const int nFitPars = 3;
 	      //double fitPar[nFitPars], fitParErr[nFitPars];
	      //for (int ifitpar=0; ifitpar<nFitPars; ifitpar++) {
	      //	fitPar[ifitpar] = f_pol1->GetParameter(ifitpar);
	      //	fitParErr[ifitpar] = f_pol1->GetParError(ifitpar);
	      //}

	      double histLowEdge = h_par[ipar]->GetBinLowEdge(1);
	      double histHighEdge = h_par[ipar]->GetBinLowEdge(h_par[ipar]->GetNbinsX()+1);

	      TLine *l_mc = new TLine(histLowEdge, MCparVal[ipar], histHighEdge, MCparVal[ipar]);
	      l_mc->SetLineWidth(4);
	      l_mc->SetLineColor(kRed);
	      l_mc->Draw("SAME");
	      if (outfile) l_mc->Write(Form("l_%s",h_par[ipar]->GetName()));
	      
	      TLine *l_apr_2009 = new TLine(refYear, h_par[ipar]->GetMinimum(), refYear, h_par[ipar]->GetMaximum());
	      l_apr_2009->SetLineWidth(7);
	      l_apr_2009->SetLineColor(kRed);
	      l_apr_2009->Draw("SAME");
	      if (outfile) l_apr_2009->Write(Form("l_apr_2009_%s",h_par[ipar]->GetName()));

	      double xBaseline = 2008.33;
	      TLine *l_baseline = new TLine(xBaseline, h_par[ipar]->GetMinimum(), xBaseline, h_par[ipar]->GetMaximum());
	      l_baseline->SetLineWidth(3);
	      l_baseline->SetLineColor(kRed);
	      //l_baseline->Draw("SAME");
	      if (outfile) l_baseline->Write(Form("l_baseline_%s",h_par[ipar]->GetName()));
	      
	      double dataval = f_pol1_p[0];
	      double dataerr = f_pol1_p_err[0];
	      //double dataval = h_par[ipar]->GetBinContent(apr2009bin);
	      //double dataerr = h_par[ipar]->GetBinError(apr2009bin);

	      double datamc = (dataval/MCparVal[ipar]);
	      double datamcerr = 100*datamc*sqrt(pow(MCparErr[ipar]/MCparVal[ipar],2) + pow(dataerr/dataval,2));

	      datamc = 100*(datamc-1);
						 

	      TLatex *t_datamc = new TLatex(0.34,0.3,Form("Data/MC: %.02f #pm %.02f%%",datamc, datamcerr));
	      t_datamc->SetTextSize(0.06);
	      t_datamc->SetTextFont(font);
	      t_datamc->SetTextColor(kRed);
	      t_datamc->SetNDC(true);
	      t_datamc->Draw("SAME");
	      	      
	      TLatex *t_slope = new TLatex(0.34,0.2,Form("%.02f #pm %.02f%% per year",100*f_pol1_p[1]/f_pol1_p[0], 100*f_pol1_p_err[1]/f_pol1_p[0]));
	      t_slope->SetTextSize(0.06);
	      t_slope->SetTextFont(font);
	      t_slope->SetTextColor(kOrange+2);
	      t_slope->SetNDC(true);
	      t_slope->Draw("SAME");
	      
	      if (doPrint) {
		h_par[ipar]->Draw(drawOptions+"sames");
		c_proj_year->Print(Form("%s/%s_projyear.pdf",outputdir.c_str(),tempH1[i][ireco]->GetName()));
		c_proj_year->Print(Form("%s/%s_projyear_%s.png",outputdir.c_str(),tempH1[i][ireco]->GetName(),parNames[ipar].Data()));
		c_proj_year->Print(Form("%s/%s_projyear_%s.pdf",outputdir.c_str(),tempH1[i][ireco]->GetName(),parNames[ipar].Data()));
	      }
	    	    
	      if (outfile) h_par[ipar]->Write();
	      if (outfile) f_pol1->Write();
	    }

	    if (doPrint) c_proj_year->Print(Form("%s/%s_projyear.pdf]",outputdir.c_str(),tempH1[i][ireco]->GetName()));

	    if (outfile) file[i][ireco]->cd();
	  }
	  
	  if (strNameOrig.find("_mom_res_vs_trumom_")!=string::npos) {
	    TLine *l_zero;
	    l_zero = new TLine(0,0,tempH1[i][ireco]->GetXaxis()->GetBinLowEdge(tempH1[i][ireco]->GetNbinsX()+1),0);
	    l_zero->SetLineWidth(2);
	    //l_zero->SetLineColor(1);

	    l_zero->Draw("SAME");
	  }


	  if (doPrint) {
	    c_reco->Print(Form("%s/%s.png",outputdir.c_str(),tempH1[i][ireco]->GetName()));
	    c_reco->Print(Form("%s/%s.pdf",outputdir.c_str(),tempH1[i][ireco]->GetName()));
	  }
	}

	// 2D ratio
	for (int i=0; i<2; i++) {
	  //((TH2D*)tempH1[i][ireco])->Rebin2D(2,2);

	  //if (ise==0) ((TH2D*)tempH1[i][ireco])->Rebin2D(4,4);

	  if (strNameOrig.find("_lnLrat1R2R_")!=string::npos) 
	    ((TH2D*)tempH1[i][ireco])->Rebin2D(4,4);
	}

	// Set binary color scale
	//set2Dcolor(0);
	//
	//TH2D *tempH1_rat = (TH2D*)(tempH1[0][ireco])->Clone();
	//tempH1_rat->Reset();
	//tempH1_rat->Divide(tempH1[1][ireco],tempH1[0][ireco]);
	//
	//TCanvas *c_reco_rat = new TCanvas(1);
	//tempH1_rat->SetMinimum(0);
	//tempH1_rat->SetMaximum(2);
	//tempH1_rat->Draw("COLZ");
	//
	//if (drawCut) f_cut->Draw("SAME");
	//  
	//if (doPrint) {
	//  c_reco_rat->Print(Form("%s/%s_rat.png",outputdir.c_str(),tempH1[0][ireco]->GetName()));
	//  //c_reco_rat->Print(Form("%s/%s_rat.pdf",outputdir.c_str(),tempH1[0][ireco]->GetName()));
	//}
	//
	//// Reset color scale
	//set2Dcolor(1);

      }
      
    }
    
  }


  if (outfile) outfile->Close();

}


double get_res1sigma(TH1 * h1tmp){
  // return histogram only integrated from 0 upto if the integral reaches to 1-sigma. 

  //TH1D * h2sigma;
  double nentry=0;
  double integral = 0.;
  //set_edge(-100.); // initialize

  //h2sigma = (TH1D*)h1tmp->Clone(); h2sigma->Reset(); // resets the bin contents and errors of an histogram       
  // Since histogram added up 1 event with weight 1, GetEntries() is fine here. otherwise need to think about overflow.
  nentry = h1tmp->Integral();
  //cout << " nentry is " << nentry <<  ", integral= " << h1tmp->Integral(1, h1tmp->GetNbinsX()) << endl;
 
  for (int i=1; i<= h1tmp->GetNbinsX(); i++){
    integral = h1tmp->Integral(1,i) / (double)nentry; // integrate 
    //cout << i <<  " integral = " <<  h1tmp->Integral(1,i ) << " fraction is " << integral << endl;
    
    if ( integral >= 0.6827)  {
      //if ( integral > get_integrallimit() )  { // get_integrallimit() is the value which we integrate the histogram up to.
      //cout << " edge ! " << endl;
      //set_edge( (i-1) * h1tmp->GetBinWidth(i) ); // return (X-axis value of the left edge of the bin)  * (bin width)
      //return h2sigma;
      double sigVal = h1tmp->GetBinLowEdge(i)+h1tmp->GetBinWidth(i);
      //cout << "sigval = " << sigVal << endl;
      return sigVal;
    }
    //} else {
    //        h2sigma->SetBinContent(i, h1tmp->GetBinContent(i));
    //        h2sigma->SetBinError(i, h1tmp->GetBinError(i));
    //}
  }
	
  //set_edge( h1tmp->GetNbinsX()*h1tmp->GetBinWidth(0) ); // assume bin width is always the same.
  //return h2sigma; // this is the case if the Ingegral() did not reach to 1-sigma.
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

