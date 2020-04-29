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

double get_res1sigma(TH1 * h1tmp);
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

  string outputdir = Form("images_oldnew_Elt%d/norm%d_log%d/nrg%d",maxE,drawNormalized,doLog,nrgSepType);

  gROOT->ProcessLine(Form(".! mkdir -p %s",outputdir.c_str())); 

  using namespace std;
  
  const int nMCs = 2;
  const int nRecos = 1;
  const int nSE = 2;
  
  TFile *file[nMCs][nRecos] = {{0}};
  
  char *h_name[nMCs] = {"h","hmc"};
  
  ofstream outtext(Form("%s/outtext.txt",outputdir.c_str()));

  double nEntriesForLifetimeNorm[nMCs][nRecos];

  string indir[nRecos] = {
    "/Users/pdeperio/T2K/SK/SoftwareAna/Data/escale/stopmu11d_fitqun_v2r1"
    //"/Users/pdeperio/T2K/SK/SoftwareAna/Data/escale/stopmu11d_tarek_fitqun_v2r1",
    //"/Users/pdeperio/T2K/SK/SoftwareAna/Data/escale/stopmu11d_tarek_all_fitqun_v2r1"
  };

  for (int ireco=0; ireco<nRecos; ireco++) {
    for (int i=0; i<2; i++) {
      string filename = Form("%s/my_cosmics_reco1_nrg%d_Elt%d_%s.root",indir[ireco].c_str(), nrgSepType, maxE, h_name[i]);
      cout << filename.c_str() << endl;
      file[i][ireco] = new TFile(filename.c_str());
      
      //nEntriesForLifetimeNorm[i][ireco] = ((TH1D*)file[i][ireco]->Get(Form("%s_cosmic_nmue_reco%d",h_name[i],ireco)))->Integral();
      nEntriesForLifetimeNorm[i][ireco] = ((TH1D*)file[i][ireco]->Get(Form("%s_cosmic_nmue_reco1",h_name[i])))->Integral();
    }
  }

  TH1 *tempH1[nMCs][nRecos] = {{0}};

  TText *text;
  TLine *line;

  TIter next(file[0][0]->GetListOfKeys());
  TKey *key;
  int nPlots=0;

  bool canvasOpened = false;

  TCanvas *c = new TCanvas(1);

  while ((key=(TKey*)next())) {
    
    string strNameOrig  = key->GetName();
    string strClass = key->GetClassName();

    if (strNameOrig.find(doHistoName.c_str())==string::npos)
      continue;

    int ise = -1;
    if (strNameOrig.find("_ise0_")!=string::npos) ise=0;
    else if (strNameOrig.find("_ise1_")!=string::npos) ise=1;

    for (int ireco=0; ireco<nRecos; ireco++) {
      
      string strName = strNameOrig;
      
      //strName.replace(strName.find("reco0"), std::string("reco0").length(), Form("reco%d",ireco));
                  
      for (int i=0; i<2; i++) {
	
	strName.replace(strName.find("h_"), std::string("h_").length(), Form("%s_",h_name[i]));
	
	file[i][ireco]->cd();
	tempH1[i][ireco] = (TH1*)(((TH1*)file[i][ireco]->Get(strName.c_str()))->Clone());
	
	if (ireco!=1) {
	  strName.replace(strName.find("reco1"), std::string("reco1").length(), Form("reco%d",ireco));
	  tempH1[i][ireco]->SetName(strName.c_str());
	  strName.replace(strName.find(Form("reco%d",ireco)), std::string(Form("reco%d",ireco)).length(), Form("reco1"));
	}

	//if (i==1) cout << file[i][ireco]->GetName() << " " << strName.c_str() << " " << tempH1[i][ireco]->GetName() << " " << tempH1[i][ireco]->GetMean() << " " << ((TH1*)file[i][ireco]->Get(strName.c_str()))->GetMean() << endl;


	tempH1[i][ireco]->SetTitle("");

	//cout << i << " " << ireco << " " << file[i][ireco]->GetName() << " " << tempH1[i][ireco]->GetMean() << endl;
      }
    }

    if (doLog) c->SetLogy(1);

    if ( !strClass.compare("TH1D") )  {
      
      
      //TCanvas *c = new TCanvas(Form("c%d",nPlots),Form("c%d",nPlots),0,0,800,600);

      double maxY = 0;

      TText textVal[nMCs][nRecos];
      double mean[nMCs][nRecos];
      double err[nMCs][nRecos];
      
      TLine vertLine;
      vertLine.SetLineColor(kAzure-8);
      vertLine.SetLineWidth(3);

      for (int ireco=0; ireco<nRecos; ireco++) {
		
	for (int i=0; i<2; i++) {

	  textVal[i][ireco].SetNDC(true);
	  textVal[i][ireco].SetTextSize(0.05);
	  
	  if (!doLog)
	    if ( (strNameOrig.find("zenith_ise1_reco0")==string::npos) )
	      tempH1[i][ireco]->SetMinimum(0);
	  else
	    tempH1[i][ireco]->SetMinimum(1);
	  
	  // Control binning
	  if ( (strNameOrig.find("z_ise0")!=string::npos) ) {
	    tempH1[i][ireco]->Rebin(2);
	  } 
	  else if ( (strNameOrig.find("wall")!=string::npos) ) {
	    tempH1[i][ireco]->Rebin(4);
	  }
	  else if (strNameOrig.find("_momsep")!=string::npos) {
	    tempH1[i][ireco]->Rebin(4);
	  }
	  else if ( (strNameOrig.find("zenith_ise1")!=string::npos) ) {
	    tempH1[i][ireco]->Rebin(4);
	  }

		

	  // Control x-axis
	  if ( (strNameOrig.find("r2")!=string::npos) ) {
	    tempH1[i][ireco]->GetXaxis()->SetRangeUser(2200000, 3500000);
	  } 
	  if ( (strNameOrig.find("z_ise0")!=string::npos) ) {
	    tempH1[i][ireco]->GetXaxis()->SetRangeUser(1650, 2000);
	  } 
	  if ( (strNameOrig.find("cosmic_mom_ise0_reco")!=string::npos) ) {
	    tempH1[i][ireco]->GetXaxis()->SetRangeUser(100, 10000);
	  } 


	  // Resolution (1-sigma from 0)
	  if ( (i && strNameOrig.find("_res_")!=string::npos) ||
	       (strNameOrig.find("_ddir_")!=string::npos) ) {
	    
	    mean[i][ireco] = get_res1sigma(tempH1[i][ireco]);
	    err[i][ireco] = 0;
	    textVal[i][ireco].SetText(0.6, 0.8-i*0.1, Form("%.2f", mean[i][ireco]));
	    tempH1[i][ireco]->Scale(1/tempH1[i][ireco]->Integral());

	    outtext << tempH1[i][ireco]->GetName() << " res = " << mean[i][ireco] << endl;

	  }

	  // (Mean and RMS values)
	  if  ( (i && strNameOrig.find("_lvtx_res_")!=string::npos) ) {
	    
	    mean[i][ireco] = tempH1[i][ireco]->GetMean();
	    err[i][ireco] = tempH1[i][ireco]->GetRMS();
	    
	    outtext << tempH1[i][ireco]->GetName() << " mean/rms = " << mean[i][ireco] << " " << err[i][ireco] << endl;

	    textVal[i][ireco].SetText(0.6, 0.8-i*0.1, Form("%.2f %.2f", mean[i][ireco], err[i][ireco]));
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
	  
	  
	  tempH1[i][ireco]->SetLineColor(ireco+1);
	  tempH1[i][ireco]->SetMarkerColor(ireco+1);
	  if (i==0) { // && ireco==1) {
	    tempH1[i][ireco]->SetLineColor(4);
	    tempH1[i][ireco]->SetMarkerColor(4);
	  }
	  tempH1[i][ireco]->SetMarkerSize(0.8);


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
	  c_reco->SetLogx(1);
	} 
	
	TCanvas *c_mcdata[2] = {0};

	bool isDrawn = 0;
	for (int i=0; i<2; i++) {
	  	  
	  if (!tempH1[i][ireco]->GetEntries()) continue;
	 
	  tempH1[i][ireco]->SetMaximum(TMath::Max(tempH1[0][ireco]->GetMaximum()*maxScale,tempH1[1][ireco]->GetMaximum())*maxScale);

	  // Lifetime
	  if (strNameOrig.find("lifetime")!=string::npos) {

	    if (!doLog)
	      tempH1[i][ireco]->GetXaxis()->SetRangeUser(0,8);
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
	      
	      if (!isDrawn) {
		drawTH1D(tempH1[i][ireco],"");
		isDrawn = 1;
	      }
	      else {
		drawTH1D(tempH1[i][ireco],"SAME");
	      }
	      
	      //f[i][ireco]->Draw("SAME");
	      //if (doPrint) textVal[i][ireco].Draw("SAME");
	      //}
	    }
	    
	    
	    else {
	  
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
	      tempH1[i][ireco]->Draw("hist");
	      f[i][ireco]->Draw("SAME");
	      if (doPrint) textVal[i][ireco].Draw("SAME");

	    }
	    
	  } 

	  else {
	    
	    if (!isDrawn) {
	      drawTH1D(tempH1[i][ireco],"");
	      isDrawn = 1;
	    }
	    else {
	      drawTH1D(tempH1[i][ireco],"SAME");
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
	  if (strNameOrig.find("mom_over_range")!=string::npos) gausDist = 1;
	  
	  if (gausDist>=0) {
	    double gausMean = tempH1[i][ireco]->GetMean();
	    double gausRMS = tempH1[i][ireco]->GetRMS();
	    gausRange[0]=gausMean-gausRMS*(1+gausDist); 
	    gausRange[1]=gausMean+gausRMS*(1+gausDist/4);

	    //cout << "range = " << i << " " << ireco << " " << gausRange[0] << " " << gausRange[1] << endl;
	    
	    f[i][ireco] = new TF1(Form("f_reco%d_i%d",ireco,i),"gaus",gausRange[0],gausRange[1]);
	    tempH1[i][ireco]->Fit(f[i][ireco],"NRQ");
	    
	    mean[i][ireco] = f[i][ireco]->GetParameter(1);
	    if (gausDist==0) { // momentum resolution
	      err[i][ireco] = f[i][ireco]->GetParameter(2); 
	      outtext << tempH1[i][ireco]->GetName() << " mean/sig = " << mean[i][ireco] << " " <<  err[i][ireco]  << endl; 

	    }
	    else err[i][ireco] = f[i][ireco]->GetParError(1);   // momentum/range
	      
	    //if (f[i][ireco]->GetParameter(1) && f[i][ireco]->GetParError(1)==f[i][ireco]->GetParError(1)) {
	    //cout << tempH1[i][ireco]->GetName() << " " << mean[i][ireco] << " ± " << f[i][ireco]->GetParameter(2) << " " << err[i][ireco]/f[i][ireco]->GetParameter(1)*100 << endl;
	    
	    //cout << tempH1[i][ireco]->GetName() << " " << mean[i][ireco] << " ± " << err[i][ireco] << endl; 
	    
	    textVal[i][ireco].SetText(0.6, 0.8-i*0.1, Form("%.3f %.3f", mean[i][ireco], err[i][ireco]));
	    
	    if (doPrint) {
	      //f[i][ireco]->Draw("SAME");
	      if (doPrint) textVal[i][ireco].Draw("SAME");
	    }
	  }

	  
	  // Resolution
	  else if ( (i && strNameOrig.find("_res_")!=string::npos) ||
		    (strNameOrig.find("_ddir_")!=string::npos) ||
		    (strNameOrig.find("_mom_ise1_")!=string::npos) ) {
	    if (doPrint) textVal[i][ireco].Draw("SAME");
	  }
	  
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
	
	// Calculate differences
	TText textValDiff;
	textValDiff.SetNDC(true);
	textValDiff.SetTextSize(0.05);
	int i=0;
	if ( (strNameOrig.find("mom_over_range")!=string::npos) ||
	     (strNameOrig.find("_mom_ise1_")!=string::npos) ||
	     (strNameOrig.find("_ddir_")!=string::npos) ) {
	
	  double diff = 100*(mean[1][ireco]-mean[0][ireco])/mean[0][ireco];
	  double errdiff = 100*sqrt(pow(err[0][ireco],2)+pow(err[1][ireco],2))/mean[0][ireco];
	  textValDiff.SetText(0.6, 0.6, Form("%.3f %.3f", diff, errdiff));
	  if (doPrint) textValDiff.Draw("SAME");
	  
	  outtext << tempH1[i][ireco]->GetName() << " " << mean[1][ireco] << " ± " << err[1][ireco] << " " <<  mean[0][ireco] << " ± " << err[0][ireco] << " " << diff << " ± " << errdiff << endl; 
		      
	} 	
	
	
	
	if (doPrint) {
	  string strNamePrint = strNameOrig;
	  strNamePrint.replace(strNamePrint.find("h_"), 2, "");
	  strNamePrint.replace(strNamePrint.find("reco1"), 5, Form("reco%d",ireco));

	  if (doLog)
	    strNamePrint.append("_log");
	  
	  if (strNameOrig.find("lifetime")!=string::npos && doLifeTimeEff) {
	    for (int imc=0; imc<nMCs; imc++) {
	      if (c_mcdata[imc])
		c_mcdata[imc]->Print(Form("%s/%s_imc%d.png",outputdir.c_str(),strNamePrint.c_str(),imc));
	    }
	  }
	  else {
	    c_reco->Print(Form("%s/%s.png",outputdir.c_str(),strNamePrint.c_str()));
	    //c_reco->Print(Form("%s/%s.pdf",outputdir.c_str(),strNamePrint.c_str()));
	  }
	}
      }

      // Overlay all data/MC
      c->cd();
      bool isDrawn = 0;
      for (int ireco=0; ireco<nRecos; ireco++) {

	for (int i=0; i<2; i++) {
 	  tempH1[i][ireco]->SetMaximum(maxY*maxScale);

	  tempH1[i][ireco]->SetTitle(tempH1[i][ireco]->GetName());

	  if (!tempH1[i][ireco]->GetEntries()) continue;
	  
	  if (!isDrawn) {
	    drawTH1D(tempH1[i][ireco],"");
	    isDrawn = 1;
	  }
	  else {
	    drawTH1D(tempH1[i][ireco],"SAME");
	  }

	  tempH1[i][ireco]->SetTitle("");
		  
	}
	
	
  
	if (doPrint) {
	  string strNameAll = strNameOrig;
	  strNameAll.replace(strNameAll.find("h_"), 2, "");   
	  strNameAll.replace(strNameAll.find("reco1"), 5,"");
	  if (doLog) strNameAll.append("_log");
	  c->Print(Form("%s/%sall.png",outputdir.c_str(),strNameAll.c_str()));
	  //c->Print(Form("%s/%sall.pdf",outputdir.c_str(),strNameAll.c_str()));

	  if (!canvasOpened) {
	    c->Print(Form("%s/all.pdf(",outputdir.c_str()));
	    canvasOpened = true;
	  } else {
	    c->Print(Form("%s/all.pdf",outputdir.c_str()));
	  }
	 
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
	if (!isDrawn) {
	  tempH1_rat[ireco]->Draw("hist");
	  isDrawn = 1;
	} else {
	  tempH1_rat[ireco]->Draw("hist same");
	}
	
      }

      if (doPrint) {
	string strNameAll = strNameOrig;
	strNameAll.replace(strNameAll.find("h_"), 2, "");   
	strNameAll.replace(strNameAll.find("reco1"), 5,"");
	c_reco_rat->Print(Form("%s/%s_rat.png",outputdir.c_str(),strNameAll.c_str()));
	//c_reco_rat->Print(Form("%s/%s_rat.pdf",outputdir.c_str(),strNameAll.c_str()));
      }
      
    }
      
    if ( !strClass.compare("TH2D") )  {
      
      // Cut lines
      TF1 *f_cut;
      bool drawCut = 0;
      if (strNameOrig.find("_lnLrat1R_vs_mom_ise1")!=string::npos) {
	f_cut = new TF1("f_cut","x*0.2",0,tempH1[0][0]->GetBinLowEdge(tempH1[0][0]->GetNbinsX()+1));
	f_cut->SetLineColor(1);
	drawCut = 1;
      }
      //else if (strNameOrig.find("_lnLrat1R2R_")!=string::npos) {
      //	f_cut = new TF1("f_cut","150",0,1000);
      //	drawCut = 1;
      //}

      
      // Fit slices
      TH1D *h_gauspars[2][2] = {{0}};
            
      for (int ireco=0; ireco<nRecos; ireco++) {
	
	for (int i=0; i<2; i++) {
	 
	  if (drawNormalized) {
	    tempH1[i][ireco]->Scale(1/tempH1[i][ireco]->Integral());
	    
	    
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
	  
	  
	  if ( (i && strNameOrig.find("_mom_res_vs_trumom_")!=string::npos) ) {
	    
	    tempH1[i][ireco]->GetXaxis()->SetTitle("Momentum (MeV/c)");
	    tempH1[i][ireco]->GetYaxis()->SetTitle("Momentum Resolution (%)");
	    if (ise==0) ((TH2D*)tempH1[i][ireco])->Rebin2D(2,1);
	    ((TH2D*)tempH1[i][ireco])->FitSlicesY();
	    
	    for (int ipar=0; ipar<2; ipar++) {
	      h_gauspars[ireco][ipar] = (TH1D*)gDirectory->Get(Form("%s_%d",tempH1[i][ireco]->GetName(),ipar+1))->Clone();
	      h_gauspars[ireco][ipar]->SetMarkerColor(ireco+1);
	      h_gauspars[ireco][ipar]->SetLineColor(ireco+1);
	      if (ipar) h_gauspars[ireco][ipar]->SetMarkerStyle(ipar*22);
	      h_gauspars[ireco][ipar]->SetTitle("Mean (Circles), #sigma (Triangles)");
	      h_gauspars[ireco][ipar]->GetXaxis()->SetTitle("Momentum (MeV/c)");
	      h_gauspars[ireco][ipar]->GetYaxis()->SetTitle("Fitted Mean and #sigma (%)");
	      if (ise==0) h_gauspars[ireco][ipar]->GetYaxis()->SetRangeUser(-30,10);
	      else h_gauspars[ireco][ipar]->GetYaxis()->SetRangeUser(-5,25);
	    }
	    
	  }
	  /**/
	
	  
	  //else if (strNameOrig.find("cost_vs_mom")!=string::npos) {
	  //
	  //  for (int ibin=1; ibin<=tempH1[i][ireco]->GetNbinsX(); ibin++) {
	  //
	  //    double col_integral = 0;      
	  //    for (int jbin=1; jbin<=tempH1[i][ireco]->GetNbinsY(); jbin++) {
	  //	col_integral += tempH1[i][ireco]->GetBinContent(ibin,jbin);
	  //    }
	  //    
	  //    // Column normalize
	  //    for (int jbin=1; jbin<=tempH1[i][ireco]->GetNbinsY(); jbin++) {
	  //	if (col_integral)
	  //	  tempH1[i][ireco]->SetBinContent(ibin,jbin, tempH1[i][ireco]->GetBinContent(ibin,jbin)/col_integral);
	  //    }
	  //  }
	  //  
	  //  tempH1[i][ireco]->SetMaximum(0.05);
	  //}      	  
	  
	  TCanvas *c_reco = new TCanvas(1);
	  if (doLog) c_reco->SetLogz(1);

	  tempH1[i][ireco]->Draw("COLZ");

	  if (drawCut) f_cut->Draw("SAME");
	  
	  if (doPrint) {
	    c_reco->Print(Form("%s/%s.png",outputdir.c_str(),tempH1[i][ireco]->GetName()));
	    //c_reco->Print(Form("%s/%s.pdf",outputdir.c_str(),tempH1[i][ireco]->GetName()));
	  }
	}

	
	// 2D ratio
	for (int i=0; i<2; i++) {
	  ((TH2D*)tempH1[i][ireco])->Rebin2D(2,2);

	  if (ise==0) ((TH2D*)tempH1[i][ireco])->Rebin2D(4,4);

	  if (strNameOrig.find("_lnLrat1R2R_")!=string::npos) 
	    ((TH2D*)tempH1[i][ireco])->Rebin2D(4,4);
	}

	// Set binary color scale
	set2Dcolor(0);

	TH2D *tempH1_rat = (TH2D*)(tempH1[0][ireco])->Clone();
	tempH1_rat->Reset();
	tempH1_rat->Divide(tempH1[1][ireco],tempH1[0][ireco]);
	
	TCanvas *c_reco_rat = new TCanvas(1);
	tempH1_rat->SetMinimum(0);
	tempH1_rat->SetMaximum(2);
	tempH1_rat->Draw("COLZ");

	if (drawCut) f_cut->Draw("SAME");
	  
	if (doPrint) {
	  c_reco_rat->Print(Form("%s/%s_rat.png",outputdir.c_str(),tempH1[0][ireco]->GetName()));
	  //c_reco_rat->Print(Form("%s/%s_rat.pdf",outputdir.c_str(),tempH1[0][ireco]->GetName()));
	}

	// Reset color scale
	set2Dcolor(1);

      }
     

      if (h_gauspars[0][0]) {
	TCanvas *c_gauspars = new TCanvas(1);
	int isDrawn=0;
	for (int ireco=0; ireco<nRecos; ireco++) {
	  for (int ipar=0; ipar<2; ipar++) {
	   	      
	    if (!isDrawn) {
	      h_gauspars[ireco][ipar]->Draw();
	      isDrawn++;
	    } else {
	      h_gauspars[ireco][ipar]->Draw("same");
	      isDrawn++;
	    }
	  }
	}
	if (doPrint) {
	  c_gauspars->Print(Form("%s/%s_gasupars.png",outputdir.c_str(),tempH1[0][0]->GetName()));
	  //c_gauspars->Print(Form("%s/%s_gauspars.pdf",outputdir.c_str(),tempH1[0][0]->GetName()));
	}


	// Ratio
	isDrawn=0;
	int ireco=1;
	h_gauspars[ireco][0]->Add(h_gauspars[0][0],-1);  // mean
	h_gauspars[ireco][1]->Divide(h_gauspars[0][1]);  // sigma
	for (int ipar=0; ipar<2; ipar++) {

	  h_gauspars[ireco][ipar]->SetLineColor(ipar+1);
	  h_gauspars[ireco][ipar]->SetMarkerColor(ipar+1);
	  h_gauspars[ireco][ipar]->GetYaxis()->SetRangeUser(-1,2);
	  h_gauspars[ireco][ipar]->GetYaxis()->SetTitle("Mean_{New}-Mean_{Old} (%); #sigma_{New}/#sigma_{Old} Ratio");

	  if (!isDrawn) {
	    h_gauspars[ireco][ipar]->Draw();
	    isDrawn++;
	  } else {
	    h_gauspars[ireco][ipar]->Draw("same");
	    isDrawn++;
	  }
	  
	}
	if (doPrint) {
	  c_gauspars->Print(Form("%s/%s_gasupars_rat.png",outputdir.c_str(),tempH1[0][0]->GetName()));
	  //c_gauspars->Print(Form("%s/%s_gauspars_rat.pdf",outputdir.c_str(),tempH1[0][0]->GetName()));
	}
      }
    }
    
  }

  c->Print(Form("%s/all.pdf]",outputdir.c_str()));


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
