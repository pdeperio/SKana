#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>


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

using namespace std;

int main(int argc, char *argv[]){

  //gROOT->ProcessLine(".x ~/.rootlogon.C");  
  
  bool doPrint = 1;
  bool doLog = 0;
  if (argc>1)
    doPrint = atoi(argv[1]);
  if (argc>2)
    doLog = atoi(argv[2]);
    
  //cout << doPrint << " " << doLog << endl;

  gROOT->ProcessLine(Form(".! mkdir -p images")); // save eps files in figure/


  using namespace std;

  const int nRecos = 2;;
  
  TFile *file[2][nRecos];
  
  char *h_name[2] = {"h","hmc"};
  
  for (int ireco=0; ireco<nRecos; ireco++) {
    for (int i=0; i<2; i++) {
      string filename = Form("my_sk4_cosmics_reco%d_%s.root",ireco, h_name[i]);
      cout << filename.c_str() << endl;
      file[i][ireco] = new TFile(filename.c_str());
    }
  }

  
  TH1 *tempH1[2][nRecos] = {{0}};

  TText *text;
  TLine *line;

  TIter next(file[0][0]->GetListOfKeys());
  TKey *key;
  int nPlots=0;
  while ((key=(TKey*)next())) {
    
    string strNameOrig  = key->GetName();
    string strClass = key->GetClassName();

    int do1r2RmomsepFill = 2;
    
    //if (strNameOrig.find("momsep")==string::npos)
    //continue;
    
    //cout << strNameOrig.c_str() << endl;

    for (int ireco=0; ireco<nRecos; ireco++) {
      
      string strName = strNameOrig;
      
      strName.replace(strName.find("reco0"), std::string("reco0").length(), Form("reco%d",ireco));
      
      for (int i=0; i<2; i++) {
	
	strName.replace(strName.find("h_"), std::string("h_").length(), Form("%s_",h_name[i]));
	
	//cout << strName.c_str() << endl;
	tempH1[i][ireco] = (TH1*)(((TH1*)file[i][ireco]->Get(strName.c_str()))->Clone());
      }
    }

    TCanvas *c = new TCanvas(1);

    if ( !strClass.compare("TH1D") )  {
      
      //cout << strName.c_str() << endl;
      
      //TCanvas *c = new TCanvas(Form("c%d",nPlots),Form("c%d",nPlots),0,0,800,600);

      double maxY = 0;
      for (int ireco=0; ireco<nRecos; ireco++) {
		
	for (int i=0; i<2; i++) {

	  // Resolution
	  if (i && strNameOrig.find("_res_")!=string::npos)
	    tempH1[i][ireco]->Scale(1/tempH1[i][ireco]->Integral());

	  // Likelihood Ratios
	  else if (i && strNameOrig.find("_lnLrat1R")!=string::npos)
	    tempH1[i][ireco]->Scale(tempH1[0][ireco]->Integral()/tempH1[i][ireco]->Integral());


	  maxY = tempH1[i][ireco]->GetMaximum() > maxY ? tempH1[i][ireco]->GetMaximum() : maxY;
	  
	  tempH1[i][ireco]->SetLineColor(ireco+1);
	  tempH1[i][ireco]->SetMarkerColor(ireco+1);
	  
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

      TF1 *f[2][nRecos];
      for (int ireco=0; ireco<nRecos; ireco++) {
	TCanvas *c_reco = new TCanvas(1);

	if (doLog) {
	  if (strNameOrig.find("_lnLrat1R")!=string::npos) 
	    c_reco->SetLogy(1);
	}

	for (int i=0; i<2; i++) {
	  
	  if (strNameOrig.find("eff")!=string::npos) 
	    tempH1[i][ireco]->SetMinimum(0);


	  tempH1[i][ireco]->SetMaximum(TMath::Max(tempH1[0][ireco]->GetMaximum(),tempH1[1][ireco]->GetMaximum())*1.1);
	  
	  if (!i) {
	    tempH1[i][ireco]->Draw("");
	  }
	  else {
	    tempH1[i][ireco]->Draw("hist SAME");
	  }
	  
	  // Livetime
	  if (strNameOrig.find("lifetime")!=string::npos) {
	    f[i][ireco] = new TF1(Form("f_reco%d_i%d",ireco,i),"[0]*exp(-[1]*x)",1.25,20);
	    tempH1[i][ireco]->Fit(f[i][ireco],"NR");
	    
	    if (f[i][ireco]->GetParameter(1) && f[i][ireco]->GetParError(1)==f[i][ireco]->GetParError(1)) {
	      cout << tempH1[i][ireco]->GetName() << " " << 1/f[i][ireco]->GetParameter(1) << " +/- " << (1/f[i][ireco]->GetParameter(1)) * (f[i][ireco]->GetParError(1)/f[i][ireco]->GetParameter(1)) << endl;
	      f[i][ireco]->Draw("SAME");
	    }
	  }

	  // Range/Momentum
	  else if (strNameOrig.find("range_over_mom")!=string::npos) {
	    f[i][ireco] = new TF1(Form("f_reco%d_i%d",ireco,i),"gaus",0.36,0.45);
	    tempH1[i][ireco]->Fit(f[i][ireco],"NR");
	    
	    //if (f[i][ireco]->GetParameter(1) && f[i][ireco]->GetParError(1)==f[i][ireco]->GetParError(1)) {
	    cout << tempH1[i][ireco]->GetName() << " " << f[i][ireco]->GetParameter(1) << " +/- " << f[i][ireco]->GetParameter(2) << " " << f[i][ireco]->GetParError(1)/f[i][ireco]->GetParameter(1)*100 << endl;
	    f[i][ireco]->Draw("SAME");
	  }
	  

	}
	
	
	if (doPrint) {
	  string strNamePrint = strNameOrig;
	  strNamePrint.replace(strNamePrint.find("h_"), 2, "");
	  strNamePrint.replace(strNamePrint.find("reco0"), 5, Form("reco%d",ireco));

	  if (doLog)
	    strNamePrint.append("_log");
	  
	  c_reco->Print(Form("images/%s.png",strNamePrint.c_str()));
	  c_reco->Print(Form("images/%s.pdf",strNamePrint.c_str()));
	}
      }

      c->cd();
      for (int ireco=0; ireco<nRecos; ireco++) {

	for (int i=0; i<2; i++) {
 	  tempH1[i][ireco]->SetMaximum(maxY*1.1);
	  
	  if (!ireco && !i) 
	    tempH1[i][ireco]->Draw();
	  else {
	    if (i==1)
	      tempH1[i][ireco]->Draw("hist SAME");
	    else 
	      tempH1[i][ireco]->Draw("SAME");
	  }
	}
      }

      if (doPrint) {
	strNameOrig.replace(strNameOrig.find("h_"), 2, "");   
	strNameOrig.replace(strNameOrig.find("reco0"), 5,"");
	c->Print(Form("images/%sall.png",strNameOrig.c_str()));
	c->Print(Form("images/%sall.pdf",strNameOrig.c_str()));
      }
    }
 
    if ( !strClass.compare("TH2D") )  {
      
      // Cut lines
      TF1 *f_cut;
      bool drawCut = 0;
      if (strNameOrig.find("_lnLrat1R_")!=string::npos) {
	f_cut = new TF1("f_cut","x*0.2",0,2000);
	drawCut = 1;
      }
      else if (strNameOrig.find("_lnLrat1R2R_")!=string::npos) {
      	f_cut = new TF1("f_cut","150",0,1000);
	drawCut = 1;
      }

      for (int ireco=0; ireco<nRecos; ireco++) {
	for (int i=0; i<2; i++) {

	  // Likelihood Ratios
	  //if (i && strNameOrig.find("_lnLrat1R")!=string::npos)
	  //tempH1[i][ireco]->Scale(1/tempH1[i][ireco]->Integral());

	  
	  TCanvas *c_reco = new TCanvas(1);

	  tempH1[i][ireco]->Draw("COLZ");

	  if (drawCut) f_cut->Draw("SAME");

	  if (doPrint) {
	    c_reco->Print(Form("images/%s.png",tempH1[i][ireco]->GetName()));
	    c_reco->Print(Form("images/%s.pdf",tempH1[i][ireco]->GetName()));
	  }
	}
      }

    }
      
  }
}
