{
  gStyle->SetOptStat(0);

  TString rawfilename = "pc2pe_tst080884_and_6.root";
  //TString rawfilename = "pc2pe_tst080871_to_5.root";
  //TString rawfilename = "pc2pe_tst061889.root";
  //TString rawfilename = "";
  
  const int nFiles = 2;
  TString filenames[nFiles] = {
    //"figures_may22_fixtimestability/pc2pe_output.root",
    //"figures_may21_tightwindow/pc2pe_output.root"

    "../output_apr24/"+rawfilename,
    "../output_may22_fixtimestability/"+rawfilename
  };

  TString filetitles[nFiles] = {
    "Long",
    "Short"
  };

  const int nVars = 2;
  TString variables[nVars] = {
    //"qisk_ton",
    //"rhit_occu"
 
    "qisk_ton",
    "nhit_ton"
  };
  
  TString VarTitles[nVars] = {
    "Mean Charge",
    "Hit Rate"
  };

  TString DrawOpts = "hist p";

  TH1F *h_var[nFiles][nVars];
  TCanvas *c_var[nVars];
  for (int ivar=0; ivar<nVars; ivar++) 
    c_var[ivar] = new TCanvas(1);

  TString dataset = "";
  if (!rawfilename.CompareTo("")) dataset = "_sk4";
  
  for (int ifile=0; ifile<nFiles; ifile++) {

    TFile *infile = new TFile(filenames[ifile]);
    
    for (int ivar=0; ivar<nVars; ivar++) {
      
      h_var[ifile][ivar] = (TH1F*)infile->Get("h_"+variables[ivar]+dataset)->Clone();
      h_var[ifile][ivar]->SetLineColor(ifile+1);
      h_var[ifile][ivar]->SetMarkerColor(ifile+1);

      c_var[ivar]->cd();
      if (!ifile) {
	h_var[ifile][ivar]->Draw(DrawOpts);
      }
      else h_var[ifile][ivar]->Draw(DrawOpts+"SAME");

    }
  }

  for (int ivar=0; ivar<nVars; ivar++) 
    c_var[ivar]->Print("figures/comparewindow_"+variables[ivar]+dataset+".png");

  bool bIsDrawn = 0;
  for (int ivar=0; ivar<nVars; ivar++) {
    TCanvas *c_rat = new TCanvas(1);
    
    for (int ifile=1; ifile<nFiles; ifile++) {
      TH1F *h_rat = (TH1F*)h_var[ifile][ivar]->Clone();
      h_rat->Divide(h_var[0][ivar]);
      
      if (!bIsDrawn) {
	h_rat->Draw(DrawOpts);
	bIsDrawn = 1;
      }
      else h_rat->Draw(DrawOpts+"SAME");

      if (!rawfilename.CompareTo("")) h_rat->GetYaxis()->SetRangeUser(0.8, 1.1);
      else h_rat->GetYaxis()->SetRangeUser(0.16, 0.26);

      h_rat->GetYaxis()->SetTitle(VarTitles[ivar]);
      TF1 *f_fit = new TF1(filetitles[ifile], "pol0", 0, 11147);
      h_rat->Fit(filetitles[ifile]);
      f_fit->Draw("same");
    }
    c_rat->Print("figures/comparewindow_"+variables[ivar]+"_ratio"+dataset+".png");

  }
}
