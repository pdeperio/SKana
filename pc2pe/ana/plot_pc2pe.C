{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.15);

  const int nFiles = 3;

  TString TreeVarNames[nFiles] = {
    "_sk4", "_sk5", "_sk5i"
  }
  int Colors[nFiles] = {kBlack, kRed, kBlue};

  TString FileTitles[nFiles] = {
    "SK4", "SK5", "SK5 Inv."
  };

  TFile *infile = new TFile("pc2pe_output.root");

  TCanvas *c_pc2pe = new TCanvas(1);
  c_pc2pe->SetLogy(1);
  TLegend *leg = new TLegend(0.65, 0.7, 1, 1);
  leg->SetHeader("Dataset (RMS)");

  TString DrawOpts = "HIST";

  // 1D distributions
  TH1D *h_pc2pe[nFiles];
  
  for (int ifile=0; ifile<nFiles; ifile++) {

    TString histname = "pc2pe"+TreeVarNames[ifile];
    h_pc2pe[ifile] = new TH1D(histname, ";pc2pe Ratio;Number of Channels", 50, 0.1, 2);
    
    pc2pe->Project(histname, "rationorm"+TreeVarNames[ifile]);

    h_pc2pe[ifile]->SetLineColor(Colors[ifile]);
    h_pc2pe[ifile]->SetLineWidth(2);

    if (!ifile) 
      h_pc2pe[ifile]->Draw(); //DrawOpts += "same";
    else 
      h_pc2pe[ifile]->Draw("same");
    
    leg->AddEntry(h_pc2pe[ifile], FileTitles[ifile]+Form(" (%.3f)", h_pc2pe[ifile]->GetRMS()), "l");
  }  

  leg->Draw();
  c_pc2pe->Print("figures/pc2pe.png");

  
  // 2D Distributions
  TH2D *h_pc2pe_2d[nFiles];
  
  for (int ifile=0; ifile<nFiles-1; ifile++) {

    TCanvas *c_pc2pe_2d = new TCanvas(1);

    TString histname = "pc2pe_2d"+TreeVarNames[ifile];
    h_pc2pe_2d[ifile] = new TH2D(histname, ";pc2pe"+TreeVarNames[ifile]+" Ratio;pc2pe"+TreeVarNames[ifile+1]+" Ratio;Number of Channels", 50, 0.1, 2, 50, 0.1, 2);
    
    pc2pe->Project(histname, "rationorm"+TreeVarNames[ifile]+":rationorm"+TreeVarNames[ifile+1]);

    h_pc2pe_2d[ifile]->Draw("colz");

    c_pc2pe_2d->Print("figures/pc2pe"+TreeVarNames[ifile]+"_vs"+TreeVarNames[ifile+1]+".png");
  }
}
