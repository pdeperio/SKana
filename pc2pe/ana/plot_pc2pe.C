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

  TString DrawOpts = "HIST";

  // PMT Type separated
  const int nPMTtypes = 4;
  int PMTflags[nPMTtypes] = {3, 4, 6, -1};
  TString PMTTypeNames[nPMTtypes] = {"SK2", "SK3", "HK", "Other"};

  // 1D distributions

  // All

  TCanvas *c_pc2pe = new TCanvas(1);
  c_pc2pe->SetLogy(1);
  TLegend *leg = new TLegend(0.65, 0.7, 1, 1);
  leg->SetHeader("Dataset (RMS)");

  TH1D *h_pc2pe[nFiles];

  for (int ifile=0; ifile<nFiles; ifile++) {

    TString histname = "pc2pe"+TreeVarNames[ifile];
    h_pc2pe[ifile] = new TH1D(histname, ";pc2pe Ratio;Number of Channels", 50, 0.1, 2);

    // Select good PMT flags
    TString cut_pmttype = "";
    for (int jpmttype=0; jpmttype<nPMTtypes-1; jpmttype++) {
      cut_pmttype += "pmtflag"+TreeVarNames[ifile]+Form("==%d",PMTflags[jpmttype]);
      if (jpmttype<nPMTtypes-2) cut_pmttype += " || ";
    }
  
    pc2pe->Project(histname, "rationorm"+TreeVarNames[ifile], cut_pmttype);

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


  
  TH1D *h_pc2pe_pmttypes[nFiles][nPMTtypes];
  for (int ifile=0; ifile<nFiles; ifile++) {

    TCanvas *c_pc2pe_pmttype = new TCanvas(1);
    c_pc2pe_pmttype->SetLogy(1);
    TLegend *leg_pc2pe_pmttype = new TLegend(0.6, 0.7, 1, 1);
    leg_pc2pe_pmttype->SetHeader("PMT Type (Channels, Mean, RMS)");
    
    h_pc2pe[ifile]->SetLineColor(1);
    h_pc2pe[ifile]->Draw();
    h_pc2pe[ifile]->SetTitle(FileTitles[ifile]);
    
    leg_pc2pe_pmttype->AddEntry(h_pc2pe[ifile], Form("All (%d, %.2f, %.3f)", (int)h_pc2pe[ifile]->GetEntries(), h_pc2pe[ifile]->GetMean(), h_pc2pe[ifile]->GetRMS()), "l");

    for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++) {
      
      TString histname = "pc2pe"+TreeVarNames[ifile]+"_pmt"+PMTTypeNames[ipmttype];
      h_pc2pe_pmttypes[ifile][ipmttype] = new TH1D(histname, ";pc2pe Ratio;Number of Channels", 50, 0.1, 2);

      TString cut_pmttype = "pmtflag"+TreeVarNames[ifile]+Form("==%d",PMTflags[ipmttype]);
      if (ipmttype==nPMTtypes-1) {
	cut_pmttype = "";
	for (int jpmttype=0; jpmttype<nPMTtypes-1; jpmttype++) {
	  cut_pmttype += "pmtflag"+TreeVarNames[ifile]+Form("!=%d",PMTflags[jpmttype]);
	  if (jpmttype<nPMTtypes-2) cut_pmttype += " && ";
	}
      }
      
      pc2pe->Project(histname, "rationorm"+TreeVarNames[ifile], cut_pmttype);

      int nentries = h_pc2pe_pmttypes[ifile][ipmttype]->GetEntries();
      if (!nentries) continue;
      
      h_pc2pe_pmttypes[ifile][ipmttype]->SetLineColor(ipmttype+2);
      h_pc2pe_pmttypes[ifile][ipmttype]->SetLineWidth(2);

      if (ipmttype>=nPMTtypes-1) continue;
      
      //if (!ipmttype) 
      //h_pc2pe_pmttypes[ifile][ipmttype]->Draw(); //DrawOpts += "same";
      //else
      h_pc2pe_pmttypes[ifile][ipmttype]->Draw("same");
    
      leg_pc2pe_pmttype->AddEntry(h_pc2pe_pmttypes[ifile][ipmttype], PMTTypeNames[ipmttype]+Form(" (%d, %.2f, %.3f)", nentries, h_pc2pe_pmttypes[ifile][ipmttype]->GetMean(), h_pc2pe_pmttypes[ifile][ipmttype]->GetRMS()), "l");
    }
    
    leg_pc2pe_pmttype->Draw();
    c_pc2pe_pmttype->Print("figures/pc2pe_pmttype"+TreeVarNames[ifile]+".png");
    
  }  


  
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
