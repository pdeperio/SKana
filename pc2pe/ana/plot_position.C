{
  gStyle->SetOptStat(0);
  //gStyle->SetPadRightMargin(0.12);

  const int nFiles = 3;

  TString TreeVarNames[nFiles] = {
    "_sk4", "_sk5", "_sk5i"
  }
  int Colors[nFiles] = {kBlack, kRed, kBlue};

  TString FileTitles[nFiles] = {
    "SK4", "SK5", "SK5 Inv."
  };

  // PMT Type separated
  const int nPMTtypes = 4;
  int PMTflags[nPMTtypes] = {3, 4, 6, -1};
  TString PMTTypeNames[nPMTtypes] = {"SK2", "SK3", "HK", "Other"};


  TFile *infile = new TFile("pc2pe_output.root");

  TString DrawOpts = "][ P HIST";

  TCanvas *c_pc2pe = new TCanvas(1);
  TLegend *leg = new TLegend(0.2, 0.78, 0.88, 0.88);
  leg->SetNColumns(3);
  
  TCanvas *c_pc2pe_rms = new TCanvas(1);
  TLegend *leg_rms = new TLegend(0.2, 0.2, 0.8, 0.3);
  leg_rms->SetNColumns(3);

  TH1D *h_pc2pe[nFiles];
  TH1D *h_pc2pe_count[nFiles];

  float minY = 0.95;
  float maxY = 1.07;
  
  for (int ifile=0; ifile<nFiles; ifile++) {

    // Select good PMT flags
    //TString cut_pmttype = "(";
    //for (int jpmttype=0; jpmttype<nPMTtypes-1; jpmttype++) {
    //  cut_pmttype += "pmtflag"+TreeVarNames[ifile]+Form("==%d",PMTflags[jpmttype]);
    //  if (jpmttype<nPMTtypes-2) cut_pmttype += " || ";
    //}
    //cut_pmttype += ")";
    //cut_pmttype = "1";  // Ignore cut
    
    // 1D
    h_pc2pe[ifile] = (TH1D*)h_group_qisk_ton->Clone();
    TString histname = "group_pc2pe"+TreeVarNames[ifile];
    h_pc2pe[ifile]->SetName(histname);
    h_pc2pe[ifile]->SetTitle("");
    h_pc2pe[ifile]->Reset();

    h_pc2pe_count[ifile] = (TH1D*)h_pc2pe[ifile]->Clone();
    h_pc2pe_count[ifile]->SetName(histname+"_count");
    
    pc2pe->Project(histname, "group", "rationorm"+TreeVarNames[ifile]);
    pc2pe->Project(histname+"_count", "group");
    
    h_pc2pe[ifile]->Divide(h_pc2pe_count[ifile]);
    
    h_pc2pe[ifile]->SetMarkerColor(Colors[ifile]);
    h_pc2pe[ifile]->SetMarkerStyle(ifile+2);
    h_pc2pe[ifile]->SetMarkerSize(1.5);

    h_pc2pe[ifile]->GetYaxis()->SetRangeUser(minY, maxY);

    c_pc2pe->cd();
    
    if (!ifile) 
      h_pc2pe[ifile]->Draw(DrawOpts); //DrawOpts += "same";
    else 
      h_pc2pe[ifile]->Draw(DrawOpts+" same");
    
    leg->AddEntry(h_pc2pe[ifile], FileTitles[ifile], "p");

    // Geometry separation
    for (int isep=18; isep<=26; isep+=8) {
      TLine *l_Sep = new TLine(isep, minY, isep, maxY);
      l_Sep->SetLineColor(kGreen-2);
      l_Sep->SetLineWidth(2);
      l_Sep->Draw();
    }

    TLine *l_unity = new TLine(0, 1, 35, 1);
    l_unity->SetLineColor(kGray+1);
    l_unity->Draw();


    // 2D
    histname = "2d_"+histname;
    TH2F *h_pc2pe_vs_group = new TH2F(histname, FileTitles[ifile]+";PMT Group; Relative Gain (pc2pe ratio)", 35, 0, 35, 50, 0.4, 1.6); 
    pc2pe->Project(histname, "rationorm"+TreeVarNames[ifile]+":group");

    TCanvas *c_pc2pe_vs_group = new TCanvas(1);
    
    TLegend *leg_pc2pe_vs_group = new TLegend(0.2, 0.2, 0.7, 0.3);
    leg_pc2pe_vs_group->SetNColumns(2);
    
    h_pc2pe_vs_group->Draw("COLZ");
    h_pc2pe_vs_group->ProfileX()->Draw("SAME");
    
    TH1F *h_pc2pe_copy = (TH1F*)h_pc2pe[ifile]->Clone();
    h_pc2pe_copy->Draw(DrawOpts+"same");
    h_pc2pe_copy->SetMarkerStyle(4);
    h_pc2pe_copy->SetMarkerColor(kGreen);
    leg_pc2pe_vs_group->AddEntry(h_pc2pe_copy, "Mean", "p");

    TH1F *h_pc2pe_rms = (TH1F*)h_pc2pe_copy->Clone();
    h_pc2pe_rms->Reset();

    TH1F *h_pc2pe_rms_plus_mean = (TH1F*)h_pc2pe_rms->Clone();
    
    for (int ibin=2; ibin<=h_pc2pe_vs_group->GetNbinsX()-1; ibin++) {
      TH1F *h_proj = (TH1F*)h_pc2pe_vs_group->ProjectionY(Form("_px%d",ibin), ibin, ibin+1);
      h_pc2pe_rms->SetBinContent(ibin, h_proj->GetRMS());
      h_pc2pe_rms->SetBinError(ibin, h_proj->GetRMSError());

      h_pc2pe_rms_plus_mean->SetBinContent(ibin, h_proj->GetMean()+h_proj->GetRMS());
      h_pc2pe_rms_plus_mean->SetBinError(ibin, h_proj->GetRMSError());
    }
    h_pc2pe_rms_plus_mean->SetMarkerStyle(2);
    h_pc2pe_rms_plus_mean->SetMarkerColor(kRed);
    h_pc2pe_rms_plus_mean->Draw("same");
    leg_pc2pe_vs_group->AddEntry(h_pc2pe_rms_plus_mean, "Mean+RMS", "p");
    
    for (int isep=18; isep<=26; isep+=8) {
      TLine *l_Sep = new TLine(isep, 0.4, isep, 1.6);
      l_Sep->SetLineColor(kGreen-2);
      l_Sep->SetLineWidth(2);
      l_Sep->Draw();
    }

    TLine *l_unity = new TLine(0, 1, 35, 1);
    l_unity->SetLineColor(kGray+1);
    l_unity->Draw();

    leg_pc2pe_vs_group->Draw();

    c_pc2pe_vs_group->Print("figures/pc2pe_vs_group_rmsoverlay"+TreeVarNames[ifile]+".png");
    
    // Compare RMSs
    c_pc2pe_rms->cd();

    h_pc2pe_rms->GetYaxis()->SetTitle("pc2pe RMS");
    h_pc2pe_rms->SetMarkerColor(Colors[ifile]);
    h_pc2pe_rms->SetMarkerStyle(ifile+2);
    h_pc2pe_rms->SetMarkerSize(1.5);

    h_pc2pe_rms->GetYaxis()->SetRangeUser(0, 0.09);
      
    if (!ifile) 
      h_pc2pe_rms->Draw(DrawOpts); //DrawOpts += "same";
    else 
      h_pc2pe_rms->Draw(DrawOpts+" same");
    
    leg_rms->AddEntry(h_pc2pe_rms, FileTitles[ifile], "p");

    // Geometry separation
    for (int isep=18; isep<=26; isep+=8) {
      TLine *l_Sep = new TLine(isep, 0, isep, 0.09);
      l_Sep->SetLineColor(kGreen-2);
      l_Sep->SetLineWidth(2);
      l_Sep->Draw();
    }

    TLine *l_unity = new TLine(0, 1, 35, 1);
    l_unity->SetLineColor(kGray+1);
    l_unity->Draw();

  }  

  c_pc2pe->cd();
  leg->Draw();
  c_pc2pe->Print("figures/pc2pe_grouped.png");

  c_pc2pe_rms->cd();
  leg_rms->Draw();
  c_pc2pe_rms->Print("figures/pc2pe_rms.png");
}
