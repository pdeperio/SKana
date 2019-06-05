{
  gStyle->SetOptStat(0);
  //gStyle->SetPadRightMargin(0.12);

  /*
  const int nFiles = 6;

  TString FileDatasets[nFiles] = {
    "_sk4",
    "_sk5",
    "_sk5i",
    "_sk5n",
    "_sk5avg",
    "_sk4official"
  };
  enum config_enum {sk4, sk5, sk5i, sk5n, sk5avg, sk4official};

  TString TreeVarNames[nFiles] = {
    "pc2pe",
    "pc2pe",
    "pc2pe",
    "pc2pe",
    "pc2pe",
    "pc2pe"
  }
  
  int Colors[nFiles] = {kBlack, kOrange, kGreen, kMagenta, kRed, kGray+2};
  int Markers[nFiles] = {25, 23, 26, 23, 20, 24};
  int StartingMarker = 24;
  
  TString FileTitles[nFiles] = {
    "SK4 Reana",
    "SK5 Mar.",
    "SK5 Inv.",
    "SK5 Apr.",
    "SK5 Avg.",
    "SK3/4 off."
  };

  TString CanvasAppend = "_all";

  /**/

  /*
  const int nFiles = 2;

  TString FileDatasets[nFiles] = {
    "_sk4official",
    "_sk3"
  };
  enum config_enum {sk4official, sk4_qe, sk4};

  TString TreeVarNames[nFiles] = {
    "pc2pe",
    "qe"
  }
  
  int Colors[nFiles] = {kBlack, kRed,};
  int Markers[nFiles] = {24, 25};
  int StartingMarker = 24;
  
  TString FileTitles[nFiles] = {
    "SK3/4 off.",
    "SK3 QE"
  };

  TString CanvasAppend = "_sk3only";

  /**/
  /*
  
  const int nFiles = 4;
  TString TreeVarNames[nFiles] = {
    "pc2pe",
    "pc2pe",
    "pc2pe",    
    "pc2pe"
  };

  TString FileDatasets[nFiles] = {
    "_sk5",
    "_sk5i",
    "_sk5n",
    "_sk5avg"    
  };
  enum config_enum {sk5, sk5i, sk5n, sk5avg, sk4};
  
  int Colors[nFiles] = {kGreen+2, kRed, kBlue, kBlack};
  int Markers[nFiles] = {32, 26, 32, 20};
  int StartingMarker = 24;
  
  TString FileTitles[nFiles] = {
    "March",
    "Upside Down",
    "April",
    "Average"
  };

  TString CanvasAppend = "_sk5avg";

  /**/

  /*
  const int nFiles = 3;
  TString TreeVarNames[nFiles] = {
    "pc2pe",
    "pc2pe",
    "pc2pe"
  };

  TString FileDatasets[nFiles] = {
    "_sk5",
    "_sk5i",
    "_sk5n"
  };
  enum config_enum {sk5, sk5i, sk5n, sk4};
  
  int Colors[nFiles] = {kGreen+2, kRed, kBlue};
  int Markers[nFiles] = {23, 26, 23};
  int StartingMarker = 24;
  
  TString FileTitles[nFiles] = {
    "March",
    "Upside Down",
    "April"
  };

  TString CanvasAppend = "_sk5only";

  /**/
  /*
  const int nFiles = 3;

  TString TreeVarNames[nFiles] = {
    "pc2pe",
    "pc2pe",
    "pc2pe"
  };

  TString FileDatasets[nFiles] = {
    "_sk5",
    "_sk5n",
    "_sk5avg"
  };

  enum config_enum {sk5, sk5n, sk5avg, sk4};
  
  int Colors[nFiles] = {kBlue, kRed, kBlack};
  int Markers[nFiles] = {24, 25, 20};
  int StartingMarker = 24;
  
  TString FileTitles[nFiles] = {
    "March",
    "April",
    "Average"
  };

  TString CanvasAppend = "_sk5avg";
  /**/

  /*
  const int nFiles = 3;

  TString FileDatasets[nFiles] = {
    "_sk5avg",
    "_sk4official",
    "_sk4"
  };

  TString TreeVarNames[nFiles] = {
    "pc2pe",
    "pc2pe",
    "pc2pe"
  };
  
  enum config_enum {sk5avg, sk4official, sk4};
  
  int Colors[nFiles] = {kBlack, kRed, kBlue};
  int Markers[nFiles] = {20, 4, 21};
  int StartingMarker = 20;
  
  TString FileTitles[nFiles] = {
    "SK5",
    "SK3/4 official",
    "SK4"
  };

  TString CanvasAppend = "_final";

  /**/

  
  const int nFiles = 2;

  TString FileDatasets[nFiles] = {
    "_sk5avg",
    "_sk4official"
  };

  TString TreeVarNames[nFiles] = {
    "pc2pe",
    "pc2pe"
  };
  
  enum config_enum {sk5avg, sk4official, sk4};
  
  int Colors[nFiles] = {kBlack, kRed};
  int Markers[nFiles] = {20, 4};
  int StartingMarker = 20;
  
  TString FileTitles[nFiles] = {
    "SK5",
    "SK3/4 official"
  };

  TString CanvasAppend = "_final";

  /*
  const int nFiles = 2;

  TString FileDatasets[nFiles] = {
    "_sk4",
    "_sk4official"
  };

  TString TreeVarNames[nFiles] = {
    "pc2pe",
    "pc2pe"
  };
  
  enum config_enum {sk4, sk4official};
  
  int Colors[nFiles] = {kGreen-2, kRed-2};
  int Markers[nFiles] = {20, 20};
  int StartingMarker = 20;
  
  TString FileTitles[nFiles] = {
    "SK4",
    "SK3/4 official"
  };

  TString CanvasAppend = "_sk4only";

  /**/

  
  // PMT Type separated
  const int nPMTtypes = 4;
  int PMTflags[nPMTtypes] = {3, 4, 6, -1};
  TString PMTTypeNames[nPMTtypes] = {"SK2", "SK3", "HK", "Other"};

  const int nGroups = 35;
  
  TFile *infile = new TFile("pc2pe_output.root");

  TString DrawOpts = "][ P HIST";

  TCanvas *c_pc2pe = new TCanvas(1);
  TLegend *leg = new TLegend(0.2, 0.85, 0.9, 0.95);
  leg->SetNColumns(3);
  
  TCanvas *c_pc2pe_rms = new TCanvas(1);
  TLegend *leg_rms = new TLegend(0.2, 0.2, 0.8, 0.3);
  leg_rms->SetNColumns(3);

  TH1D *h_pc2pe[nFiles];
  TH1D *h_pc2pe_count[nFiles];

  float minY = 0.95;
  float maxY = 1.07;

  TCanvas *c_projections = new TCanvas("c_projections", "c_projections", 0, 0, 2000, 1600);
  c_projections->Divide(6, 6);
  TLegend *leg_projections = new TLegend(0, 0, 1, 1);
  TLegend *leg_lines = new TLegend(0, 0, 1, 1);

  TFile *outfile = new TFile("pc2pe_hists.root", "RECREATE");

  for (int ifile=0; ifile<nFiles; ifile++) {

    //if (TreeVarNames[ifile].Contains("sk5i")) continue;
      
    infile->cd();
    
    // Reject bad PMTs
    TString cut_all = "";//TreeVarNames[ifile]+"_bad"+FileDatasets[ifile]+" == 0";
    
    // 1D
    TString varname = TreeVarNames[ifile]+FileDatasets[ifile];
    TString histname = varname+"_group";
    TString histtitle = FileTitles[ifile];//+" "+TreeVarNames[ifile];
    
    h_pc2pe[ifile] = new TH1D(histname, ";PMT Group;"+TreeVarNames[ifile]+" Averaged", nGroups, 0, nGroups);

    h_pc2pe_count[ifile] = (TH1D*)h_pc2pe[ifile]->Clone();
    h_pc2pe_count[ifile]->SetName(histname+"_count");

    // Be careful here, this ruins the pc2pe average
    // This is for plot_spe.C study without HK PMTs, remove for making full plots
    cut_all += "";//"&& pmtflag"+TreeVarNames[ifile]+"!=6";
    
    pc2pe->Project(histname, "group", varname, cut_all);
    pc2pe->Project(histname+"_count", "group", cut_all);
    
    h_pc2pe[ifile]->Divide(h_pc2pe_count[ifile]);
    
    h_pc2pe[ifile]->SetMarkerColor(Colors[ifile]);
    if (nFiles>=6) h_pc2pe[ifile]->SetMarkerStyle(ifile + StartingMarker);
    else h_pc2pe[ifile]->SetMarkerStyle(Markers[ifile]);
    h_pc2pe[ifile]->SetMarkerSize(1.5);

    h_pc2pe[ifile]->GetYaxis()->SetRangeUser(minY, maxY);

    c_pc2pe->cd();
    
    if (!ifile) 
      h_pc2pe[ifile]->Draw(DrawOpts); //DrawOpts += "same";
    else 
      h_pc2pe[ifile]->Draw(DrawOpts+" same");
    
    leg->AddEntry(h_pc2pe[ifile], histtitle, "p");

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
    histname += "_2d";
    TH2F *h_pc2pe_vs_group = new TH2F(histname, FileTitles[ifile]+";PMT Group;"+TreeVarNames[ifile], nGroups, 0, nGroups, 50, 0.4, 1.6); 
    pc2pe->Project(histname, varname+":group", cut_all);

    TCanvas *c_pc2pe_vs_group = new TCanvas(1);
    
    TLegend *leg_pc2pe_vs_group = new TLegend(0.2, 0.2, 0.7, 0.3);
    leg_pc2pe_vs_group->SetNColumns(2);
    
    h_pc2pe_vs_group->Draw("COLZ");
    TH1D *h_pc2pe_vs_group_projx = h_pc2pe_vs_group->ProfileX();
    h_pc2pe_vs_group_projx->Draw("SAME");
    outfile->cd();
    h_pc2pe_vs_group_projx->Write();
    infile->cd();
    
    TH1F *h_pc2pe_copy = (TH1F*)h_pc2pe[ifile]->Clone();
    h_pc2pe_copy->Draw(DrawOpts+"same");
    h_pc2pe_copy->SetMarkerStyle(4);
    h_pc2pe_copy->SetMarkerColor(kGreen);
    leg_pc2pe_vs_group->AddEntry(h_pc2pe_copy, "Mean", "p");

    TH1F *h_pc2pe_rms = (TH1F*)h_pc2pe_copy->Clone();
    h_pc2pe_rms->Reset();

    TH1F *h_pc2pe_rms_plus_mean = (TH1F*)h_pc2pe_rms->Clone();
    
    for (int ibin=2; ibin<=h_pc2pe_vs_group->GetNbinsX()-1; ibin++) {
      TH1F *h_proj = (TH1F*)h_pc2pe_vs_group->ProjectionY(histname+Form("_px%d",ibin), ibin, ibin+1);

      double RMS = h_proj->GetRMS();
      double RMSError = h_proj->GetRMSError();
      double Mean = h_proj->GetMean();
      
      h_pc2pe_rms->SetBinContent(ibin, RMS);
      h_pc2pe_rms->SetBinError(ibin, RMSError);
      

      h_pc2pe_rms_plus_mean->SetBinContent(ibin, Mean+RMS);
      h_pc2pe_rms_plus_mean->SetBinError(ibin, RMSError); // Wrong

      // Draw all projections
      c_projections->cd(ibin+2)->SetLogy(1);
      h_proj->SetTitle(Form("Group %d; pc2pe Relative Gain; Number of Channels", ibin-1));
      //h_proj->SetLineWidth(2);
      h_proj->SetLineColor(Colors[ifile]);
      if (!ifile) h_proj->Draw();
      else h_proj->Draw("same");
      if (ibin==2) leg_projections->AddEntry(h_proj, histtitle, "l");

      if (ifile==sk4) {
	for (int iline=-2; iline<=2; iline++) {
	  TLine *l_proj = new TLine(Mean+RMS*iline, 0, Mean+RMS*iline, h_proj->GetMaximum()*5);
	  l_proj->SetLineStyle(abs(iline)+2);
	  l_proj->SetLineColor(Colors[ifile]);
	  l_proj->Draw();
	  TString line_leg = "Mean";
	  if (iline) line_leg += Form(" #pm %d*RMS", iline);
	  if (ibin==2 && iline>=0) leg_lines->AddEntry(l_proj, line_leg, "l");
	}
      }       
    }

    c_pc2pe_vs_group->cd();
    
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

    c_pc2pe_vs_group->Print("figures/"+varname+"_vs_group_rmsoverlay.png");
    
    // Compare RMSs
    c_pc2pe_rms->cd();

    h_pc2pe_rms->GetYaxis()->SetTitle("pc2pe RMS");
    h_pc2pe_rms->SetMarkerColor(Colors[ifile]);
    if (nFiles>=6) h_pc2pe_rms->SetMarkerStyle(ifile+StartingMarker);
    else h_pc2pe_rms->SetMarkerStyle(Markers[ifile]);
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

    TLine *l_unity = new TLine(0, 1, nGroups, 1);
    l_unity->SetLineColor(kGray+1);
    l_unity->Draw();

    outfile->cd();
    h_pc2pe_rms->Write();
  }

  c_pc2pe->cd();
  leg->Draw();
  c_pc2pe->Print("figures/pc2pe_grouped"+CanvasAppend+".png");
  c_pc2pe->Print("figures/pc2pe_grouped"+CanvasAppend+".pdf");

  c_pc2pe_rms->cd();
  leg_rms->Draw();
  c_pc2pe_rms->Print("figures/pc2pe_rms"+CanvasAppend+".png");
  c_pc2pe_rms->Print("figures/pc2pe_rms"+CanvasAppend+".pdf");

  c_projections->cd(2);
  leg_projections->Draw();
  c_projections->cd(3);
  leg_lines->Draw();
  c_projections->Print("figures/pc2pe_group_projections.png");
}
