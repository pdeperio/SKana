{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.15);
  
  const int nSPEfiles = 2;
  TString SPEruns[nSPEfiles] = {
    "80329",
    "80467"
  };
  
  TString SPEfiledir = "../spe_ana_junjie/";

  const int nFiles = 3;

  TString TreeVarNames[nFiles] = {
    "_sk4", "_sk5", "_sk5i"
  };
  int Colors[nFiles] = {kBlack, kBlue, kRed};

  TString FileTitles[nFiles] = {
    "SK4", "SK5", "SK5 Inv."
  };

  
  TTree *t_pc2pe;

  bool MakeFile = 0;
  TString pc2pe_and_spe_file = "pc2pe_output.root";
  if (!MakeFile) {

    TFile *infile = new TFile(pc2pe_and_spe_file);
    t_pc2pe = pc2pe;

  }
  else {
    TFile *infile = new TFile("pc2pe_output.root","update");
    t_pc2pe = (TTree*)infile->Get("pc2pe");
        
    int channel;
    t_pc2pe->SetBranchAddress("channel", &channel);

    for (int ispe=0; ispe<nSPEfiles; ispe++) {

      TFile *spefile = new TFile(SPEfiledir+"fit_result_"+SPEruns[ispe]+".root");
      TTree *spe = spe;

      int Channel;
      double spe_Peak;
      double spe_Peakerr;

      spe->SetBranchAddress("Channel", &Channel);
      spe->SetBranchAddress("Peak", &spe_Peak);
      spe->SetBranchAddress("Peakerr", &spe_Peakerr);

      TBranch *b_spe_Peak = t_pc2pe->Branch("spe_Peak_"+SPEruns[ispe], &spe_Peak, "spe_Peak_"+SPEruns[ispe]+"/D");
      TBranch *b_spe_Peakerr = t_pc2pe->Branch("spe_Peakerr_"+SPEruns[ispe], &spe_Peakerr, "spe_Peakerr_"+SPEruns[ispe]+"/D");

      for (int ichannel=0; ichannel<t_pc2pe->GetEntries(); ichannel++) {
	t_pc2pe->GetEntry(ichannel);

	int ientry;      
	for (ientry=0; ientry<t_pc2pe->GetEntries(); ientry++) {
	  spe->GetEntry(ientry);

	  if (Channel < channel) continue;
	  else if (channel == Channel) {
	    b_spe_Peak->Fill();
	    b_spe_Peakerr->Fill();
	    break;
	  } else {
	    spe_Peak = -1;      
	    spe_Peakerr = -1;
	    b_spe_Peak->Fill();
	    b_spe_Peakerr->Fill();
	    break;
	  }
	}
      }
    }

    //TFile *out = new TFile(pc2pe_and_spe_file,"RECREATE");
    infile->cd();
    t_pc2pe->Write();
  }


  TH1D *h_pc2pe[nFiles];
  TH1D *h_pc2pe_count[nFiles];


  TH1D *h_spe[nFiles];
  TH1D *h_spe_count[nFiles];
  TString DrawOpts = "][ P HIST";

  TCanvas *c_pc2pe = new TCanvas(1);
  TLegend *leg = new TLegend(0.2, 0.78, 0.88, 0.88);
  leg->SetNColumns(3);

  float minY = 0.95;
  float maxY = 1.06;

  float normpar = 2.71;

  int ifile = 1;
  h_pc2pe[ifile] = (TH1D*)h_group_qisk_ton->Clone();
  TString histname = "group_pc2pe"+TreeVarNames[ifile];
  h_pc2pe[ifile]->SetName(histname);
  h_pc2pe[ifile]->SetTitle(Form(";PMT Group;Avg. Rel.Gain or SPE Peak/%.2f", normpar));

  h_pc2pe[ifile]->Reset();
  
  h_pc2pe_count[ifile] = (TH1D*)h_pc2pe[ifile]->Clone();
  h_pc2pe_count[ifile]->SetName(histname+"_count");
    
  pc2pe->Project(histname, "group", "rationorm"+TreeVarNames[ifile]);
  pc2pe->Project(histname+"_count", "group");
    
  h_pc2pe[ifile]->Divide(h_pc2pe_count[ifile]);
    
  h_pc2pe[ifile]->SetMarkerColor(Colors[ifile]);
  h_pc2pe[ifile]->SetMarkerStyle(4);
  h_pc2pe[ifile]->SetMarkerSize(1.5);

  h_pc2pe[ifile]->GetYaxis()->SetRangeUser(minY, maxY);

  c_pc2pe->cd();
    
  h_pc2pe[ifile]->Draw(DrawOpts); //DrawOpts += "same";
    
  leg->AddEntry(h_pc2pe[ifile], FileTitles[ifile]+" pc2pe", "p");


  for (int ispe=0; ispe<nSPEfiles; ispe++) {
    h_spe[ispe] = (TH1D*)h_group_qisk_ton->Clone();
    TString histname = "group_spe_"+SPEruns[ispe];
    h_spe[ispe]->SetName(histname);
    h_spe[ispe]->Reset();
  
    h_spe_count[ispe] = (TH1D*)h_spe[ispe]->Clone();
    h_spe_count[ispe]->SetName(histname+"_count");

    TString GoodPeakCut = "spe_Peak_"+SPEruns[ispe]+">0";
    pc2pe->Project(histname, "group", "spe_Peak_"+SPEruns[ispe], GoodPeakCut);
    pc2pe->Project(histname+"_count", "group", GoodPeakCut);

    h_spe[ispe]->Divide(h_spe_count[ispe]);
    h_spe[ispe]->Scale(1/normpar); // Hardcode normalize
    
    h_spe[ispe]->SetMarkerColor(ispe+1);
    h_spe[ispe]->SetMarkerStyle(ispe+2);
    h_spe[ispe]->SetMarkerSize(1.5);

    h_spe[ispe]->Draw(DrawOpts+ "same");

    leg->AddEntry(h_spe[ispe], "SPE (Run "+SPEruns[ispe]+")", "p");

  }
  
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

  leg->Draw();

  c_pc2pe->Print("figures/pc2pe_spe_overlay.png");
  
  // 2D
  
  TH2D *h_pc2pe_vs_spe[nSPEfiles];

  for (int ispe=0; ispe<nSPEfiles; ispe++) {

    TCanvas *c_pc2pe_vs_spe = new TCanvas(1);
    
    TString histname = "pc2pe_vs_spe"+SPEruns[ispe];
    h_pc2pe_vs_spe[ispe] = new TH2D(histname, "SPE Run " + SPEruns[ispe]+";SPE Peak;pc2pe"+TreeVarNames[ifile]+" Ratio;Number of Channels", 50, 1, 5, 50, 0.1, 2);
    
    pc2pe->Project(histname, "rationorm"+TreeVarNames[ifile]+":spe_Peak_"+SPEruns[ispe]);
    
    h_pc2pe_vs_spe[ispe]->Draw("colz");
    c_pc2pe_vs_spe->Print("figures/pc2pe_vs_spe_"+SPEruns[ispe]+".png");
  }
}
