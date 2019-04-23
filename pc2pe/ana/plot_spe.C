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
  int channel, group, pc2pe_bad_sk5;
  float rationorm_sk5;


  bool MakeFile = 0;
  TString pc2pe_and_spe_file = "pc2pe_output.root";
  if (!MakeFile) {

    TFile *infile = new TFile(pc2pe_and_spe_file);
    t_pc2pe = pc2pe;
    t_pc2pe->SetBranchAddress("channel", &channel);
    t_pc2pe->SetBranchAddress("group", &group);
    t_pc2pe->SetBranchAddress("rationorm_sk5", &rationorm_sk5);
    t_pc2pe->SetBranchAddress("pc2pe_bad_sk5", &pc2pe_bad_sk5);

  }
  else {
    TFile *infile = new TFile("pc2pe_output.root","update");
    t_pc2pe = (TTree*)infile->Get("pc2pe");
        
    t_pc2pe->SetBranchAddress("channel", &channel);

    for (int ispe=0; ispe<nSPEfiles; ispe++) {

      TFile *spefile = new TFile(SPEfiledir+"fit_result_"+SPEruns[ispe]+".root");
      TTree *spe = spe;

      int Channel;
      double spe_Peak;
      double spe_Peakerr;
      double spe_Mean;
      
      spe->SetBranchAddress("Channel", &Channel);
      spe->SetBranchAddress("Peak", &spe_Peak);
      spe->SetBranchAddress("Peakerr", &spe_Peakerr);

      TBranch *b_spe_Peak = t_pc2pe->Branch("spe_Peak_"+SPEruns[ispe], &spe_Peak, "spe_Peak_"+SPEruns[ispe]+"/D");
      TBranch *b_spe_Peakerr = t_pc2pe->Branch("spe_Peakerr_"+SPEruns[ispe], &spe_Peakerr, "spe_Peakerr_"+SPEruns[ispe]+"/D");

      TBranch *b_spe_Mean = t_pc2pe->Branch("spe_Mean_"+SPEruns[ispe], &spe_Mean, "spe_Mean_"+SPEruns[ispe]+"/D");

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

	spe_Mean = -1;
	TH1D *h_spe_tmp = (TH1D*)spefile->Get(Form("h_spe_onoff_%d", ichannel+1));
	if (h_spe_tmp && !h_spe_tmp->IsZombie())
	  spe_Mean = h_spe_tmp->GetMean();
	b_spe_Mean->Fill();
	  
      }
    }

    //TFile *out = new TFile(pc2pe_and_spe_file,"RECREATE");
    infile->cd();
    t_pc2pe->Write();
  }


  TH1D *h_pc2pe[nFiles];
  TH1D *h_pc2pe_count[nFiles];

  TString DrawOpts = "][ P HIST";

  TCanvas *c_pc2pe = new TCanvas("c_pc2pe", "c_pc2pe", 600, 600, 1400, 800);
  TLegend *leg = new TLegend(0.2, 0.78, 0.88, 0.88);
  leg->SetNColumns(3);

  float minY = 0.95;
  float maxY = 1.06;

  int ifile = 1;
  h_pc2pe[ifile] = (TH1D*)h_group_qisk_ton->Clone();
  TString histname = "group_pc2pe"+TreeVarNames[ifile];
  h_pc2pe[ifile]->SetName(histname);
  h_pc2pe[ifile]->SetTitle(Form(";PMT Group;Avg. Rel.Gain, SPE Peak,Mean"));

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

  const int nVars = 2;
  TString spe_var[nVars] = {"Peak", "Mean"};
  
  for (int ispe=0; ispe<nSPEfiles; ispe++) {
    for (int ivar=0; ivar<2; ivar++) {
      TH1D* h_spe = (TH1D*)h_group_qisk_ton->Clone();
      TString histname = "group_spe_"+spe_var[ivar]+"_"+SPEruns[ispe];
      h_spe->SetName(histname);
      h_spe->Reset();
  
      h_spe_count = (TH1D*)h_spe->Clone();
      h_spe_count->SetName(histname+"_count");

      TString GoodPeakCut = "spe_"+spe_var[ivar]+"_"+SPEruns[ispe]+">0";
      pc2pe->Project(histname, "group", "spe_"+spe_var[ivar]+"_"+SPEruns[ispe], GoodPeakCut);
      pc2pe->Project(histname+"_count", "group", GoodPeakCut);

      h_spe->Divide(h_spe_count);

      // Norm by average
      float avg = 0;
      for (int ibin=2; ibin<h_spe->GetNbinsX(); ibin++)
	avg += h_spe->GetBinContent(ibin);
      avg /= (h_spe->GetNbinsX() - 2);
      h_spe->Scale(1/avg);

      h_spe->SetMarkerColor(ispe+1);
      h_spe->SetMarkerStyle(ispe+2 + ivar*23);
      h_spe->SetMarkerSize(1.5);

      h_spe->Draw(DrawOpts+ "same");

      leg->AddEntry(h_spe, "SPE "+spe_var[ivar]+Form("/%.2f",avg)+" (Run "+SPEruns[ispe]+")", "p");
    }
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
  
  for (int ispe=0; ispe<nSPEfiles; ispe++) {
    for (int ivar=0; ivar<2; ivar++) {
      
      TCanvas *c_pc2pe_vs_spe = new TCanvas(1);
    
      TString histname = "pc2pe_vs_spe_"+spe_var[ivar]+"_"+SPEruns[ispe];
      TH2D *h_pc2pe_vs_spe = new TH2D(histname, ";SPE "+spe_var[ivar]+";pc2pe"+TreeVarNames[ifile]+" Ratio;Number of Channels", 50, 1, 5, 50, 0.1, 2);
    
      pc2pe->Project(histname, "rationorm"+TreeVarNames[ifile]+":spe_"+spe_var[ivar]+"_"+SPEruns[ispe]);
    
      h_pc2pe_vs_spe->Draw("colz");

      h_pc2pe_vs_spe->SetTitle("SPE Run " + SPEruns[ispe]+Form(" (Correlation = %.2f)", h_pc2pe_vs_spe->GetCorrelationFactor()));
    
      c_pc2pe_vs_spe->Print("figures/pc2pe_vs_spe_"+spe_var[ivar]+"_"+SPEruns[ispe]+".png");
    }
  }

  // Draw averaged SPE charge distributions
  TH1D *h_spe_charge_avg[3] = {0};
  int nChannels[3] = {0};
  TString RMSrangeString[3] = {"pc2pe < mean-2*RMS", "|pc2pe - mean| < 0.5*RMS", "pc2pe > mean+2*RMS"};
  
  TFile *inhistfile = new TFile("pc2pe_hists.root");
  TH1D *h_pc2pe_group_mean = (TH1D*)inhistfile->Get("group_pc2pe_sk5_2d_pfx");
  TH1D *h_pc2pe_group_rms = (TH1D*)inhistfile->Get("group_pc2pe_sk5");


  for (int ispe=0; ispe<nSPEfiles; ispe++) {
    TFile *spefile = new TFile(SPEfiledir+"fit_result_"+SPEruns[ispe]+".root");


    for (int iline=0; iline<3; iline++) {
      h_spe_charge_avg[iline] = (TH1D*)spefile->Get("h_spe_onoff_1")->Clone();
      h_spe_charge_avg[iline]->Reset();
      TString histname = h_spe_charge_avg[iline]->GetName();
      histname += "_"+SPEruns[ispe]+Form("_avg_%d", iline);
      h_spe_charge_avg[iline]->SetName(histname);
      nChannels[iline] = 0;
    }
    
    cout << spefile->GetName() << endl;
    
    for (int ichannel=0; ichannel<t_pc2pe->GetEntries(); ichannel++) {
      t_pc2pe->GetEntry(ichannel);

      if (ichannel%1000==0) cout << "Channel " << channel << endl;
      
      if (pc2pe_bad_sk5) continue;

      TH1D *h_spe_charge = (TH1D*)spefile->Get(Form("h_spe_onoff_%d", channel));
      if (!h_spe_charge || h_spe_charge->IsZombie()) continue;
      
      double pc2pe_mean = h_pc2pe_group_mean->GetBinContent(group+1);
      double pc2pe_rms = h_pc2pe_group_rms->GetBinContent(group+1);

      int iline = -1;
      if (rationorm_sk5 < pc2pe_mean - 2*pc2pe_rms) iline = 0;
      else if (fabs(rationorm_sk5-pc2pe_mean) < 0.5*pc2pe_rms) iline = 1;
      else if (rationorm_sk5 > pc2pe_mean + 2*pc2pe_rms) iline = 2;

      if (iline<0) continue;

      int nBaseBins = h_spe_charge_avg[iline]->GetNbinsX();
      int nBinsThis = h_spe_charge->GetNbinsX();
      if (nBinsThis > nBaseBins)
	h_spe_charge->Rebin(nBinsThis/nBaseBins);

      h_spe_charge->Scale(1/h_spe_charge->Integral());
      
      h_spe_charge_avg[iline]->Add(h_spe_charge);
      nChannels[iline]++;

    }

    TCanvas *c_spe_charge_avg = new TCanvas(1);
    TLegend *leg_spe_charge = new TLegend(0.4, 0.6, 1, 0.88);
    leg_spe_charge->SetHeader("Group Selection (# of Channels)");
    
    for (int iline=0; iline<3; iline++) {
      h_spe_charge_avg[iline]->SetTitle("Normalized SPE Charge Distributions (Run "+SPEruns[ispe]+")");
      h_spe_charge_avg[iline]->GetYaxis()->SetTitle(Form("Average Rate"));
      h_spe_charge_avg[iline]->Scale(1./nChannels[iline]);
      h_spe_charge_avg[iline]->SetLineColor(iline+1);
      h_spe_charge_avg[iline]->SetMarkerColor(iline+1);
      h_spe_charge_avg[iline]->SetLineWidth(2);
      if (!iline) h_spe_charge_avg[iline]->Draw();
      else h_spe_charge_avg[iline]->Draw("same");
      leg_spe_charge->AddEntry(h_spe_charge_avg[iline], RMSrangeString[iline]+Form(" (%d)", nChannels[iline]), "lp");
    }
    leg_spe_charge->Draw();
    c_spe_charge_avg->Print("figures/spe_avg_charge_"+SPEruns[ispe]+".png");
  }
}
