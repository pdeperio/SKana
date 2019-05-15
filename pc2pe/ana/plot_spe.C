void plot_spe(bool MakeFile = 0)
{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.15);
  gErrorIgnoreLevel = kWarning;  // Suppress "Info in <TCanvas::Print>" messages

  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);

  const int nSPEfiles = 1;
  TString SPEruns[nSPEfiles] = {
    //"80329",
    "80467"
  };
  
  TString SPEfiledir = "../spe_ana_junjie/";

  const int nFiles = 5;

  TString TreeVarNames[nFiles] = {
    "_sk4", "_sk5", "_sk5i", "_sk5n", "_sk5avg"
  };
  int Colors[nFiles] = {
    kGray+2,
    kRed,
    kBlue,
    kGreen-2,
    kBlack
  };

  TString FileTitles[nFiles] = {
    "SK4",
    "SK5",
    "SK5 Inv.",
    "SK5 New",
    "SK5 Avg."
  };
  enum config_enum {sk4, sk5, sk5i, sk5n, sk5avg};

  const int nGroups = 35;

  // PMT Type separated
  const int nPMTtypes = 4;
  int PMTflags[nPMTtypes] = {3, 4, 6, -1};
  TString PMTTypeNames[nPMTtypes] = {"SK2", "SK3", "HK", "Other"};
  int PMTtypeColors[nPMTtypes] = {kBlue, kRed, kBlack, kGray};
  
  TTree *t_pc2pe;
  int channel, group, pc2pe_bad, pmtflag;
  float rationorm, rationorm_err;

  TString pc2pe_and_spe_file = "pc2pe_output.root";
  if (!MakeFile) {

    TFile *infile = new TFile(pc2pe_and_spe_file);
    t_pc2pe = pc2pe;
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
      double spe_Meanerr;
      
      spe->SetBranchAddress("Channel", &Channel);
      spe->SetBranchAddress("Peak", &spe_Peak);
      spe->SetBranchAddress("Peakerr", &spe_Peakerr);

      TBranch *b_spe_Peak = t_pc2pe->Branch("spe_Peak_"+SPEruns[ispe], &spe_Peak, "spe_Peak_"+SPEruns[ispe]+"/D");
      TBranch *b_spe_Peakerr = t_pc2pe->Branch("spe_Peakerr_"+SPEruns[ispe], &spe_Peakerr, "spe_Peakerr_"+SPEruns[ispe]+"/D");

      TBranch *b_spe_Mean = t_pc2pe->Branch("spe_Mean_"+SPEruns[ispe], &spe_Mean, "spe_Mean_"+SPEruns[ispe]+"/D");
      TBranch *b_spe_Meanerr = t_pc2pe->Branch("spe_Meanerr_"+SPEruns[ispe], &spe_Meanerr, "spe_Meanerr_"+SPEruns[ispe]+"/D");

      int last_entry=0;

      for (int ichannel=0; ichannel<t_pc2pe->GetEntries(); ichannel++) {
	t_pc2pe->GetEntry(ichannel);

	if (ichannel%1000==0) cout << "ichannel " << ichannel << endl;

	int ientry;
	for (ientry=last_entry; ientry<t_pc2pe->GetEntries(); ientry++) {
	  spe->GetEntry(ientry);

	  if (Channel < channel) continue;
	  else if (channel == Channel) {
	    b_spe_Peak->Fill();
	    b_spe_Peakerr->Fill();
	    last_entry = ientry;
	    break;
	  } else {
	    spe_Peak = -1;      
	    spe_Peakerr = -1;
	    b_spe_Peak->Fill();
	    b_spe_Peakerr->Fill();
	    last_entry = ientry;
	    break;
	  }
	}

	spe_Mean = -1;
	spe_Meanerr = -1;
	TH1D *h_spe_tmp = (TH1D*)spefile->Get(Form("h_spe_onoff_%d", ichannel+1));
	if (h_spe_tmp && !h_spe_tmp->IsZombie()) {
	  spe_Mean = h_spe_tmp->GetMean();
	  spe_Meanerr = h_spe_tmp->GetMeanError();
	}
	b_spe_Mean->Fill();
	b_spe_Meanerr->Fill();
      }
    }

    //TFile *out = new TFile(pc2pe_and_spe_file,"RECREATE");
    infile->cd();
    t_pc2pe->Write();
  }

  int ifile = sk5avg;

  t_pc2pe->SetBranchAddress("channel", &channel);
  t_pc2pe->SetBranchAddress("group", &group);
  t_pc2pe->SetBranchAddress("rationorm"+TreeVarNames[ifile], &rationorm);
  t_pc2pe->SetBranchAddress("rationorm_err"+TreeVarNames[ifile], &rationorm_err);
  t_pc2pe->SetBranchAddress("pc2pe_bad"+TreeVarNames[ifile], &pc2pe_bad);
  t_pc2pe->SetBranchAddress("pmtflag"+TreeVarNames[ifile], &pmtflag);

  double spe_means[nSPEfiles];
  double spe_meanerrs[nSPEfiles];
  for (int ispe=0; ispe<nSPEfiles; ispe++) {
    t_pc2pe->SetBranchAddress("spe_Mean_"+SPEruns[ispe], &spe_means[ispe]);
    t_pc2pe->SetBranchAddress("spe_Meanerr_"+SPEruns[ispe], &spe_meanerrs[ispe]);
  }

  TH1D *h_pc2pe[nFiles];
  TH1D *h_pc2pe_count[nFiles];

  TString DrawOpts = "][ P HIST";

  TCanvas *c_pc2pe = new TCanvas("c_pc2pe", "c_pc2pe", 600, 600, 1400, 800);
  TLegend *leg = new TLegend(0.2, 0.78, 0.88, 0.88);
  leg->SetNColumns(3);

  float minY = 0.95;
  float maxY = 1.06;

  TString histname = "group_pc2pe"+TreeVarNames[ifile];
  h_pc2pe[ifile] = new TH1D(histname, ";PMT Group;Avg. Rel.Gain, SPE Peak,Mean", nGroups, 0, nGroups);
  
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
      TString histname = "group_spe_"+spe_var[ivar]+"_"+SPEruns[ispe];
      TH1D* h_spe =  new TH1D(histname, "Run "+SPEruns[ispe]+";PMT Group;"+spe_var[ivar], nGroups, 0, nGroups);
  
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

      // Plot only one SPE run)
      if (!ispe) {
	h_spe->Draw(DrawOpts+ "same");
	leg->AddEntry(h_spe, "SPE "+spe_var[ivar]+Form("/%.2f",avg)+" (Run "+SPEruns[ispe]+")", "p");
      }
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
  TF1 *f_prof;
  
  for (int ispe=0; ispe<nSPEfiles; ispe++) {
    for (int ivar=0; ivar<2; ivar++) {
      
      TCanvas *c_pc2pe_vs_spe = new TCanvas(1);
    
      TString histname = "pc2pe_vs_spe_"+spe_var[ivar]+"_"+SPEruns[ispe];
      TH2D *h_pc2pe_vs_spe = new TH2D(histname, ";SPE "+spe_var[ivar]+";"+FileTitles[ifile]+" pc2pe;Number of Channels", 120, 2, 4.4, 120, 0.6, 1.7);

      TString cut_all = "pc2pe_bad"+TreeVarNames[ifile]+"==0";  // No pc2pe bad channels
      cut_all += " && badchannel==0";  // No official bad channels
      cut_all += " && pmtflag"+TreeVarNames[ifile]+"!=6";  // No HK PMTs

      TString xvar = "spe_"+spe_var[ivar]+"_"+SPEruns[ispe];
      TString yvar = "rationorm"+TreeVarNames[ifile];
      
      pc2pe->Project(histname, yvar+":"+xvar, cut_all);
    
      h_pc2pe_vs_spe->Draw("colz");

      if (ivar) {
	TProfile *prof = h_pc2pe_vs_spe->ProfileX();
	prof->SetMarkerColor(kGray);
	prof->SetLineColor(kGray);
	prof->Draw("SAME");
	f_prof = new TF1(Form("f_%s", h_pc2pe_vs_spe->GetName()), "pol1", 2, 4.4);
	prof->Fit(f_prof, "N");
	f_prof->SetLineColor(kGray+2);
	f_prof->SetLineWidth(1);
	f_prof->Draw("SAME");

	TLatex *fLineText = new TLatex(3.25, 1.53, Form("pc2pe = %.2f #times SPEmean + %.2f", f_prof->GetParameter(1), f_prof->GetParameter(0)));
	fLineText->SetTextSize(0.04);
	fLineText->Draw();
      }
      
      h_pc2pe_vs_spe->SetTitle("SPE Run " + SPEruns[ispe]+Form(" (Correlation = %.2f)", h_pc2pe_vs_spe->GetCorrelationFactor()));
            
      //TPaveText *tRMSoverMeans = new TPaveText(0.2, 0.7, 0.5, 0.88);
      TPaveText *tRMSoverMeans;
      if (ivar) tRMSoverMeans = new TPaveText(2.1, 1.4, 3, 1.6);
      else tRMSoverMeans = new TPaveText(3.4, 1.4, 4.3, 1.6);
      tRMSoverMeans->SetFillColorAlpha(0, 0);
      tRMSoverMeans->SetBorderSize(1);
      tRMSoverMeans->AddText("RMS/Mean (%)");
      
      double RMS2d, Mean2d, RMSErr2d, MeanErr2d;
      for (int iaxis=0; iaxis<=1; iaxis++) {
	Mean2d = h_pc2pe_vs_spe->GetMean(iaxis+1);
	MeanErr2d = h_pc2pe_vs_spe->GetMeanError(iaxis+1);
	RMS2d = h_pc2pe_vs_spe->GetRMS(iaxis+1);
	RMSErr2d = h_pc2pe_vs_spe->GetRMSError(iaxis+1);
	
	double RMSoverMean = RMS2d/Mean2d*100;
	double RMSoverMeanErr = RMSoverMean*sqrt(pow(MeanErr2d/Mean2d, 2)+pow(RMSErr2d/RMS2d, 2));

	TString AxisTitle = "pc2pe: ";
	if (!iaxis) AxisTitle = "SPE "+spe_var[ivar]+": ";
	tRMSoverMeans->AddText(AxisTitle+Form("%.2f #pm %.2f", RMSoverMean, RMSoverMeanErr));
      }
      tRMSoverMeans->Draw();

      c_pc2pe_vs_spe->Print("figures/pc2pe_vs_spe_"+spe_var[ivar]+"_"+SPEruns[ispe]+".png");

    }
  }

  // PMT Type separated
  TF1 *f_anomalous;
  
  for (int ispe=0; ispe<nSPEfiles; ispe++) {

    cout << "Run " << SPEruns[ispe] << endl;
    
    for (int ivar=0; ivar<2; ivar++) {
      
      cout << "SPE " << spe_var[ivar] << endl;

      TString axis_vars[2];
      axis_vars[0] = "spe_"+spe_var[ivar]+"_"+SPEruns[ispe];
      axis_vars[1] = "rationorm"+TreeVarNames[ifile];
      
      cout << "    " << axis_vars[0].Data() << " " << axis_vars[1].Data() <<endl;
      
      TCanvas *c_pc2pe_vs_spe = new TCanvas(1);
      TLegend *leg_pc2pe_pmttype = new TLegend(0.2, 0.65, 0.5, 0.88);
      leg_pc2pe_pmttype->SetHeader("PMT Type (Corr., x #sigma/#mu (%), y #sigma/#mu (%))");

      for (int ipmttype=0; ipmttype<nPMTtypes-1; ipmttype++) {

	TString histname = "pc2pe_vs_spe_"+spe_var[ivar]+"_"+SPEruns[ispe]+"_pmt"+PMTTypeNames[ipmttype];;
	TH2D *h_pc2pe_vs_spe_pmttype = new TH2D(histname, ";SPE "+spe_var[ivar]+";"+FileTitles[ifile]+" pc2pe;Number of Channels", 120, 2, 4.4, 120, 0.6, 1.7);

	TString cutstring = "pc2pe_bad"+TreeVarNames[ifile]+"==0 && badchannel==0 && pmtflag"+TreeVarNames[ifile]+Form("==%d", PMTflags[ipmttype])+" && rationorm"+TreeVarNames[ifile]+" < 4";

	pc2pe->Project(histname, axis_vars[1]+":"+axis_vars[0], cutstring);
	pc2pe->Draw(axis_vars[1]+":"+axis_vars[0], cutstring, "goff");

	TGraph *gr = new TGraph(pc2pe->GetSelectedRows(), pc2pe->GetV2(), pc2pe->GetV1());
	
	if (!gr->GetN()) continue;
	
	gr->SetTitle(";SPE "+spe_var[ivar]+";"+FileTitles[ifile]+" pc2pe");
	gr->SetMarkerSize(0.06);
	//gr->SetMarkerColor(PMTtypeColors[ipmttype]);
	gr->SetMarkerColorAlpha(PMTtypeColors[ipmttype], 1);
	gr->SetFillColor(PMTtypeColors[ipmttype]);
	
	if (!ipmttype) gr->Draw("AP");
	else gr->Draw("same P");

	if (ivar) {
	  f_anomalous = new TF1("f_anomalous", "[0]+[1]*(x-[2])", 2, 4.4);
	  f_anomalous->SetParameters(f_prof->GetParameter(0), f_prof->GetParameter(1), 0.3);
	  f_anomalous->SetLineColor(kGray);
	  f_anomalous->Draw("SAME");
	}

	TString RMSoverMeans;

	cout << PMTTypeNames[ipmttype] << ": ";

	for (int iaxis=0; iaxis<=1; iaxis++) {
	  double RMS2d, Mean2d, RMSErr2d, MeanErr2d;
	  Mean2d = h_pc2pe_vs_spe_pmttype->GetMean(iaxis+1);
	  MeanErr2d = h_pc2pe_vs_spe_pmttype->GetMeanError(iaxis+1);
	  RMS2d = h_pc2pe_vs_spe_pmttype->GetRMS(iaxis+1);
	  RMSErr2d = h_pc2pe_vs_spe_pmttype->GetRMSError(iaxis+1);
	  
	  double RMSoverMean = RMS2d/Mean2d*100;
	  double RMSoverMeanErr = RMSoverMean*sqrt(pow(MeanErr2d/Mean2d, 2)+pow(RMSErr2d/RMS2d, 2));

	  cout << Form("%.2f ± %.2f", RMSoverMean, RMSoverMeanErr);
	  RMSoverMeans += Form("%.2f", RMSoverMean);
	  if (!iaxis) {
	    cout  << ", ";
	    RMSoverMeans += ", ";
	  }
	}
	cout << endl;
	
	leg_pc2pe_pmttype->AddEntry(gr, PMTTypeNames[ipmttype]+Form(" (%.2f, %s)",
								    gr->GetCorrelationFactor(),
								    RMSoverMeans.Data()),
				    "f");

	//if (ipmttype==2) gr->Print();

      }
      cout << endl;
      
      leg_pc2pe_pmttype->Draw();
      c_pc2pe_vs_spe->Print("figures/pc2pe_vs_spe_"+spe_var[ivar]+"_"+SPEruns[ispe]+"_pmttypes.png");
    }
    cout << endl;

  }

  break;
    
  // Draw averaged SPE charge distributions
  const int nRMSranges = 3;
  TH1D *h_spe_charge_avg[nRMSranges] = {0};
  int nChannels[nRMSranges] = {0};
  TString RMSrangeString[nRMSranges] = {"pc2pe < mean-2*RMS", "|pc2pe - mean| < 0.5*RMS", "pc2pe > mean+2*RMS"};
  TString RMSrangeName[nRMSranges] = {"minus2rms", "mean", "plus2rms"};

  ofstream fPMTlists[nRMSranges];
  TCanvas *c_spe_all[nRMSranges] = {0};
  TString c_spe_all_names[nRMSranges];

  TH1D *h_spe_charge_avg_pmttype[nPMTtypes] = {0};
  int nChannels_pmttype[nPMTtypes] = {0};
  for (int iline=0; iline<nRMSranges; iline++) {
    fPMTlists[iline].open("channel_list_pc2pe_"+RMSrangeName[iline]+".txt");
    if (iline!=1) {
      c_spe_all[iline] = new TCanvas("c_"+RMSrangeName[iline], "c_"+RMSrangeName[iline], 800, 600);
      c_spe_all_names[iline] = "figures/spe_dists_"+RMSrangeName[iline]+".pdf";
      c_spe_all[iline]->Print(c_spe_all_names[iline]+"[");
    }
  }

  TString c_spe_anomalous_name = "spe_anomalous";
  TCanvas *c_spe_anomalous = new TCanvas(c_spe_anomalous_name, c_spe_anomalous_name, 800, 600);
  c_spe_anomalous_name = "figures/"+c_spe_anomalous_name+".pdf";
  c_spe_anomalous->Print(c_spe_anomalous_name+"[");


  TFile *inhistfile = new TFile("pc2pe_hists.root");
  TH1D *h_pc2pe_group_mean = (TH1D*)inhistfile->Get("group_pc2pe"+TreeVarNames[ifile]+"_2d_pfx");
  TH1D *h_pc2pe_group_rms = (TH1D*)inhistfile->Get("group_pc2pe"+TreeVarNames[ifile]);

  for (int ispe=0; ispe<nSPEfiles; ispe++) {
    
    TFile *spefile = new TFile(SPEfiledir+"fit_result_"+SPEruns[ispe]+".root");

    // Initialize SPE average histo for extreme pc2pe's
    for (int iline=0; iline<nRMSranges; iline++) {
      h_spe_charge_avg[iline] = (TH1D*)spefile->Get("h_spe_onoff_1")->Clone();
      h_spe_charge_avg[iline]->Reset();
      TString histname = h_spe_charge_avg[iline]->GetName();
      histname += "_"+SPEruns[ispe]+Form("_avg_%d", iline);
      h_spe_charge_avg[iline]->SetName(histname);
      nChannels[iline] = 0;
    }

    // Initialize SPE average histo for different PMT types
    for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++) {
      h_spe_charge_avg_pmttype[ipmttype] = (TH1D*)spefile->Get("h_spe_onoff_1")->Clone();
      h_spe_charge_avg_pmttype[ipmttype]->Reset();
      TString histname = h_spe_charge_avg_pmttype[ipmttype]->GetName();
      histname += "_" + SPEruns[ispe] + Form("_avg_pmttype_%d", ipmttype);
      h_spe_charge_avg_pmttype[ipmttype]->SetName(histname);
      nChannels_pmttype[ipmttype] = 0;
    }

    cout << spefile->GetName() << endl;
    
    for (int ichannel=0; ichannel<t_pc2pe->GetEntries(); ichannel++) {
      t_pc2pe->GetEntry(ichannel);

      if (ichannel%1000==0) cout << "Channel " << channel << endl;

      // Skip bad
      if (pc2pe_bad) continue;

      TH1D *h_spe_charge = (TH1D*)spefile->Get(Form("h_spe_onoff_%d", channel));
      if (!h_spe_charge || h_spe_charge->IsZombie()) continue;

      int ipmttype;
      for (ipmttype=0; ipmttype<nPMTtypes; ipmttype++)
	if (pmtflag==PMTflags[ipmttype]) break;

      double pc2pe_mean = h_pc2pe_group_mean->GetBinContent(group+1);
      double pc2pe_rms = h_pc2pe_group_rms->GetBinContent(group+1);

      int iline = -1;
      if (rationorm < pc2pe_mean - 2*pc2pe_rms) iline = 0;
      else if (fabs(rationorm-pc2pe_mean) < 0.5*pc2pe_rms) iline = 1;
      else if (rationorm > pc2pe_mean + 2*pc2pe_rms) iline = 2;

      TString histname = h_spe_charge->GetTitle();
      histname += " ("+PMTTypeNames[ipmttype]+"), ";
      histname += Form("pc2pe = %.2f #pm %.2f, ", rationorm, rationorm_err);
      histname += Form("<SPE> = %.2f #pm %.2f", spe_means[ispe], spe_meanerrs[ispe]);
      h_spe_charge->SetTitle(histname);

      // Save individual plots for one SPE run only
      if (!ispe) {

	// Skip HK PMTs for RMS study (must correspond to plot_position.C)
	if (iline>=0 && pmtflag!=6) {

	  fPMTlists[iline] << channel << endl;
	  
	  if (c_spe_all[iline]) {
	    c_spe_all[iline]->cd();
	    h_spe_charge->Draw();
	    c_spe_all[iline]->Print(c_spe_all_names[iline]);
	  }
	}

	// Anomalous population
	if (rationorm < f_anomalous->Eval(spe_means[ispe])) {
	  c_spe_anomalous->cd();
	  h_spe_charge->Draw();
	  c_spe_anomalous->Print(c_spe_anomalous_name);
	}
      }

      // Make SPE charge histo has consistent binning before summing
      int nBaseBins = h_spe_charge_avg[ipmttype]->GetNbinsX();
      int nBinsThis = h_spe_charge->GetNbinsX();
      if (nBinsThis > nBaseBins)
	h_spe_charge->Rebin(nBinsThis/nBaseBins);
      
      // Normalize SPE charge histo
      h_spe_charge->Scale(1/h_spe_charge->Integral());

      if (iline>=0 && pmtflag!=6) {      
	h_spe_charge_avg[iline]->Add(h_spe_charge);
	nChannels[iline]++;
      }
      
      h_spe_charge_avg_pmttype[ipmttype]->Add(h_spe_charge);
      nChannels_pmttype[ipmttype]++;
    }

    TCanvas *c_spe_charge_avg = new TCanvas(1);
    TLegend *leg_spe_charge = new TLegend(0.5, 0.6, 0.95, 0.88);
    leg_spe_charge->SetHeader("Group Selection (# of Channels)");
    
    for (int iline=0; iline<nRMSranges; iline++) {
      h_spe_charge_avg[iline]->Scale(1./nChannels[iline]);

      h_spe_charge_avg[iline]->SetTitle("Normalized SPE Charge Distributions (Run "+SPEruns[ispe]+")");
      h_spe_charge_avg[iline]->GetYaxis()->SetTitle(Form("Average Rate"));

      h_spe_charge_avg[iline]->SetLineColor(iline+1);
      h_spe_charge_avg[iline]->SetMarkerColor(iline+1);
      h_spe_charge_avg[iline]->SetLineWidth(2);

      if (!iline) h_spe_charge_avg[iline]->Draw();
      else h_spe_charge_avg[iline]->Draw("same");

      leg_spe_charge->AddEntry(h_spe_charge_avg[iline], RMSrangeString[iline]+Form(" (%d)", nChannels[iline]), "lp");
    }
    leg_spe_charge->Draw();
    c_spe_charge_avg->Print("figures/spe_avg_charge_"+SPEruns[ispe]+".png");

    c_spe_charge_avg->SetLogy(1);
    c_spe_charge_avg->Print("figures/spe_avg_charge_"+SPEruns[ispe]+"_logy.png");

    
    TCanvas *c_spe_charge_avg_pmttype = new TCanvas(1);
    TLegend *leg_spe_charge_pmttype = new TLegend(0.5, 0.6, 0.95, 0.88);
    leg_spe_charge_pmttype->SetHeader("PMT Type (# of Channels)");
    
    for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++) {

      if (!nChannels_pmttype[ipmttype]) continue;
      
      h_spe_charge_avg_pmttype[ipmttype]->Scale(1./nChannels_pmttype[ipmttype]);

      h_spe_charge_avg_pmttype[ipmttype]->SetTitle("Normalized SPE Charge Distributions (Run "+SPEruns[ispe]+")");
      h_spe_charge_avg_pmttype[ipmttype]->GetYaxis()->SetTitle(Form("Average Rate"));

      h_spe_charge_avg_pmttype[ipmttype]->SetLineColor(ipmttype+2);
      h_spe_charge_avg_pmttype[ipmttype]->SetMarkerColor(ipmttype+2);
      h_spe_charge_avg_pmttype[ipmttype]->SetLineWidth(2);

      if (!ipmttype) h_spe_charge_avg_pmttype[ipmttype]->Draw();
      else h_spe_charge_avg_pmttype[ipmttype]->Draw("same");

      leg_spe_charge_pmttype->AddEntry(h_spe_charge_avg_pmttype[ipmttype], PMTTypeNames[ipmttype]+Form(" (%d)", nChannels_pmttype[ipmttype]), "lp");
    }
    leg_spe_charge_pmttype->Draw();
    c_spe_charge_avg_pmttype->Print("figures/spe_avg_charge_pmttype_"+SPEruns[ispe]+".png");

    c_spe_charge_avg_pmttype->SetLogy(1);
    c_spe_charge_avg_pmttype->Print("figures/spe_avg_charge_pmttype_"+SPEruns[ispe]+"_logy.png");

  }

  for (int iline=0; iline<nRMSranges; iline++)
    if (c_spe_all[iline])
      c_spe_all[iline]->Print(c_spe_all_names[iline]+"]");

  c_spe_anomalous->Print(c_spe_anomalous_name+"]");

}
