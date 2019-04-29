void print_avgerr(TH1D *h_in) {
  double avgErr = 0;
  int nBins = 0;

  // PMT Channel binning
  for (int ibin=2; ibin<=h_in->GetNbinsX(); ibin++) {
    if (h_in->GetBinContent(ibin)>0 && h_in->GetBinError(ibin)>0) {
      avgErr += h_in->GetBinError(ibin)/h_in->GetBinContent(ibin);
      nBins++;
    }
  }
  avgErr /= nBins;
  cout << h_in->GetName() << " avg. rel. err. = " << avgErr << endl;
}

// https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html
TH1D* calc_weighted_average(TH1D *h_in1, TH1D *h_in2) {
  
  TH1D *h_avg = (TH1D*)h_in1->Clone();
  h_avg->SetName(Form("%savg", h_in1->GetName()));
  h_avg->Reset();
  
  for (int ibin=0; ibin<=h_in1->GetNbinsX()+1; ibin++) {

    double err1 = h_in1->GetBinError(ibin);
    double err2 = h_in2->GetBinError(ibin);
    if (err1<=0 || err2<=0 || err1!=err1 || err2!=err2) {
      h_avg->SetBinContent(ibin, 1);
      h_avg->SetBinError(ibin, 0);
      continue; 
    }

    double val1 = h_in1->GetBinContent(ibin);
    double val2 = h_in2->GetBinContent(ibin);
    
    double sum_mean_over_err2 = val1/(err1*err1) + val2/(err2*err2);
    double sum_err2 = 1/((err1*err1) + (err2*err2));
    
    h_avg->SetBinContent(ibin, sum_mean_over_err2/sum_err2);
    h_avg->SetBinError(ibin, sqrt(1/sum_err2));

    cout << val1 << " " << val2 << " " << err1 << " " << err2 << " " << sum_mean_over_err2/sum_err2 << " " << sqrt(1/sum_err2) << endl;
  }

  print_avgerr(h_avg);
  return h_avg;
}

void plot() {  
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;  // Suppress "Info in <TCanvas::Print>" messages

  TString datadir = "../output/";

  const int nChannels = 11147;
  
  const int nFiles = 8;
  
  int Colors[nFiles] = {kBlack, kRed, kBlue, kGreen-2, kGray+1, kMagenta, kCyan+1, kGreen};

  TString FileNames[nFiles] = {
    "pc2pe_tst061892_to_5.root",
    "pc2pe_tst080871_to_5.root",
    "pc2pe_tst080885.root", // Inv
    "pc2pe_tst081028.root", // New
    "pc2pe_tst061889.root",
    "pc2pe_tst080877.root",
    "pc2pe_tst080884_and_6.root", // Inv
    "pc2pe_tst081030.root" // New
  };

  TString FileTitles[nFiles] = {
    "SK4 Low", // (61892-61895)",
    "SK5 Low", // (80871-80875)",
    "SK5 Low Inv.", // (80885)",
    "SK5 Low New", // (80885)",
    "SK4 High", // (61889)",
    "SK5 High", // (80877)",
    "SK5 High Inv.", // (80884, 80886)"
    "SK5 High New" // (80884, 80886)"
  };

  const int nConfigs = nFiles/2 + 1; // +1 for SK5 weighted average
  enum config_enum {sk4, sk5, sk5i, sk5n, sk5avg};
  
  TString TreeVarNames[nConfigs] = {
    "_sk4", "_sk5", "_sk5i", "_sk5n", "_sk5avg"
  };

  float TimeDuration[nFiles] = {
    38201,
    8050,
    18229,
    728,
    1725,
    692,  // Inv
    1722  // New
  };
  
  // Must be propagated to llaser_qb_c.cc
  float ontime_window[nFiles][2] = {
    1180, 1400,
    1000, 1300,
    1000, 1300,
    1000, 1300,
    1140, 1200,
    975, 1038,
    975, 1038,
    975, 1038
  };

  // Must have same length as ontime_window above?
  float offtime_window[nFiles][2] = {
    480, 700,
    450, 750,
    450, 750,
    420, 480,
    420, 483,
    420, 483,
    420, 483
    //400, 800,
    //410, 900,
    //410, 900,
    //400, 750,
    //410, 600,
    //410, 600    
  };

  float rHitThresholds[nFiles/2] = {0.009, 0.005, 0.005, 0.005};
  //float QMeanThresholds[nFiles/2] = {22, 22, 22};
  float QMeanThresholds[nFiles/2] = {50, 50, 50, 50};

  TCanvas *c_ttof = new TCanvas("c_ttof","c_ttof", 0, 0, 1200, 600);
  c_ttof->SetLogy(1);
  TLegend *leg_ttof = new TLegend(0.2, 0.5, 0.5, 0.85);

  TCanvas *c_nqisk = new TCanvas("c_nqisk","c_nqisk", 20, 500, 600, 400);
  //c_nqisk->SetLogy(1);
  TLegend *leg_nqisk = new TLegend(0.2, 0.58, 0.75, 0.87);
  leg_nqisk->SetHeader("Dataset (# events, Mean, RMS)");

  TCanvas *c_nHitsOnTime = new TCanvas("c_nHitsOnTime","c_nHitsOnTime", 20, 500, 600, 400);
  //c_nHitsOnTime->SetLogy(1);
  TLegend *leg_nHitsOnTime = new TLegend(0.2, 0.58, 0.75, 0.87);
  leg_nHitsOnTime->SetHeader("Dataset (# events, Mean, RMS)");

  
  TCanvas *c_QOnTime = new TCanvas("c_QOnTime","c_QOnTime", 500, 500, 600, 400);
  TLegend *leg_QOnTime = new TLegend(0.2, 0.5, 0.75, 0.87);
  leg_QOnTime->SetHeader("Dataset (# events, Mean, RMS)");
  
  //TCanvas *c_qisk_ton = new TCanvas("c_qisk_ton","c_qisk_ton", 500, 500, 600, 400);
  //TLegend *leg_qisk_ton = new TLegend(0.2, 0.5, 0.75, 0.87);
  //leg_qisk_ton->SetHeader("Dataset (# events, Mean, RMS)");

  //TH1D *h_nhit_occu[nFiles/2];
  TH1D *h_qmean[nFiles/2];

  TH1D *h_rhit_occu[nFiles/2];
  TH1D *h_rhit_occu_wtf[nFiles/2]; // Needed to save previous histogram, which doesn't survive later
  
  TTree *ConnectionTable[nConfigs];

  bool bNormByTime = 0;

  for (int ifile=0; ifile<nFiles; ifile++) {
    
    TFile *infile = new TFile(datadir+FileNames[ifile]);
    
    ////////////////////////////////////////////
    // Plot TToF
    c_ttof->cd();

    //ttof->Sumw2();

    ttof->SetLineColor(Colors[ifile]);
    ttof->SetLineWidth(2);
    if (ifile%(nFiles/2)==nFiles/2-1) ttof->SetLineStyle(2);

    // Normalize to total number of laser trigger events
    int nevents = hnqisk->Integral();

    double normalization = nevents;
    TString axis_norm = "event";
    if (bNormByTime) {
      normalization = TimeDuration[ifile];
      axis_norm = "second";
    }
    
    ttof->Scale(1./normalization);
    ttof->GetYaxis()->SetTitle("Hits/"+axis_norm);
    ttof->SetTitle("Hit times and On/Off-time Windows");
    
    if (!ifile) ttof->Draw();
    else {
      ttof->Draw("sames");
      TPaveStats *st = (TPaveStats*)ttof->FindObject("stats");
    }
    
    leg_ttof->AddEntry(ttof, FileTitles[ifile], "l");

    // Set axis ranges
    float ttof_miny = 0.04;
    float ttof_maxy = 7000;
    ttof->GetXaxis()->SetRangeUser(300, 1700);
    ttof->GetYaxis()->SetRangeUser(ttof_miny, ttof_maxy);

    // Time window lines
    for (int iwindow=0; iwindow<2; iwindow++) {      
      TLine *l_OnTimeWindow = new TLine(ontime_window[ifile][iwindow], ttof_miny, ontime_window[ifile][iwindow], ttof_maxy);
      l_OnTimeWindow->SetLineColor(Colors[ifile]);
      l_OnTimeWindow->SetLineStyle(2);
      l_OnTimeWindow->SetLineWidth(2);
      l_OnTimeWindow->Draw();

      TLine *l_OffTimeWindow = new TLine(offtime_window[ifile][iwindow], ttof_miny, offtime_window[ifile][iwindow], ttof_maxy);
      l_OffTimeWindow->SetLineColor(Colors[ifile]);
      l_OffTimeWindow->SetLineStyle(3);
      l_OffTimeWindow->SetLineWidth(2);
      l_OffTimeWindow->Draw();
    }

    ////////////////////////////////////////////
    // Plot total nqisk for low intensity only
    if (ifile<nFiles/2) {
      c_nqisk->cd();
      //hnqisk->Sumw2();
      hnqisk->SetLineColor(Colors[ifile]);
      hnqisk->SetLineWidth(2);

      hnqisk->Scale(1./normalization);
      hnqisk->SetTitle("Total Number of Hits;Nqisk;/"+axis_norm);

      if (!ifile) hnqisk->Draw();
      else hnqisk->Draw("sames");

      hnqisk->GetXaxis()->SetRangeUser(0, 500);
      hnqisk->GetYaxis()->SetRangeUser(0, 0.41);
      
      leg_nqisk->AddEntry(hnqisk, FileTitles[ifile]+Form(" (%d, %.0f, %.0f)", nevents, hnqisk->GetMean(), hnqisk->GetRMS()), "l");
    }

    
    ////////////////////////////////////////////
    // Plot total on-time "nqisk" for low intensity only
    if (ifile<nFiles/2) {
      c_nHitsOnTime->cd();
      //hnHitsOnTime->Sumw2();
      hnHitsOnTime->SetLineColor(Colors[ifile]);
      hnHitsOnTime->SetLineWidth(2);

      hnHitsOnTime->Scale(1./normalization);
      hnHitsOnTime->SetTitle("Total Number of On-Time Hits;On-time Hits;/"+axis_norm);

      if (!ifile) hnHitsOnTime->Draw();
      else hnHitsOnTime->Draw("sames");

      hnHitsOnTime->GetXaxis()->SetRangeUser(0, 450);
      hnHitsOnTime->GetYaxis()->SetRangeUser(0, 0.39);
      
      leg_nHitsOnTime->AddEntry(hnHitsOnTime, FileTitles[ifile]+Form(" (%d, %.0f, %.0f)", nevents, hnHitsOnTime->GetMean(), hnHitsOnTime->GetRMS()), "l");
    }



    ////////////////////////////////////////////
    // Plot On-time charge for high intensity only
    if (ifile>=nFiles/2) {
      c_QOnTime->cd();

      hQOnTime->Rebin(5);
      hQOnTime->SetLineColor(Colors[ifile]);
      hQOnTime->SetLineWidth(2);

      hQOnTime->Scale(1./normalization);
      hQOnTime->SetTitle("On-time Charge;Charge (pe);/"+axis_norm);

      if (ifile==nFiles/2) hQOnTime->Draw();
      else hQOnTime->Draw("sames");

      //hQOnTime->GetXaxis()->SetRangeUser(4e5, 1.2e6);
      hQOnTime->GetYaxis()->SetRangeUser(0, 0.17);

      leg_QOnTime->AddEntry(hQOnTime, FileTitles[ifile]+Form(" (%d, %.0f, %.0f)", nevents, hQOnTime->GetMean(), hQOnTime->GetRMS()), "l");
    }

    ////////////////////////////////////////////
    // Plot On-time charge per channel for high intensity only
    if (ifile>=nFiles/2) {

      h_qisk_ton->Sumw2();
      
      h_qisk_ton->Scale(1./normalization);
      
      h_qisk_ton->SetTitle(FileTitles[ifile]+";PMT Cable;Q_{mean} (pe)");

      h_qmean[ifile-nFiles/2] = (TH1D*)h_qisk_ton->Clone();

      h_qmean[ifile-nFiles/2]->SetName("h_qisk_ton"+TreeVarNames[ifile-nFiles/2]);
      //leg_qisk_ton->AddEntry(h_qisk_ton, FileTitles[ifile], "l");
    }
    
    ////////////////////////////////////////////
    // Plot On-time hits per channel for low intensity only
    if (ifile<nFiles/2) {
      //c_nhit_ton->cd();
      TCanvas *c_nhit_ton = new TCanvas(Form("c_nhit_ton_%d",ifile),Form("c_nhit_ton_%d",ifile), 500, 500, 600, 400);

      h_nhit_ton->Sumw2();
      
      h_nhit_ton->SetLineColor(Colors[ifile]);
      h_nhit_ton->SetMarkerColor(Colors[ifile]);
      h_nhit_ton->SetMarkerSize(0.1);
      h_nhit_ton->SetLineWidth(2);

      h_nhit_ton->Scale(1./normalization);
      h_nhit_ton->SetTitle(FileTitles[ifile]+";PMT Cable;On-time Mean Hit Rate (/"+axis_norm+")");
      //TString histname = h_nhit_ton->GetName()+TreeVarNames[ifile];
      //h_nhit_ton->SetName(histname);
      
      //if (ifile==nFiles/2)
      h_nhit_ton->Draw();
      //else h_nhit_ton->Draw("sames");

      //h_nhit_ton->GetXaxis()->SetRangeUser(4e5, 1e6);
      float maxY = h_nhit_ton->GetMaximum();
      h_nhit_ton->GetYaxis()->SetRangeUser(-0.01, maxY);
      
      //leg_nhit_ton->AddEntry(h_nhit_ton, FileTitles[ifile], "l");
      c_nhit_ton->Print("figures/nhit_ton"+TreeVarNames[ifile]+".png");

      //print_avgerr(h_nhit_ton);
    }

    ////////////////////////////////////////////
    // Plot Off-time hits per channel for low intensity only
    if (ifile<nFiles/2) {
      //c_nhit_toff->cd();
      TCanvas *c_nhit_toff = new TCanvas(Form("c_nhit_toff_%d",ifile),Form("c_nhit_toff_%d",ifile), 500, 500, 600, 400);

      h_nhit_toff->Sumw2();
      
      h_nhit_toff->SetLineColor(Colors[ifile]);
      h_nhit_toff->SetMarkerColor(Colors[ifile]);
      h_nhit_toff->SetMarkerSize(0.1);
      h_nhit_toff->SetLineWidth(2);

      h_nhit_toff->Scale(1./normalization);
      h_nhit_toff->SetTitle(FileTitles[ifile]+";PMT Cable;Off-time Mean Hit Rate (/"+axis_norm+")");
      //TString histname = h_nhit_toff->GetName()+TreeVarNames[ifile];
      //h_nhit_toff->SetName(histname);

      //if (ifile==nFiles/2)
      h_nhit_toff->Draw();
      //else h_nhit_toff->Draw("sames");

      //h_nhit_toff->GetXaxis()->SetRangeUser(4e5, 1e6);
      float maxY = h_nhit_toff->GetMaximum();
      h_nhit_toff->GetYaxis()->SetRangeUser(-0.01, maxY);
      
      //leg_nhit_toff->AddEntry(h_nhit_toff, FileTitles[ifile], "l");
      c_nhit_toff->Print("figures/nhit_toff"+TreeVarNames[ifile]+".png");
      //print_avgerr(h_nhit_toff);
    }

    ////////////////////////////////////////////
    // Plot On-time minus Off-time hits per channel for low intensity only
    if (ifile<nFiles/2) {

      TCanvas *c_nhit_ton_minus_toff = new TCanvas(Form("c_nhit_ton_minus_toff_%d",ifile),Form("c_nhit_ton_minus_toff_%d",ifile), 500, 500, 600, 400);

      TH1D *h_nhit_ton_minus_toff = (TH1D*)h_nhit_ton->Clone();
      h_nhit_ton_minus_toff->Add(h_nhit_toff, -1);
      
      h_nhit_ton_minus_toff->SetLineColor(Colors[ifile]);
      h_nhit_ton_minus_toff->SetMarkerColor(Colors[ifile]);
      h_nhit_ton_minus_toff->SetMarkerSize(0.1);
      h_nhit_ton_minus_toff->SetLineWidth(2);

      h_nhit_ton_minus_toff->SetTitle(FileTitles[ifile]+";PMT Cable;(Ontime - Offtime) Hit Rate (/"+axis_norm+")");
      TString histname = h_nhit_ton_minus_toff->GetName()+TreeVarNames[ifile];
      h_nhit_ton_minus_toff->SetName(histname);

      //if (ifile==nFiles/2)
      h_nhit_ton_minus_toff->Draw();
      //else h_nhit_ton_minus_toff->Draw("sames");

      //h_nhit_ton_minus_toff->GetXaxis()->SetRangeUser(4e5, 1e6);
      float maxY = h_nhit_ton_minus_toff->GetMaximum();
      h_nhit_ton_minus_toff->GetYaxis()->SetRangeUser(-0.005, maxY);
      
      //leg_nhit_ton_minus_toff->AddEntry(h_nhit_ton_minus_toff, FileTitles[ifile], "l");

      //print_avgerr(h_nhit_ton_minus_toff);
    }

    ///////////////////////////////
    // Occupancy correction
    if (ifile<nFiles/2) {

      TH1D *h_rhit_occu_tmp = (TH1D*)h_nhit_ton_minus_toff->Clone();
      for(int ibin=0 ; ibin<=h_rhit_occu_tmp->GetNbinsX() ; ++ibin) {

	double ValBeforeCorr = h_rhit_occu_tmp->GetBinContent(ibin);
	double ValAfterCorr = -log(1 - ValBeforeCorr);
	h_rhit_occu_tmp->SetBinContent(ibin, ValAfterCorr);

	// Not perfect error propogation
	h_rhit_occu_tmp->SetBinError(ibin, h_rhit_occu_tmp->GetBinError(ibin) * ValAfterCorr/ValBeforeCorr);
      }

      h_rhit_occu_tmp->SetTitle(FileTitles[ifile]+" (Occupancy Corrected);PMT Cable;(Ontime - Offtime) Hit Rate (/"+axis_norm+")");

      h_rhit_occu_wtf[ifile] = (TH1D*)h_rhit_occu_tmp->Clone();
      h_rhit_occu_wtf[ifile]->SetName("h_rhit_occu"+TreeVarNames[ifile]);
      //h_nhit_occu[ifile]->Scale(float(normalization));
      //h_nhit_occu[ifile]->GetYaxis()->SetTitle("Number of Events");
      h_rhit_occu[ifile] = h_rhit_occu_wtf[ifile];

      //print_avgerr(h_rhit_occu[ifile]);
    }
  }

  c_ttof->cd();
  leg_ttof->Draw();
  c_ttof->Print("figures/ttof.png");

  c_nqisk->cd();
  leg_nqisk->Draw();
  c_nqisk->Print("figures/nqisk.png");

  c_nHitsOnTime->cd();
  leg_nHitsOnTime->Draw();
  c_nHitsOnTime->Print("figures/nHitsOnTime.png");
  
  c_QOnTime->cd();
  leg_QOnTime->Draw();
  c_QOnTime->Print("figures/QOnTime.png");

  //c_qisk_ton->cd();
  //leg_qisk_ton->Draw();

  TH1D *h_qmean_over_rhit[nConfigs];
  TH1D *h_pc2pe[nConfigs];

  for (int ifile=0; ifile<nFiles/2; ifile++) {

    // Hit Rate
    TCanvas *c_rhit_occu = new TCanvas(Form("c_rhit_occu_%d",ifile),Form("c_rhit_occu_%d",ifile), 500, 500, 600, 400);
    
    h_rhit_occu[ifile]->SetLineColor(Colors[ifile]);
    h_rhit_occu[ifile]->SetMarkerColor(Colors[ifile]);
    h_rhit_occu[ifile]->SetMarkerSize(0.1);
    h_rhit_occu[ifile]->SetLineWidth(2);

    h_rhit_occu[ifile]->Draw();

    TLine *l_rHitThresh = new TLine(0, rHitThresholds[ifile], nChannels, rHitThresholds[ifile]);
    l_rHitThresh->SetLineColor(Colors[ifile+nFiles/2]);
    l_rHitThresh->SetLineWidth(2);
    l_rHitThresh->Draw();

    c_rhit_occu->Print("figures/rhit"+TreeVarNames[ifile]+".png");

    print_avgerr(h_rhit_occu[ifile]);
    
    // Number of Hits
    //TCanvas *c_nhit_occu = new TCanvas(Form("c_nhit_occu_%d",ifile),Form("c_nhit_occu_%d",ifile), 500, 500, 600, 400);
    //  
    //h_nhit_occu[ifile]->SetLineColor(Colors[ifile]);
    //h_nhit_occu[ifile]->SetMarkerColor(Colors[ifile]);
    //h_nhit_occu[ifile]->SetMarkerSize(0.1);
    //h_nhit_occu[ifile]->SetLineWidth(2);
    //if (ifile%(nFiles/2)==nFiles/2-1) h_nhit_occu[ifile]->SetLineStyle(2);
    //h_nhit_occu[ifile]->Draw();

    // Mean Charge    
    TCanvas *c_qmean = new TCanvas(Form("c_qmean_%d",ifile),Form("c_qmean_%d",ifile), 500, 500, 600, 400);
      
    h_qmean[ifile]->SetLineColor(Colors[ifile]);
    h_qmean[ifile]->SetMarkerColor(Colors[ifile]);
    h_qmean[ifile]->SetMarkerSize(0.1);
    h_qmean[ifile]->SetLineWidth(2);
    h_qmean[ifile]->Draw();

    float maxY = h_qmean[ifile]->GetMaximum();
    h_qmean[ifile]->GetYaxis()->SetRangeUser(-10, maxY);
    if (ifile==nFiles/2) h_qmean[ifile]->GetYaxis()->SetRangeUser(-10, 120);

    TLine *l_QMeanThresh = new TLine(0, QMeanThresholds[ifile], nChannels, QMeanThresholds[ifile]);
    l_QMeanThresh->SetLineColor(Colors[ifile+nFiles/2]);
    l_QMeanThresh->SetLineWidth(2);
    l_QMeanThresh->Draw();

    c_qmean->Print("figures/qmean"+TreeVarNames[ifile]+".png");

    print_avgerr(h_qmean[ifile]);

    // Take temporary division
    h_qmean_over_rhit[ifile] = (TH1D*)h_qmean[ifile]->Clone();
    h_qmean_over_rhit[ifile]->SetName("h_qmean_over_rhit"+TreeVarNames[ifile]);
    h_qmean_over_rhit[ifile]->Divide(h_rhit_occu[ifile]);
    print_avgerr(h_qmean_over_rhit[ifile]);

  }


  // Read official bad channel list
  vector<int> badchannels;
  ifstream fbadchannel;
  fbadchannel.open("const/badchannels_sk5.txt");
  while (!fbadchannel.eof()){
    int badchannel;
    fbadchannel >> badchannel;
    badchannels.push_back(badchannel);
  }

  // Read old official SK4 table
  vector<double> sk4_pgain;
  ifstream fsk4pgain;
  fsk4pgain.open("const/pgain_30.00");

  string line;
  getline(fsk4pgain, line);

  while (!fsk4pgain.eof()){
    int channel;
    double pgain;
    fsk4pgain >> channel >> pgain;
    sk4_pgain.push_back(pgain);
  }
  
  // Read QE table
  vector<double> sk4_qe;
  ifstream fsk4qe;
  fsk4qe.open("const/qetable4_1.dat");

  while (!fsk4qe.eof()){
    int channel;
    double qe;
    fsk4qe >> channel >> qe;
    sk4_qe.push_back(qe);
  }


  if (bNormByTime) break;
  
  TFile *outfile = new TFile("pc2pe_output.root", "RECREATE");

  TTree *t_pc2pe = new TTree("pc2pe","SK pc2pe Ratios");//[nFiles/2];
  
  int group = -1;

  int pmtflag[nConfigs] = {0};

  int channel;
  float rhit[nConfigs], qmean[nConfigs], ratio[nConfigs], ratio_norm[nConfigs], ratio_norm_sk4;
  float rhit_err[nConfigs], qmean_err[nConfigs], ratio_err[nConfigs], ratio_norm_err[nConfigs];
  int pc2pe_bad[nConfigs]={0}, badchannel;
  int pc2pe_bad_sk4;
  float qe_sk4;
  int qe_bad_sk4;
  
  float rmean[nConfigs] = {0};
  
  for (int ifile=0; ifile<nConfigs; ifile++) {

    TFile *infile = new TFile(datadir+FileNames[ifile]);

    TTree *ConnectionTableTmp = (TTree*)infile->Get("ConnectionTable");

    if (ifile<nConfigs) {
      if (ifile==sk5avg) {
	TFile *infile_sk5avg = new TFile(datadir+"pc2pe_tst081028_kludge4avg.root");
	TTree *ConnectionTableTmp = (TTree*)infile_sk5avg->Get("ConnectionTable");
      } 
      ConnectionTable[ifile] = (TTree*)ConnectionTableTmp->Clone();
    }

    ConnectionTable[ifile]->SetBranchAddress("group", &group);
    ConnectionTable[ifile]->SetBranchAddress("pmtflag", &pmtflag[ifile]);

    if (!ifile) {
      t_pc2pe->Branch("channel", &channel, "channel/I");
      t_pc2pe->Branch("group", &group, "group/I");
      t_pc2pe->Branch("badchannel", &badchannel, "badchannel/I");
      t_pc2pe->Branch("rationorm_sk4official", &ratio_norm_sk4, "rationorm_sk4official/F");
      t_pc2pe->Branch("pmtflag_sk4official", &pmtflag[ifile], "pmtflag_sk4official/I");
      t_pc2pe->Branch("pc2pe_bad_sk4official", &pc2pe_bad_sk4, "pc2pe_bad_sk4official/I");
      t_pc2pe->Branch("qe_sk4", &qe_sk4, "qe_sk4/F");
      t_pc2pe->Branch("qe_bad_sk4", &qe_bad_sk4, "qe_bad_sk4/I");
    }
    //t_pc2pe->Branch("rhit"+TreeVarNames[ifile], &rhit[ifile], "rhit"+TreeVarNames[ifile]+"/F");
    //t_pc2pe->Branch("qmean"+TreeVarNames[ifile], &qmean[ifile], "qmean"+TreeVarNames[ifile]+"/F");
    //t_pc2pe->Branch("ratio"+TreeVarNames[ifile], &ratio[ifile], "ratio"+TreeVarNames[ifile]+"/F");
    t_pc2pe->Branch("rationorm"+TreeVarNames[ifile], &ratio_norm[ifile], "rationorm"+TreeVarNames[ifile]+"/F");
    t_pc2pe->Branch("pc2pe_bad"+TreeVarNames[ifile], &pc2pe_bad[ifile], "pc2pe_bad"+TreeVarNames[ifile]+"/I");

    //t_pc2pe->Branch("rhit_err"+TreeVarNames[ifile], &rhit_err[ifile], "rhit_err"+TreeVarNames[ifile]+"/F");
    //t_pc2pe->Branch("qmean_err"+TreeVarNames[ifile], &qmean_err[ifile], "qmean_err"+TreeVarNames[ifile]+"/F");
    //t_pc2pe->Branch("ratio_err"+TreeVarNames[ifile], &ratio_err[ifile], "ratio_err"+TreeVarNames[ifile]+"/F");
    t_pc2pe->Branch("rationorm_err"+TreeVarNames[ifile], &ratio_norm_err[ifile], "rationorm_err"+TreeVarNames[ifile]+"/F");

    t_pc2pe->Branch("pmtflag"+TreeVarNames[ifile], &pmtflag[ifile], "pmtflag"+TreeVarNames[ifile]+"/I");

    //cout << "Bad Channels in Dataset: " << TreeVarNames[ifile].Data() << endl;

    // First get ratio normalization

    // Except for later weighted average
    if (ifile == sk5avg) {
      h_pc2pe[ifile] = (TH1D*)h_qmean_over_rhit[sk5]->Clone();
      h_pc2pe[ifile]->Reset();
      h_pc2pe[ifile]->SetName("h_pc2pe"+TreeVarNames[ifile]);
      continue;
    }
    
    int nGoodChannels = 0;
    
    for (int ibin=2; ibin<=h_rhit_occu[ifile]->GetNbinsX(); ibin++) {

      ConnectionTable[ifile]->GetEntry(ibin-2);

      rhit[ifile] = h_rhit_occu[ifile]->GetBinContent(ibin);
      qmean[ifile] = h_qmean[ifile]->GetBinContent(ibin);

      if (rhit[ifile] > rHitThresholds[ifile] && qmean[ifile] > QMeanThresholds[ifile]
	  //&& (pmtflag[ifile]==3 || pmtflag[ifile]==4)  // For ignoring HK PMTs
	  ) {	
	ratio[ifile] = h_qmean_over_rhit[ifile]->GetBinContent(ibin);
	rmean[ifile] += ratio[ifile];

	nGoodChannels++;
      }
      //else cout << ibin-1 << endl;
    }
    rmean[ifile] /= (float)nGoodChannels;

    h_pc2pe[ifile] = (TH1D*)h_qmean_over_rhit[ifile]->Clone();
    h_pc2pe[ifile]->SetName("h_pc2pe"+TreeVarNames[ifile]);
    h_pc2pe[ifile]->Scale(1/rmean[ifile]);
    
    //cout << TreeVarNames[ifile] << " rmean = " << rmean[ifile] << endl;
  }
  
  for (int ibin=2; ibin<=h_rhit_occu[0]->GetNbinsX(); ibin++) {

    channel = ibin-1;
    
    for (int ifile=0; ifile<nConfigs; ifile++) {

      ConnectionTable[ifile]->GetEntry(channel-1);

      ratio_norm[ifile] = ratio[ifile] = 1;
      ratio_norm_err[ifile] = ratio_err[ifile] = 0;

      if (ifile != sk5avg) {

	rhit[ifile] = h_rhit_occu[ifile]->GetBinContent(ibin);
	qmean[ifile] = h_qmean[ifile]->GetBinContent(ibin);
	
	// Good channel selection by threshold
	if (rhit[ifile] > rHitThresholds[ifile] && qmean[ifile] > QMeanThresholds[ifile]
	    //&& (pmtflag[ifile]==3 || pmtflag[ifile]==4)  // For ignoring HK PMTs
	    ) {

	  pc2pe_bad[ifile] = 0;

	  ratio[ifile] = h_qmean_over_rhit[ifile]->GetBinContent(ibin);
	  ratio_norm[ifile] = h_pc2pe[ifile]->GetBinContent(ibin);
	  
	  ratio_err[ifile] = h_qmean_over_rhit[ifile]->GetBinError(ibin);
	  ratio_norm_err[ifile] = h_pc2pe[ifile]->GetBinError(ibin);

	  if (fabs(ratio_norm[ifile])>1.5 || fabs(ratio_norm[ifile])<0.5) 
	    cout << "Outlier: " << TreeVarNames[ifile].Data() << " " << channel << " " << ratio_norm[ifile] << " " << ratio[ifile] << " " << rhit[ifile] << " " << qmean[ifile] << endl;
	}

	else {
	  pc2pe_bad[ifile] = 1;
	  h_pc2pe[ifile]->SetBinContent(ibin, ratio_norm[ifile]);
	  h_pc2pe[ifile]->SetBinError(ibin, ratio_norm_err[ifile]);
	}
      }

      // SK5 Weighted average
      else {

	pc2pe_bad[ifile] = pc2pe_bad[sk5] && pc2pe_bad[sk5n];
	
	double err1 = ratio_norm_err[sk5];
	double err2 = ratio_norm_err[sk5n];
	if (!pc2pe_bad[ifile] && err1>0 && err2>0) {
	
	  double val1 = ratio_norm[sk5];
	  double val2 = ratio_norm[sk5n];
	
	  double sum_mean_over_err2 = val1/(err1*err1) + val2/(err2*err2);
	  double sum_err2 = 1/(err1*err1) + 1/(err2*err2);
	
	  ratio_norm[ifile] = sum_mean_over_err2/sum_err2;
	  ratio_norm_err[ifile] = sqrt(1/sum_err2);

	  h_pc2pe[ifile]->SetBinContent(ibin, ratio_norm[ifile]);
	  h_pc2pe[ifile]->SetBinError(ibin, ratio_norm_err[ifile]);
	}
      }
    }

    // Official bad channels
    //if(std::find(badchannels.begin(), badchannels.end(), channel) != badchannels.end())
    badchannel = 0;
    for (int ib=0; ib<badchannels.size(); ib++) {
      if (badchannels[ib] < channel) continue;
      else if (channel == badchannels[ib]) {
	badchannel = 1;
	break;
      } else {
	break;
      }
    }

    // Official SK3/SK4 table
    ratio_norm_sk4 = sk4_pgain[channel-1];
    pc2pe_bad_sk4 = 0;
    if (sk4_pgain[channel-1] == 1.00000000000) pc2pe_bad_sk4 = 1;    

    qe_sk4 = sk4_qe[channel-1];
    qe_bad_sk4 = 0;
    if (qe_sk4 == 1) qe_bad_sk4 = 1;    

    t_pc2pe->Fill();
  }
  
  for (int ifile=0; ifile<nConfigs; ifile++)
    print_avgerr(h_pc2pe[ifile]);

  outfile->Write();
}
