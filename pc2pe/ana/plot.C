{
  gStyle->SetOptStat(0);
  
  TString datadir = "../output/";

  const int nChannels = 11147;
  
  const int nFiles = 6;

  int Colors[nFiles] = {kBlack, kRed, kBlue, kGray+1, kMagenta, kCyan+1};

  TString FileNames[nFiles] = {
    "pc2pe_tst061892_to_5.root",
    "pc2pe_tst080871_to_5.root",
    "pc2pe_tst080885.root",
    "pc2pe_tst061889.root",
    "pc2pe_tst080877.root",
    "pc2pe_tst080884_and_6.root"
  };

  TString FileTitles[nFiles] = {
    "SK4 Low", // (61892-61895)",
    "SK5 Low", // (80871-80875)",
    "SK5 Low Inv.", // (80885)",
    "SK4 High", // (61889)",
    "SK5 High", // (80877)",
    "SK5 High Inv." // (80884, 80886)"
  };

  TString TreeVarNames[nFiles/2] = {
    "_sk4", "_sk5", "_sk5i"
  };

  // Must be propagated to llaser_qb_c.cc
  float ontime_window[nFiles][2] = {
    1180, 1400,
    1000, 1300,
    1000, 1300,
    1140, 1200,
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
    420, 483
    //400, 800,
    //410, 900,
    //410, 900,
    //400, 750,
    //410, 600,
    //410, 600    
  };

  float rHitThresholds[nFiles/2] = {0.009, 0.005, 0.005};
  //float QMeanThresholds[nFiles/2] = {22, 22, 22};
  float QMeanThresholds[nFiles/2] = {50, 50, 50};

  TCanvas *c_ttof = new TCanvas("c_ttof","c_ttof", 0, 0, 1200, 600);
  c_ttof->SetLogy(1);
  TLegend *leg_ttof = new TLegend(0.2, 0.5, 0.5, 0.85);

  TCanvas *c_nqisk = new TCanvas("c_nqisk","c_nqisk", 20, 500, 600, 400);
  //c_nqisk->SetLogy(1);
  TLegend *leg_nqisk = new TLegend(0.2, 0.58, 0.75, 0.87);
  leg_nqisk->SetHeader("Dataset (# events, Mean, RMS)");
  
  TCanvas *c_QOnTime = new TCanvas("c_QOnTime","c_QOnTime", 500, 500, 600, 400);
  TLegend *leg_QOnTime = new TLegend(0.2, 0.5, 0.75, 0.87);
  leg_QOnTime->SetHeader("Dataset (# events, Mean, RMS)");
  
  //TCanvas *c_qisk_ton = new TCanvas("c_qisk_ton","c_qisk_ton", 500, 500, 600, 400);
  //TLegend *leg_qisk_ton = new TLegend(0.2, 0.5, 0.75, 0.87);
  //leg_qisk_ton->SetHeader("Dataset (# events, Mean, RMS)");

  //TH1D *h_nhit_occu[nFiles/2];
  TH1D *h_qmean[nFiles/2];
  TH1D *h_rhit_occu[nFiles/2];
  
  TH1D *h_group_qmean[nFiles/2];
  TH1D *h_group_rhit_occu[nFiles/2];
  TH1D *h_group_pc2pe[nFiles/2];

  TTree *ConnectionTable[nFiles/2];

  for (int ifile=0; ifile<nFiles; ifile++) {
    
    TFile *infile = new TFile(datadir+FileNames[ifile]);

    TTree *ConnectionTableTmp = (TTree*)infile->Get("ConnectionTable");
    if (ifile<nFiles/2) 
      ConnectionTable[ifile] = (TTree*)ConnectionTableTmp->Clone();
    
    ////////////////////////////////////////////
    // Plot TToF
    c_ttof->cd();

    //ttof->Sumw2();

    ttof->SetLineColor(Colors[ifile]);
    ttof->SetLineWidth(2);
    if (ifile%3==2) ttof->SetLineStyle(2);

    // Normalize to total number of laser trigger events
    int nevents = hnqisk->Integral();
    
    ttof->Scale(1./nevents);
    ttof->GetYaxis()->SetTitle("Hits/event");
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
      if (ifile%3==2) hnqisk->SetLineStyle(2);

      hnqisk->Scale(1./nevents);
      hnqisk->SetTitle("Total Number of Hits;Nqisk;Area Normalized");

      if (!ifile) hnqisk->Draw();
      else hnqisk->Draw("sames");

      hnqisk->GetXaxis()->SetRangeUser(0, 500);
      hnqisk->GetYaxis()->SetRangeUser(0, 0.30);
      
      leg_nqisk->AddEntry(hnqisk, FileTitles[ifile]+Form(" (%d, %.0f, %.0f)", nevents, hnqisk->GetMean(), hnqisk->GetRMS()), "l");
    }


    ////////////////////////////////////////////
    // Plot On-time charge for high intensity only
    if (ifile>=nFiles/2) {
      c_QOnTime->cd();

      hQOnTime->Rebin(5);
      hQOnTime->SetLineColor(Colors[ifile]);
      hQOnTime->SetLineWidth(2);
      if (ifile%3==2) hQOnTime->SetLineStyle(2);

      hQOnTime->Scale(1./nevents);
      hQOnTime->SetTitle("On-time Charge;Charge (pe);Area Normalized");

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
      
      h_qisk_ton->Scale(1./nevents);
      
      h_qisk_ton->SetTitle(FileTitles[ifile]+";PMT Cable;Q_{mean} (pe)");

      h_qmean[ifile-3] = (TH1D*)h_qisk_ton->Clone();
      
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
      if (ifile%3==2) h_nhit_ton->SetLineStyle(2);

      h_nhit_ton->Scale(1./nevents);
      h_nhit_ton->SetTitle(FileTitles[ifile]+";PMT Cable;On-time Mean Hit Rate (/event)");

      //if (ifile==nFiles/2)
      h_nhit_ton->Draw();
      //else h_nhit_ton->Draw("sames");

      //h_nhit_ton->GetXaxis()->SetRangeUser(4e5, 1e6);
      float maxY = h_nhit_ton->GetMaximum();
      h_nhit_ton->GetYaxis()->SetRangeUser(-0.01, maxY);
      
      //leg_nhit_ton->AddEntry(h_nhit_ton, FileTitles[ifile], "l");
      c_nhit_ton->Print("figures/nhit_ton"+TreeVarNames[ifile]+".png");
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
      if (ifile%3==2) h_nhit_toff->SetLineStyle(2);

      h_nhit_toff->Scale(1./nevents);
      h_nhit_toff->SetTitle(FileTitles[ifile]+";PMT Cable;Off-time Mean Hit Rate (/event)");

      //if (ifile==nFiles/2)
      h_nhit_toff->Draw();
      //else h_nhit_toff->Draw("sames");

      //h_nhit_toff->GetXaxis()->SetRangeUser(4e5, 1e6);
      float maxY = h_nhit_toff->GetMaximum();
      h_nhit_toff->GetYaxis()->SetRangeUser(-0.01, maxY);
      
      //leg_nhit_toff->AddEntry(h_nhit_toff, FileTitles[ifile], "l");
      c_nhit_toff->Print("figures/nhit_toff"+TreeVarNames[ifile]+".png");
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
      if (ifile%3==2) h_nhit_ton_minus_toff->SetLineStyle(2);

      h_nhit_ton_minus_toff->SetTitle(FileTitles[ifile]+";PMT Cable;(Ontime - Offtime) Hit Rate (/event)");

      //if (ifile==nFiles/2)
      h_nhit_ton_minus_toff->Draw();
      //else h_nhit_ton_minus_toff->Draw("sames");

      //h_nhit_ton_minus_toff->GetXaxis()->SetRangeUser(4e5, 1e6);
      float maxY = h_nhit_ton_minus_toff->GetMaximum();
      h_nhit_ton_minus_toff->GetYaxis()->SetRangeUser(-0.005, maxY);
      
      //leg_nhit_ton_minus_toff->AddEntry(h_nhit_ton_minus_toff, FileTitles[ifile], "l");
    }

    ///////////////////////////////
    // Occupancy correction
    if (ifile<nFiles/2) {

      TH1D *h_rhit_occu_tmp = (TH1D*)h_nhit_ton_minus_toff->Clone();
      for(int ibin=0 ; ibin<=h_rhit_occu_tmp->GetNbinsX() ; ++ibin) 
	h_rhit_occu_tmp->SetBinContent(ibin, -log(1-h_rhit_occu_tmp->GetBinContent(ibin)));

      h_rhit_occu_tmp->SetTitle(FileTitles[ifile]+" (Occupancy Corrected);PMT Cable;(Ontime - Offtime) Hit Rate (/event)");

      h_rhit_occu[ifile] = (TH1D*)h_rhit_occu_tmp->Clone();
      //h_nhit_occu[ifile] = (TH1D*)h_rhit_occu_tmp->Clone();
      //h_nhit_occu[ifile]->Scale(float(nevents));
      //h_nhit_occu[ifile]->GetYaxis()->SetTitle("Number of Events");
    }    

    //////////////////////////////////////////////////////////////////
    /////////////// Below same but for grouped PMTs //////////////////
    
    ////////////////////////////////////////////
    // Plot On-time charge per group for high intensity only
    if (ifile>=nFiles/2) {

      h_group_qisk_ton->Sumw2();
      
      h_group_qisk_ton->Scale(1./nevents);
      
      h_group_qisk_ton->SetTitle(FileTitles[ifile]+";PMT Group;Q_{mean} (pe)");

      h_group_qmean[ifile-3] = (TH1D*)h_group_qisk_ton->Clone();
      
      //leg_qisk_ton->AddEntry(h_group_qisk_ton, FileTitles[ifile], "l");
    }
    
    ////////////////////////////////////////////
    // Plot On-time hits per group for low intensity only
    if (ifile<nFiles/2) {
      //c_nhit_ton->cd();
      TCanvas *c_nhit_ton = new TCanvas(Form("c_group_nhit_ton_%d",ifile),Form("c_group_nhit_ton_%d",ifile), 500, 500, 600, 400);

      h_group_nhit_ton->Sumw2();
      
      h_group_nhit_ton->SetLineColor(Colors[ifile]);
      h_group_nhit_ton->SetMarkerColor(Colors[ifile]);
      h_group_nhit_ton->SetMarkerSize(0.1);
      h_group_nhit_ton->SetLineWidth(2);
      if (ifile%3==2) h_group_nhit_ton->SetLineStyle(2);

      h_group_nhit_ton->Scale(1./nevents);
      h_group_nhit_ton->SetTitle(FileTitles[ifile]+";PMT Group;On-time Mean Hit Rate (/event)");

      //if (ifile==nFiles/2)
      h_group_nhit_ton->Draw();
      //else h_group_nhit_ton->Draw("sames");

      //h_group_nhit_ton->GetXaxis()->SetRangeUser(4e5, 1e6);
      float maxY = h_group_nhit_ton->GetMaximum();
      h_group_nhit_ton->GetYaxis()->SetRangeUser(-0.01, maxY);
      
      //leg_nhit_ton->AddEntry(h_group_nhit_ton, FileTitles[ifile], "l");
      c_nhit_ton->Print("figures/nhit_ton"+TreeVarNames[ifile]+".png");
    }

    ////////////////////////////////////////////
    // Plot Off-time hits per group for low intensity only
    if (ifile<nFiles/2) {
      //c_nhit_toff->cd();
      TCanvas *c_nhit_toff = new TCanvas(Form("c_group_nhit_toff_%d",ifile),Form("c_group_nhit_toff_%d",ifile), 500, 500, 600, 400);

      h_group_nhit_toff->Sumw2();
      
      h_group_nhit_toff->SetLineColor(Colors[ifile]);
      h_group_nhit_toff->SetMarkerColor(Colors[ifile]);
      h_group_nhit_toff->SetMarkerSize(0.1);
      h_group_nhit_toff->SetLineWidth(2);
      if (ifile%3==2) h_group_nhit_toff->SetLineStyle(2);

      h_group_nhit_toff->Scale(1./nevents);
      h_group_nhit_toff->SetTitle(FileTitles[ifile]+";PMT Group;Off-time Mean Hit Rate (/event)");

      //if (ifile==nFiles/2)
      h_group_nhit_toff->Draw();
      //else h_group_nhit_toff->Draw("sames");

      //h_group_nhit_toff->GetXaxis()->SetRangeUser(4e5, 1e6);
      float maxY = h_group_nhit_toff->GetMaximum();
      h_group_nhit_toff->GetYaxis()->SetRangeUser(-0.01, maxY);
      
      //leg_nhit_toff->AddEntry(h_group_nhit_toff, FileTitles[ifile], "l");
      c_nhit_toff->Print("figures/nhit_toff"+TreeVarNames[ifile]+".png");
    }

    ////////////////////////////////////////////
    // Plot On-time minus Off-time hits per group for low intensity only
    if (ifile<nFiles/2) {

      TCanvas *c_nhit_ton_minus_toff = new TCanvas(Form("c_group_nhit_ton_minus_toff_%d",ifile),Form("c_group_nhit_ton_minus_toff_%d",ifile), 500, 500, 600, 400);

      TH1D *h_group_nhit_ton_minus_toff = (TH1D*)h_group_nhit_ton->Clone();
      h_group_nhit_ton_minus_toff->Add(h_group_nhit_toff, -1);
      
      h_group_nhit_ton_minus_toff->SetLineColor(Colors[ifile]);
      h_group_nhit_ton_minus_toff->SetMarkerColor(Colors[ifile]);
      h_group_nhit_ton_minus_toff->SetMarkerSize(0.1);
      h_group_nhit_ton_minus_toff->SetLineWidth(2);
      if (ifile%3==2) h_group_nhit_ton_minus_toff->SetLineStyle(2);

      h_group_nhit_ton_minus_toff->SetTitle(FileTitles[ifile]+";PMT Group;(Ontime - Offtime) Hit Rate (/event)");

      //if (ifile==nFiles/2)
      h_group_nhit_ton_minus_toff->Draw();
      //else h_group_nhit_ton_minus_toff->Draw("sames");

      //h_group_nhit_ton_minus_toff->GetXaxis()->SetRangeUser(4e5, 1e6);
      float maxY = h_group_nhit_ton_minus_toff->GetMaximum();
      h_group_nhit_ton_minus_toff->GetYaxis()->SetRangeUser(-0.005, maxY);
      
      //leg_nhit_ton_minus_toff->AddEntry(h_group_nhit_ton_minus_toff, FileTitles[ifile], "l");
    }

    ///////////////////////////////
    // Occupancy correction
    if (ifile<nFiles/2) {

      TH1D *h_group_rhit_occu_tmp = (TH1D*)h_group_nhit_ton_minus_toff->Clone();
      //for(int ibin=0 ; ibin<=h_group_rhit_occu_tmp->GetNbinsX() ; ++ibin)
      //h_group_rhit_occu_tmp->SetBinContent(ibin, -log(1-h_group_rhit_occu_tmp->GetBinContent(ibin)));
      
      //h_group_rhit_occu_tmp->SetTitle(FileTitles[ifile]+" (Occupancy Corrected);PMT Group;(Ontime - Offtime) Hit Rate (/event)");

      h_group_rhit_occu[ifile] = (TH1D*)h_group_rhit_occu_tmp->Clone();
      //h_group_nhit_occu[ifile] = (TH1D*)h_group_rhit_occu_tmp->Clone();
      //h_group_nhit_occu[ifile]->Scale(float(nevents));
      //h_group_nhit_occu[ifile]->GetYaxis()->SetTitle("Number of Events");
    }    
  }

  c_ttof->cd();
  leg_ttof->Draw();
  c_ttof->Print("figures/ttof.png");

  c_nqisk->cd();
  leg_nqisk->Draw();
  c_nqisk->Print("figures/nqisk.png");

  c_QOnTime->cd();
  leg_QOnTime->Draw();
  c_QOnTime->Print("figures/QOnTime.png");

  //c_qisk_ton->cd();
  //leg_qisk_ton->Draw();

  for (int ifile=0; ifile<nFiles/2; ifile++) {

    // Hit Rate
    TCanvas *c_rhit_occu = new TCanvas(Form("c_rhit_occu_%d",ifile),Form("c_rhit_occu_%d",ifile), 500, 500, 600, 400);
      
    h_rhit_occu[ifile]->SetLineColor(Colors[ifile]);
    h_rhit_occu[ifile]->SetMarkerColor(Colors[ifile]);
    h_rhit_occu[ifile]->SetMarkerSize(0.1);
    h_rhit_occu[ifile]->SetLineWidth(2);
    if (ifile%3==2) h_rhit_occu[ifile]->SetLineStyle(2);
    h_rhit_occu[ifile]->Draw();

    TLine *l_rHitThresh = new TLine(0, rHitThresholds[ifile], nChannels, rHitThresholds[ifile]);
    l_rHitThresh->SetLineColor(Colors[ifile+3]);
    l_rHitThresh->SetLineWidth(2);
    l_rHitThresh->Draw();

    c_rhit_occu->Print("figures/rhit"+TreeVarNames[ifile]+".png");
    // Number of Hits
    //TCanvas *c_nhit_occu = new TCanvas(Form("c_nhit_occu_%d",ifile),Form("c_nhit_occu_%d",ifile), 500, 500, 600, 400);
    //  
    //h_nhit_occu[ifile]->SetLineColor(Colors[ifile]);
    //h_nhit_occu[ifile]->SetMarkerColor(Colors[ifile]);
    //h_nhit_occu[ifile]->SetMarkerSize(0.1);
    //h_nhit_occu[ifile]->SetLineWidth(2);
    //if (ifile%3==2) h_nhit_occu[ifile]->SetLineStyle(2);
    //h_nhit_occu[ifile]->Draw();

    // Mean Charge    
    TCanvas *c_qmean = new TCanvas(Form("c_qmean_%d",ifile),Form("c_qmean_%d",ifile), 500, 500, 600, 400);
      
    h_qmean[ifile]->SetLineColor(Colors[ifile]);
    h_qmean[ifile]->SetMarkerColor(Colors[ifile]);
    h_qmean[ifile]->SetMarkerSize(0.1);
    h_qmean[ifile]->SetLineWidth(2);
    if (ifile%3==2) h_qmean[ifile]->SetLineStyle(2);
    h_qmean[ifile]->Draw();

    float maxY = h_qmean[ifile]->GetMaximum();
    h_qmean[ifile]->GetYaxis()->SetRangeUser(-10, maxY);
    if (ifile==nFiles/2) h_qmean[ifile]->GetYaxis()->SetRangeUser(-10, 120);

    TLine *l_QMeanThresh = new TLine(0, QMeanThresholds[ifile], nChannels, QMeanThresholds[ifile]);
    l_QMeanThresh->SetLineColor(Colors[ifile+3]);
    l_QMeanThresh->SetLineWidth(2);
    l_QMeanThresh->Draw();

    c_qmean->Print("figures/qmean"+TreeVarNames[ifile]+".png");

    //////////////////////////////////////////////////////////////////
    /////////////// Below same but for grouped PMTs //////////////////
    
    // Hit Rate
    TCanvas *c_group_rhit_occu = new TCanvas(Form("c_group_rhit_occu_%d",ifile),Form("c_group_rhit_occu_%d",ifile), 500, 500, 600, 400);
      
    h_group_rhit_occu[ifile]->SetLineColor(Colors[ifile]);
    h_group_rhit_occu[ifile]->SetMarkerColor(Colors[ifile]);
    h_group_rhit_occu[ifile]->SetMarkerSize(0.1);
    h_group_rhit_occu[ifile]->SetLineWidth(2);
    if (ifile%3==2) h_group_rhit_occu[ifile]->SetLineStyle(2);
    h_group_rhit_occu[ifile]->Draw();

    c_group_rhit_occu->Print("figures/group_rhit"+TreeVarNames[ifile]+".png");
    // Number of Hits
    //TCanvas *c_nhit_occu = new TCanvas(Form("c_nhit_occu_%d",ifile),Form("c_nhit_occu_%d",ifile), 500, 500, 600, 400);
    //  
    //h_group_nhit_occu[ifile]->SetLineColor(Colors[ifile]);
    //h_group_nhit_occu[ifile]->SetMarkerColor(Colors[ifile]);
    //h_group_nhit_occu[ifile]->SetMarkerSize(0.1);
    //h_group_nhit_occu[ifile]->SetLineWidth(2);
    //if (ifile%3==2) h_group_nhit_occu[ifile]->SetLineStyle(2);
    //h_group_nhit_occu[ifile]->Draw();

    // Mean Charge    
    TCanvas *c_group_qmean = new TCanvas(Form("c_group_qmean_%d",ifile),Form("c_group_qmean_%d",ifile), 500, 500, 600, 400);
      
    h_group_qmean[ifile]->SetLineColor(Colors[ifile]);
    h_group_qmean[ifile]->SetMarkerColor(Colors[ifile]);
    h_group_qmean[ifile]->SetMarkerSize(0.1);
    h_group_qmean[ifile]->SetLineWidth(2);
    if (ifile%3==2) h_group_qmean[ifile]->SetLineStyle(2);
    h_group_qmean[ifile]->Draw();

    float maxY = h_group_qmean[ifile]->GetMaximum();
    h_group_qmean[ifile]->GetYaxis()->SetRangeUser(-10, maxY);
    if (ifile==nFiles/2) h_group_qmean[ifile]->GetYaxis()->SetRangeUser(-10, 120);

    c_group_qmean->Print("figures/group_qmean"+TreeVarNames[ifile]+".png");    
  }
  
  TFile *outfile = new TFile("pc2pe_output.root", "RECREATE");

  TTree *t_pc2pe = new TTree("pc2pe","SK pc2pe Ratios");//[nFiles/2];
  
  int group = -1;

  int pmtflag[nFiles/2] = {0};

  int channel;
  float rhit[nFiles/2], qmean[nFiles/2], ratio[nFiles/2], ratio_norm[nFiles/2];

  float rmean[nFiles/2] = {0};

  for (int ifile=0; ifile<nFiles/2; ifile++) {

    ConnectionTable[ifile]->SetBranchAddress("group", &group);
    ConnectionTable[ifile]->SetBranchAddress("pmtflag", &pmtflag[ifile]);

    if (!ifile) {
      t_pc2pe->Branch("channel", &channel, "channel/I");
      t_pc2pe->Branch("group", &group, "group/I");
    }
    t_pc2pe->Branch("rhit"+TreeVarNames[ifile], &rhit[ifile], "rhit"+TreeVarNames[ifile]+"/F");
    t_pc2pe->Branch("qmean"+TreeVarNames[ifile], &qmean[ifile], "qmean"+TreeVarNames[ifile]+"/F");
    t_pc2pe->Branch("ratio"+TreeVarNames[ifile], &ratio[ifile], "ratio"+TreeVarNames[ifile]+"/F");
    t_pc2pe->Branch("rationorm"+TreeVarNames[ifile], &ratio_norm[ifile], "rationorm"+TreeVarNames[ifile]+"/F");

    t_pc2pe->Branch("pmtflag"+TreeVarNames[ifile], &pmtflag[ifile], "pmtflag"+TreeVarNames[ifile]+"/I");

    // First get ratio normalization
    cout << "Bad Channels in Dataset: " << TreeVarNames[ifile].Data() << endl;

    int nGoodChannels = 0;

    for (channel=1; channel<h_rhit_occu[ifile]->GetNbinsX(); channel++) {
      rhit[ifile] = h_rhit_occu[ifile]->GetBinContent(channel);
      qmean[ifile] = h_qmean[ifile]->GetBinContent(channel);

      ratio[ifile] = 1;
      if (rhit[ifile] > rHitThresholds[ifile] && qmean[ifile] > QMeanThresholds[ifile]) {	
	ratio[ifile] = qmean[ifile]/rhit[ifile];
	rmean[ifile] += ratio[ifile];
	nGoodChannels++;
      }

      else cout << channel << endl;
    }
    
    rmean[ifile] /= (float)nGoodChannels;

    cout << "Total bad channels = " << nChannels - 1 - nGoodChannels << endl << endl;


    ////////////////////////////////////////////////////////////
    /////////////// Below for grouped PMTs //////////////////
    h_group_pc2pe[ifile] = (TH1D*)h_group_qmean[ifile]->Clone();
    h_group_pc2pe[ifile]->Divide(h_group_rhit_occu[ifile]);
    int nGroups = 0;
    float rmean_group = 0;

    for (int group=1; group<=h_group_pc2pe[ifile]->GetNbinsX(); group++) {

      double group_ratio = h_group_pc2pe[ifile]->GetBinContent(group);

      if (group_ratio<=0) continue;

      rmean_group += group_ratio;
      nGroups++;	
    }
    rmean_group /= (float)nGroups;
    h_group_pc2pe[ifile]->Scale(1/rmean_group); 
  }

  TCanvas *c_group_pc2pe = new TCanvas("c_group_pc2pe","c_group_pc2pe", 500, 500, 600, 400);
  TLegend *leg_group_pc2pe = new TLegend(0.2, 0.75, 0.88, 0.85);
  leg_group_pc2pe->SetNColumns(3);
  
  for (int ifile=0; ifile<nFiles/2; ifile++) {
    h_group_pc2pe[ifile]->SetLineColor(Colors[ifile]);
    h_group_pc2pe[ifile]->SetLineWidth(2);
    h_group_pc2pe[ifile]->GetYaxis()->SetTitle("Relative Gain (pc2pe ratio)");
    h_group_pc2pe[ifile]->GetYaxis()->SetRangeUser(0.95, 1.07);
    if (!ifile) h_group_pc2pe[ifile]->Draw();
    else h_group_pc2pe[ifile]->Draw("sames");
      
    leg_group_pc2pe->AddEntry(h_group_pc2pe[ifile], FileTitles[ifile].ReplaceAll(" Low",""), "l");
  }

  leg_group_pc2pe->Draw();
  c_group_pc2pe->Print("figures/group_pc2pe.png");    

  //t_pc2pe[ifile] = new TTree("pc2pe", FileTitles[ifile].ReplaceAll(" Low", ""));
  
  // Then fill tree with normalized ratios
  for (channel=1; channel<h_rhit_occu[0]->GetNbinsX(); channel++) {
			      
    for (int ifile=0; ifile<nFiles/2; ifile++) {

      ConnectionTable[ifile]->GetEntry(channel-1);

      rhit[ifile] = h_rhit_occu[ifile]->GetBinContent(channel);
      qmean[ifile] = h_qmean[ifile]->GetBinContent(channel);

      ratio_norm[ifile] = ratio[ifile] = 1;
      
      if (rhit[ifile] > rHitThresholds[ifile] && qmean[ifile] > QMeanThresholds[ifile]) {

	ratio[ifile] = qmean[ifile]/rhit[ifile];
	ratio_norm[ifile] = ratio[ifile]/rmean[ifile];

	if (fabs(ratio_norm[ifile])>1.5 || fabs(ratio_norm[ifile])<0.5) 
	  cout << "WTF " << TreeVarNames[ifile].Data() << " " << channel << " " << ratio_norm[ifile] << " " << ratio[ifile] << " " << rhit[ifile] << " " << qmean[ifile] << endl;

      }
    }
    t_pc2pe->Fill();
  }

  outfile->Write();
}
