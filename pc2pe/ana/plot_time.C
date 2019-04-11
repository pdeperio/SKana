{
  gStyle->SetOptStat(0);
  gStyle->SetTitleX(0.25);
  gStyle->SetPadRightMargin(0.1);

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

  TString TimeAxis = "#splitline{%Y-%m-%d}{%H:%M}";

  // Plot total nqisk for low intensity only
  for (int ifile=0; ifile<nFiles/2; ifile++) {
    
    TFile *infile = new TFile(datadir+FileNames[ifile]);

    TCanvas *c_nqisk = new TCanvas(Form("c_nqisk_vs_time_%d",ifile),Form("c_nqisk_vs_time_%d",ifile), 500, 500, 1000, 400);

    hnqisk_vs_time->SetTitle(FileTitles[ifile]+";GMT;Total Number of Hits");
    hnqisk_vs_time->Draw("colz");
    hnqisk_vs_time->ProfileX()->Draw("SAME");
    
    if (TimeAxis.CompareTo("")) {
      TAxis *xax = hnqisk_vs_time->GetXaxis();
      xax->SetTimeDisplay(1);
      xax->SetLabelOffset(0.02);
      xax->SetTimeFormat(TimeAxis);
      xax->SetTimeOffset(0, "gmt");
    }
    //c_nqisk->Print(Form("figures/nqisk_vs_time_%d.png",ifile));

    // Timing cut on SK5 Low Inv. for laser stability
    if (ifile==2) {
      float TimeCut = 1553583000;
      TLine *l_TimeCut = new TLine(TimeCut, 0, TimeCut, 500);
      l_TimeCut->SetLineColor(kGreen-2);
      l_TimeCut->SetLineStyle(2);
      l_TimeCut->SetLineWidth(2);
      l_TimeCut->Draw();

    }
  }

  
  // Plot On-time charge for high intensity only
  for (int ifile=nFiles/2; ifile<nFiles; ifile++) {
    
    TFile *infile = new TFile(datadir+FileNames[ifile]);

    TCanvas *c_QOnTime = new TCanvas(Form("c_QOnTime_vs_time_%d",ifile),Form("c_QOnTime_vs_time_%d",ifile), 500, 500, 1000, 400);

    hQOnTime_vs_time->SetTitle(FileTitles[ifile]+";GMT;Charge [pC]");
    hQOnTime_vs_time->Draw("colz");
    hQOnTime_vs_time->ProfileX()->Draw("SAME");
    
    if (TimeAxis.CompareTo("")) {
      TAxis *xax = hQOnTime_vs_time->GetXaxis();
      xax->SetTimeDisplay(1);
      xax->SetLabelOffset(0.02);
      xax->SetTimeFormat(TimeAxis);
      xax->SetTimeOffset(0, "gmt");
    }
    //c_QOnTime->Print(Form("figures/QOnTime_vs_time_%d.png",ifile));
  }
}
