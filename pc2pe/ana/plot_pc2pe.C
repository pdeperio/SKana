{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.15);

  const int nFiles = 6;

  TString TreeVarNames[nFiles] = {
    "_sk4",
    "_sk5",
    "_sk5i",
    "_sk5n",
    "_sk5avg",
    "_qe_sk4" // QE
    //rationorm_sk4official
  }
  int Colors[nFiles] = {
    kGray+2,
    kRed,
    kBlue,
    kGreen-2,
    kBlack,
    kMagenta
  };

  TString FileTitles[nFiles] = {
    "SK4",
    "SK5",
    "SK5 Inv.",
    "SK5 New",
    "SK5 Avg.",
    "SK4 QE"
    //"SK4 official"
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
    TString histtitle = ";pc2pe Ratio";
    if (ifile==nFiles-1) histtitle = ";SK4 QE";
    histtitle += ";Number of Channels";
    
    h_pc2pe[ifile] = new TH1D(histname, histtitle, 50, 0.1, 2);

    // Select good PMT flags
    int useFile = ifile;
    if (ifile==nFiles-1) useFile = 0;
    TString cut_pmttype = "(";
    for (int jpmttype=0; jpmttype<nPMTtypes-1; jpmttype++) {
      cut_pmttype += "pmtflag"+TreeVarNames[useFile]+Form("==%d",PMTflags[jpmttype]);
      if (jpmttype<nPMTtypes-2) cut_pmttype += " || ";
    }
    cut_pmttype += ")";

    if (ifile==nFiles-1) cut_pmttype += " && qe_bad_sk4 == 0";
    else cut_pmttype += " && pc2pe_bad"+TreeVarNames[ifile]+" == 0";

    TString plotVar = "rationorm"+TreeVarNames[ifile];
    if (ifile==nFiles-1) plotVar = "qe_sk4";
    
    pc2pe->Project(histname, plotVar, cut_pmttype);

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
      TString histtitle = ";pc2pe Ratio";
      if (ifile==nFiles-1) histtitle = ";SK4 QE";
      histtitle += ";Number of Channels";
      
      h_pc2pe_pmttypes[ifile][ipmttype] = new TH1D(histname, histtitle, 50, 0.1, 2);
      
      int useFile = ifile;
      if (ifile==nFiles-1) useFile = 0;
      TString cut_pmttype = "pmtflag"+TreeVarNames[useFile]+Form("==%d",PMTflags[ipmttype]);
      if (ipmttype==nPMTtypes-1) {
	cut_pmttype = "";
	for (int jpmttype=0; jpmttype<nPMTtypes-1; jpmttype++) {
	  cut_pmttype += "pmtflag"+TreeVarNames[useFile]+Form("!=%d",PMTflags[jpmttype]);
	  if (jpmttype<nPMTtypes-2) cut_pmttype += " && ";
	}
      }
      if (ifile==nFiles-1) cut_pmttype += " && qe_bad_sk4 == 0";
      else cut_pmttype += " && pc2pe_bad"+TreeVarNames[ifile]+" == 0";

      
      TString plotVar = "rationorm"+TreeVarNames[ifile];
      if (ifile==nFiles-1) plotVar = "qe_sk4";
      pc2pe->Project(histname, plotVar, cut_pmttype);

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
  
  for (int ifile=0; ifile<nFiles; ifile++) {
    for (int jfile=ifile+1; jfile<nFiles; jfile++) {

      TCanvas *c_pc2pe_2d = new TCanvas(1);

      TString histname = "pc2pe_2d"+TreeVarNames[ifile]+"_vs"+TreeVarNames[jfile];
      TString xTitle = "pc2pe "+FileTitles[ifile];
      if (ifile==nFiles-1) xTitle = "QE SK4";
      TString yTitle = "pc2pe "+FileTitles[jfile];
      if (jfile==nFiles-1) yTitle = "QE SK4";

      h_pc2pe_2d[ifile] = new TH2D(histname, ";"+xTitle+";"+yTitle+";Number of Channels", 120, 0.1, 2, 120, 0.1, 2);

      TString plotVarX = "rationorm"+TreeVarNames[ifile];
      if (ifile==nFiles-1) plotVarX = "qe_sk4";
      TString plotVarY = "rationorm"+TreeVarNames[jfile];
      if (jfile==nFiles-1) plotVarY = "qe_sk4";

      TString cut_badpmt = "";
      if (ifile==nFiles-1) cut_badpmt = "qe_bad_sk4 == 0";
      else cut_badpmt = "pc2pe_bad"+TreeVarNames[ifile]+" == 0";
      if (jfile==nFiles-1) cut_badpmt += "&& qe_bad_sk4 == 0";
      else cut_badpmt += " && pc2pe_bad"+TreeVarNames[ifile]+" == 0";

      
      pc2pe->Project(histname, plotVarY+":"+plotVarX, cut_badpmt);

      h_pc2pe_2d[ifile]->Draw("colz");

      h_pc2pe_2d[ifile]->SetTitle(Form("Correlation = %.2f", h_pc2pe_2d[ifile]->GetCorrelationFactor()));

      TPaveText *tRMSoverMeans = new TPaveText(0.2, 1.38, 1, 1.9);
      tRMSoverMeans->SetFillColorAlpha(0, 0);
      tRMSoverMeans->SetBorderSize(1);
      tRMSoverMeans->AddText("RMS/Mean (%)");
      
      double RMS2d, Mean2d, RMSErr2d, MeanErr2d;
      for (int iaxis=0; iaxis<=1; iaxis++) {
	Mean2d = h_pc2pe_2d[ifile]->GetMean(iaxis+1);
	MeanErr2d = h_pc2pe_2d[ifile]->GetMeanError(iaxis+1);
	RMS2d = h_pc2pe_2d[ifile]->GetRMS(iaxis+1);
	RMSErr2d = h_pc2pe_2d[ifile]->GetRMSError(iaxis+1);
	
	double RMSoverMean = RMS2d/Mean2d*100;
	double RMSoverMeanErr = RMSoverMean*sqrt(pow(MeanErr2d/Mean2d, 2)+pow(RMSErr2d/RMS2d, 2));
	
	TString AxisTitle = yTitle;
	if (!iaxis) AxisTitle = xTitle;
	tRMSoverMeans->AddText(AxisTitle+Form(": %.2f #pm %.2f", RMSoverMean, RMSoverMeanErr));
      }
      tRMSoverMeans->Draw();

      c_pc2pe_2d->Print("figures/pc2pe"+TreeVarNames[ifile]+"_vs"+TreeVarNames[jfile]+".png");
    }
  }

  // Position Separated Scatter Plots
  const int nGroups = 4;
  enum GroupsEnum {Barrel, Top, Bottom, HK};
  
  TString GroupSelection[nGroups] = {
    "group <= 17",
    "18 <= group && group <= 25",
    "group >= 26",
    ""
  };

  TString GroupName[nGroups] = {
    "Barrel",
    "Top",
    "Bottom",
    "HK"
  };
  
  int GroupColors[nGroups] = {kBlack, kRed, kGreen, kGray+1};

  for (int ifile=0; ifile<nFiles; ifile++) {    
    for (int jfile=ifile+1; jfile<nFiles; jfile++) {

      TCanvas *c_pc2pe_2d_sep = new TCanvas(1);
      TLegend *leg_pc2pe_2d_sep = new TLegend(0.2, 0.65, 0.45, 0.88);
      leg_pc2pe_2d_sep->SetHeader("PMT Type (Correlation)");

      int iuseFile = ifile;
      if (ifile==nFiles-1) iuseFile = 0;
      int juseFile = jfile;
      if (jfile==nFiles-1) juseFile = 0;

      for (int igroup=0; igroup<nGroups; igroup++) {

	TString cut_ihkpmt, cut_jhkpmt;
	if (igroup == HK) {
	  cut_ihkpmt = "pmtflag"+TreeVarNames[iuseFile]+"==6";
	  cut_jhkpmt = "pmtflag"+TreeVarNames[juseFile]+"==6";
	} else {
	  cut_ihkpmt = "pmtflag"+TreeVarNames[iuseFile]+"!=6";
	  cut_jhkpmt = "pmtflag"+TreeVarNames[juseFile]+"!=6";
	}
	
	TString cut_all = cut_ihkpmt + "&&" + cut_jhkpmt;
	if (igroup != HK) cut_all += "&&" + GroupSelection[igroup];

	TString cut_badpmt = "&& ";
	if (ifile==nFiles-1) cut_badpmt+ = "qe_bad_sk4 == 0";
	else cut_badpmt += "pc2pe_bad"+TreeVarNames[ifile]+" == 0";
	if (jfile==nFiles-1) cut_badpmt += "&& qe_bad_sk4 == 0";
	else cut_badpmt += " && pc2pe_bad"+TreeVarNames[ifile]+" == 0";
	cut_all += cut_badpmt;
	
	TString histname = "pc2pe_2d_sep"+TreeVarNames[ifile]+"_vs"+TreeVarNames[jfile]+"_"+GroupName[igroup];
    
	TString plotVarX = "rationorm"+TreeVarNames[ifile];
	if (ifile==nFiles-1) plotVarX = "qe_sk4";
	TString plotVarY = "rationorm"+TreeVarNames[jfile];
	if (jfile==nFiles-1) plotVarY = "qe_sk4";
	
	pc2pe->Draw(plotVarY+":"+plotVarX, cut_all, "goff");
	
	//h_pc2pe_2d[ifile]->SetTitle(Form("Correlation = %.2f", h_pc2pe_2d[ifile]->GetCorrelationFactor()));

	TGraph *gr = new TGraph(pc2pe->GetSelectedRows(), pc2pe->GetV2(), pc2pe->GetV1());

	if (!gr->GetN()) continue;
	
	gr->SetTitle(";pc2pe "+FileTitles[ifile]+";pc2pe "+FileTitles[jfile]);
	gr->SetMarkerSize(0.06);
	gr->SetMarkerColor(GroupColors[igroup]);
	gr->SetFillColor(GroupColors[igroup]);
	
	if (!igroup) gr->Draw("AP");
	else gr->Draw("same P");

	leg_pc2pe_2d_sep->AddEntry(gr, GroupName[igroup]+Form(" (%.2f)",gr->GetCorrelationFactor()), "f");
	
      }

      leg_pc2pe_2d_sep->Draw();
      c_pc2pe_2d_sep->Print("figures/pc2pe_2d_sep"+TreeVarNames[ifile]+"_vs"+TreeVarNames[jfile]+".png");


    }
  }
}
