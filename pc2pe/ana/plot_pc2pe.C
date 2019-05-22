{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.15);
  gErrorIgnoreLevel = kWarning;  // Suppress "Info in <TCanvas::Print>" messages

  
  const int nFiles = 9;

  enum file_enum {sk4, sk5, sk5i, sk5n, sk5avg, sk4official, sk3_qe, sk4_qe, sk5_qe};
  
  TString TreeVarNames[nFiles] = {
    "pc2pe",
    "pc2pe",
    "pc2pe",
    "pc2pe",
    "pc2pe",
    "pc2pe",
    "qe",
    "qe",
    "qe"
  };

  TString FileDatasets[nFiles] = {
    "_sk4",
    "_sk5",
    "_sk5i",
    "_sk5n",
    "_sk5avg",
    "_sk4official",
    "_sk3", // QE
    "_sk4", // QE    
    "_sk5" // QE
  }
  
  int Colors[nFiles] = {
    kGray+2,
    kRed,
    kBlue,
    kGreen-2,
    kBlack,
    kGray+1,
    kMagenta-7,
    kMagenta-9,
    kMagenta
  };

  TString FileTitles[nFiles] = {
    "SK4 Reana.",
    "SK5 Mar.",
    "SK5 Inv.",
    "SK5 Apr.",
    "SK5 Avg.",
    "SK3/4 off.",
    "SK3",
    "SK4",    
    "SK5"
  };

  TString CanvasAppend = "_all";

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

  
  int Colors[nFiles] = {
    kBlack,
    kRed,
    kBlue
  };

  TString FileTitles[nFiles] = {
    "March",
    "Upside Down",
    "April"
  };

  TString CanvasAppend = "_sk5only";

  /**/
  
  /*
  const int nFiles = 2;

  TString FileDatasets[nFiles] = {
    "_sk5avg",
    "_sk4official"
  };

  TString TreeVarNames[nFiles] = {
    "pc2pe",
    "pc2pe"
  };
  
  int Colors[nFiles] = {
    kBlack,
    kGray+2
  };

  TString FileTitles[nFiles] = {
    "SK5",
    "SK3/4 official"
  };
  
  TString CanvasAppend = "";
  /**/
  
  TFile *infile = new TFile("pc2pe_output.root");

  TString DrawOpts = "HIST";

  // PMT Type separated
  const int nPMTtypes = 4;
  int PMTflags[nPMTtypes] = {3, 4, 6, -1};
  TString PMTTypeNames[nPMTtypes] = {"SK2", "SK3", "HK", "Other"};
  int PMTtypeColors[nPMTtypes] = {kBlue, kRed, kBlack, kGray};

  // 1D distributions

  // All

  TCanvas *c_pc2pe = new TCanvas(1);
  c_pc2pe->SetLogy(1);
  TLegend *leg = new TLegend(0.65, 0.7, 1, 1);
  leg->SetHeader("Dataset (RMS)");

  TH1D *h_pc2pe[nFiles];

  for (int ifile=0; ifile<nFiles; ifile++) {

    TString histname = TreeVarNames[ifile]+FileDatasets[ifile];
    TString histtitle = ";"+TreeVarNames[ifile]; //+" "+FileTitles[ifile];
    histtitle += ";Number of Channels";
    
    h_pc2pe[ifile] = new TH1D(histname, histtitle, 50, 0.1, 2);

    // Select good PMT flags
    TString cut_pmttype = "(";
    for (int jpmttype=0; jpmttype<nPMTtypes-1; jpmttype++) {
      cut_pmttype += "pmtflag"+FileDatasets[ifile]+Form("==%d",PMTflags[jpmttype]);
      if (jpmttype<nPMTtypes-2) cut_pmttype += " || ";
    }
    cut_pmttype += ")";

    cut_pmttype += " && "+TreeVarNames[ifile]+"_bad"+FileDatasets[ifile]+" == 0";
    
    TString plotVar = TreeVarNames[ifile]+FileDatasets[ifile];
    
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
  c_pc2pe->Print("figures/pc2pe"+CanvasAppend+".png");

  
  TH1D *h_pc2pe_pmttypes[nFiles][nPMTtypes];
  for (int ifile=0; ifile<nFiles; ifile++) {

    TCanvas *c_pc2pe_pmttype = new TCanvas(1);
    c_pc2pe_pmttype->SetLogy(1);
    TLegend *leg_pc2pe_pmttype = new TLegend(0.6, 0.7, 1, 1);
    leg_pc2pe_pmttype->SetHeader("PMT Type (Channels, Mean, RMS)");
    
    h_pc2pe[ifile]->SetLineColor(kGray+2);
    h_pc2pe[ifile]->Draw();
    h_pc2pe[ifile]->SetTitle(FileTitles[ifile]);

    cout << "1D " << FileTitles[ifile] << " " << TreeVarNames[ifile] << endl;
    cout << "     Mean,  RMS" << endl;
    
    leg_pc2pe_pmttype->AddEntry(h_pc2pe[ifile], Form("All (%d, %.2f, %.3f)", (int)h_pc2pe[ifile]->GetEntries(), h_pc2pe[ifile]->GetMean(), h_pc2pe[ifile]->GetRMS()), "l");

    cout << "All: " <<  Form("%.2f ± %.2f, ", h_pc2pe[ifile]->GetMean(), h_pc2pe[ifile]->GetMeanError());
    cout << Form("%.2f ± %.2f", 100*h_pc2pe[ifile]->GetRMS(), 100*h_pc2pe[ifile]->GetRMSError()) << endl;
    
    for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++) {
      
      TString histname = TreeVarNames[ifile]+FileDatasets[ifile]+"_pmt"+PMTTypeNames[ipmttype];
      TString histtitle = ";"+TreeVarNames[ifile];
      histtitle += ";Number of Channels";
      
      h_pc2pe_pmttypes[ifile][ipmttype] = new TH1D(histname, histtitle, 50, 0.1, 2);

      // Use SK4 PMT type list
      TString cut_pmttype = "pmtflag"+FileDatasets[ifile]+Form("==%d",PMTflags[ipmttype]);
      if (ipmttype==nPMTtypes-1) {
	cut_pmttype = "";
	for (int jpmttype=0; jpmttype<nPMTtypes-1; jpmttype++) {
	  cut_pmttype += "pmtflag"+FileDatasets[ifile]+Form("!=%d",PMTflags[jpmttype]);
	  if (jpmttype<nPMTtypes-2) cut_pmttype += " && ";
	}
      }
      cut_pmttype += " && "+TreeVarNames[ifile]+"_bad"+FileDatasets[ifile]+" == 0";
      
      TString plotVar = TreeVarNames[ifile]+FileDatasets[ifile];

      pc2pe->Project(histname, plotVar, cut_pmttype);
      
      int nentries = h_pc2pe_pmttypes[ifile][ipmttype]->GetEntries();
      if (!nentries) continue;
      
      h_pc2pe_pmttypes[ifile][ipmttype]->SetLineColor(PMTtypeColors[ipmttype]);
      h_pc2pe_pmttypes[ifile][ipmttype]->SetLineWidth(2);

      if (ipmttype>=nPMTtypes-1) continue;
      
      //if (!ipmttype) 
      //h_pc2pe_pmttypes[ifile][ipmttype]->Draw(); //DrawOpts += "same";
      //else
      h_pc2pe_pmttypes[ifile][ipmttype]->Draw("same");
    
      leg_pc2pe_pmttype->AddEntry(h_pc2pe_pmttypes[ifile][ipmttype], PMTTypeNames[ipmttype]+Form(" (%d, %.2f, %.3f)", nentries, h_pc2pe_pmttypes[ifile][ipmttype]->GetMean(), h_pc2pe_pmttypes[ifile][ipmttype]->GetRMS()), "l");

      cout << PMTTypeNames[ipmttype] << ": " << Form("%.2f ± %.2f, ", h_pc2pe_pmttypes[ifile][ipmttype]->GetMean(), h_pc2pe_pmttypes[ifile][ipmttype]->GetMeanError());
      cout <<  Form("%.2f ± %.2f", 100*h_pc2pe_pmttypes[ifile][ipmttype]->GetRMS(), 100*h_pc2pe_pmttypes[ifile][ipmttype]->GetRMSError()) << endl;

    }
    cout << endl;
    
    leg_pc2pe_pmttype->Draw();
    c_pc2pe_pmttype->Print("figures/"+TreeVarNames[ifile]+FileDatasets[ifile]+"_pmttype.png");
    
  }

  // 2D Distributions
  TH2D *h_pc2pe_2d[nFiles];

  TString var_title[2];
  for (int ifile=0; ifile<nFiles; ifile++) {

    var_title[0] = FileTitles[ifile]+" "+TreeVarNames[ifile];
    TString i_var = TreeVarNames[ifile]+FileDatasets[ifile];
      
    for (int jfile=ifile+1; jfile<nFiles; jfile++) {

      var_title[1] = FileTitles[jfile]+" "+TreeVarNames[jfile];
      TString j_var = TreeVarNames[jfile]+FileDatasets[jfile];
      
      cout << var_title[1].Data() << " vs " << var_title[0].Data() << endl;

      TCanvas *c_pc2pe_2d = new TCanvas(1);
      
      TString histname = j_var+"_vs_"+i_var;

      h_pc2pe_2d[ifile] = new TH2D(histname, ";"+var_title[0]+";"+var_title[1]+";Number of Channels", 120, 0.1, 2, 120, 0.1, 2);

      TString cut_badpmt = TreeVarNames[ifile]+"_bad"+FileDatasets[ifile]+" == 0";
      cut_badpmt += " && "+ TreeVarNames[jfile]+"_bad"+FileDatasets[jfile]+" == 0";
      
      pc2pe->Project(histname, j_var+":"+i_var, cut_badpmt);

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
	
	tRMSoverMeans->AddText(var_title[iaxis]+Form(": %.2f #pm %.2f", RMSoverMean, RMSoverMeanErr));

	cout << Form("%.2f ± %.2f", RMSoverMean, RMSoverMeanErr);
	if (!iaxis) {
	  cout  << ", ";
	}
      }
      cout << endl << endl;
      
      tRMSoverMeans->Draw();
      
      c_pc2pe_2d->Print("figures/"+histname+".png");
    }
  }
  /**/

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

    var_title[0] = FileTitles[ifile]+" "+TreeVarNames[ifile];
    TString i_var = TreeVarNames[ifile]+FileDatasets[ifile];

    for (int jfile=ifile+1; jfile<nFiles; jfile++) {

      var_title[1] = FileTitles[jfile]+" "+TreeVarNames[jfile];
      TString j_var = TreeVarNames[jfile]+FileDatasets[jfile];
      
      TCanvas *c_pc2pe_2d_sep = new TCanvas(1);
      TLegend *leg_pc2pe_2d_sep = new TLegend(0.2, 0.65, 0.45, 0.88);
      leg_pc2pe_2d_sep->SetHeader("PMT Type (Corr., x #sigma/#mu (%), y #sigma/#mu (%))");

      TString histname_base = j_var+"_vs_"+i_var;

      for (int igroup=0; igroup<nGroups; igroup++) {

	TString cut_ihkpmt, cut_jhkpmt;
	if (igroup == HK) {
	  cut_ihkpmt = "pmtflag"+FileDatasets[ifile]+"==6";
	  cut_jhkpmt = "pmtflag"+FileDatasets[jfile]+"==6";
	} else {
	  cut_ihkpmt = "pmtflag"+FileDatasets[ifile]+"!=6";
	  cut_jhkpmt = "pmtflag"+FileDatasets[jfile]+"!=6";
	}
	
	TString cut_all = cut_ihkpmt + " && " + cut_jhkpmt;
	if (igroup != HK) cut_all += " && " + GroupSelection[igroup];

	TString cut_badpmt = TreeVarNames[ifile]+"_bad"+FileDatasets[ifile]+" == 0";
	cut_badpmt += " && "+ TreeVarNames[jfile]+"_bad"+FileDatasets[jfile]+" == 0";
	cut_all += " && " + cut_badpmt;
	
	TString histname = histname_base+"_"+GroupName[igroup];
    	
	pc2pe->Draw(j_var+":"+i_var, cut_all, "goff");
	
	TGraph *gr = new TGraph(pc2pe->GetSelectedRows(), pc2pe->GetV2(), pc2pe->GetV1());

	if (!gr->GetN()) continue;
	
	gr->SetTitle(";"+var_title[0]+";"+var_title[1]);
	gr->SetMarkerSize(0.06);
	gr->SetMarkerColor(GroupColors[igroup]);
	gr->SetFillColor(GroupColors[igroup]);
	
	if (!igroup) gr->Draw("AP");
	else gr->Draw("same P");

	leg_pc2pe_2d_sep->AddEntry(gr, GroupName[igroup]+Form(" (%.2f, %.2f, %.2f)",
							      gr->GetCorrelationFactor(),
							      100*gr->GetRMS(1)/gr->GetMean(1),
							      100*gr->GetRMS(2)/gr->GetMean(2)
							      ), "f");
	
      }

      leg_pc2pe_2d_sep->Draw();
      c_pc2pe_2d_sep->Print("figures/"+histname_base+"_groupsep.png");

    }
  }
  /**/

  // PMT type separated 2D dists
  for (int ifile=0; ifile<nFiles; ifile++) {    

    var_title[0] = FileTitles[ifile]+" "+TreeVarNames[ifile];
    TString i_var = TreeVarNames[ifile]+FileDatasets[ifile];

    for (int jfile=ifile+1; jfile<nFiles; jfile++) {

      var_title[1] = FileTitles[jfile]+" "+TreeVarNames[jfile];
      TString j_var = TreeVarNames[jfile]+FileDatasets[jfile];
      
      TCanvas *c_pc2pe_2d_pmtsep = new TCanvas(1);
      TLegend *leg_pc2pe_2d_pmtsep = new TLegend(0.6, 0.7, 0.99, 0.99);
      leg_pc2pe_2d_pmtsep->SetHeader("PMT Type (Corr., x #sigma/#mu (%), y #sigma/#mu (%))");

      cout << var_title[1].Data() << " vs " << var_title[0].Data() << endl;

      TString histname_base = j_var+"_vs_"+i_var;

      for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++) {
      
	TString histname = histname_base+"_pmtsep_"+PMTTypeNames[ipmttype];
	TString histtitle = ";"+var_title[0]+";"+var_title[1];
	
	TH2D *h_pc2pe_pmttypes_tmp = new TH2D(histname, histtitle, 150, 0.7, 1.6, 150, 0.7, 1.6);

	TString cut_pmttype = "pmtflag"+FileDatasets[ifile]+Form("==%d",PMTflags[ipmttype]);
	cut_pmttype += " && pmtflag"+FileDatasets[jfile]+Form("==%d ",PMTflags[ipmttype]);
	
	TString cut_all = cut_pmttype;
	
	TString cut_badpmt = TreeVarNames[ifile]+"_bad"+FileDatasets[ifile]+" == 0";
	cut_badpmt += " && "+ TreeVarNames[jfile]+"_bad"+FileDatasets[jfile]+" == 0";
	cut_all += " && " + cut_badpmt;
	
	pc2pe->Project(histname, j_var+":"+i_var, cut_all);
	pc2pe->Draw(j_var+":"+i_var, cut_all, "goff");
	
	TGraph *gr = new TGraph(pc2pe->GetSelectedRows(), pc2pe->GetV2(), pc2pe->GetV1());

	if (!gr->GetN()) continue;
	
	gr->SetTitle(histtitle);
	gr->SetMarkerSize(0.06);
	gr->SetMarkerColor(PMTtypeColors[ipmttype]);
	gr->SetFillColor(PMTtypeColors[ipmttype]);
	
	if (!ipmttype) gr->Draw("AP");
	else gr->Draw("same P");

	TString RMSoverMeans;

	cout << PMTTypeNames[ipmttype] << ": ";

	for (int iaxis=0; iaxis<=1; iaxis++) {
	  double RMS2d, Mean2d, RMSErr2d, MeanErr2d;
	  Mean2d = h_pc2pe_pmttypes_tmp->GetMean(iaxis+1);
	  MeanErr2d = h_pc2pe_pmttypes_tmp->GetMeanError(iaxis+1);
	  RMS2d = h_pc2pe_pmttypes_tmp->GetRMS(iaxis+1);
	  RMSErr2d = h_pc2pe_pmttypes_tmp->GetRMSError(iaxis+1);
	  
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
	
	leg_pc2pe_2d_pmtsep->AddEntry(gr, PMTTypeNames[ipmttype]+Form(" (%.2f, %s)",
								      gr->GetCorrelationFactor(),
								      RMSoverMeans.Data(), "f"));
	
      }
      cout << endl;
      
      leg_pc2pe_2d_pmtsep->Draw();
      c_pc2pe_2d_pmtsep->Print("figures/"+histname_base+"_pmtsep.png");

    }
  }
}
