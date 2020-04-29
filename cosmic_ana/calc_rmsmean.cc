{

  bool doReal = 1;

  int iPlot = 2;
  TString plotTitle = "RMS/Mean of ";
  if (iPlot==1) plotTitle = "#frac{Data-MC}{MC} for ";
  else if (iPlot==2) plotTitle = "Slope for ";
  else if (iPlot==3) plotTitle = "";

  const int nPars = 5;
  const TString parNames[nPars] = {"gausmean", "sigma", "mean", "rms", "mispid"};
  const TString parTitles[nPars] = {"Gaus. Mean of ", "#sigma of ", "Mean of ", "RMS of ", "Mis-PID Rate (%)"};
  enum parEnum {igausmean, isigma, imean, irms, imispid};

  int ipar = 2;
  TString yAxisTitle;
  //if (ipar==2) yAxisTitle = ">";
  
  if (iPlot==1) yAxisTitle += " (%)";
  else if (iPlot==2) yAxisTitle += " (%/year) ";
  else if (iPlot==3) yAxisTitle += " (MeV/c/cm)";
  


  const int nRanges = 5; //7
  const double rangeLow = 0;
  const double rangeHigh = 25; //35 // cm
  const double rangeStep = 5;  // (rangeHigh-rangeLow)/nRanges
 
  TFile *infile = new TFile("pars_vs_mom.root");

  TH1D *h_dist[nRanges+1];
  
  TH1D *h_rms_vs_range[2], *h_rms_vs_range_mc[2];
  for (int ireco=0; ireco<2; ireco++) {
    h_rms_vs_range[ireco] = new TH1D(Form("h_rms_vs_range_reco%d",ireco),";Range (m); "+plotTitle+parTitles[ipar]+"Mom."+yAxisTitle, nRanges+1, rangeLow-rangeStep*1,rangeHigh);
    h_rms_vs_range_mc[ireco] = new TH1D(Form("hmc_rms_vs_range_reco%d",ireco),";Range (m); "+plotTitle+parTitles[ipar]+"Mom."+yAxisTitle, nRanges+1, rangeLow-rangeStep*1,rangeHigh);

    h_rms_vs_range[ireco]->SetLineWidth(2);
    h_rms_vs_range[ireco]->SetLineColor(ireco+1);
    h_rms_vs_range[ireco]->SetMarkerColor(ireco+1);

    h_rms_vs_range_mc[ireco]->SetLineWidth(2);
    h_rms_vs_range_mc[ireco]->SetLineColor(kRed);
    h_rms_vs_range_mc[ireco]->SetMarkerColor(kRed);
  }

  TString histName;

  int icut=1;
  
  for (int ireco=0; ireco<2; ireco++) {
    for (int idist=0; idist<nRanges+1; idist++) {
    
      if (idist==0) 
	histName = Form("h_cosmic_emom_vs_year_reco%d_%s_vs_year",ireco, parNames[ipar].Data());
      else {
	double rangeLowEdge = (rangeLow+(idist-1)*rangeStep);
	double rangeHighEdge = (rangeLow+(idist)*rangeStep);
	histName = Form("h_cosmic_mom_over_range_vs_year_reco%d_rangesep%03d_%03d_icut%d_%s_vs_year",ireco,(int)rangeLowEdge,(int)rangeHighEdge, icut, parNames[ipar].Data());
      }
      
      h_dist[idist] = (TH1D*)infile->Get(histName);

      if (iPlot==0) {
	
	TH1D *h_rms = new TH1D(histName+"_rmsmean",histName+"_rmsmean",100,h_dist[idist]->GetMinimum(),h_dist[idist]->GetMaximum());
	h_rms->Sumw2();
	
	int nBins = h_dist[idist]->GetNbinsX();
	for (int ibin=1; ibin<=nBins; ibin++) {
	  
	  double val = h_dist[idist]->GetBinContent(ibin);
	  if ( val<=0.001 || val!=val) continue;
	  
	  h_rms->Fill(val);
	  
	}
      
	double rmsmean = 100*h_rms->GetRMS()/h_rms->GetMean();
      
	//cout << histName << " " << h_rms->GetRMS() << " " << h_rms->GetMean() << " " <<  rmsmean << endl;
      
	if (!doReal && ireco==0) {
	  double rmsmean_apfit = 0;
	  if (idist==0) rmsmean_apfit = 0.41;
	  else if (idist>1) rmsmean_apfit = 0.39;
	  rmsmean = rmsmean_apfit;
	}
      
	h_rms_vs_range[ireco]->SetBinContent(idist+1, rmsmean);
	
      } else if (iPlot>=1) {
	TString f_name = Form("f_pol1_%s",histName.Data());

	TF1 *f_pol1 = (TF1*)infile->Get(f_name);
	
	double p0 = f_pol1->GetParameter(0);
	double p1 = f_pol1->GetParameter(1);
	double p0err = f_pol1->GetParError(0);
	double p1err = f_pol1->GetParError(1);
	
	//int apr2009bin = h_dist[idist]->FindBin(2009+3./12);
	//p0 = h_dist[idist]->GetBinContent(apr2009bin);
	//p0err = h_dist[idist]->GetBinError(apr2009bin);

	//TLine *l_horiz = (TLine*)infile->Get("l_"+histName.ReplaceAll("h_","hmc_"));
	TLine *l_horiz = (TLine*)infile->Get("l_"+histName);
	double MCparVal = l_horiz->GetY1();

	h_rms_vs_range_mc[ireco]->SetBinContent(idist+1, MCparVal);

	if (iPlot==1) {
	  double datamc = p0/MCparVal;
	  double datamcerr = 100*datamc*p0err/p0;
	  datamc = 100*(datamc-1);
	  
	  // From: hignight_20130510_escale_update.pdf
	  if (!doReal && ireco==1) {
	    if (idist==0) { datamc = 0.7; datamcerr = 0.2; }
	    if (idist==1) { datamc = 0; datamcerr = 0; }
	    if (idist==2) { datamc = -1.3; datamcerr = 0.3; }
	    if (idist==3) { datamc = -1.1; datamcerr = 0.2; }
	    if (idist==4) { datamc = -1.3; datamcerr = 0.2; }
	    if (idist==5) { datamc = -0.6; datamcerr = 0.2; }
	    if (idist==6) { datamc = 0.2; datamcerr = 0.3; }
	    if (idist==7) { datamc = -0.8; datamcerr = 0.3; }
	  }

	  h_rms_vs_range[ireco]->SetBinContent(idist+1, datamc);
	  h_rms_vs_range[ireco]->SetBinError(idist+1, datamcerr);
	} else if (iPlot==2) {
	  double slope = 100*p1/p0;
	  double slopeerr = 100*p1err/p0;
	  
	  h_rms_vs_range[ireco]->SetBinContent(idist+1, slope);
	  h_rms_vs_range[ireco]->SetBinError(idist+1, slopeerr);
	} else if (iPlot==3) {
	  double datamc = p0;
	  double datamcerr = p0err;

	  h_rms_vs_range[ireco]->SetBinContent(idist+1, datamc);
	  h_rms_vs_range[ireco]->SetBinError(idist+1, datamcerr);
	}
	
	

      }
      /**/
    }
  }
  
  double ymin, ymax;
  if (ipar==imean || ipar==igausmean) {
    if (iPlot==0) h_rms_vs_range[1]->GetYaxis()->SetRangeUser(0,1);
    else if (iPlot==1) {
      ymin = -6; ymax = 3;
      h_rms_vs_range[1]->GetYaxis()->SetRangeUser(ymin,ymax);
    }
    else if (iPlot==2) {
      ymin = -0.53; ymax = 0.53;
      h_rms_vs_range[1]->GetYaxis()->SetRangeUser(ymin,ymax);
    }
    else if (iPlot==3) h_rms_vs_range[1]->GetYaxis()->SetRangeUser(2,2.6);
  } else  if (ipar==isigma || ipar==irms) {
    if (iPlot==0) h_rms_vs_range[1]->GetYaxis()->SetRangeUser(0,1);
    else if (iPlot==1) h_rms_vs_range[1]->GetYaxis()->SetRangeUser(-15,30);
    else if (iPlot==2) h_rms_vs_range[1]->GetYaxis()->SetRangeUser(-1.5,1.5);
  }


  h_rms_vs_range[1]->Draw("");
  h_rms_vs_range[1]->Draw("hist same");
  //h_rms_vs_range[0]->Draw("same");

  if (iPlot==3) {
    h_rms_vs_range[1]->GetXaxis()->SetRangeUser(0,25);
    h_rms_vs_range[1]->SetLineColor(kBlue);
    h_rms_vs_range[1]->SetMarkerColor(kBlue);
    
    h_rms_vs_range_mc[1]->Draw("SAME");
    h_rms_vs_range[1]->Draw("SAME");

    TLegend *leg = new TLegend(0.6, 0.6, 0.80, 0.85);
    leg->SetFillStyle(0);
    leg->SetTextFont(132);    
    leg->AddEntry(h_rms_vs_range_mc[1], "MC", "l");
    leg->AddEntry(h_rms_vs_range[1], "Data", "p");
    leg->Draw("SAME");
  }

  // Decay-e label
  TText *t_decaye = new TText(-2, ymin+(ymax-ymin)*0.08, "Decay-e");
  t_decaye->SetTextAngle(90);
  t_decaye->SetTextSize(0.06);
  t_decaye->Draw("SAME");
  
  TLine *l_decaye;
  l_decaye = new TLine(0,ymin,0,ymax);
  l_decaye->SetLineWidth(2);
  l_decaye->SetLineColor(kBlack);
  l_decaye->Draw("SAME");

  // Muon label
  TArrow *arw = new TArrow();
  arw->SetLineWidth(2);
  arw->DrawArrow(0, ymin+(ymax-ymin)*0.15, 6,  ymin+(ymax-ymin)*0.15, 0.05,">");

  TText *t_muon = new TText(6.3, ymin+(ymax-ymin)*0.13, "Muons");
  t_muon->SetTextSize(0.06);
  t_muon->Draw("SAME");
  
  if (iPlot==0) {
    //for (int ireco=1; ireco<2; ireco++) {
    for (int ireco=0; ireco<2; ireco++) {
      int nBins = h_rms_vs_range[ireco]->GetNbinsX();
      for (int ibin=1; ibin<=nBins; ibin++) {
	
	if (!doReal && ireco==0 && (ibin!=1 && ibin!=5)) continue;
      
	double x = h_rms_vs_range[ireco]->GetBinLowEdge(ibin);
	double y = h_rms_vs_range[ireco]->GetBinContent(ibin);
      
	cout << ireco << " " << x << " " << y << endl;

	TText *t_val = new TText(x+0.8,y+(ireco*2-1)*0.05-0.02, Form("%.02f", y));
	t_val->SetTextSize(0.05);
	t_val->SetTextColor(ireco+1);
	t_val->Draw();

      }
    }  

    //c1->Print(Form("images/rmsmean_real_%s_icut%d.png",parNames[ipar].Data(),icut);
    
  } else if (iPlot>=1) {

    TLine *l_zero;
    l_zero = new TLine(-5,0,rangeHigh,0);
    l_zero->SetLineWidth(2);
    l_zero->SetLineColor(kGray);

    l_zero->Draw("SAME");

    // Energy scale error
    if (iPlot>=1) {
      double escaleerr = 2.4 / (1+(iPlot-1)*4);
      TLine *l_escaleerr;
      l_escaleerr = new TLine(-5,-escaleerr,rangeHigh,-escaleerr);
      l_escaleerr->SetLineWidth(2);
      l_escaleerr->SetLineStyle(2);
      l_escaleerr->SetLineColor(kGreen-2);
      l_escaleerr->Draw("SAME");
      l_escaleerr = new TLine(-5,escaleerr,rangeHigh,escaleerr);
      l_escaleerr->SetLineWidth(2);
      l_escaleerr->SetLineStyle(2);
      l_escaleerr->SetLineColor(kGreen-2);
      l_escaleerr->Draw("SAME");

      //TLatex t_escaleerr(rangeHigh+2.5, 0, "#deltaE_{e-scale} ~ #pm 2.4%");
      //t_escaleerr->SetTextAngle(90);
      //t_escaleerr->SetTextAlign(22);
      //t_escaleerr->SetTextSize(0.06);
      //t_escaleerr->SetTextColor(kGreen-2);
      //t_escaleerr->Draw("SAME");
      TLegend *leg = new TLegend(0.4, 0.7-(0.07*(iPlot-1)), 0.8, 0.84);
      leg->SetFillStyle(0);
      leg->AddEntry(h_rms_vs_range[1], "fiTQun");
      if (iPlot==1) leg->AddEntry(l_escaleerr, "E-scale error (2.4%)","l");
      else if (iPlot==2) {
	leg->SetMargin(.3);
	leg->AddEntry(l_escaleerr, "#frac{E-scale error (2.4%)}{5 years}","l");
      }
      leg->Draw();

    }

    //h_rms_vs_range[1]->Draw("hist same");
    //h_rms_vs_range[0]->Draw("hist same");
    if (iPlot==1) c1->Print(Form("images/momrange_datamc_real_%s_icut%d.pdf",parNames[ipar].Data(),icut));
    else if (iPlot==2) c1->Print(Form("images/momrange_slope_real_%s_icut%d.pdf",parNames[ipar].Data(),icut));
    else if (iPlot==3) c1->Print(Form("images/momrange_absval_real_%s_icut%d.pdf",parNames[ipar].Data(),icut));
  }

}
