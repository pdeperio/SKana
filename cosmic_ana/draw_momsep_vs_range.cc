{
  gStyle->SetPadRightMargin(0.08);

  gStyle->SetOptFit(0);
  gStyle->SetStatX(.95);
  gStyle->SetStatY(1);
  gStyle->SetStatW(0.4);
  gStyle->SetStatH(0.2);

  const TString pidTitle[2] = {"1st Sub-Event (#mu)", "2nd Sub-Event (Decay-e)"};

  const int nRanges = 7; //5; //7;
  const double rangeLow = 0;
  const double rangeHigh = 35; //25; //35; // cm
  const double rangeStep = 5;  // (rangeHigh-rangeLow)/nRanges
 
 
  const int nPars = 2;
  const TString parNames[nPars] = {"mean", "sigma"};
  const TString parTitles[nPars] = {"Mean", "#sigma"};
  enum parEnum {imean, isigma};

  TFile *infile = new TFile("pars_vs_mom.root");

  TGraphErrors *g_mean_vs_momsep = 0;
  TGraphErrors *g_rms_vs_momsep = 0;
  int ipoint = 0;

  int icut=3;
  
  for (int ireco=0; ireco<2; ireco++) {

    for (int ipar=0; ipar<nPars; ipar++) {
            
      TLegend *leg = new TLegend(0.3, 0.65, 0.53, 0.88);
      leg->SetFillStyle(0);
      
      TCanvas *c_par = new TCanvas(Form("c_par%d_ireco%d",ipar,ireco),Form("c_par%d_ireco%d",ipar,ireco),0,0,1800,1000);
      c_par->Divide((nRanges+1)/2,2);
      //TCanvas *c_par = new TCanvas(Form("c_momrange_par%d_ireco%d",ipar,ireco),Form("c_momrange_par%d_ireco%d",ipar,ireco),0,0,700,1000);
      //c_par->Divide(2,(nRanges+1)/2);
	
      int idraw = 2;

      for (int imom=0; imom<nRanges; imom++) {
      
	double rangeLowEdge = (rangeLow+imom*rangeStep);
	double rangeHighEdge = (rangeLow+(imom+1)*rangeStep);
	TString h_name = Form("h_cosmic_mom_over_range_vs_year_reco%d_rangesep%03d_%03d_icut%d",ireco,(int)rangeLowEdge,(int)rangeHighEdge,icut);

	TString h_name_par = Form(Form("%s_%s_vs_year",h_name.Data(),parNames[ipar].Data()));

	cout << h_name.Data() << endl;
	TH1D *h_tmp = (TH1D*)infile->Get(h_name_par.Data());
	h_tmp->SetTitle(Form("%d to %d m", (int)rangeLowEdge, (int)rangeHighEdge)); 
	  
	h_tmp->SetLineColor(kBlue);
	h_tmp->SetMarkerColor(kBlue);
	
	double yrange[2];

	//if (ipar==0) {
	//  if (imom) {
	//    yrange[0] = 2.1;
	//    yrange[1] = 2.4;
	//  }
	//  else {
	//    yrange[0] = 2.4;
	//    yrange[1] = 2.8;
	//  }
	//  h_tmp->GetYaxis()->SetRangeUser(yrange[0],yrange[1]);
	//
	//
	//}
	  
	c_par->cd(idraw);
	h_tmp->Draw();

	TString f_name = Form("f_pol1_%s",h_name_par.Data());
	TF1 *f_pol1 = (TF1*)infile->Get(f_name);
	f_pol1->Draw("SAME");
	  
	TLine *l_horiz = (TLine*)infile->Get("l_"+h_name_par);
	l_horiz->Draw("SAME");
	  
	TLine *l_apr_2009 = (TLine*)infile->Get("l_apr_2009_"+h_name_par);
	//l_apr_2009->SetY1(yrange[0]);
	//l_apr_2009->SetY2(yrange[1]);
	l_apr_2009->Draw("SAME");
	  
	h_tmp->Draw("same");

	double MCparVal = l_horiz->GetY1();
      
	double f_pol1_p[2];
	double f_pol1_p_err[2];
	for (int ipar2=0; ipar2<2; ipar2++) {
	  f_pol1_p[ipar2] = f_pol1->GetParameter(ipar2);
	  f_pol1_p_err[ipar2] = f_pol1->GetParError(ipar2);
	}

	//double datamc = (f_pol1_p[0]/MCparVal);
	//double datamcerr = 100*datamc*sqrt(pow(f_pol1_p_err[0]/f_pol1_p[0],2));
	double dataval = f_pol1_p[0];
	double dataerr = f_pol1_p_err[0];
	//int apr2009bin = h_tmp->FindBin(2009+3./12);
	//double dataval = h_tmp->GetBinContent(apr2009bin);
	//double dataerr = h_tmp->GetBinError(apr2009bin);
	double datamc = (dataval/MCparVal);
	double datamcerr = 100*datamc*sqrt(pow(dataerr/dataval,2));

	datamc = 100*(datamc-1);
      
	//TLatex *t_datamc = new TLatex(0.34,0.28,Form("Data/MC: %.02f%%",datamc));
	TLatex *t_datamc = new TLatex(0.34,0.28,Form("Data/MC: %.02f #pm %.02f%%",datamc,datamcerr));
	t_datamc->SetTextSize(0.06);
	t_datamc->SetTextFont(132);
	t_datamc->SetTextColor(kRed);
	t_datamc->SetNDC(true);
	t_datamc->Draw("SAME");
	      	      
	TLatex *t_slope = new TLatex(0.34,0.2,Form("%.02f #pm %.02f%% per year",100*f_pol1_p[1]/f_pol1_p[0], 100*f_pol1_p_err[1]/f_pol1_p[0]));
	t_slope->SetTextSize(0.06);
	t_slope->SetTextFont(132);
	t_slope->SetTextColor(kOrange+2);
	t_slope->SetNDC(true);
	t_slope->Draw("SAME");
	  

	if (ipar==0) {
	    
	  TString h_rms = Form("%s_%s_vs_year",h_name.Data(),parNames[isigma].Data());
	
	  TLine *l_rms = (TLine*)infile->Get("l_"+h_rms);
	  double rmsVal = l_rms->GetY1();

	  //TText *rmsText = new TText(l_baseline->GetY2(), l_horiz->GetX1()+0.2, Form("RMS_{MC} = %f.0",rmsVal));
	  TLatex *rmsText = new TLatex(.1, .1, Form("RMS_{MC} = %.0f",rmsVal));
	  rmsText->SetTextColor(kRed);
	  rmsText->SetTextSize(0.05);
	  //rmsText->DrawTextNDC(.22,.18,Form("RMS of MC = %.0f",rmsVal) );


	}

	idraw++;
	  
      }
	
      c_par->Print(Form("images/%s_cut%d.png",c_par->GetName(),icut));
      //c_par->Print(Form("images/%s.pdf",c_par->GetName()));
    
    }
  }
}    
