{
  gStyle->SetPadRightMargin(0.08);

  //gStyle->SetOptFit(111);
  gStyle->SetStatX(.95);
  gStyle->SetStatY(1);
  gStyle->SetStatW(0.4);
  gStyle->SetStatH(0.2);

  const TString pidTitle[2] = {"1st Sub-Event (#mu)", "2nd Sub-Event (Decay-e)"};

  // Muon momentum binning for likelihoods
  const int nMuonMomBinsLike = 7; //10;
  const double muonMomBinsLike[nMuonMomBinsLike+1] = {0, 500, 1000, 1500, 2000, 3000, 4000, 5000}; //, 6000, 8000, 12000};

  const int nMuonLikes = 3;
  const TString muonLikeString[nMuonLikes] = {"Le_Lmu","Le_Lpi0","Lmu_Lpi0"};
  const TString muonLikeTitle[nMuonLikes] = {"ln(L_{e}/L_{#mu})","ln(L_{e}/L_{#pi^{0}})","ln(L_{#mu}/L_{#pi^{0}})"};
  // Electron momentum binning for likelihoods
  const int nDcyeMomBinsLike = 8;
  const double dcyeMomBinsLike[nDcyeMomBinsLike+1] = {0, 15, 25, 30, 35, 40, 45, 50, 100};

  const int nDcyeLikes = 1;
  const TString dcyeLikeString[nDcyeLikes] = {"Le_Lmu"};
  const TString dcyeLikeTitle[nDcyeLikes] = {"ln(L_{e}/L_{#mu})"};
 
  const int nPars = 3;
  const TString parNames[nPars] = {"mean", "sigma", "mispid"};
  const TString parTitles[nPars] = {"Mean", "RMS","Mis-PID Rate (%)"};
  enum parEnum {imean, isigma, imispid};

  TFile *infile = new TFile("pars_vs_mom.root");

  TGraphErrors *g_mean_vs_momsep = 0;
  TGraphErrors *g_rms_vs_momsep = 0;
  int ipoint = 0;
  
  //for (int ipar=0; ipar<1; ipar++) {
  int ipar=0;

    for (int ise=0; ise<2; ise++) {
      
      int nMomBins;
      double *momBins;
      
      int nLikes;
      TString *likeString, *likeTitle;
      
      double momFactor = 1;
      TString momUnit = "MeV/c";
      //if (ise==0) {
      //	momFactor = 1000;
      //	momUnit = "GeV/c";
      //}

      if (ise==0) {
	nMomBins = nMuonMomBinsLike;
	momBins = muonMomBinsLike;
	nLikes = nMuonLikes;
	likeString = muonLikeString;
	likeTitle = muonLikeTitle;
      } else {
	nMomBins = nDcyeMomBinsLike;
	momBins = dcyeMomBinsLike;
	nLikes = nDcyeLikes;
	likeString = dcyeLikeString;
	likeTitle = dcyeLikeTitle;
      }
      
      TLegend *leg = new TLegend(0.3, 0.65, 0.53, 0.88);
      leg->SetFillStyle(0);
      
      for (int ilike=0; ilike<nLikes; ilike++) {
	
	//TCanvas *c_par = new TCanvas(Form("c_par%d_ise%d_ilike%d",ipar,ise,ilike),Form("c_par%d_ise%d_ilike%d",ipar,ise,ilike),0,0,1800,1000);
	//c_par->Divide(nMomBins/2,2);

	TCanvas *c_par = new TCanvas(Form("c_likelihood_par%d_ise%d_ilike%d",ipar,ise,ilike),Form("c_likelihood_par%d_ise%d_ilike%d",ipar,ise,ilike),0,0,700,1000);
	c_par->Divide(2, TMath::Ceil(nMomBins/2.));
	
	int idraw = 1;

	for (int imom=0; imom<nMomBins; imom++) {
	  
	  TString h_name = Form("h_cosmic_ise%d_%s_vs_year_reco1_momsep%d_%d_%s_vs_year",ise,likeString[ilike].Data(),(int)momBins[imom],(int)momBins[imom+1],parNames[ipar].Data());

	  cout << h_name.Data() << endl;
	  TH1D *h_tmp = (TH1D*)infile->Get(h_name);
	  h_tmp->SetTitle(Form("%d to %d %s", (int)(momBins[imom]/momFactor), (int)(momBins[imom+1]/momFactor), momUnit.Data())); 
	  h_tmp->GetYaxis()->SetTitle(Form("%s of %s",parTitles[ipar].Data(),muonLikeTitle[ilike].Data()));

	  h_tmp->SetLineColor(kBlue);
	  h_tmp->SetMarkerColor(kBlue);
	  c_par->cd(idraw);
	  h_tmp->Draw();

	  TString f_name = Form("f_pol1_h_cosmic_ise%d_%s_vs_year_reco1_momsep%d_%d_%s_vs_year",ise,likeString[ilike].Data(),(int)momBins[imom],(int)momBins[imom+1],parNames[ipar].Data());
	  TF1 *f_pol1 = (TF1*)infile->Get(f_name);
	  f_pol1->Draw("SAME");

	  TLine *l_horiz = (TLine*)infile->Get("l_"+h_name);
	  l_horiz->Draw("SAME");
	  
	  TLine *l_apr_2009 = (TLine*)infile->Get("l_apr_2009_"+h_name);
	  l_apr_2009->Draw("SAME");
	  

	  //TLine *l_baseline = (TLine*)infile->Get("l_baseline_"+h_name);
	  //l_baseline->Draw("SAME");

	  double MCparVal = l_horiz->GetY1();
      
	  double f_pol1_p[2];
	  double f_pol1_p_err[2];
	  for (int ipar2=0; ipar2<2; ipar2++) {
	    f_pol1_p[ipar2] = f_pol1->GetParameter(ipar2);
	    f_pol1_p_err[ipar2] = f_pol1->GetParError(ipar2);
	  }

	  double datamc = (f_pol1_p[0]/MCparVal);
	  double datamcerr = 100*datamc*sqrt(pow(f_pol1_p_err[0]/f_pol1_p[0],2));
      
	  datamc = 100*(datamc-1);
      
	  TLatex *t_datamc = new TLatex(0.34,0.28,Form("Data/MC: %.02f%%",datamc));
	  t_datamc->SetTextSize(0.06);
	  t_datamc->SetTextFont(132);
	  t_datamc->SetTextColor(kRed);
	  t_datamc->SetNDC(true);
	  if (ipar==2) t_datamc->Draw("SAME");
	      	      
	  TLatex *t_slope = new TLatex(0.34,0.22,Form("%.02f #pm %.02f%% per year",100*f_pol1_p[1]/f_pol1_p[0], 100*f_pol1_p_err[1]/f_pol1_p[0]));
	  t_slope->SetTextSize(0.06);
	  t_slope->SetTextFont(132);
	  t_slope->SetTextColor(kOrange+2);
	  t_slope->SetNDC(true);
	  if (ipar==2) t_slope->Draw("SAME");
	  

	  if (ipar==0) {
	    
	    TString h_rms = Form("h_cosmic_ise%d_%s_vs_year_reco1_momsep%d_%d_%s_vs_year",ise,likeString[ilike].Data(),(int)momBins[imom],(int)momBins[imom+1],parNames[1].Data());
	  
	    TLine *l_rms = (TLine*)infile->Get("l_"+h_rms);
	    double rmsVal = l_rms->GetY1();

	    //TText *rmsText = new TText(l_baseline->GetY2(), l_horiz->GetX1()+0.2, Form("RMS_{MC} = %f.0",rmsVal));
	    TLatex *rmsText = new TLatex(.1, .1, Form("RMS_{MC} = %.0f",rmsVal));
	    rmsText->SetTextColor(kRed);
	    rmsText->SetTextSize(0.05);
	    rmsText->SetTextFont(132);
	    rmsText->DrawTextNDC(.35,.2,Form("RMS of MC = %.0f",rmsVal) );


	  }

	  h_tmp->Draw("same");
	  

	  idraw++;
	  
	}
	
	//c_par->Print(Form("images/%s.png",c_par->GetName()));
	//c_par->Print(Form("images/%s.pdf",c_par->GetName()));
    
      }
      
    }    
  }
}
