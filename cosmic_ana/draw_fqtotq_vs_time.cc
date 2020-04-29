{
  gStyle->SetPadRightMargin(0.08);

  gStyle->SetOptFit(0);
  gStyle->SetStatX(.95);
  gStyle->SetStatY(1);
  gStyle->SetStatW(0.4);
  gStyle->SetStatH(0.2);

  const int nPars = 3;
  const TString parNames[nPars] = {"mean", "sigma", "mispid"};
  const TString parTitles[nPars] = {"Mean", "#sigma","Mis-PID Rate (%)"};
  enum parEnum {imean, isigma, imispid};

  TFile *infile = new TFile("pars_vs_mom.root");

  double refYear = 2009+3./12;
  
  for (int ipar=0; ipar<2; ipar++) {
      
      int nMomBins;
      double *momBins;
      
      int nLikes;
      TString *likeString, *likeTitle;
      
      TCanvas *c_par = new TCanvas(1);
      
      int idraw = 1;

      TString h_name = Form("h_cosmic_fqtotqnorm_vs_year_reco1_%s_vs_year",parNames[ipar].Data());
      cout << h_name.Data() << endl;
      TH1D *h_tmp = (TH1D*)infile->Get(h_name);

      TString f_name = Form("f_pol1_h_cosmic_fqtotqnorm_vs_year_reco1_%s_vs_year",parNames[ipar].Data());
      TF1 *f_pol1 = (TF1*)infile->Get(f_name);

      TLine *l_horiz = (TLine*)infile->Get("l_"+h_name);
      double MCval = l_horiz->GetY1();

      TF1 *f2 = new TF1("f2",Form("1.023 + 0.02341*(x-%f) - 0.001287*pow(x-%f,2)",refYear,refYear),2008,2014);
      f2->SetTitle(";Year;Charge Rate Data/MC Ratio");
      f2->Draw();

      double xBaseline = f2->GetX(1,2008.1,2008.5,1e-10,500);

      TLine *l_baseline_mean = new TLine(xBaseline, 0.985, xBaseline, 1.11);
      l_baseline_mean->SetLineWidth(3);
      //l_baseline_mean->SetLineColor(kGreen-2);
      l_baseline_mean->SetLineColor(kRed);
      l_baseline_mean->Draw("SAME");
      
      TLine *l_unity = new TLine(2008, 1, 2014, 1);
      l_unity->SetLineWidth(3);
      //l_unity->SetLineColor(kGreen-2);
      l_unity->SetLineColor(kRed);
      l_unity->Draw("SAME");

      h_tmp->Scale(1/MCval);
      h_tmp->GetYaxis()->SetTitle("Charge Rate Data/MC Ratio");
      h_tmp->SetLineColor(kBlue);
      h_tmp->SetMarkerColor(kBlue);
      h_tmp->Draw("Same");


      //f_pol1->Draw("SAME");
      //l_horiz->Draw("SAME");
      
      //TLine *l_apr_2009 = (TLine*)infile->Get("l_apr_2009_"+h_name);
      //l_apr_2009->Draw("SAME");
      
      
      //TLine *l_baseline = (TLine*)infile->Get("l_baseline_"+h_name);
      //l_baseline->Draw("SAME");
      
      //h_tmp->Draw("same");
      
      c_par->Print(Form("images/h_fqtotq_%s.png",c_par->GetName()));
  

  }
}
