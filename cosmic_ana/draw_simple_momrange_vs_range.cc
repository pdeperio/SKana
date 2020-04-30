{
  const int nRanges = 7;
  const double rangeLow = 0;
  const double rangeHigh = 35; // cm
  const double rangeStep = 5;  // (rangeHigh-rangeLow)/nRanges
  
  const int nRangeCuts = 4;
  const TString rangeTitle[nRangeCuts] = {
    "Official E-Scale",
    "Top & Side Entering",
    "Top Entering",
    "Side Entering"
  };

  char *h_name[2] = {"h", "hmc"};
  int ireco = 0;
  int _nrgSepType = 0;
  int _maxE = 999999;
  int _fileNum = -1;

  TH1D *h_cosmic_mom_over_range_rangesep[2][nRanges][nRangeCuts];

  TH1D *h_mom_over_range_vs_range[2][nRangeCuts];
  TH1D *h_mom_over_range_vs_range_rat[nRangeCuts];

  TGraphErrors *g_mom_over_range_vs_range[2][nRangeCuts];
  TGraphErrors *g_mom_over_range_vs_range_rat[nRangeCuts];

  TCanvas *c1 = new TCanvas(1);
  bool isDrawn = 0;
  
  for (int imc=0; imc<2; imc++) {

    TString infilename = Form("%s_reco%d_nrg%d_Elt%d_file%d.root",h_name[imc],ireco,_nrgSepType,_maxE,_fileNum);

    TFile *infile = new TFile(infilename);

    for (int icut=0; icut<nRangeCuts; icut++) {

      h_mom_over_range_vs_range[imc][icut] = new TH1D(Form("%s_mom_over_range_vs_range_icut%d", h_name[imc], icut), ";Range (m); Mean Mom./Range (GeV/c / m)", nRanges, rangeLow, rangeHigh);
      h_mom_over_range_vs_range[imc][icut]->Sumw2();

      g_mom_over_range_vs_range[imc][icut] = new TGraphErrors();
      
      for (int irange=0; irange<nRanges; irange++) {
	
	double rangeLowEdge = (rangeLow+irange*rangeStep);
	double rangeHighEdge = (rangeLow+(irange+1)*rangeStep);
	
	TString histName = Form("%s_cosmic_mom_over_range_reco%d_rangesep%03d_%03d_icut%d",h_name[imc],ireco,(int)rangeLowEdge,(int)rangeHighEdge,icut);

	h_cosmic_mom_over_range_rangesep[imc][irange][icut] = (TH1D*)infile->Get(histName);

	TH1D *h_this = h_cosmic_mom_over_range_rangesep[imc][irange][icut];

	h_mom_over_range_vs_range[imc][icut]->SetBinContent(irange+1, h_this->GetMean());
	h_mom_over_range_vs_range[imc][icut]->SetBinError(irange+1, h_this->GetMeanError());

	g_mom_over_range_vs_range[imc][icut]->SetPoint(irange, rangeLowEdge+icut+1+(imc-1)*.2, h_this->GetMean());
	g_mom_over_range_vs_range[imc][icut]->SetPointError(irange, 0, h_this->GetMeanError());
      }
      
      g_mom_over_range_vs_range[imc][icut]->SetLineWidth(2);
      g_mom_over_range_vs_range[imc][icut]->SetLineColor(icut+1);
      g_mom_over_range_vs_range[imc][icut]->SetMarkerColor(icut+1);
      g_mom_over_range_vs_range[imc][icut]->SetMarkerSize(1);
      
      //g_mom_over_range_vs_range[imc][icut]->SetLineStyle(imc+1);
      g_mom_over_range_vs_range[imc][icut]->SetMarkerStyle((imc+1)*4);
      
      c1->cd();
      if (!isDrawn) {
	h_mom_over_range_vs_range[imc][icut]->SetLineColor(0);
	h_mom_over_range_vs_range[imc][icut]->GetYaxis()->SetRangeUser(2, 2.8);
	h_mom_over_range_vs_range[imc][icut]->Draw("HIST");
	g_mom_over_range_vs_range[imc][icut]->Draw("SAME p");
	isDrawn = 1;
      }
      else
	g_mom_over_range_vs_range[imc][icut]->Draw("same p");
      
    }
    
    
  }

  TLegend *leg = new TLegend(0.5, 0.5, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->AddEntry(g_mom_over_range_vs_range[1][0], "Nominal", "lp");
  leg->AddEntry(g_mom_over_range_vs_range[0][0], "PMT Shift", "lp");

  int imc=0;
  for (int icut=0; icut<nRangeCuts; icut++)
    leg->AddEntry(g_mom_over_range_vs_range[imc][icut], rangeTitle[icut], "lp");
  
  leg->Draw();


  

  
  for (int icut=0; icut<nRangeCuts; icut++) {
    h_mom_over_range_vs_range_rat[icut] = (TH1D*)h_mom_over_range_vs_range[0][icut]->Clone();
    h_mom_over_range_vs_range_rat[icut]->Divide(h_mom_over_range_vs_range[1][icut]);
    
    g_mom_over_range_vs_range_rat[icut] = new TGraphErrors();
    for (int irange=0; irange<nRanges; irange++) {
      double rangeLowEdge = (rangeLow+irange*rangeStep);
      double rangeHighEdge = (rangeLow+(irange+1)*rangeStep);
      g_mom_over_range_vs_range_rat[icut]->SetPoint(irange, rangeLowEdge+icut+1, h_mom_over_range_vs_range_rat[icut]->GetBinContent(irange+1));
      g_mom_over_range_vs_range_rat[icut]->SetPointError(irange, 0, h_mom_over_range_vs_range_rat[icut]->GetBinError(irange+1));
    }
  }



  TCanvas *crat = new TCanvas(1);

  h_mom_over_range_vs_range_rat[0]->SetLineColor(0);
  h_mom_over_range_vs_range_rat[0]->GetYaxis()->SetRangeUser(0.97, 1.03);
  h_mom_over_range_vs_range_rat[0]->Draw("HIST");
  h_mom_over_range_vs_range_rat[0]->GetYaxis()->SetTitle("PMT Shift / Nominal Ratio");

  TLine *l_unity = new TLine(rangeLow, 1., rangeHigh, 1.);
  l_unity->Draw();

  for (int icut=0; icut<nRangeCuts; icut++) {
    g_mom_over_range_vs_range_rat[icut]->SetLineWidth(2);
    g_mom_over_range_vs_range_rat[icut]->SetLineColor(icut+1);
    g_mom_over_range_vs_range_rat[icut]->SetMarkerColor(icut+1);
    g_mom_over_range_vs_range_rat[icut]->SetMarkerSize(1);

    g_mom_over_range_vs_range_rat[icut]->Draw("SAME p");	
  }
  
  
}
