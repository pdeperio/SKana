{

  const TString pidTitle[2] = {"1st Sub-Event (#mu)", "2nd Sub-Event (Decay-e)"};

  // Muon momentum binning for likelihoods
  const int nMuonMomBinsLike = 10;
  const double muonMomBinsLike[nMuonMomBinsLike+1] = {0, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 12000};

  // Electron momentum binning for likelihoods
  const int nDcyeMomBinsLike = 8;
  const double dcyeMomBinsLike[nDcyeMomBinsLike+1] = {0, 15, 25, 30, 35, 40, 45, 50, 100};

  TFile *infile = new TFile("/Users/pdeperio/T2K/SK/SoftwareAna/Data/escale/stopmu11d_fitqun_v3r0/output/stopmu_mc_1380_13a.fq_v3r0.root");

  TGraphErrors *g_mean_vs_momsep = 0;
  TGraphErrors *g_rms_vs_momsep = 0;
  int ipoint = 0;
  
  for (int ipar=0; ipar<nPars; ipar++) {

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
	
	TCanvas *c_par = new TCanvas(Form("c_par%d_ise%d_ilike%d",ipar,ise,ilike),Form("c_par%d_ise%d_ilike%d",ipar,ise,ilike),0,0,1800,1000);
	c_par->Divide(nMomBins/2,2);
	
	int idraw = 1;

	for (int imom=0; imom<nMomBins; imom++) {
	  
	  TString h_name = Form("h_cosmic_ise%d_%s_vs_year_reco1_momsep%d_%d_%s_vs_year",ise,likeString[ilike].Data(),(int)momBins[imom],(int)momBins[imom+1],parNames[ipar].Data());
	  TH1D *h_tmp = (TH1D*)infile->Get(h_name);
	  h_tmp->SetTitle(Form("%d to %d %s", (int)(momBins[imom]/momFactor), (int)(momBins[imom+1]/momFactor), momUnit.Data())); 

	  h_tmp->SetLineColor(kBlue);
	  h_tmp->SetMarkerColor(kBlue);
	  c_par->cd(idraw);
	  h_tmp->Draw();

	  TLine *l_horiz = (TLine*)infile->Get("l_"+h_name);
	  l_horiz->Draw("SAME");
	  
	  TLine *l_vert = (TLine*)infile->Get("l_apr_2009_"+h_name);
	  l_vert->Draw("SAME");
	  
	  h_tmp->Draw("same");
	  

	  idraw++;
	  
	}
	
	c_par->Print(Form("images/%s.png",c_par->GetName()));
	//c_par->Print(Form("images/%s.pdf",c_par->GetName()));
    
      }
      
    }    
  }
}
