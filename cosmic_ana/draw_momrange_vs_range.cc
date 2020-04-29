{
  gStyle->SetPadRightMargin(0.08);

  const TString pidTitle[2] = {"1st Sub-Event (#mu)", "2nd Sub-Event (Decay-e)"};

  const int nRanges = 5; //7;
  const double rangeLow = 0;
  const double rangeHigh = 25; //35; // cm
  const double rangeStep = 5;  // (rangeHigh-rangeLow)/nRanges
 
  const int nPars = 2;
  const TString parNames[nPars] = {"mean", "sigma"};
  const TString parTitles[nPars] = {"Mean", "#sigma"};
  enum parEnum {imean, isigma};

  const double yRange[2][nPars] = {300, 100,
				   2, 2};

  const int nFiles = 2;
  TFile *infile[nFiles];
  infile[1] = new TFile("images_v3r0/pars_vs_mom.root");
  //infile[0] = new TFile("images_v3r1/pars_vs_mom.root");
  infile[0] = new TFile("pars_vs_mom.root");

  TGraphErrors *g_mean_vs_momsep = 0;
  TGraphErrors *g_rms_vs_momsep = 0;
  int ipoint = 0;

  bool doZoom = 1;
  
  
  for (int ise=1; ise<2; ise++) {
      
    //if (doZoom && ise) continue;
          
    int nMomBins;
    double *momBins;
      
    int nLikes;
    TString *likeString, *likeTitle;
      
    double momFactor = 1;
    TString momUnit = "MeV/c";
    if (ise==0) {
      momFactor = 1000;
      momUnit = "GeV/c";
    }

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

    TCanvas *c_par = new TCanvas(Form("c_ise%d",ise),Form("c_ise%d",ise),0,0,1200,1000);
    c_par->Divide(2,2);

    int idraw = 1;

    for (int ipar=0; ipar<nPars; ipar++) {

      //TCanvas *c_par;

      //if (!doZoom) {
      //  c_par = new TCanvas("c_"+parNames[ipar],"c_"+parNames[ipar],0,0,1200,500);
      //  c_par->Divide(2,1);
      //} else 
      //c_par = new TCanvas("c_"+parNames[ipar],"c_"+parNames[ipar],0,0,600,500);
          
      TLegend *leg = new TLegend(0.3, 0.65, 0.53, 0.88);
      leg->SetFillStyle(0);

      for (int ilike=0; ilike<nLikes; ilike++) {
	
	if (ipar==imispid && ilike>0) continue;

	for (int ifile=0; ifile<nFiles; ifile++) {	       
	  
	  TGraphErrors *g_par = new TGraphErrors(nMomBins);
	  g_par->SetName(Form("g_ise%d_ipar%d_ilike%d_ifile%d",ise,ipar,ilike,ifile));

	  for (int imom=0; imom<nMomBins; imom++) {
	    
	    double midMom = (momBins[imom+1]+momBins[imom])/2/momFactor;
	    double midMomErr = (momBins[imom+1]-momBins[imom])/2/momFactor;

	    TString f_name = Form("f_pol1_h_cosmic_ise%d_%s_vs_year_reco1_momsep%d_%d_%s_vs_year",ise,likeString[ilike].Data(),(int)momBins[imom],(int)momBins[imom+1],parNames[ipar].Data());
	    cout << f_name.Data() << endl;

	    
	    TF1 *f_pol1 = (TF1*)infile[ifile]->Get(f_name);
	    
	    double p0 = f_pol1->GetParameter(0);
	    double p1 = f_pol1->GetParameter(1);
	    double p0err = f_pol1->GetParError(0);
	    double p1err = f_pol1->GetParError(1);
	    
	    //g_par->SetPoint(imom, midMom, p1/fabs(p0)*100);
	    //g_par->SetPointError(imom, midMomErr, fabs((p1/p0*100)*sqrt(pow(p1err/p1,2)+pow(p0err/p0,2))));
	    
	    g_par->SetPoint(imom, midMom, p1);
	    g_par->SetPointError(imom, midMomErr, p1err);
	    
	  }
	  
	  c_par->cd(idraw);//->SetLogx(1);
	  if (!ilike) g_par->Draw("AP");
	  else g_par->Draw("same P");
	  //g_par->GetYaxis()->SetRangeUser(-yRange[ise][ipar],yRange[ise][ipar]);

	  if (ise==0 && doZoom) {
	    double zoomRangeY = 20;
	    g_par->GetXaxis()->SetRangeUser(0,2);
	    g_par->GetYaxis()->SetRangeUser(-zoomRangeY/(ipar+1),zoomRangeY/(ipar+1));
	  }
	
	  if (ipar==imispid) {
	    if (ise==0) g_par->GetYaxis()->SetRangeUser(-0.1,0.1);
	    else if (ise==1) g_par->GetYaxis()->SetRangeUser(-0.4,0.4);
	  }
	
	  g_par->SetLineStyle(ifile+1);
	  g_par->SetLineColor(ilike+1);
	  g_par->SetMarkerStyle(ifile*6+20);
	  g_par->SetMarkerColor(ilike+1);
	  g_par->SetTitle(pidTitle[ise]);
	  g_par->GetXaxis()->SetTitle(Form("Momentum (%s)",momUnit.Data()));

	  leg->AddEntry(g_par, likeTitle[ilike], "pl");

	  if (ise==1) g_par->GetYaxis()->SetTitle(Form("d[%s of %s] / dt  (1/Year)",parTitles[ipar].Data(), likeTitle[ilike].Data()));
	  else {
	    g_par->GetYaxis()->SetTitle(Form("d[%s of ln(L_{a}/L_{b}] / dt  (1/Year)",parTitles[ipar].Data()));
	  }
	}
      }
      //if (ise==0) leg->Draw("SAME");
      
      idraw++;
      //if (!doZoom) 
      //	c_par->Print(Form("images/%s_ise%d.png",c_par->GetName(),ise));
      //else
      //	c_par->Print(Form("images/%s_ise%d_zoom.png",c_par->GetName(),ise));
    }

    if (!doZoom) {
      c_par->Print(Form("images/%s.png",c_par->GetName()));
      //c_par->Print(Form("images/%s.pdf",c_par->GetName()));
    } else {
      c_par->Print(Form("images/%s_zoom.png",c_par->GetName()));
    }
    
  }
}
