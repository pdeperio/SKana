{
  gStyle->SetPadRightMargin(0.08);

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
 
  const int nPars = 5;
  const TString parNames[nPars] = {"gausmean", "sigma", "mean", "rms", "mispid"};
  const TString parTitles[nPars] = {"Gaus. Mean", "#sigma", "Mean", "RMS", "Mis-PID Rate (%)"};
  enum parEnum {igausmean, isigma, imean, irms, imispid};

  const double yRange[2][nPars] = {150, 80, 150, 80, 1,
				   2, 2, 2, 2, 1};


  int iPlot = 1;
  TString plotName = "slope";
  if (iPlot==1) plotName = "datamc";
  else if (iPlot==2) plotName = "absval";
  

  const int nFiles = 1;
  TFile *infile[nFiles];
  //infile[1] = new TFile("images_v3r0/pars_vs_mom.root");
  //infile[0] = new TFile("images_v3r1/pars_vs_mom.root");
  infile[0] = new TFile("pars_vs_mom.root");

  TGraphErrors *g_mean_vs_momsep = 0;
  TGraphErrors *g_rms_vs_momsep = 0;
  int ipoint = 0;

  bool doZoom = 1;
  
  
  for (int ise=0; ise<2; ise++) {
      
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

    /*
    TCanvas *c_par = new TCanvas(Form("c_ise%d",ise),Form("c_ise%d",ise),0,0,1200,500);
    c_par->Divide(2,1);
    /**/

    int idraw = 1;

    for (int ipar=0; ipar<nPars; ipar++) {

      //if (ipar>=imispid) continue;

      
      TCanvas *c_par;

      if (!doZoom) {
        //c_par = new TCanvas(Form("c_likevsmom_%s_ise%d",parNames[ipar].Data(),ise),Form("c_likevsmom_%s_ise%d",parNames[ipar].Data(),ise),0,0,1200,500);
	c_par = new TCanvas(Form("c_likevsmom_%s_ise%d",parNames[ipar].Data(),ise),Form("c_likevsmom_%s_ise%d",parNames[ipar].Data(),ise),0,0,600,500);
        //c_par->Divide(2,1);
      } else 
	c_par = new TCanvas(Form("c_likevsmom_%s_ise%d",parNames[ipar].Data(),ise),Form("c_likevsmom_%s_ise%d",parNames[ipar].Data(),ise),0,0,600,500);
      /**/
    
      //TLegend *leg = new TLegend(0.15, 0.65, 0.33, 0.88);
      TLegend *leg = new TLegend(0.2, 0.18, 0.45, 0.45);
      leg->SetFillStyle(0);
      leg->SetTextFont(132);

      for (int ilike=0; ilike<nLikes; ilike++) {
	
	if (ipar==imispid && ilike>0) continue;

	//if (ilike>1) continue;

	for (int ifile=0; ifile<nFiles; ifile++) {	       
	  
	  TGraphErrors *g_par = new TGraphErrors(nMomBins);
	  g_par->SetName(Form("g_ise%d_ipar%d_ilike%d_ifile%d",ise,ipar,ilike,ifile));

	  TGraphErrors *g_par_mc = new TGraphErrors(nMomBins);
	  g_par_mc->SetName(Form("gmc_ise%d_ipar%d_ilike%d_ifile%d",ise,ipar,ilike,ifile));

	  for (int imom=0; imom<nMomBins; imom++) {
	    
	    double midMom = (momBins[imom+1]+momBins[imom])/2/momFactor;
	    double midMomErr = (momBins[imom+1]-momBins[imom])/2/momFactor;
	    
	    TString h_name = Form("h_cosmic_ise%d_%s_vs_year_reco1_momsep%d_%d_%s_vs_year",ise,likeString[ilike].Data(),(int)momBins[imom],(int)momBins[imom+1],parNames[ipar].Data());

	    TString f_name = Form("f_pol1_%s",h_name.Data());
	    //cout << f_name.Data() << endl;
	    
	    TF1 *f_pol1 = (TF1*)infile[ifile]->Get(f_name);
	    
	    double p0 = f_pol1->GetParameter(0);
	    double p1 = f_pol1->GetParameter(1);
	    double p0err = f_pol1->GetParError(0);
	    double p1err = f_pol1->GetParError(1);

	    //int apr2009bin = h_dist[idist]->FindBin(2009+3./12);
	    //p0 = h_dist[idist]->GetBinContent(apr2009bin);
	    //p0err = h_dist[idist]->GetBinError(apr2009bin);
	    	    
	    TLine *l_horiz = (TLine*)infile[ifile]->Get("l_"+h_name);
	    double MCparVal = l_horiz->GetY1();

	    //g_par->SetPoint(imom, midMom, p1/fabs(p0)*100);
	    //g_par->SetPointError(imom, midMomErr, fabs((p1/p0*100)*sqrt(pow(p1err/p1,2)+pow(p0err/p0,2))));
	    
	    if (iPlot==0) {
	      g_par->SetPoint(imom, midMom, p1);
	      g_par->SetPointError(imom, midMomErr, p1err);
	    } else if (iPlot==1) {

	      //cout << ise << " " << ilike << " " << imom << " " << p0 << " " << MCparVal << " " << p0err << endl;

	      if (MCparVal && MCparVal==MCparVal) {
		//double datamc = p0/MCparVal;
		//double datamcerr = 100*datamc*p0err/p0;
		//datamc = 100*(datamc-1);
		
 		double datamc = p0-MCparVal;
		double datamcerr = p0err;
		
		g_par->SetPoint(imom, midMom, datamc);
		g_par->SetPointError(imom, midMomErr, datamcerr);
	      }
	    } else if (iPlot==2) {
	      g_par->SetPoint(imom, midMom, p0);
	      g_par->SetPointError(imom, midMomErr, p0err);

	      g_par_mc->SetPoint(imom, midMom, MCparVal);
	    }
	    
	    
	  }
	  
	  c_par->cd(idraw);//->SetLogx(1);
	  if (ilike==0 && ifile==0) {
	    g_par->Draw("AP");

	    TLine *l_zero;
	    if (ise==0) l_zero = new TLine(0,0,5,0);
	    else l_zero = new TLine(0,0,110,0);
	    l_zero->SetLineWidth(2);
	    l_zero->SetLineColor(kGray);

	    if (iPlot!=2) l_zero->Draw("SAME");
	    g_par->Draw("same P");
	    
	  }
	  else g_par->Draw("same P");

	  if (iPlot==2) {
 
	    if (ise==0 && ipar!=imispid) g_par_mc->SetLineColor(ilike+1);
	    else g_par_mc->SetLineColor(kRed);

	    g_par_mc->SetLineWidth(2);
	    g_par_mc->Draw("same L");
	      
	  }

	  if (iPlot==0) {
	    g_par->GetYaxis()->SetRangeUser(-yRange[ise][ipar],yRange[ise][ipar]);
	    
	    if (ise==0 && doZoom) {
	      double zoomRangeY = 10;
	      g_par->GetXaxis()->SetRangeUser(0,5);
	      //g_par->GetYaxis()->SetRangeUser(-zoomRangeY/(ipar+1),zoomRangeY/(ipar+1));
	      //g_par->GetYaxis()->SetRangeUser(-zoomRangeY,zoomRangeY);
	      g_par->GetYaxis()->SetRangeUser(-30,20);
	    }
	    
	    if (ipar==imispid) {
	      if (ise==0) g_par->GetYaxis()->SetRangeUser(-0.1,0.1);
	      else if (ise==1) g_par->GetYaxis()->SetRangeUser(-0.4,0.4);
	    }
	  } else if (iPlot==1) {
	    if (ise==0) {
	      if (ipar==imean) g_par->GetYaxis()->SetRangeUser(-800, 300);
	      else if (ipar==isigma) g_par->GetYaxis()->SetRangeUser(-400, 200);

	      if (doZoom) {
		g_par->GetXaxis()->SetRangeUser(0, 5);
		if (ipar==imean) g_par->GetYaxis()->SetRangeUser(-300, 200);
		//g_par->GetXaxis()->SetRangeUser(0, 2);
		//if (ipar==imean) g_par->GetYaxis()->SetRangeUser(-60, 60);
		else if (ipar==isigma) g_par->GetYaxis()->SetRangeUser(-50, 50);


	      }
	    }
	    else {
	      if (ipar==imean) g_par->GetYaxis()->SetRangeUser(-15, 1);
	      else if (ipar==isigma) g_par->GetYaxis()->SetRangeUser(-4, 12);
	    }
	  } else if (iPlot==2) {
	    if (ise==0) {
	      if (ipar==imispid) g_par->GetYaxis()->SetRangeUser(0,1);
	    }
	    else
	      if (ipar==imispid) g_par->GetYaxis()->SetRangeUser(0,5);
	  }

	  
	  g_par->SetLineStyle(ifile+1);
	  g_par->SetLineColor(ilike+1);
	  g_par->SetMarkerStyle(ifile*6+20);
	  g_par->SetMarkerColor(ilike+1);
	  g_par->SetTitle(pidTitle[ise]);
	  g_par->GetXaxis()->SetTitle(Form("Momentum (%s)",momUnit.Data()));
	  
	  leg->AddEntry(g_par, likeTitle[ilike], "pl");
	  
	  if (iPlot==0) {
	    if (ise==1) g_par->GetYaxis()->SetTitle(Form("d[%s of %s] / dt  (1/Year)",parTitles[ipar].Data(), likeTitle[ilike].Data()));
	    else {
	      g_par->GetYaxis()->SetTitle(Form("d[%s of ln(L_{a}/L_{b}] / dt  (1/Year)",parTitles[ipar].Data()));
	    }
	  } else if (iPlot==1) {
	    if (ise==1) g_par->GetYaxis()->SetTitle(Form("Data-MC [%s of %s]",parTitles[ipar].Data(), likeTitle[ilike].Data()));
	    else {
	      g_par->GetYaxis()->SetTitle(Form("Data-MC [%s of ln(L_{a}/L_{b})]",parTitles[ipar].Data()));
	    }
	  } else if (iPlot==2) {
	    if (ipar==imispid) {
	      g_par->GetYaxis()->SetTitle(Form("%s",parTitles[ipar].Data()));
	    } else {
	      if (ise==1) g_par->GetYaxis()->SetTitle(Form("%s of %s",parTitles[ipar].Data(), likeTitle[ilike].Data()));
	      else {
		g_par->GetYaxis()->SetTitle(Form("%s of ln(L_{a}/L_{b})",parTitles[ipar].Data()));
	      }
	    }
	  }

	}
      }
      //if (ise==0) leg->Draw("SAME");
      
      
      idraw++;
      /*
      if (!doZoom) 
      	c_par->Print(Form("images/%s_%s.pdf",c_par->GetName(),plotName.Data()));
      else
      	c_par->Print(Form("images/%s_%s_zoom.pdf",c_par->GetName(),plotName.Data()));
      /**/
    }

    /*
    if (!doZoom) {
      c_par->Print(Form("images/%s_%s.png",c_par->GetName(),plotName.Data()));
      //c_par->Print(Form("images/%s.pdf",c_par->GetName()));
    } else {
      c_par->Print(Form("images/%s_%s_zoom.png",c_par->GetName(),plotName.Data()));
    }
    /**/
    
  }
}
