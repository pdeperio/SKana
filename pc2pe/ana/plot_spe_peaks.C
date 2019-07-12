void plot_spe_peaks(TString spe_var = "Mean")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  // PMT Type separated
  const int nPMTtypes = 4;
  int PMTflags[nPMTtypes] = {3, 4, 6, -1};
  TString PMTTypeNames[nPMTtypes] = {"SK2", "SK3", "HK", "All"};
  int PMTtypeColors[nPMTtypes] = {kBlue, kRed, kBlack, kGray+1};

  TFile *infile = new TFile("pc2pe_output.root","update");
  t_pc2pe = (TTree*)infile->Get("pc2pe");


  TH1D *h_spepeaks[nPMTtypes];
  for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++) {
    h_spepeaks[ipmttype] = new TH1D("h_spepeaks"+PMTTypeNames[ipmttype], ";SPE "+spe_var+";Number of Channels", 100, 2.3, 4.1);
    h_spepeaks[ipmttype]->Sumw2();
  }  
  
  TString cutbase = "spe_"+spe_var+"_80467>0 && badchannel==0";

  int ipmttype = nPMTtypes-1;
  t_pc2pe->Project("h_spepeaks"+PMTTypeNames[ipmttype], "spe_"+spe_var+"_80467", cutbase);
  h_spepeaks[ipmttype]->SetLineColor(PMTtypeColors[ipmttype]);
  h_spepeaks[ipmttype]->SetMarkerColor(PMTtypeColors[ipmttype]);
  h_spepeaks[ipmttype]->SetLineWidth(2);
  h_spepeaks[ipmttype]->Draw();

  c1->SetLogy(1);

  for (int ipmttype=0; ipmttype<nPMTtypes-1; ipmttype++) {
    
    t_pc2pe->Project("h_spepeaks"+PMTTypeNames[ipmttype], "spe_"+spe_var+"_80467", cutbase+Form("&& pmtflag_sk5==%d", PMTflags[ipmttype]), "same");

    h_spepeaks[ipmttype]->SetLineWidth(2);
    h_spepeaks[ipmttype]->SetLineColor(PMTtypeColors[ipmttype]);
    h_spepeaks[ipmttype]->SetMarkerColor(PMTtypeColors[ipmttype]);
    h_spepeaks[ipmttype]->Draw("SAME");
  }

  TLegend *leg = new TLegend(0.67, 0.67, 1, 1);
  leg->SetHeader("PMT type (#mu, #sigma)");
		
  for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++) {

    double mean = h_spepeaks[ipmttype]->GetMean();
    double rms = h_spepeaks[ipmttype]->GetRMS();
    
    TF1 *f_gaus = new TF1("f_gaus"+PMTTypeNames[ipmttype], "gaus", mean-rms, mean+rms);
    h_spepeaks[ipmttype]->Fit(f_gaus, "RN");
    f_gaus->SetLineColor(PMTtypeColors[ipmttype]);
    f_gaus->SetLineWidth(1.7);
    f_gaus->Draw("SAME");

    leg->AddEntry(h_spepeaks[ipmttype], PMTTypeNames[ipmttype] + Form(" (%.3f, %.3f)", f_gaus->GetParameter(1),  f_gaus->GetParameter(2), "lp"));
    
  }
  leg->Draw();

  c1->Print("figures/spe_"+spe_var+".pdf");
  
}
