{

  //const int nRuns = 7;
  //TString RunNumbers[nRuns] = {
  //  "080116",
  //  "080122",
  //  "080124",
  //  "080134",
  //  "080136",
  //  "080138",
  //  "080140",
  //  "080142"
  //};
  
  //const int nRuns = 7;
  //TString RunNumbers[nRuns] = {
  //  "080148",
  //  "080152",
  //  "080157",
  //  "080159",
  //  "080161",
  //  "080163",
  //  "080167"
  //};
  
  //float Livetimes[nRuns] = {
  //  958,
  //  1400,
  //  961,
  //  1143,
  //  1069,
  //  1064,
  //  919
  //};
   
  const int nRuns = 8;
  TString RunNumbers[nRuns] = {
    "080249",
    "080254",
    "080263",
    "080265",
    "080269",
    "080275",
    "080278",
    "080282"
  };

  float Livetimes[nRuns] = {
    1393,
    994,
    1221,
    1291,
    1087,
    1349,
    1075,
    1082
  };

  TString RunNames[nRuns] = {
    "Nominal (Start)",
    "+75 V",
    "+50 V",
    "+25 V",
    "Nominal",
    "-25 V",
    "-50 V",
    "-75 V"
    //"Nominal (After)"
  };


  TLegend *leg = new TLegend(0.5, 0.4, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  leg->SetHeader("#DeltaVoltage, Total Hit Rate");
  
  bool bNorm = 0;
  
  for (int irun=0; irun<nRuns; irun++) {

    TFile *infile = new TFile("output/spe_individual_"+RunNumbers[irun]+".root");

    hqisk->SetLineWidth(2);
    hqisk->SetLineColor(irun+1);

    double integral = hqisk->Integral();

    if (bNorm) {
      hqisk->Scale(1/integral);
      hqisk->GetYaxis()->SetTitle("Hit rate (Normalized)");
      hqisk->GetYaxis()->SetRangeUser(0, 0.039);
    }
    
    else {
      hqisk->Scale(1./Livetimes[irun]);
      hqisk->GetYaxis()->SetTitle("Hit rate (/second)");
      hqisk->GetYaxis()->SetRangeUser(0, 11000);
    }
    
    hqisk->GetXaxis()->SetTitle("Charge (pC)");
    hqisk->SetTitle("Charge Distribution Across All PMTs");
    
    if (!irun)
      hqisk->Draw();

    else 
      hqisk->Draw("SAME");

 
    leg->AddEntry(hqisk, RunNames[irun]+Form(", %.0f kHz", integral/Livetimes[irun]/1000), "l");
    
  }
  leg->Draw();

  TString cname = "spe_ana/spe_hv_scan_190208";
  if (bNorm) cname += "_norm";
  
  c1->Print(cname+".pdf");
}
