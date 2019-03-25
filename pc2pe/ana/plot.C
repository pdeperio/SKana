{
  gStyle->SetOptStat(0);
  
  const int nFiles = 2;

  TString datadir = "../output/";
  
  TString FileNames[nFiles] = {
    "pc2pe_tst061892_to_5.root",
    "pc2pe_tst080871.root"
  };

  TString FileTitles[nFiles] = {
    "SK4 (61892-61895)",
    "SK5 (80871)"
  };

  TCanvas *c_ttof = new TCanvas(1);
  c_ttof->SetLogy(1);

  TCanvas *c_nquisk = new TCanvas(1);
  c1->SetLogy(1);

  TLegend *leg = new TLegend(0.2, 0.7, 0.6, 0.85);
  
  
  for (int ifile=0; ifile<nFiles; ifile++) {
    
    TFile *infile = new TFile(datadir+FileNames[ifile]);

    c_ttof->cd();
    ttof->GetXaxis()->SetRangeUser(0, 2000);
    ttof->SetLineColor(ifile+1);
    ttof->SetLineWidth(2);
    ttof->GetYaxis()->SetTitle("Hits/event");
    
    int nevents = hnqisk->Integral();
    ttof->Scale(1./nevents);
    
    if (!ifile) ttof->Draw();
    else ttof->Draw("same");

    leg->AddEntry(ttof, FileTitles[ifile], "l");
    
    c_nquisk->cd();
    hnqisk->GetXaxis()->SetRangeUser(0, 600);
    hnqisk->SetLineColor(ifile+1);
    hnqisk->SetLineWidth(2);
    hnqisk->Scale(1./nevents);
    if (!ifile) hnqisk->Draw();
    else hnqisk->Draw("same");
  }

  c_ttof->cd();
  leg->Draw();
  c_ttof->Print("ttof_sk4_vs_sk5.png");

  //TCanvas *c1 = new TCanvas(1);
  //hQOnTime->Draw();
}
