{

  gStyle->SetMarkerSize(0.1);
  //gStyle->SetOptFit(0);
  
  const int nFiles = 3;
  TString FileNames[nFiles] = {
    "hvtable_1.8e7_sorted.txt",
    "sk5GainHV_80249-80282_gausspeak.txt",
    "targethv.txt"
  }

  TGraph *gr_target[nFiles];

  TCanvas *c1 = new TCanvas(1);

  // Plot absolute values
  for (int ifile=0; ifile<nFiles; ifile++) {

    gr_target[ifile] = new TGraph(FileNames[ifile]);

    gr_target[ifile]->SetMarkerColor(ifile+1);

    gr_target[ifile]->SetTitle("Target HV Comparison @ 1.8e7 gain;Channel;High Voltage (V)");
  
    if (!ifile) 
      gr_target[ifile]->Draw("AP");
    
    else 
      gr_target[ifile]->Draw("sameP");
  
  }

  // Plot differences
  TGraph *gr_diff[nFiles];
  for (int ifile=1; ifile<nFiles; ifile++) {
    gr_diff[ifile] = new TGraph();
    gr_diff[ifile]->SetTitle("Difference in Target HVs;Channel;#DeltaHV (V)");
  }
  
  for (int ipoint=0; ipoint<gr_target[0]->GetN(); ipoint++) {
    
    double x0, y0;
    gr_target[0]->GetPoint(ipoint, x0, y0);

    for (int ifile=1; ifile<nFiles; ifile++) {

      double x1, y1;

      gr_target[ifile]->GetPoint(ipoint, x1, y1);

      if (fabs(x0-x1) > 0.0001) {
	cout << "Error: Channel numbers not equal " << x0 << " " << x1 << endl;
	exit (-1);
      }

      gr_diff[ifile]->SetPoint(ipoint, x0, y1-y0);
    }
  }

  

  TCanvas *c1 = new TCanvas(1);
  bool isDrawn = 0;
  for (int ifile=1; ifile<nFiles; ifile++) {
    gr_diff[ifile]->SetMarkerColor(ifile+1);
    if (!isDrawn) {
      gr_diff[ifile]->Draw("AP");
      isDrawn = 1;
    }
    gr_diff[ifile]->Draw("sameP");
  }


  // Analyze differences
  TH1F *h_diff[nFiles];
  
  for (int ifile=1; ifile<nFiles; ifile++) {

    h_diff[ifile] = new TH1F(Form("h_diff_%d", ifile), Form("Difference in Target HVs; #DeltaHV (V); Number of Channels", ifile), 500, -2600, 800);
    h_diff[ifile]->SetLineColor(ifile+1);

    for (int ipoint=0; ipoint<gr_diff[ifile]->GetN(); ipoint++) {
      
      double x1, y1;
      gr_diff[ifile]->GetPoint(ipoint, x1, y1);

      h_diff[ifile]->Fill(y1);
    }    
  }
	
 
  // Fit bulk
  TF1 *f_diff[nFiles];

  TCanvas *c1 = new TCanvas(1);
  c1->SetLogy(1);
  
  bool isDrawn = 0;
  for (int ifile=1; ifile<nFiles; ifile++) {

    f_diff[ifile] = new TF1(Form("f_diff_%d", ifile), "gaus", 0, 40);
    f_diff[ifile]->SetLineColor(ifile+1);

    h_diff[ifile]->Fit(f_diff[ifile], "LNR");

    if (!isDrawn) {
      h_diff[ifile]->Draw();
      isDrawn = 1;
    }

    h_diff[ifile]->Draw("same");  
    f_diff[ifile]->Draw("SAME");
  }

  float sigma_threshold = 14;
  cout << endl << "Checking for larger than " << sigma_threshold << " sigma from baseline file: " << FileNames[0].Data() << endl;
  
  for (int ifile=1; ifile<nFiles; ifile++) {

    cout << endl << "Differences in file " << FileNames[ifile].Data() << endl;
    
    for (int ipoint=0; ipoint<gr_diff[ifile]->GetN(); ipoint++) {
      
      double x1, y1;
      gr_diff[ifile]->GetPoint(ipoint, x1, y1);

      float diff_from_mean = y1 - f_diff[ifile]->GetParameter(1);
      float sigma = f_diff[ifile]->GetParameter(2);
      if ( fabs(diff_from_mean) > sigma_threshold*sigma )
	cout << x1 << " " << y1 << endl;
    }
    cout << endl;
  }
 
}
