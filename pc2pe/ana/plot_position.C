{
  gStyle->SetOptStat(0);

  const int nFiles = 3;

  TString TreeVarNames[nFiles] = {
    "_sk4", "_sk5", "_sk5i"
  }
  int Colors[nFiles] = {kBlack, kRed, kBlue};

  TString FileTitles[nFiles] = {
    "SK4", "SK5", "SK5 Inv."
  };

  TFile *infile = new TFile("pc2pe_output.root");

  TString DrawOpts = "][ P HIST";

  TCanvas *c_pc2pe = new TCanvas(1);
  TLegend *leg = new TLegend(0.2, 0.78, 0.88, 0.88);
  leg->SetNColumns(3);
  
  TH1D *h_pc2pe[nFiles];
  TH1D *h_pc2pe_count[nFiles];

  float minY = 0.95;
  float maxY = 1.07;
  
  for (int ifile=0; ifile<nFiles; ifile++) {

    h_pc2pe[ifile] = (TH1D*)h_group_qisk_ton->Clone();
    TString histname = "group_pc2pe"+TreeVarNames[ifile];
    h_pc2pe[ifile]->SetName(histname);
    h_pc2pe[ifile]->SetTitle("");
    h_pc2pe[ifile]->Reset();

    h_pc2pe_count[ifile] = (TH1D*)h_pc2pe[ifile]->Clone();
    h_pc2pe_count[ifile]->SetName(histname+"_count");
    
    pc2pe->Project(histname, "group", "rationorm"+TreeVarNames[ifile]);
    pc2pe->Project(histname+"_count", "group");
    
    h_pc2pe[ifile]->Divide(h_pc2pe_count[ifile]);
    
    h_pc2pe[ifile]->SetMarkerColor(Colors[ifile]);
    h_pc2pe[ifile]->SetMarkerStyle(ifile+2);
    h_pc2pe[ifile]->SetMarkerSize(1.5);

    h_pc2pe[ifile]->GetYaxis()->SetRangeUser(minY, maxY);
    
    if (!ifile) 
      h_pc2pe[ifile]->Draw(DrawOpts); //DrawOpts += "same";
    else 
      h_pc2pe[ifile]->Draw(DrawOpts+" same");
    
    leg->AddEntry(h_pc2pe[ifile], FileTitles[ifile], "p");

    // Geometry separation
    for (int isep=18; isep<=26; isep+=8) {
      TLine *l_Sep = new TLine(isep, minY, isep, maxY);
      l_Sep->SetLineColor(kGreen-2);
      l_Sep->SetLineWidth(2);
      l_Sep->Draw();
    }

    TLine *l_unity = new TLine(0, 1, 35, 1);
    l_unity->SetLineColor(kGray+1);
    l_unity->Draw();
	  
    
    
  }  

  leg->Draw();
  c_pc2pe->Print("figures/pc2pe_grouped.png");
}
