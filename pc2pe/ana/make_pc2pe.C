{
  const int nFiles = 3;
  
  TString TreeVarNames[nFiles] = {
    "_sk4", "_sk5", "_sk5i"
  }

  int channel;
  float rationorm;
  
  TFile *infile = new TFile("pc2pe_output.root");
  TTree *pc2pe = (TTree*)infile->Get("pc2pe");
  
  for (int ifile=0; ifile<nFiles; ifile++) {

    ofstream pgain_file;
    pgain_file.open("pgain"+TreeVarNames[ifile]+"_19apr10");
    pgain_file << "   0  0  11146" << endl;  // Header (Version ? NPMTs)

    pc2pe->SetBranchAddress("rationorm"+TreeVarNames[ifile], &rationorm);
    pc2pe->SetBranchAddress("channel", &channel);

    for (int ipmt=0; ipmt<pc2pe->GetEntries(); ipmt++) {

      pc2pe->GetEntry(ipmt);
      
      pgain_file << setw(5);
      pgain_file << channel;
      pgain_file << " " << rationorm << endl;
    }
  }

}
