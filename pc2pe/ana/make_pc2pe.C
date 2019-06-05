{
  const int nFiles = 5;
  
  TString TreeVarNames[nFiles] = {
    "_sk4", "_sk5", "_sk5i", "_sk5n", "_sk5avg"
  }

  int channel;
  float pc2pe;

  TString folder = "figures_june06_hksep/";
  //TString folder = "";
  
  TFile *infile = new TFile(folder+"pc2pe_output.root");
  TTree *t_pc2pe = (TTree*)infile->Get("pc2pe");
  
  for (int ifile=0; ifile<nFiles; ifile++) {

    ofstream pgain_file;
    pgain_file.open(folder+"pgain"+TreeVarNames[ifile]+"_19jun06");
    pgain_file << "   0  0  11146" << endl;  // Header (Version ? NPMTs)

    t_pc2pe->SetBranchAddress("pc2pe"+TreeVarNames[ifile], &pc2pe);
    t_pc2pe->SetBranchAddress("channel", &channel);

    for (int ipmt=0; ipmt<t_pc2pe->GetEntries(); ipmt++) {

      t_pc2pe->GetEntry(ipmt);
      
      pgain_file << setw(5);
      pgain_file << channel;
      pgain_file << " " << pc2pe << endl;
    }
  }

}
