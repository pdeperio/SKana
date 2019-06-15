#include "ConnectionTable.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "TTreeIndex.h"

ConnectionTable::ConnectionTable(TString ConnectionFile) {

  ifstream connect;
  connect.open(ConnectionFile.Data());

  if (!connect.is_open()){
    cerr << "Error: " << ConnectionFile.Data() << " not opened." << endl;
    exit (-1);
  }
  cout << "Opened: " << ConnectionFile.Data() << endl;

  tConnection = new TTree("ConnectionTableRaw", "PMT Information");

  tConnection->Branch("cableid", &cableid, "cableid/I");
  //std::string supadd;
  //std::string supsubadd;
  tConnection->Branch("supserial", &supserial, "supserial/I");
  tConnection->Branch("modserial", &modserial, "modserial/I");
  tConnection->Branch("hutnum", &hutnum, "hutnum/I");
  tConnection->Branch("tkobnum", &tkobnum, "tkobnum/I");
  tConnection->Branch("tkomodadd", &tkomodadd, "tkomodadd/I");
  tConnection->Branch("qbch", &qbch, "qbch/I");
  tConnection->Branch("hvcrate", &hvcrate, "hvcrate/I");
  tConnection->Branch("hvmodadd", &hvmodadd, "hvmodadd/I");
  tConnection->Branch("hvch", &hvch, "hvch/I");
  tConnection->Branch("oldhv", &oldhv, "oldhv/F");
  //std::string pmtserial; 
  tConnection->Branch("pmtflag", &pmtflag, "pmtflag/I");
  //std::string qbip;
  tConnection->Branch("pmtx", &pmtx, "pmtx/I");
  tConnection->Branch("pmty", &pmty, "pmty/I");
  tConnection->Branch("pmtz", &pmtz, "pmtz/I");
  //Int_t odpadnum_hut, odpadnum_crate, odpadnum_mod, odpadnum_ch;

  string line;
  for (Int_t head = 0; head < 52; head++){
      getline(connect, line);
  }//reading header of the connection table

  while (!connect.eof()){

    connect >> cableid >> supadd >> supsubadd >> supserial >> modserial >> hutnum >> tkobnum >> tkomodadd >> qbch >> hvcrate >> hvmodadd >> hvch >> oldhv >> pmtserial >> pmtflag >> qbip >> pmtx >> pmty >> pmtz >> odpadnum_hut >> odpadnum_crate >> odpadnum_mod >> odpadnum_ch;
   
    // ID PMTs only
    if (cableid <1 || cableid > 11146) continue;

    tConnection->Fill();

  }
  connect.close();

  SortTree();
  AddGroupings();
  AddProductionYear();

  tConnection->SetBranchAddress("group", &group);
  tConnection->SetBranchAddress("cableid", &cableid);
}

ConnectionTable::~ConnectionTable () {}

void ConnectionTable::SortTree() {

  cout << "Sorting PMT information tree" << endl;

  // Create index on tunix branch (time)
  Int_t nb_idx = tConnection->BuildIndex("cableid");
  TTreeIndex *att_index = (TTreeIndex*)tConnection->GetTreeIndex();

  TTree* tTmp = (TTree*)tConnection->CloneTree(0);
  tTmp->SetName("ConnectionTable");
  
  for( Long64_t i = 0; i < att_index->GetN(); i++ ) {
    tConnection->GetEntry( att_index->GetIndex()[i] );
    tTmp->Fill();
  }

  tConnection->Delete();
  tConnection = tTmp;

}

void ConnectionTable::AddGroupings() {

  TString GroupingFile = "pmt_groups.txt";

  ifstream groupings;
  groupings.open(GroupingFile.Data());

  if (!groupings.is_open()){
    cerr << "Error: " << GroupingFile.Data() << " not opened." << endl;
    exit (-1);
  }
  cout << "Opened: " << GroupingFile.Data() << endl;

  TBranch *bGroup = tConnection->Branch("group", &group, "group/I"); 
  tConnection->SetBranchAddress("group", &group); 
  tConnection->SetBranchAddress("cableid", &cableid);

  Long64_t nentries = tConnection->GetEntries(); 

  int ientry=0;
  while (!groupings.eof()){

    tConnection->GetEntry(ientry);

    if (ientry >= nentries) break;

    int cableidtmp;

    groupings >> cableidtmp >> group;

    if (cableidtmp != cableid) {
      cerr << "Error " << ientry << ": cableidtmp (" << cableidtmp << ") != cableid (" << cableid << ")" << endl;
      exit (-2);
    }

    group_arr[ientry] = group;
    bGroup->Fill();
    ientry++;
  }
}

void ConnectionTable::AddProductionYear() {

  const int nSKs = 2;

  const TString SKname[nSKs] = {"_sk4", "_sk5"};
  const TString filename[nSKs] = {"ana/const/pmt_prod_year.dat", "ana/const/pmtinf.dat"};

  for (int isk=0; isk<nSKs; isk++) {
    
    ifstream f_prodyear;
    f_prodyear.open(filename[isk]);
    
    if (!f_prodyear.is_open()) {
      cerr << "Error: " << filename[isk].Data() << " not opened." << endl;
      exit (-1);
    }
    cout << "Opened: " << filename[isk].Data() << endl;
    
    TBranch *bProdyear = tConnection->Branch("prodyear"+SKname[isk], &prodyear[isk], "prodyear"+SKname[isk]+"/F");
    tConnection->SetBranchAddress("cableid", &cableid);
    
    Long64_t nentries = tConnection->GetEntries();
    
    int ientry=0;
    while (!f_prodyear.eof()){
      
      tConnection->GetEntry(ientry);
      
      if (ientry >= nentries) break;
      
      int cableidtmp;
      
      f_prodyear >> cableidtmp;
      
      if (cableidtmp != cableid) {
	cerr << "Error " << ientry << ": cableidtmp (" << cableidtmp << ") != cableid (" << cableid << ")" << endl;
	exit (-2);
      }
      
      int pmt_type, pmtx, pmty, pmtz, something;
      float prod_year;
      
      if (!isk)
	f_prodyear >> prod_year >> pmtx >> pmty >> pmtz >> something;
      
      else
	f_prodyear >> pmt_type >> prod_year;
      
      prodyear[isk] = prod_year;
      bProdyear->Fill();
      ientry++;
    }
  }
}
