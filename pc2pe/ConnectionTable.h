#ifndef ConnectionTable_hh
#define ConnectionTable_hh

#include "TTree.h"
#include "TString.h"

using namespace std;

class ConnectionTable
{
 public:
  
  ConnectionTable(TString Filename);
  ~ConnectionTable();

  Int_t cableid;
  string supadd;
  string supsubadd;
  Int_t supserial;
  Int_t modserial;
  Int_t hutnum;
  Int_t tkobnum;
  Int_t tkomodadd;
  Int_t qbch;
  Int_t hvcrate;
  Int_t hvmodadd;
  Int_t hvch;
  Float_t oldhv; // applied hv (to change)
  string pmtserial; // serial (check criteria)
  Int_t pmtflag;
  string qbip;
  Int_t pmtx, pmty, pmtz;
  Int_t odpadnum_hut, odpadnum_crate, odpadnum_mod, odpadnum_ch;

  Int_t group;
  Int_t group_arr[11146];

  int LoadPMT(int ipmt) {return tConnection->GetEntry(ipmt);};
  void WriteTree() {tConnection->Write();};

 private:

  TTree *tConnection;

  void SortTree();
  void AddGroupings();
};


#endif

