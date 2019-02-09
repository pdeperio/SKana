#ifndef _SKRootRead_
#define _SKRootRead_

#include <sstream>

//root stuff
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"

// skroot stuff
/*
#include "DataDefinition.h"
*/
#include "TreeManager.h"
#include "tqrealroot.h"

// sk fortran commons
#include "skheadC.h"
#include "geopmtC.h"
#include "skbadcC.h"

#include "sktqC.h"

// sk fortran code
extern "C" void geoset_();
extern "C" void skbadch_(int*, int*, int*);
extern "C" void skbadopt_(int*);

class SKRootRead {

 public:

  SKRootRead( TreeManager * );
  ~SKRootRead();

  //
  // for skread option (see also skoptn.F)
  //
  int GetSkOptionBit();
  inline void SetSkOption(string str){skoptn=str;}
  //
  // for bad channel masking
  //
  inline void SetSkBadchRun(int i){nrunbadch=i;}
  inline void SetSkBadchOption(int i){nbadopt=i;}

  void skread();

  //
  // hit information within 1.3usec for class
  //
  Int_t nqisk;
  Float_t qismsk;
  Float_t qimxsk;
  Int_t mxqisk;
  vector<int> ihcab, kidsk, iidsk, jidsk;
  vector<float> qidsk, tidsk, xidsk, yidsk, zidsk, dxidsk, dyidsk, dzidsk;

 private:
  Header    *HEAD;
  TQReal    *TQREAL;
  MCInfo    *MC;

  int nrunsk_prv;
  int nsubsk_prv;

  string skoptn;
  int nrunbadch, nbadopt;

  //ClassDef(SKRootRead,1) // increase version number when structure changed. 
};


#endif

