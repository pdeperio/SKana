//
//  SKRootRead.cc
//
//  First version by Y.Koshio 07-FEB-2012
//
//

#include <iostream>

/*
#include <TBranch.h>
#include <TObjArray.h>
*/

using namespace std;

#include "SKRootRead.h"

SKRootRead::SKRootRead( TreeManager * lmgr )
{
  HEAD    = lmgr->GetHEAD();
  TQREAL  = lmgr->GetTQREALINFO();
  MC      = lmgr->GetMC();

  nrunsk_prv = -9999;
  nsubsk_prv = -9999;

  nrunbadch = -9999;
  nbadopt = -1;

  skoptn = "0";
}

void SKRootRead::skread()
{

  //
  // Set header variable
  //

  skhead_.nrunsk = HEAD->nrunsk;
  skhead_.nsubsk = HEAD->nsubsk;
  skhead_.nevsk  = HEAD->nevsk;

  int nrun = HEAD->nrunsk;
  int nsub = HEAD->nsubsk;
  int ierr;

  //
  // Set PMT geometry
  //

  if(HEAD->mdrnsk == 0 || HEAD->mdrnsk == 999999) { // Monte Carlo
    if(MC->ivmcp == 1) {
      skheadg_.sk_geometry = 1;
    }
    else if(MC->ivmcp == 1001) {
      skheadg_.sk_geometry = 1;
    }
    else if(MC->ivmcp == 1002) {
      skheadg_.sk_geometry = 2;
    }
    else if(MC->ivmcp == 1003) {
      skheadg_.sk_geometry = 3;
    }
    else if(MC->ivmcp == 1004) {
      skheadg_.sk_geometry = 4;
    }
    else{
      cout << " Geometry version " << MC->ivmcp
	   << " is not supported... "
	   << endl;
      exit(-1);
    }
  }
  else { // real data
    if(nrunsk_prv != nrun){

      if(nrunsk_prv != -9999) {
	cout << "Run change from " << nrunsk_prv 
	     << " to " << nrun
	     << endl;
      }
      else {cout << "Start Run " << nrun << endl;}

      if( HEAD->sk_geometry < 1 && HEAD->sk_geometry > 5) {
	cout << "SK_GEOMETRY "  << HEAD->sk_geometry << " is not supported" << endl;
	exit(-1);
      }
      skheadg_.sk_geometry = HEAD->sk_geometry;
    }
  }
  geoset_();

  //
  // Set skread option
  //
  //cout << "skoptn : " << skoptn << endl;

  int ibitopt = GetSkOptionBit();
  //cout << ibitopt << endl;

  //
  // Set bad channel
  //

  if(nrunsk_prv != nrun || nsubsk_prv != nsub){

    if(ibitopt & (1<<(31-25))){
      skbadopt_(&nbadopt);
      if(nrunbadch == -9999) {
	skbadch_(&nrun, &nsub, &ierr);
      }
      else {
	int nsubbadch = 1;
	skbadch_(&nrunbadch, &nsubbadch, &ierr);
      }
      cout << "Number of bad channel " << combad_.nbad
	/* << " isqbad " << combad_.isqbad[0]
	   << " ibad-0 " << combad_.ibad[0]
	   << " ibad-4 " << combad_.ibad[4]
	   << " ibad-5 " << combad_.ibad[5] */
	   << endl;
    }
    else {
      cout << "No bad channel mask, Run:" << HEAD->nrunsk << " Sub:" << HEAD->nsubsk << endl;
    }

    nrunsk_prv = HEAD->nrunsk;
    nsubsk_prv = HEAD->nsubsk;
  }

  //
  // class
  //
  nqisk=0;
  qismsk=0;
  qimxsk=0;
  mxqisk=0;
  ihcab.clear();
  tidsk.clear();
  qidsk.clear();

  // PMT geometry
  xidsk.clear();
  yidsk.clear();
  zidsk.clear();
  dxidsk.clear();
  dyidsk.clear();
  dzidsk.clear();
  kidsk.clear();
  iidsk.clear();
  jidsk.clear();

  //
  // fortran
  //
  if(ibitopt & (1<<(31-30))){
    skq_.nqisk = 0;
    skq_.qismsk = 0;
    skq_.qimxsk = 0;
    skq_.mxqisk = 0;
    for(int i=0; i<MAXPM; i++) skq_.qisk[i]=0;

    skt_.timnsk = 100000.;
    skt_.timxsk = -100000.;
    for(int i=0; i<MAXPM; i++) skt_.tisk[i]=0;
  }
  //

  int ifcab[MAXPM]={0};

  int icabbit=0xffff;
  for (int j = 0; j < TQREAL->nhits; j++) {
    int icab = TQREAL->cables[j] & icabbit;
    int iflag = TQREAL->cables[j] >> 16;
    if(iflag & 1) {                   // 1.3usec window
      if(icab > 0 && icab < 11147) { // inner PMT
	if(combad_.ibad[icab-1] == 0) { // good channel
	  if(ifcab[icab-1] == 0) {      // avoid multi hits
	    // class
	    nqisk++;
	    qismsk+=TQREAL->Q[j];
	    if(TQREAL->Q[j] > qimxsk) {
	      qimxsk=TQREAL->Q[j];
	      mxqisk=icab;
	    }
	    ihcab.push_back(icab);
	    tidsk.push_back(TQREAL->T[j]);
	    qidsk.push_back(TQREAL->Q[j]);

	    // fortran
	    if(ibitopt & (1<<(31-30))){
	      skchnl_.ihcab[skq_.nqisk] = icab;
	      skt_.tisk[icab-1] = TQREAL->T[j];
	      skq_.qisk[icab-1] = TQREAL->Q[j];
	      if(TQREAL->Q[j] > qimxsk) {
		skq_.qimxsk=TQREAL->Q[j];
		skq_.mxqisk=icab;
	      }
	      skq_.nqisk++;
	      skq_.qismsk+=TQREAL->Q[j];
	    }

	    // class
	    xidsk.push_back(geopmt_.xyzpm[icab-1][0]);
	    yidsk.push_back(geopmt_.xyzpm[icab-1][1]);
	    zidsk.push_back(geopmt_.xyzpm[icab-1][2]);
	    dxidsk.push_back(geopmt_.dxyzpm[icab-1][0]);
	    dyidsk.push_back(geopmt_.dxyzpm[icab-1][1]);
	    dzidsk.push_back(geopmt_.dxyzpm[icab-1][2]);
	    kidsk.push_back(geopmt_.kijpm[icab-1][0]);
	    iidsk.push_back(geopmt_.kijpm[icab-1][1]);
	    jidsk.push_back(geopmt_.kijpm[icab-1][2]);

	    ifcab[icab-1]=1;
	  }
	  else {
	    cout << " Warning : multi-hit ; Run " << HEAD->nrunsk
		 << "  Sub " << HEAD->nsubsk
		 << "  Ev  " << HEAD->nevsk
		 << "  Cable " << icab
		 << endl;
	  }
	}
      }
    }
  }
}

int SKRootRead::GetSkOptionBit(){

  int ibitopt = 0, iflag=0;
  istringstream iss(skoptn);
  string str;
  while( getline(iss, str, ',') ){
    int num = atoi(str.c_str());
    if(num == 31) iflag =1;
    ibitopt |= (1<<(31-num));
    //cout << num << " " << ibitopt << endl;
  }
  if(iflag != 1) {
    cout << "please set 31 in SkOption" << endl;
    exit(-1);
  }

  return ibitopt;
}

