//
// sample.cc     27-DEC-2011     Y.Koshio
//
//

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

#include <math.h>
#include "TH1D.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"

#include "SuperManager.h"
#include "SKRootRead.h"

float xball = 35.35;
float yball = -70.7;
float zball = 0.;

#include "geotnkC.h"

//#define CWATER 22.305986  // From sample_snld.cc
#define CWATER 21.6438  // From llaser_qb.F


float ontimemin = 1130, ontimemax = 1250, offtimemin = 850., offtimemax = 1100;

int BasicReduction(Header *HEAD){

  //
  // basic reduction (see also skhead.h)
  //

  int ireduc=0;

  // require TQ
  if(!(HEAD->ifevsk & 1)){
    cout << " No QBEE TQ " << endl;
    ireduc++;
  }

  // require TRG
  if(!(HEAD->ifevsk & (1<<1))){
    cout << " No HARD TRG " << endl;
    ireduc++;
  }

  // skip pedestal
  if(HEAD->ifevsk & (1<<9) || HEAD->idtgsk & (1<<30)){
    cout << " Pedestal " << endl;
    ireduc++;
  }

  // skip slow data
  if(HEAD->ifevsk & (1<<14)){
    cout << " slow data " << endl;
    ireduc++;
  }

  // skip run info
  if(HEAD->ifevsk & (1<<15)){
    cout << " run info " << endl;
    ireduc++;
  }

  // spacer
  if(HEAD->ifevsk & (1<<19)){
    cout << " Remove spacer " << endl;
    ireduc++;
  }

  // incomplete
  if(HEAD->ifevsk & (1<<20)){
    cout << " Incomplete " << endl;
    ireduc++;
  }

  // LED burst 
  if(HEAD->ifevsk & (1<<26)){
    cout << " LED burst " << endl;
    ireduc++;
  }

  // skip inner off
  if(HEAD->ifevsk & (1<<28)){
    cout << " inner off " << endl;
    ireduc++;
  }

  // require EVNT_HDR
  if(!(HEAD->ifevsk & (1<<31))){
    cout << " No EVNT_HDR " << endl;
    ireduc++;
  }

  // skip event number zero
  if(HEAD->nevsk == 0){
    cout << " event number zero " << endl;
    ireduc++;
  }

  if(ireduc != 0) cout << " event #" << HEAD->nevsk << " is skipped because of the above " << ireduc << " reasons" << endl;

  return ireduc;

}

int main(int argc, char *argv[])
{
  if(argc < 3) {
    printf("Usage: sample input.root output.root\n");
    exit(0);
  }

  // open skroot file
  int id = 10;

  //
  SuperManager* Smgr = SuperManager::GetManager(); 

  // No output file
  Smgr->CreateTreeManager(id,"\0","\0",2);  // mode=2

  // Set input file
  TreeManager* mgr = Smgr->GetTreeManager(id);
  //for(int ifile=1; ifile < argc; ifile++)
  mgr->SetInputFile(argv[1]);
  
  TString OutputFile = argv[2];


  // Root initialize
  mgr->Initialize();

  // obtain tree
  TTree* tree = mgr->GetTree();

  // obtain Header information
  Header    *HEAD    = mgr->GetHEAD();

  // obtain all the TQ information 
  TQReal    *TQREAL  = mgr->GetTQREALINFO();
  TQReal    *TQAREAL = mgr->GetTQAREALINFO();

  // obtain each analysis infomation
  LoweInfo  *LOWE    = mgr->GetLOWE();
  AtmpdInfo *ATMPD   = mgr->GetATMPD();
  UpmuInfo  *UPMU    = mgr->GetUPMU();
  MuInfo    *MU      = mgr->GetMU();
  SLEInfo   *SLE     = mgr->GetSLE();

  // obtain TQ information within 1.3u sec
  SKRootRead *rootread = new SKRootRead( mgr );

  // skread option
  rootread->SetSkOption("31"); // look at src/skrd/skoptn.F

  //
  // For calibration run, you have to set run number by yourself
  //
  //rootread->SetSkBadchRun(nrunb);
  //rootread->SetSkBadchRun(38764);
  //rootread->SetSkBadchOption(1);

  TFile *fout = new TFile(OutputFile,"RECREATE");
  
  // histogram
  TH1D *hhitmap = new TH1D ("hhitmap","HITMAP",11146, 0.5, 11146.5);
  TH1D *htisk = new TH1D ("htisk","TISK",3000, 0, 3000);
  TH1D *hqisk = new TH1D("hqisk", "QISK", 500, -5, 20);
  
  TH1D *ttof = new TH1D ("ttof","Hit Times;T-ToF [ns]",3000, 0., 3000.);
  TH1D *nhitsub = new TH1D ("nhitsub"," ",100, 0., 25000.);
  TH1D *htdiff = new TH1D ("htdiff", "TIMEDiff", 200, 0, 2);

  TH1D *hnHitsOnTime = new TH1D ("hnHitsOnTime", "Number of On-time Hits;nHits", 500, 0, 12000);
  TH1D *hQOnTime = new TH1D ("hQOnTime", "Total On-time Charge;Charge [pC]", 500, 0, 1000000);


  // initialize of the previous event clock
  int prev_nt48sk[3] = {0,0,0};

  // total number of events
  int ntotal = tree->GetEntries();

  // after selection
  int nsel = 0;
  int nOnTime = 0;

  float nHitsAvg = 0;
  float Qavg = 0;

  // number of hit
  int nhiton[11147] = {0};
  int nhitoff[11147] = {0};

  int nrunsk;

  // main loop
  for (int i = 0; i < ntotal; i++) {

    // read all branches
    tree->GetEntry(i);

    if(i==0) nrunsk = HEAD->nrunsk;

    int ireduc = BasicReduction(HEAD);
    if(ireduc != 0) continue;

    double curclk = (double(HEAD->nt48sk[2]) + double(HEAD->nt48sk[1])*6.5536e4 + double(HEAD->nt48sk[0])*4.294967296e9)*2.0e1;
    double tdiff = 0.;
    if(i > 0) {
      float preclk = (double(prev_nt48sk[2]) + double(prev_nt48sk[1])*6.5536e4 + double(prev_nt48sk[0])*4.294967296e9)*2.0e1;
      tdiff = (curclk - preclk)/1e6;
    }
    htdiff->Fill(tdiff);
    prev_nt48sk[0] = HEAD->nt48sk[0];
    prev_nt48sk[1] = HEAD->nt48sk[1];
    prev_nt48sk[2] = HEAD->nt48sk[2];

    int icabbit=0xffff;
    bool hit15012=0;
    for (int j = 0; j < TQREAL->nhits; j++) {
      int icab = TQREAL->cables[j] & icabbit;
      //if(icab == 15012) if(TQREAL->Q[j] > 100) hit15012 = 1; // for run 78381
      if(icab == 15012) hit15012 = 1; 
    }
    if (!hit15012) continue;
    
    rootread->skread(); // fill hit information
    
    //if (rootread->nqisk < 200 || rootread->nqisk > 400) continue;  // From sample_snld

    //cout << "NEVSK: " << HEAD->nevsk << " NQISK: " << rootread->nqisk << endl;
    nsel++;

    int nHitsOnTime = 0;
    float QOnTime = 0;

    for(int j=0; j<rootread->nqisk; j++){

      int icab=rootread->ihcab[j]; // cable number
      if(icab > 11146) cout << " icab: " << icab << endl;

      float xpmt = rootread->xidsk[j]; // PMT position
      float ypmt = rootread->yidsk[j];
      float zpmt = rootread->zidsk[j];
      
      hhitmap->Fill(icab);

      float tisk = rootread->tidsk[j]; // t
      htisk -> Fill(tisk);
      float qisk = rootread->qidsk[j]; // q
      hqisk -> Fill(qisk);
      //cout << "Q: " << qisk << endl;    
      
      float distan=sqrt((xball-xpmt)*(xball-xpmt) + (yball-ypmt)*(yball-ypmt) + (zball-zpmt)*(zball-zpmt));
      float tof=rootread->tidsk[j]-distan/CWATER;
      ttof->Fill(tof);
      
      if(tof > ontimemin && tof < ontimemax) {
	nhiton[icab]++;
        nHitsOnTime++;
        QOnTime += qisk;
      }
      if(tof > offtimemin && tof < offtimemax) {
	nhitoff[icab]++;
      }
    }

    hnHitsOnTime->Fill(nHitsOnTime);
    hQOnTime->Fill(QOnTime);
    if (nHitsOnTime) nOnTime++;

  }
  
  cout << "Total number of events: " << ntotal << " Selected: " << nsel;
  cout << "On-time: " << nOnTime << endl;
  
  for(int i=1; i<11147;i++){
    int n = nhiton[i] - nhitoff[i];
    nhitsub->Fill(n);
  }
  
  fout->Write();
  fout->Close();

  // close file
  Smgr->DeleteTreeManager(id);

  // end
  SuperManager::DestroyManager();

}
