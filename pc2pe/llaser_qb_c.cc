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
#include "ConnectionTable.h"

float xball = 35.35;
float yball = -70.7;
float zball = 0.;

#include "geotnkC.h"

//#define CWATER 22.305986  // From sample_snld.cc
#define CWATER 21.6438  // From llaser_qb.F

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
  if(argc < 5) {
    printf("Usage: sample input.root output.root RUN_NUMBER ConnectionTableFile\n");
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

  int RUN_NUMBER = atoi(argv[3]);

  // Determined by eye with ttof plot from ana/plot.C
  float ontimemin, ontimemax, offtimemin, offtimemax;

  // SK4 Low intensity
  if (61892<=RUN_NUMBER && RUN_NUMBER<=61895) {
    ontimemin = 1180;
    ontimemax = 1400;
    offtimemin = 480;//400;
    offtimemax = 700;//800;
  } 

  // SK5 Low intensity
  else if (80871<=RUN_NUMBER && RUN_NUMBER<=80875) {
    ontimemin = 1000;
    ontimemax = 1300;
    offtimemin = 450;//410;
    offtimemax = 750;//900;
  } 

  // SK5 Low intensity (source inverted)
  else if (RUN_NUMBER==80885) {
    ontimemin = 1000;
    ontimemax = 1300;
    offtimemin = 450;//410;
    offtimemax = 750;//900;
  } 

  // SK4 High intensity
  else if (RUN_NUMBER==61889) {
    ontimemin = 1140;
    ontimemax = 1200;
    offtimemin = 420;//400;
    offtimemax = 480;//750;
  }

  // SK5 High intensity
  else if (RUN_NUMBER==80877) {
    ontimemin = 975;
    ontimemax = 1038;
    offtimemin = 420;//410;
    offtimemax = 483;//600;
  } 

  // SK5 High intensity (source inverted)
  else if (RUN_NUMBER==80884 || RUN_NUMBER==80886) {
    ontimemin = 975;
    ontimemax = 1038;
    offtimemin = 420;//410;
    offtimemax = 483;//600;
  } 

  else {
    cerr << "Error: RUN_NUMBER=" << RUN_NUMBER << " not implemented" <<endl;
    exit (-1);
  }

  // Root initialize
  mgr->Initialize();

  // obtain tree
  TTree* tree = mgr->GetTree();

  // obtain Header information
  Header    *HEAD    = mgr->GetHEAD();

  // obtain all the TQ information 
  TQReal    *TQREAL  = mgr->GetTQREALINFO();
  //TQReal    *TQAREAL = mgr->GetTQAREALINFO();

  // obtain each analysis infomation
  //LoweInfo  *LOWE    = mgr->GetLOWE();
  //AtmpdInfo *ATMPD   = mgr->GetATMPD();
  //UpmuInfo  *UPMU    = mgr->GetUPMU();
  //MuInfo    *MU      = mgr->GetMU();
  //SLEInfo   *SLE     = mgr->GetSLE();

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
  
  ConnectionTable PMTTable(argv[4]);
  PMTTable.WriteTree();

  const int nCables = 11147;

  // histogram
  TH1D *hhitmap = new TH1D ("hhitmap","HITMAP",11146, 0.5, 11146.5);
  TH1D *htisk = new TH1D ("htisk","TISK",3000, 0, 3000);
  TH1D *hqisk = new TH1D("hqisk", "QISK", 500, -5, 20);
  TH1D *hnqisk = new TH1D("hnqisk", "Nqisk;Nqisk", 1200, 0, 12000);
  
  TH1D *ttof = new TH1D ("ttof","Hit Times;T-ToF [ns]",3000, 0, 3000.);

  TH1D *hnHitsOnTime = new TH1D ("hnHitsOnTime", "Number of On-time Hits;nHits", 1200, 0, 12000);
  TH1D *hQOnTime = new TH1D ("hQOnTime", "Total On-time Charge;Charge [pC]", 500, 400000, 2200000);

  TH1D *h_nhit_ton = new TH1D("h_nhit_ton", "Per Channel Event Rate (On-time);Channel;Events", nCables, 0, nCables);
  TH1D *h_nhit_toff = new TH1D("h_nhit_toff", "Per Channel Event Rate (Off-time);Channel;Events", nCables, 0, nCables);
  TH1D *h_qisk_ton = new TH1D("h_qisk_ton", "Per Channel Charge (On-time);Channel;Total Charge (pe)", nCables, 0, nCables);
  TH1D *h_qisk_toff = new TH1D("h_qisk_toff", "Per Channel Charge (Off-time);Channel;Total Charge (pe)", nCables, 0, nCables);

  // total number of events
  int ntotal = tree->GetEntries();

  // after selection
  int nsel = 0;
  int nOnTime = 0;

  int nrunsk;

  // main loop
  for (int i = 0; i < ntotal; i++) {

    // read all branches
    tree->GetEntry(i);

    if(i==0) nrunsk = HEAD->nrunsk;

    

    int ireduc = BasicReduction(HEAD);
    if(ireduc != 0) continue;

    int icabbit=0xffff;
    bool hit15012=0;
    for (int j = 0; j < TQREAL->nhits; j++) {
      int icab = TQREAL->cables[j] & icabbit;
      //if(icab == 15012) if(TQREAL->Q[j] > 100) hit15012 = 1; // for run 78381
      if(icab == 15012) hit15012 = 1; 
    }
    if (!hit15012) continue;
    
    rootread->skread(); // fill hit information
    
    //cout << "NEVSK: " << HEAD->nevsk << " NQISK: " << rootread->nqisk << endl;
    nsel++;

    int nHitsOnTime = 0;
    float QOnTime = 0;

    for(int j=0; j<rootread->nqisk; j++){

      int icab=rootread->ihcab[j]; // cable number
      if(!icab || icab >= nCables) {
        cout << "Error: icab = " << icab << endl;
        exit (-1);
      }

      float xpmt = rootread->xidsk[j]; // PMT position
      float ypmt = rootread->yidsk[j];
      float zpmt = rootread->zidsk[j];
      
      hhitmap->Fill(icab);

      float tisk = rootread->tidsk[j]; // t
      htisk -> Fill(tisk);
      float qisk = rootread->qidsk[j]; // q
      hqisk -> Fill(qisk);
      //cout << "Q: " << qisk << endl;    
      
      // This shows up in high intensity analysis but not low. Needed?
      //if (qisk<0) continue;

      float distan=sqrt((xball-xpmt)*(xball-xpmt) + (yball-ypmt)*(yball-ypmt) + (zball-zpmt)*(zball-zpmt));
      float tof=rootread->tidsk[j]-distan/CWATER;
      ttof->Fill(tof);
      
      if(tof > ontimemin && tof < ontimemax) {
        h_nhit_ton->Fill(icab);
        h_qisk_ton->Fill(icab, qisk);

        nHitsOnTime++;
        QOnTime += qisk;
      }
      if(tof > offtimemin && tof < offtimemax) {
        h_nhit_toff->Fill(icab);
        h_qisk_toff->Fill(icab, qisk);
      }
    }

    hnqisk->Fill(rootread->nqisk);  // Integral of this gives "nlaser" in llaser_qb.F
    hnHitsOnTime->Fill(nHitsOnTime);
    hQOnTime->Fill(QOnTime);
    if (nHitsOnTime) nOnTime++;

  }
  
  cout << "Total number of events: " << ntotal << " Selected: " << nsel;
  cout << ", On-time: " << nOnTime << endl;
  
  fout->Write();
  fout->Close();

  // close file
  Smgr->DeleteTreeManager(id);

  // end
  SuperManager::DestroyManager();

}
