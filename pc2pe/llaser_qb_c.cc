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
#include "TH1F.h"
#include "TH2F.h"
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
  if(argc < 7) {
    printf("Usage: sample input.root output.root RUN_NUMBER ConnectionTableFile StartTime EndTime\n");
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

  float low_intensity_window = -1;//100;

  // SK4 Low intensity
  if (61892<=RUN_NUMBER && RUN_NUMBER<=61895) {
    ontimemin = 1180;
    ontimemax = 1400;
    offtimemin = 480;//400;
    offtimemax = 700;//800;
    if (low_intensity_window>0) {
      ontimemax = ontimemin+low_intensity_window;
      offtimemax = offtimemin+low_intensity_window;
    }
  } 

  // SK5 Low intensity
  else if ((80871<=RUN_NUMBER && RUN_NUMBER<=80875) || RUN_NUMBER==81028) {
    ontimemin = 1000;
    ontimemax = 1300;
    offtimemin = 450;//410;
    offtimemax = 750;//900;
    if (low_intensity_window>0) {
      ontimemax = ontimemin+low_intensity_window;
      offtimemax = offtimemin+low_intensity_window;
    }
  } 

  // SK5 Low intensity (source inverted)
  else if (RUN_NUMBER==80885) {
    ontimemin = 1000;
    ontimemax = 1300;
    offtimemin = 450;//410;
    offtimemax = 750;//900;
    if (low_intensity_window>0) {
      ontimemax = ontimemin+low_intensity_window;
      offtimemax = offtimemin+low_intensity_window;
    }
  } 

  // SK4 High intensity
  else if (RUN_NUMBER==61889) {
    ontimemin = 1140;
    ontimemax = 1200;
    offtimemin = 420;//400;
    offtimemax = 480;//750;
  }

  // SK5 High intensity
  else if (RUN_NUMBER==80877 || RUN_NUMBER==81030) {
    ontimemin = 975;
    ontimemax = 1060; //1038; misses HK PMTs
    offtimemin = 420;//410;
    offtimemax = 505; //483;//600;
  } 

  // SK5 High intensity (source inverted)
  else if (RUN_NUMBER==80884 || RUN_NUMBER==80886) {
    ontimemin = 975;
    ontimemax = 1060; //1038; misses HK PMTs
    offtimemin = 420;//410;
    offtimemax = 505; //483;//600;
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

  float StartTime = atof(argv[5]);
  float EndTime = atof(argv[6]);
  float PadTime = 60;

  const int nCables = 11147;

  // histogram
  TH1F *hhitmap = new TH1F ("hhitmap","HITMAP",11146, 0.5, 11146.5);
  TH1F *htisk = new TH1F ("htisk","TISK",3000, 0, 3000);
  TH1F *hqisk = new TH1F("hqisk", "QISK", 500, -5, 20);
  TH1F *hnqisk = new TH1F("hnqisk", "Nqisk;Nqisk", 1200, 0, 12000);
  TH2F *hnqisk_vs_time = new TH2F("hnqisk_vs_time", "Nqisk Time Dependence;Epoch Time (s);Nqisk", 500, StartTime-PadTime, EndTime+PadTime, 100, 0, 500);
  
  TH1F *ttof = new TH1F ("ttof","Hit Times;T-ToF [ns]",3000, 0, 3000.);
  //TH1F *ttof_channel[nCables-1];
  //for (int ipmt=0; ipmt<nCables-1; ipmt++)
  //   ttof_channel[ipmt]= new TH1F (Form("ttof_channel_%d",ipmt+1),Form("PMT %d Hit Times;T-ToF [ns]",ipmt+1),3000, 0, 3000.);

  TH1F *hnHitsOnTime = new TH1F ("hnHitsOnTime", "Number of On-time Hits;nHits", 1200, 0, 12000);
  TH1F *hQOnTime = new TH1F ("hQOnTime", "Total On-time Charge;Charge [pC]", 500, 400000, 2400000);
  TH2F *hQOnTime_vs_time = new TH2F ("hQOnTime_vs_time", "Total On-time Charge Time Dependence;Epoch Time (s);Charge [pC]", 500, StartTime-PadTime, EndTime+PadTime, 100, 400000, 2400000);

  TH1F *h_nhit_ton = new TH1F("h_nhit_ton", "Per Channel Event Rate (On-time);Channel;Events", nCables, 0, nCables);
  TH1F *h_nhit_toff = new TH1F("h_nhit_toff", "Per Channel Event Rate (Off-time);Channel;Events", nCables, 0, nCables);
  TH1F *h_qisk_ton = new TH1F("h_qisk_ton", "Per Channel Charge (On-time);Channel;Total Charge (pe)", nCables, 0, nCables);
  TH1F *h_qisk_toff = new TH1F("h_qisk_toff", "Per Channel Charge (Off-time);Channel;Total Charge (pe)", nCables, 0, nCables);

  const int nGroups = 35;
  TH1F *h_group_nhit_ton = new TH1F("h_group_nhit_ton", "Per Group Event Rate (On-time);Group;Events", nGroups, 0, nGroups);
  TH1F *h_group_nhit_toff = new TH1F("h_group_nhit_toff", "Per Group Event Rate (Off-time);Group;Events", nGroups, 0, nGroups);
  TH1F *h_group_qisk_ton = new TH1F("h_group_qisk_ton", "Per Group Charge (On-time);Group;Total Charge (pe)", nGroups, 0, nGroups);
  TH1F *h_group_qisk_toff = new TH1F("h_group_qisk_toff", "Per Group Charge (Off-time);Group;Total Charge (pe)", nGroups, 0, nGroups);

  float monitorQrange[2] = {2850, 2950};
  if (RUN_NUMBER>80000) {
    monitorQrange[0] = 570;
    monitorQrange[1] = 670;
  }

  TH2F *h_monitor_q = new TH2F ("h_monitor_q_vs_time","Monitor PMT Charge;Epoch Time (s);Q [pC]", 500, StartTime-PadTime, EndTime+PadTime, 100, monitorQrange[0], monitorQrange[1]);
  TH1F *h_monitor_t = new TH1F ("h_monitor_t","Monitor PMT Hit Time;Hit Time [ns]", 1000, 0, 30000);


  // total number of events
  int ntotal = tree->GetEntries();

  // after selection
  int nsel = 0;
  int nOnTime = 0;

  int nrunsk;

  struct tm t;
  time_t t_of_day;

  // main loop
  for (int i = 0; i < ntotal; i++) {

    // read all branches
    tree->GetEntry(i);

    if(i==0) nrunsk = HEAD->nrunsk;

    t.tm_year = HEAD->ndaysk[0];
    t.tm_mon = HEAD->ndaysk[1]-1;
    t.tm_mday = HEAD->ndaysk[2];
    t.tm_hour = HEAD->ntimsk[0];
    t.tm_min = HEAD->ntimsk[1];
    t.tm_sec = HEAD->ntimsk[2];
    t.tm_isdst = 0;        // Is DST on? 1 = yes, 0 = no, -1 = unknown
    t_of_day = mktime(&t);

    // Laser stability cut
    if (RUN_NUMBER==80885 && t_of_day < 1553583000) continue;
    else if (RUN_NUMBER==81028 && t_of_day < 1555993848) continue;

    int ireduc = BasicReduction(HEAD);
    if(ireduc != 0) continue;

    int icabbit=0xffff;
    bool hit15012=0;
    float q15012;
    float t15012;
    for (int j = 0; j < TQREAL->nhits; j++) {
      int icab = TQREAL->cables[j] & icabbit;
      //if(icab == 15012) if(TQREAL->Q[j] > 100) hit15012 = 1; // for run 78381
      if(icab == 15012) {
        hit15012 = 1; 
        q15012 = TQREAL->Q[j];
        t15012 = TQREAL->T[j];
      }
    }
    if (!hit15012) continue;
    h_monitor_q->Fill(t_of_day, q15012);
    h_monitor_t->Fill(t15012);
   
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

      int group = PMTTable.group_arr[icab-1];

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
      //ttof_channel[icab-1]->Fill(tof);

      if(tof > ontimemin && tof < ontimemax) {
        h_nhit_ton->Fill(icab);
        h_qisk_ton->Fill(icab, qisk);
        h_group_nhit_ton->Fill(group);
        h_group_qisk_ton->Fill(group, qisk);

        nHitsOnTime++;
        QOnTime += qisk;
      }
      if(tof > offtimemin && tof < offtimemax) {
        h_nhit_toff->Fill(icab);
        h_qisk_toff->Fill(icab, qisk);
        h_group_nhit_toff->Fill(group);
        h_group_qisk_toff->Fill(group, qisk);
      }
    }

    hnqisk->Fill(rootread->nqisk);  // Integral of this gives "nlaser" in llaser_qb.F
    hnHitsOnTime->Fill(nHitsOnTime);
    hQOnTime->Fill(QOnTime);

    hnqisk_vs_time->Fill(t_of_day, rootread->nqisk); 
    hQOnTime_vs_time->Fill(t_of_day, QOnTime);

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
