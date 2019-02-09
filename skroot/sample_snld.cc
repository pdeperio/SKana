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

//float xball = -176.8;
//float yball = -70.7;
float xball = 35.35;
float yball = -70.7;
float zball = 0.;

#include "geotnkC.h"
//float water_hight = 2247.4; //(cm) at 9:00 Nov. 20, 2018
//float water_hight = 3105.4; //(cm) at 8:00 Dec. 3, 2018
//float water_hight = 3500.6;//(cm) at 8:00 Dec. 10, 2018
//float water_hight = 3850.0;

#define CWATER 22.305986


//float ontimemin = 1100., ontimemax = 1500., offtimemin = 300., offtimemax = 1100.; // for run 61893
float ontimemin = -400., ontimemax = 0., offtimemin = -800., offtimemax = -400.; //for run 78559

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

  if(ireduc != 0) cout << " event #" << HEAD->nevsk << " is skipped because of the above reasons" << endl;

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
  TH1D *htisk = new TH1D ("htisk","TISK",500, -900., 600.);
  //TH1D *hqisk = new TH1D ("hqisk","QISK",168, -2., 82.);
  TH1D *hqisk = new TH1D("hqisk", "QISK", 200, -4, 18);
  //TH1D *hqisk = new TH1D("hqisk", "QISK", 500, -5, 20);
  
  TH1D *ttof = new TH1D ("ttof","T-Tof",500, -900., 600.);
  TH1D *nhitsub = new TH1D ("nhitsub"," ",100, 0., 25000.);
  //TH1D *nhitsub = new TH1D ("nhitsub"," ",100, 0., 800.);
  TH1D *htdiff = new TH1D ("htdiff", "TIMEDiff", 200, 0, 2);

  //TH1D *h_spe_all_on = new TH1D ("h_spe_all_on","h_spe_all_on",168, -2, 82);
  //TH1D *h_spe_all_off = new TH1D ("h_spe_all_off","h_spe_all_off",168, -2, 82);
  TH1D *h_spe_all_on = new TH1D ("h_spe_all_on","h_spe_all_on",200, -4, 18);
  TH1D *h_spe_all_off = new TH1D ("h_spe_all_off","h_spe_all_off",200, -4, 18);
  //TH1D *h_spe_all_on = new TH1D ("h_spe_all_on","h_spe_all_on",500, -5, 20);
  //TH1D *h_spe_all_off = new TH1D ("h_spe_all_off","h_spe_all_off",500, -5, 20);
  
  TH1D *h_spe_on[MAXPM];
  TH1D *h_spe_off[MAXPM];

  float onhit550, offhit550;//extract the hit information of a specific channel 550;
  TTree * ton550 = new TTree("cable550_on", "cable550 ontime hits");
  ton550->Branch("onhit550", &onhit550, "onhit550/F");
  TTree * toff550 = new TTree("cable550_off", "cable550 offtime hits");
  toff550->Branch("offhit550", &offhit550, "offhit550/F");
  

  for (Int_t i=0; i<MAXPM; i++){
    //h_spe_on[i] = new TH1D (Form("h_spe_on_%d",i+1),"",168, -2, 82);
    //h_spe_off[i] = new TH1D (Form("h_spe_off_%d",i+1),"",168, -2, 82);             
    h_spe_on[i] = new TH1D (Form("h_spe_on_%d",i+1),"",200, -4, 18);
    h_spe_off[i] = new TH1D (Form("h_spe_off_%d",i+1),"",200, -4, 18);
    //h_spe_on[i] = new TH1D (Form("h_spe_on_%d",i+1),"",500, -5, 20);
    //h_spe_off[i] = new TH1D (Form("h_spe_off_%d",i+1),"",500, -5, 20);
  }

  // initialize of the previous event clock
  int prev_nt48sk[3] = {0,0,0};

  // total number of events
  int ntotal = tree->GetEntries();

  // after selection
  int nsel = 0;

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

    if(tdiff > 0.36 && tdiff < 0.64) {
    //    if (tdiff > 0.84 && tdiff < 1.16) {
      int icabbit=0xffff;
      int hit15024=0;
      for (int j = 0; j < TQREAL->nhits; j++) {
	int icab = TQREAL->cables[j] & icabbit;
	//if(icab == 15024) if(TQREAL->Q[j] > 100) hit15024 = 1; // for run 78381
	if(icab == 15024) hit15024 = 1; // for run 78559
      }

      if(hit15024 == 1) {
	rootread->skread(); // fill hit information

	if(rootread->nqisk > 200 && rootread->nqisk < 400){
	  //cout << "NEVSK: " << HEAD->nevsk << " NQISK: " << rootread->nqisk << endl;
	  nsel++;
	  for(int j=0; j<rootread->nqisk; j++){

	    int icab=rootread->ihcab[j]; // cable number
	    if(icab > 11146) cout << " icab: " << icab << endl;

	    float xpmt = rootread->xidsk[j]; // PMT position
	    float ypmt = rootread->yidsk[j];
	    float zpmt = rootread->zidsk[j];

	    //float water_surface = water_hight - ZPTKTK;

	    //if(zpmt < water_surface) {

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
		h_spe_all_on->Fill(qisk);
		h_spe_on[icab-1]->Fill(qisk);
		if (icab == 550) {
		  onhit550 = qisk;
		  ton550->Fill();
		}
		nhiton[icab]++;
	      }
	      if(tof > offtimemin && tof < offtimemax) {
		h_spe_all_off->Fill(qisk);
		h_spe_off[icab-1]->Fill(qisk);
		if (icab == 550) {
		  offhit550 = qisk;
		  toff550->Fill();
		}
		nhitoff[icab]++;
	      }
	      //}
	  }
	}
      }
    }
  }
  cout << "Total number of events: " << ntotal << " Selected: " << nsel;

  TH1D *h_spe_all_sub = (TH1D*)h_spe_all_on->Clone();
  h_spe_all_sub->SetName("h_spe_all_sub");
  h_spe_all_sub->Add(h_spe_all_off,-1);
  //q00028sub->Add(q00028off,-1);

  for(int i=1; i<11147;i++){
    int n = nhiton[i] - nhitoff[i];
    nhitsub->Fill(n);
  }

  //TFile *fout = new TFile(Form("spe_individual_%d.root",nrunsk),"RECREATE");
  hhitmap->Write();
  htisk->Write();
  hqisk->Write();
  ttof->Write();
  nhitsub->Write();
  htdiff->Write();

  h_spe_all_on->Write();
  h_spe_all_off->Write();
  h_spe_all_sub->Write();

  for(Int_t i=0; i<MAXPM; i++){
    h_spe_on[i]->Write();
    h_spe_off[i]->Write();
  }
  ton550->Write();
  toff550->Write();
  fout->Close();

  // close file
  Smgr->DeleteTreeManager(id);

  // end
  SuperManager::DestroyManager();

}
