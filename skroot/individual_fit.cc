#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <TLegend.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TColor.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLatex.h>
#include <TMath.h>
#include <TList.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TRandom3.h>

using namespace std;

TString outdir = "hv_ana";

#ifdef __CINT__
void individual_fit(TString PMTtype = "", TString InputDir = "hv_ana"){
#else
TString PMTtype = "";
TString InputDir = outdir;

int getArgs(int argc, char* argv[]);

int main(int argc, char *argv[]) {

    int args = getArgs(argc, argv);
    if(args != 0){
        std::cerr << "Usage " << std::endl;
        return 0;
    }
#endif

    outdir += PMTtype + "/";

    bool AnalyzeHK = 0;
    if (PMTtype.Contains("hk")) AnalyzeHK = 1;

    gErrorIgnoreLevel = kWarning; // For removing TCanvas::Print msgs
    
    gStyle->SetFrameBorderMode(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetPadColor(0);

    gStyle->SetOptTitle(kTRUE);
    gStyle->SetOptStat("e");
    gStyle->SetOptFit(1111);
    //gStyle->SetPalette(kCool);
    gStyle->SetStatX(0.43);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.18);
    gStyle->SetStatH(0.1);
    //TCanvas * c1 = new TCanvas("c1","c1",1250,1250);
    //TCanvas * c2 = new TCanvas("c2","c2",1000,1250);
    //c1->Divide(5,5);
    //c2->Divide(5,4);
    
    Float_t st0_lower_left_x = 0.1;
    Float_t st0_lower_left_y = 0.65;
    Float_t st_Width = 0.45;
    Float_t st_Height = 0.2;

    const int nfile = 7;
    Int_t runno[] = {80282, 80278, 80275, 80269, 80265, 80263, 80254};
    Double_t hvshift[] = {-75, -50, -25, 0, 25, 50, 75};
    //Double_t threshold[] = {-0.69, -0.69, -0.69, -0.69, -0.69, -0.69, -0.69};
    
    const int nPMTtypes = 3;
    enum PMTtypeEnum {hk, sk2, sk3};
    TString PMTtypeNames[nPMTtypes] = {"HK", "SK2", "SK3"};

    TFile *fin1 = new TFile(InputDir+"hvscan_parameter"+PMTtype+".root", "Read");
    //TFile *fin2 = new TFile("hvscan_parameter_offset.root", "Read");
    
    TFile *fitfin[7];
    for (Int_t i = 0; i < nfile; i++){
        fitfin[i] = new TFile(InputDir+Form("fit_result_%d"+PMTtype+".root", runno[i]), "Read");
    }
    
    Int_t PMTinfo[11146] = {0};
    
    ifstream connect;
    connect.open("connection.super.sk-5.dat");
    std::string line;
    
    for (Int_t head = 0; head < 51; head++){
        getline(connect, line);
    }//reading header of the connection table
    
    if (connect.is_open()){
        while (!connect.eof()){
            Int_t cableid;
            std::string supadd;
            std::string supsubadd;
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
            std::string pmtserial; // serial (check criteria)
            Int_t pmtflag;
            std::string qbip;
            Int_t pmtx, pmty, pmtz;
            Int_t odpadnum_hut, odpadnum_crate, odpadnum_mod, odpadnum_ch;
            
            connect >> cableid >> supadd >> supsubadd >> supserial >> modserial >> hutnum >> tkobnum >> tkomodadd >> qbch >> hvcrate >> hvmodadd >> hvch >> oldhv >> pmtserial >> pmtflag >> qbip >> pmtx >> pmty >> pmtz >> odpadnum_hut >> odpadnum_crate >> odpadnum_mod >> odpadnum_ch;
            
            if (cableid <1 || cableid > 11146) continue;
            PMTinfo[cableid-1]=pmtflag;
            //PMTinfo[cableid-1][1]=pmtz;
            //PMTinfo[cableid-1][2]=oldhv;
        }
	
	connect.close();
	
    } else {
      cout << "Error: Connection file not open" << endl;
      exit (-1);
    }
    
    
    std::vector<Int_t> badchsk;
    std::vector<Int_t> largechsk;
    std::vector<Int_t> zerochsk;
    std::vector<Int_t> okaych;
    
    TTree *tr[nPMTtypes] = {0};
    Int_t nentry[nPMTtypes] = {0};
    Int_t skch[nPMTtypes] = {0};
    
    for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++) {
      if (AnalyzeHK && ipmttype!=hk) continue;
      else if (!AnalyzeHK && ipmttype==hk) continue;

      TString PMTtypeNameLower = PMTtypeNames[ipmttype];
      PMTtypeNameLower.ToLower();

      tr[ipmttype] = (TTree*)fin1->Get("hvscan_"+PMTtypeNameLower);

      nentry[ipmttype] = tr[ipmttype]->GetEntries();
      
      tr[ipmttype]->SetBranchAddress("Channel", &skch[ipmttype]);

      for (Int_t itry = 0; itry < nentry[ipmttype]; itry++){
        tr[ipmttype]->GetEntry(itry);
        okaych.push_back(skch[ipmttype]);
      }
    }
    
    
    ifstream intxt1;
    intxt1.open(InputDir+"badfitting.txt");
    if (intxt1.fail()){
        std::cerr << "Error opening file intxt1" << std::endl;
        exit(1);
    }
   
    while(getline(intxt1,line)){
        if (!intxt1.eof()){
            Int_t channel;
            std::string pmt, norm, beta;
            Double_t lhv, hhv, chi2, prob;
            std::string comment;

            intxt1 >> channel >> pmt >> lhv >> hhv >> norm >> beta >> chi2 >> prob >> comment;
            
            //std::cout << channel << " " << pmt << std::endl;
            okaych.erase(std::remove(okaych.begin(), okaych.end(), channel), okaych.end());

	    if (AnalyzeHK && pmt != PMTtypeNames[hk]) {
	      continue;
	    }
	    else if (!AnalyzeHK && 
		     !(pmt == PMTtypeNames[sk2] || pmt == PMTtypeNames[sk3])) {
	      continue;
	    }
	    
            if (comment == "Failed_to_find_minimum"){
                badchsk.push_back(channel);
                //std::cout << "Add 1 badch for SK PMT" << std::endl;
            }
            else if (comment == "Chi2_>>_1"){
                largechsk.push_back(channel);
                //std::cout << "Add 1 badch for SK PMT" << std::endl;
            }
            else if (comment == "Chi2_=_0"){
                zerochsk.push_back(channel);
            }
        }
    }
    
    Int_t badchsksize = badchsk.size();
    Int_t largechsksize = largechsk.size();
    Int_t zerochsksize = zerochsk.size();
    Int_t okaychsize = okaych.size();
    
    
    TCanvas * c1 = new TCanvas("c1","c1",250*((nfile+1)/2),500);
    c1->Divide((nfile+1)/2,2);

    TString PMTtypeName = "SK";
    if (AnalyzeHK) PMTtypeName = "HK";

    TString CanvasName = outdir+"Failed_HV_Fit_"+PMTtypeName+".pdf";
    
    c1->Print(CanvasName+"[");

    for (Int_t isk =0; isk < badchsksize; isk++){

        Int_t c1divide = isk * (nfile+1) % (8) + 1;
        Int_t c2divide = 0;
	
        //std:: cout << "c1divide: " << c1divide << " isk: " << isk+1 << std::endl;

	TGraphErrors * gr1sk = (TGraphErrors*)fin1->Get(PMTtypeName+Form("_PMT_HVscan_Cable_%06d", badchsk[isk]));
	
        c1->cd(c1divide);
        gr1sk->Draw("AP");
        c1->Update();

	TH1D * hsk1[nfile];
	
        for (Int_t i = 0; i < nfile; i++){

	    c2divide = c1divide + i + 1;
            hsk1[i] = (TH1D*)fitfin[i]->Get(Form("h_spe_onoff_%d",badchsk[isk]));

	    std:: cout << "c2divide: " << c2divide << std::endl;
            c1->cd(c2divide);
            hsk1[i]->Draw();
            c1->Update();

	    if (c2divide < (nfile+1)) {
	      //std::cout << "Drawing c1 " << c2divide << std::endl;
            }
            else {
                c1->Print(CanvasName);
                c1->Clear();
                c1->Divide((nfile+1)/2,2);
            }
        }
	
        if (isk == badchsksize - 1 && c2divide != 0){
            c1->Print(CanvasName);
        }
        //c1count++;
    }
    c1->Modified();
    c1->Print(CanvasName+"]");
    
    c1->Clear();
    c1->Divide((nfile+1)/2,2);

    CanvasName = outdir+"Large_HV_Chi2_"+PMTtypeName+".pdf";
    c1->Print(CanvasName+"[");

    for (Int_t isk =0; isk < largechsksize; isk++){

        Int_t c1divide = isk * (nfile+1) % 8 + 1;
        Int_t c2divide = 0;
        //std:: cout << "c1divide: " << c1divide << " isk: " << isk+1 << std::endl;

	TGraphErrors * gr1sk = (TGraphErrors*)fin1->Get(PMTtypeName+Form("_PMT_HVscan_Cable_%06d", largechsk[isk]));
	
        c1->cd(c1divide);
        gr1sk->Draw("AP");
        c1->Update();
	
        TH1D * hsk1[nfile];
	
        for (Int_t i = 0; i < nfile; i++){
	  
            c2divide = c1divide + i + 1;

	    hsk1[i] = (TH1D*)fitfin[i]->Get(Form("h_spe_onoff_%d",largechsk[isk]));

	    //std:: cout << "c2divide: " << c2divide << std::endl;
            c1->cd(c2divide);
            hsk1[i]->Draw();
            c1->Update();

	    if (c2divide < (nfile+1)) {
	        //std::cout << "Drawing c1 " << c2divide << std::endl;
            }
            else {
                c1->Print(CanvasName);
                c1->Clear();
                c1->Divide((nfile+1)/2,2);
            }
        }
        if (isk == largechsksize - 1 && c2divide != 0){
            c1->Print(CanvasName);
        }
        //c1count++;
    }
    c1->Modified();
    c1->Print(CanvasName+"]");
    
    c1->Clear();
    c1->Divide((nfile+1)/2,2);

    CanvasName = outdir+"Zero_HV_Chi2_"+PMTtypeName+".pdf";
    c1->Print(CanvasName+"[");

    for (Int_t isk =0; isk < zerochsksize; isk++){

        Int_t c1divide = isk * (nfile+1) % 8 + 1;
        Int_t c2divide = 0;

	//std:: cout << "c1divide: " << c1divide << " isk: " << isk+1 << std::endl;

	TGraphErrors * gr1sk = (TGraphErrors*)fin1->Get(PMTtypeName+Form("_PMT_HVscan_Cable_%06d", zerochsk[isk]));
	
        c1->cd(c1divide);
        gr1sk->Draw("AP");
        c1->Update();

	TH1D * hsk1[nfile];

	for (Int_t i = 0; i < nfile; i++){
	  
            c2divide = c1divide + i + 1;

	    hsk1[i] = (TH1D*)fitfin[i]->Get(Form("h_spe_onoff_%d",largechsk[isk]));

	    //std:: cout << "c2divide: " << c2divide << std::endl;
            c1->cd(c2divide);
            hsk1[i]->Draw();
            c1->Update();

	    if (c2divide < (nfile+1)) {
	      //std::cout << "Drawing c1 " << c2divide << std::endl;
            }
            else {
                c1->Print(CanvasName);
                c1->Clear();
                c1->Divide((nfile+1)/2,2);
                
            }
        }
        if (isk == zerochsksize - 1 && c2divide != 0){
            c1->Print(CanvasName);
        }
        //c1count++;
    }
    c1->Modified();
    c1->Print(CanvasName+"]");
    
    c1->Clear();
    c1->Divide((nfile+1)/2,2);

    CanvasName = outdir+"Okay_fit_"+PMTtypeName+".pdf";
    c1->Print(CanvasName+"[");

    for (Int_t iok =0; iok < okaychsize; iok++){

        if (iok%100 == 0) cout << "Plotting #" << iok << endl;
      
        Int_t c1divide = iok * (nfile+1) % 8 + 1;
        Int_t c2divide = 0;

	//std:: cout << "c1divide: " << c1divide << " iok: " << iok+1 << "okchannel: " << okaych[iok] << std::endl;

	TGraphErrors * gr1ok;
        Int_t okchannel = okaych[iok];

        gr1ok = (TGraphErrors*)fin1->Get(PMTtypeName+Form("_PMT_HVscan_Cable_%06d", okchannel));

	c1->cd(c1divide);
        gr1ok->Draw("AP");
        c1->Update();

	TH1D * hsk1[nfile];

	for (Int_t i = 0; i < nfile; i++){
	  
            c2divide = c1divide + i + 1;
            hsk1[i] = (TH1D*)fitfin[i]->Get(Form("h_spe_onoff_%d",okchannel));

	    //std:: cout << "c2divide: " << c2divide << std::endl;
            c1->cd(c2divide);
            hsk1[i]->Draw();
            c1->Update();

	    if (c2divide < (nfile+1)) {
	      //std::cout << "Drawing c1 " << c2divide << std::endl;
            }
            else {
                c1->Print(CanvasName);
                c1->Clear();
                c1->Divide((nfile+1)/2,2);
                
            }
        }
        if (iok == largechsksize - 1 && c2divide != 0){
            c1->Print(CanvasName);
        }
        //c1count++;
    }
    c1->Modified();
    c1->Print(CanvasName+"]");
    
    fin1->Close();

    for (Int_t i = 0; i < nfile; i++){
        fitfin[i]->Close();
    }
}

#ifndef __CINT__
    int getArgs(int argc, char* argv[]){
        
        while( (argc > 1) && (argv[1][0] == '-') ){
            switch(argv[1][1]){
	      
	        case 't':
                    PMTtype = argv[2];
		    InputDir += PMTtype + "/";
                   ++argv; --argc;
                    break;
		    
	        case 'i':
                    InputDir = argv[2];
                    ++argv; --argc;
                    break;
            }
            
            ++argv; --argc;
        }
        
        return 0;
        
    }
#endif
