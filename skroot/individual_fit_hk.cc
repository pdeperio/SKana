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

void individual_fit_hk(){
    gStyle->SetOptTitle(kTRUE);
    gStyle->SetOptStat("e");
    gStyle->SetOptFit(1111);
    //gStyle->SetPalette(kCool);
    gStyle->SetStatX(0.43);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.18);
    gStyle->SetStatH(0.1);
    //TCanvas * c2 = new TCanvas("c2","c2",1000,1250);
    //c2->Divide(5,4);
    
    Float_t st0_lower_left_x = 0.1;
    Float_t st0_lower_left_y = 0.65;
    Float_t st_Width = 0.45;
    Float_t st_Height = 0.2;
   
    //Int_t runno[] = {78565, 78561, 78559, 78568};
    //Double_t hvshift[] = {-100, -50, 0, 50};
    //Double_t threshold[] = {-0.69, -0.69, -0.69, -0.72};
    const int nfile = 7;
    
    Int_t runno[] = {80282, 80278, 80275, 80269, 80265, 80263, 80254};
    Double_t hvshift[] = {-75, -50, -25, 0, 25, 50, 75};
    Double_t threshold[] = {-0.69, -0.69, -0.69, -0.69, -0.69, -0.69, -0.69};
    
    TFile *fin1 = new TFile("hvscan_parameter.root", "Read");
    //TFile *fin2 = new TFile("hvscan_parameter_offset.root", "Read");
    TFile *fitfin[7];
    for (Int_t i = 0; i < nfile; i++){
        fitfin[i] = new TFile(Form("fit_result_%d.root", runno[i]), "Read");
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
    }
    
    std::vector<Int_t> badchhk;
    std::vector<Int_t> largechhk;
    std::vector<Int_t> okaych;
    
    TTree * trhk = (TTree*)fin1->Get("hvscan_hk");
    Int_t nentry_hk = trhk->GetEntries();
    Int_t hkch;
    trhk->SetBranchAddress("Channel",&hkch);
    
    for (Int_t itryhk = 0; itryhk < nentry_hk; itryhk++){
        trhk->GetEntry(itryhk);
        okaych.push_back(hkch);
    }
    
    ifstream intxt1;
    intxt1.open("badfitting.txt");
    if (intxt1.fail()){
        std::cerr << "Error opening file intxt1" << std::endl;
        exit(1);
    }
    //std::string line;    
    while(getline(intxt1,line)){
        if (!intxt1.eof()){
            Int_t channel;
            std::string pmt, norm, beta;
            Double_t lhv, hhv, chi2, prob;
            std::string comment;
            
            intxt1 >> channel >> pmt >> lhv >> hhv >> norm >> beta >> chi2 >> prob >> comment;
            
            std::cout << channel << " " << pmt << std::endl;
            okaych.erase(std::remove(okaych.begin(), okaych.end(), channel), okaych.end());
            
            if (pmt == "HK" && comment == "Failed_to_find_minimum"){
                badchhk.push_back(channel);
                //std::cout << "Add 1 badch for HK PMT" << std::endl;
            }
            else if (pmt == "HK" && comment == "Chi2_>>_1"){
                largechhk.push_back(channel);
                //std::cout << "Add 1 badch for HK PMT" << std::endl;
            }
        }
    }
    
    Int_t badchhksize = badchhk.size();
    Int_t largechhksize = largechhk.size();
    Int_t okaychsize = okaych.size();
    
    
    TCanvas * c1 = new TCanvas("c1","c1",250*((nfile+1)/2),500);
    c1->Divide((nfile+1)/2,2);
    c1->Print("Failed_HV_Fit_HK.pdf[");
    for (Int_t ihk =0; ihk < badchhksize; ihk++){
        Int_t c1divide = ihk * (nfile+1) % (8) + 1;
        Int_t c2divide;
        std:: cout << "c1divide: " << c1divide << " ihk: " << ihk+1 << std::endl;
        TGraphErrors * gr1hk = (TGraphErrors*)fin1->Get(Form("HK_PMT_HVscan_Cable_%06d", badchhk[ihk]));
        c1->cd(c1divide);
        gr1hk->Draw("AP");
        c1->Update();
        TH1D * hhk1[nfile];
        for (Int_t i = 0; i < nfile; i++){
            c2divide = c1divide + i + 1;
            hhk1[i] = (TH1D*)fitfin[i]->Get(Form("h_spe_onoff_%d",badchhk[ihk]));
            std:: cout << "c2divide: " << c2divide << std::endl;
            c1->cd(c2divide);
            hhk1[i]->Draw();
            c1->Update();
            if (c2divide < (nfile+1)) {
                std::cout << "Drawing c1 " << c2divide << std::endl;
            }
            else {
                c1->Print("Failed_HV_Fit_HK.pdf");
                c1->Clear();
                c1->Divide((nfile+1)/2,2);
                
            }
        }
        if (ihk == badchhksize - 1 && c2divide != 0){
            c1->Print("Failed_HV_Fit_HK.pdf");
            //c1->Clear();
        }
        //c1count++;
    }
    c1->Modified();
    c1->Print("Failed_HV_Fit_HK.pdf]");
    
    c1->Clear();
    c1->Divide((nfile+1)/2,2);
    c1->Print("Large_HV_Chi2_HK.pdf[");
    for (Int_t ihk =0; ihk < largechhksize; ihk++){
        Int_t c1divide = ihk * (nfile+1) % 8 + 1;
        Int_t c2divide;
        std:: cout << "c1divide: " << c1divide << " ihk: " << ihk+1 << std::endl;
        TGraphErrors * gr1hk = (TGraphErrors*)fin1->Get(Form("HK_PMT_HVscan_Cable_%06d", largechhk[ihk]));
        c1->cd(c1divide);
        gr1hk->Draw("AP");
        c1->Update();
        TH1D * hhk1[nfile];
        for (Int_t i = 0; i < nfile; i++){
            c2divide = c1divide + i + 1;
            hhk1[i] = (TH1D*)fitfin[i]->Get(Form("h_spe_onoff_%d",largechhk[ihk]));
            std:: cout << "c2divide: " << c2divide << std::endl;
            c1->cd(c2divide);
            hhk1[i]->Draw();
            c1->Update();
            if (c2divide < (nfile+1)) {
                std::cout << "Drawing c1 " << c2divide << std::endl;
            }
            else {
                c1->Print("Large_HV_Chi2_HK.pdf");
                c1->Clear();
                c1->Divide((nfile+1)/2,2);
                
            }
        }
        if (ihk == largechhksize - 1 && c2divide != 0){
            c1->Print("Large_HV_Chi2_HK.pdf");
            //c1->Clear();
        }
        //c1count++;
    }
    c1->Modified();
    c1->Print("Large_HV_Chi2_HK.pdf]");
    
    
   
    c1->Clear();
    c1->Divide((nfile+1)/2,2);
    c1->Print("Okay_fit_HK.pdf[");
    for (Int_t iok =0; iok < okaychsize; iok++){
        Int_t c1divide = iok * (nfile+1) % 8 + 1;
        Int_t c2divide = 0;
        std:: cout << "c1divide: " << c1divide << " iok: " << iok+1 << "okchannel: " << okaych[iok] << std::endl;
        TGraphErrors * gr1ok;
        Int_t okchannel = okaych[iok];
        gr1ok = (TGraphErrors*)fin1->Get(Form("HK_PMT_HVscan_Cable_%06d", okchannel));
        c1->cd(c1divide);
        gr1ok->Draw("AP");
        c1->Update();
        TH1D * hsk1[nfile];
        for (Int_t i = 0; i < nfile; i++){
            c2divide = c1divide + i + 1;
            hsk1[i] = (TH1D*)fitfin[i]->Get(Form("h_spe_onoff_%d",okchannel));
            std:: cout << "c2divide: " << c2divide << std::endl;
            c1->cd(c2divide);
            hsk1[i]->Draw();
            c1->Update();
            if (c2divide < (nfile+1)) {
                std::cout << "Drawing c1 " << c2divide << std::endl;
            }
            else {
                c1->Print("Okay_fit_HK.pdf");
                c1->Clear();
                c1->Divide((nfile+1)/2,2);
                
            }
        }
        if (iok == okaychsize - 1 && c2divide != 0){
            c1->Print("Okay_fit_HK.pdf");
        }
        //c1count++;
    }
    c1->Modified();
    c1->Print("Okay_fit_HK.pdf]");
    
    /*
     c1->Clear();
    c1->Print("Failed_HV_Fit_SK_off.pdf[");
    for (Int_t isk =0; isk < badchskoffsize; isk++){
        Int_t c1divide = (isk*5+1) % 25;
        Int_t c2divide = 0;
        std:: cout << "c1divide: " << c1divide << " isk: " << isk+1 << std::endl;
        TGraphErrors * gr2sk = (TGraphErrors*)fin1->Get(Form("SK_PMT_HVscan_Cable_%06d", badchsk_off[isk]));
        c1->cd(c1divide);
        gr2sk->Draw("AP");
        c1->Update();
        
        TH1D * hsk2[4];
        for (Int_t i = 0; i < 4; i++){
            c2divide = c1divide + i + 1;
            hsk2[i] = (TH1D*)fitfin[i]->Get(Form("h_spe_onoff_%d",badchsk_off[isk]));
            std:: cout << "c2divide: " << c2divide << std::endl;
            c1->cd(c2divide);
            hsk2[i]->Draw();
            c1->Update();
            if (c2divide < 25) {
                std::cout << "Drawing c1 " << c2divide << std::endl;
            }
            else {
                c1->Print("Failed_HV_Fit_SK_off.pdf");
                
            }
        }
        if (isk == badchskoffsize - 1 && c2divide != 25){
            c1->Print("Failed_HV_Fit_SK_off.pdf");
            c1->Clear();
        }
        //c1count++;
    }
    c1->Modified();
    c1->Print("Failed_HV_Fit_SK_off.pdf]");
    */
    fin1->Close();
    for (Int_t i = 0; i < 7; i++){
        fitfin[i]->Close();
    }
}
