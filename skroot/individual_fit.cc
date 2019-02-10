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
//#include <TRandom3.h>

using namespace std;

#ifdef __CINT__
void individual_fit(){
#else
int main(int argc, char *argv[]) {
#endif


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
   
    Int_t runno[] = {80282, 80278, 80275, 80269, 80265, 80263, 80254};
    Double_t hvshift[] = {-75, -50, -25, 0, 25, 50, 75};
    Double_t threshold[] = {-0.69, -0.69, -0.69, -0.69, -0.69, -0.69, -0.69};
    TString outdir = "hv_ana/";
    
    //TFile *fin1 = new TFile("hvscan_parameter_sk.root", "Read");
    //TFile *fin2 = new TFile("hvscan_parameter_offset.root", "Read");
    //TFile *fitfin[7];
    //for (Int_t i = 0; i < nfile; i++){
    //fitfin[i] = new TFile(outdir+Form("fit_result_%d_sk.root", runno[i]), "Read");
    //}
    
    const int nfile = 7;
    const int MAXPM = 11146;
    Int_t PMTinfo[MAXPM] = {0};
    TFile *fin1 = new TFile("hvscan_parameter_sk.root", "Read");
    //TFile *fin2 = new TFile("hvscan_parameter_offset.root", "Read");
    TFile *fitfin[nfile];
    for (Int_t i = 0; i < nfile; i++){
      fitfin[i] = new TFile(outdir+Form("fit_result_%d_sk.root", runno[i]), "Read");
    }
    

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
    }
    
    std::vector<Int_t> badchsk;
    std::vector<Int_t> largechsk;
    std::vector<Int_t> zerochsk;
    std::vector<Int_t> okaych;
    
    TTree * trsk2 = (TTree*)fin1->Get("hvscan_sk2");
    TTree * trsk3 = (TTree*)fin1->Get("hvscan_sk3");
    
    Int_t nentry_sk2 = trsk2->GetEntries();
    Int_t nentry_sk3 = trsk3->GetEntries();
    //Int_t nentry_hk = trhk->GetEntries();
    Int_t sk2ch, sk3ch;
    trsk2->SetBranchAddress("Channel", &sk2ch);
    trsk3->SetBranchAddress("Channel", &sk3ch);
    
    for (Int_t itry2 = 0; itry2 < nentry_sk2; itry2++){
        trsk2->GetEntry(itry2);
        okaych.push_back(sk2ch);
    }
    for (Int_t itry3 = 0; itry3 < nentry_sk3; itry3++){
        trsk3->GetEntry(itry3);
        okaych.push_back(sk3ch);
    }
    

    ofstream outxt;
    outxt.open("strangenumbers.txt");
    outxt << setw(8) << " Channel" << setw(6) << "   PMT" << setw(5) << "  lhv" << setw(5) << "  hhv" << setw(8) << "    norm" << setw(8) << "   index" << setw(8) << "    chi2" << setw(8) << "    prob" << "\n";
      
    ifstream intxt1;
    intxt1.open("badfitting_sk.txt");
    if (intxt1.fail()){
        std::cerr << "Error opening file intxt1" << std::endl;
        exit(1);
    }
    getline(intxt1, line);
   
    if(intxt1.is_open()){
        while(!intxt1.eof()){
            Int_t channel;
            std::string pmt, norm, beta;
	    std::string chi2, prob;
            std::string lhv, hhv, comment;
            
            intxt1 >> channel >> pmt >> lhv >> hhv >> norm >> beta >> chi2 >> prob >> comment;
            
            std::cout << channel << " " << pmt << std::endl;
            okaych.erase(std::remove(okaych.begin(), okaych.end(), channel), okaych.end());
	    if ((lhv == "-nan") || (hhv == "-nan") || (chi2 == "inf"))
	      outxt << setw(8) << channel << setw(6) << pmt.c_str() << setw(5) << lhv.c_str() << setw(5) << hhv.c_str() << setw(8) << norm.c_str() << setw(8) << beta.c_str() << setw(8) << chi2.c_str() << setw(8) << prob.c_str() << "\n";
	    
	    
            if ((pmt == "SK2" || pmt == "SK3") && comment == "Failed_to_find_minimum"){
                badchsk.push_back(channel);
                std::cout << "Add 1 badch for SK PMT" << std::endl;
            }
            else if ((pmt == "SK2" || pmt == "SK3") && comment == "Chi2_>>_1"){
                largechsk.push_back(channel);
                std::cout << "Add 1 badch for SK PMT" << std::endl;
            }
            else if ((pmt == "SK2" || pmt == "SK3") && comment == "Chi2_=_0"){
                zerochsk.push_back(channel);
		std::cout << "Add 1 badch for SK PMT" << std::endl;
	    }
        }
        intxt1.close();
    }
    outxt.close();
    
    //intxt1.close();
    Int_t badchsksize = badchsk.size();
    Int_t largechsksize = largechsk.size();
    Int_t zerochsksize = zerochsk.size();
    Int_t okaychsize = okaych.size();
    cout << "In total bad channel: " << badchsksize << " okay channel: " << okaychsize << " largechi2 channel: " << largechsksize << endl << endl;
    
    TCanvas * c1 = new TCanvas("c1","c1",250*((nfile+1)/2),500);
    c1->Divide((nfile+1)/2,2);
    c1->Print("Failed_HV_Fit_SK.pdf(");
    for (Int_t isk =0; isk < badchsksize; isk++){
      //if (badchsk[isk] == 2740 || badchsk[isk] == 5543 || badchsk[isk] == 6066) continue;
        Int_t c1divide = 1;
        Int_t c2divide = 0;
        std:: cout << " isk: " << isk << " bad cable " << badchsk[isk] << " total bad cable# " << badchsksize <<  std::endl;
        TGraphErrors * gr1sk = (TGraphErrors*)fin1->Get(Form("SK%d_PMT_HVscan_Cable_%06d", PMTinfo[badchsk[isk]-1]-1, badchsk[isk]));
        c1->cd(c1divide);
        gr1sk->Draw("AP");
        c1->Update();
        TH1D * hsk1[nfile];
        for (Int_t i = 0; i < nfile; i++){
            c2divide = c1divide + i + 1;
	    //std:: cout << "c2divide: " << c2divide << " file: " << runno[i] << " badcable " << badchsk[isk] << std::endl;
	    //if (!fitfin[i]->GetListOfKeys()->Contains(Form("h_spe_onoff_%d",badchsk[isk]))) continue;
	    //if (c2divide < (nfile+1)) {
	    //std::cout << "Drawing c1 " << c2divide << " nfile+1: " << nfile+1 << std::endl;
	    //}
	    hsk1[i] = (TH1D*)fitfin[i]->Get(Form("h_spe_onoff_%d",badchsk[isk]));
	    //if (hsk1[i] == (TH1D*)0) continue;
            //std:: cout << "c2divide: " << c2divide << " file: " << fitfin[i]->Data() << std::endl;
            c1->cd(c2divide);
            hsk1[i]->Draw();
            c1->Update();
            //if (c2divide < (nfile+1)) {
	    //std::cout << "Drawing c1 " << c2divide << std::endl;
            //}
            if (c2divide == (nfile+1)) {
                c1->Print("Failed_HV_Fit_SK.pdf");
                c1->Clear();
                c1->Divide((nfile+1)/2,2);
                
            }
        }
        //if (isk == badchsksize - 1 && c2divide != 0){
	//c1->Print("Failed_HV_Fit_SK.pdf");
        //}
        //c1count++;
    }
    c1->Modified();
    c1->Print("Failed_HV_Fit_SK.pdf)");
    
    c1->Clear();
    c1->Divide((nfile+1)/2,2);
    c1->Print("Large_HV_Chi2_SK.pdf(");
    for (Int_t isk =0; isk < largechsksize; isk++){
      if (largechsk[isk] == 2740 || largechsk[isk] == 5543 || largechsk[isk] == 6066) continue;
        Int_t c1divide = isk * (nfile+1) % 8 + 1;
        Int_t c2divide = 0;
	std:: cout << " isk: " << isk << " lchi2 cable " << largechsk[isk] << " total lchi2 cable# " << largechsksize <<  std::endl;
	//std:: cout << "c1divide: " << c1divide << " isk: " << isk+1 << std::endl;
        TGraphErrors * gr1sk = (TGraphErrors*)fin1->Get(Form("SK%d_PMT_HVscan_Cable_%06d", PMTinfo[largechsk[isk]-1]-1, largechsk[isk]));
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
            //if (c2divide < (nfile+1)) {
	      //std::cout << "Drawing c1 " << c2divide << std::endl;
	      //}
            if (c2divide == nfile+1) {
                c1->Print("Large_HV_Chi2_SK.pdf");
                c1->Clear();
                c1->Divide((nfile+1)/2,2);
                
            }
        }
        //if (isk == largechsksize - 1 && c2divide != 0){
	//c1->Print("Large_HV_Chi2_SK.pdf");
        //}
        //c1count++;
    }
    c1->Modified();
    c1->Print("Large_HV_Chi2_SK.pdf)");
    
    /*c1->Clear();
    c1->Divide((nfile+1)/2,2);
    c1->Print("Zero_HV_Chi2_SK.pdf(");
    for (Int_t isk =0; isk < zerochsksize; isk++){
        Int_t c1divide = isk * (nfile+1) % 8 + 1;
        Int_t c2divide = 0;
        std:: cout << "c1divide: " << c1divide << " isk: " << isk+1 << std::endl;
        TGraphErrors * gr1sk = (TGraphErrors*)fin1->Get(Form("SK%d_PMT_HVscan_Cable_%06d", PMTinfo[zerochsk[isk]-1]-1, zerochsk[isk]));
        c1->cd(c1divide);
        gr1sk->Draw("AP");
        c1->Update();
        TH1D * hsk1[nfile];
        for (Int_t i = 0; i < nfile; i++){
            c2divide = c1divide + i + 1;
            hsk1[i] = (TH1D*)fitfin[i]->Get(Form("h_spe_onoff_%d",largechsk[isk]));
            std:: cout << "c2divide: " << c2divide << std::endl;
            c1->cd(c2divide);
            hsk1[i]->Draw();
            c1->Update();
            //if (c2divide < (nfile+1)) {
                //std::cout << "Drawing c1 " << c2divide << std::endl;
            //}
            if (c2divide == nfile+1)  {
                c1->Print("Zero_HV_Chi2_SK.pdf");
                c1->Clear();
                c1->Divide((nfile+1)/2,2);
                
            }
        }
        //if (isk == zerochsksize - 1 && c2divide != 0){
	//  c1->Print("Zero_HV_Chi2_SK.pdf");
        //}
        //c1count++;
    }
    c1->Modified();
    c1->Print("Zero_HV_Chi2_SK.pdf)");*/
    
    c1->Clear();
    c1->Divide((nfile+1)/2,2);
    c1->Print("Okay_fit.pdf(");
    for (Int_t iok =0; iok < okaychsize; iok++){
        Int_t c1divide = iok * (nfile+1) % 8 + 1;
        Int_t c2divide = 0;
        std:: cout  << " iok: " << iok+1 << "okchannel: " << okaych[iok] << " in total okay# " << okaychsize << std::endl;
        Int_t okchannel = okaych[iok];
        TGraphErrors * gr1ok = (TGraphErrors*)fin1->Get(Form("SK%d_PMT_HVscan_Cable_%06d", PMTinfo[okchannel-1]-1, okchannel));;
        c1->cd(c1divide);
        gr1ok->Draw("AP");
        c1->Update();
        TH1D * hsk1[nfile];
        for (Int_t i = 0; i < nfile; i++){
            c2divide = c1divide + i + 1;
            hsk1[i] = (TH1D*)fitfin[i]->Get(Form("h_spe_onoff_%d",okchannel));
	    //            std:: cout << "c2divide: " << c2divide << std::endl;
            c1->cd(c2divide);
            hsk1[i]->Draw();
            c1->Update();
            //if (c2divide < 25) {
	    //std::cout << "Drawing c1 " << c2divide << std::endl;
            //}
            if (c2divide == nfile+1) {
                c1->Print("Okay_fit.pdf");
                c1->Clear();
                c1->Divide((nfile+1)/2,2);
                
            }
        }
        //if (iok == largechsksize - 1 && c2divide != 0){
        //    c1->Print("Okay_fit.pdf");
        //}
        //c1count++;
    }
    c1->Modified();
    c1->Print("Okay_fit.pdf)");
    

    
    fin1->Close();

    for (Int_t i = 0; i < nfile; i++){
        fitfin[i]->Close();
    }
}
