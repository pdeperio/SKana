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
#include <TLatex.h>
#include <TPaveStats.h>
#include <TList.h>
#include <TText.h>
#include <TSpectrum.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TString.h>

using namespace std;

#ifndef __CINT__
int getArgs(int argc, char* argv[]);
int runno = -1;
TString ConnectionFile = "";
Float_t hv_shift = 0;
int SelectPMT = -1;
TString PMTtype = "";
TString InputDir = "output";
#endif

TString outdir = "hv_ana";

TF1* FitHK(TH1D* h, Int_t rebin = 2, Float_t hv = 0, Float_t PMTflag = 6){
    Int_t npx = 1000;
    Int_t entry = h->GetEntries();
    Int_t nbin = h->GetNbinsX();
    if (rebin > 1) h->Rebin(rebin);
    
    
    TF1 *func1pe = new TF1("fgaus","gaus",0, 12);
    TF1 *func1pe_new = new TF1("fgaus_new", "gaus",0, 12);
    TF1 *func2peak_bs = new TF1("f2peakbs","gaus(0)+0.5*[3]*[0]*((TMath::Erf((x/[2])))+(1-TMath::Erf((x-[1]+0.5*[2])/[2]))-1)+0.25*[5]*[0]*((TMath::Erf((x-[1]-0.5*[2])/[2])+(1-TMath::Erf((x-[4]*[1])/[2])))-1)",0,12);
    //TF1 *func0peak = new TF1("fgaus0","gaus",-1, 12);
    TF1 *func1peak = new TF1("fgaus1","gaus",0, 12);
    TF1 *funcbs = new TF1("fbs","0.5*[3]*[0]*((TMath::Erf((x/[2])))+(1-TMath::Erf((x-[1]+0.5*[2])/[2]))-1)+0.25*[5]*[0]*((TMath::Erf((x-[1]-0.5*[2])/[2])+(1-TMath::Erf((x-[4]*[1])/[2])))-1)",0, 12);
    
    
    func1pe->SetLineColor(2);
    func1pe_new->SetLineColor(46);
    func2peak_bs->SetLineColor(38);
    //func0peak->SetLineColor(30);
    func1peak->SetLineColor(9);
    funcbs->SetLineColor(7);
    
    func1pe->SetLineStyle(3);
    //func0peak->SetLineStyle(3);
    func1peak->SetLineStyle(3);
    funcbs->SetLineStyle(3);
    
    func1pe->SetNpx(npx);
    func1pe_new->SetNpx(npx);
    func2peak_bs->SetNpx(npx);
    func1peak->SetNpx(npx);
    //func0peak->SetNpx(npx);
    funcbs->SetNpx(npx);
    
    Double_t peakx_pre[2] = {2.9, 1.5};
    Double_t peaky_pre[2] = {h->GetBinContent(h->FindBin(peakx_pre[0])), h->GetBinContent(h->FindBin(peakx_pre[1]))};
    
    //preliminarily fit the main peak
    func1pe->SetParameters(h->GetMaximum(), peakx_pre[0], 1);
    func1pe->SetParLimits(0, 0, h->GetMaximum());
    func1pe->SetParLimits(1, 1.5, 10);
    func1pe->SetParLimits(2, 0, 3);
    h->Fit(func1pe, "BNQ0", "", 1.5, 10);
    
    //refit the main peak with tuned range, set to the fit of the hv=-100 data
    Double_t func1pe_low = func1pe->GetParameter(1) - func1pe->GetParameter(2);
    Double_t func1pe_high = func1pe->GetParameter(1) + 2*func1pe->GetParameter(2);
    for (Int_t ipar = 0; ipar < 3; ipar++){
        func1pe_new->SetParameter(ipar, func1pe->GetParameter(ipar));
    }
    h->Fit(func1pe_new, (hv==-100)?"BQ":"BNQ0", "", func1pe_low, func1pe_high);
    
    if (hv >= -100){//fit 2peak+erf
        for(int iSet=0; iSet<2; iSet++){
            func2peak_bs->SetParName(iSet,Form("Scale_{%d}",iSet));
            func2peak_bs->SetParName(iSet*3+1,Form("Peak_{%d}",iSet));
            func2peak_bs->SetParName(iSet*3+2,Form("#sigma_{%d}",iSet));
        }
        
        for(int iSet=0; iSet<3;iSet++){
            func2peak_bs->SetParameter(iSet, func1pe_new->GetParameter(iSet));
            func2peak_bs->SetParLimits(iSet, func1pe_new->GetParameter(iSet)-func1pe_new->GetParError(iSet), func1pe_new->GetParameter(iSet)+func1pe_new->GetParError(iSet));
            //func1peak->SetParameter(iSet, func1pe_new->GetParameter(iSet));
        }
        
        //func2peak_bs->SetParameter(3, func1pe_new->GetParameter(0)*0.25);
        //func2peak_bs->SetParameter(4, func1pe_new->GetParameter(1)-1.5*func1pe_new->GetParameter(2));
        //func2peak_bs->SetParameter(5, 0.7*func1pe_new->GetParameter(2));
        //func2peak_bs->SetParLimits(3, 0, 0.4*func1pe_new->GetParameter(0));
        //func2peak_bs->SetParLimits(4, 0, func1pe_new->GetParameter(1)-func1pe_new->GetParameter(2));
        //func2peak_bs->SetParLimits(5, 0.3*func1pe_new->GetParameter(2), func1pe_new->GetParameter(2));
        func2peak_bs->SetParName(3, "1pe BG ratio");
        func2peak_bs->SetParName(5, "2pe BG ratio");
        func2peak_bs->SetParName(4, "2pe peak region");
        func2peak_bs->SetParameter(3, 0.5);
        func2peak_bs->SetParLimits(3, 0, 1);
        func2peak_bs->SetParameter(4, 2);
        func2peak_bs->SetParLimits(4, 1, 3);
        func2peak_bs->SetParameter(5, 0.2);
        func2peak_bs->SetParLimits(5, 0, func2peak_bs->GetParameter(3));
        
       
        h->Fit(func2peak_bs,"BQ+","",0, 12);
        
        for(int iSet = 0; iSet < 3; iSet++){
            func1peak->FixParameter(iSet, func2peak_bs->GetParameter(iSet));
            funcbs->FixParameter(iSet, func2peak_bs->GetParameter(iSet));
            funcbs->FixParameter(iSet+3, func2peak_bs->GetParameter(iSet+3));
        }
        //funcbs->FixParameter(6, func2peak_bs->GetParameter(6));
        
        //h->Fit(func0peak,"BQ+","", -1, 10);
        h->Fit(func1peak,"BQ+","", 0, 10);
        h->Fit(funcbs,"BQ+","", 0, 10);
        
    }
    
    TF1 *func;
    if (hv >= -100) func = func2peak_bs;
    else func = func1pe_new;
    
    gPad->cd();
    func1pe->Delete();
    func1pe_new->Delete();
    
    h->Draw();
    gPad->Update();
    
    Float_t st0_lower_left_x = 0.5;
    Float_t st0_lower_left_y = 0.5;
    Float_t st_Width = 0.4;
    Float_t st_Height = 0.4;
    
    TPaveStats *st0 = (TPaveStats*)gPad->GetPrimitive("stats");
    st0->SetName("result");
    TList *listOfLines = st0->GetListOfLines();
    if (st0) {
        st0->SetX1NDC(st0_lower_left_x);
        st0->SetY1NDC(st0_lower_left_y);
        st0->SetX2NDC(st0_lower_left_x + st_Width);
        st0->SetY2NDC(st0_lower_left_y + st_Height);
        TLatex *myt[6];
        myt[0] = new TLatex(0,0,Form("Single PE RESULT"));
        //if (func->GetNpar()>6) myt[1] = new TLatex(0,0,Form("SPE Peak  = %1.2e [pC]",func->GetParameter(1)>func->GetParameter(4)?func->GetParameter(1):func->GetParameter(4)));
        myt[1] = new TLatex(0,0,Form("SPE Peak  = %1.2e [pC]",func->GetParameter(1)));
        myt[2] = new TLatex(0,0,Form("SPE Peak #sigma = %1.2e [pC]",func->GetParameter(2)));
        myt[3] = new TLatex(0,0,Form("Func Max = %1.2e [pC]",func->GetMaximumX(0,10)));
        myt[4] = new TLatex(0,0,Form("Chi2ndf  = %1.2e",func->GetChisquare()/func->GetNDF()));
     
        for (Int_t i = 0; i < 5; i++) {
            myt[i] ->SetTextFont(42);
            myt[i] ->SetTextSize(0.03);
            myt[i] ->SetTextColor(kBlack);
            listOfLines->Add(myt[i]);
        }
        //if (func->GetNpar() > 6){
            //myt[5] ->SetTextFont(42);
            //myt[5] ->SetTextSize(0.03);
            //myt[5] ->SetTextColor(kBlack);
            //listOfLines->Add(myt[4]);
        //}
        h->SetStats(0);
        
        gPad->Modified();
    }
    
    /*if (func->GetParameter(1)<func->GetParameter(4)){
        Double_t temp, temperr, tempamp, tempsig;
        tempamp = func->GetParameter(0);
        temp = func->GetParameter(1);
        tempsig = func->GetParameter(2);
        temperr = func->GetParError(1);
        func->SetParameter(0, func->GetParameter(3));
        func->SetParameter(1, func->GetParameter(4));
        func->SetParameter(2, func->GetParameter(5));
        func->SetParError(1, func->GetParError(4));
        func->SetParameter(3, tempamp);
        func->SetParameter(4, temp);
        func->SetParameter(5, tempsig);
        func->SetParError(4, temperr);
    }*/
    return func;
}

TF1*/*std::vector<Double_t>*/ FitSK(TH1D* h, Int_t rebin = 4, Float_t hv = 0, Float_t PMTflag = 3){
    Int_t npx = 1000;
    Int_t entry = h->GetEntries();
    Int_t nbin = h->GetNbinsX();
    if (rebin > 1) h->Rebin(rebin);
    TF1 *func = 0;
    //Int_t binmx = h->GetMaximumBin();
    //Double_t xmax = h->FindBin(binmx);

    double fit_range[2] = {0, 12+hv*0.01};
    double hv_scale = 0.01;
    
    TF1 *func1pe = new TF1("fgaus","gaus", fit_range[0], fit_range[1]);
    TF1 *func1pe_new = new TF1("fgaus_new", "gaus", fit_range[0], fit_range[1]);
    TF1 *func2peak = new TF1("func2peak", "gaus(0)+[3]*(TMath::Exp(-(x-[4])/[5]))", fit_range[0], fit_range[1]);
    TF1 *func1peak = new TF1("fgaus1","gaus", fit_range[0], fit_range[1]);
    TF1 *funcexbs = new TF1("fexbs","[3]*(TMath::Exp(-(x-[4])/[5]))", fit_range[0], fit_range[1]);
    TF1 *funcskgaus = new TF1("funcskgaus","[3]*gaus(0)*(1+[4]*(TMath::Erf(((x-[1])/[2]))))", fit_range[0], fit_range[1]);
    
    
    func1pe->SetLineColor(2);
    func1pe_new->SetLineColor(46);
    func2peak->SetLineColor(38);
    funcskgaus->SetLineColor(8);
    func1peak->SetLineColor(9);
    funcexbs->SetLineColor(7);
    
    func1pe->SetLineStyle(3);
    func1peak->SetLineStyle(3);
    funcexbs->SetLineStyle(3);
    
    func1pe->SetNpx(npx);
    func1pe_new->SetNpx(npx);
    func2peak->SetNpx(npx);
    func1peak->SetNpx(npx);
    funcskgaus->SetNpx(npx);
    funcexbs->SetNpx(npx);
    
    Double_t peakx_pre[2] = {3+hv*hv_scale, 1.5+hv*hv_scale};
    Double_t peaky_pre[2] = {h->GetBinContent(h->FindBin(peakx_pre[0])), h->GetBinContent(h->FindBin(peakx_pre[1]))};
    
    //preliminarily fit the main peak
    func1pe->SetParameters(h->GetMaximum(), peakx_pre[0], 2);
    func1pe->SetParLimits(0, 0, h->GetMaximum());
    func1pe->SetParLimits(1, 1+hv*hv_scale, 10+hv*hv_scale);
    func1pe->SetParLimits(2, 1, 5);
    h->Fit(func1pe, "BNQ0", "", 1+hv*hv_scale, 10+hv*hv_scale);
    
    //refit the main peak with tuned range
    Double_t func1pe_low = func1pe->GetParameter(1) - func1pe->GetParameter(2);
    Double_t func1pe_high = func1pe->GetParameter(1) + 2*func1pe->GetParameter(2);
    for (Int_t ipar = 0; ipar < 3; ipar++){
        func1pe_new->SetParameter(ipar, func1pe->GetParameter(ipar));
        func1pe_new->SetParLimits(1, func1pe->GetParameter(1) - 0.5*func1pe->GetParameter(2), func1pe->GetParameter(1) + 0.5*func1pe->GetParameter(2));
        func1pe->SetParLimits(2, 1, 5);
    }
    h->Fit(func1pe_new, (func1pe_new->GetParameter(2)<0.3&&func1pe_new->GetParameter(1)<1&&func1pe_new->GetParameter(2)/func1pe_new->GetParameter(1)>0.01)?"BQ+":"BNQ0", "", func1pe_low, func1pe_high);
    //h->Fit(func1pe_new, (hv < -70)?"BQ":"BNQ0", "", func1pe_low, func1pe_high);
    func = func1pe_new;
    
    if (func1pe_new->GetParameter(2)>0.3&&func1pe_new->GetParameter(1)>1){//fit gaus + expo
        func2peak->SetParName(0,"Scale");
        func2peak->SetParName(1,"Peak");
        func2peak->SetParName(2,"#sigma");
        
        for(int iSet=0; iSet<3;iSet++){
	  func2peak->SetParameter(iSet, func1pe_new->GetParameter(iSet));
	  func2peak->SetParLimits(iSet, func1pe_new->GetParameter(iSet)-func1pe_new->GetParError(iSet), func1pe_new->GetParameter(iSet)+func1pe_new->GetParError(iSet));
            //func1peak->FixParameter(iSet, func1pe_new->GetParameter(iSet));
        }
        
        func2peak->SetParameter(3, func1pe_new->GetParameter(0));
        func2peak->SetParameter(4, func1pe_new->GetParameter(1)-2*func1pe_new->GetParameter(2)>0?func1pe_new->GetParameter(1)-2*func1pe_new->GetParameter(2):0);
        func2peak->SetParameter(5, func1pe_new->GetParameter(2));
        //func2peak->SetParameter(5, 0.5);
        //func2peak->SetParameter(6, func1pe_new->GetParameter(1)-func1pe_new->GetParameter(2));
        func2peak->SetParLimits(3, 0, 1.5 * func1pe_new->GetParameter(0));
        func2peak->SetParLimits(4, 0, func1pe_new->GetParameter(1));
        //func2peak->SetParLimits(5, 0, func1pe_new->GetParameter(2));
        func2peak->SetParLimits(5, 0, 2*func1pe_new->GetParameter(2));
        //func2peak->SetParLimits(6, func1pe_new->GetParameter(1)-1.5*func1pe_new->GetParameter(2), func1pe_new->GetParameter(1));
        
        h->Fit(func2peak,"BQ+","",0, 12+hv*hv_scale);
        //h->Fit(func2peak_bs,"BQ","",func1pe_new->GetParameter(1) - func1pe_new->GetParameter(2), func1pe_new->GetParameter(1) + 2 * func1pe_new->GetParameter(2));
        /*for(int iSet = 3; iSet < 5; iSet++){
         funcexpo->FixParameter(iSet, func2peak->GetParameter(iSet));
         }
         
         for(int iSet = 5; iSet < 7; iSet++){
         funcpos->FixParameter(iSet, func2peak->GetParameter(iSet));
         }*/
        
        func = func2peak;
        for(int iSet = 0; iSet < 3; iSet++){
            func1peak->FixParameter(iSet, func2peak->GetParameter(iSet));
            //funcexbs->FixParameter(iSet, func2peak->GetParameter(iSet));
            funcexbs->FixParameter(iSet+3, func2peak->GetParameter(iSet+3));
        }
        
        //h->Fit(funcexpo, "BQ+","", 0+hv*0.008, 10+hv*0.008);
        h->Fit(func1peak,"BQ+","", 0, 10+hv*hv_scale);
        h->Fit(funcexbs,"BQ+","", 0, 10+hv*hv_scale);
        //h->Fit(funcpos,"BQ+","", 0+hv*0.008, 10+hv*0.008);

        funcskgaus->Delete();
	func1pe_new->Delete();
    }
    
    else if (func1pe_new->GetParameter(2)/func1pe_new->GetParameter(1) < 0.01){//fit assym gaus
        funcskgaus->SetParameter(0, h->GetMaximum());
        funcskgaus->SetParameter(1, 1);
        funcskgaus->SetParameter(2, 1);
        funcskgaus->SetParameter(3, 1);
        funcskgaus->SetParameter(4, 5);
        funcskgaus->SetParLimits(0, 0.5*h->GetMaximum(), 1.5*h->GetMaximum());
        funcskgaus->SetParLimits(1, 0, 2);
        funcskgaus->SetParLimits(2, 0.5, 4);
        funcskgaus->SetParLimits(3, 0.5, 2);
        funcskgaus->SetParLimits(4, 0, 10);
        
        h->Fit(funcskgaus, "BV+","", 0, 12+hv*hv_scale);
        func = funcskgaus;

	func2peak->Delete();
	func1pe_new->Delete();
    }

    if (!func) {
      cout << "Error: func not set" << endl;
      exit(-2);
    }    
    
    //if (hv >= -100) func = func2peak;
    //else func = func1pe_new;
    
    gPad->cd();
    func1pe->Delete();
    func1peak->Delete();
    funcexbs->Delete();
    
    h->Draw();
    gPad->Update();
    
    Float_t st0_lower_left_x = 0.5;
    Float_t st0_lower_left_y = 0.55;
    Float_t st_Width = 0.4;
    Float_t st_Height = 0.35;
    
    TPaveStats *st0 = (TPaveStats*)gPad->GetPrimitive("stats");
    st0->SetName("result");
    TList *listOfLines = st0->GetListOfLines();
    if (st0) {
        st0->SetX1NDC(st0_lower_left_x);
        st0->SetY1NDC(st0_lower_left_y);
        st0->SetX2NDC(st0_lower_left_x + st_Width);
        st0->SetY2NDC(st0_lower_left_y + st_Height);
        TLatex *myt[6];
        myt[0] = new TLatex(0,0,Form("Single PE RESULT"));
        //myt[1] = new TLatex(0,0,Form("# of Hits = %d", h->GetEntries()));
        myt[1] = new TLatex(0,0,Form("SPE Peak  = %1.2e [pC]",func->GetParameter(1)>func->GetMaximumX()-1?func->GetParameter(1):func->GetMaximumX()));
        myt[2] = new TLatex(0,0,Form("SPE Peak #sigma = %1.2e [pC]",func->GetParameter(2)));
        myt[3] = new TLatex(0,0,Form("Func Max = %1.2e [pC]",func->GetMaximumX(0,10)));
        myt[4] = new TLatex(0,0,Form("Chi2ndf  = %1.2e",func->GetChisquare()/func->GetNDF()));
        //if (func->GetNpar() > 3){
            //myt[4] = new TLatex(0,0,Form("Poisson Peak = %1.2e [pC]",func->GetParameter(6)));
            //myt[5] = new TLatex(0,0,Form("Poisson Peak 2 = %1.2e [pC]",func->GetParameter(6)));
        //}
        if (func->GetNpar()>5){
            myt[5] = new TLatex(0,0,Form("Decay = %1.2e",func->GetParameter(5)));
        }
        else myt[5] = new TLatex(0,0,Form("Skew = %1.2e",func->GetParameter(4)));
        //myt[6] = new TLatex(0,0,Form("FWHM = %3.1f [%%]", (result.peakx!=0?result.FWHM/result.peakx*100.:-1)));
        //myt[6] = new TLatex(0,0,Form("Peak (1pe only) = %2.2f [pC]",    result.peak1pex));
        //myt[7] = new TLatex(0,0,Form("Sigma (1pe only) = %2.1f [%%]",     result.peak1pesigma));
        //myt[8] = new TLatex(0,0,Form("Gain (1pe only) = %1.3e",  result.gainpeak1pex));
        for (Int_t i = 0; i < 6; i++) {
            myt[i] ->SetTextFont(42);
            myt[i] ->SetTextSize(0.03);
            myt[i] ->SetTextColor(kBlack);
            listOfLines->Add(myt[i]);
        }
        
        h->SetStats(0);
        
        gPad->Modified();
    }
    return func;
}

#ifdef __CINT__
void fit_hvscan_spe(int runno = -1, TString ConnectionFile = "", Float_t hv_shift=0, TString PMTtype = "", TString InputDir = "output", int SelectPMT=-1){
    
#else
int main(int argc, char *argv[]) {// process the arguments
    int args = getArgs(argc, argv);
    if(args != 0){
        std::cerr << "Usage " << std::endl;
        return 0;
    }
#endif

    outdir += PMTtype + "/";

    gStyle->SetFrameBorderMode(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetPadColor(0);

    gStyle->SetOptTitle(kTRUE);
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(109);
    
    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    
    const int MAXPM = 11146;
    const int nPMTtype = 3;
    
    Float_t PMTinfo[MAXPM][3];//flag, place, hv
    
    ifstream connect;
    connect.open(ConnectionFile.Data());
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
            
            PMTinfo[cableid-1][0] = pmtflag;
            PMTinfo[cableid-1][1] = pmtz;
            PMTinfo[cableid-1][2] = oldhv;
        }
        connect.close();
    }
    
    else {
        cout << "Error: Connection file (" << ConnectionFile.Data() << ") not open" << endl;
        exit (-1);
    }
    
    TFile *f;
    TFile *fout;
    TH1D * h_on[MAXPM];
    TH1D * h_off[MAXPM];
    
    int negbad = 0;
    int smallbad = 0;
    int largebad =0;
    int posbad = 0;
    int nofit = 0;//counters
        
    Int_t chid;
    Double_t peak;
    Double_t peakerr;
    //Double_t perr;
    Double_t chi2;
    Double_t ndf;
    //Double_t pvalue;
    Double_t highv;
    Double_t rchi2;
    Double_t sigma;
    Int_t nhit;
        
    TTree *tree = new TTree("spe", "SPE Peaks");
    tree->Branch("Channel", &chid, "chid/I");
    tree->Branch("HV", &highv, "highv/D");
    tree->Branch("Peak", &peak, "peak/D");
    tree->Branch("Peakerr", &peakerr, "peakerr/D");
    tree->Branch("Chi2", &rchi2, "rchi2/D");
    tree->Branch("Sigma", &sigma, "sigma/D");
    tree->Branch("Nhits", &nhit, "nhit/I");
        
    ofstream outfile;
    outfile.open(outdir+Form("badpeak_%d%s.txt", runno, PMTtype.Data()));
    outfile << setw(10) << "     Cable" << setw(10) << "    HV [V]" << setw(10) << " Peak [pC]" << setw(20) << "         Err Message\n";
        
    ofstream outfile2;
    outfile2.open(outdir+Form("nofit_%d%s.txt", runno, PMTtype.Data()));

    TString neg_bad_name = outdir+Form("neg_bad_fit_%d%s.pdf", runno, PMTtype.Data());
    TString small_bad_name = outdir+Form("small_bad_fit_%d%s.pdf", runno, PMTtype.Data());
    TString big_bad_name = outdir+Form("big_bad_fit_%d%s.pdf", runno, PMTtype.Data());
    TString pos_bad_name = outdir+Form("pos_bad_fit_%d%s.pdf", runno, PMTtype.Data());
    TString okay_fit_name = outdir+Form("okay_fit_sample_%d%s.pdf[", runno, PMTtype.Data());

    c1->Print(neg_bad_name+"[");
    c1->Print(small_bad_name+"[");
    c1->Print(big_bad_name+"[");
    c1->Print(pos_bad_name+"[");
    c1->Print(okay_fit_name+"[");
    
    TH1D * h_spe_peak = new TH1D("h_spe_peak","spe_peak",120,0,6);
    TH1D * h_spe_chi2 = new TH1D("h_spe_chi2","spe_chi2",200,0,40);
    h_spe_peak->GetXaxis()->SetTitle("Peak [pC]");
    h_spe_chi2->GetXaxis()->SetTitle("Fit Chi2/ndf");
    
    TString infile = Form("%s/spe_individual_%d_C.root", InputDir.Data(), runno, PMTtype.Data());
    f = new TFile(infile,"read");
    fout = new TFile(outdir+Form("fit_result_%d%s.root",runno, PMTtype.Data()),"recreate");
    fout->cd();
    
    cout << "Opening: " << infile.Data() << endl;
    
    if (!f->IsOpen()) {
        cout << "Error opening file" << endl;
        exit (-1);
    }
    
    cout << "Inputted HV shift = " << hv_shift << endl;

    int PMTRange[2] = {0, MAXPM};
    if (SelectPMT>=0) {
      PMTRange[0] = SelectPMT;
      PMTRange[1] = SelectPMT+1;
    }

    bool AnalyzeHK = PMTtype.Contains("hk");
    
    for (Int_t iPMT = PMTRange[0]; iPMT < PMTRange[1]; iPMT++){

        // Select only HK or SK PMTs 
        if (AnalyzeHK) {

	  // Flag 6: HK PMT
	  if (PMTinfo[iPMT][0] != 6) continue;
	}
       
	else {

	  // Flag 3: SK-2 PMT, Flag 4: SK-3 PMT, 
	  if (!(PMTinfo[iPMT][0] == 3 || PMTinfo[iPMT][0] == 4))
	    continue;
        }
	
        h_on[iPMT] = (TH1D*)f->Get(Form("h_spe_on_%d",iPMT+1));
        h_off[iPMT] = (TH1D*)f->Get(Form("h_spe_off_%d",iPMT+1));
        
        h_on[iPMT]->Sumw2();
        h_off[iPMT]->Sumw2();
        
        h_on[iPMT]->Add(h_off[iPMT],-1);
        
        h_on[iPMT]->SetTitle(Form("PMT %d",iPMT+1));
        h_on[iPMT]->SetName(Form("h_spe_onoff_%d", iPMT+1));
        h_on[iPMT]->GetXaxis()->SetTitle("Charge [pC]");
        
        h_on[iPMT]->GetXaxis()->SetRangeUser(-3,15);
        h_on[iPMT]->Draw();
	
        // Check for dead PMTs
        if (h_on[iPMT]->Integral() < 20){
            nofit++;
            outfile2 << iPMT+1 << std::endl;
            continue;
        }
        
	highv = PMTinfo[iPMT][2];
               
	TF1 *fge;
	if (AnalyzeHK) fge = FitHK((TH1D*)h_on[iPMT], 2, hv_shift, PMTinfo[iPMT][0]);
	else fge = FitSK((TH1D*)h_on[iPMT], 4, hv_shift, PMTinfo[iPMT][0]);
            
	c1->Update();
            
	h_on[iPMT]->Write();
                
	// Manual peak shifting in to fix strange cases

	// For HK
	if (AnalyzeHK && (fge->GetParameter(1) < 0 || fge->GetParameter(1) > 10)) {
	  peak = fge->GetMaximumX(0,6);
	  Double_t FWHMlow  = peak- fge->GetX(fge->Eval(peak)*0.5, 0, peak);
	  Double_t FWHMhigh = fge->GetX(fge->Eval(peak)*0.5, peak, 12) - peak;
	  sigma = (FWHMlow+FWHMhigh)/(2.*TMath::Sqrt(TMath::Log(2)*2));
	  peakerr = sigma * 1e-4;
	}

	// For SK
	else if (!AnalyzeHK && (fge->GetParameter(1) < 0 || fge->GetParameter(1) < fge->GetMaximumX() - 1)) {
	  peak = fge->GetMaximumX(0, 8);
	  //Double_t FWHMlow  = peak- fge->GetX(fge->Eval(peak)*0.5, 0, peak);
	  Double_t FWHMhigh = fge->GetX(fge->Eval(peak)*0.5, peak, 12) - peak;
	  sigma = 2*FWHMhigh/(2.*TMath::Sqrt(TMath::Log(2)*2));
	  peakerr = 0.01 * sigma;
	}

	// Good fits
	else {
	  peak = fge->GetParameter(1);
	  peakerr = fge->GetParError(1);
	  sigma = fge->GetParameter(2);
	}
        
	chi2 = fge->GetChisquare();
	ndf = fge->GetNDF();
	rchi2 = chi2/ndf;
	chid = iPMT+1;
	nhit = h_on[iPMT]->GetEntries();
	//sigma = fge->GetParameter(2);
	//peakerr = fge->GetParError(1);
	//peakerr = fge->GetParError(4);
	tree->Fill();
                
	std::cout << "Handling PMT " << iPMT+1 << std::endl;
	if (peak <= 0){
	  negbad++;
	  h_on[iPMT]->Draw();
	  c1->Update();
	  c1->Print(neg_bad_name);
	  outfile << setw(10) << Form("%06d",iPMT+1) << setw(10) << Form("%0.2f", PMTinfo[iPMT][2])  << setw(10) << Form("%0.2f",peak) << setw(20) << "       negative peak" << "\n";
	}
                
	else if (peak <= 2&&peak > 0){
	  smallbad++;
	  h_on[iPMT]->Draw();
	  c1->Update();
	  c1->Print(small_bad_name);
	  outfile << setw(10) << Form("%06d",iPMT+1) << setw(10) << Form("%0.2f", PMTinfo[iPMT][2])  << setw(10) << Form("%0.2f",peak) << setw(20) << "          small peak" << "\n";
	}
                
	else if (peak > 4.5&&peak <= 6){
	  largebad++;
	  h_on[iPMT]->Draw();
	  c1->Update();
	  c1->Print(big_bad_name);
	  outfile << setw(10) << Form("%06d",iPMT+1) << setw(10) << Form("%0.2f", PMTinfo[iPMT][2])  << setw(10) << Form("%0.2f",peak) << setw(20) << "          large peak" << "\n";
	}
                
	else if (peak > 6){
	  posbad++;
	  h_on[iPMT]->Draw();
	  c1->Update();
	  c1->Print(pos_bad_name);
	  outfile << setw(10) << Form("%06d",iPMT+1) << setw(10) << Form("%0.2f", PMTinfo[iPMT][2])  << setw(10) << Form("%0.2f",peak) << setw(20) << "     larger 6pC peak" << "\n";
	}
                
	else if (peak > 2 && peak < 4.5){
	  //if (h_on[iPMT]->GetEntries()==0) continue;
	  h_on[iPMT]->Draw();
	  c1->Update();
	  c1->Print(okay_fit_name);
	  //outfile << setw(10) << Form("%06d",iPMT+1) << setw(10) << Form("%.2f", PMTinfo[iPMT][2]+hvshift)  << setw(10) << Form("%.2f",peak) << setw(20) << "     larger 6pC peak" << "\n";
	}
	    
	h_spe_peak->Fill(peak);
	h_spe_chi2->Fill(chi2/ndf);
            
	fge->Delete();
    }
            
    h_spe_peak->Write();
    h_spe_chi2->Write();
    tree->Write();
    
    c1->Print(neg_bad_name+"]");
    c1->Print(small_bad_name+"]");
    c1->Print(big_bad_name+"]");
    c1->Print(pos_bad_name+"]");
    c1->Print(okay_fit_name+"]");

    outfile << "\n" << "\n" << "================================================\n";
    outfile << "In total " << negbad << " PMTs failed the fit to negative peak." << std::endl;
    outfile << "In total " << smallbad << " PMTs failed the fit to small not okay peak." << std::endl;
    outfile << "In total " << largebad << " PMTs failed the fit to large but okay peak." << std::endl;
    outfile << "In total " << posbad << " PMTs failed the fit to too large peak." << std::endl;
    outfile << "In total " << nofit << " PMTs have no data." << std::endl;
        
        
    f->Close();
    fout->Close();
    outfile.close();
    outfile2.close();

    h_spe_peak->Clear();
    h_spe_chi2->Clear();
    tree->Clear();
    
}
    
#ifndef __CINT__
    int getArgs(int argc, char* argv[]){
        
        while( (argc > 1) && (argv[1][0] == '-') ){
            switch(argv[1][1]){
                    
                case 'r':
                    runno    = atoi(argv[2]);
                    ++argv; --argc;
                    break;
                    
                case 'c':
                    ConnectionFile = argv[2];
                    ++argv; --argc;
                    break;
                    
                case 'h':
                    hv_shift = atof(argv[2]);
                    ++argv; --argc;
                    break;
                    
                case 'p':
                    SelectPMT = atoi(argv[2]);
                    ++argv; --argc;
                    break;

	        case 't':
                    PMTtype = argv[2];
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

