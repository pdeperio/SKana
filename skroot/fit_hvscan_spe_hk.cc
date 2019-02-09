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

using namespace std;

TF1* FitSK(TH1D* h, Int_t rebin = 1, Float_t hv = 0, Float_t PMTflag = 6){
    Int_t npx = 1000;
    Int_t entry = h->GetEntries();
    Int_t nbin = h->GetNbinsX();
    if (rebin > 1) h->Rebin(rebin);
    
    
    TF1 *func1pe = new TF1("fgaus","gaus",-1+hv*0.008, 12+hv*0.008);
    TF1 *func1pe_new = new TF1("fgaus_new", "gaus",-1+hv*0.008, 12+hv*0.008);
    TF1 *func2peak_bs = new TF1("f2peakbs","gaus(0)+gaus(3)+(0.5*(([6]*[0]*((TMath::Erf(((x-[4])/[5])))+(1-TMath::Erf((x-[7]*[1])/[2]))-1))))",-1+hv*0.008,12+hv*0.008);
    TF1 *func0peak = new TF1("fgaus0","gaus",-1+hv*0.008, 12+hv*0.008);
    TF1 *func1peak = new TF1("fgaus1","gaus",-1+hv*0.008, 12+hv*0.008);
    TF1 *funcbs = new TF1("fbs","(0.5*(([6]*[0]*((TMath::Erf(((x-[4])/[5])))+(1-TMath::Erf((x-[7]*[1])/[2]))-1))))",-1+hv*0.008, 12+hv*0.008);
    
    
    func1pe->SetLineColor(2);
    func1pe_new->SetLineColor(46);
    func2peak_bs->SetLineColor(38);
    func0peak->SetLineColor(30);
    func1peak->SetLineColor(9);
    funcbs->SetLineColor(7);
    
    func1pe->SetLineStyle(3);
    func0peak->SetLineStyle(3);
    func1peak->SetLineStyle(3);
    funcbs->SetLineStyle(3);
    
    func1pe->SetNpx(npx);
    func1pe_new->SetNpx(npx);
    func2peak_bs->SetNpx(npx);
    func1peak->SetNpx(npx);
    func0peak->SetNpx(npx);
    funcbs->SetNpx(npx);
    
    Double_t peakx_pre[2] = {2.9+hv*0.008, 1.5+hv*0.008};
    Double_t peaky_pre[2] = {h->GetBinContent(h->FindBin(peakx_pre[0])), h->GetBinContent(h->FindBin(peakx_pre[1]))};
    
    //preliminarily fit the main peak
    func1pe->SetParameters(h->GetMaximum(), peakx_pre[0], 1);
    func1pe->SetParLimits(0, 0, h->GetMaximum());
    func1pe->SetParLimits(1, 1.5+hv*0.008, 10+hv*0.008);
    func1pe->SetParLimits(2, 0, 3);
    h->Fit(func1pe, "BNQ0", "", 1.5+hv*0.008, 10+hv*0.008);
    
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
            func1peak->SetParameter(iSet, func1pe_new->GetParameter(iSet));
        }
        
        func2peak_bs->SetParameter(3, func1pe_new->GetParameter(0)*0.25);
        func2peak_bs->SetParameter(4, func1pe_new->GetParameter(1)-1.5*func1pe_new->GetParameter(2));
        func2peak_bs->SetParameter(5, 0.7*func1pe_new->GetParameter(2));
        func2peak_bs->SetParLimits(3, 0, 0.5*func1pe_new->GetParameter(0));
        func2peak_bs->SetParLimits(4, 0, func1pe_new->GetParameter(1)-func1pe_new->GetParameter(2));
        func2peak_bs->SetParLimits(5, 0.3*func1pe_new->GetParameter(2), func1pe_new->GetParameter(2));
        func2peak_bs->SetParName(6, "1pe BG ratio");
        func2peak_bs->SetParameter(6, 0.5);
        func2peak_bs->SetParLimits(6, 0, 1);
        func2peak_bs->SetParameter(7, 2);
        func2peak_bs->SetParLimits(7, 1, 3);
        
         h->Fit(func2peak_bs,"Q+","",-1+hv*0.008, 12+hv*0.008);
        
        for(int iSet = 0; iSet < 3; iSet++){
            func0peak->SetParameter(iSet, func2peak_bs->GetParameter(iSet+3));
            funcbs->SetParameter(iSet, func2peak_bs->GetParameter(iSet));
            funcbs->SetParameter(iSet+3, func2peak_bs->GetParameter(iSet+3));
        }
        funcbs->SetParameter(6, func2peak_bs->GetParameter(6));
        
        h->Fit(func0peak,"Q+","", -1+hv*0.008, 10+hv*0.008);
        h->Fit(func1peak,"Q+","", -1+hv*0.008, 10+hv*0.008);
        h->Fit(funcbs,"Q+","", -1+hv*0.008, 10+hv*0.008);
        
    }
    
    TF1 *func;
    if (hv >= -100) func = func2peak_bs;
    else func = func1pe_new;
    
    gPad->cd();
    
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
        if (func->GetNpar()>6) myt[1] = new TLatex(0,0,Form("SPE Peak  = %1.2e [pC]",func->GetParameter(1)>func->GetParameter(4)?func->GetParameter(1):func->GetParameter(4)));
        else myt[1] = new TLatex(0,0,Form("SPE Peak  = %1.2e [pC]",func->GetParameter(1)));
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
    func->SetParError(1, func2peak_bs->GetParError(1));
    return func;
}


void fit_hvscan_spe_hk(){
    
    gStyle->SetOptTitle(kTRUE);
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(kCool);
    
    //Int_t runno[] = {78559};
    //Double_t hvshift[] = {0};
    //Double_t threshold = {-0.69};
    Int_t runno[] = {80282, 80278, 80275, 80269, 80265, 80263, 80254};
    Double_t hvshift[] = {-75, -50, -25, 0, 25, 50, 75};
    Double_t threshold[] = {-0.69, -0.69, -0.69, -0.69, -0.69, -0.69, -0.69};
    TCanvas * c1 = new TCanvas("c1","c1",800,800);

    const int nfile = 7;//number of files to read in
    const int MAXPM = 11146;
    const int nPMTtype = 6;
    
    Float_t PMTinfo[MAXPM][3];//flag, place, hv
    //int npmt;
    //int pmtid;
    
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
            PMTinfo[cableid-1][0]=pmtflag;
            PMTinfo[cableid-1][1]=pmtz;
            PMTinfo[cableid-1][2]=oldhv;
        }
    }
    
    
    connect.close();
    
    TFile *f[nfile];
    TFile *fout[nfile];
    TH1D * h_on[MAXPM];
    TH1D * h_off[MAXPM];
    
    
   
    for (Int_t ifile = 0; ifile < nfile; ifile++){
        int negbad = 0;
        int smallbad = 0;
        int largebad =0;
        int posbad = 0;
        int nofit = 0;//counters
        
        Int_t chid;
        Double_t peak;
        //Double_t perr;
        Double_t chi2;
        Double_t peakerr;
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
        outfile.open(Form("badpeak_%d.txt", runno[ifile]));
        outfile << setw(10) << "     Cable" << setw(10) << "    HV [V]" << setw(10) << " Peak [pC]" << setw(20) << "         Err Message\n";
        
        ofstream outfile2;
        outfile2.open(Form("nofit_%d.txt", runno[ifile]));
        
        
        c1->Print(Form("neg_bad_fit_%d.pdf[", runno[ifile]));
        c1->Print(Form("small_bad_fit_%d.pdf[", runno[ifile]));
        c1->Print(Form("big_bad_fit_%d.pdf[", runno[ifile]));
        c1->Print(Form("pos_bad_fit_%d.pdf[", runno[ifile]));
        c1->Print(Form("okay_fit_sample_%d.pdf[", runno[ifile]));
        
        
        TH1D * h_spe_peak_hk = new TH1D("h_spe_peak_hk","spe_peak_hk",120,0,6);
        TH1D * h_spe_chi2_hk = new TH1D("h_spe_chi2_hk","spe_chi2_hk",200,0,40);
        
        
        h_spe_peak_hk->GetXaxis()->SetTitle("Peak [pC]");
        h_spe_chi2_hk->GetXaxis()->SetTitle("Fit Chi2/ndf");
        
        f[ifile] = new TFile(Form("spe_individual_%d_C.root",runno[ifile]),"read");
        fout[ifile] = new TFile(Form("fit_result_%d.root",runno[ifile]),"recreate");
        fout[ifile]->cd();
        
        for (Int_t iPMT = 0; iPMT < MAXPM; iPMT++){
            if (PMTinfo[iPMT][0] != 6) continue;
            h_on[iPMT] = (TH1D*)f[ifile]->Get(Form("h_spe_on_%d",iPMT+1));
            h_off[iPMT] = (TH1D*)f[ifile]->Get(Form("h_spe_off_%d",iPMT+1));
            h_on[iPMT]->Sumw2();
            h_off[iPMT]->Sumw2();
            h_on[iPMT]->Add(h_off[iPMT],-1);
            h_on[iPMT]->SetTitle(Form("PMT %d",iPMT+1));
            h_on[iPMT]->SetName(Form("h_spe_onoff_%d", iPMT+1));
            h_on[iPMT]->GetXaxis()->SetTitle("Charge [pC]");
            h_on[iPMT]->GetXaxis()->SetRangeUser(-3,15);
            h_on[iPMT]->Draw();
            //TH1D *h = (TH1D*)h_on[iPMT]->Clone();
            if (h_on[iPMT]->GetEntries() < 10 && PMTinfo[iPMT][0] == 6){
                nofit++;
                outfile2 << iPMT+1 << std::endl;
                continue;
            }
         
            if (PMTinfo[iPMT][0] == 6){
                highv = PMTinfo[iPMT][2] + hvshift[ifile];
                
                TF1 *fge = FitSK((TH1D*)h_on[iPMT], 2, hvshift[ifile], PMTinfo[iPMT][0]);
               
                c1->Update();
                h_on[iPMT]->Write();
                if (fge->GetNpar() > 6) {
                    if ((TMath::Abs(fge->GetMaximumX(0,6) - fge->GetParameter(1)) >= 1) || fge->GetParameter(1) > 10){
                        peak = fge->GetMaximumX(0, 6);
                        Double_t FWHMlow  = peak-fge->GetX(fge->Eval(peak)*0.5, 0, peak);
                        Double_t FWHMhigh = fge->GetX(fge->Eval(peak)*0.5, peak, 12) - peak;
                        sigma = (FWHMlow+FWHMhigh)/(2.*TMath::Sqrt(TMath::Log(2)*2));
                        peakerr = sigma * 1e-4;
                    }
                    else {
                        peak = fge->GetParameter(1)>fge->GetParameter(4)?fge->GetParameter(1):fge->GetParameter(4);
                        sigma = fge->GetParameter(1)>fge->GetParameter(4)?fge->GetParameter(2):fge->GetParameter(5);
                        peakerr = fge->GetParameter(1)>fge->GetParameter(4)?fge->GetParError(1):fge->GetParError(4);
                    }
                }
                else {
                    if (fge->GetParameter(1) < 0 || fge->GetParameter(1) > 10){
                        peak = fge->GetMaximumX(0,6);
                        Double_t FWHMlow  = peak- fge->GetX(fge->Eval(peak)*0.5, 0, peak);
                        Double_t FWHMhigh = fge->GetX(fge->Eval(peak)*0.5, peak, 12) - peak;
                        sigma = (FWHMlow+FWHMhigh)/(2.*TMath::Sqrt(TMath::Log(2)*2));
                        peakerr = sigma * 1e-4;
                    }
                    else {
                        peak = fge->GetParameter(1);
                        sigma = fge->GetParameter(2);
                        peakerr = fge->GetParError(1);
                    }
                }
                chi2 = fge->GetChisquare();
                ndf = fge->GetNDF();
                rchi2 = chi2/ndf;
                chid = iPMT+1;
                nhit = h_on[iPMT]->GetEntries();
                //peakerr = fge->GetParError(1);
                //peakerr = fge->GetParError(4);
                tree->Fill();
            }
            
            std::cout << "Handling PMT " << iPMT+1 << std::endl;
            if (peak <= 0){
                negbad++;
                //h_on[iPMT]->Draw();
                c1->Update();
                c1->Print(Form("neg_bad_fit_%d.pdf", runno[ifile]));
                outfile << setw(10) << Form("%06d",iPMT+1) << setw(10) << Form("%0.2f", PMTinfo[iPMT][2]+hvshift[ifile])  << setw(10) << Form("%0.2f",peak) << setw(20) << "       negative peak" << "\n";
            }
            
            else if (peak <= 2&&peak > 0){
                smallbad++;
                //h_on[iPMT]->Draw();
                c1->Update();
                c1->Print(Form("small_bad_fit_%d.pdf", runno[ifile]));
                outfile << setw(10) << Form("%06d",iPMT+1) << setw(10) << Form("%0.2f", PMTinfo[iPMT][2]+hvshift[ifile])  << setw(10) << Form("%0.2f",peak) << setw(20) << "          small peak" << "\n";
            }
            
            else if (peak > 4.5&&peak <= 6){
                largebad++;
                //h_on[iPMT]->Draw();
                c1->Update();
                c1->Print(Form("big_bad_fit_%d.pdf", runno[ifile]));
                outfile << setw(10) << Form("%06d",iPMT+1) << setw(10) << Form("%0.2f", PMTinfo[iPMT][2]+hvshift[ifile])  << setw(10) << Form("%0.2f",peak) << setw(20) << "          large peak" << "\n";
            }
            
            else if (peak > 6){
                posbad++;
                //h_on[iPMT]->Draw();
                c1->Update();
                c1->Print(Form("pos_bad_fit_%d.pdf", runno[ifile]));
                outfile << setw(10) << Form("%06d",iPMT+1) << setw(10) << Form("%0.2f", PMTinfo[iPMT][2]+hvshift[ifile])  << setw(10) << Form("%0.2f",peak) << setw(20) << "     larger 6pC peak" << "\n";
            }
            
            else if (peak > 2&&peak < 4.5){
                //if (h_on[iPMT]->GetEntries()==0) continue;
                //h_on[iPMT]->Draw();
                c1->Update();
                c1->Print(Form("okay_fit_sample_%d.pdf", runno[ifile]));
                //outfile << setw(10) << Form("%06d",iPMT+1) << setw(10) << Form("%.2f", PMTinfo[iPMT][2]+hvshift[ifile])  << setw(10) << Form("%.2f",peak) << setw(20) << "     larger 6pC peak" << "\n";
            }
    
            else if (PMTinfo[iPMT][0]==6 && PMTinfo[iPMT][1]<=1656){
                //hkpeak[ifile].push_back(peak);
                //hkchi2[ifile].push_back(chi2/ndf);
                //skperr[ifile].push_back(peakerr);
                //skpvalue[ifile].push_back(pvalue);
                h_spe_peak_hk->Fill(peak);
                h_spe_chi2_hk->Fill(chi2/ndf);
                //h_spe_p_hk->Fill(pvalue);
                //h_spe_perr_hk->Fill(peakerr);
            }
            
        }
        h_spe_peak_hk->Write();
        h_spe_chi2_hk->Write();
        //h_spe_p_hk->Write();
        //h_spe_perr_hk->Write();
        tree->Write();
        
        auto *legend = new TLegend(0.1,0.85,0.7,0.9);
        //legend->SetNColumns(2);
        //legend->AddEntry(h_spe_peak, "SK PMT", "l");
        legend->AddEntry(h_spe_peak_hk, "HK PMT", "l");
        
        
        c1->Print(Form("hk_spe_peak_%06d.pdf[",runno[ifile]));
        //c1->SetLogy();
        //gStyle->SetOptStats(0);
        //h_spe_peak->SetLineStyle(2);
        //h_spe_peak->SetLineColor(kRed);
        //h_spe_peak->Draw();
        h_spe_peak_hk->SetLineColor(kBlue);
        h_spe_peak_hk->SetFillColor(kBlue);
        h_spe_peak_hk->Draw();
        legend->Draw("same");
        c1->Update();
        c1->Print(Form("hk_spe_peak_%06d.pdf",runno[ifile]));
        c1->Print(Form("hk_spe_peak_%06d.pdf]",runno[ifile]));
        
        c1->Print(Form("hk_spe_chi2_%06d.pdf[",runno[ifile]));
        //h_spe_chi2->SetLineStyle(2);
        //h_spe_chi2->SetLineColor(kRed);
        //h_spe_chi2->Draw();
        h_spe_chi2_hk->SetLineColor(kBlue);
        h_spe_chi2_hk->SetFillColor(kBlue);
        h_spe_chi2_hk->Draw();
        legend->Draw("same");
        c1->Update();
        c1->Print(Form("hk_spe_chi2_%06d.pdf",runno[ifile]));
        c1->Print(Form("hk_spe_chi2_%06d.pdf]",runno[ifile]));
        
    
        c1->Print(Form("neg_bad_fit_%d.pdf]", runno[ifile]));
        c1->Print(Form("small_bad_fit_%d.pdf]", runno[ifile]));
        c1->Print(Form("big_bad_fit_%d.pdf]", runno[ifile]));
        c1->Print(Form("pos_bad_fit_%d.pdf]", runno[ifile]));
        c1->Print(Form("okay_fit_sample_%d.pdf]", runno[ifile]));
        
        outfile << "\n" << "\n" << "================================================\n";
        outfile << "In total " << negbad << " PMTs failed the fit to negative peak." << std::endl;
        outfile << "In total " << smallbad << " PMTs failed the fit to small not okay peak." << std::endl;
        outfile << "In total " << largebad << " PMTs failed the fit to large but okay peak." << std::endl;
        outfile << "In total " << posbad << " PMTs failed the fit to too large peak." << std::endl;
        outfile << "In total " << nofit << " PMTs have no data." << std::endl;
        
        
        f[ifile]->Close();
        fout[ifile]->Close();
        outfile.close();
        outfile2.close();
        
        
        h_spe_peak_hk->Clear();
        h_spe_chi2_hk->Clear();
        //h_spe_p_hk->Clear();
        //h_spe_perr_hk->Clear();
        tree->Clear();
        c1->Clear();
        
        
    }
   
    
    
}
