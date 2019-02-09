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
#include <TFitResultPtr.h>
#include <TVirtualFitter.h>
#include <TPaveStats.h>

using namespace std;

void hv_curve_hk(){
    gStyle->SetOptTitle(kTRUE);
    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptFit(1111);
    gStyle->SetPalette(kCool);
    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    
    
    const int nfile = 7;//number of files to read in
    const int MAXPM = 11146;
    const int MAXHKPM = 136;
    const int nPMTtype = 6;
    
    Double_t PMTinfo[MAXPM][2] = {0};//flags
    std::vector<Double_t> HKPMThv(MAXHKPM,0);
    Int_t hksize;
    //double year;
    //char skip[100];
    
    int pmtid;
    
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
            //PMTinfo[cableid-1][1]=pmtz;
            PMTinfo[cableid-1][1]=oldhv;
        }
    }
    
    
    connect.close();
    
    ofstream outtxt1;
    outtxt1.open("badcables.txt");
    outtxt1 << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(10) << "    HV [V]" << setw(10) << "  Shift[V]" << setw(10) << " Peak [pC]" << "\n";
    
    ofstream outtxt2;
    outtxt2.open("badfitting.txt");
    outtxt2 << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(8) << " LHV [V]" << setw(8) << " HHV [V]" << setw(25) << "    Norm Factor" << setw(20) << "     Index" << setw(15) << "      Chi2" << setw(15) << "    Prob" << setw(25) << "                  Comment" << "\n";
    
    ofstream outtxt3;
    outtxt3.open("badokfitting.txt");
    outtxt3 << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(15) << "    Chi2"  << setw(15) << "    Prob" << setw(25) << "    Norm Factor" << setw(20) << "     Index" << "\n";
    
    //Int_t runno[] = {78565, 78561, 78559, 78568};
    //Double_t hvshift[] = {-100, -50, 0, 50};
    //Double_t threshold[] = {-0.69, -0.69, -0.69, -0.72};
    Int_t runno[] = {80282, 80278, 80275, 80269, 80265, 80263, 80254};
    Double_t hvshift[] = {-75, -50, -25, 0, 25, 50, 75};
    Double_t threshold[] = {-0.69, -0.69, -0.69, -0.69, -0.69, -0.69, -0.69};
    
    TFile *f[nfile];
    //std::vector<std::vector<Double_t>> sk2peak(nfile);
    //std::vector<std::vector<Double_t>> sk3peak(nfile);
    std::vector<std::vector<Double_t>> hkpeak(nfile,vector<Double_t>(MAXHKPM,0));
    //std::vector<std::vector<Double_t>> sk2peakerr(nfile);
    //std::vector<std::vector<Double_t>> sk3peakerr(nfile);
    std::vector<std::vector<Double_t>> hkpeakerr(nfile,vector<Double_t>(MAXHKPM,0));
    //std::vector<std::vector<Double_t>> sk2hv(nfile);
    //std::vector<std::vector<Double_t>> sk3hv(nfile);
    std::vector<std::vector<Double_t>> hkhv(nfile,vector<Double_t>(MAXHKPM,0));
    //std::vector<std::vector<Int_t>> sk2cable(nfile);
    //std::vector<std::vector<Int_t>> sk3cable(nfile);
    std::vector<std::vector<Int_t>> hkcable(nfile,vector<Int_t>(MAXHKPM,0));
    std::vector<std::vector<Double_t>> hksigma(nfile,vector<Double_t>(MAXHKPM,0));
    
    for (Int_t ifile = 0; ifile < nfile; ifile++){
        std::cout << Form("Reading fit result file %d", runno[ifile]) << std::endl;
        f[ifile] = new TFile(Form("fit_result_%d.root",runno[ifile]),"read");
        TTree *tree = (TTree*)(f[ifile]->Get("spe"));
        Int_t entry = tree->GetEntries();
        std::cout << entry << std::endl;
        //c1->Print(Form("skcable_%d.pdf[", ifile));
        //TH1D *hskcable = new TH1D("hskcable", "hskcable", 100, 0, 11200);
        
        Double_t peak;
        Double_t peakerr;
        Double_t highv;
        Double_t sigma;
        Int_t chid;
        
        tree->SetBranchAddress("HV", &highv);
        tree->SetBranchAddress("Peak", &peak);
        tree->SetBranchAddress("Peakerr",&peakerr);
        tree->SetBranchAddress("Sigma", &sigma);
        tree->SetBranchAddress("Channel", &chid);
        
        for (Int_t iPMT = 0; iPMT < entry; iPMT++){
            tree->GetEntry(iPMT);
            if (PMTinfo[chid-1][0] != 6) continue;
            //if (chid == 3761 || chid == 69) continue;//the two cables don't show up in the nominal hv run
            
            std::cout << highv << " " << peak << " " << chid << std::endl;
            hkpeak[ifile][chid-1] = peak;
            hkhv[ifile][chid-1] = highv;
            hkcable[ifile][chid-1] = chid;
            //TH1D *h2 = (TH1D*)f[ifile]->Get(Form("h_spe_onoff_%d", chid));
            //hkpeakerr[ifile].push_back(h2->GetFunction("f1peakall")->GetParError(1));
            hkpeakerr[ifile][chid-1] = peakerr;
            //hkpeakerr[ifile].push_back(sqrt(pow(peakerr,2.0)+pow(0.055,2.0)));
            //cout << "Fit peak err " << sqrt(pow(peakerr,2.0)+pow(0.055,2.0)) << endl;
            cout << "Fit peak err " << Form("%1.2e", peakerr) << endl;
            //hkpeakerr[ifile].push_back(5e-3*sigma);
            hksigma[ifile][chid-1] = sigma;
            HKPMThv[chid-1] = PMTinfo[iPMT][1];
        }
        
        f[ifile]->Close();
        std::cout << Form("Closed fit result file %d", runno[ifile]) << std::endl;
    }
    
    
    
    
    TFile *fout = new TFile("hvscan_parameter.root", "recreate");
    //TTree *trsk2 = new TTree("hvscan_sk2", "SK2 PMT HV Scan Parameter");
    //TTree *trsk3 = new TTree("hvscan_sk3", "SK3 PMT HV Scan Parameter");
    TTree *trhk = new TTree("hvscan_hk", "HK PMT HV Scan Parameter");
    
    
    Int_t channelid_hk;
    Double_t geighthv_hk;
    Double_t gfourhv_hk;
    Double_t gshifthv_hk;
    Double_t geighthverr_hk;
    Double_t gfourhverr_hk;
    Double_t rchi2_hk;
    Double_t prob_hk;
    Double_t norm_hk;
    Double_t beta_hk;
    Double_t offset_hk = 0;
    Double_t normerr_hk;
    Double_t betaerr_hk;
    Double_t offseterr_hk = 0;
    Double_t resolution_hk;
    
    trhk->Branch("Channel", &channelid_hk, "channelid_hk/I");
    trhk->Branch("HighGainHV", &geighthv_hk, "geighthv_hk/D");
    trhk->Branch("HighGainHVerr", &geighthverr_hk, "geighthverr_hk/D");
    trhk->Branch("LowGainHV", &gfourhv_hk, "gfourhv_hk/D");
    trhk->Branch("LowGainHVerr", &gfourhverr_hk, "gfourhverr_hk/D");
    trhk->Branch("HVshiftHighGain", &gshifthv_hk, "gshifthv_hk/D");
    trhk->Branch("Chi2", &rchi2_hk, "rchi2_hk/D");
    trhk->Branch("Prob", &prob_hk, "prob_hk/D");
    trhk->Branch("Normalization", &norm_hk, "norm_hk/D");
    trhk->Branch("Index", &beta_hk, "beta_hk/D");
    trhk->Branch("Resolution", &resolution_hk, "resolution_hk/D");
    //trhk->Branch("Offset", &offset_hk, "offset_hk/D");
    trhk->Branch("NormErr", &normerr_hk, "normerr_hk/D");
    trhk->Branch("IndexErr", &betaerr_hk, "betaerr_hk/D");
    //trhk->Branch("OffsetErr", &offseterr_hk, "offseterr_hk/D");
    
    hksize = hkcable[0].size();
    TGraphErrors *ghv_hk[hksize];
    //TGraphErrors *gthr_sk[MAXPM];
    //TGraphErrors *gthr_hk[MAXHKPM];
   
    TF1 *fHVhk[hksize];
    TF1 *fHVinvhk[hksize];
    TF1 *fHVinvhkerr[hksize];
    
    
    for (Int_t q = 0; q < hksize; q++){
        ghv_hk[q] = new TGraphErrors();
        ghv_hk[q]->SetMarkerStyle(8);
        ghv_hk[q]->SetMarkerSize(1);
        ghv_hk[q]->SetMarkerColor(4);
        ghv_hk[q]->SetLineColor(15);
        ghv_hk[q]->SetLineStyle(3);
        ghv_hk[q]->GetXaxis()->SetTitle("HV [V]");
        //ghv_hk[q]->GetXaxis()->SetRangeUser(1500,2500);
        ghv_hk[q]->GetYaxis()->SetTitle("Peak [pC]");
        ghv_hk[q]->SetName(Form("HK_PMT_HVscan_Cable_%06d", hkcable[0][q]));
        ghv_hk[q]->SetTitle(Form("HK_PMT_HVscan_Cable_%06d", hkcable[0][q]));
        
        fHVhk[q] = new TF1(Form("fHVhk%d", q), "[0]*pow(x,[1])", 1500, 2500);
        fHVinvhk[q] = new TF1(Form("fHVinvhk%d", q), "pow(x/[0],1./[1])", 1500, 2500);
        fHVinvhkerr[q] = new TF1(Form("fHVinvhkerr%d", q), "pow(x/[0],1./[1])", 1500, 2500);
        fHVhk[q]->SetLineColor(4);
        fHVhk[q]->SetLineStyle(3);
        fHVhk[q]->SetLineWidth(1);
        fHVhk[q]->SetParName(0, "Factor");
        fHVhk[q]->SetParName(1, "Index #beta");
        fHVhk[q]->SetParameter(0, 4e-20);
        fHVhk[q]->SetParameter(1, 6);
    }
    
    for (Int_t ifile = 0; ifile < nfile; ifile++){
        for (Int_t q = 0; q < hksize; q++){
            if (hkpeak[ifile][q] <= 0 || hkpeak[ifile][q] > 10){
                outtxt1 << setw(10) << hkcable[ifile][q] << setw(4) << "  HK" << setw(10) << hkhv[ifile][q] << setw(10) << hvshift[ifile] << std::endl;
                continue;
            }
            ghv_hk[q]->SetPoint(ifile, hkhv[ifile][q], hkpeak[ifile][q]);
            ghv_hk[q]->SetPointError(ifile, 0.2, hkpeakerr[ifile][q]);
            if (hksigma[ifile][q] < 0.1) ghv_hk[q]->RemovePoint(ifile);
        }
    }
    
    Double_t targetQ = 2.884; // 1.8e7 gain
    // 1.8e7/(1e-12/1.60217657e-19) = 2.884
    Double_t targetQHPK = 2.243; // 1.4e7 gain
    //1.4e7/(1e-12/1.60217657e-19) = 2.243
    
    for (Int_t q = 0; q < hksize; q++){
        //ghv_hk[q]->GetXaxis()->SetLimits(1500,2500);
        ghv_hk[q]->Fit(fHVhk[q], "BQN0");
        TVirtualFitter *gfitter = TVirtualFitter::Fitter(ghv_hk[q]);
        gfitter->SetPrecision(0.05);
        TFitResultPtr fitr = ghv_hk[q]->Fit(fHVhk[q], "SB+");
        int status = int(fitr);
        std::cout << Form("Fitting HK cable %06d", q) << std::endl;
        //ghv_hk[q]->Draw();
        //c1->Update();
        ghv_hk[q]->Write();
        channelid_hk = hkcable[0][q];
        norm_hk = fHVhk[q]->GetParameter(0);
        beta_hk = fHVhk[q]->GetParameter(1);
        //offset_hk = fHVhk[q]->GetParameter(2);
        normerr_hk = fHVhk[q]->GetParError(0);
        betaerr_hk = fHVhk[q]->GetParError(1);
        //offseterr_hk = fHVhk[q]->GetParError(2);
        rchi2_hk = fHVhk[q]->GetChisquare()/fHVhk[q]->GetNDF();
        prob_hk = fHVhk[q]->GetProb();
        fHVinvhk[q]->SetParameter(0, fHVhk[q]->GetParameter(0));
        fHVinvhk[q]->SetParameter(1, fHVhk[q]->GetParameter(1));
        fHVinvhkerr[q]->SetParameter(0, fHVhk[q]->GetParameter(0)-fHVhk[q]->GetParError(0));
        fHVinvhkerr[q]->SetParameter(1, fHVhk[q]->GetParameter(1)-fHVhk[q]->GetParError(1));
        //fHVinvhk[q]->SetParameter(2, fHVhk[q]->GetParameter(2));
        geighthv_hk = fHVinvhk[q]->Eval(targetQ);
        geighthverr_hk = fHVinvhkerr[q]->Eval(targetQ) - geighthv_hk;
        gfourhv_hk = fHVinvhk[q]->Eval(targetQHPK);
        gfourhverr_hk = fHVinvhkerr[q]->Eval(targetQHPK) - gfourhv_hk;
        gshifthv_hk = geighthv_hk - HKPMThv[q];
        resolution_hk = hksigma[2][q]/hkpeak[2][q];
        trhk->Fill();
        
        if (status == 4){
            fHVhk[q]->SetLineColor(2);
            outtxt2 << setw(10) << hkcable[2][q] << setw(4) << "  HK" << setw(8) << Form("%1.2f", gfourhv_hk) << setw(8) << Form("%1.2f", geighthv_hk) << setw(25) << Form("%1.2e(%1.2e)", norm_hk, normerr_hk) << setw(20) << Form("%1.2f(%1.2f)", beta_hk, betaerr_hk) << setw(15) << Form("%1.3f", rchi2_hk) << setw(15) << Form("%1.2e", fHVhk[q]->GetProb()) << setw(25) << "   Failed_to_find_minimum" << "\n";
        }
        /*else if (rchi2_hk > 20 && prob_hk < 0.0001){
            fHVhk[q]->SetLineColor(6);
            outtxt2 << setw(10) << hkcable[2][q] << setw(4) << "  HK" << setw(8) << Form("%1.2f", gfourhv_hk) << setw(8) << Form("%1.2f", geighthv_hk) << setw(25) << Form("%1.2e(%1.2e)", norm_hk, normerr_hk) << setw(20) << Form("%1.2f(%1.2f)", beta_hk, betaerr_hk) << setw(15) << Form("%1.3f", rchi2_hk) << setw(15) << Form("%1.2e", fHVhk[q]->GetProb()) << setw(25) << "                Chi2_>>_1" << "\n";
        }*/
        
        //mghv_hk->Add(ghv_hk[q]);
        //mghv_hk->Add(fHVhk[q]);
        //mgthr_hk->Add(gthr_hk[id]);
        //fHV->Clear();
    }
    
    c1->Print("HK_SPE_HV.pdf[");
    fHVhk[0]->Draw();
    for (Int_t q = 1; q < hksize; q++){
        //ghv_hk[q]->Draw("same");
        //if (fHVhk[q]->GetLineColor() != 4) continue;
        fHVhk[q]->Draw("same");
        c1->Update();
    }
    //mghv_hk->Fit(fHV, "BQ");
    //mghv_hk->GetXaxis()->SetLimits(1500,2500);
    //legend2->Draw("same");
    //c1->Update();
    c1->Modified();
    c1->Print("HK_SPE_HV.pdf");
    c1->Print("HK_SPE_HV.pdf]");
    
    
    TPaveText *t = new TPaveText(.4,.15,.9,.25,"NDC");
    TPaveText *trun = new TPaveText(.4,.25,.9,.35,"NDC");
    t->SetFillStyle(0);
    t->SetBorderSize(0);
    t->SetTextColor(kRed);
    trun->SetFillStyle(0);
    trun->SetBorderSize(0);
    trun->SetTextColor(kRed);
    //trun->AddText(Form("Run %s", runlist.Data()));
    
    c1->Print("HK_HV_Curves.pdf[");
    for (Int_t q = 0; q < hksize; q++){
        ghv_hk[q]->Draw("AP");
        t->Clear();
        t->AddText(Form("HV at 1.4e7 gain : %4.2f#pm%4.2f [V]", fHVinvhk[q]->Eval(targetQHPK), fHVinvhkerr[q]->Eval(targetQHPK) - fHVinvhk[q]->Eval(targetQHPK)));
        ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlue);
        t->AddText(Form("HV at 1.8e7 gain   : %4.2f#pm%4.2f [V]", fHVinvhk[q]->Eval(targetQ),fHVinvhkerr[q]->Eval(targetQ) - fHVinvhk[q]->Eval(targetQ)));
        ((TText*)t->GetListOfLines()->Last())->SetTextColor(kMagenta);
        t->Draw("same");
        trun->Draw("same");
        c1->Update();
        
        TPaveStats *ps = (TPaveStats*)ghv_hk[q]->GetListOfFunctions()->FindObject("stats");
        ps->SetX1NDC(0.1);
        ps->SetX2NDC(0.4);
        c1->Modified();
        c1->Update();
        c1->Print("HK_HV_Curves.pdf");
    }
    c1->Print("HK_HV_Curves.pdf]");
    
    trhk->Write();
    //mghv_sk->Write();
    //mghv_hk->Write();
    fout->Close();
    outtxt1.close();
    outtxt2.close();
    outtxt3.close();
}
