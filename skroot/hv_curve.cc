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

void hv_curve(){
    gStyle->SetOptTitle(kTRUE);
    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptFit(1111);
    gStyle->SetPalette(kCool);
    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    
    
    const int nfile = 7;//number of files to read in
    const int MAXPM = 11146;
    const int MAXHKPM = 136;
    const int nPMTtype = 6;
    
    Int_t PMTinfo[MAXPM][2] = {0};//flags
    std::vector<Double_t> SKPMThv(MAXPM,0);
    int npmt;
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
    outtxt1 << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(10) << "    HV [V]" << setw(10) << "  Shift[V]" << setw(10) << " Peak [pC]" << setw(10) << "      Qisk" << "\n";
    
    ofstream outtxt2;
    outtxt2.open("badfitting.txt");
    outtxt2 << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(8) << " LHV [V]" << setw(8) << " HHV [V]" << setw(25) << "    Norm Factor" << setw(20) << "     Index" << setw(15) << "      Chi2" << setw(15) << "    Prob" << setw(25) << "                  Comment" << "\n";
    
    ofstream outtxt3;
    outtxt3.open("badokfitting.txt");
    outtxt3 << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(15) << "    Chi2"  << setw(15) << "    Prob" << setw(25) << "    Norm Factor" << setw(20) << "     Index" << "\n";
    
    //Int_t runno[] = {78565, 78561, 78559, 78568};
    //Double_t hvshift[] = {-100, -50, 0, 50};
    //Double_t threshold[] = {-0.69, -0.69, -0.69, -0.72};
    Int_t runno[] = {80152, 80157, 80159, 80148, 80161, 80163, 80167};
    Double_t hvshift[] = {-75, -50, -25, 0, 25, 50, 75};
    Double_t threshold[] = {-0.69, -0.69, -0.69, -0.69, -0.69, -0.69, -0.69};
    
    TFile *f[nfile];
    std::vector<std::vector<Double_t>> skpeak(nfile, std::vector<Double_t>(MAXPM,0));
    //std::vector<std::vector<Double_t>> sk3peak(nfile, std::vector<Double_t>(MAXPM,0));
    //std::vector<std::vector<Double_t>> hkpeak(nfile);
    std::vector<std::vector<Double_t>> skpeakerr(nfile, std::vector<Double_t>(MAXPM,0));
    //std::vector<std::vector<Double_t>> sk3peakerr(nfile, std::vector<Double_t>(MAXPM,0));
    //std::vector<std::vector<Double_t>> hkpeakerr(nfile);
    std::vector<std::vector<Double_t>> skhv(nfile, std::vector<Double_t>(MAXPM,0));
    //std::vector<std::vector<Double_t>> sk3hv(nfile, std::vector<Double_t>(MAXPM,0));
    //std::vector<std::vector<Double_t>> hkhv(nfile);
    std::vector<std::vector<Int_t>> skcable(nfile, std::vector<Int_t>(MAXPM,0));
    //std::vector<std::vector<Int_t>> sk3cable(nfile, std::vector<Double_t>(MAXPM,0));
    //std::vector<std::vector<Int_t>> hkcable(nfile);
    std::vector<std::vector<Int_t>> sknhit(nfile, std::vector<Int_t>(MAXPM,0));
    //std::vector<std::vector<Int_t>> sk3nhit(nfile, std::vector<Double_t>(MAXPM,0));
    std::vector<std::vector<Double_t>> sksigma(nfile, std::vector<Double_t>(MAXPM,0));
    
    for (Int_t ifile = 0; ifile < nfile; ifile++){
        std::cout << Form("Reading fit result file %d", runno[ifile]) << std::endl;
        f[ifile] = new TFile(Form("fit_result_%d_sk.root",runno[ifile]),"read");
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
        Int_t nhits;
        
        tree->SetBranchAddress("HV", &highv);
        tree->SetBranchAddress("Peak", &peak);
        tree->SetBranchAddress("Peakerr", &peakerr);
        tree->SetBranchAddress("Sigma", &sigma);
        tree->SetBranchAddress("Channel", &chid);
        tree->SetBranchAddress("Nhits", &nhits);
        
        for (Int_t iPMT = 0; iPMT < entry; iPMT++){
            tree->GetEntry(iPMT);
            //if (PMTinfo[chid-1][0]!=3 || PMTinfo[chid-1][0]!=4) continue;
            //if (chid == 3761 || chid == 69) continue;//the two cables don't show up in the nominal hv run
            std::cout << highv << " " << peak << " " << chid << std::endl;
            if (PMTinfo[chid-1][0]==3){
                
                if (nhits < 10){
                    outtxt1 << setw(10) << skcable[ifile][chid-1] << setw(4) << " SK2" << setw(10) << skhv[ifile][chid-1] << setw(10) << hvshift[ifile] << setw(10) << nhits << std::endl;
                    continue;
                }
                skpeak[ifile][chid-1] = peak;
                skhv[ifile][chid-1] = highv;
                skcable[ifile][chid-1] = chid;
                sksigma[ifile][chid-1] = sigma;
                //TH1D *h1 = (TH1D*)f[ifile]->Get(Form("h_spe_onoff_%d", chid));
                skpeakerr[ifile][chid-1] = pow((pow(peakerr,2)+pow(0.22,2)), 1/2);
                SKPMThv[chid-1] = PMTinfo[chid-1][1];
                sknhit[ifile][chid-1] = nhits;
                
            }
            
            else if (PMTinfo[chid-1][0]==4){
                if (nhits < 10) {
                    outtxt1 << setw(10) << skcable[ifile][chid-1] << setw(4) << " SK3" << setw(10) << skhv[ifile][chid-1] << setw(10) << hvshift[ifile] << setw(10) << nhits << std::endl;
                    continue;
                }
                skpeak[ifile][chid-1] = peak;
                skhv[ifile][chid-1] = highv;
                skcable[ifile][chid-1] = chid;
                sksigma[ifile][chid-1] = sigma;
                //TH1D *h1 = (TH1D*)f[ifile]->Get(Form("h_spe_onoff_%d", chid));
                //sk3peakerr[ifile].push_back(h1->GetFunction("f1peakall")->GetParError(1));
                skpeakerr[ifile][chid-1] = pow((pow(peakerr,2)+pow(0.22,2)), 1/2);
                SKPMThv[chid-1] = PMTinfo[chid-1][1];
                sknhit[ifile][chid-1] = nhits;
            }
            
        }
        
        f[ifile]->Close();
        std::cout << Form("Closed fit result file %d", runno[ifile]) << std::endl;
    }
    
    
    
    
    TFile *fout = new TFile("hvscan_parameter_sk.root", "recreate");
    TTree *trsk2 = new TTree("hvscan_sk2", "SK2 PMT HV Scan Parameter");
    TTree *trsk3 = new TTree("hvscan_sk3", "SK3 PMT HV Scan Parameter");
    //TTree *trhk = new TTree("hvscan_hk", "HK PMT HV Scan Parameter");
    
    Int_t channelid_sk2;
    Double_t geighthv_sk2;
    Double_t geighthverr_sk2;
    Double_t gfourhv_sk2;
    Double_t gfourhverr_sk2;
    Double_t gshifthv_sk2;
    Double_t rchi2_sk2;
    Double_t prob_sk2;
    Double_t norm_sk2;
    Double_t beta_sk2;
    Double_t offset_sk2 = 0;
    Double_t normerr_sk2;
    Double_t betaerr_sk2;
    Double_t offseterr_sk2 = 0;
    Int_t channelid_sk3;
    Double_t geighthv_sk3;
    Double_t geighthverr_sk3;
    Double_t gfourhv_sk3;
    Double_t gfourhverr_sk3;
    Double_t gshifthv_sk3;
    Double_t rchi2_sk3;
    Double_t prob_sk3;
    Double_t norm_sk3;
    Double_t beta_sk3;
    Double_t offset_sk3 = 0;
    Double_t normerr_sk3;
    Double_t betaerr_sk3;
    Double_t offseterr_sk3 = 0;
    //TFile *fscan = new TFile("ID_Scan.root","recreate");
    
    trsk2->Branch("Channel", &channelid_sk2, "channelid_sk2/I");
    trsk2->Branch("HighGainHV", &geighthv_sk2, "geighthv_sk2/D");
    trsk2->Branch("HighGainHVerr", &geighthverr_sk2, "geighthverr_sk2/D");
    trsk2->Branch("LowGainHV", &gfourhv_sk2, "gfourhv_sk2/D");
    trsk2->Branch("LowGainHVerr", &gfourhverr_sk2, "gfourhverr_sk2/D");
    trsk2->Branch("HVshiftsk2", &gshifthv_sk2, "gshifthv_sk2/D");
    trsk2->Branch("Chi2", &rchi2_sk2, "rchi2_sk2/D");
    trsk2->Branch("Prob", &prob_sk2, "prob_sk2/D");
    trsk2->Branch("Normalization", &norm_sk2, "norm_sk2/D");
    trsk2->Branch("Index", &beta_sk2, "beta_sk2/D");
    //trsk->Branch("Offset", &offset_sk, "offset_sk/D");
    trsk2->Branch("NormErr", &normerr_sk2, "normerr_sk2/D");
    trsk2->Branch("IndexErr", &betaerr_sk2, "betaerr_sk2/D");
    //trsk->Branch("OffsetErr", &offseterr_sk, "offseterr_sk/D");
   
    trsk3->Branch("Channel", &channelid_sk3, "channelid_sk3/I");
    trsk3->Branch("HighGainHV", &geighthv_sk3, "geighthv_sk3/D");
    trsk3->Branch("HighGainHVerr", &geighthverr_sk3, "geighthv_sk3/D");
    trsk3->Branch("LowGainHV", &gfourhv_sk3, "gfourhv_sk3/D");
    trsk3->Branch("LowGainHVerr", &gfourhverr_sk3, "gfourhv_sk3/D");
    trsk3->Branch("HVshiftsk3", &gshifthv_sk3, "gshifthv_sk3/D");
    trsk3->Branch("Chi2", &rchi2_sk3, "rchi2_sk3/D");
    trsk3->Branch("Prob", &prob_sk3, "prob_sk3/D");
    trsk3->Branch("Normalization", &norm_sk3, "norm_sk3/D");
    trsk3->Branch("Index", &beta_sk3, "beta_sk3/D");
    //trsk->Branch("Offset", &offset_sk, "offset_sk/D");
    trsk3->Branch("NormErr", &normerr_sk3, "normerr_sk3/D");
    trsk3->Branch("IndexErr", &betaerr_sk3, "betaerr_sk3/D");
    //trsk->Branch("OffsetErr", &offseterr_sk, "offseterr_sk/D");
    
    
    //Int_t sk2size = skcable[0].size();
    //Int_t sk3size = skcable[0].size();
    TGraphErrors *ghv_sk[MAXPM];
    //TGraphErrors *ghv_sk[sk3size];
   
    TF1 *fHVsk[MAXPM];
    //TF1 *fHVsk[sk3size];
    TF1 *fHVinvsk[MAXPM];
    TF1 *fHVinvskerr[MAXPM];
    //TF1 *fHVinvsk[sk3size];
    
    //fHV = new TF1("fHV", "[0]*pow(x,[1])", 1500, 2500);
    //fHVinv = new TF1("fHVinv", "pow(x/[0],1./[1])", 1500, 2500);
    //fHV = new TF1("fHV", "[0]*pow((x-[2]),[1])");
    //fHVinv = new TF1("fHVinv", "pow(x/[0],1./[1])+[2]");
    //fHV->SetTitle(Form("ID%03d %s", ipmt+1, PMTserial[ipmt].Data()));
    //fHV->SetLineColor(kMagenta);
    //fHV->SetLineStyle(2);
    //fHV->SetLineWidth(1);
    //fHV->SetParName(0, "Factor");
    //fHV->SetParName(1, "Index #beta");
    //fHV->SetParName(2, "Voltage Offset");
    //fHV->SetParameter(0, 6e-29);
    //fHV->SetParameter(1, 8.5);
    //fHV->SetParLimits(0, 0, 1e-28);
    //fHV->SetParLimits(1, 3, 10);
    //fHV->SetParameter(2, -500);
    
    for (Int_t p = 0; p < MAXPM; p++){
        ghv_sk[p] = new TGraphErrors();
        ghv_sk[p]->SetMarkerStyle(8);
        ghv_sk[p]->SetMarkerColor(2);
        ghv_sk[p]->SetMarkerSize(1);
        ghv_sk[p]->SetLineColor(15);
        ghv_sk[p]->SetLineStyle(3);
        ghv_sk[p]->GetXaxis()->SetTitle("HV [V]");
        ghv_sk[p]->GetYaxis()->SetTitle("Peak [pC]");
        if (PMTinfo[p][0] == 3 && skpeak[3][p] != 0){
            ghv_sk[p]->SetName(Form("SK2_PMT_HVscan_Cable_%06d", skcable[3][p]));
            ghv_sk[p]->SetTitle(Form("SK2_PMT_HVscan_Cable_%06d", skcable[3][p]));
        }
        else if (PMTinfo[p][0] == 4 && skpeak[3][p] != 0){
            ghv_sk[p]->SetName(Form("SK3_PMT_HVscan_Cable_%06d", skcable[3][p]));
            ghv_sk[p]->SetTitle(Form("SK3_PMT_HVscan_Cable_%06d", skcable[3][p]));
        }
        
        //fHVsk[p] = new TF1(Form("fHVsk%d", p), "[0]*pow(x-[2],[1])", 1500, 2500);
        //fHVinvsk[p] = new TF1(Form("fHVinvsk%d", p), "pow(x/[0],1./[1])+[2]", 1500, 2500);
        fHVsk[p] = new TF1(Form("fHVsk%d", p), "[0]*pow(x,[1])", 1500, 2500);
        fHVinvsk[p] = new TF1(Form("fHVinvsk%d", p), "pow(x/[0],1./[1])", 1500, 2500);
        fHVinvskerr[p] = new TF1(Form("fHVinvskerr%d", p), "pow(x/[0],1./[1])", 1500, 2500);
        //fHVsk[p]->SetLineColorAlpha(4, 0.5);
        fHVsk[p]->SetLineColor(4);
        fHVsk[p]->SetLineStyle(3);
        fHVsk[p]->SetLineWidth(1);
        fHVsk[p]->SetParName(0, "Factor");
        fHVsk[p]->SetParName(1, "Index #beta");
        //fHVsk[p]->SetParName(2, "Voltage Offset");
        fHVsk[p]->SetParameter(0, 5e-18);
        fHVsk[p]->SetParameter(1, 6);
        //fHVsk[p]->SetParLimits(0, 0, 2e-17);
        fHVsk[p]->SetParLimits(1, 3, 11);
        //fHVsk[p]->SetParameter(2, -500);
        //ghv_sk[p]->SetPoint(ifile, skhv[ifile][p], hkpeak[ifile][p]);
        //gthr_[p] = new TGraphErrors();
        //g_thr[p]->SetMarkerStyle(8);
        //g_thr[p]->SetMarkerColor(2);
        //g_thr[p]->SetLineColor(15);
        //g_thr[p]->GetXaxis()->SetTitle("Threshold [mV]");
        //g_thr[p]->GetYaxis()->SetTitle("Peak [pC]");
    }
    
    /*for (Int_t pp = 0; pp < sk3size; pp++){
        ghv_sk[p] = new TGraphErrors();
        ghv_sk[p]->SetMarkerStyle(8);
        ghv_sk[p]->SetMarkerColor(2);
        ghv_sk[p]->SetMarkerSize(1);
        ghv_sk[p]->SetLineColor(15);
        ghv_sk[p]->SetLineStyle(3);
        ghv_sk[p]->GetXaxis()->SetTitle("HV [V]");
        ghv_sk[p]->GetYaxis()->SetTitle("Peak [pC]");
        ghv_sk[p]->SetName(Form("SK_PMT_HVscan_Cable_%06d", skcable[0][p]));
        ghv_sk[p]->SetTitle(Form("SK_PMT_HVscan_Cable_%06d", skcable[0][p]));
        
        //fHVsk[p] = new TF1(Form("fHVsk%d", p), "[0]*pow(x-[2],[1])", 1500, 2500);
        //fHVinvsk[p] = new TF1(Form("fHVinvsk%d", p), "pow(x/[0],1./[1])+[2]", 1500, 2500);
        fHVsk[p] = new TF1(Form("fHVsk3%d", pp), "[0]*pow(x,[1])", 1500, 2500);
        fHVinvsk[p] = new TF1(Form("fHVinvsk3%d", pp), "pow(x/[0],1./[1])", 1500, 2500);
        //fHVsk[p]->SetLineColorAlpha(4, 0.5);
        fHVsk[p]->SetLineColor(4);
        fHVsk[p]->SetLineStyle(3);
        fHVsk[p]->SetLineWidth(1);
        fHVsk[p]->SetParName(0, "Factor");
        fHVsk[p]->SetParName(1, "Index #beta");
        //fHVsk[p]->SetParName(2, "Voltage Offset");
        fHVsk[p]->SetParameter(0, 1e-18);
        fHVsk[p]->SetParameter(1, 8);
        fHVsk[p]->SetParLimits(0, 0, 5e-18);
        fHVsk[p]->SetParLimits(1, 5, 11);
        //fHVsk[p]->SetParameter(2, -500);
        //ghv_sk[p]->SetPoint(ifile, skhv[ifile][p], hkpeak[ifile][p]);
        //gthr_[p] = new TGraphErrors();
        //g_thr[p]->SetMarkerStyle(8);
        //g_thr[p]->SetMarkerColor(2);
        //g_thr[p]->SetLineColor(15);
        //g_thr[p]->GetXaxis()->SetTitle("Threshold [mV]");
        //g_thr[p]->GetYaxis()->SetTitle("Peak [pC]");
     }*/
    
    vector<Int_t> removepoint(MAXPM, 0);
    for (Int_t ifile = 0; ifile < nfile; ifile++){
        for (Int_t p = 0; p < MAXPM; p++){
            /*if (sk2peak[ifile][p] <=0 || skhv[ifile][p] <= 0){
                outtxt1 << setw(10) << skcable[ifile][p] << setw(4) << " SK2" << setw(10) << skhv[ifile][p] << setw(10) << hvshift[ifile] << std::endl;
                continue;
            }*/
            if (PMTinfo[p][0] != 3 && PMTinfo[p][0] != 4) continue;
            ghv_sk[p]->SetPoint(ifile, skhv[ifile][p], skpeak[ifile][p]);
            ghv_sk[p]->SetPointError(ifile, 0.2, skpeakerr[ifile][p]);
            if (sksigma[ifile][p] < 0.3 || skhv[ifile][p] <= 0 || skpeak[ifile][p]<=0)
                ghv_sk[p]->RemovePoint(ifile);
            /*if (skhv[ifile][p] > 0 && skpeak[ifile][p]>0){
                ghv_sk[p]->SetPoint(ifile, skhv[ifile][p], skpeak[ifile][p]);
                ghv_sk[p]->SetPointError(ifile, 0.2, skpeakerr[ifile][p]);
            }*/
            
        }
        
        /*for (Int_t pp = 0; pp < sk3size; pp++){
            if (sk3peak[ifile][p] <=0 || skhv[ifile][p] <= 0){
                outtxt1 << setw(10) << skcable[ifile][p] << setw(4) << " SK3" << setw(10) << skhv[ifile][p] << setw(10) << hvshift[ifile] << std::endl;
                continue;
            }
            ghv_sk[p]->SetPoint(ifile, skhv[ifile][p], sk3peak[ifile][p]);
            ghv_sk[p]->SetPointError(ifile, 0.2, sk3peakerr[ifile][p]);
        }*/
        
    }
    
    Double_t targetQ = 2.884; // 1.8e7 gain
    // 1.8e7/(1e-12/1.60217657e-19) = 2.884
    Double_t targetQHPK = 2.243; // 1.4e7 gain
    //1.4e7/(1e-12/1.60217657e-19) = 2.243
    
    
    for (Int_t p = 0; p < MAXPM; p++){
        //ghv_sk[p]->GetXaxis()->SetLimits(1500,2500);
        if (PMTinfo[p][0] != 3 && PMTinfo[p][0] !=4) continue;
        if (ghv_sk[p]->GetN() < 3) {
            cout << "PMT " << p+1 << " has fewer than 3 points." << endl;
            continue;
            
        }
        ghv_sk[p]->Fit(fHVsk[p], "BQN0");
        TVirtualFitter * gfitter = TVirtualFitter::Fitter(ghv_sk[p]);
        gfitter->SetPrecision(0.5);
        TFitResultPtr fitr = ghv_sk[p]->Fit(fHVsk[p],"SB+");
        int status = int(fitr);
        std::cout << Form("Fitting SK%d cable %06d", PMTinfo[p][0] - 1, p) << std::endl;
        //ghv_sk[p]->Draw();
        //c1->Update();
        ghv_sk[p]->Write();
        
        if (PMTinfo[p][0] == 3){
            channelid_sk2 = skcable[2][p];
            norm_sk2 = fHVsk[p]->GetParameter(0);
            beta_sk2 = fHVsk[p]->GetParameter(1);
            //offset_sk = fHVsk[p]->GetParameter(2);
            normerr_sk2 = fHVsk[p]->GetParError(0);
            betaerr_sk2 = fHVsk[p]->GetParError(1);
            //offseterr_sk = fHVsk[p]->GetParError(2);
            rchi2_sk2 = fHVsk[p]->GetChisquare()/fHVsk[p]->GetNDF();
            prob_sk2 = fHVsk[p]->GetProb();
            fHVinvsk[p]->SetParameter(0, fHVsk[p]->GetParameter(0));
            fHVinvsk[p]->SetParameter(1, fHVsk[p]->GetParameter(1));
            fHVinvskerr[p]->SetParameter(0, fHVsk[p]->GetParameter(0)-fHVsk[p]->GetParError(0));
            fHVinvskerr[p]->SetParameter(1, fHVsk[p]->GetParameter(1)-fHVsk[p]->GetParError(1));
            //fHVinvsk[p]->SetParameter(2, fHVsk[p]->GetParameter(2));
            geighthv_sk2 = fHVinvsk[p]->Eval(targetQ);
            gfourhv_sk2 = fHVinvsk[p]->Eval(targetQHPK);
            geighthverr_sk2 = fHVinvskerr[p]->Eval(targetQ) - geighthv_sk2;
            gfourhverr_sk2 = fHVinvskerr[p]->Eval(targetQHPK) - gfourhv_sk2;
            gshifthv_sk2 = geighthv_sk2 - SKPMThv[p];
            trsk2->Fill();
            
            if (status == 4){
                fHVsk[p]->SetLineColor(2);
                outtxt2 << setw(10) << skcable[2][p] << setw(4) << " SK2" << setw(8) << Form("%1.2f", gfourhv_sk2) << setw(8) << Form("%1.2f", geighthv_sk2) << setw(25) << Form("%1.2e(%1.2e)", norm_sk2, normerr_sk2) << setw(20) << Form("%1.2f(%1.2f)", beta_sk2, betaerr_sk2) << setw(15) << Form("%1.3f", rchi2_sk2 ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "   Failed_to_find_minimum" << "\n";
                
                if (rchi2_sk2 < 10 && rchi2_sk2 > 0 && prob_sk3 > 0.0001){
                    outtxt3 << setw(10) << skcable[2][p] << setw(4) << " SK2" << setw(15) << Form("%1.3f", rchi2_sk2 ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << Form("%1.2e(%1.2e)", norm_sk2, normerr_sk2) << setw(20) << Form("%1.2f(%1.2f)", beta_sk2, betaerr_sk2) << "\n";
                }
            }
            else if (rchi2_sk2 > 10 && prob_sk2 < 0.0001){
                fHVsk[p]->SetLineColor(6);
                outtxt2 << setw(10) << skcable[2][p] << setw(4) << " SK2" << setw(8) << Form("%1.2f", gfourhv_sk2) << setw(8) << Form("%1.2f", geighthv_sk2) << setw(25) << Form("%1.2e(%1.2e)", norm_sk2, normerr_sk2) << setw(20) << Form("%1.2f(%1.2f)", beta_sk2, betaerr_sk2) << setw(15) << Form("%1.3f", rchi2_sk2 ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                Chi2_>>_1" << "\n";
            }
            else if (rchi2_sk2 == 0 && prob_sk2 == 1){
                fHVsk[p]->SetLineColor(3);
                outtxt2 << setw(10) << skcable[2][p] << setw(4) << " SK2" << setw(8) << Form("%1.2f", gfourhv_sk2) << setw(8) << Form("%1.2f", geighthv_sk2) << setw(25) << Form("%1.2e(%1.2e)", norm_sk2, normerr_sk2) << setw(20) << Form("%1.2f(%1.2f)", beta_sk2, betaerr_sk2) << setw(15) << Form("%1.3f", rchi2_sk2 ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                 Chi2_=_0" << "\n";
            }
        }
        
        
        else if (PMTinfo[p][0] == 4){
            channelid_sk3 = skcable[2][p];
            norm_sk3 = fHVsk[p]->GetParameter(0);
            beta_sk3 = fHVsk[p]->GetParameter(1);
            //offset_sk = fHVsk[p]->GetParameter(2);
            normerr_sk3 = fHVsk[p]->GetParError(0);
            betaerr_sk3 = fHVsk[p]->GetParError(1);
            //offseterr_sk = fHVsk[p]->GetParError(2);
            rchi2_sk3 = fHVsk[p]->GetChisquare()/fHVsk[p]->GetNDF();
            prob_sk3 = fHVsk[p]->GetProb();
            fHVinvsk[p]->SetParameter(0, fHVsk[p]->GetParameter(0));
            fHVinvsk[p]->SetParameter(1, fHVsk[p]->GetParameter(1));
            //fHVinvsk[p]->SetParameter(2, fHVsk[p]->GetParameter(2));
            geighthv_sk3 = fHVinvsk[p]->Eval(targetQ);
            gfourhv_sk3 = fHVinvsk[p]->Eval(targetQHPK);
            geighthverr_sk3 = fHVinvskerr[p]->Eval(targetQ) - geighthv_sk3;
            gfourhverr_sk3 = fHVinvskerr[p]->Eval(targetQHPK) - gfourhv_sk3;
            gshifthv_sk3 = geighthv_sk3 - SKPMThv[p];
            trsk3->Fill();
            
            if (status == 4){
                fHVsk[p]->SetLineColor(2);
                outtxt2 << setw(10) << skcable[2][p] << setw(4) << " SK3" << setw(8) << Form("%1.2f", gfourhv_sk3) << setw(8) << Form("%1.2f", geighthv_sk3) << setw(25) << Form("%1.2e(%1.2e)", norm_sk3, normerr_sk3) << setw(20) << Form("%1.2f(%1.2f)", beta_sk3, betaerr_sk3) << setw(15) << Form("%1.3f", rchi2_sk3) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "   Failed_to_find_minimum" << "\n";
                
                if (rchi2_sk3 < 10 && rchi2_sk3 >0 && prob_sk3 > 0.0001){
                    outtxt3 << setw(10) << skcable[2][p] << setw(4) << " SK3" << setw(15) << Form("%1.3f", rchi2_sk3) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << Form("%1.2e(%1.2e)", norm_sk3, normerr_sk3) << setw(20) << Form("%1.2f(%1.2f)", beta_sk3, betaerr_sk3) << "\n";
                }
            }
            else if (rchi2_sk3 > 10 && prob_sk3 < 0.0001){
                fHVsk[p]->SetLineColor(6);
                outtxt2 << setw(10) << skcable[2][p] << setw(4) << " SK3" << setw(8) << Form("%1.2f", gfourhv_sk3) << setw(8) << Form("%1.2f", geighthv_sk3) << setw(25) << Form("%1.2e(%1.2e)", norm_sk3, normerr_sk3) << setw(20) << Form("%1.2f(%1.2f)", beta_sk3, betaerr_sk3) << setw(15) << Form("%1.3f", rchi2_sk3) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                Chi2_>>_1" << "\n";
            }
            
            else if (rchi2_sk3 == 0 && prob_sk3 == 1){
                fHVsk[p]->SetLineColor(3);
                outtxt2 << setw(10) << skcable[2][p] << setw(4) << " SK3" << setw(8) << Form("%1.2f", gfourhv_sk3) << setw(8) << Form("%1.2f", geighthv_sk3) << setw(25) << Form("%1.2e(%1.2e)", norm_sk3, normerr_sk3) << setw(20) << Form("%1.2f(%1.2f)", beta_sk3, betaerr_sk3) << setw(15) << Form("%1.3f", rchi2_sk3) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                 Chi2_=_0" << "\n";
            }
        }
        
       
        //mghv_sk->Add(ghv_sk[p]);
        //mghv_sk->Add(fHVsk[p]);
        //mgthr_sk->Add(gthr_sk[id]);
        //fHV->Clear();
    }
    
    /*for (Int_t pp = 0; pp < sk3size; pp++){
        //ghv_sk[p]->GetXaxis()->SetLimits(1500,2500);
        ghv_sk[p]->Fit(fHVsk[p], "BQN0");
        TVirtualFitter *gfitter = TVirtualFitter::Fitter(ghv_sk[p]);
        gfitter->SetPrecision(0.1);
        TFitResultPtr fitr = ghv_sk[p]->Fit(fHVsk[p],"SB+");
        int status = int(fitr);
        std::cout << Form("Fitting SK 3 cable %06d", pp) << std::endl;
        //ghv_sk[p]->Draw();
        //c1->Update();
        ghv_sk[p]->Write();
        channelid_sk3 = skcable[2][p];
        norm_sk3 = fHVsk[p]->GetParameter(0);
        beta_sk3 = fHVsk[p]->GetParameter(1);
        //offset_sk = fHVsk[p]->GetParameter(2);
        normerr_sk3 = fHVsk[p]->GetParError(0);
        betaerr_sk3 = fHVsk[p]->GetParError(1);
        //offseterr_sk = fHVsk[p]->GetParError(2);
        rchi2_sk3 = fHVsk[p]->GetChisquare()/fHVsk[p]->GetNDF();
        prob_sk3 = fHVsk[p]->GetProb();
        fHVinvsk[p]->SetParameter(0, fHVsk[p]->GetParameter(0));
        fHVinvsk[p]->SetParameter(1, fHVsk[p]->GetParameter(1));
        //fHVinvsk[p]->SetParameter(2, fHVsk[p]->GetParameter(2));
        geighthv_sk3 = fHVinvsk[p]->Eval(targetQ);
        gfourhv_sk3 = fHVinvsk[p]->Eval(targetQHPK);
        gshifthv_sk3 = geighthv_sk3 - SKPMThv3[p];
        trsk3->Fill();
        
        if (status == 4){
            fHVsk[p]->SetLineColor(2);
            outtxt2 << setw(10) << skcable[2][p] << setw(4) << " SK3" << setw(8) << Form("%1.2f", gfourhv_sk3) << setw(8) << Form("%1.2f", geighthv_sk3) << setw(25) << Form("%1.2e(%1.2e)", norm_sk3, normerr_sk3) << setw(20) << Form("%1.2f(%1.2f)", beta_sk3, betaerr_sk3) << setw(15) << Form("%1.3f", rchi2_sk3) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "   Failed_to_find_minimum" << "\n";
            
            if (rchi2_sk3 < 10 && rchi2_sk3 >0 && prob_sk3 > 0.0001){
                outtxt3 << setw(10) << skcable[2][p] << setw(4) << " SK3" << setw(15) << Form("%1.3f", rchi2_sk3) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << Form("%1.2e(%1.2e)", norm_sk3, normerr_sk3) << setw(20) << Form("%1.2f(%1.2f)", beta_sk3, betaerr_sk3) << "\n";
            }
        }
        else if (rchi2_sk3 > 10 && prob_sk3 < 0.0001){
            fHVsk[p]->SetLineColor(6);
            outtxt2 << setw(10) << skcable[2][p] << setw(4) << " SK3" << setw(8) << Form("%1.2f", gfourhv_sk3) << setw(8) << Form("%1.2f", geighthv_sk3) << setw(25) << Form("%1.2e(%1.2e)", norm_sk3, normerr_sk3) << setw(20) << Form("%1.2f(%1.2f)", beta_sk3, betaerr_sk3) << setw(15) << Form("%1.3f", rchi2_sk3) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                Chi2_>>_1" << "\n";
        }
        
        else if (rchi2_sk3 == 0 && prob_sk3 == 1){
            fHVsk[p]->SetLineColor(3);
            outtxt2 << setw(10) << skcable[2][p] << setw(4) << " SK3" << setw(8) << Form("%1.2f", gfourhv_sk3) << setw(8) << Form("%1.2f", geighthv_sk3) << setw(25) << Form("%1.2e(%1.2e)", norm_sk3, normerr_sk3) << setw(20) << Form("%1.2f(%1.2f)", beta_sk3, betaerr_sk3) << setw(15) << Form("%1.3f", rchi2_sk3) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                 Chi2_=_0" << "\n";
        }
        //mghv_sk->Add(ghv_sk[p]);
        //mghv_sk->Add(fHVsk[p]);
        //mgthr_sk->Add(gthr_sk[id]);
        //fHV->Clear();
    }*/

    
    c1->Print("SK_SPE_HV.pdf[");
    //mghv_sk->Draw();
    fHVsk[0]->Draw();
    for (Int_t p = 1; p < MAXPM; p++){
        if (PMTinfo[p][0] != 3 && PMTinfo[p][0] != 4) continue;
        fHVsk[p]->Draw("same");
    }
    
    c1->Update();
    c1->Modified();
    c1->Print("SK_SPE_HV.pdf");
    c1->Print("SK_SPE_HV.pdf]");

    TPaveText *t = new TPaveText(.4,.15,.9,.25,"NDC");
    TPaveText *trun = new TPaveText(.4,.25,.9,.35,"NDC");
    t->SetFillStyle(0);
    t->SetBorderSize(0);
    t->SetTextColor(kRed);
    trun->SetFillStyle(0);
    trun->SetBorderSize(0);
    trun->SetTextColor(kRed);
    
    c1->Print("SK_HV_Curves.pdf[");
    //mghv_sk->Draw();
    for (Int_t p = 0; p < MAXPM; p++){
        if (PMTinfo[p][0] != 3 && PMTinfo[p][0] != 4) continue;
        ghv_sk[p]->Draw("AP");
        /*std::cout << skcable[0][p] << std::endl;
        t->Clear();
        if (fHVinvsk[p]->Eval(targetQ) < 2500 && fHVinvskerr[p]->Eval(targetQ) < 2500){
            t->AddText(Form("HV at 1.4e7 gain : %4.2f#pm%4.2f [V]", fHVinvsk[p]->Eval(targetQHPK), fHVinvskerr[p]->Eval(targetQHPK) - fHVinvsk[p]->Eval(targetQHPK)));
            ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlue);
            t->AddText(Form("HV at 1.8e7 gain   : %4.2f#pm%4.2f [V]", fHVinvsk[p]->Eval(targetQ),fHVinvskerr[p]->Eval(targetQ) - fHVinvsk[p]->Eval(targetQ)));
            ((TText*)t->GetListOfLines()->Last())->SetTextColor(kMagenta);
            t->Draw("same");
            trun->Draw("same");
            c1->Update();
        }
        
        TPaveStats *ps = (TPaveStats*)ghv_sk[p]->GetListOfFunctions()->FindObject("stats");
        ps->SetX1NDC(0.1);
        ps->SetX2NDC(0.4);*/
        
        c1->Modified();
        c1->Update();
        c1->Print("SK_HV_Curves.pdf");
    }
    c1->Print("SK_HV_Curves.pdf]");

   
    trsk2->Write();
    trsk3->Write();
   
    fout->Close();
    outtxt1.close();
    outtxt2.close();
    outtxt3.close();
}
