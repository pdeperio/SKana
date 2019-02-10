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

TString outdir = "hv_ana";

#ifdef __CINT__
void hv_curve(TString PMTtype = "", TString InputDir = "hv_ana"){

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
    
    gStyle->SetFrameBorderMode(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetPadColor(0);

    gStyle->SetOptTitle(kTRUE);
    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptFit(1111);
    //gStyle->SetPalette(kCool);
    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    
    
    const int nfile = 7;//number of files to read in
    const int MAXPM = 11146;
    const int MAXHKPM = 136;
    const int nPMTtypes = 3;
    enum PMTtypeEnum {hk, sk2, sk3};
    TString PMTtypeNames[nPMTtypes] = {"HK", "SK2", "SK3"};
 
    Int_t runno[] = {80282, 80278, 80275, 80269, 80265, 80263, 80254};
    Double_t hvshift[] = {-75, -50, -25, 0, 25, 50, 75};
    //Double_t threshold[] = {-0.69, -0.69, -0.69, -0.69, -0.69, -0.69, -0.69};

    Int_t PMTinfo[MAXPM][2] = {0};//flags
    std::vector<Double_t> SKPMThv(MAXPM,0);
    
    int pmtid;
    
    ifstream connect;
    connect.open("/disk02/usr6/pdeperio/precalib_2018/connection.super.sk-5.dat");
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
    outtxt1.open(outdir+"badcables.txt");
    outtxt1 << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(10) << "    HV [V]" << setw(10) << "  Shift[V]" << setw(10) << " Peak [pC]" << setw(10) << "      Qisk" << "\n";
    
    ofstream outtxt2;
    outtxt2.open(outdir+"badfitting.txt");
    outtxt2 << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(8) << " LHV [V]" << setw(8) << " HHV [V]" << setw(25) << "    Norm Factor" << setw(20) << "     Index" << setw(15) << "      Chi2" << setw(15) << "    Prob" << setw(25) << "                  Comment" << "\n";
    
    ofstream outtxt3;
    outtxt3.open(outdir+"badokfitting.txt");
    outtxt3 << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(15) << "    Chi2"  << setw(15) << "    Prob" << setw(25) << "    Norm Factor" << setw(20) << "     Index" << "\n";
    
    TFile *f[nfile];
    
    std::vector<std::vector<Double_t> > skpeak(nfile, std::vector<Double_t>(MAXPM,0));
    std::vector<std::vector<Double_t> > skpeakerr(nfile, std::vector<Double_t>(MAXPM,0));
    std::vector<std::vector<Double_t> > skhv(nfile, std::vector<Double_t>(MAXPM,0));
    std::vector<std::vector<Int_t> > skcable(nfile, std::vector<Int_t>(MAXPM,0));
    std::vector<std::vector<Int_t> > sknhit(nfile, std::vector<Int_t>(MAXPM,0));
    std::vector<std::vector<Double_t> > sksigma(nfile, std::vector<Double_t>(MAXPM,0));
    std::vector<std::vector<Double_t> > goodchannel(nfile, std::vector<Double_t>(MAXPM,0));

    for (Int_t ifile = 0; ifile < nfile; ifile++){

        TString InputFile = Form("%s/fit_result_%d%s.root",InputDir.Data(),runno[ifile],PMTtype.Data());
        std::cout << "Reading fit result file: " << InputFile.Data() << endl;
	f[ifile] = new TFile(InputFile,"read");
	
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

	    //std::cout << "WTF " << highv << " " << peak << " " << chid << std::endl;

	    TString PMTtypeName = "UnknownPMTtype";
	    if (AnalyzeHK) {
	      if (PMTinfo[chid-1][0]==6) PMTtypeName = PMTtypeNames[hk];
	      else continue;
	      
	    } else {
	      
	      if (PMTinfo[chid-1][0]==3) PMTtypeName = PMTtypeNames[sk2];
	      else if (PMTinfo[chid-1][0]==4) PMTtypeName = PMTtypeNames[sk3];
	      else continue;
	    }
	    
	    if ( (AnalyzeHK && (peak <= 0 || peak > 10)) ||
		 (!AnalyzeHK && nhits < 10) ) {
	      outtxt1 << setw(10) << skcable[ifile][chid-1] << setw(4) << " " << PMTtypeName.Data() << setw(10) << skhv[ifile][chid-1] << setw(10) << hvshift[ifile] << setw(10) << nhits << std::endl;

	      continue;
	    }
	    
	    skpeak[ifile][chid-1] = peak;
	    skhv[ifile][chid-1] = highv;
	    skcable[ifile][chid-1] = chid;
	    sksigma[ifile][chid-1] = sigma;
	    skpeakerr[ifile][chid-1] = peakerr; //pow((pow(peakerr,2)+pow(0.22,2)), 1/2);
	    SKPMThv[chid-1] = PMTinfo[chid-1][1];
	    sknhit[ifile][chid-1] = nhits;
	    goodchannel[ifile][chid-1] = 1;
            
        }

	f[ifile]->Close();
        std::cout << "Closed fit result file" << std::endl;
    }
    
    TFile *fout = new TFile(outdir+Form("hvscan_parameter%s.root", PMTtype.Data()), "recreate");
	    
    TTree *tr[nPMTtypes] = {0};
    if (AnalyzeHK)
      tr[hk] = new TTree("hvscan_hk", PMTtypeNames[hk]+" PMT HV Scan Parameter");

    else {
      tr[sk2] = new TTree("hvscan_sk2", PMTtypeNames[sk2]+" PMT HV Scan Parameter");
      tr[sk3] = new TTree("hvscan_sk3", PMTtypeNames[sk3]+" PMT HV Scan Parameter");
    }
    
    Int_t channelid[nPMTtypes];
    Double_t geighthv[nPMTtypes];
    Double_t geighthverr[nPMTtypes];
    Double_t gfourhv[nPMTtypes];
    Double_t gfourhverr[nPMTtypes];
    Double_t gshifthv[nPMTtypes];
    Double_t rchi2[nPMTtypes];
    Double_t prob[nPMTtypes];
    Double_t norm[nPMTtypes];
    Double_t beta[nPMTtypes];
    Double_t offset[nPMTtypes] = {0};
    Double_t normerr[nPMTtypes];
    Double_t betaerr[nPMTtypes];
    Double_t offseterr[nPMTtypes] = {0};
    Double_t resolution[nPMTtypes] = {0};

    for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++) {

      if (!tr[ipmttype]) continue;
      
      tr[ipmttype]->Branch("Channel", &channelid[ipmttype], "channelid/I");
      tr[ipmttype]->Branch("HighGainHV", &geighthv[ipmttype], "geighthv/D");
      tr[ipmttype]->Branch("HighGainHVerr", &geighthverr[ipmttype], "geighthverr/D");
      tr[ipmttype]->Branch("LowGainHV", &gfourhv[ipmttype], "gfourhv/D");
      tr[ipmttype]->Branch("LowGainHVerr", &gfourhverr[ipmttype], "gfourhverr/D");
      tr[ipmttype]->Branch("HVshift", &gshifthv[ipmttype], "gshifthv/D");
      tr[ipmttype]->Branch("Chi2", &rchi2[ipmttype], "rchi2/D");
      tr[ipmttype]->Branch("Prob", &prob[ipmttype], "prob/D");
      tr[ipmttype]->Branch("Normalization", &norm[ipmttype], "norm/D");
      tr[ipmttype]->Branch("Resolution", &resolution[ipmttype], "resolution/D");
      tr[ipmttype]->Branch("Index", &beta[ipmttype], "beta/D");
      //trsk->Branch("Offset", &offset_sk, "offset_sk/D");
      tr[ipmttype]->Branch("NormErr", &normerr[ipmttype], "normerr/D");
      tr[ipmttype]->Branch("IndexErr", &betaerr[ipmttype], "betaerr/D");
      //trsk->Branch("OffsetErr", &offseterr[ipmttype], "offseterr/D");
    }
    
    //Int_t sk2size = skcable[0].size();
    //Int_t sk3size = skcable[0].size();
    TGraphErrors *ghv_sk[MAXPM] = {0};
    //TGraphErrors *ghv_sk[sk3size];
   
    TF1 *fHVsk[MAXPM] = {0};
    //TF1 *fHVsk[sk3size] = {0};
    TF1 *fHVinvsk[MAXPM] = {0};
    TF1 *fHVinvskerr[MAXPM] = {0};
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

        int ipmttype = -1;
        if (AnalyzeHK) {
	  if (PMTinfo[p][0]==6) ipmttype = hk;
	  else continue;
	      
	} else {
	  if (PMTinfo[p][0]==3) ipmttype = sk2;
	  else if (PMTinfo[p][0]==4) ipmttype = sk3;
	  else continue;
	}

	if (!goodchannel[0][p]) continue;

	ghv_sk[p] = new TGraphErrors();
        ghv_sk[p]->SetMarkerStyle(8);
        ghv_sk[p]->SetMarkerColor(2);
        ghv_sk[p]->SetMarkerSize(1);
        ghv_sk[p]->SetLineColor(15);
        ghv_sk[p]->SetLineStyle(3);
        ghv_sk[p]->GetXaxis()->SetTitle("HV [V]");
        ghv_sk[p]->GetYaxis()->SetTitle("Peak [pC]");

	ghv_sk[p]->SetName(PMTtypeNames[ipmttype]+Form("_PMT_HVscan_Cable_%06d", skcable[0][p]));
	ghv_sk[p]->SetTitle(PMTtypeNames[ipmttype]+Form(" PMT HV scan (Cable %06d)", skcable[0][p]));
        
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
        if (!AnalyzeHK) fHVsk[p]->SetParLimits(1, 3, 11);
        //fHVsk[p]->SetParameter(2, -500);
        //ghv_sk[p]->SetPoint(ifile, skhv[ifile][p], hkpeak[ifile][p]);
        //gthr_[p] = new TGraphErrors();
        //g_thr[p]->SetMarkerStyle(8);
        //g_thr[p]->SetMarkerColor(2);
        //g_thr[p]->SetLineColor(15);
        //g_thr[p]->GetXaxis()->SetTitle("Threshold [mV]");
        //g_thr[p]->GetYaxis()->SetTitle("Peak [pC]");
    }

    vector<Int_t> removepoint(MAXPM, 0);
    for (Int_t ifile = 0; ifile < nfile; ifile++){
        for (Int_t p = 0; p < MAXPM; p++){
            /*if (sk2peak[ifile][p] <=0 || skhv[ifile][p] <= 0){
                outtxt1 << setw(10) << skcable[ifile][p] << setw(4) << " SK2" << setw(10) << skhv[ifile][p] << setw(10) << hvshift[ifile] << std::endl;
                continue;
            }*/

	    if (!ghv_sk[p]) continue;
	    
	    ghv_sk[p]->SetPoint(ifile, skhv[ifile][p], skpeak[ifile][p]);

            if (AnalyzeHK) {
	      ghv_sk[p]->SetPointError(ifile, 0.5, skpeakerr[ifile][p]);
	      if (sksigma[ifile][p] < 0.1)
		ghv_sk[p]->RemovePoint(ifile);
	      
	    } else {
	      ghv_sk[p]->SetPointError(ifile, 0.2, skpeakerr[ifile][p]);
	      if (sksigma[ifile][p] < 0.3 || skhv[ifile][p] <= 0 || skpeak[ifile][p]<=0)
                ghv_sk[p]->RemovePoint(ifile);
	    }
        }
    }
    
    Double_t targetQ = 2.884; // 1.8e7 gain
    // 1.8e7/(1e-12/1.60217657e-19) = 2.884
    Double_t targetQHPK = 2.243; // 1.4e7 gain
    //1.4e7/(1e-12/1.60217657e-19) = 2.243
    
    TString fitOpts = "SMBE+";

    for (Int_t p = 0; p < MAXPM; p++){
        if (!ghv_sk[p]) continue;
	//ghv_sk[p]->GetXaxis()->SetLimits(1500,2500);

        if (ghv_sk[p]->GetN() < 3) {
            cout << "PMT " << p+1 << " has fewer than 3 points." << endl;
            continue;
        }

        int ipmttype = -1;
	if (PMTinfo[p][0]==6) ipmttype = hk;
	else if (PMTinfo[p][0]==3) ipmttype = sk2;
	else if (PMTinfo[p][0]==4) ipmttype = sk3;

        std::cout << endl << Form("Fitting %s cable %06d", PMTtypeNames[ipmttype].Data(), p) << std::endl;

	ghv_sk[p]->Fit(fHVsk[p], "BQN0");
        TVirtualFitter * gfitter = TVirtualFitter::Fitter(ghv_sk[p]);

	if (AnalyzeHK) gfitter->SetPrecision(1);
	else gfitter->SetPrecision(0.5);
	
        TFitResultPtr fitr = ghv_sk[p]->Fit(fHVsk[p], fitOpts);
        int status = (int)fitr;
        for (Int_t pointsrm = 0; pointsrm < 2; pointsrm++){
	  status = (int)fitr;
	    cout << "Status = " << status << endl;
	    if (status == 4) {
                ghv_sk[p]->RemovePoint(pointsrm==0?pointsrm:(nfile-pointsrm-1));
                fitr = ghv_sk[p]->Fit(fHVsk[p], fitOpts);
            }
	    else break;
        }
	
        //ghv_sk[p]->Draw();
        //c1->Update();
        ghv_sk[p]->Write();
        
	channelid[ipmttype] = skcable[0][p];
	norm[ipmttype] = fHVsk[p]->GetParameter(0);
	beta[ipmttype] = fHVsk[p]->GetParameter(1);
	//offset_sk = fHVsk[p]->GetParameter(2);
	normerr[ipmttype] = fHVsk[p]->GetParError(0);
	betaerr[ipmttype] = fHVsk[p]->GetParError(1);
	//offseterr_sk = fHVsk[p]->GetParError(2);
	rchi2[ipmttype] = fHVsk[p]->GetChisquare()/fHVsk[p]->GetNDF();
	prob[ipmttype] = fHVsk[p]->GetProb();
	fHVinvsk[p]->SetParameter(0, fHVsk[p]->GetParameter(0));
	fHVinvsk[p]->SetParameter(1, fHVsk[p]->GetParameter(1));
	fHVinvskerr[p]->SetParameter(0, fHVsk[p]->GetParameter(0)-fHVsk[p]->GetParError(0));
	fHVinvskerr[p]->SetParameter(1, fHVsk[p]->GetParameter(1)-fHVsk[p]->GetParError(1));
	//fHVinvsk[p]->SetParameter(2, fHVsk[p]->GetParameter(2));
	geighthv[ipmttype] = fHVinvsk[p]->Eval(targetQ);
	gfourhv[ipmttype] = fHVinvsk[p]->Eval(targetQHPK);
	geighthverr[ipmttype] = fHVinvskerr[p]->Eval(targetQ) - geighthv[ipmttype];
	gfourhverr[ipmttype] = fHVinvskerr[p]->Eval(targetQHPK) - gfourhv[ipmttype];
	gshifthv[ipmttype] = geighthv[ipmttype] - SKPMThv[p];
        resolution[ipmttype] = sksigma[2][p]/skpeak[2][p];
	tr[ipmttype]->Fill();
            
	if (status == 4){
	  fHVsk[p]->SetLineColor(2);
	  outtxt2 << setw(10) << skcable[0][p] << setw(4) << " " << PMTtypeNames[ipmttype] << setw(8) << Form("%1.2f", gfourhv[ipmttype]) << setw(8) << Form("%1.2f", geighthv[ipmttype]) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "   Failed_to_find_minimum" << "\n";

	  if (!AnalyzeHK) {
	    if (rchi2[ipmttype] < 10 && rchi2[ipmttype] > 0 && prob[ipmttype] > 0.0001){
	      outtxt3 << setw(10) << skcable[0][p] << setw(4) << " SK2" << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << "\n";
	    }
	  }
	}

	else if (!AnalyzeHK) {
	  if (rchi2[ipmttype] > 10 && prob[ipmttype] < 0.0001){
	    fHVsk[p]->SetLineColor(6);
	    outtxt2 << setw(10) << skcable[0][p] << setw(4) << " " << PMTtypeNames[ipmttype] << setw(8) << Form("%1.2f", gfourhv[ipmttype]) << setw(8) << Form("%1.2f", geighthv[ipmttype]) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                Chi2_>>_1" << "\n";
	  }
	  else if (rchi2[ipmttype] == 0 && prob[ipmttype] == 1){
	    fHVsk[p]->SetLineColor(3);
	    outtxt2 << setw(10) << skcable[0][p] << setw(4) << " " << PMTtypeNames[ipmttype] << setw(8) << Form("%1.2f", gfourhv[ipmttype]) << setw(8) << Form("%1.2f", geighthv[ipmttype]) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                 Chi2_=_0" << "\n";
	  }
	}
       
        //mghv_sk->Add(ghv_sk[p]);
        //mghv_sk->Add(fHVsk[p]);
        //mgthr_sk->Add(gthr_sk[id]);
        //fHV->Clear();
    }
    
    bool isDrawn = 0;
    for (Int_t p = 0; p < MAXPM; p++){

      if (!fHVsk[p]) continue;

      if (!isDrawn)
	fHVsk[p]->Draw();
      
      else {
	fHVsk[p]->Draw("same");
	isDrawn = 1;
      }
    }
    
    TString CanvasName = "_SPE_HV.pdf";
    if (AnalyzeHK) CanvasName = outdir+"HK"+CanvasName;
    else CanvasName = outdir+"SK"+CanvasName;

    c1->Update();
    c1->Modified();
    c1->Print(CanvasName);

    TPaveText *t = new TPaveText(.4,.15,.9,.25,"NDC");
    TPaveText *trun = new TPaveText(.4,.25,.9,.35,"NDC");
    t->SetFillStyle(0);
    t->SetBorderSize(0);
    t->SetTextColor(kRed);
    trun->SetFillStyle(0);
    trun->SetBorderSize(0);
    trun->SetTextColor(kRed);
    
    CanvasName = "_HV_Curves.pdf";
    if (AnalyzeHK) CanvasName = outdir+"HK"+CanvasName;
    else CanvasName = outdir+"SK"+CanvasName;
    c1->Print(CanvasName+"[");
    
    //mghv_sk->Draw();
    for (Int_t p = 0; p < MAXPM; p++){
      
        if (!ghv_sk[p]) continue;

        ghv_sk[p]->Draw("AP");
	
        //std::cout << skcable[0][p] << std::endl;
        t->Clear();

        t->AddText(Form("HV at 1.4e7 gain : %4.2f#pm%4.2f [V]", fHVinvsk[p]->Eval(targetQHPK), fHVinvskerr[p]->Eval(targetQHPK) - fHVinvsk[p]->Eval(targetQHPK)));
        ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlue);
        t->AddText(Form("HV at 1.8e7 gain   : %4.2f#pm%4.2f [V]", fHVinvsk[p]->Eval(targetQ),fHVinvskerr[p]->Eval(targetQ) - fHVinvsk[p]->Eval(targetQ)));
        ((TText*)t->GetListOfLines()->Last())->SetTextColor(kMagenta);
        t->Draw("same");
        trun->Draw("same");
        c1->Update();

	/*
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
        
        TPaveStats *ps = (TPaveStats*)ghv_sk[p]->GetListOfFunctions()->FindObject("stats");
        ps->SetX1NDC(0.1);
        ps->SetX2NDC(0.4);
        c1->Modified();
        c1->Update();
        c1->Print(CanvasName);
    }
    c1->Print(CanvasName+"]");

    for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++)
      if (tr[ipmttype])
	tr[ipmttype]->Write();
    
    fout->Close();
    outtxt1.close();
    outtxt2.close();
    outtxt3.close();
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
