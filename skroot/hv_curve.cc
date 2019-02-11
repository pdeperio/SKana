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
#include <TROOT.h>
#include <algorithm>
#include <iterator>
using namespace std;

TString outdir = "hv_ana";

#ifdef __CINT__
void hv_curve(TString PMTtype = "", TString InputDir = "hv_ana"){

#else

TString PMTtype = "";
TString InputDir = outdir;
const int MAXRANGE = 99999999;
int PlotRange[2] = {0, MAXRANGE};

int getArgs(int argc, char* argv[]);

int main(int argc, char *argv[]) {

    int args = getArgs(argc, argv);
    if(args != 0){
        std::cerr << "Usage " << std::endl;
        return 0;
    }
#endif
  
    outdir += PMTtype + "/";

    enum AnalyzeEnum {all, sk, hk};
    int AnalyzeWhat = all;
    if (PMTtype.Contains("sk")) AnalyzeWhat = sk;
    if (PMTtype.Contains("hk")) AnalyzeWhat = hk;

    gStyle->SetFrameBorderMode(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetPadColor(0);
  
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptTitle(kTRUE);
    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptFit(1111);
    //gStyle->SetPalette(kCool);
    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    
    
    const int nfile = 7;//number of files to read in
    const int MAXPM = 11146;
    const int MAXHKPM = 136;
    const int nPMTtypes = 3;
    enum PMTtypeEnum {hkpmt, sk2pmt, sk3pmt};
    TString PMTtypeNames[nPMTtypes] = {"HK", "SK2", "SK3"};
 
    Int_t runno[] = {80282, 80278, 80275, 80269, 80265, 80263, 80254};
    Double_t hvshift[] = {-75, -50, -25, 0, 25, 50, 75};
    //Double_t threshold[] = {-0.69, -0.69, -0.69, -0.69, -0.69, -0.69, -0.69};

    Int_t PMTinfo[MAXPM][2] = {0};//flags
    Double_t SKPMThv[MAXPM] = {0};
    
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

	connect.close();
	
    } else {
      cout << "Error: Connection file not open" << endl;
      exit (-1);
    }

    TString filename = outdir+"badcables";
    if (PlotRange[1] != MAXRANGE) filename += Form("_%05d", PlotRange[0]);
    ofstream outtxt_badcables;
    outtxt_badcables.open(filename+".txt");
    outtxt_badcables << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(10) << "    HV [V]" << setw(10) << "  Shift[V]" << setw(10) << " Peak [pC]" << setw(10) << "      Qisk" << "\n";

    filename = outdir+"badfittings";
    if (PlotRange[1] != MAXRANGE) filename += Form("_%05d", PlotRange[0]);
    ofstream outtxt_badfit;
    outtxt_badfit.open(filename+".txt");
    outtxt_badfit << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(8) << " LHV [V]" << setw(8) << " HHV [V]" << setw(25) << "    Norm Factor" << setw(20) << "     Index" << setw(15) << "      Chi2" << setw(15) << "    Prob" << setw(25) << "                  Comment" << "\n";
    
    filename = outdir+"badokfitting";
    if (PlotRange[1] != MAXRANGE) filename += Form("_%05d", PlotRange[0]);
    ofstream outtxt_badokfit;
    outtxt_badokfit.open(filename+".txt");
    outtxt_badokfit << setw(10) << "    Cable#" << setw(4) << " PMT" << setw(15) << "    Chi2"  << setw(15) << "    Prob" << setw(25) << "    Norm Factor" << setw(20) << "     Index" << "\n";

    filename = outdir+"deadchannels";
    if (PlotRange[1] != MAXRANGE) filename += Form("_%05d", PlotRange[0]);
    ofstream outtxt_dead;
    outtxt_dead.open(filename+".txt");
    
    TFile *f[nfile];
    Double_t skpeak[nfile][MAXPM] = {0};
    Double_t skpeakerr[nfile][MAXPM] = {0};
    Double_t skhv[nfile][MAXPM] = {0};
    Int_t skcable[nfile][MAXPM] = {0};
    Int_t sknhit[nfile][MAXPM] = {0};
    Double_t sksigma[nfile][MAXPM] = {0};
    Double_t goodchannel[nfile][MAXPM] = {0};

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
        
        for (Int_t iPMT =  PlotRange[0]; iPMT < min(entry, PlotRange[1]); iPMT++){
            tree->GetEntry(iPMT);

		   //if (PMTinfo[chid-1][0]!=3 || PMTinfo[chid-1][0]!=4) continue;
            //if (chid == 3761 || chid == 69) continue;//the two cables don't show up in the nominal hv run

	    //std::cout << "WTF " << highv << " " << peak << " " << chid << std::endl;
      //std::cout << highv << " " << peak << " " << chid << " " << PMTinfo[chid-1] << " " << std::endl;
  
	    TString PMTtypeName = "UnknownPMTtype";
	    int iPMTtype = -1;

	    // Flag 6: HK PMT
	    if (PMTinfo[iPMT][0] == 6) {
	      iPMTtype = hk;
	      PMTtypeName = PMTtypeNames[hkpmt];
	    }
	    
	    // Flag 3: SK-2 PMT, Flag 4: SK-3 PMT, 
	    else if (PMTinfo[iPMT][0] == 3) {
	      iPMTtype = sk;
	      PMTtypeName = PMTtypeNames[sk2pmt];
	    }
	    
	    else if (PMTinfo[iPMT][0] == 4) {
	      iPMTtype = sk;
	      PMTtypeName = PMTtypeNames[sk3pmt];
	    }
	  
	    else {
	      //cout << "Unknown PMT type: " << PMTinfo[iPMT][0] << endl;
	      continue;
	    }


	    // Select only HK or SK PMTs if specified above
	    if (AnalyzeWhat == hk && iPMTtype != hk) continue;
	    else if (AnalyzeWhat == sk && iPMTtype != sk) continue;
	    
	    if (peak <= 0 || peak > 10) {
	      outtxt_badcables << setw(10) << skcable[ifile][chid-1] << setw(4) << " " << PMTtypeName.Data() << setw(10) << skhv[ifile][chid-1] << setw(10) << hvshift[ifile] << setw(10) << nhits << std::endl;

	      continue;
	    }
	    
	    double ErrScaling = 2;
	    skpeak[ifile][chid-1] = peak;
	    skhv[ifile][chid-1] = highv;
	    skcable[ifile][chid-1] = chid;
	    sksigma[ifile][chid-1] = sigma;
	    skpeakerr[ifile][chid-1] = ErrScaling*peakerr; //pow((pow(peakerr,2)+pow(0.22,2)), 1/2);
	    SKPMThv[chid-1] = PMTinfo[chid-1][1];
	    sknhit[ifile][chid-1] = nhits;
	    goodchannel[ifile][chid-1] = 1;
            
	    //cout << "Add " << PMTtypeName.Data() << " element " << nfile << " " << chid-1 << " Peak err is " << peakerr << endl;
	    //cout << skpeak[ifile][chid-1] << " " << skhv[ifile][chid-1] << " " << sksigma[ifile][chid-1] << endl << endl;            
        }
	
	f[ifile]->Close();
        std::cout << "Closed fit result file" << std::endl;
    }
    
    filename = outdir+"hvscan_parameter";
    if (PlotRange[1] != MAXRANGE) filename += Form("_%05d", PlotRange[0]);
    TFile *fout = new TFile(filename+".root", "recreate");
	    
    TTree *tr[nPMTtypes] = {0};
    for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++) {
      if (AnalyzeWhat == hk && ipmttype!=hk) continue;
      else if (AnalyzeWhat == sk && ipmttype==hk) continue;

      tr[ipmttype] = new TTree("hvscan_"+PMTtypeNames[ipmttype], PMTtypeNames[ipmttype]+" PMT HV Scan Parameter");
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
    
    TGraphErrors *ghv_sk[MAXPM] = {0};
   
    TF1 *fHVsk[MAXPM] = {0};
    TF1 *fHVinvsk[MAXPM] = {0};
    TF1 *fHVinvskerr[MAXPM] = {0};

    for (Int_t p = 0; p < MAXPM; p++){

        if (p%1000==0) cout << "Making graph and functions for channel: " << p << endl;
      
        int ipmttype = -1;
	if (PMTinfo[p][0] == 6) ipmttype = hkpmt;
	else if (PMTinfo[p][0] == 3) ipmttype = sk2pmt;
	else if (PMTinfo[p][0] == 4) ipmttype = sk3pmt;
	else {
	  //cout << "Unknown PMT type: " << PMTinfo[iPMT][0] << endl;
	  continue;
	}

	if (skcable[0][p] <= 0) {
	  outtxt_dead << p+1 << endl;
	  continue;
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
        fHVsk[p] = new TF1(Form("fHVsk%d", p), "[0]*1e-20*pow(x,[1])", 1500, 2500);
        fHVinvsk[p] = new TF1(Form("fHVinvsk%d", p), "pow(x/([0]*1e-20),1./[1])", 1500, 2500);
        fHVinvskerr[p] = new TF1(Form("fHVinvskerr%d", p), "pow(x/([0]*1e-20),1./[1])", 1500, 2500);
        //fHVsk[p]->SetLineColorAlpha(4, 0.5);
        fHVsk[p]->SetLineColor(4);
        fHVsk[p]->SetLineStyle(3);
        fHVsk[p]->SetLineWidth(1);
        fHVsk[p]->SetParName(0, "Factor");
        fHVsk[p]->SetParName(1, "Index #beta");
        //fHVsk[p]->SetParName(2, "Voltage Offset");

        if (ipmttype == hkpmt){
	          fHVsk[p]->SetParameter(0, 0.5);
	          fHVsk[p]->SetParameter(1, 6);
	      }
	
     	  else {
	          fHVsk[p]->SetParameter(0, 5);
	          fHVsk[p]->SetParameter(1, 6);
	          //fHVsk[p]->SetParLimits(0, 0, 2e-17);
	          //fHVsk[p]->SetParLimits(1, 0, 15);
	      }
	      fHVsk[p]->SetParLimits(1, 0, 15);

        //fHVsk[p]->SetParameter(2, -500);
        //ghv_sk[p]->SetPoint(ifile, skhv[ifile][p], hkpeak[ifile][p]);
        //gthr_[p] = new TGraphErrors();
        //g_thr[p]->SetMarkerStyle(8);
        //g_thr[p]->SetMarkerColor(2);
        //g_thr[p]->SetLineColor(15);
        //g_thr[p]->GetXaxis()->SetTitle("Threshold [mV]");
        //g_thr[p]->GetYaxis()->SetTitle("Peak [pC]");
    }

    for (Int_t ifile = 0; ifile < nfile; ifile++){
        for (Int_t p = 0; p < MAXPM; p++){
            /*if (sk2peak[ifile][p] <=0 || skhv[ifile][p] <= 0){
                outtxt_badcables << setw(10) << skcable[ifile][p] << setw(4) << " SK2" << setw(10) << skhv[ifile][p] << setw(10) << hvshift[ifile] << std::endl;
                continue;
            }*/

	    if (!ghv_sk[p]) continue;
	    
	    ghv_sk[p]->SetPoint(ifile, skhv[ifile][p], skpeak[ifile][p]);
	    ghv_sk[p]->SetPointError(ifile, 0.5, skpeakerr[ifile][p]);
	    
	    if (sksigma[ifile][p] < 0.1 || skhv[ifile][p] <= 0 || skpeak[ifile][p]<=0){
	      ghv_sk[p]->RemovePoint(ifile);
	      cout << Form("Removing %dth point from PMT %d", ifile+1, skcable[ifile][p]) << endl << endl;
	      
	    }
        }
    }
    
    Double_t targetQ = 2.884; // 1.8e7 gain
    // 1.8e7/(1e-12/1.60217657e-19) = 2.884
    Double_t targetQHPK = 2.243; // 1.4e7 gain
    //1.4e7/(1e-12/1.60217657e-19) = 2.243
    
    TString fitOpts = "SMBE";
    int minpoint = 4;

    for (Int_t p = 0; p < MAXPM; p++){

        if (!ghv_sk[p]) continue;

	//ghv_sk[p]->GetXaxis()->SetLimits(1500,2500);

        int npoints = ghv_sk[p]->GetN();

        if (npoints < minpoint) {
            cout << "PMT " <<  skcable[0][p] << " has fewer than " << minpoint << " points." << endl;
            continue;
        }

        int ipmttype = -1;
	if (PMTinfo[p][0]==6) ipmttype = hkpmt;
	else if (PMTinfo[p][0]==3) ipmttype = sk2pmt;
	else if (PMTinfo[p][0]==4) ipmttype = sk3pmt;

        std::cout << endl << Form("Fitting %s cable %06d", PMTtypeNames[ipmttype].Data(), p) << std::endl;

	ghv_sk[p]->Fit(fHVsk[p], "BQN0");
        TVirtualFitter * gfitter = TVirtualFitter::Fitter(ghv_sk[p]);

	gfitter->SetPrecision(1);

        TFitResultPtr fitr = ghv_sk[p]->Fit(fHVsk[p], fitOpts);
        int status = (int)fitr;

 	Double_t chisquare = fHVsk[p]->GetChisquare()/fHVsk[p]->GetNDF();

	vector<Double_t> slopes(npoints-1);

	Double_t x[npoints];
	Double_t y[npoints];
	for (Int_t point = 0; point < npoints; point++){
	  ghv_sk[p]->GetPoint(point,x[point], y[point]);
	}
	for (Int_t point = 0; point < npoints - 1; point++){
	  slopes[point] = (y[point+1]-y[point])/(x[point+1]-x[point]);
	}
	vector<int> speakpoints(2);
	speakpoints[1] = distance(slopes.begin(), max_element(slopes.begin(), slopes.end()));
	speakpoints[0] = distance(slopes.begin(), min_element(slopes.begin(), slopes.end()));

	int error = 4;
	float chi2_max = 20;
	
	if (status == error || chisquare > chi2_max) {
	  for (int ifile = 0; ifile < ghv_sk[p]->GetN(); ifile++){
	    ghv_sk[p]->SetPointError(ifile, 0.5, skpeakerr[ifile][p]*1.5);
	  }
	}

        for (Int_t pointsrm = 0; pointsrm < 2; pointsrm++){
	  if (pointsrm > npoints - minpoint) break;
          
	  cout << "  status = " << status << endl;
	  
	  if (status == error || chisquare > chi2_max) {

	    if (status == error) gfitter->SetPrecision((pointsrm+1)*15);

	    //if (ipmttype == hk) ghv_sk[p]->RemovePoint(pointsrm==0?pointsrm:(nfile-pointsrm-1));
	    ghv_sk[p]->RemovePoint(speakpoints[pointsrm]);
	    
	    fitr = ghv_sk[p]->Fit(fHVsk[p], fitOpts);
	  }
	  status = int(fitr);	  
	  chisquare = fHVsk[p]->GetChisquare()/fHVsk[p]->GetNDF();
	  
	  if (status != error && chisquare <= chi2_max) break;
        }
	
	if ((status == error || chisquare > chi2_max) && ghv_sk[p]->GetN()>minpoint){
	  if (status == error) gfitter->SetPrecision(100);

	  ghv_sk[p]->RemovePoint(0);
	  fitr = ghv_sk[p]->Fit(fHVsk[p], fitOpts);
	  status = (int)fitr;
	  chisquare = fHVsk[p]->GetChisquare()/fHVsk[p]->GetNDF();
	}
        //std::cout << Form("Fitting SK%d cable %06d", PMTinfo[p][0] - 1, p) << std::endl;
        
        ghv_sk[p]->Draw("AP");
        //c1->Update();
	fout->cd();
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
            
	if (status == error){
	  fHVsk[p]->SetLineColor(2);
	  outtxt_badfit << setw(10) << skcable[0][p] << setw(4) << " " << PMTtypeNames[ipmttype] << setw(8) << Form("%1.2f", gfourhv[ipmttype]) << setw(8) << Form("%1.2f", geighthv[ipmttype]) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "   Failed_to_find_minimum" << "\n";

	  if (rchi2[ipmttype] < chi2_max && rchi2[ipmttype] > 0 && prob[ipmttype] > 0.0001){
	    outtxt_badokfit << setw(10) << skcable[0][p] << setw(4) << " " << PMTtypeNames[ipmttype] << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << "\n";
	  }
	}

	if (rchi2[ipmttype] > chi2_max && prob[ipmttype] < 0.0001){
	    fHVsk[p]->SetLineColor(6);
	    outtxt_badfit << setw(10) << skcable[0][p] << setw(4) << " " << PMTtypeNames[ipmttype] << setw(8) << Form("%1.2f", gfourhv[ipmttype]) << setw(8) << Form("%1.2f", geighthv[ipmttype]) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                Chi2_>>_1" << "\n";
	}
	else if (rchi2[ipmttype] == 0 && prob[ipmttype] == 1){
	  fHVsk[p]->SetLineColor(3);
	  outtxt_badfit << setw(10) << skcable[0][p] << setw(4) << " " << PMTtypeNames[ipmttype] << setw(8) << Form("%1.2f", gfourhv[ipmttype]) << setw(8) << Form("%1.2f", geighthv[ipmttype]) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                 Chi2_=_0" << "\n";
	}
       
        //mghv_sk->Add(ghv_sk[p]);
        //mghv_sk->Add(fHVsk[p]);
        //mgthr_sk->Add(gthr_sk[id]);
        //fHV->Clear();
    }
  
    for (int ipmttype=0; ipmttype<nPMTtypes; ipmttype++)
      if (tr[ipmttype])
	tr[ipmttype]->Write();

    fout->Close();

    bool isDrawn = 0;
    for (Int_t p = 0; p < MAXPM; p++){

      if (!fHVsk[p]) continue;

      if (!isDrawn) {
	fHVsk[p]->Draw();
      	isDrawn = 1;
      }
      
      else 
	fHVsk[p]->Draw("same");
      
    }
    
    TString CanvasName = outdir+"SPE_HV";
    if (AnalyzeWhat==hk) CanvasName += "_HK";
    else if (AnalyzeWhat==sk) CanvasName += "_SK";
    if (PlotRange[1] != MAXRANGE) CanvasName += Form("_%05d", PlotRange[0]);
    CanvasName += ".pdf";

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
    
    CanvasName = outdir+"HV_Curves";
    if (AnalyzeWhat==hk) CanvasName += "_HK";
    else if (AnalyzeWhat==sk) CanvasName += "_SK";
    if (PlotRange[1] != MAXRANGE) CanvasName += Form("_%05d", PlotRange[0]);
    CanvasName += ".pdf";

    c1->Print(CanvasName+"[");
    
    for (Int_t p = 0; p < MAXPM; p++){
      
        if (!ghv_sk[p]) continue;
        if (ghv_sk[p]->GetN() < minpoint) continue;

        ghv_sk[p]->Draw("AP");
	
        cout << "Drawing cable " << skcable[0][p] << endl;

        t->Clear();

        //if (fHVinvsk[p]->Eval(targetQ) < 2500 && fHVinvskerr[p]->Eval(targetQ) < 2500){

	  t->AddText(Form("HV at 1.4e7 gain : %4.2f [V]", fHVinvsk[p]->Eval(targetQHPK)));
          //t->AddText(Form("HV at 1.4e7 gain : %4.2f#pm%4.2f [V]", fHVinvsk[p]->Eval(targetQHPK), fHVinvskerr[p]->Eval(targetQHPK) - fHVinvsk[p]->Eval(targetQHPK)));
	  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlue);
	  t->AddText(Form("HV at 1.8e7 gain : %4.2f [V]", fHVinvsk[p]->Eval(targetQ)));      
          //t->AddText(Form("HV at 1.8e7 gain   : %4.2f#pm%4.2f [V]", fHVinvsk[p]->Eval(targetQ),fHVinvskerr[p]->Eval(targetQ) - fHVinvsk[p]->Eval(targetQ)));
	  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kMagenta);
          t->Draw("same");
          trun->Draw("same");
          c1->Update();
       //}
        
	TPaveStats *ps = (TPaveStats*)ghv_sk[p]->GetListOfFunctions()->FindObject("stats");

	if (ps != (TPaveStats*)0){
	  //ps = (TPaveStats*)fHVsk[p]->FindObject("stats");
	  ps->SetX1NDC(0.1);
	  ps->SetX2NDC(0.5);
	  ps->SetY1NDC(0.65);
	  ps->SetY2NDC(0.9);
	}
	else {
	  cout << "No function in the graph for cable " << skcable[0][p] << " " << fHVsk[p]->GetNDF() << endl;
	}
	c1->Modified();
	c1->Update();
	c1->Print(CanvasName);
    }
    c1->Print(CanvasName+"]");
    
    outtxt_badcables.close();
    outtxt_badfit.close();
    outtxt_badokfit.close();
    outtxt_dead.close();
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

	        case 'l':
	   	    PlotRange[0] = atoi(argv[2]);
                    ++argv; --argc;
                    break;		    
	        case 'u':
	  	    PlotRange[1] = atoi(argv[2]);
                    ++argv; --argc;
                    break;
            }
            
            ++argv; --argc;
        }
        
        return 0;
        
    }
#endif
