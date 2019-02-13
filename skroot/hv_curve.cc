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

    gErrorIgnoreLevel = kWarning; // For removing TCanvas::Print msgs

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

    /*TString KEKdir = "const";
    ifstream intxt;
    std::vector<double> KEKcable;
    intxt.open(KEKdir+"/pmt_prod_year.dat");
    if (intxt.is_open()){
      while (!intxt.eof()){
	Int_t cableid, w, x, y, z;
	double KEK; 
	
	intxt >> cableid >> KEK >> w >> x >> y >> z;
	if (KEK > 0) continue;
	KEKcable.push_back(cableid);
      }

      intxt.close();

    } else {
      cout << "Error: KEK PMT file not open" << endl;
      exit (-1);
      }*/
    

    
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
    outtxt_dead << setw(10) << "   Channel" << setw(6) << "   PMT"<< "\n";

    filename = outdir+"targethv";
    if (PlotRange[1] != MAXRANGE) filename += Form("_%05d", PlotRange[0]);
    ofstream outtxt_target;
    outtxt_target.open(filename+".txt");
    outtxt_target << "    Cable#" << setw(4) << " PMT" << "\n";

    //TFile *fsyscheck = new TFile(Form("%s/fit_result_80249%s.root",InputDir.Data(),PMTtype.Data()));
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
	//TTree *systree = (TTree*)(fsyscheck->Get("spe"));
	//Int_t sysentry = systree->GetEntries();
        TTree *tree = (TTree*)(f[ifile]->Get("spe"));
        Int_t entry = tree->GetEntries();
        //std::cout << entry << " " << sysentry << std::endl;
        //c1->Print(Form("skcable_%d.pdf[", ifile));
        //TH1D *hskcable = new TH1D("hskcable", "hskcable", 100, 0, 11200);
        
        Double_t peak;
        Double_t peakerr;
        Double_t highv;
        Double_t sigma;
        Int_t chid;
        Int_t nhits;

	//systree->SetBranchAddress("Peak", &syspeak);
	//systree->SetBranchAddress("Channel", &syschid);

	//vector<Double_t> syspeakvec;
	//vector<Int_t> syschidvec;
	
	//for (int isys = PlotRange[0]; isys < min(sysentry, PlotRange[1]); isys++){
	//systree->GetEntry(isys);
	//syspeakvec.push_back(syspeak);
	//syschidvec.push_back(syschid);
	//}
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
	    if (PMTinfo[chid-1][0] == 6) {
	      iPMTtype = hk;
	      PMTtypeName = PMTtypeNames[hkpmt];
	    }
	    
	    // Flag 3: SK-2 PMT, Flag 4: SK-3 PMT, 
	    else if (PMTinfo[chid-1][0] == 3) {
	      iPMTtype = sk;
	      PMTtypeName = PMTtypeNames[sk2pmt];
	    }
	    
	    else if (PMTinfo[chid-1][0] == 4) {
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

	    double ErrScaling = 1;
	    //vector<int>::iterator isysvec = (find(syschidvec.begin(), syschidvec.end(), chid));
	    //if (isysvec != syschidvec.end()){
	    //int sysindex = distance(syschidvec.begin(), isysvec);
	    //skpeakerr[ifile][chid-1] = sqrt(pow(peakerr,2)+pow(0.11,2));
	    //}
	    //skpeakerr[ifile][chid-1] = ErrScaling*peakerr;
	    skpeak[ifile][chid-1] = peak;
	    ErrScaling = 1./peak; 
	    //if (iPMTtype == hk)	ErrScaling = 1./peak;
	    //else if (iPMTtype == sk) ErrScaling = 1./peak;
	    //skpeakerr[ifile][chid-1] = ErrScaling*peakerr;       
	    skhv[ifile][chid-1] = highv;
	    skcable[ifile][chid-1] = chid;
	    sksigma[ifile][chid-1] = sigma;
	    SKPMThv[chid-1] = PMTinfo[chid-1][1];
	    sknhit[ifile][chid-1] = nhits;
	    skpeakerr[ifile][chid-1] = ErrScaling*(sqrt(pow(peakerr,2)+pow(0.055,2)));
	    //skpeakerr[ifile][chid-1] = ErrScaling*peakerr;
	    goodchannel[ifile][chid-1] = 1;
            
	    //cout << "Add " << PMTtypeName.Data() << " element " << nfile << " " << chid-1 << " Peak err is " << peakerr << endl;
	    //cout << skpeak[ifile][chid-1] << " " << skhv[ifile][chid-1] << " " << sksigma[ifile][chid-1] << endl << endl;            
        }
	
	f[ifile]->Close();
        std::cout << "Closed fit result file" << std::endl;
    }
    //fsyscheck->Close();
    //std::cout << "Closed Systematic Err Reference File" << std::endl;
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
    TGraphErrors *ghv_orig[MAXPM] = {0};

    TF1 *fHVsk[MAXPM] = {0};
    TF1 *fHVinvsk[MAXPM] = {0};
    TF1 *fHVinvskerr[MAXPM] = {0};

    for (Int_t p = PlotRange[0]; p < min(MAXPM, PlotRange[1]); p++){

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
	  outtxt_dead << setw(10) << p+1 << setw(6) << PMTtypeNames[ipmttype] << endl;
	  cout << p+1 << " " << skcable[0][p] << " " << skhv[0][p] << " " << skpeak[0][p] << endl;
	  continue;
	}

	if (!goodchannel[0][p]) continue;
	
	ghv_sk[p] = new TGraphErrors();

        ghv_sk[p]->SetMarkerStyle(8);
        ghv_sk[p]->SetMarkerColor(2);
        ghv_sk[p]->SetMarkerSize(1);
        ghv_sk[p]->SetLineColor(15);
        ghv_sk[p]->SetLineStyle(3);
        ghv_sk[p]->GetXaxis()->SetTitle("ln(HV[V])");
        ghv_sk[p]->GetYaxis()->SetTitle("ln(Peak[pC])");

	ghv_sk[p]->SetName(PMTtypeNames[ipmttype]+Form("_PMT_HVscan_Cable_%06d", skcable[0][p]));
	ghv_sk[p]->SetTitle(PMTtypeNames[ipmttype]+Form(" PMT HV scan (Cable %06d)", skcable[0][p]));

        //fHVsk[p] = new TF1(Form("fHVsk%d", p), "[0]*pow(x-[2],[1])", 1500, 2500);
        //fHVinvsk[p] = new TF1(Form("fHVinvsk%d", p), "pow(x/[0],1./[1])+[2]", 1500, 2500);
        fHVsk[p] = new TF1(Form("fHVsk%d", p), "[0]+x*[1]", 7, 8);
        fHVinvsk[p] = new TF1(Form("fHVinvsk%d", p), "(x-[0])*1.0/[1]", 0, 2.7);
        fHVinvskerr[p] = new TF1(Form("fHVinvskerr%d", p), "(x-[0])*1.0/[1]", 0, 2.7);
        //fHVsk[p]->SetLineColorAlpha(4, 0.5);
        fHVsk[p]->SetLineColor(4);
        fHVsk[p]->SetLineStyle(3);
        fHVsk[p]->SetLineWidth(1);
        fHVsk[p]->SetParName(0, "Factor");
        fHVsk[p]->SetParName(1, "Index #beta");
        //fHVsk[p]->SetParName(2, "Voltage Offset");

        if (ipmttype == hkpmt){
	          fHVsk[p]->SetParameter(0, -45);
	          fHVsk[p]->SetParameter(1, 6);
	      }
	
     	  else {
	          fHVsk[p]->SetParameter(0, -45);
	          fHVsk[p]->SetParameter(1, 8);
	          //fHVsk[p]->SetParLimits(0, 0, 2e-17);
	          //fHVsk[p]->SetParLimits(1, 0, 15);
	      }
	//fHVsk[p]->SetParLimits(0, -10, 10);
	//fHVsk[p]->SetParLimits(1, 3, 15);
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
      for (Int_t p = PlotRange[0]; p < min(MAXPM, PlotRange[1]); p++){
            /*if (sk2peak[ifile][p] <=0 || skhv[ifile][p] <= 0){
                outtxt_badcables << setw(10) << skcable[ifile][p] << setw(4) << " SK2" << setw(10) << skhv[ifile][p] << setw(10) << hvshift[ifile] << std::endl;
                continue;
            }*/

	    if (skcable[ifile][p] <= 0) continue; 

	    if (!ghv_sk[p]) continue;
	    
	    ghv_sk[p]->SetPoint(ifile, TMath::Log(skhv[ifile][p]), TMath::Log(skpeak[ifile][p]));
	    ghv_sk[p]->SetPointError(ifile, TMath::Log(0.5), skpeakerr[ifile][p]);	    
	    
	    if (sksigma[ifile][p] < 0.1 || skhv[ifile][p] <= 10 || skpeak[ifile][p]<=0){
	      ghv_sk[p]->RemovePoint(ifile);
	      cout << Form("Removing %dth point from PMT %d", ifile+1, skcable[ifile][p]) << endl << endl;
	      
	    }
        }
    }
    
    Double_t targetQ = TMath::Log(2.884); // 1.8e7 gain
    // 1.8e7/(1e-12/1.60217657e-19) = 2.884
    Double_t targetQHPK = TMath::Log(2.243); // 1.4e7 gain
    //1.4e7/(1e-12/1.60217657e-19) = 2.243
    
    TString fitOpts = "CSMRE";
    int minpoint = 4;
    int indexmin = 7;
    int indexmax = 10;   
    
    for (Int_t p = PlotRange[0]; p < min(MAXPM, PlotRange[1]); p++){
      //if (skcable[0][p] <= 0) continue;
        if (!ghv_sk[p]) {
	  
	  if (skcable[0][p]>0) outtxt_target <<  skcable[0][p] << setw(6) << " -1" << endl;
	  else if (skhv[3][p]>=0) outtxt_target << p+1 << setw(6) << " " << skhv[3][p] << endl;
	  
	  continue;
	}
	
	//fHVsk[p]->SetParLimits(1, indexmin, indexmax);
	//ghv_sk[p]->GetXaxis()->SetLimits(1500,2500);

        int npoints = ghv_sk[p]->GetN();
	//ghv_orig[p] = (TGraphErrors*)ghv_sk[p]->Clone();
        //ghv_orig[p]->SetMarkerColor(kGreen);
	//ghv_orig[p]->SetName(Form("%s_orig", ghv_sk[p]->GetName()));

        if (npoints < minpoint) {
            cout << "PMT " <<  skcable[0][p] << " has fewer than " << minpoint << " points." << endl;
	    
	    if (skcable[0][p]>0) outtxt_target <<  skcable[0][p] << setw(6) << " -1" << endl;
	    else if (skhv[3][p]>=0) outtxt_target << p+1 << setw(6) << " " << skhv[3][p] << endl;
	    
            continue;
        }

	//ghv_orig[p] = (TGraphErrors*)ghv_sk[p]->Clone();
	//ghv_orig[p]->SetMarkerColor(kGreen);
	//ghv_orig[p]->SetName(Form("%s_orig", ghv_sk[p]->GetName()));
	

        int ipmttype = -1;
	if (PMTinfo[p][0]==6) ipmttype = hkpmt;
	else if (PMTinfo[p][0]==3) ipmttype = sk2pmt;
	else if (PMTinfo[p][0]==4) ipmttype = sk3pmt;

        std::cout << endl << Form("Fitting %s cable %06d", PMTtypeNames[ipmttype].Data(), p) << std::endl;
	
	ghv_sk[p]->GetXaxis()->SetLimits(TMath::Log(skhv[3][p]-100), TMath::Log(skhv[3][p]+100));
	ghv_sk[p]->Fit(fHVsk[p], "BQN0", "", TMath::Log(skhv[3][p]-100), TMath::Log(skhv[3][p]+100));
        TVirtualFitter * gfitter = TVirtualFitter::Fitter(ghv_sk[p]);

	gfitter->SetPrecision(1);

        TFitResultPtr fitr = ghv_sk[p]->Fit(fHVsk[p], fitOpts, "", TMath::Log(skhv[3][p]-100), TMath::Log(skhv[3][p]+100));
        int status = (int)fitr;

 	Double_t chisquare = fHVsk[p]->GetChisquare()/fHVsk[p]->GetNDF();
	//int error = 4;
	//float chi2_max = 20;
	//float chi2_min = 0.001;  
	/*vector<Double_t> slopes(npoints-1);

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

	//slopes.erase( remove(slopes.begin(), slopes.end(), max_element(slopes.begin(), slopes.end()) ), slopes.end() );
	//speakpoints[2] = distance(slopes.begin(), max_element(slopes.begin(), slopes.end()));
	
	int error = 4;
	float chi2_max = 20;
	float chi2_min = 0.001;
	//int indexmin = 3;
	//int indexmax = 15;*/
	 /*if (chisquare < chi2_min){
	  fitr = ghv_sk[p]->Fit(fHVsk[p], fitOpts, "", TMath::Log(skhv[3][p]-100), TMath::Log(skhv[3][p]+100));
	  status = int(fitr);
	  chisquare = fHVsk[p]->GetChisquare()/fHVsk[p]->GetNDF();
	  }*/

        for (Int_t pointsrm = 0; pointsrm < 8; pointsrm++){
	  //if (pointsrm > ghv_sk[p]->GetN() - minpoint - 1) break;
          
	  //cout << "  status = " << status << endl;
	  
	  if ((fHVsk[p]->GetParameter(1)-7) < 0.1 || fHVsk[p]->GetParameter(1) >= indexmin || fHVsk[p]->GetParameter(1) <= indexmax) {

	    //gfitter->SetPrecision((pointsrm+1)*5);
	    //if (ipmttype == hk) ghv_sk[p]->RemovePoint(pointsrm==0?pointsrm:(nfile-pointsrm-1));
	    //ghv_sk[p]->RemovePoint(speakpoints[pointsrm]);
	    for (int ifile = 0; ifile < ghv_sk[p]->GetN(); ifile++){
	      ghv_sk[p]->SetPointError(ifile, TMath::Log(0.5), skpeakerr[ifile][p]*(1.2+pointsrm*0.1));
	    }
	    fitr = ghv_sk[p]->Fit(fHVsk[p], fitOpts, "", TMath::Log(skhv[3][p]-100), TMath::Log(skhv[3][p]+100));
	  }
	  //status = int(fitr);	  
	  //chisquare = fHVsk[p]->GetChisquare()/fHVsk[p]->GetNDF();
	  
	  if (fHVsk[p]->GetParameter(1) > indexmin && fHVsk[p]->GetParameter(1) < indexmax) break;
        }
	
	/*if ((status == error || chisquare > chi2_max || fHVsk[p]->GetParameter(1) > indexmin || fHVsk[p]->GetParameter(1) < indexmax) && ghv_sk[p]->GetN()>minpoint){
	  //gfitter->SetPrecision(30);
	  for (int extrarm = 0; extrarm < 2; extrarm++){
	    if (extrarm > ghv_sk[p]->GetN() - minpoint - 1) break;
	    gfitter->SetPrecision(30+extrarm*10);
	    if (extrarm > 0)
	      ghv_sk[p]->RemovePoint(0);
	      //ghv_sk[p]->RemovePoint(TMath::Abs(speakpoints[1]-speakpoints[0])>0?speakpoints[1]-1:0);
	    fitr = ghv_sk[p]->Fit(fHVsk[p], fitOpts, "", TMath::Log(skhv[3][p]-100), TMath::Log(skhv[3][p]+100));
	  }
	    status = (int)fitr;
	    chisquare = fHVsk[p]->GetChisquare()/fHVsk[p]->GetNDF();
	    if (status != error && chisquare <= chi2_max && fHVsk[p]->GetParameter(1) > indexmin && fHVsk[p]->GetParameter(1) < indexmax) break;
	}
        //std::cout << Form("Fitting SK%d cable %06d", PMTinfo[p][0] - 1, p) << std::endl;*/
        
	//ghv_orig[p]->Draw("AP");
	ghv_sk[p]->Draw("sameP");
        //c1->Update();
	fout->cd();
        ghv_sk[p]->Write();
        //ghv_orig[p]->Write();
	
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

	geighthv[ipmttype] = TMath::Exp(fHVinvsk[p]->Eval(targetQ));
	gfourhv[ipmttype] = TMath::Exp(fHVinvsk[p]->Eval(targetQHPK));
	if (isinf(geighthv[ipmttype])) geighthv[ipmttype] = -3;
	else if (isnan(geighthv[ipmttype])) geighthv[ipmttype] = -4;

	/*if (beta[ipmttype] >= indexmax || beta[ipmttype] <= indexmin) {
	  //norm[ipmttype] = 0;
	  //normerr[ipmttype] = 0;
	  //beta[ipmttype] = 0;
	  //beta[ipmttype] = 0;
	  gfourhv[ipmttype] = skhv[1][p];
	  geighthv[ipmttype] = skhv[3][p];
	  rchi2[ipmttype] = 100;
	  prob[ipmttype] = 0;
	}*/

	//else if (geighthv[ipmttype] >= 2500 || gfourhv[ipmttype] <= 100){
	//gfourhv[ipmttype] = skhv[1][p];
	//geighthv[ipmttype] = skhv[3][p];
	//}

	//gfourhv[ipmttype] = TMath::Exp(fHVinvsk[p]->Eval(targetQHPK));
	geighthverr[ipmttype] = TMath::Exp(fHVinvskerr[p]->Eval(targetQ)) - geighthv[ipmttype];
	gfourhverr[ipmttype] = TMath::Exp(fHVinvskerr[p]->Eval(targetQHPK)) - gfourhv[ipmttype];
	gshifthv[ipmttype] = geighthv[ipmttype] - SKPMThv[p];
        resolution[ipmttype] = sksigma[2][p]/skpeak[2][p];
	tr[ipmttype]->Fill();

	outtxt_target <<  skcable[0][p] << setw(6) << " " << geighthv[ipmttype] << endl;
           
	if (beta[ipmttype] <= indexmin){

	  fHVsk[p]->SetLineColor(2);
	  outtxt_badfit  << setw(10) << skcable[0][p] << setw(4) << " "+PMTtypeNames[ipmttype] << setw(8) << Form("%1.2f", gfourhv[ipmttype]) << setw(8) << Form("%1.2f", geighthv[ipmttype]) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(15) << "     Flat_Curve" << "\n";
	}

	else if (beta[ipmttype] >= indexmax){
	  fHVsk[p]->SetLineColor(46);
	  outtxt_badfit << setw(10) << skcable[0][p] << setw(4) << " "+PMTtypeNames[ipmttype] << setw(8) << Form("%1.2f", gfourhv[ipmttype]) << setw(8) << Form("%1.2f", geighthv[ipmttype]) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(15) << "    Steep_Curve" << "\n";
	  
	}

	/*else if (rchi2[ipmttype] > chi2_max){
	    fHVsk[p]->SetLineColor(6);
	    if (std::find(KEKcable.begin(), KEKcable.end(), skcable[0][p]) != KEKcable.end())
	      outtxt_badfit << setw(6) <<"   KEK" << setw(10) << skcable[0][p] << setw(4) << " "+PMTtypeNames[ipmttype] << setw(8) << Form("%1.2f", gfourhv[ipmttype]) << setw(8) << Form("%1.2f", geighthv[ipmttype]) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                Chi2_>>_1" << "\n";
	    
	    else
	      outtxt_badfit << setw(6) <<"    no" << setw(10) << skcable[0][p] << setw(4) << " "+PMTtypeNames[ipmttype] << setw(8) << Form("%1.2f", gfourhv[ipmttype]) << setw(8) << Form("%1.2f", geighthv[ipmttype]) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                Chi2_>>_1" << "\n";
	}
	if (rchi2[ipmttype] == 0 && prob[ipmttype] == 1){
	  fHVsk[p]->SetLineColor(3);
	  
	  if (std::find(KEKcable.begin(), KEKcable.end(), skcable[0][p]) != KEKcable.end()) outtxt_badfit << setw(6) <<"   KEK" outtxt_badokfit << setw(10) << skcable[0][p] << setw(4) << " " << PMTtypeNames[ipmttype] << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << "\n";
	  else outtxt_badfit << setw(6) <<"    no" << setw(10) << skcable[0][p] << setw(4) << " " << PMTtypeNames[ipmttype] << setw(8) << Form("%1.2f", gfourhv[ipmttype]) << setw(8) << Form("%1.2f", geighthv[ipmttype]) << setw(25) << Form("%1.2e(%1.2e)", norm[ipmttype], normerr[ipmttype]) << setw(20) << Form("%1.2f(%1.2f)", beta[ipmttype], betaerr[ipmttype]) << setw(15) << Form("%1.3f", rchi2[ipmttype] ) << setw(15) << Form("%1.2e", fHVsk[p]->GetProb()) << setw(25) << "                 Chi2_=_0" << "\n";
	  }*/
       
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
    
    for (Int_t p = PlotRange[0]; p < min(MAXPM, PlotRange[1]); p++){
      
        if (!ghv_sk[p]) continue;
        if (ghv_sk[p]->GetN() < minpoint) continue;

        //ghv_orig[p]->Draw("AP");
        ghv_sk[p]->Draw("AP");
	
        if (p%1000==0) cout << "Drawing cable " << skcable[0][p] << endl;

        t->Clear();

        if (TMath::Exp(fHVinvsk[p]->Eval(targetQ)) < 2500 && TMath::Exp(fHVinvskerr[p]->Eval(targetQ) < 2500)){

	t->AddText(Form("HV at 1.4e7 gain : %4.2f [V]", TMath::Exp(fHVinvsk[p]->Eval(targetQHPK))));
          //t->AddText(Form("HV at 1.4e7 gain : %4.2f#pm%4.2f [V]", fHVinvsk[p]->Eval(targetQHPK), fHVinvskerr[p]->Eval(targetQHPK) - fHVinvsk[p]->Eval(targetQHPK)));
	  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlue);
	  t->AddText(Form("HV at 1.8e7 gain : %4.2f [V]", TMath::Exp(fHVinvsk[p]->Eval(targetQ))));      
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
	}

	else cout << "Target HV larger than 2500V for cable: " << skcable[0][p+1] << endl;
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
