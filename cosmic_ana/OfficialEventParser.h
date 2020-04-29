#ifndef OfficialEventParser_h
#define OfficialEventParser_h

#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "Event.h"
#include <iostream>
#include <fstream>
#include <string>
#include <TH3.h>
#include <TVector3.h>
#include "global.h"
//#include "pdg_codes.h"
#include "geotank.h"

using namespace std;
using namespace global;

namespace fiTQun_parameters{
  enum particleIndices {ig, ie, imu, ipip, ikp, ipr, icg};
  enum TrackIndices {itrkg, itrke, itrkmu, itrkpip, itrkkp, itrkp, itrkcg, itrkpi0, itrkpipknk};

  const int nPartHyps = 2; // imu, ipip (particleIndices+2)
  
  const int fq_pi0_mode = 0;
}
using namespace fiTQun_parameters;

// Peakfinding method (Use prefit vertex)
const int imethod=0;

enum selectReco {bothreco, allapfit, allfitqun};
const int recoSelect = bothreco;

  
const int nSE = 2;
const int nSample = 2;
const int nCuts = 1;
const int nNrgs = 2;
const int nMRhyps = 3;

const int nGates = 3;

const int ise_muon = 0;    
const int ise_dcye = 1;


const int nMomBins[nSE] = {100, 50};
const double momLow[nSE] = {0 ,0};
const double momHigh[nSE] = {10000, 100};


const int nRanges = 7;
const double rangeLow = 0;
const double rangeHigh = 35; // cm
const double rangeStep = 5;  // (rangeHigh-rangeLow)/nRanges



const int yearRange[2] = {2008, 2014};
const int nMonthsPerYear = 12;
const int nYearBins = nMonthsPerYear*(yearRange[1]-yearRange[0]);
const TString monthNames[nMonthsPerYear] = {"Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

const int runRange[2] = {61500, 70900};
const int nRunBins = (runRange[1]-runRange[0]);
      
const int nComps = 2;

const int nLikeBins = 100;

// Muon momentum binning for likelihoods
const int nMuonMomBinsLike = 10;
const double muonMomBinsLike[nMuonMomBinsLike+1] = {0, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 12000};
// Muon mom vs potot: 
/*
Minimizer is Linear
Chi2                      =       100492
NDf                       =           11
p0                        =      166.457   +/-   0.202932    
p1                        =    0.0927583   +/-   6.66376e-06 
*/

const int nMuonLikes = 5;
const TString muonLikeString[nMuonLikes] = {"Le_Lmu","Le_Lpi0","Lmu_Lpi0","L2r_L1r","L2r_L1r_cut"};
//const TString muonLikeTitle[nMuonLikes] = {"ln(L_{e}/L_{#mu})","ln(L_{e}/L_{#pi^{0}})","ln(L_{#mu}/L_{#pi^{0}})"};
const TString muonLikeTitle[nMuonLikes] = {"Distance from Cut Line in ln(L_{e}/L_{#mu}) vs p_{e} Space","ln(L_{e}/L_{#pi^{0}})","ln(L_{#mu}/L_{#pi^{0}})","ln(L_{2R}/L_{1R})","Distance from 2R/1R Cut Line"};
const double muonLikeTopLeft[nMuonLikes] = {3000, 600, 5000, 2000, 2000};
const double muonLikeTopRight[nMuonLikes] = {7000, 24000, 30000, 25000, 4000};
const double muonLikeBotLeft[nMuonLikes] = {-3000, -800, -5000, -2000, -2000};
const double muonLikeBotRight[nMuonLikes] = {-45000, -55000, -45000, -23000, -35000};

// Electron momentum binning for likelihoods
const int nDcyeMomBinsLike = 8;
const double dcyeMomBinsLike[nDcyeMomBinsLike+1] = {0, 15, 25, 30, 35, 40, 45, 50, 100};

const int nDcyeLikes = 1;
const TString dcyeLikeString[nDcyeLikes] = {"Le_Lmu"};//,"L2r_L1r"};
const TString dcyeLikeTitle[nDcyeLikes] = {"Distance from Cut Line in ln(L_{e}/L_{#mu}) vs p_{e} Space"};//,"ln(L_{2R}/L_{1R})"};
const double dcyeLikeTopLeft[nDcyeLikes] = {400};//, 500};
const double dcyeLikeTopRight[nDcyeLikes] = {600};//, 1500};
const double dcyeLikeBotLeft[nDcyeLikes] = {-500};//, -500};
const double dcyeLikeBotRight[nDcyeLikes] = {-1500};//, -600};

const int nRangeCuts = 4;
const TString rangeTitle[nRangeCuts] = {
  "Official E-Scale",
  "Top & Side Entering",
  "Top Entering",
  "Side Entering"
};

class OfficialEventParser{
private:
        Event *E;

	TFile * efile[nSample][nSE];

	vector<TH1 *> h_vect_hists[nSample][nRecos];

	TH2D * h_cosmic_emom_vs_year[nSample][nRecos];
	TH2D * h_cosmic_mumom_vs_year[nSample][nRecos];
	TH2D * h_cosmic_lifetime_vs_year[nSample][nRecos];
	TH3D * h_cosmic_emom_vs_lifetime_vs_year[nSample][nRecos];
	TH3D * h_cosmic_range_vs_mumom_vs_year[nSample][nRecos];

	TH2D * h_cosmic_nmue_vs_year[nSample][nRecos];
	TH2D * h_cosmic_ingate_vs_year[nSample][nRecos];

	TH2D * h_cosmic_wall_vs_year[nSample][nSE][nRecos];
	TH2D * h_cosmic_z_vs_year[nSample][nSE][nRecos][nCuts];
	TH2D * h_cosmic_r2_vs_year[nSample][nSE][nRecos][nCuts];
	TH2D * h_cosmic_zzoom_vs_year[nSample][nSE][nRecos][nCuts];
	TH2D * h_cosmic_r2zoom_vs_year[nSample][nSE][nRecos][nCuts];

	TH1D * h_cosmic_z_nearwall[nSample][nSE][nRecos][nCuts];
	TH1D * h_cosmic_r2_nearwall[nSample][nSE][nRecos][nCuts];

	TH2D * h_cosmic_pcflg_vs_year[nSample][nSE][nRecos];
	
	TH2D * h_cosmic_vtx_prfvtx_diff_vs_year[nSample][nSE][nRecos];

	TH2D * h_cosmic_mom_over_range_rangesep_vs_year[nSample][nRecos][nRanges][nRangeCuts];

	TH1D * h_cosmic_mom_binning[nSample][nRecos][nSE];
	TH2D * h_cosmic_Lmuon_vs_year[nSample][nRecos][nMuonLikes][nMuonMomBinsLike];
	TH2D * h_cosmic_Ldcye_vs_year[nSample][nRecos][nDcyeLikes][nDcyeMomBinsLike];

	TH2D * h_cosmic_muonmisid_vs_year[nSample][nRecos][nMuonMomBinsLike];
	TH2D * h_cosmic_dcyemisid_vs_year[nSample][nRecos][nDcyeMomBinsLike];

	TH2D * h_cosmic_Lmuon_vs_zenith[nSample][nRecos][nMuonLikes][nMuonMomBinsLike];
	TH2D * h_cosmic_Ldcye_vs_zenith[nSample][nRecos][nDcyeLikes][nDcyeMomBinsLike];

	TH2D * h_cosmic_lnLrat1R_vs_mom[nSample][nSE][nRecos];
	TH2D * h_cosmic_lnLrat2R1R_vs_mom[nSample][nSE][nRecos];

	TH2D * h_cosmic_n50_vs_year[nSample][nRecos];
	TH2D * h_cosmic_q50_vs_year[nSample][nRecos];
	TH2D * h_cosmic_clnhits_vs_year[nSample][nSE][nRecos];
	TH2D * h_cosmic_cltotq_vs_year[nSample][nSE][nRecos];

	TH2D * h_cosmic_cltstart_vs_year[nSample][nSE][nRecos];
	TH2D * h_cosmic_cltend_vs_year[nSample][nSE][nRecos];
	TH2D * h_cosmic_cllength_vs_year[nSample][nSE][nRecos];
	TH2D * h_cosmic_cllength_vs_cltstart[nSample][nSE][nRecos];

	TH2D * h_cosmic_cltotqnorm_vs_year[nSample][nSE][nRecos];
	TH2D * h_cosmic_clnhitsnorm_vs_year[nSample][nSE][nRecos];

	TH2D * h_cosmic_fqtotqnorm_vs_year[nSample][nSE][nRecos];
	TH2D * h_cosmic_fqnhitsnorm_vs_year[nSample][nSE][nRecos];

	TH2D * h_cosmic_fqnhits_vs_year[nSample][nSE][nRecos];
	TH2D * h_cosmic_fqtotq_vs_year[nSample][nSE][nRecos];
	TH2D * h_cosmic_fqtotmu_vs_year[nSample][nSE][nRecos];
	TH2D * h_cosmic_fq1rtotmu_vs_year[nSample][nSE][nRecos];

	TH1D * h_cosmic_mom_res[nSample][nSE][nRecos][nCuts];
	TH2D * h_cosmic_mom_res_vs_trumom[nSample][nSE][nRecos][nCuts];
	TH1D * h_cosmic_vtx_res[nSample][nSE][nRecos][nCuts];
	TH1D * h_cosmic_tvtx_res[nSample][nSE][nRecos][nCuts];
	TH1D * h_cosmic_lvtx_res[nSample][nSE][nRecos][nCuts];
	TH1D * h_cosmic_range_res[nSample][nRecos][nCuts];
	TH1D * h_cosmic_dir_res[nSample][nSE][nRecos][nCuts]; 
	TH1D * h_cosmic_time_res[nSample][nSE][nRecos][nCuts];
	TH1D * h_cosmic_ddir_res[nSample][nRecos][nCuts];

	bool fcmc;      // true : MC, false : data
	int  _version;  // SK detector version
	int  _nrgSepType;  // nrgSepType (0: all, 1: sub-GeV, 2: multi-GeV)
	int  _maxE;  
	int  _fileNum;  
	bool _osc;      // flag to oscillate MC or not
	int _datatype; 
	int _towall_dir; // towall_dir for this event is filled.
	bool _from_wall_flag; // if this flag is true, calculate from_wall distance instead of to_wall.
	double _datalivetime; // livetime of data
	double nentries_norm[nSample][nRecos];

	enum towall_dir {out_of_tank=0, to_side, to_top, to_bottom};
	enum find_true_particle { find_e, find_mu};

	// Computed variables
	

	// Truth variables
	int found_true_particle[nSE];
	int pid_tru[nSE];
     	float mom_tru[nSE], etimev;
	TVector3 v_dir_tru[nSE], v_vtx_tru[nSE];
	
	// Reconstructed variables
	int fillPID;
	float lnLrat1R2R_muon;
	float etime;

	float lnLrat1R[nSE];
	float mom_rec[nSE];
	int imom_rec_idx[nSE];
	TVector3 v_dir_rec[nSE], v_vtx_rec[nSE];
	float wall[nSE], r2[nSE], z[nSE], abs_vtx_prfvtx_diff[nSE];

	double likelihood_mom[nSE];

	double MRnll[nMRhyps+1];
	int MrPID[nMRhyps+1][nMRhyps+1];

	// Resolutions
	float mom_res[nSE], vtx_res[nSE], tvtx_res[nSE], lvtx_res[nSE], dir_res[nSE], time_res, ddir_res, range_res;
	
	// Histo filling
	void FillResolutionHistos(int imc, int ireco, int ise, int icut);
	void FillLikelihoodHistos(int imc, int ireco, int ise, int icut);
	
	TH1D* MakeNewHisto(const char* name, const char* title, Int_t nbinsx, const Double_t *xbins, vector<TH1 *>&h_vect);
	TH1D* MakeNewHisto(const char* name, const char* title, Int_t nbinsx, const Double_t xlow, const Double_t xhigh, vector<TH1 *>&h_vect);

	TH2D* MakeNewHisto(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, vector<TH1 *>&h_vect);
	
	TH3D* MakeNewHisto(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, Int_t nbinsz, Double_t zlow, Double_t zup, vector<TH1 *>&h_vect);
	
	TH3D* MakeNewHisto(const char* name, const char* title, Int_t nbinsx, const Double_t* xbins, Int_t nbinsy, const Double_t* ybins, Int_t nbinsz, const Double_t* zbins, vector<TH1 *>&h_vect);

public :

	OfficialEventParser();
	void InitHistograms(int ireco, int imc);
	void InitEvent();
	void Parse(int ireco=0);
	void Print();
	void Terminate();
	void operator()(  TChain *fChain, int ireco );
	
	void FillCosmic(int ireco);

	double weight();

	bool isCosmic(int ireco);
	bool is1R(int ireco);
	int   isInGate(int ireco, int ise);
	bool isSubGeV(int ireco);
	bool RvtxInFid(float *vertex);
	bool ZvtxInFid(float *vertex);
	bool hasRvtxInFid(int ireco);
	bool hasZvtxInFid(int ireco);

        void set_version( int x)  { _version = x; }
	void set_nrgSepType( int x)  { _nrgSepType = x; }
	void set_maxE( int x)  { _maxE = x; }
	void set_fileNum( int x)  { _fileNum = x; }
        int  get_version( void)   { return _version;}
	bool get_osc(void)        { return _osc;}
	void set_osc(bool x)      { _osc = x;}
	int  get_datatype(void)   { return _datatype;}
	void set_datatype(int x)  { _datatype = x;}

	float cut_nhitac(); // return cut value of nhitac for sk1,sk2, and sk3
	float cut_potot();
	int cut_ehit( int i );
	
	int get_truth_info(int imc, int ireco, int ise, int &pid_tru, float &mom_tru, TVector3 &v_dir_tru, TVector3 &v_vtx_tru, float &etimev);
	
	double get_towall(float *vertex, float *dir);
	double get_fromwall(float *vertex, float *dir);

	double getLnLrat1R(int ireco, int ise);
	double getLnLrat1R2R(int ireco, int &fillPID);
	int get1rType(int ireco,int ise);

	double getMin1Rnll(int ireco, int &PID1);
	double getMin2Rnll(int ireco, int &PID1, int &PID2);
	double getMin3Rnll(int ireco, int &PID1, int &PID2, int &PID3);
	double getMin4Rnll(int ireco, int &PID1, int &PID2, int &PID3, int &PID4);
	double getRingCountingLikelihood(int ireco);

	
	void  set_towall_dir( int x){ _towall_dir = x;}
	int   get_towall_dir( void ){ return _towall_dir ;}
	double calc_towall(double x, double y, double z, float *vertex);
	
	void  set_from_wall_flag ( bool x ) { _from_wall_flag = x;}
	bool  get_from_wall_flag ( void ) { return _from_wall_flag;}
	
	void  set_datalivetime ( double x ) { _datalivetime = x; cout << _datalivetime << endl;}
	double get_datalivetime( void ) { return _datalivetime;} 

	int get_num_ring(int ireco);
	
	void *get_eststop_vertex(int ireco, int ise, float *stop_vertex);
	float *get_vertex(int ireco, int ise=0);
	float *get_1Rdir(int ireco, int ise=0);
	float get_1Rmom(int ireco, int ise=0);
	float get_dcye_time(int ireco, int ise=0);

	void get_2Rdir(int ireco, float *dir1, float *dir2);

	float get_nhits_over_time(int ireco);
	int get_n50(int ireco, int ise);
	int get_q50(int ireco, int ise);

	float get_evis(int ireco);
	float get_potot(int ireco);
	float get_wall(float *vertex);

	int   muedcy(int ireco);

	double get_year(int ireco);
	double get_run(int ireco);

	void rotateVector(float *vector, double theta, double phi);


	double OfficialEventParser::minimum_distance(TVector2 v, TVector2 w, TVector2 p) {
	  // Return minimum distance between line segment vw and point p
	  const double l2 = (v-w).Mod2();  // i.e. |w-v|^2 -  avoid a sqrt
	  //if (l2 == 0.0) return distance(p, v);   // v == w case
	  // Consider the line extending the segment, parameterized as v + t (w - v).
	  // We find projection of point p onto the line. 
	  // It falls where t = [(p-v) . (w-v)] / |w-v|^2
	  const double t = (p - v)*(w - v) / l2;
	  //if (t < 0.0) return distance(p, v);       // Beyond the 'v' end of the segment
	  //else if (t > 1.0) return distance(p, w);  // Beyond the 'w' end of the segment
	  const TVector2 projection = v + t * (w - v);  // Projection falls on the segment
	  
	  double distance = (p-projection).Mod();
	  return distance;
	};


};
#endif

#ifdef __OFFICIALEVENTPARSER_CXX
void OfficialEventParser::operator()(  TChain *fChain, int ireco ){

  fChain->SetBranchStatus("*",0);
  
  if ((ireco==fitqun&&recoSelect==bothreco) || recoSelect==allfitqun) {
  
    //cout << "Setting fiTQun branch addresses..." << endl;

    fChain->SetBranchAddress("nring", &E->nring);
    fChain->SetBranchAddress("nrun", &E->nrun);

    
    //fChain->SetBranchAddress("nev", &E->nev);
    //fChain->SetBranchAddress("nsub", &E->nsub);
    //fChain->SetBranchAddress("cate", &E->cate);
    
    fChain->SetBranchAddress("potot", &E->potot);
    fChain->SetBranchAddress("nhit", &E->nhit);
    //fChain->SetBranchAddress("pomax", &E->pomax);
    //fChain->SetBranchAddress("potota", &E->potota);
    //fChain->SetBranchAddress("nhita", &E->nhita);
    fChain->SetBranchAddress("nhitac", &E->nhitac);
    //fChain->SetBranchAddress("pomaxa", &E->pomaxa);
    fChain->SetBranchAddress("wall", &E->wall);
    fChain->SetBranchAddress("evis", &E->evis);
    /*
      unused in this OfficialEventParser
      fChain->SetBranchAddress("rtsum", &E->rtsum);
      fChain->SetBranchAddress("rtmax", &E->rtmax);
      fChain->SetBranchAddress("wlen", &E->wlen);
    */

    //fChain->SetBranchAddress("ip", E->ip);
    //fChain->SetBranchAddress("pos", E->pos);
    //fChain->SetBranchAddress("dir", E->dir);
    //fChain->SetBranchAddress("dirtot", E->dirtot);
    //fChain->SetBranchAddress("ang", E->ang);
    //fChain->SetBranchAddress("rtot", E->rtot);
    //fChain->SetBranchAddress("amom", E->amom);
    //fChain->SetBranchAddress("rtote", E->rtote);
    //fChain->SetBranchAddress("amome", E->amome);
    //fChain->SetBranchAddress("rtotm", E->rtotm);
    //fChain->SetBranchAddress("amomm", E->amomm);
    //fChain->SetBranchAddress("nsube", &E->nsube);
    //fChain->SetBranchAddress("ndcy", &E->ndcy);
    //fChain->SetBranchAddress("ngate", &E->ngate);
    //fChain->SetBranchAddress("nbye", &E->nbye);
    //fChain->SetBranchAddress("probms", E->probms);
    //fChain->SetBranchAddress("prmslg", E->prmslg);
      fChain->SetBranchAddress("date", E->date);
    /*
      unused in this OfficialEventParser
      fChain->SetBranchAddress("time", E->time);
      fChain->SetBranchAddress("elpsday", &E->elpsday);
      fChain->SetBranchAddress("numpo", E->numpo);
      fChain->SetBranchAddress("apos", E->apos);
      fChain->SetBranchAddress("adir", E->adir);
      fChain->SetBranchAddress("aang", &E->aang);
      fChain->SetBranchAddress("agood", &E->agood);
      fChain->SetBranchAddress("wgain", &E->wgain);
      fChain->SetBranchAddress("nbad", &E->nbad);
      fChain->SetBranchAddress("nbada", &E->nbada);
    */
    //fChain->SetBranchAddress("msdir", E->msdir);
    /*
      fChain->SetBranchAddress("amomp", E->amomp);
      fChain->SetBranchAddress("ange", E->ange);
      fChain->SetBranchAddress("angm", E->angm);
      fChain->SetBranchAddress("angp", E->angp);
      fChain->SetBranchAddress("ntot", E->ntot);
      fChain->SetBranchAddress("probth", E->probth);
      fChain->SetBranchAddress("probpt", E->probpt);
    */
    //fChain->SetBranchAddress("pi0like", E->pi0like);
    //fChain->SetBranchAddress("pi0_e", E->pi0_e);
    //fChain->SetBranchAddress("pi0_dir", E->pi0_dir);
    //fChain->SetBranchAddress("pi0mass", E->pi0mass);
    /*
      fChain->SetBranchAddress("evisold", E->evisold);
      fChain->SetBranchAddress("evisoldxe", E->evisoldxe);
      fChain->SetBranchAddress("evisnew", E->evisnew);
    */
    //fChain->SetBranchAddress("nmue", &E->nmue);
    //fChain->SetBranchAddress("etype", E->etype);
    //fChain->SetBranchAddress("etime", E->etime);
    //fChain->SetBranchAddress("epos", E->epos);
    //fChain->SetBranchAddress("edir", E->edir);
    //fChain->SetBranchAddress("egood", E->egood);
    //fChain->SetBranchAddress("ehit", E->ehit);
    //fChain->SetBranchAddress("mueprob", E->mueprob);

    /*
      unused in this OfficialEventParser

      fChain->SetBranchAddress("Rnring", &E->Rnring);
      fChain->SetBranchAddress("Rdir", E->Rdir);
      fChain->SetBranchAddress("Rang", E->Rang);
      fChain->SetBranchAddress("Riring", &E->Riring);
      fChain->SetBranchAddress("Rtwout", E->Rtwout);
      fChain->SetBranchAddress("Rtwith", E->Rtwith);
      fChain->SetBranchAddress("Alwout", &E->Alwout);
      fChain->SetBranchAddress("Alwith", &E->Alwith);
      fChain->SetBranchAddress("Qsmi", &E->Qsmi);
      fChain->SetBranchAddress("Qsmo", &E->Qsmo);
      fChain->SetBranchAddress("Qexi", &E->Qexi);
      fChain->SetBranchAddress("Qexo", &E->Qexo);
      fChain->SetBranchAddress("Pe5d", &E->Pe5d);
      fChain->SetBranchAddress("En5d", &E->En5d);
      fChain->SetBranchAddress("Eh5d", &E->Eh5d);
      fChain->SetBranchAddress("Pe5do", &E->Pe5do);
      fChain->SetBranchAddress("En5do", &E->En5do);
      fChain->SetBranchAddress("Eh5do", &E->Eh5do);
      fChain->SetBranchAddress("Rtadd", &E->Rtadd);
      fChain->SetBranchAddress("Pdgeta", &E->Pdgeta);
      fChain->SetBranchAddress("Pd5d", &E->Pd5d);
      fChain->SetBranchAddress("Pdthre", &E->Pdthre);
      fChain->SetBranchAddress("Pd5do", &E->Pd5do);
      fChain->SetBranchAddress("Delpd", &E->Delpd);
      fChain->SetBranchAddress("Ropena", E->Ropena);
      fChain->SetBranchAddress("Maxth", &E->Maxth);
      fChain->SetBranchAddress("Pkang", &E->Pkang);
      fChain->SetBranchAddress("Qrfct", &E->Qrfct);
      fChain->SetBranchAddress("Pdfct", &E->Pdfct);
      fChain->SetBranchAddress("Pkfct", &E->Pkfct);
      fChain->SetBranchAddress("Agfct", &E->Agfct);
    */


    //fChain->SetBranchAddress("Dlfct", &E->Dlfct);

    /*
      unused in this OfficialEventParser

      fChain->SetBranchAddress("Iflag", &E->Iflag);
      fChain->SetBranchAddress("Pmfct", &E->Pmfct);
      fChain->SetBranchAddress("Imfct", &E->Imfct);
      fChain->SetBranchAddress("Rilike", &E->Rilike);
      fChain->SetBranchAddress("ri_ver", &E->ri_ver);
      fChain->SetBranchAddress("ri_pid", &E->ri_pid);
      fChain->SetBranchAddress("ri_nring", &E->ri_nring);
      fChain->SetBranchAddress("ri_flag", E->ri_flag);
      fChain->SetBranchAddress("ri_dlfct", E->ri_dlfct);
      fChain->SetBranchAddress("ri_pdfct", E->ri_pdfct);
      fChain->SetBranchAddress("ri_pkfct", E->ri_pkfct);
      fChain->SetBranchAddress("ri_vafct", E->ri_vafct);
      fChain->SetBranchAddress("ri_total", E->ri_total);
      fChain->SetBranchAddress("ri_dir", E->ri_dir);
      fChain->SetBranchAddress("ri_imfct", E->ri_imfct);
      fChain->SetBranchAddress("ri_pmfct", E->ri_pmfct);
    */

    fChain->SetBranchAddress("npar", &E->npar);
    fChain->SetBranchAddress("wallv", &E->wallv);
    fChain->SetBranchAddress("ipv", E->ipv);
    fChain->SetBranchAddress("posv", E->posv);
    fChain->SetBranchAddress("dirv", E->dirv);
    fChain->SetBranchAddress("pmomv", E->pmomv);
    //   fChain->SetBranchAddress("light_flag", &E->light_flag);
    //fChain->SetBranchAddress("npar2", &E->npar2);
    //fChain->SetBranchAddress("wallv2", E->wallv2);
    //fChain->SetBranchAddress("ipv2", E->ipv2);
    //fChain->SetBranchAddress("iorg", E->iorg);
    //fChain->SetBranchAddress("posv2", E->posv2);
    //fChain->SetBranchAddress("dirv2", E->dirv2);
    //fChain->SetBranchAddress("pmomv2", E->pmomv2);
    //fChain->SetBranchAddress("numnu", &E->numnu);
    //fChain->SetBranchAddress("mode", &E->mode);
    //fChain->SetBranchAddress("ipnu", E->ipnu);
    //fChain->SetBranchAddress("pnu", E->pnu);
    //fChain->SetBranchAddress("dirnu", E->dirnu);
    /*
      fChain->SetBranchAddress("flxg", E->flxg);
      fChain->SetBranchAddress("flxh01", E->flxh01);
      fChain->SetBranchAddress("kflux", E->kflux);
      fChain->SetBranchAddress("bs71", E->bs71);
      fChain->SetBranchAddress("bs74", E->bs74);
      fChain->SetBranchAddress("flxf", E->flxf);
      fChain->SetBranchAddress("flxh1d", E->flxh1d);
      fChain->SetBranchAddress("flxb03", E->flxb03);
      fChain->SetBranchAddress("flxf03", E->flxf03);
    */
    //fChain->SetBranchAddress("flxh06", E->flxh06);
    //fChain->SetBranchAddress("live", &E->live);
    /*
      unused in this OfficialEventParser
      fChain->SetBranchAddress("sacth", &E->sacth);
      fChain->SetBranchAddress("sactg", &E->sactg);
      fChain->SetBranchAddress("sacth1d", &E->sacth1d);
      fChain->SetBranchAddress("scan", E->scan);
      fChain->SetBranchAddress("dirtotepi", E->dirtotepi);
      fChain->SetBranchAddress("dirtotenpi", E->dirtotenpi);
      fChain->SetBranchAddress("dirtotmue", E->dirtotmue);
      fChain->SetBranchAddress("dirsum", E->dirsum);
      fChain->SetBranchAddress("etot", &E->etot);
      fChain->SetBranchAddress("etotepi", &E->etotepi);
      fChain->SetBranchAddress("etotenpi", &E->etotenpi);
      fChain->SetBranchAddress("etotmue", &E->etotmue);
      fChain->SetBranchAddress("oscweight", E->oscweight);
    */
    //fChain->SetBranchAddress("oscwgt", &E->oscwgt);

    /*
      unused in this OfficialEventParser
      fChain->SetBranchAddress("ent_pos", E->ent_pos);
      fChain->SetBranchAddress("ent_dir", E->ent_dir);
      fChain->SetBranchAddress("length", &E->length);
      fChain->SetBranchAddress("tr_mom1", &E->tr_mom1);
      fChain->SetBranchAddress("A_ent_mom", &E->A_ent_mom);
      fChain->SetBranchAddress("A_ent_pos", E->A_ent_pos);
      fChain->SetBranchAddress("A_ent_dir", E->A_ent_dir);
      fChain->SetBranchAddress("A_ext_mom", &E->A_ext_mom);
      fChain->SetBranchAddress("A_ext_pos", E->A_ext_pos);
      fChain->SetBranchAddress("A_ext_dir", E->A_ext_dir);
      fChain->SetBranchAddress("Fit_pos", E->Fit_pos);
      fChain->SetBranchAddress("Fit_dir", E->Fit_dir);
      fChain->SetBranchAddress("Fit_len", &E->Fit_len);
      fChain->SetBranchAddress("Fit_mom", &E->Fit_mom);
      fChain->SetBranchAddress("Fit_pid", &E->Fit_pid);
      fChain->SetBranchAddress("Um_ehit8m", &E->Um_ehit8m);
      fChain->SetBranchAddress("Um_ohit8m", &E->Um_ohit8m);
      fChain->SetBranchAddress("Um_qent", &E->Um_qent);
      fChain->SetBranchAddress("Sh_chi1p", &E->Sh_chi1p);
      fChain->SetBranchAddress("Sh_delta", &E->Sh_delta);
      fChain->SetBranchAddress("Sh_mean", &E->Sh_mean);
      fChain->SetBranchAddress("Sh_meanq", &E->Sh_meanq);
      fChain->SetBranchAddress("Tr_stop", E->Tr_stop);
      fChain->SetBranchAddress("Tr_mom", &E->Tr_mom);
      fChain->SetBranchAddress("Tr_len", &E->Tr_len);
      fChain->SetBranchAddress("Tr_len1", &E->Tr_len1);
      fChain->SetBranchAddress("Pid_flg", &E->Pid_flg);
      fChain->SetBranchAddress("Crs1", &E->Crs1);
      fChain->SetBranchAddress("Crs2", &E->Crs2);
      fChain->SetBranchAddress("iclass", &E->iclass);
      fChain->SetBranchAddress("mu_class", &E->mu_class);
      fChain->SetBranchAddress("mu_dec", &E->mu_dec);
      fChain->SetBranchAddress("mu_dir", E->mu_dir);
      fChain->SetBranchAddress("mu_pos", E->mu_pos);
      fChain->SetBranchAddress("mu_good", &mu_good);
      fChain->SetBranchAddress("history", &history);
      fChain->SetBranchAddress("Pdst", &Pdst, &b_Pdst);
      fChain->SetBranchAddress("idoff", &idoff, &b_idoff);
      fChain->SetBranchAddress("anthit", &anthit, &b_anthit);
      fChain->SetBranchAddress("idseq", &idseq, &b_idseq);
      fChain->SetBranchAddress("tstfrac", &tstfrac, &b_tstfrac);
      fChain->SetBranchAddress("judge", &judge, &b_judge);
      fChain->SetBranchAddress("Upcrs1", &Upcrs1, &b_Upcrs1);
      fChain->SetBranchAddress("Upcrs2", &Upcrs2, &b_Upcrs2);
      fChain->SetBranchAddress("lst", &E->lst);
      fChain->SetBranchAddress("jd", &E->jd);
      fChain->SetBranchAddress("fjd", &E->fjd);
      fChain->SetBranchAddress("alt", &E->alt);
      fChain->SetBranchAddress("azi", &E->azi);
      fChain->SetBranchAddress("ra", &E->ra);
      fChain->SetBranchAddress("dec", &E->dec);
      fChain->SetBranchAddress("glat", &E->glat);
      fChain->SetBranchAddress("glong", &E->glong);
      fChain->SetBranchAddress("nuceff_version", &E->nuceff_version);
      fChain->SetBranchAddress("charge_exchange", &E->charge_exchange);
      fChain->SetBranchAddress("absorbed", &E->absorbed);
      fChain->SetBranchAddress("multipi_gen", &E->multipi_gen);
      fChain->SetBranchAddress("scattering", &E->scattering);
      fChain->SetBranchAddress("piless_dcy", &E->piless_dcy);
    */


   //fChain->SetBranchAddress("cluster_ncand", &E->cluster_ncand);
   fChain->SetBranchAddress("cluster_tstart", E->cluster_tstart);
   fChain->SetBranchAddress("cluster_tend", E->cluster_tend);
   fChain->SetBranchAddress("cluster_nhits", E->cluster_nhits);
   fChain->SetBranchAddress("cluster_totq", E->cluster_totq);
   //fChain->SetBranchAddress("cluster_goodflag", E->cluster_goodflag);
   //fChain->SetBranchAddress("cluster_npeaks", E->cluster_npeaks);
   //fChain->SetBranchAddress("cluster_ipeak", E->cluster_ipeak);
   //fChain->SetBranchAddress("cluster_timeofpeak", E->cluster_timeofpeak);
   //fChain->SetBranchAddress("muechk_ncand", E->muechk_ncand);
   //fChain->SetBranchAddress("muechk_tpeak", E->muechk_tpeak);
   //fChain->SetBranchAddress("muechk_bg", E->muechk_bg);
   //fChain->SetBranchAddress("muechk_mean", E->muechk_mean);
   //fChain->SetBranchAddress("muechk_excess", E->muechk_excess);
   //fChain->SetBranchAddress("muechk_signif", E->muechk_signif);
   //fChain->SetBranchAddress("muechk_icluster", E->muechk_icluster);
   //fChain->SetBranchAddress("trgoff", &E->trgoff);
   //fChain->SetBranchAddress("fqntwnd", &E->fqntwnd);
   fChain->SetBranchAddress("fqtwnd_iclstr", E->fqtwnd_iclstr);
   //fChain->SetBranchAddress("fqtwnd_npeak", E->fqtwnd_npeak);
   //fChain->SetBranchAddress("fqtwnd_prftt0", E->fqtwnd_prftt0);
   fChain->SetBranchAddress("fqtwnd_prftvtx", E->fqtwnd_prftvtx);
   //fChain->SetBranchAddress("fqtwnd", E->fqtwnd);
   //fChain->SetBranchAddress("fqtwnd_peakt0", E->fqtwnd_peakt0);
   //fChain->SetBranchAddress("fqtwnd_peakiness", E->fqtwnd_peakiness);
   fChain->SetBranchAddress("fqnse", &E->fqnse);
   fChain->SetBranchAddress("fqitwnd", E->fqitwnd);
   fChain->SetBranchAddress("fqipeak", E->fqipeak);
   fChain->SetBranchAddress("fqnhitpmt", E->fqnhitpmt);
   fChain->SetBranchAddress("fqtotq", E->fqtotq);
   //fChain->SetBranchAddress("fq0rtotmu", E->fq0rtotmu);
   //fChain->SetBranchAddress("fq0rnll", E->fq0rnll);
   fChain->SetBranchAddress("fqn50", E->fqn50);
   fChain->SetBranchAddress("fqq50", E->fqq50);
   fChain->SetBranchAddress("fq1rpcflg", E->fq1rpcflg);
   fChain->SetBranchAddress("fq1rmom", E->fq1rmom);
   fChain->SetBranchAddress("fq1rt0", E->fq1rt0);
   fChain->SetBranchAddress("fq1rtotmu", E->fq1rtotmu);
   fChain->SetBranchAddress("fq1rnll", E->fq1rnll);
   fChain->SetBranchAddress("fq1rvtx", E->fq1rvtx);
   fChain->SetBranchAddress("fq1rdir", E->fq1rdir);
   //fChain->SetBranchAddress("fq1rpar7", E->fq1rpar7);
   //fChain->SetBranchAddress("fqpi0pcflg", E->fqpi0pcflg);
   //fChain->SetBranchAddress("fqpi0mom1", E->fqpi0mom1);
   //fChain->SetBranchAddress("fqpi0mom2", E->fqpi0mom2);
   //fChain->SetBranchAddress("fqpi0momtot", E->fqpi0momtot);
   //fChain->SetBranchAddress("fqpi0dconv1", E->fqpi0dconv1);
   //fChain->SetBranchAddress("fqpi0dconv2", E->fqpi0dconv2);
   //fChain->SetBranchAddress("fqpi0t0", E->fqpi0t0);
   //fChain->SetBranchAddress("fqpi0totmu", E->fqpi0totmu);
   fChain->SetBranchAddress("fqpi0nll", E->fqpi0nll);
   //fChain->SetBranchAddress("fqpi0mass", E->fqpi0mass);
   //fChain->SetBranchAddress("fqpi0photangle", E->fqpi0photangle);
   //fChain->SetBranchAddress("fqpi0vtx", E->fqpi0vtx);
   //fChain->SetBranchAddress("fqpi0dir1", E->fqpi0dir1);
   //fChain->SetBranchAddress("fqpi0dir2", E->fqpi0dir2);
   //fChain->SetBranchAddress("fqpi0dirtot", E->fqpi0dirtot);
   //fChain->SetBranchAddress("fq2rpcflg", E->fq2rpcflg);
   //fChain->SetBranchAddress("fq2rmom1", E->fq2rmom1);
   //fChain->SetBranchAddress("fq2rmom2", E->fq2rmom2);
   //fChain->SetBranchAddress("fq2rpar71", E->fq2rpar71);
   //fChain->SetBranchAddress("fq2rpar72", E->fq2rpar72);
   //fChain->SetBranchAddress("fq2rt0", E->fq2rt0);
   //fChain->SetBranchAddress("fq2rtotmu", E->fq2rtotmu);
   fChain->SetBranchAddress("fq2rnll", E->fq2rnll);
   //fChain->SetBranchAddress("fq2rvtx", E->fq2rvtx);
   //fChain->SetBranchAddress("fq2rdir1", E->fq2rdir1);
   //fChain->SetBranchAddress("fq2rdir2", E->fq2rdir2);
   //fChain->SetBranchAddress("fq3rpcflg", E->fq3rpcflg);
   //fChain->SetBranchAddress("fq3rmom3", E->fq3rmom3);
   //fChain->SetBranchAddress("fq3rpar73", E->fq3rpar73);
   //fChain->SetBranchAddress("fq3rtotmu", E->fq3rtotmu);
   //fChain->SetBranchAddress("fq3rnll", E->fq3rnll);
   //fChain->SetBranchAddress("fq3rdir3", E->fq3rdir3);
   //fChain->SetBranchAddress("fq4rpcflg", E->fq4rpcflg);
   //fChain->SetBranchAddress("fq4rmom4", E->fq4rmom4);
   //fChain->SetBranchAddress("fq4rpar74", E->fq4rpar74);
   //fChain->SetBranchAddress("fq4rtotmu", E->fq4rtotmu);
   //fChain->SetBranchAddress("fq4rnll", E->fq4rnll);
   //fChain->SetBranchAddress("fq4rdir4", E->fq4rdir4);

    //fChain->SetBranchAddress("Npvc", &E->Npvc);
    //fChain->SetBranchAddress("Ipvc", E->Ipvc); 
    //fChain->SetBranchAddress("Ichvc", E->Ichvc); 
    //fChain->SetBranchAddress("Iorgvc", E->Iorgvc);
    //fChain->SetBranchAddress("Iflvc", E->Iflvc); 
    //fChain->SetBranchAddress("Abspvc", E->Abspvc); 
    //fChain->SetBranchAddress("Pvc", E->Pvc); 

    fChain->SetBranchAddress("nscndprt", &E->nscndprt);
    //fChain->SetBranchAddress("itrkscnd", E->itrkscnd); 
    fChain->SetBranchAddress("vtxscnd", E->vtxscnd); 
    fChain->SetBranchAddress("pscnd", E->pscnd);
    fChain->SetBranchAddress("iprtscnd", E->iprtscnd);
    fChain->SetBranchAddress("tscnd", E->tscnd); 
    //fChain->SetBranchAddress("iprntprt", E->iprntprt);
    fChain->SetBranchAddress("lmecscnd", E->lmecscnd);
    //fChain->SetBranchAddress("iprnttrk", E->iprnttrk);
    //fChain->SetBranchAddress("iorgprt", E->iorgprt); 
    //fChain->SetBranchAddress("pprnt", E->pprnt); 
    //fChain->SetBranchAddress("iflgscnd", E->iflgscnd);
  } 

  // STMU
  else if ((ireco==apfit&&recoSelect==bothreco) || recoSelect==allapfit) {
    
    //cout << "Setting APFIT branch addresses..." << endl;
    
   fChain->SetBranchAddress("nrun",      	 &E->stnrun);		 
   fChain->SetBranchAddress("nev",		 &E->stnev);		 
   fChain->SetBranchAddress("nsub",		 &E->stnsub);		 
   fChain->SetBranchAddress("potot",		 &E->stpotot);	 
   fChain->SetBranchAddress("range",		 &E->strange);	 
   fChain->SetBranchAddress("bgood",		 &E->stbgood);	 
   fChain->SetBranchAddress("egood",		 &E->stegood);	 
   fChain->SetBranchAddress("amom",		 &E->stamom);		 
   fChain->SetBranchAddress("vmom",		 &E->stvmom);		 
   fChain->SetBranchAddress("msrtot",            &E->stmsrtot);	 
   fChain->SetBranchAddress("rtot",		 &E->strtot);		 
   fChain->SetBranchAddress("pomax",		 &E->stpomax);	 
   fChain->SetBranchAddress("prmslg",            &E->stprmslg);	 
   fChain->SetBranchAddress("wtot",		 &E->stwtot);		 
   fChain->SetBranchAddress("angwall",		 &E->stangwall);	 
   fChain->SetBranchAddress("posmu",		  E->stposmu);		 
   fChain->SetBranchAddress("pose",		  E->stpose);		 
   fChain->SetBranchAddress("dirmu",		  E->stdirmu);		 
   fChain->SetBranchAddress("dire",		  E->stdire);		 
   fChain->SetBranchAddress("rangemc",		 &E->strangemc);	 
   fChain->SetBranchAddress("t",		 &E->stt);		 
   fChain->SetBranchAddress("n50",		 &E->stn50);		 
   fChain->SetBranchAddress("amomdcye",		 &E->stamomdcye);	 
   fChain->SetBranchAddress("q50",		 &E->stq50);		 
   fChain->SetBranchAddress("nall",		 &E->stnall);		 
   fChain->SetBranchAddress("nmue",		 &E->stnmue);		 
   fChain->SetBranchAddress("muetype",		 &E->stmuetype);	 
   fChain->SetBranchAddress("nday",		  E->stnday);		 
   fChain->SetBranchAddress("ntim",		  E->stntim);		 
   fChain->SetBranchAddress("rangev",		 &E->strangev);	 
   fChain->SetBranchAddress("amomv",		 &E->stamomv);	 
   fChain->SetBranchAddress("amomve",		 &E->stamomve);	 
   fChain->SetBranchAddress("posvmu",		  E->stposvmu);	 
   fChain->SetBranchAddress("posve",		  E->stposve);		 
   fChain->SetBranchAddress("dirvmu",		  E->stdirvmu);	 
   fChain->SetBranchAddress("dirve",		  E->stdirve);		 
   fChain->SetBranchAddress("dcyepidv",		 &E->stdcyepidv);	 
   fChain->SetBranchAddress("posemc",		  E->stposemc);	 
   fChain->SetBranchAddress("decaytv",		 &E->stdecaytv);  	 

    
    

  }
  
  //cout << "Done setting branch addresses" << endl;
  
}
#endif
