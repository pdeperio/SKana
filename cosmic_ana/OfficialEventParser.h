#ifndef OfficialEventParser_h
#define OfficialEventParser_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "Event.h"
#include <iostream>
#include <fstream>
#include <string>
#include <TH2.h>
#include "global.h"
//#include "pdg_codes.h"

using namespace std;

const int nRecos = 2;
enum reco_enum {apfit, fitqun};

namespace fiTQun_parameters{
  enum particleIndices {ig, ie, imu, ipip, ikp, ipr, icg};
  enum TrackIndices {itrkg, itrke, itrkmu, itrkpip, itrkkp, itrkp, itrkcg, itrkpi0, itrkpipknk};

}
using namespace fiTQun_parameters;


//
// SK tank parameter from geotnk.h 
//
const Double_t DITKTK = 3930. ; 
const Double_t HITKTK=4140. ; 
const Double_t DIINTK=3380. ; 
const Double_t HIINTK=3620. ; 
const Double_t HIWAL=3605.7 ; 
const Double_t TTATTK=260. ; 
const Double_t TBATTK=260. ; 
const Double_t TWATTK=275. ; 
const Double_t RTKTK = DITKTK/2. ; 
const Double_t ZPTKTK=HITKTK/2. ; 
const Double_t ZMTKTK=-HITKTK/2. ; 
const Double_t RINTK=DIINTK/2. ; 
const Double_t ZPINTK=HIINTK/2. ; 
const Double_t ZMINTK=-HIINTK/2. ; 
const Double_t ZPWAL=HIWAL/2. ; 
const Double_t ZMWAL=-HIWAL/2. ; 
const Double_t ZCNTTK=HITKTK/2. ;
const Double_t RMED  = 55. ;
const Double_t ZMED  = 55. ;
  
const int nSE = 2;
const int nSample = 2;
const int nCuts = 2;

const int nMoms = 20;  //  mu, e
const double momLow[nSE] = {0, 0};
const double momHigh[nSE] = {1000, 200};
const double momStep[nSE] = {50, 10}; // (momHigh-momLow)/nMoms

class OfficialEventParser{
private:
        Event *E;

	TFile * efile[nSample][nSE];

	TH1D * h_cosmic_evis[nSample][nRecos];
	TH1D * h_cosmic_nhit_over_time[nSample][nRecos];
	TH1D * h_cosmic_nmue[nSample][nRecos];
	TH1D * h_cosmic_z[nSample][nRecos];
	TH1D * h_cosmic_wall[nSample][nSE][nRecos];
	TH1D * h_cosmic_pcflag[nSample][nSE][nRecos];
	TH1D * h_cosmic_zenith[nSample][nSE][nRecos];
	TH1D * h_cosmic_mom[nSample][nSE][nRecos][nCuts];
	TH1D * h_cosmic_mom_res[nSample][nSE][nRecos][nCuts];
	TH1D * h_cosmic_tvtx_res[nSample][nSE][nRecos][nCuts];
	TH1D * h_cosmic_lvtx_res[nSample][nSE][nRecos][nCuts];
	TH1D * h_cosmic_dir_res[nSample][nSE][nRecos][nCuts]; 
	TH2D * h_cosmic_range_vs_mom[nSample][nRecos][nCuts];
	TH1D * h_cosmic_range_over_mom[nSample][nRecos][nCuts];
	TH1D * h_cosmic_r2[nSample][nRecos];
	TH2D * h_cosmic_lnLrat1R_vs_mom[nSample][nSE][nRecos];
	TH1D * h_cosmic_lnLrat1R[nSample][nSE][nRecos];
	TH1D * h_cosmic_lnLrat1R_momsep[nSample][nSE][nRecos][nMoms];
	TH1D * h_cosmic_lnLrat1R2R[nSample][nRecos];
	TH1D * h_cosmic_lnLrat1R2R_momsep[nSample][nRecos][nMoms];
	TH2D * h_cosmic_lnLrat1R2R_vs_mom[nSample][nRecos];
	TH2D * h_cosmic_lnLrat1R2R_vs_pid[nSample][nRecos];
	TH1D * h_cosmic_ddir[nSample][nRecos];
	TH1D * h_cosmic_lifetime[nSample][nRecos][nCuts];
	TH1D * h_eff[nSample][nRecos];


	bool fcmc;      // true : MC, false : data
	int  _version;  // SK detector version
	bool _osc;      // flag to oscillate MC or not
	int _datatype; 
	int _towall_dir; // towall_dir for this event is filled.
	bool _from_wall_flag; // if this flag is true, calculate from_wall distance instead of to_wall.
	double _datalivetime; // livetime of data

	enum towall_dir {out_of_tank=0, to_side, to_top, to_bottom};
	enum find_true_particle { find_e, find_mu};

public :

	OfficialEventParser();
	void InitHistograms();
	void InitEvent();
	void Parse(int ireco=0);
	void Print();
	void Terminate();
	void operator()(  TChain *fChain );
	
	void FillCosmic(int ireco);

	double weight();

	bool isCosmic(int ireco);
	bool is1R(int ireco);
	bool hasRvtxInFid(int ireco);
	bool hasZvtxInFid(int ireco);

        void set_version( int x)  { _version = x; }
        int  get_version( void)   { return _version;}
	bool get_osc(void)        { return _osc;}
	void set_osc(bool x)      { _osc = x;}
	int  get_datatype(void)   { return _datatype;}
	void set_datatype(int x)  { _datatype = x;}

	float cut_nhitac(); // return cut value of nhitac for sk1,sk2, and sk3
	float cut_potot();
	int cut_ehit( int i );
	
	double get_towall(int nth_ring);
	double get_fromwall(int nth_ring);

	double getLnLrat1R(int ireco, int ise);
	double getLnLrat1R2R(int ireco, int &fillPID);
	int get1rType(int ireco,int ise);
	double getMin2Rnll(int ireco, int &fillPID);
	
	
	void  set_towall_dir( int x){ _towall_dir = x;}
	int   get_towall_dir( void ){ return _towall_dir ;}
	double calc_towall(double x, double y, double z, int ireco=0);
	
	void  set_from_wall_flag ( bool x ) { _from_wall_flag = x;}
	bool  get_from_wall_flag ( void ) { return _from_wall_flag;}
	
	void  set_datalivetime ( double x ) { _datalivetime = x; cout << _datalivetime << endl;}
	double get_datalivetime( void ) { return _datalivetime;} 

	int get_num_ring(int ireco);
	
	float *get_vertex(int ireco, int ise=0);
	float *get_1Rdir(int ireco, int ise=0);
	float get_1Rmom(int ireco, int ise=0);
	float get_dcye_time(int ireco, int ise=0);

	float get_nhits_over_time(int ireco);

	float get_evis(int ireco);
	float get_wall(float *vertex);

	int   muedcy(int ireco);

};
#endif

#ifdef __OFFICIALEVENTPARSER_CXX
void OfficialEventParser::operator()(  TChain *fChain ){

/*
	I listed up the only common parameters between data and MC. ... Naho 
*/
   fChain->SetBranchAddress("nring", &E->nring);
   fChain->SetBranchAddress("nrun", &E->nrun);
   fChain->SetBranchAddress("nev", &E->nev);
   fChain->SetBranchAddress("nsub", &E->nsub);
   fChain->SetBranchAddress("cate", &E->cate);
   fChain->SetBranchAddress("potot", &E->potot);
   fChain->SetBranchAddress("nhit", &E->nhit);
   fChain->SetBranchAddress("pomax", &E->pomax);
   fChain->SetBranchAddress("potota", &E->potota);
   fChain->SetBranchAddress("nhita", &E->nhita);
   fChain->SetBranchAddress("nhitac", &E->nhitac);
   fChain->SetBranchAddress("pomaxa", &E->pomaxa);
   fChain->SetBranchAddress("wall", &E->wall);
   fChain->SetBranchAddress("evis", &E->evis);
/*
unused in this OfficialEventParser
   fChain->SetBranchAddress("rtsum", &E->rtsum);
   fChain->SetBranchAddress("rtmax", &E->rtmax);
   fChain->SetBranchAddress("wlen", &E->wlen);
*/

   fChain->SetBranchAddress("ip", E->ip);
   fChain->SetBranchAddress("pos", E->pos);
   fChain->SetBranchAddress("dir", E->dir);
   fChain->SetBranchAddress("dirtot", E->dirtot);
   fChain->SetBranchAddress("ang", E->ang);
   fChain->SetBranchAddress("rtot", E->rtot);
   fChain->SetBranchAddress("amom", E->amom);
   fChain->SetBranchAddress("rtote", E->rtote);
   fChain->SetBranchAddress("amome", E->amome);
   fChain->SetBranchAddress("rtotm", E->rtotm);
   fChain->SetBranchAddress("amomm", E->amomm);
   fChain->SetBranchAddress("nsube", &E->nsube);
   fChain->SetBranchAddress("ndcy", &E->ndcy);
   fChain->SetBranchAddress("ngate", &E->ngate);
   fChain->SetBranchAddress("nbye", &E->nbye);
   fChain->SetBranchAddress("probms", E->probms);
   fChain->SetBranchAddress("prmslg", E->prmslg);
/*
unused in this OfficialEventParser
   fChain->SetBranchAddress("date", E->date);
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
   fChain->SetBranchAddress("msdir", E->msdir);
/*
   fChain->SetBranchAddress("amomp", E->amomp);
   fChain->SetBranchAddress("ange", E->ange);
   fChain->SetBranchAddress("angm", E->angm);
   fChain->SetBranchAddress("angp", E->angp);
   fChain->SetBranchAddress("ntot", E->ntot);
   fChain->SetBranchAddress("probth", E->probth);
   fChain->SetBranchAddress("probpt", E->probpt);
*/
   fChain->SetBranchAddress("pi0like", E->pi0like);
   fChain->SetBranchAddress("pi0_e", E->pi0_e);
   fChain->SetBranchAddress("pi0_dir", E->pi0_dir);
   fChain->SetBranchAddress("pi0mass", E->pi0mass);
/*
   fChain->SetBranchAddress("evisold", E->evisold);
   fChain->SetBranchAddress("evisoldxe", E->evisoldxe);
   fChain->SetBranchAddress("evisnew", E->evisnew);
*/
   fChain->SetBranchAddress("nmue", &E->nmue);
   fChain->SetBranchAddress("etype", E->etype);
   fChain->SetBranchAddress("etime", E->etime);
   fChain->SetBranchAddress("epos", E->epos);
   fChain->SetBranchAddress("edir", E->edir);
   fChain->SetBranchAddress("egood", E->egood);
   fChain->SetBranchAddress("ehit", E->ehit);
   fChain->SetBranchAddress("mueprob", E->mueprob);

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


   fChain->SetBranchAddress("Dlfct", &E->Dlfct);

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
   fChain->SetBranchAddress("npar2", &E->npar2);
   fChain->SetBranchAddress("wallv2", E->wallv2);
   fChain->SetBranchAddress("ipv2", E->ipv2);
   fChain->SetBranchAddress("iorg", E->iorg);
   fChain->SetBranchAddress("posv2", E->posv2);
   fChain->SetBranchAddress("dirv2", E->dirv2);
   fChain->SetBranchAddress("pmomv2", E->pmomv2);
   fChain->SetBranchAddress("numnu", &E->numnu);
   fChain->SetBranchAddress("mode", &E->mode);
   fChain->SetBranchAddress("ipnu", E->ipnu);
   fChain->SetBranchAddress("pnu", E->pnu);
   fChain->SetBranchAddress("dirnu", E->dirnu);
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
   fChain->SetBranchAddress("flxh06", E->flxh06);
   fChain->SetBranchAddress("live", &E->live);
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
   fChain->SetBranchAddress("oscwgt", &E->oscwgt);

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
   fChain->SetBranchAddress("cluster_ncand", &E->cluster_ncand);
   fChain->SetBranchAddress("cluster_tstart", E->cluster_tstart);
   fChain->SetBranchAddress("cluster_tend", E->cluster_tend);
   fChain->SetBranchAddress("cluster_nhits", E->cluster_nhits);
   fChain->SetBranchAddress("cluster_totq", E->cluster_totq);
   fChain->SetBranchAddress("cluster_goodflag", E->cluster_goodflag);
   fChain->SetBranchAddress("cluster_npeaks", E->cluster_npeaks);
   fChain->SetBranchAddress("cluster_ipeak", E->cluster_ipeak);
   fChain->SetBranchAddress("cluster_timeofpeak", E->cluster_timeofpeak);

   fChain->SetBranchAddress("nse", &E->nse);
   fChain->SetBranchAddress("trgoff", &E->trgoff);
   fChain->SetBranchAddress("fqnhitpmt", E->fqnhitpmt);
   fChain->SetBranchAddress("fqtotq", E->fqtotq);
   fChain->SetBranchAddress("fq0rtotmu", E->fq0rtotmu);
   fChain->SetBranchAddress("fq0rnll", E->fq0rnll);
   fChain->SetBranchAddress("fq1rpcflg", E->fq1rpcflg);
   fChain->SetBranchAddress("fq1rmom", E->fq1rmom);
   fChain->SetBranchAddress("fq1rt0", E->fq1rt0);
   fChain->SetBranchAddress("fq1rtotmu", E->fq1rtotmu);
   fChain->SetBranchAddress("fq1rnll", E->fq1rnll);
   fChain->SetBranchAddress("fq1rvtx", E->fq1rvtx);
   fChain->SetBranchAddress("fq1rdir", E->fq1rdir);
   fChain->SetBranchAddress("fq1rpar7", E->fq1rpar7);
   fChain->SetBranchAddress("fqpi0pcflg", E->fqpi0pcflg);
   fChain->SetBranchAddress("fqpi0mom1", E->fqpi0mom1);
   fChain->SetBranchAddress("fqpi0mom2", E->fqpi0mom2);
   fChain->SetBranchAddress("fqpi0momtot", E->fqpi0momtot);
   fChain->SetBranchAddress("fqpi0dconv1", E->fqpi0dconv1);
   fChain->SetBranchAddress("fqpi0dconv2", E->fqpi0dconv2);
   fChain->SetBranchAddress("fqpi0t0", E->fqpi0t0);
   fChain->SetBranchAddress("fqpi0totmu", E->fqpi0totmu);
   fChain->SetBranchAddress("fqpi0nll", E->fqpi0nll);
   fChain->SetBranchAddress("fqpi0mass", E->fqpi0mass);
   fChain->SetBranchAddress("fqpi0photangle", E->fqpi0photangle);
   fChain->SetBranchAddress("fqpi0vtx", E->fqpi0vtx);
   fChain->SetBranchAddress("fqpi0dir1", E->fqpi0dir1);
   fChain->SetBranchAddress("fqpi0dir2", E->fqpi0dir2);
   fChain->SetBranchAddress("fqpi0dirtot", E->fqpi0dirtot);
   fChain->SetBranchAddress("fq2rpcflg", E->fq2rpcflg);
   fChain->SetBranchAddress("fq2rmom1", E->fq2rmom1);
   fChain->SetBranchAddress("fq2rmom2", E->fq2rmom2);
   fChain->SetBranchAddress("fq2rpar71", E->fq2rpar71);
   fChain->SetBranchAddress("fq2rpar72", E->fq2rpar72);
   fChain->SetBranchAddress("fq2rt0", E->fq2rt0);
   fChain->SetBranchAddress("fq2rtotmu", E->fq2rtotmu);
   fChain->SetBranchAddress("fq2rnll", E->fq2rnll);
   fChain->SetBranchAddress("fq2rvtx", E->fq2rvtx);
   fChain->SetBranchAddress("fq2rdir1", E->fq2rdir1);
   fChain->SetBranchAddress("fq2rdir2", E->fq2rdir2);
   fChain->SetBranchAddress("fq3rpcflg", E->fq3rpcflg);
   fChain->SetBranchAddress("fq3rmom3", E->fq3rmom3);
   fChain->SetBranchAddress("fq3rpar73", E->fq3rpar73);
   fChain->SetBranchAddress("fq3rtotmu", E->fq3rtotmu);
   fChain->SetBranchAddress("fq3rnll", E->fq3rnll);
   fChain->SetBranchAddress("fq3rdir3", E->fq3rdir3); 
   fChain->SetBranchAddress("fq4rpcflg", E->fq4rpcflg);
   fChain->SetBranchAddress("fq4rmom4", E->fq4rmom4); 
   fChain->SetBranchAddress("fq4rpar74", E->fq4rpar74); 
   fChain->SetBranchAddress("fq4rtotmu", E->fq4rtotmu); 
   fChain->SetBranchAddress("fq4rnll", E->fq4rnll);
   fChain->SetBranchAddress("fq4rdir4", E->fq4rdir4); 
   fChain->SetBranchAddress("fqtestn1r", &E->fqtestn1r); 
   fChain->SetBranchAddress("fqtest1rstage", E->fqtest1rstage); 
   fChain->SetBranchAddress("fqtest1rse", E->fqtest1rse);
   fChain->SetBranchAddress("fqtest1rpid", E->fqtest1rpid); 
   fChain->SetBranchAddress("fqtest1rpcflg", E->fqtest1rpcflg); 
   fChain->SetBranchAddress("fqtest1rmom", E->fqtest1rmom);
   fChain->SetBranchAddress("fqtest1rt0", E->fqtest1rt0); 
   fChain->SetBranchAddress("fqtest1rtotmu", E->fqtest1rtotmu);
   fChain->SetBranchAddress("fqtest1rnll", E->fqtest1rnll);
   fChain->SetBranchAddress("fqtest1rvtx", E->fqtest1rvtx); 
   fChain->SetBranchAddress("fqtest1rdir", E->fqtest1rdir);
   fChain->SetBranchAddress("fqtest1rpar7", E->fqtest1rpar7);
   fChain->SetBranchAddress("fqtestnpi0", &E->fqtestnpi0); 
   fChain->SetBranchAddress("fqtestpi0stage", E->fqtestpi0stage); 
   fChain->SetBranchAddress("fqtestpi0pcflg", E->fqtestpi0pcflg); 
   fChain->SetBranchAddress("fqtestpi0mom1", E->fqtestpi0mom1); 
   fChain->SetBranchAddress("fqtestpi0mom2", E->fqtestpi0mom2);
   fChain->SetBranchAddress("fqtestpi0momtot", E->fqtestpi0momtot); 
   fChain->SetBranchAddress("fqtestpi0dconv1", E->fqtestpi0dconv1); 
   fChain->SetBranchAddress("fqtestpi0dconv2", E->fqtestpi0dconv2); 
   fChain->SetBranchAddress("fqtestpi0t0", E->fqtestpi0t0);
   fChain->SetBranchAddress("fqtestpi0totmu", E->fqtestpi0totmu); 
   fChain->SetBranchAddress("fqtestpi0nll", E->fqtestpi0nll);
   fChain->SetBranchAddress("fqtestpi0mass", E->fqtestpi0mass);
   fChain->SetBranchAddress("fqtestpi0photangle", E->fqtestpi0photangle);
   fChain->SetBranchAddress("fqtestpi0vtx", E->fqtestpi0vtx);
   fChain->SetBranchAddress("fqtestpi0dir1", E->fqtestpi0dir1);
   fChain->SetBranchAddress("fqtestpi0dir2", E->fqtestpi0dir2);
   fChain->SetBranchAddress("fqtestpi0dirtot", E->fqtestpi0dirtot); 

   fChain->SetBranchAddress("stmuent", E->stmuent);
   fChain->SetBranchAddress("stmudir", E->stmudir);
   fChain->SetBranchAddress("stmugood", &E->stmugood);
   fChain->SetBranchAddress("stmuqent", &E->stmuqent);


   fChain->SetBranchAddress("Npvc", &E->Npvc);
   fChain->SetBranchAddress("Ipvc", E->Ipvc); 
   fChain->SetBranchAddress("Ichvc", E->Ichvc); 
   fChain->SetBranchAddress("Iorgvc", E->Iorgvc);
   fChain->SetBranchAddress("Iflvc", E->Iflvc); 
   fChain->SetBranchAddress("Abspvc", E->Abspvc); 
   fChain->SetBranchAddress("Pvc", E->Pvc); 

   fChain->SetBranchAddress("nscndprt", &E->nscndprt);
   fChain->SetBranchAddress("itrkscnd", E->itrkscnd); 
   fChain->SetBranchAddress("vtxscnd", E->vtxscnd); 
   fChain->SetBranchAddress("pscnd", E->pscnd);
   fChain->SetBranchAddress("iprtscnd", E->iprtscnd);
   fChain->SetBranchAddress("tscnd", E->tscnd); 
   fChain->SetBranchAddress("iprntprt", E->iprntprt);
   fChain->SetBranchAddress("lmecscnd", E->lmecscnd);
   fChain->SetBranchAddress("iprnttrk", E->iprnttrk);
   fChain->SetBranchAddress("iorgprt", E->iorgprt); 
   fChain->SetBranchAddress("pprnt", E->pprnt); 
   fChain->SetBranchAddress("iflgscnd", E->iflgscnd);
}
#endif
