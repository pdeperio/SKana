#ifndef _EVENT_
#define _EVENT_

#include "TROOT.h"

struct Event{

/* Data Block */

   // Declaration of leave types
   Int_t           nring;
   UInt_t          nrun;
   Int_t           nev;
   Int_t           nsub;
   UInt_t          cate;
   Float_t         potot;
   UInt_t          nhit;
   Float_t         pomax;
   Float_t         potota;
   UInt_t          nhita;
   UInt_t          nhitac;
   Float_t         pomaxa;
   Float_t         wall;
   Float_t         evis;
/*
unused in this OfficialEventParser
   Float_t         rtsum;
   Float_t         rtmax;
   Float_t         wlen;
*/
   UInt_t          ip[10];   //[nring]
   Float_t         pos[3];
   Float_t         dir[10][3];   //[nring]
   Float_t         dirtot[3];
   Float_t         ang[10];   //[nring]
   Float_t         rtot[10];   //[nring]
   Float_t         amom[10];   //[nring]
   Float_t         rtote[10];   //[nring]
   Float_t         amome[10];   //[nring]
   Float_t         rtotm[10];   //[nring]
   Float_t         amomm[10];   //[nring]
   UInt_t          nsube;
   UInt_t          ndcy;
   UInt_t          ngate;
   UInt_t          nbye;
   Float_t         probms[10][6];   //[nring]
   Float_t         prmslg[10][6];   //[nring]
/*
unused in this OfficialEventParser
   Int_t           date[3];
   Int_t           time[4];
   Float_t         elpsday;
   Int_t           numpo[10];   //[nring]
   Float_t         apos[3];
   Float_t         adir[3];
   Float_t         aang;
   Float_t         agood;
   Float_t         wgain;
   Int_t           nbad;
   Int_t           nbada;
*/
   Float_t         msdir[10][3][6];   //[nring]
/*
   Float_t         amomp[10];
   Float_t         ange[10];
   Float_t         angm[10];
   Float_t         angp[10];
   Int_t           ntot[10];
   Float_t         probth[10][6];
   Float_t         probpt[10][6];
*/
   Float_t         pi0like[2];
   Float_t         pi0_e[2][2];
   Float_t         pi0_dir[2][2][3];
   Float_t         pi0mass[2];
/*
   Float_t         evisold;
   Float_t         evisoldxe;
   Float_t         evisnew;

*/

   Int_t           nmue;
   UInt_t          etype[10];   //[nmue]
   Float_t         etime[10];   //[nmue]
   Float_t         epos[10][3];   //[nmue]
   Float_t         edir[10][3];   //[nmue]
   Float_t         egood[10];   //[nmue]
   Float_t         ehit[10];   //[nmue]
   Float_t         mueprob[2];
/*
unused in this OfficialEventParser

   Int_t           Rnring;
   Float_t         Rdir[30][3];   //[Rnring]
   Float_t         Rang[30];   //[Rnring]
   UInt_t          Riring;
   Float_t         Rtwout[30];   //[Rnring]
   Float_t         Rtwith[30];   //[Rnring]
   Float_t         Alwout;
   Float_t         Alwith;
   Float_t         Qsmi;
   Float_t         Qsmo;
   Float_t         Qexi;
   Float_t         Qexo;
   Float_t         Pe5d;
   Float_t         En5d;
   Float_t         Eh5d;
   Float_t         Pe5do;
   Float_t         En5do;
   Float_t         Eh5do;
   Float_t         Rtadd;
   Float_t         Pdgeta;
   Float_t         Pd5d;
   Float_t         Pdthre;
   Float_t         Pd5do;
   Float_t         Delpd;
   Float_t         Ropena[30];   //[Rnring]
   Int_t           Maxth;
   Float_t         Pkang;
   Float_t         Qrfct;
   Float_t         Pdfct;
   Float_t         Pkfct;
   Float_t         Agfct;
*/

   Float_t         Dlfct;
/*
unused in this OfficialEventParser

   Int_t           Iflag;
   Float_t         Pmfct;
   Float_t         Imfct;
   Float_t         Rilike;
   Int_t           ri_ver;
   Float_t         ri_pid;
   Int_t           ri_nring;
   Float_t         ri_flag[10];   //[ri_nring]
   Float_t         ri_dlfct[10];   //[ri_nring]
   Float_t         ri_pdfct[10];   //[ri_nring]
   Float_t         ri_pkfct[10];   //[ri_nring]
   Float_t         ri_vafct[10];   //[ri_nring]
   Float_t         ri_total[10];   //[ri_nring]
   Float_t         ri_dir[10][3];   //[ri_nring]
   Float_t         ri_imfct[10];   //[ri_nring]
   Float_t         ri_pmfct[10];   //[ri_nring]
*/

   Int_t           npar;
   Float_t         wallv;
   UInt_t          ipv[50];   //[npar]
   Float_t         posv[3];
   Float_t         dirv[50][3];   //[npar]
   Float_t         pmomv[50];   //[npar]
/*
   Int_t           light_flag[50];   //[npar]
*/ 
   Int_t           npar2;
   Float_t         wallv2[50];   //[npar2]
   UInt_t          ipv2[50];   //[npar2]
   UInt_t          iorg[50];   //[npar2]
   Float_t         posv2[50][3];   //[npar2]
   Float_t         dirv2[50][3];   //[npar2]
   Float_t         pmomv2[50];   //[npar2]

   Int_t           numnu;
   Int_t           mode;
   Int_t           ipnu[50];   //[numnu]
   Float_t         pnu[50];   //[numnu]
   Float_t         dirnu[50][3];   //[numnu]
/*
unused in this OfficialEventParser
   Float_t         flxg[3];
   Float_t         flxh01[3];
   Float_t         kflux[4];
   Float_t         bs71[3];
   Float_t         bs74[3];
   Float_t         flxf[3];
   Float_t         flxh1d[3];
   Float_t         flxb03[3];
   Float_t         flxf03[3];
*/
   Float_t         flxh06[3];

   Float_t         live;
/* 
Unused in this OfficialEventParser
   Float_t         sacth;
   Float_t         sactg;
   Float_t         sacth1d;
   Int_t           scan[2];
   Float_t         dirtotepi[3];
   Float_t         dirtotenpi[3];
   Float_t         dirtotmue[3];
   Float_t         dirsum[3];
   Float_t         etot;
   Float_t         etotepi;
   Float_t         etotenpi;
   Float_t         etotmue;
   Float_t         oscweight[2][4];
*/
   Float_t         oscwgt;
/*
Unused in this OfficialEventParser
   Float_t         ent_pos[3];
   Float_t         ent_dir[3];
   Float_t         length;
   Float_t         tr_mom1;
   Float_t         A_ent_mom;
   Float_t         A_ent_pos[3];
   Float_t         A_ent_dir[3];
   Float_t         A_ext_mom;
   Float_t         A_ext_pos[3];
   Float_t         A_ext_dir[3];
   Float_t         Fit_pos[3];
   Float_t         Fit_dir[3];
   Float_t         Fit_len;
   Float_t         Fit_mom;
   Int_t           Fit_pid;
   Int_t           Um_ehit8m;
   Int_t           Um_ohit8m;
   Float_t         Um_qent;
   Float_t         Sh_chi1p;
   Float_t         Sh_delta;
   Float_t         Sh_mean;
   Float_t         Sh_meanq;
   Float_t         Tr_stop[3];
   Float_t         Tr_mom;
   Float_t         Tr_len;
   Float_t         Tr_len1;
   Int_t           Pid_flg;
   Float_t         Crs1;
   Float_t         Crs2;
   Int_t           iclass;
   Int_t           mu_class;
   Int_t           mu_dec;
   Float_t         mu_dir[3];
   Float_t         mu_pos[3];
   Float_t         mu_good;
   Int_t           history;
   Int_t           Pdst;
   Int_t           idoff;
   Float_t         anthit;
   Int_t           idseq;
   Float_t         tstfrac;
   Int_t           judge;
   Float_t         Upcrs1;
   Float_t         Upcrs2;
   Float_t         lst;
   Int_t           jd;
   Float_t         fjd;
   Float_t         alt;
   Float_t         azi;
   Float_t         ra;
   Float_t         dec;
   Float_t         glat;
   Float_t         glong;
   Int_t           nuceff_version;
   Int_t           charge_exchange;
   Int_t           absorbed;
   Int_t           multipi_gen;
   Int_t           scattering;
   Int_t           piless_dcy;
*/
   Int_t           cluster_ncand;
   Float_t         cluster_tstart[10];   //[cluster_ncand]
   Float_t         cluster_tend[10];   //[cluster_ncand]
   Int_t           cluster_nhits[10];   //[cluster_ncand]
   Float_t         cluster_totq[10];   //[cluster_ncand]
   Int_t           cluster_goodflag[10];   //[cluster_ncand]
   Int_t           cluster_npeaks[10][4];   //[cluster_ncand]
   Int_t           cluster_ipeak[10][4][10];   //[cluster_ncand]
   Float_t         cluster_timeofpeak[10][4][10];   //[cluster_ncand]
   Int_t           muechk_ncand[4];
   Float_t         muechk_tpeak[4][10];
   Int_t           muechk_bg[4][10];
   Float_t         muechk_mean[4][10];
   Float_t         muechk_excess[4][10];
   Float_t         muechk_signif[4][10];
   Int_t           muechk_icluster[4][10];
   Int_t           nse;
   Float_t         trgoff;
   Int_t           fqnhitpmt[10];   //[nse]
   Float_t         fqtotq[10];   //[nse]
   Float_t         fq0rtotmu[10];   //[nse]
   Float_t         fq0rnll[10];   //[nse]
   Int_t           fq1rpcflg[10][7];   //[nse]
   Float_t         fq1rmom[10][7];   //[nse]
   Float_t         fq1rt0[10][7];   //[nse]
   Float_t         fq1rtotmu[10][7];   //[nse]
   Float_t         fq1rnll[10][7];   //[nse]
   Float_t         fq1rvtx[10][7][3];   //[nse]
   Float_t         fq1rdir[10][7][3];   //[nse]
   Float_t         fq1rpar7[10][7];   //[nse]
   Int_t           fqpi0pcflg[2];
   Float_t         fqpi0mom1[2];
   Float_t         fqpi0mom2[2];
   Float_t         fqpi0momtot[2];
   Float_t         fqpi0dconv1[2];
   Float_t         fqpi0dconv2[2];
   Float_t         fqpi0t0[2];
   Float_t         fqpi0totmu[2];
   Float_t         fqpi0nll[2];
   Float_t         fqpi0mass[2];
   Float_t         fqpi0photangle[2];
   Float_t         fqpi0vtx[2][3];
   Float_t         fqpi0dir1[2][3];
   Float_t         fqpi0dir2[2][3];
   Float_t         fqpi0dirtot[2][3];
   Int_t           fq2rpcflg[2][2];
   Float_t         fq2rmom1[2][2];
   Float_t         fq2rmom2[2][2];
   Float_t         fq2rpar71[2][2];
   Float_t         fq2rpar72[2][2];
   Float_t         fq2rt0[2][2];
   Float_t         fq2rtotmu[2][2];
   Float_t         fq2rnll[2][2];
   Float_t         fq2rvtx[2][2][3];
   Float_t         fq2rdir1[2][2][3];
   Float_t         fq2rdir2[2][2][3];
   Int_t           fq3rpcflg[2][2][2];
   Float_t         fq3rmom3[2][2][2];
   Float_t         fq3rpar73[2][2][2];
   Float_t         fq3rtotmu[2][2][2];
   Float_t         fq3rnll[2][2][2];
   Float_t         fq3rdir3[2][2][2][3];
   Int_t           fq4rpcflg[2][2][2][2];
   Float_t         fq4rmom4[2][2][2][2];
   Float_t         fq4rpar74[2][2][2][2];
   Float_t         fq4rtotmu[2][2][2][2];
   Float_t         fq4rnll[2][2][2][2];
   Float_t         fq4rdir4[2][2][2][2][3];
   Int_t           fqtestn1r;
   Int_t           fqtest1rstage[30];   //[fqtestn1r]
   Int_t           fqtest1rse[30];   //[fqtestn1r]
   Int_t           fqtest1rpid[30];   //[fqtestn1r]
   Int_t           fqtest1rpcflg[30];   //[fqtestn1r]
   Float_t         fqtest1rmom[30];   //[fqtestn1r]
   Float_t         fqtest1rt0[30];   //[fqtestn1r]
   Float_t         fqtest1rtotmu[30];   //[fqtestn1r]
   Float_t         fqtest1rnll[30];   //[fqtestn1r]
   Float_t         fqtest1rvtx[30][3];   //[fqtestn1r]
   Float_t         fqtest1rdir[30][3];   //[fqtestn1r]
   Float_t         fqtest1rpar7[30];   //[fqtestn1r]
   Int_t           fqtestnpi0;
   Int_t           fqtestpi0stage[30];   //[fqtestnpi0]
   Int_t           fqtestpi0pcflg[30];   //[fqtestnpi0]
   Float_t         fqtestpi0mom1[30];   //[fqtestnpi0]
   Float_t         fqtestpi0mom2[30];   //[fqtestnpi0]
   Float_t         fqtestpi0momtot[30];   //[fqtestnpi0]
   Float_t         fqtestpi0dconv1[30];   //[fqtestnpi0]
   Float_t         fqtestpi0dconv2[30];   //[fqtestnpi0]
   Float_t         fqtestpi0t0[30];   //[fqtestnpi0]
   Float_t         fqtestpi0totmu[30];   //[fqtestnpi0]
   Float_t         fqtestpi0nll[30];   //[fqtestnpi0]
   Float_t         fqtestpi0mass[30];   //[fqtestnpi0]
   Float_t         fqtestpi0photangle[30];   //[fqtestnpi0]
   Float_t         fqtestpi0vtx[30][3];   //[fqtestnpi0]
   Float_t         fqtestpi0dir1[30][3];   //[fqtestnpi0]
   Float_t         fqtestpi0dir2[30][3];   //[fqtestnpi0]
   Float_t         fqtestpi0dirtot[30][3];   //[fqtestnpi0]

   Float_t         stmuent[3];
   Float_t         stmudir[3];
   Float_t         stmugood;
   Float_t         stmuqent;

   Int_t           Npvc;
   Int_t           Ipvc[100];   //[Npvc]
   Int_t           Ichvc[100];   //[Npvc]
   Int_t           Iorgvc[100];   //[Npvc]
   Int_t           Iflvc[100];   //[Npvc]
   Float_t         Abspvc[100];   //[Npvc]
   Float_t         Pvc[100][3];   //[Npvc]


   Int_t           nscndprt;
   Int_t           itrkscnd[200];   //[nscndprt]
   Float_t         vtxscnd[200][3];   //[nscndprt]
   Float_t         pscnd[200][3];   //[nscndprt]
   Int_t           iprtscnd[200];   //[nscndprt]
   Float_t         tscnd[200];   //[nscndprt]
   Int_t           iprntprt[200];   //[nscndprt]
   Int_t           lmecscnd[200];   //[nscndprt]
   Int_t           iprnttrk[200];   //[nscndprt]
   Int_t           iorgprt[200];   //[nscndprt]
   Float_t         pprnt[200][3];   //[nscndprt]
   Int_t           iflgscnd[200];   //[nscndprt]
};

#endif
