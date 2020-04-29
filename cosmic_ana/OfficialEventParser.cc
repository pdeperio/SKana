#define __OFFICIALEVENTPARSER_CXX

#include "OfficialEventParser.h"
#include "Event.h"
#include <string>
#include "global.h"
#include <cmath>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <map>
#include <TMath.h>
#include <TVector3.h>
//#include "pdg_codes.h"

#define __FC_RUN_SELECTION__  true 
//#define __FC_RUN_SELECTION__  ( E->nrun >= 66251 && E->nrun <= 66945 )
//#define __FC_RUN_SELECTION__  ( E->nrun <= 66955  ) 

//    This needs to be included later

/*
$SKOFL_ROOT/inc/geotnkC.h
*/

using namespace std;
using namespace global;

OfficialEventParser::OfficialEventParser(){
  E = new Event();
  set_version(-1); // initialize
  set_osc(false); // not oscillate
  set_datatype(fcData); // fcdata, fcpc, etc. written in global.cc

}

void OfficialEventParser::InitEvent(){
  set_from_wall_flag(false);
}

void OfficialEventParser::InitHistograms(){
  // Fill histograms in a root file.
  string outputnameroot="my_";
  if (get_version()==1) outputnameroot+="sk1";
  else if (get_version()==2) outputnameroot+="sk2";
  else if (get_version()==3) outputnameroot+="sk3";
  else if (get_version()==4) outputnameroot+="sk4";
  outputnameroot+="_cosmics";


  for (int ireco=0; ireco<nRecos; ireco++) {
	
    char *h_name[2] = {"h","hmc"};
    for (int i=0; i<nSample; i++) {

      string outputname=outputnameroot+Form("_reco%d_",ireco);
      outputname+=h_name[i];
      outputname+=".root";
      efile[i][ireco] = new TFile(outputname.c_str(),"recreate");

      for (int ise=0; ise<nSE; ise++) {
	if (ise==0) {
	  h_cosmic_lnLrat1R_vs_mom[i][ise][ireco]     = new TH2D(Form("%s_cosmic_lnLrat1R_vs_mom_ise_%d_reco%d",h_name[i],ise,ireco), Form("%s_cosmic_lnLrat1R_vs_mom_ise_%d_reco%d;p_{e} (MeV/c); ln(L_{e}/L_{#mu})",h_name[i],ise,ireco), 250, 0, 1000, 1000, -2000, 2000);
	  h_cosmic_lnLrat1R[i][ise][ireco]     = new TH1D(Form("%s_cosmic_lnLrat1R_ise%d_reco%d",h_name[i],ise,ireco), Form("%s_cosmic_lnLrat1R_ise%d_reco%d; ln(L_{e}/L_{#mu})",h_name[i],ise,ireco), 200, -2000, 2000);
	} else {
	  h_cosmic_lnLrat1R_vs_mom[i][ise][ireco]     = new TH2D(Form("%s_cosmic_lnLrat1R_vs_mom_ise%d_reco%d",h_name[i],ise,ireco), Form("%s_cosmic_lnLrat1R_vs_mom_ise%d_reco%d;p_{e} (MeV/c); ln(L_{e}/L_{#mu})",h_name[i],ise,ireco), 100, 0, 200, 500, -1000, 1000);
	  h_cosmic_lnLrat1R[i][ise][ireco]     = new TH1D(Form("%s_cosmic_lnLrat1R_ise%d_reco%d",h_name[i],ise,ireco), Form("%s_cosmic_lnLrat1R_ise%d_reco%d;ln(L_{e}/L_{#mu})",h_name[i],ise,ireco), 200, -1000, 1000);
	}
	
	h_cosmic_lnLrat1R_vs_mom[i][ise][ireco]->Sumw2();
	h_cosmic_lnLrat1R[i][ise][ireco]->Sumw2();

	if (ise==0) h_cosmic_zenith[i][ise][ireco]     = new TH1D(Form("%s_cosmic_zenith_ise%d_reco%d",h_name[i],ise,ireco),  Form("%s_cosmic_zenith_ise%d_reco%d;cos(#theta_{zenith})",h_name[i],ise,ireco), 50, 0, 1);
	else
	  h_cosmic_zenith[i][ise][ireco]     = new TH1D(Form("%s_cosmic_zenith_ise%d_reco%d",h_name[i],ise,ireco),  Form("%s_cosmic_zenith_ise%d_reco%d;cos(#theta_{zenith})",h_name[i],ise,ireco), 20, -1, 1);
	h_cosmic_zenith[i][ise][ireco]	       ->Sumw2();


	h_cosmic_wall[i][ise][ireco]     = new TH1D(Form("%s_cosmic_wall_ise%d_reco%d",h_name[i],ise,ireco), Form("%s_cosmic_wall_ise%d_reco%d;Wall (cm)",h_name[i],ise,ireco)  , 100, 0, 500);
	h_cosmic_wall[i][ise][ireco]  	       ->Sumw2();

	h_cosmic_pcflag[i][ise][ireco]     = new TH1D(Form("%s_cosmic_pcflag_ise%d_reco%d",h_name[i],ise,ireco),Form("%s_cosmic_pcflag_ise%d_reco%d",h_name[i],ise,ireco) , 4, -2, 2);
	h_cosmic_pcflag[i][ise][ireco]         ->Sumw2();

	for (int imom=0; imom<nMoms; imom++) {
	  
	  double momLowEdge = (momLow[ise]+imom*momStep[ise]);
	  double momHighEdge = (momLow[ise]+(imom+1)*momStep[ise]);

	  if (ise==0)
	    h_cosmic_lnLrat1R_momsep[i][ise][ireco][imom] = new TH1D(Form("%s_cosmic_lnLrat1R_ise%d_reco%d_momsep%03d_%03d",h_name[i],ise,ireco,(int)momLowEdge,(int)momHighEdge), Form("%s_cosmic_lnLrat1R_ise%d_reco%d_momsep%03d_%03d; ln(L_{e}/L_{#mu})",h_name[i],ise,ireco,(int)momLowEdge,(int)momHighEdge), 200, -2000, 2000);		  
	  else 
	    h_cosmic_lnLrat1R_momsep[i][ise][ireco][imom] = new TH1D(Form("%s_cosmic_lnLrat1R_ise%d_reco%d_momsep%03d_%03d",h_name[i],ise,ireco,(int)momLowEdge,(int)momHighEdge), Form("%s_cosmic_lnLrat1R_ise%d_reco%d_momsep%03d_%03d; ln(L_{e}/L_{#mu})",h_name[i],ise,ireco,(int)momLowEdge,(int)momHighEdge), 200, -1000, 1000);

	  h_cosmic_lnLrat1R_momsep[i][ise][ireco][imom]->Sumw2();
	  
	}
	
	for (int icut=0; icut<nCuts; icut++) {
	  
	  if (ise==0)
	    h_cosmic_mom[i][ise][ireco][icut]     = new TH1D(Form("%s_cosmic_mom_ise%d_reco%d_cut%d",h_name[i],ise,ireco,icut), Form("%s_cosmic_mom_ise%d_reco%d_cut%d; Momentum (MeV/c)",h_name[i],ise,ireco,icut)   , 100, 0, 1000);
	  else 
	    h_cosmic_mom[i][ise][ireco][icut]     = new TH1D(Form("%s_cosmic_mom_ise%d_reco%d_cut%d",h_name[i],ise,ireco,icut), Form("%s_cosmic_mom_ise%d_reco%d_cut%d; Momentum (MeV/c)",h_name[i],ise,ireco,icut), 50, 0, 100);
	  h_cosmic_mom[i][ise][ireco][icut] ->Sumw2();

	  h_cosmic_mom_res[i][ise][ireco][icut]     = new TH1D(Form("%s_cosmic_mom_res_ise%d_reco%d_cut%d",h_name[i],ise,ireco,icut), Form("%s_cosmic_mom_res_ise%d_reco%d_cut%d; Momentum Resolution (%)",h_name[i],ise,ireco,icut), 200, -20, 20);
	  h_cosmic_mom_res[i][ise][ireco][icut] ->Sumw2();
	  
	  h_cosmic_tvtx_res[i][ise][ireco][icut]     = new TH1D(Form("%s_cosmic_tvtx_res_ise%d_reco%d_cut%d",h_name[i],ise,ireco,icut), Form("%s_cosmic_tvtx_res_ise%d_reco%d_cut%d; Tvtx Resolution (cm)",h_name[i],ise,ireco,icut), 200, 0, 100);
	  h_cosmic_tvtx_res[i][ise][ireco][icut] ->Sumw2();

	  h_cosmic_lvtx_res[i][ise][ireco][icut]     = new TH1D(Form("%s_cosmic_lvtx_res_ise%d_reco%d_cut%d",h_name[i],ise,ireco,icut), Form("%s_cosmic_lvtx_res_ise%d_reco%d_cut%d; Lvtx Resolution (cm)",h_name[i],ise,ireco,icut), 200, -100, 100);
	  h_cosmic_lvtx_res[i][ise][ireco][icut] ->Sumw2();

	  h_cosmic_dir_res[i][ise][ireco][icut]     = new TH1D(Form("%s_cosmic_dir_res_ise%d_reco%d_cut%d",h_name[i],ise,ireco,icut), Form("%s_cosmic_dir_res_ise%d_reco%d_cut%d; Direction Resolution (%)",h_name[i],ise,ireco,icut), 100, 0, 30);
	  h_cosmic_dir_res[i][ise][ireco][icut] ->Sumw2();

	}

      }
      
      h_cosmic_evis[i][ireco]     = new TH1D(Form("%s_cosmic_evis_reco%d",h_name[i],ireco),Form("%s_cosmic_evis_reco%d; Evis (MeV)",h_name[i],ireco) , 50, 0, 1000);
      h_cosmic_nhit_over_time[i][ireco]     = new TH1D(Form("%s_cosmic_nhit_over_time_reco%d",h_name[i],ireco), Form("%s_cosmic_nhit_over_time_reco%d; Hits/Time Window (/ns)",h_name[i],ireco), 250, 0, 5);
      h_cosmic_nmue[i][ireco]     = new TH1D(Form("%s_cosmic_nmue_reco%d",h_name[i],ireco), Form("%s_cosmic_nmue_reco%d; Number of Decay-e",h_name[i],ireco), 10, 0, 10);
      h_cosmic_z[i][ireco]     = new TH1D(Form("%s_cosmic_z_reco%d",h_name[i],ireco),  Form("%s_cosmic_z_reco%d;z (cm)",h_name[i],ireco), 200, 1600, 2000);
      h_cosmic_r2[i][ireco]     = new TH1D(Form("%s_cosmic_r2_reco%d",h_name[i],ireco),  Form("%s_cosmic_r2_reco%d; R^{2} (cm^{2})",h_name[i],ireco), 50, 0., 2856100.);
      h_cosmic_ddir[i][ireco] = new TH1D(Form("%s_cosmic_ddir_reco%d",h_name[i],ireco), Form("%s_cosmic_ddir_reco%d;Degrees",h_name[i],ireco),  50, 0., 30.);
      
      h_eff[i][ireco] = new TH1D(Form("%s_eff_reco%d",h_name[i],ireco), Form("%s_eff_reco%d",h_name[i],ireco), 9, 0, 9);

      
      h_cosmic_evis[i][ireco]     	  ->Sumw2();
      h_cosmic_nhit_over_time[i][ireco] ->Sumw2();
      h_cosmic_nmue[i][ireco]     	  ->Sumw2();
      h_cosmic_z[i][ireco]    	  ->Sumw2();
      h_cosmic_r2[i][ireco]   	  ->Sumw2();
      h_cosmic_ddir[i][ireco] 	  ->Sumw2();

      h_cosmic_lnLrat1R2R_vs_mom[i][ireco]     = new TH2D(Form("%s_cosmic_lnLrat1R2R_vs_mom_reco%d",h_name[i],ireco),  Form("%s_cosmic_lnLrat1R2R_vs_mom_reco%d;p_{#mu} (MeV/c);ln(L_{2R}/L_{1R})",h_name[i],ireco), 250, 0, 1000, 500, -1000, 1000);
      h_cosmic_lnLrat1R2R_vs_mom[i][ireco]     ->Sumw2();

      h_cosmic_lnLrat1R2R_vs_pid[i][ireco]     = new TH2D(Form("%s_cosmic_lnLrat1R2R_vs_pid_reco%d",h_name[i],ireco),  Form("%s_cosmic_lnLrat1R2R_vs_pid_reco%d;ln(L_{2R}/L_{1R});PID",h_name[i],ireco), 200, -1000, 1000, 4, 0, 4);
      h_cosmic_lnLrat1R2R_vs_pid[i][ireco]     ->Sumw2();
	
      h_cosmic_lnLrat1R2R[i][ireco]     = new TH1D(Form("%s_cosmic_lnLrat1R2R_reco%d",h_name[i],ireco),  Form("%s_cosmic_lnLrat1R2R_reco%d;ln(L_{2R}/L_{1R})",h_name[i],ireco), 200, -1000, 1000);
      h_cosmic_lnLrat1R2R[i][ireco]     ->Sumw2();

      for (int imom=0; imom<nMoms; imom++) {
	
	int ise=0;
	double momLowEdge = (momLow[ise]+imom*momStep[ise]);
	double momHighEdge = (momLow[ise]+(imom+1)*momStep[ise]);
	
	h_cosmic_lnLrat1R2R_momsep[i][ireco][imom] = new TH1D(Form("%s_cosmic_lnLrat1R2R_reco%d_momsep%03d_%03d",h_name[i],ireco,(int)momLowEdge,(int)momHighEdge), Form("%s_cosmic_lnLrat1R2R_reco%d_momsep%03d_%03d; ln(L_{2R}/L_{1R})",h_name[i],ireco,(int)momLowEdge,(int)momHighEdge), 200, -1000, 1000);
	h_cosmic_lnLrat1R2R_momsep[i][ireco][imom]->Sumw2();
	
      }
      
      for (int icut=0; icut<nCuts; icut++) {
	h_cosmic_lifetime[i][ireco][icut] = new TH1D(Form("%s_cosmic_lifetime_reco%d_cut%d",h_name[i],ireco,icut),Form("%s_cosmic_lifetime_reco%d_cut%d; Decay-e Lifetime (#mus)",h_name[i],ireco,icut), 100, 0, 25);
	h_cosmic_lifetime[i][ireco][icut] ->Sumw2();


	h_cosmic_range_vs_mom[i][ireco][icut]     = new TH2D(Form("%s_cosmic_range_vs_mom_reco%d_cut%d",h_name[i],ireco, icut), Form("%s_cosmic_range_vs_mom_reco%d_cut%d;Momentum (MeV/c);Range (cm)",h_name[i],ireco, icut), 250, 0, 1000, 100, 0, 700);
	h_cosmic_range_vs_mom[i][ireco][icut]   ->Sumw2();

	h_cosmic_range_over_mom[i][ireco][icut]     = new TH1D(Form("%s_cosmic_range_over_mom_reco%d_cut%d",h_name[i],ireco, icut), Form("%s_cosmic_range_over_mom_reco%d_cut%d; Range/Momentum (cm / MeV/c)",h_name[i],ireco, icut), 100, 0, 1);
	h_cosmic_range_over_mom[i][ireco][icut] ->Sumw2();
      }

    }
  }
  cout << " current data livetime is " << get_datalivetime() << endl;

}

void OfficialEventParser::Parse(int ireco){
	// This is called event by event, and filling histograms.

	// initialize flag event by event if need
	InitEvent();
	
	if (get_datatype() == fcData || get_datatype() == fcMC){

	  if ( get_datatype() == fcData && ! __FC_RUN_SELECTION__ ) {
	  } else {


	    if (isCosmic(ireco)) {
	      FillCosmic(ireco);
	    }

	  }
	}
}


bool OfficialEventParser::isCosmic(int ireco){

  bool isTrue = true;
  
  return isTrue;

}


void OfficialEventParser::FillCosmic(int ireco){

  int imc=-1;
  if (get_datatype() == fcMC) imc = 1; // mc
  else if (get_datatype() == fcData) imc = 0;

  int ise_muon = 0;    
  int ise_dcye = 1;

  float etime = get_dcye_time(ireco, ise_dcye);
  h_cosmic_lifetime[imc][ireco][0]->Fill(etime, weight());

  // Cut: Apply decay-e time cut since Josh's data was preselected
  if (etime>1.2) {

    double evis = get_evis(ireco);
    h_cosmic_evis[imc][ireco]->Fill(evis, weight());

    // Cut: evis for sufficient light in event (should probably use 
    // potot instead, since no directly equivalent variable for evis in FQ)
    if (evis > 100) {
    
      TVector3 vertex(get_vertex(ireco, ise_muon));
      float r2 = vertex[0]*vertex[0]+vertex[1]*vertex[1];
    
      float vertex_xyz[3];
      vertex.GetXYZ(vertex_xyz);
      h_cosmic_wall[imc][ise_muon][ireco]->Fill(get_wall(vertex_xyz),weight());
    
      // Cut: Vertex in fiducial R^2
      h_cosmic_r2[imc][ireco]->Fill(r2 , weight());
      if ( hasRvtxInFid(ireco) ) {

	float muon_mom; // = get_1Rmom(ireco, ise_dcye);
	if (ireco==apfit) muon_mom = E->amomm[0];
	if (ireco==fitqun) muon_mom = E->fq1rmom[ise_muon][imu];
	int muon_mom_idx = min((int)((muon_mom-momLow[ise_muon])/momStep[ise_muon]),nMoms-1);
	      
	// PID Likelihood
	float lnLrat1R_muon = getLnLrat1R(ireco, ise_muon);
	h_cosmic_lnLrat1R_vs_mom[imc][ise_muon][ireco]->Fill(E->fq1rmom[ise_muon][ie],lnLrat1R_muon, weight());
	h_cosmic_lnLrat1R[imc][ise_muon][ireco]->Fill(lnLrat1R_muon, weight());
	
	h_cosmic_lnLrat1R_momsep[imc][ise_muon][ireco][muon_mom_idx]->Fill(lnLrat1R_muon, weight());

	// Cut: PID must be muon-like cut
	if ( (ireco==apfit && get1rType(ireco,ise_muon) == 3 ) ||
	     (ireco==fitqun && get1rType(ireco,ise_muon) == imu ) ) {
	
	  // Single ring likelihood
	  int fillPID;
	  float lnLrat1R2R_muon = getLnLrat1R2R(ireco,fillPID); 
	  h_cosmic_lnLrat1R2R_vs_mom[imc][ireco]->Fill(muon_mom,lnLrat1R2R_muon, weight()); 
	  h_cosmic_lnLrat1R2R[imc][ireco]->Fill(lnLrat1R2R_muon, weight()); 
	  h_cosmic_lnLrat1R2R_vs_pid[imc][ireco]->Fill(lnLrat1R2R_muon, fillPID, weight()); 
	  
	  h_cosmic_lnLrat1R2R_momsep[imc][ireco][muon_mom_idx]->Fill(lnLrat1R2R_muon, weight()); 
	  
	  // Cut: 1-ring
	  if (is1R(ireco)) {

	    h_cosmic_z[imc][ireco]->Fill(vertex[2] , weight()); 
	  
	    h_cosmic_pcflag[imc][ise_muon][ireco]->Fill(E->fq1rpcflg[ise_muon][imu], weight());
	    
	    h_cosmic_mom[imc][ise_muon][ireco][0]->Fill(muon_mom,weight());

	    TVector3 muon_dir(get_1Rdir(ireco, ise_muon));
	    h_cosmic_zenith[imc][ise_muon][ireco]->Fill(-muon_dir[2],weight());
	    
	    int nmue = muedcy(ireco);
	    h_cosmic_nmue[imc][ireco]->Fill(nmue , weight());
	    
	    // Cut: 1 decay-e
	    if (nmue==1) {

	      float dcye_mom = -1; // = get_1Rmom(ireco, ise_dcye);
	      //if (ireco==apfit) dcye_mom = ; // How to get decay-e momentum?
	      if (ireco==fitqun) dcye_mom = E->fq1rmom[ise_dcye][ie];
	      int dcye_mom_idx = min((int)((dcye_mom-momLow[ise_dcye])/momStep[ise_dcye]),nMoms-1);

	      // Look at muon momentum again after removing in-gates
	      h_cosmic_mom[imc][ise_muon][ireco][1]->Fill(muon_mom,weight());
	      
	      // PID for sub-event is electron
	      float lnLrat1R_dcye = getLnLrat1R(ireco, ise_dcye);
	      h_cosmic_lnLrat1R_vs_mom[imc][ise_dcye][ireco]->Fill(E->fq1rmom[ise_dcye][ie],lnLrat1R_dcye, weight());
	      h_cosmic_lnLrat1R[imc][ise_dcye][ireco]->Fill(lnLrat1R_dcye, weight());
	      
	      h_cosmic_lnLrat1R_momsep[imc][ise_dcye][ireco][dcye_mom_idx]->Fill(lnLrat1R_dcye, weight());
	      
	      // Cut: Out-of-gate (APFIT) and e-like (fiTQun)
	      if ( (ireco==apfit && E->egood[0]) ||
		   (ireco==fitqun && get1rType(ireco,ise_dcye) == ie ) ) 
		{
		  
		  h_cosmic_lifetime[imc][ireco][1]->Fill(etime, weight());
		  
		  TVector3 vertex_dcye(get_vertex(ireco,ise_dcye));
		  
		  float vertex_dcye_xyz[3];
		  vertex_dcye.GetXYZ(vertex_dcye_xyz);
		  float wall_dcye = get_wall(vertex_dcye_xyz);
		  h_cosmic_wall[imc][ise_dcye][ireco]->Fill(wall_dcye,weight());

		  TVector3 vertex_diff = vertex_dcye-vertex;
		  double range = vertex_diff.Mag();

		  h_cosmic_range_vs_mom[imc][ireco][0]->Fill(muon_mom, range);
		  h_cosmic_range_over_mom[imc][ireco][0]->Fill(range/muon_mom);

		  // Cut: Decay-e vertex in FV
		  if (wall_dcye>200) {
		  
		    // PC flag of decay-e
		    int pcflg_dcye = E->fq1rpcflg[ise_dcye][ie];
		    h_cosmic_pcflag[imc][ise_dcye][ireco]->Fill(pcflg_dcye, weight());

		    // Cut: PC flag should have no effect for stopping muons, but apply anyways to check
		    if ( ireco==apfit ||
			 (ireco==fitqun && !pcflg_dcye) ) {		
		      
		      h_cosmic_range_vs_mom[imc][ireco][1]->Fill(muon_mom, range);
		      h_cosmic_range_over_mom[imc][ireco][1]->Fill(range/muon_mom);

		      h_cosmic_mom[imc][ise_dcye][ireco][1]->Fill(dcye_mom, weight());
	      	      
		      TVector3 dcye_dir(get_1Rdir(ireco, ise_dcye));
		      h_cosmic_zenith[imc][ise_dcye][ireco]->Fill(-dcye_dir[2],weight());

		      double ddir = vertex_diff.Angle(muon_dir)*180/TMath::Pi();
		      h_cosmic_ddir[imc][ireco]->Fill(ddir, weight());
		  
		      // MC Resolution
		      if (imc) {
			for (int ise=0; ise<nSE; ise++) {
			  float mom_rec, mom_tru;
			  TVector3 v_dir_rec, v_vtx_rec;
			  TVector3 v_dir_tru, v_vtx_tru;
			  float mom_res, tvtx_res, lvtx_res, dir_res;
		      
			  float time_scnd;
			  int found_true_particle = 0;
			  
			  // Muon information
			  if (ise==ise_muon) {
			    mom_rec = muon_mom;
			    mom_tru = E->Abspvc[0];
		  
			    v_dir_rec = muon_dir;
			    for (int ix=0; ix<3; ix++) 
			      v_dir_tru[ix] = E->Pvc[0][ix]/E->Abspvc[0];
		  	
			    v_vtx_rec = vertex;
			    for (int ix=0; ix<3; ix++) 
			      v_vtx_tru[ix] = E->posv[ix];

			    found_true_particle = 1;
			  }
			  
			  // Electron information
			  else if (ise==ise_dcye) {
			    mom_rec = dcye_mom;
			    v_dir_rec = dcye_dir;
			    v_vtx_rec = vertex_dcye;
			    
			    // Find true dcy-e info in particle stack
			    int ip=0;
			    for (ip=0; ip<E->nscndprt; ip++) {
			      
			      for (int ix=0; ix<3; ix++) {
				v_dir_tru[ix] = E->pscnd[ip][ix];
				v_vtx_tru[ix] = E->vtxscnd[ip][ix];
			      }
			      mom_tru = v_dir_tru.Mag();
			      v_dir_tru *= 1/mom_tru;
			      
			      if ((abs(E->iprtscnd[ip])==11 || abs(E->iprtscnd[ip])==22) && 
				  E->lmecscnd[ip]==5) {
				found_true_particle++;
				break;
			      }
			    }
			    if (found_true_particle>1) {
			      cout << "Error expected 1 true decay-e but found " << found_true_particle << " " << E->tscnd[ip-1] << " " << E->iprtscnd[ip-1] << endl;
			      //exit (1);
			    }
			  }

			  if (found_true_particle>0) {
			    
			    mom_res = (mom_rec-mom_tru)/mom_tru*100;
			    h_cosmic_mom_res[imc][ise][ireco][0]->Fill(mom_res);
			    
			    dir_res = v_dir_rec.Angle(v_dir_tru)*180/TMath::Pi();
			    h_cosmic_dir_res[imc][ise][ireco][0]->Fill(dir_res);
			    
			    lvtx_res = (v_vtx_rec-v_vtx_tru).Dot(v_dir_tru);
			    h_cosmic_lvtx_res[imc][ise][ireco][0]->Fill(lvtx_res);
			    
			    tvtx_res = ((v_vtx_rec-v_vtx_tru).Cross(v_dir_tru)).Mag();
			    h_cosmic_tvtx_res[imc][ise][ireco][0]->Fill(tvtx_res);

			  }
			}
		      }
		  
		    }
		  }
		}
	    }
	  }
	}
      }
    }
  }


}


float OfficialEventParser::get_nhits_over_time(int ireco){

  float hits_time = 0;

  for (int icluster=0; icluster<E->cluster_ncand; icluster++) {

    if (E->cluster_goodflag[icluster]==1) {
      hits_time = E->cluster_nhits[icluster]/(E->cluster_tend[icluster]-E->cluster_tstart[icluster]);
      break;
    }
  }

  return hits_time;
  
}


float OfficialEventParser::get_dcye_time(int ireco, int ise){
  
  float time = -1;
  
  if (ireco==apfit) {
    if (ise==0) time = 0;
    else if (ise>0) time = E->etime[ise-1];
  }
  else if (ireco==fitqun) {
    if (ise>=0) {
      int i1RPID0 = get1rType(ireco, 0);
      int i1RPID1 = get1rType(ireco, ise);
      
      time  = (E->fq1rt0[ise][i1RPID1] - E->fq1rt0[0][i1RPID0])/1000;
    }
  }
  
  return time;
}


int OfficialEventParser::get1rType(int ireco, int ise){
  
  if (ireco==apfit) {
    return E->ip[0];
  }
  else if (ireco==fitqun) {

    double lnLrat1R=getLnLrat1R(ireco,ise);

    if (lnLrat1R > E->fq1rmom[ise][ie]*0.2) {//PID'd as electron
      return ie;
    }
    
    else {//muon
      return imu;
    }

  }
}

double OfficialEventParser::getLnLrat1R(int ireco, int ise){

  if (ireco==apfit) return 0;
  else if (ireco==fitqun) {
    return -E->fq1rnll[ise][ie]+E->fq1rnll[ise][imu];
  }

}



double OfficialEventParser::getMin2Rnll(int ireco, int &fillPID){
  
  double minNll = 9999999;
  int minPID1, minPID2;
  for (int iPID1=0; iPID1<2; iPID1++) 
    for (int iPID2=0; iPID2<2; iPID2++) 
      if (E->fq2rnll[iPID1][iPID2]<minNll) {
	minNll = E->fq2rnll[iPID1][iPID2];
	minPID1 = iPID1; minPID2 = iPID2;
      }
  
  fillPID = -1;
  if (minPID1==0 && minPID2==0) fillPID = 0;
  else if (minPID1==0 && minPID2==1) fillPID = 1;
  else if (minPID1==1 && minPID2==0) fillPID = 2;
  else if (minPID1==1 && minPID2==1) fillPID = 3;
  
  return minNll;

}

double OfficialEventParser::getLnLrat1R2R(int ireco, int &fillPID){

  if (ireco==apfit) return 0;
  else if (ireco==fitqun) {
    int ise = 0;
    int i1RPID = get1rType(ireco,ise);
    return -getMin2Rnll(ireco,fillPID) + E->fq1rnll[ise][i1RPID];
  }

}

bool OfficialEventParser::is1R(int ireco){
  // 1 ring event
  bool isTrue = false;
  if (ireco==apfit) {
    if (E->nring == 1) isTrue = true;
	  
  }
  else if (ireco==fitqun) {
	  
    int fillPID;
    double lnLrat1R2R=getLnLrat1R2R(ireco,fillPID); //likelihood ratio: best 2R vs. best 1R
    if (lnLrat1R2R<150.) isTrue = true;
	  	    
  }
	
  return isTrue;
}

// Warning: slightly different from E->wall
float OfficialEventParser::get_wall(float *vertex) {
      
  double r = vertex[0]*vertex[0] + vertex[1]*vertex[1]; 
  r = sqrt( r ) ;
  float wall = RINTK - r ; 
  double z1 = ZPINTK - vertex[2] ;
  double z2 = vertex[2] - ZMINTK ;
  if ( z1<wall ) wall = z1 ; 
  if ( z2<wall ) wall = z2 ; 
  return wall ; 
  
}


bool OfficialEventParser::hasRvtxInFid(int ireco){
  // vertex is in the Fid in r-direction :  r=sqrt(x*x + y*y)  is less than 1490.
  bool isTrue = false;

  float *vertex = get_vertex(ireco);

  double r = sqrt(vertex[0] * vertex[0] + vertex[1] * vertex[1]);
  if (r < RINTK-200) isTrue = true;

  return isTrue;
  // R (sqrt(pos[1]**2+pos(2)**2)<1490.)
}

bool OfficialEventParser::hasZvtxInFid(int ireco){
  // vertex is in the Fid in z-direction : z=abs(pos(3))<1610.
  bool isTrue = false;

  float *vertex = get_vertex(ireco);
	
  double z = TMath::Abs(vertex[2]);
  if (z < ZPINTK-200) isTrue = true;
	
  return isTrue;
}

float OfficialEventParser::cut_nhitac(){
  if (get_version() == 1) return 10.;      // nhitac<=9 for sk1
  else if (get_version() == 2) return 16.; // nhitac<=15 for sk2
  else if (get_version() == 3) return 16.; // nhitac<=15 for sk3
  else if (get_version() == 4) return 16.; // nhitac<=15 for sk4
  else { cerr << " Need to set nhitac cut in cut_nhitac()" << endl; exit(-1);}

}

float OfficialEventParser::cut_potot(){
  if (get_version() == 1) return 3000.;      //  for sk1
  else if (get_version() == 2) return 1500.; //  for sk2
  else if (get_version() == 3) return 1500.; //  for sk3
  else if (get_version() == 4) return 1500.; //  for sk4
  else { cerr << " Need to set nhitac cut in cut_potot()" << endl; exit(-1);}

}

int OfficialEventParser::cut_ehit( int i ){
  //cut value of ehit(i) in muedcy 
  // i==1 or 2 means just first appeared ehit(i) or second in code

  if ( i==1) { 
    if (get_version() == 1) return 60;      //  for sk1
    else if (get_version() == 2) return 30; //  for sk2
    else if (get_version() == 3) return 60; //  for sk3
    else if (get_version() == 4) return 60; //  for sk4
    else { cerr << " Need to set ehit cut in cut_ehit()" << endl; exit(-1);}
  } else if (i==2) {
    if (get_version() == 1) return 40;      //  for sk1
    else if (get_version() == 2) return 20; //  for sk2
    else if (get_version() == 3) return 40; //  for sk3
    else if (get_version() == 4) return 40; //  for sk4
    else { cerr << " Need to set ehit cut in cut_ehit()" << endl; exit(-1);}
  }
}

int OfficialEventParser::muedcy(int ireco){
  int muedcy;
  int nmuemax;

  /* temporary because SK-IV is not ready  */
  if ( get_version() == 4) {
    if (ireco==apfit)
      return E->nmue;
    else if (ireco==fitqun)
      return E->nse - 1;
  }
	
  if ( E->nmue <= 0 ) return 0;

  nmuemax = E->nmue;
  muedcy = 0;
  if ( E->nmue >= 10 ) nmuemax = 10;
	
  for ( int i=0; i< nmuemax ; i++){
    if ( E->evis < 1330.  && E->etime[i] < 0.1 ) continue;
    if ( E->evis >= 1330. && E->etime[i] < 1.2 ) continue;
    if ( E->etime[i] > 0.8 && E->etime[i] < 1.2 ) continue;
    if ( E->etype[i] == 1 && E->ehit[i] >= cut_ehit(1) && E->egood[i] > 0.5  )  { muedcy++; continue; }
    if ( (E->etype[i] == 2 || E->etype[i] == 3 || E->etype[i]== 4 ) && E->ehit[i]>=cut_ehit(2) ) { muedcy++; continue; }
  }

  return muedcy;

}

int OfficialEventParser::get_num_ring(int ireco){
  
  int rings = 0;
  if (ireco==apfit) rings = E->nring;
  else if (ireco==fitqun) {
    if (is1R(ireco)) rings = 1;
    else rings = 2;
  }
  return rings;
}

float OfficialEventParser::get_evis(int ireco){
  if (ireco==apfit) return E->evis;
  else if (ireco==fitqun) return E->fq1rmom[0][ie];
}

float *OfficialEventParser::get_vertex(int ireco, int ise){

  if (ireco==apfit) {
    if (ise==0) return E->pos;
    else return E->epos[ise-1];
  }
  else if (ireco==fitqun) {
    int i1rPID = get1rType(ireco,ise);
    return E->fq1rvtx[ise][i1rPID];
  }
  
}

float *OfficialEventParser::get_1Rdir(int ireco, int ise){

  if (ireco==apfit) {
    if (ise==0) return E->dir[0];
    else return E->edir[ise-1];    
  }
  else if (ireco==fitqun) {
    int i1rPID = get1rType(ireco,ise);
    return E->fq1rdir[ise][i1rPID];
  }

}

float OfficialEventParser::get_1Rmom(int ireco, int ise){

  float mom = -1;

  if (ireco==apfit) {
    if (ise==0) {
      if (E->ip[0]==2) return E->amome[0];
      else if (E->ip[0]==3) return E->amomm[0];
    }
    
  }
  else if (ireco==fitqun) {

    int i1rPID = get1rType(ireco, ise);
    return  E->fq1rmom[ise][i1rPID];	

  }

  return mom;

}

double OfficialEventParser::get_fromwall( int nth_ring ){
  double from_wall;

  set_from_wall_flag(true);
  // inverse the ring direction, calculate "towall" instead
  from_wall = get_towall( nth_ring);
  set_from_wall_flag(false);// reset the flag

  return from_wall;
}

double OfficialEventParser::get_towall( int nth_ring){
  /*
    Calculate towall for nth ring in the event.
    Original :: $ATMPD_ROOT/src/analysis/loe/lelib/towall.F
  */
  double x1, y1, z1;
  double x, y, z;
  double dx, dy, dz;
  double a0, b0, c0;
  double cz;
  double towall;
  double  rcir;
  double  hpos[3]; // 
  bool   might_go_to_side;
  double hanbe;
  double cz1;

  // get_towall_dir() can obtain where the ring direction goes.
  //  0 : out of tank,  1 : photon goes to barrel,   2 : photon goes to top,  3 : photon goes to bottom

  // initialize
  set_towall_dir(out_of_tank); // vertex is out of the tank
  might_go_to_side = false;
  cz = 0.;
  cz1 = 0.;
  towall = -999.; // initialize
  rcir = 0.;
  hanbe = 0.;

  // get vertex
  x1=E->pos[0];	
  y1=E->pos[1];	
  z1=E->pos[2];	

  // get direction of the nth ring
  if ( get_from_wall_flag()) { // calclulate "from wall". So, inverse the ring direction,
    dx = -E->dir[nth_ring][0];
    dy = -E->dir[nth_ring][1];
    dz = -E->dir[nth_ring][2];
  } else { // default
    dx = E->dir[nth_ring][0];
    dy = E->dir[nth_ring][1];
    dz = E->dir[nth_ring][2];
  }

  //      IF(DZ) 12,20,11   --> negative : 12,   zero : 20,  plus : 11
  // the ring is down-going
  if (dz < 0.) {
    if (z1 < (double)ZMINTK) { // z vertex is at underneath of ID barrel.
      set_towall_dir(out_of_tank);
      return towall;
    } else {
      cz = (ZMINTK - z1) / (double)dz;
      x  = cz * dx + x1;
      y  = cz * dy + y1;
      z = (double)ZMINTK;
      rcir = TMath::Sqrt(x*x+y*y);

      if (rcir > (double)RINTK) {
	might_go_to_side = true;	
      } else {
	// photon goes to lower plate (bottom)
	set_towall_dir(to_bottom);
	towall = calc_towall(x,y,z);
      }
    }
  }  

  // the ring is upper-going
  if (dz > 0.){
    if ( z1 >= (double)ZPINTK){ // z vertex is higher than the top of the tank
      set_towall_dir(out_of_tank);
      return towall;
    } else {
      cz = ((double)ZPINTK - z1) / (double)dz;
      x  = cz * dx + x1;
      y  = cz * dy + y1;
      z  = (double)ZPINTK;
      rcir = TMath::Sqrt(x*x+y*y);

      if (rcir >= (double)RINTK){
	might_go_to_side = true;
      } else {
	// photon goes to upper plate (top)
	set_towall_dir(to_top);
	towall = calc_towall(x,y,z);
      }
    }
  }

  // the ring goes horizontally
  if (dz == 0. || might_go_to_side == true ) {
    double _rintk = (double)RINTK;
    a0 = dx*dx + dy*dy;
    b0 = dx*x1 + dy*y1;
    c0 = x1*x1 + y1*y1 - _rintk*_rintk; // RINTK*RINTK;
		
    hanbe = b0*b0 - c0*a0;

    if (hanbe < 0.0) {
      set_towall_dir(out_of_tank);
      return towall;
    }

    cz1 = (double)( -b0 + TMath::Sqrt(hanbe)) / (double)a0;

    if ( cz1 <= 0.0 ) {
      set_towall_dir(out_of_tank);
      return towall;
    }
		
    x = cz1 * dx + x1;
    y = cz1 * dy + y1;
    z = cz1 * dz + z1;

    if ( z >= ZPINTK || z <= ZMINTK) {
      set_towall_dir(out_of_tank);
      return towall;
    }

    set_towall_dir(to_side);
    towall = calc_towall(x,y,z);
  } 

  return towall;
}

double OfficialEventParser::calc_towall(double x, double y, double z, int ireco){
  double towall;
  double hpos[3];

  hpos[0] = x;
  hpos[1] = y;
  hpos[2] = z;

  ireco = apfit; // Currently used with multi-ring so APFIT only for now
  float *vertex = get_vertex(ireco);

  towall = TMath::Sqrt( (vertex[0] - hpos[0]) * (vertex[0] - hpos[0])
			+ (vertex[1] - hpos[1]) * (vertex[1] - hpos[1])
			+ (vertex[2] - hpos[2]) * (vertex[2] - hpos[2]));
  return towall;
}


double OfficialEventParser::weight(){
  // Flux and oscillation weght fot MC
  double w = 1;
  return w;
}

void OfficialEventParser::Print(){
}

void OfficialEventParser::Terminate(){
  
  // Efficiency
  int ise_muon = 0;
  int ise_dcye = 1;
  {
    for (int imc=0; imc<nSample; imc++) {
      for (int ireco=0; ireco<nRecos; ireco++) {
	h_eff[imc][ireco]->SetBinContent(1, h_cosmic_evis[imc][ireco]->GetEntries());
	h_eff[imc][ireco]->SetBinContent(2, h_cosmic_r2[imc][ireco]->GetEntries());
	h_eff[imc][ireco]->SetBinContent(3, h_cosmic_lnLrat1R_vs_mom[imc][ise_muon][ireco]->GetEntries());
	h_eff[imc][ireco]->SetBinContent(4, h_cosmic_lnLrat1R2R[imc][ireco]->GetEntries());
	h_eff[imc][ireco]->SetBinContent(5, h_cosmic_nmue[imc][ireco]->GetEntries());
	h_eff[imc][ireco]->SetBinContent(6, h_cosmic_lnLrat1R_vs_mom[imc][ise_dcye][ireco]->GetEntries());
	h_eff[imc][ireco]->SetBinContent(7, h_cosmic_wall[imc][ise_dcye][ireco]->GetEntries());
	h_eff[imc][ireco]->SetBinContent(8, h_cosmic_pcflag[imc][ise_dcye][ireco]->GetEntries());
	h_eff[imc][ireco]->SetBinContent(9, h_cosmic_ddir[imc][ireco]->GetEntries());
	      
	h_eff[imc][ireco]->Sumw2();
	h_eff[imc][ireco]->Scale(1/h_eff[imc][ireco]->GetBinContent(1));
      }
    }
  }

  //float scale = get_datalivetime();
  
  int i=1;
  for (int ireco=0; ireco<nRecos; ireco++) {
    
    float scale = h_cosmic_evis[0][ireco]->GetEntries()/h_cosmic_evis[1][ireco]->GetEntries();
    float scale_nmue =  h_cosmic_lnLrat1R_vs_mom[0][ise_dcye][ireco]->GetEntries()/h_cosmic_lnLrat1R_vs_mom[1][ise_dcye][ireco]->GetEntries(); 
    
    for (int ise=0; ise<2; ise++) {
      
      float scale_tmp = scale;
      if (ise>=1) scale_tmp = scale_nmue;

      h_cosmic_zenith[i][ise][ireco]         ->Scale(scale_tmp);
      h_cosmic_wall[i][ise][ireco]           ->Scale(scale_tmp);
      h_cosmic_lnLrat1R_vs_mom[i][ise][ireco]->Scale(scale_tmp);
      h_cosmic_lnLrat1R[i][ise][ireco]->Scale(scale_tmp);
      h_cosmic_pcflag[i][ise][ireco]         ->Scale(scale_tmp);

      for (int imom=0; imom<nMoms; imom++) 
	h_cosmic_lnLrat1R_momsep[i][ise][ireco][imom]->Scale(scale_tmp); 

      for (int icut=0; icut<nCuts; icut++) {
	
	float scale_tmp = scale;
	if (icut>=1) scale_tmp = scale_nmue;

	h_cosmic_mom[i][ise][ireco][icut]    ->Scale(scale_tmp);
	h_cosmic_mom_res[i][ise][ireco][icut]->Scale(scale_tmp); 
	h_cosmic_tvtx_res[i][ise][ireco][icut]->Scale(scale_tmp);
	h_cosmic_lvtx_res[i][ise][ireco][icut]->Scale(scale_tmp);
	h_cosmic_dir_res[i][ise][ireco][icut]->Scale(scale_tmp); 
	
      }
    }

    for (int imom=0; imom<nMoms; imom++) 
      h_cosmic_lnLrat1R2R_momsep[i][ireco][imom]->Scale(scale);

    h_cosmic_lnLrat1R2R[i][ireco]     ->Scale(scale);
    h_cosmic_lnLrat1R2R_vs_mom[i][ireco]     ->Scale(scale);
    h_cosmic_lnLrat1R2R_vs_pid[i][ireco]     ->Scale(scale);
    h_cosmic_evis[i][ireco]     	  ->Scale(scale);
    h_cosmic_nhit_over_time[i][ireco] ->Scale(scale);
    h_cosmic_nmue[i][ireco]     	  ->Scale(scale);
    h_cosmic_z[i][ireco]    	  ->Scale(scale);
    h_cosmic_r2[i][ireco]   	  ->Scale(scale);

    h_cosmic_ddir[i][ireco] 	  ->Scale(scale_nmue);

    for (int icut=0; icut<nCuts; icut++) {

      float scale_tmp = scale;
      if (icut>=1) scale_tmp = scale_nmue;

      h_cosmic_lifetime[i][ireco][icut] 	  ->Scale(scale_tmp);
      h_cosmic_range_vs_mom[i][ireco][icut]   ->Scale(scale_tmp);
      h_cosmic_range_over_mom[i][ireco][icut] ->Scale(scale_tmp);
    }
  }

  for (int i=0; i<nSample; i++) {
    for (int ireco=0; ireco<nRecos; ireco++) {
      efile[i][ireco]->Write();
      efile[i][ireco]->Close();
    }
  }


}
