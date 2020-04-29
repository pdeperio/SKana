#include "global.h"
#include <string>


	
double global::fidvolume = 22.458;

//
//   SK1
//
double global::sk1_data_livetime = 1489.2;
double global::sk1_fcmc_livetime = 182625.; // 500yr MC
double global::sk1_pcmc_livetime = 182625.;
std::string global::sk1fcdatafilename="/home/atmpd/ntuple/apr08_sk1/fc_dst_root/sk1.*.root";
std::string global::sk1fcmcfilename  ="/home/atmpd/ntuple/apr08_sk1/fc_mc_root/sk1.*.root";
std::string global::sk1pcdatafilename="/home/atmpd/ntuple/apr08_sk1/pc_dst_root/split.*.fit.root";
std::string global::sk1pcmcfilename="/home/atmpd/ntuple/apr08_sk1/pc_mc_root/patrd6.run000*.001.root";

//
//   SK2
//
double global::sk2_data_livetime = 798.6;
double global::sk2_fcmc_livetime = 182625.; // 500yr MC
double global::sk2_pcmc_livetime = 182625.;
std::string global::sk2fcdatafilename="/home/atmpd/ntuple/apr08_sk2/fc_dst_root/sk2.*.root";
std::string global::sk2fcmcfilename  ="/home/atmpd/ntuple/apr08_sk2/fc_mc_root/sk2.*.root";
std::string global::sk2pcdatafilename="/home/atmpd/ntuple/apr08_sk2/pc_dst_root/split.*.fit.root";
std::string global::sk2pcmcfilename="/home/atmpd/ntuple/apr08_sk2/pc_mc_root/patrd6.run000*.001.root";

//
//   SK3 
//
double global::sk3_data_livetime = 518.084;
double global::sk3_fcmc_livetime = 500.*365.25; 
double global::sk3_pcmc_livetime = 500.*365.25;
std::string global::sk3fcdatafilename="/home/atmpd/ntuple/mar09_sk3/fc_dst_root/*.root";
std::string global::sk3fcmcfilename  ="/home/atmpd/ntuple/mar09_sk3/fc_mc_root/*.root";
std::string global::sk3pcdatafilename="/home/atmpd/ntuple/mar09_sk3/pc_dst_root/pcdst.sk3.09b.root";
std::string global::sk3pcmcfilename  ="/home/atmpd/ntuple/mar09_sk3/pc_mc_root/pcmc.sk3.09b.*.root";



//
//   SK4 realtime 
//
std::string global::sk4fcdatafilename="/net/sukond1/export/data/pdeperio/cosmics/escale/lowmu11d/data/stopmu.loosered.fit.newtool_fitqun.root";

//
//   SK4 Jul10 version
//
double global::sk4_data_livetime = 1097.*(19699./26331.); // = 820.698
double global::sk4_fcmc_livetime = 11.3045*365.25; 
double global::sk4_pcdata_livetime = 448.6; 
double global::sk4_pcmc_livetime = 100.*365.25; //  Oct 21, 2009
//std::string global::sk4fcdatafilename="/home/atmpd/ntuple/jul10_sk4/fc_dst_root/*.root" ; 
std::string global::sk4fcmcfilename  ="/net/sukond1/export/data/pdeperio/cosmics/escale/lowmu11d/mc/stopmu.mc.fit.newtool_fitqun.root";
std::string global::sk4pcdatafilename="";
std::string global::sk4pcmcfilename="";


//
//   SK4 Apr10 version
//
//double global::sk4_data_livetime = 448.6;
//double global::sk4_fcmc_livetime = 100.*365.25; 
//double global::sk4_pcdata_livetime = 448.6; 
//double global::sk4_pcmc_livetime = 100.*365.25; //  Oct 21, 2009
//std::string global::sk4fcdatafilename="/home/atmpd/ntuple/apr10_sk4/fc_dst_root/*.root" ; 
//std::string global::sk4fcmcfilename  ="/home/atmpd/ntuple/apr10_sk4/fc_mc_root/*.root";
//std::string global::sk4pcdatafilename="/home/atmpd/ntuple/apr10_sk4/pc_dst_root/pcdst.sk4.10a.run61790-66624.root";
//std::string global::sk4pcmcfilename="/home/atmpd/ntuple/apr10_sk4/pc_mc_root/pcmc.sk4.10a.*.root";


	
//
//   SK4 Old
//
//std::string global::sk4fcdatafilename="/home/atmpd/ntuple/oct09_sk4/fc_dst_root/*.root";  // official area .. 289.1 days
// this has updated decay electron information.  Jul 30, 2009
//std::string global::sk4fcmcfilename  ="/disk/atmpd4/fc/2009JulySK4/fcmcred_ntuple/sk4-*.red.root";
// TEST for newNEUT std::string global::sk4fcmcfilename  ="/disk/atmpd4/sk4_dst/oct09_NewNeut/fc_mc/ntuple/sk4-*-red.root";
//std::string global::sk4fcmcfilename  ="/home/atmpd/ntuple/oct09_sk4/fc_mc_root/sk4-*-red.root";
// new NEUT test ... double global::sk4_fcmc_livetime = 19.*365.25; 
// and also  DrawHistograms::Init livetime ratio is temporary adjust to SK3PC. 
// Must fix it if SK4 PC MC is ready.                  Naho 
// fc and PC data livetime is assumed to be the same, but not for SK4, 
//std::string global::sk4pcdatafilename="/home/atmpd/ntuple/mar09_sk4/pc_dst_root/pcdst.sk4.09d.run61790-64310.root";
//std::string global::sk4pcdatafilename="/home/atmpd/ntuple/oct09_sk4/pc_dst_root/pcdst.sk4.09e.root";
//std::string global::sk4pcdatafilename=sk3pcdatafilename;
//std::string global::sk4pcmcfilename=sk3pcmcfilename;
//std::string global::sk4pcmcfilename="/home/atmpd/ntuple/oct09_sk4/pc_mc_root/pcmc.sk4.09e.*.root";



	std::string global::flux_version_name="Honda 06";


std::string globalname::datatypename[26] = {
                         "fcData",
                         "pcData",
                         "fcMC",
                         "pcMC",
                         "upmuMC",

                         "SK1fcData",
                         "SK1pcData",
                         "SK1fcMC",
                         "SK1pcMC",
                         "SK1upmuMC",

                         "SK2fcData",
                         "SK2pcData",
                         "SK2fcMC",
                         "SK2pcMC",
                         "SK2upmuMC",

                         "SK3fcData",
                         "SK3pcData",
                         "SK3fcMC",
                         "SK3pcMC",
                         "SK3upmuMC",

                         "SK4fcData",
                         "SK4pcData",
                         "SK4fcMC",
                         "SK4pcMC",
                         "SK4upmuMC",
};

std::string pi0massfitparametername::st[6]={"lin:par0", "gaus:constant", "gaus:mean", "gaus:sigma", "exp:constant","exp:slope"};

std::string globalname::eventtypename[17]={
	"Sub-GeV e 0-dcy",
	"Sub-GeV e 1-dcy",
	"Sub-GeV pi0 1-R",
	"Sub-GeV mu 0-dcy",
	"Sub-GeV mu 1-dcy",
	"Sub-GeV mu 2-dcy",
	"Sub-GeV pi0 M-R",
	"Multi-GeV e",
	"Multi-GeV mu",
	"Multi-Ring e",
	"Multi-Ring mu",
	"PC stop",
	"PC through",
	"UP stop mu",
	"Non-showering mu",
	"Showering mu"
};
