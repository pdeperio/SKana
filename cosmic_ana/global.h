#ifndef _global_
#define _global_

#include <string>

namespace global
{

extern double sk1_data_livetime;
extern double sk1_fcmc_livetime;
extern double sk1_pcmc_livetime;

extern double sk2_data_livetime;
extern double sk2_fcmc_livetime;
extern double sk2_pcmc_livetime;

extern double sk3_data_livetime;
extern double sk3_fcmc_livetime;
extern double sk3_pcmc_livetime;

extern double sk4_data_livetime;
extern double sk4_fcmc_livetime;
extern double sk4_pcmc_livetime;

/*
extern double sk1_ktonyr;
extern double sk2_ktonyr;
extern double sk3_ktonyr;
extern double sk4_ktonyr;
*/

extern double sk4_pcdata_livetime;

extern double fidvolume;

/* moved to DrawHistograms.h and DrawHistograms.cc   Apr 30, 2009 Naho
extern std::string sk1title;
extern std::string sk2title;
extern std::string sk3title;
extern std::string sk4title;
*/

extern std::string sk1fcdatafilename;
extern std::string sk1fcmcfilename;
extern std::string sk2fcdatafilename;
extern std::string sk2fcmcfilename;
extern std::string sk3fcdatafilename;
extern std::string sk3fcmcfilename;
extern std::string sk4fcdatafilename;
extern std::string sk4fcmcfilename;

extern std::string sk1pcdatafilename;
extern std::string sk1pcmcfilename;
extern std::string sk2pcdatafilename;
extern std::string sk2pcmcfilename;
extern std::string sk3pcdatafilename;
extern std::string sk3pcmcfilename;
extern std::string sk4pcdatafilename;
extern std::string sk4pcmcfilename;

extern std::string flux_version_name;


        enum DataTypes{  
                         fcData,
                         pcData,
                         fcMC,
                         pcMC,
                         upmuMC,

                         SK1fcData,
                         SK1pcData,
                         SK1fcMC,
                         SK1pcMC,
                         SK1upmuMC,

                         SK2fcData,
                         SK2pcData,
                         SK2fcMC,
                         SK2pcMC,
                         SK2upmuMC,

                         SK3fcData,
                         SK3pcData,
                         SK3fcMC,
                         SK3pcMC,
                         SK3upmuMC,

                         SK4fcData,
                         SK4pcData,
                         SK4fcMC,
                         SK4pcMC,
                         SK4upmuMC,
        };              

	enum EventTypes{
        	SubGeVE0dcy,
        	SubGeVE1dcy,
        	SubGeVPi01R,
	        SubGeVMu0dcy,
        	SubGeVMu1dcy,
	        SubGeVMu2dcy,
	        SubGeVPi0MR,
	        MulGeVE,
	        MulGeVMu,
	        MultiRingE,
        	MultiRingMu,
        	PCstop,
        	PCthrough,
        	UPstopMu,
        	NonShoweringMu,
        	ShoweringMu
        };              

}

namespace globalname{
	extern std::string datatypename[];
	extern std::string eventtypename[];
}

namespace pi0massfitparametername{
	extern std::string st[];
}


#endif
