#ifndef _global_
#define _global_

#include <string>

namespace global
{

  const int nMCs = 2;

  const int nRecos = 2;
  enum reco_enum {apfit, fitqun};

  extern std::string sk4fqdatafilename;
  extern std::string sk4fqmcfilename;
  extern std::string sk4stmudatafilename;
  extern std::string sk4stmumcfilename;

  extern std::string flux_version_name;


  enum DataTypes{  
    fcData,
    fcMC                       
  };              

  namespace globalname{
    extern std::string datatypename[];
    extern std::string eventtypename[];
  }
}

#endif
