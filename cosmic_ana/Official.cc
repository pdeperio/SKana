#include <iostream>
#include <cstdlib>
#include "TTree.h"
#include <fstream>
#include "Event.h"
#include "TChain.h"
#include "OfficialEventParser.h"
#include "global.h"
#include <vector>
#include <string>

using namespace std;
using namespace global;

int main(int argc, char *argv[]){
  int TotalEvents;
  int NEntries=0;
  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //      samples you want to loop over
  bool fcdata = true;
  bool fcmc   = true;
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // create an empty vector of strings
  vector<string> args;
  int i;
  for (i=1;i<argc;i++) 
    args.push_back(argv[i]);
  
  cout << "argc=" << argc << endl;
  
  int nrgSepType = 0; // 0: all, 1: sub-GeV, 2: multi-GeV
  if ( argc >1 ) nrgSepType = atoi(argv[1]);

  int maxE = 999999; 
  if ( argc >2 ) maxE = atoi(argv[2]);

  int fileNum = -1;
  if ( argc >3 ) fileNum = atoi(argv[3]);

  string inputfilename[nMCs][nRecos];
  
  inputfilename[0][fitqun] =  sk4fqdatafilename;
  inputfilename[0][apfit]  =  sk4stmudatafilename;
  inputfilename[1][fitqun] =  sk4fqmcfilename;
  inputfilename[1][apfit]  =  sk4stmumcfilename;
  
  TChain *_Chains[nMCs][nRecos];
  
  int version = 4; // SK4
  
  OfficialEventParser * Parser = new OfficialEventParser();
  Parser->set_version( version );
  Parser->set_nrgSepType( nrgSepType );
  Parser->set_maxE( maxE );
  Parser->set_fileNum( fileNum );
  //Parser->set_datalivetime( current_livetime );
  //Parser->InitHistograms(); // this needs livetime information
  
  for (int imc=0; imc<nMCs; imc++) {

    //if (!imc) continue;

    for (int ireco=0; ireco<nRecos; ireco++) {

      //if (ireco==fitqun) continue;

      Parser->InitHistograms(ireco, imc); // this needs livetime information

      _Chains[imc][ireco] = new TChain("h1");
      _Chains[imc][ireco]->Add(inputfilename[imc][ireco].c_str());
      
      (*Parser)(_Chains[imc][ireco],ireco);
      
      //if (ireco) _Chains[imc][0]->AddFriend(_Chains[imc][ireco]);
      //}
    
      Parser->set_datatype(imc);
    
    //for (int ireco=0; ireco<2; ireco++) {
      //int ireco=1; {
      
      NEntries = _Chains[imc][ireco]->GetEntries();
      cout << "imc = " << imc << ", ireco = " << ireco << ",  NEntries = " << NEntries << endl;
      
      TotalEvents = 0;
      
      int istart=0;
      int iend = NEntries;
      int nProcs = 8;
      if (fileNum>=0) {

	int eventsPerFile = NEntries/nProcs;

	istart = fileNum*eventsPerFile;
	iend = (fileNum+1)*eventsPerFile;

	if (iend>NEntries || fileNum>=nProcs-1) iend=NEntries;

      }

      for (int i=istart; i<iend; i++){
      //for (int i=istart; i<1500; i++){
	_Chains[imc][ireco]->GetEvent(i);
	
	if ( i % 100000 == 0 ) 
	{
	  cout << "File #" << fileNum << ", Event: " << i << endl;
	  //cout << " filename = " << _Chains[imc][ireco]->GetFile()->GetName()  << endl;
	}
	
	Parser->Parse(ireco);
	
	TotalEvents++;
      }
      
      cout <<  "File #" << fileNum << " total events : " << TotalEvents << endl;
      // Parser->Print();
      //for (int ireco=0; ireco<nRecos; ireco++) 
      delete _Chains[imc][ireco];
    }
  }
  
  //------------------------------------------------------
  
  Parser->Terminate();
}

