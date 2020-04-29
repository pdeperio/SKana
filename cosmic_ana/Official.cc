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
	bool pcdata = false;
	bool pcmc   = false;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	// create an empty vector of strings
	vector<string> args;
	int i;
	for (i=1;i<argc;i++) 
        	args.push_back(argv[i]);

	cout << "argc=" << argc << endl;

	if ( argc < 3  || argc > 4 ) {
		cerr << " format : official 1 noosc   OR   official 1 noosc test" << endl; exit(-1);
	}

	int version = atoi( argv[1] ); // 1 : SK1, 2 : SK2 , 3: SK3

	// oscillated or not
	bool isOscillated;
	if ( args[1] == "osc" || args[1] == "OSC" )    isOscillated =  true;
	else if ( args[1] == "noosc" || args[1] == "NOOSC") isOscillated =  false;
	else { cerr << "Please write osc or noosc " << endl; exit(-1);};


	bool isTestSmallSample = false; // default is fause
	if ( argc > 3 ) { 
        	if ( args[2] == "test" || args[2] == "TEST" )    isTestSmallSample =  true;
	}


        // livetime
//        cout << " arguments.size() = " << arguments.size() << endl;

        double input_livetime;
	double current_livetime;
        if ( argc>3 && isTestSmallSample == false ) {
                input_livetime = atof( args[2].c_str() );
		// this is not right way, but sort of sanity check.
		if ( input_livetime > 10000. || input_livetime < 0.) {
			cerr << " input livetime is strange, " << args[2].c_str() << endl; exit(-1);
		}

		// overwrite livetime to sk?_data_livetime
                if (version==1) sk1_data_livetime = input_livetime;
                else if (version==2) sk2_data_livetime = input_livetime;
                else if (version==3) sk3_data_livetime = input_livetime;
                else if (version==4) sk4_data_livetime = input_livetime;
        }

        if (version==1) current_livetime = sk1_data_livetime;
        else if (version==2) current_livetime = sk2_data_livetime;
        else if (version==3) current_livetime = sk3_data_livetime;
        else if (version==4) current_livetime = sk4_data_livetime;
        else {    cerr << " can not set current_livetime in  Official.cc " << endl; exit(-1);
        }


	cout << " current live time is " << current_livetime << endl;

// print out input information 
	cout << "sk version is " << version ;
	if (isOscillated) cout << ",  MC is Oscillated," ; 
	else cout << ",  MC is NOT Oscillated," ;

	if (isTestSmallSample) cout << ", and with small sample for test." << endl; 
	else cout << ", and with fill MC sample." << endl;

	string inputfcdataname;
	string inputfcmcname;
	string inputpcdataname;
	string inputpcmcname;

	if (version == 1){
		inputfcdataname    =  sk1fcdatafilename;
		inputfcmcname      =  sk1fcmcfilename;
		inputpcdataname    =  sk1pcdatafilename;
		inputpcmcname      =  sk1pcmcfilename;
	} else if (version == 2){
		inputfcdataname    =  sk2fcdatafilename;
		inputfcmcname      =  sk2fcmcfilename;
		inputpcdataname    =  sk2pcdatafilename;
		inputpcmcname      =  sk2pcmcfilename;
	} else if (version == 3){
		inputfcdataname    =  sk3fcdatafilename;
		inputfcmcname      =  sk3fcmcfilename;
		inputpcdataname    =  sk3pcdatafilename;
		inputpcmcname      =  sk3pcmcfilename;
	} else if (version == 4){
		inputfcdataname    =  sk4fcdatafilename;
		inputfcmcname      =  sk4fcmcfilename;
		inputpcdataname    =  sk4pcdatafilename;
		inputpcmcname      =  sk4pcmcfilename;
	} else if (version > 4 || version < 1){
		cout << " sk5 is not ready " << endl;
		exit(-1);
	}

	OfficialEventParser * Parser = new OfficialEventParser();
	Parser->set_version( version );
	Parser->set_datalivetime( current_livetime );
	Parser->InitHistograms(); // this needs livetime information

	// data  --------------------------------------------
	if (fcdata) {

		//Parser->set_fcmc(false);
		Parser->set_datatype(fcData);
        	TChain *_Chains = new TChain("h1");
        	_Chains->Add(inputfcdataname.c_str());

        	TotalEvents = 0;

        	NEntries = _Chains->GetEntries();
        	cerr << " NEntries = " << NEntries << endl;
        	(*Parser)(_Chains);

		// loop over data events
        	for (int i=0; i<NEntries; i++){
			if ( isTestSmallSample==true && TotalEvents > 100) continue;
                	_Chains->GetEvent(i);
                	
			for (int ireco=0; ireco<2; ireco++) Parser->Parse(ireco);

                	if ( i % 50000 == 0 ) {
                      		cerr << "Event: " << i << endl;
                      		cerr << " filename = " << _Chains->GetFile()->GetName()  << endl;
                	}

              		TotalEvents++;
        	}

        	cerr << " total data events : " << TotalEvents << endl;
//        	Parser->Print();
		delete _Chains;
	} 

	// Monte Carlo --------------------------------------------
	if (fcmc){

		// if you write any sort of counter in OfficialEventParser, be aware not to overwrite or add without initialize it.
		Parser->set_datatype(fcMC);
		Parser->set_osc(isOscillated);

        	TChain *_Chains_fcmc = new TChain("h1");
        	_Chains_fcmc->Add(inputfcmcname.c_str());

        	TotalEvents = 0;
		NEntries = 0;
        	NEntries = _Chains_fcmc->GetEntries();
        	cerr << " NEntries = " << NEntries << endl;
        	(*Parser)(_Chains_fcmc);

        	for (int i=0; i<NEntries; i++){
			if ( isTestSmallSample==true && TotalEvents > 100) continue;

                	_Chains_fcmc->GetEvent(i);
			for (int ireco=0; ireco<2; ireco++) Parser->Parse(ireco);

                	if ( i % 50000 == 0 ) {
                      		cerr << "Event: " << i << endl;
                      		cerr << " filename = " << _Chains_fcmc->GetFile()->GetName()  << endl;
                	}

              		TotalEvents++;
        	}

        	cerr << " total events : " << TotalEvents << endl;
//        	Parser->Print();
		delete _Chains_fcmc;
	}

	//------------------------------------------------------

	if (pcdata){
                Parser->set_datatype(pcData);
                TChain *_Chains_pcdata = new TChain("h1");
                _Chains_pcdata->Add(inputpcdataname.c_str());

                TotalEvents = 0;

                NEntries = _Chains_pcdata->GetEntries();
                cerr << " NEntries = " << NEntries << endl;
                (*Parser)(_Chains_pcdata);

                // loop over data events
                for (int i=0; i<NEntries; i++){
			if ( isTestSmallSample==true && TotalEvents > 100) continue;

                        _Chains_pcdata->GetEvent(i);
			for (int ireco=0; ireco<2; ireco++) Parser->Parse(ireco);

                        if ( i % 5000 == 0 ) {
                                cerr << "Event: " << i << endl;
                                cerr << " filename = " << _Chains_pcdata->GetFile()->GetName()  << endl;
                        }

                        TotalEvents++;
                }

                cerr << " total data events : " << TotalEvents << endl;
                delete _Chains_pcdata;

	}
	//------------------------------------------------------
	if (pcmc){
                Parser->set_datatype(pcMC);
                Parser->set_osc(isOscillated);

                TChain *_Chains_pcmc = new TChain("h1");
                _Chains_pcmc->Add(inputpcmcname.c_str());

                TotalEvents = 0;
                NEntries = 0; 
                NEntries = _Chains_pcmc->GetEntries();
                cerr << " NEntries = " << NEntries << endl;
                (*Parser)(_Chains_pcmc);

                for (int i=0; i<NEntries; i++){
			if ( isTestSmallSample==true && TotalEvents > 100) continue;
                        _Chains_pcmc->GetEvent(i);
			for (int ireco=0; ireco<2; ireco++) Parser->Parse(ireco);

                        if ( i % 50000 == 0 ) {
                                cerr << "Event: " << i << endl;
                                cerr << " filename = " << _Chains_pcmc->GetFile()->GetName()  << endl;
                        }

                        TotalEvents++;
                }

                cerr << " total events : " << TotalEvents << endl;
                delete _Chains_pcmc;


	}

        Parser->Terminate();


}

