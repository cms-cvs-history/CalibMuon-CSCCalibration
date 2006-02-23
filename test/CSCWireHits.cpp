//This will read AFEB data and create a root ntuple. 

#include "IORawData/CSCCommissioning/src/FileReaderDDU.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDDUEventData.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigi.h"
#include "EventFilter/CSCRawToDigi/interface/CSCEventData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCAnodeData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCALCTHeader.h"
#include "EventFilter/CSCRawToDigi/interface/CSCCLCTData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDDUTrailer.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDMBHeader.h"

#include "CalibMuon/CSCCalibration/interface/condbc.h"
#include "CalibMuon/CSCCalibration/interface/cscmap.h" 
#include <iostream>
#include <fstream>
#include <vector>
#include "string"

//root specific .h files
#include </afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TROOT.h>
#include </afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TNtuple.h>
#include </afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TFile.h>
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TH1F.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TH2F.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TH1.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TCanvas.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TTree.h"

class TCalibEvt { public:
  Int_t event;  
  Int_t chamber;
  Int_t layer;
  Int_t wire;
  Int_t tbin;
};


int main(int argc, char **argv) {

  //for debugging purposes from Alex Tumanov's code
  //set to true if you wish to debug data
  CSCAnodeData::setDebug(false);
  CSCALCTHeader::setDebug(false);
  CSCCLCTData::setDebug(false); 
  CSCEventData::setDebug(false);  
  CSCTMBData::setDebug(false);
  CSCDDUEventData::setDebug(false);

  //root ntuple
  TCalibEvt calib_evt;
  TBranch *calibevt;
  TTree *calibtree;
  TFile *calibfile;
    
  calibfile = new TFile("calibwire.root","RECREATE","Calibration Ntuple");
  calibtree = new TTree("Calibration","AFEB");
  calibevt = calibtree->Branch("EVENT", &calib_evt, "event/I:chamber/I:layer/I:wire/I:tbin/I");
  TH2F *h1 = new TH2F("h1","wire vs timebin",10,-1,8,10,0,80);

  //data file: nr.of events,chambers read
  int maxEvents = 500000;
  int misMatch = 0;
  std::string datafile=argv[1];
  std::string chamber_id; 
  int dmbID=-999,crateID=-999,chamber_num,sector,ret_code;
  int reportedChambers =0;
   
  //needed for database and mapping
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();

  if (argv[2]) maxEvents = (int) atof(argv[2]);

  FileReaderDDU ddu;
  ddu.open(datafile.c_str());
  
  //opening data files and checking data (copied from Alex Tumanov's code)  
  const unsigned short *dduBuf=0;
  int length = 1;
 
  for (int event = 1; (event < maxEvents)&&(length); ++event) {//loop over all events in file (first one has empty DMB)
    std::cout << "---------- Event: " << event << "--------------" << std::endl;

    try {
      length= ddu.next(dduBuf);    
    } catch (std::runtime_error err ){
      std::cout <<"AFEB calibration:: " << err.what()<<"  end of file?" << std::endl;
      break;
    }
    unsigned short * buf = (unsigned short *) dduBuf;
        
    CSCDDUEventData dduEvent(buf);
    const std::vector<CSCEventData> & cscData = dduEvent.cscData(); 
    reportedChambers += dduEvent.header().ncsc();
    int NChambers = cscData.size();
    int repChambers = dduEvent.header().ncsc();
    std::cout << " Reported Chambers = " << repChambers <<"   "<<NChambers<< std::endl;
    if (NChambers!=repChambers) { std::cout<< "misMatched size!!!" << std::endl; misMatch++;}

  
    for (int i_chamber=0; i_chamber<NChambers; i_chamber++) {//loop over all DMBs
      
      for(int i_layer = 1; i_layer <= 6; ++i_layer) {//loop over all layers in chambers
	
	std::vector<CSCWireDigi> wire = cscData[i_chamber].wireDigis(i_layer) ;
	const CSCDMBHeader &thisDMBheader = cscData[i_chamber].dmbHeader();
	
	dmbID = cscData[i_chamber].dmbHeader().dmbID(); //get DMB ID
	crateID = cscData[i_chamber].dmbHeader().crateID(); //get crate ID
	if(crateID == 255) continue;//255 is reserved for old crate, present only 0 and 1
	
	for (unsigned int i=0; i < wire.size(); i++){ //loop over wire digis
	  int wireGroup = wire[i].getWireGroup();
	  int wireTBin = wire[i].getBeamCrossingTag();
	  
	  calib_evt.event      = event;
	  calib_evt.wire       = wireGroup;
	  calib_evt.tbin       = wireTBin;
	  calib_evt.chamber    = i_chamber;
	  calib_evt.layer      = i_layer;
	  
	  h1->Fill(wireGroup,wireTBin);
	  
	  calibtree->Fill();	    
	  
	}//end digis loop
      }//end layer loop
    }//end chamber loop
  }//end event loop
  
    
  //get chamber ID from Igor's mapping
  cout<<"Here is crate and DMB: "<<crateID<<"  "<<dmbID<<endl;
  map->crate_chamber(crateID,dmbID,&chamber_id,&chamber_num,&sector);
  cout<<"This is from mapping: "<< chamber_id<<"  "<<chamber_num<<"  "<<sector<<endl;
  
  //info needed for database
  string test1="CSC_slice";
  string test2="afeb";
  string test3="wire";
  
  //*******to send this array to DB uncomment the next two line(s)*************
  //cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, smth_val1,1, &ret_code);
  //cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, smth_val2,1, &ret_code);
  
  //root ntuple end
  calibfile->Write();   
  calibfile->Close();

}//main

