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

#include "condbc.h"
#include "cscmap.h" 
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

using namespace std;

class TCalibEvt { public:
  Int_t event;  
  Int_t chamber;
  Int_t layer;
  Int_t wire;
  Int_t tbin;
};


int main(int argc, char **argv) {
  CSCAnodeData::setDebug(false);
  CSCALCTHeader::setDebug(false);
  CSCCLCTData::setDebug(false); 
  CSCEventData::setDebug(false);  
  CSCTMBData::setDebug(false);
  CSCDDUEventData::setDebug(false);


  TCalibEvt calib_evt;
  TBranch *calibevt;
  TTree *calibtree;
  TFile *calibfile;
  
 calibfile = new TFile("calibwire.root","RECREATE","Calibration Ntuple");
 calibtree = new TTree("Calibration","AFEB");
 calibevt = calibtree->Branch("EVENT", &calib_evt, "event/I:chamber/I:layer/I:wire/I:tbin/I");
 TH2F *h1 = new TH2F("h1","wire vs timebin",10,-1,8,10,0,80);

  int maxEvents = 500000;
  int mismatch = 0;
  std::string datafile=argv[1];
  int dmbID=-999,crateID=-999,chamber_num,sector;
  int i_chamber=0,ilayer=0;
  int reportedChambers =0;
  int ret_code;
  std::string chamber_id;
  
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();

  if (argv[2]) maxEvents = (int) atof(argv[2]);

  FileReaderDDU ddu;
  ddu.open(datafile.c_str());
  
  int evt=0;
  std::vector<int> newadc;
      
  double a_pedestal[6][80];
  double a_pedestal_sq[6][80];
  double a_ped[6][80];
  int a_dmbID[6][80];
  int a_chamberID[6][80];
  double new_ped[480];
  double new_rms[480];
  
  for(int i=0;i<6;i++){
    for (int j=0; j<80; j++){
      a_pedestal[i][j] = 0.;
      a_pedestal_sq[i][j] = 0.;
      a_ped[i][j]=0.;
      a_dmbID[i][j]=0;
      a_chamberID[i][j]=0;
    }
  }
  
  for (int i=0;i<480;i++){
    new_ped[i]=0;
    new_rms[i]=0;
  }
  
  const unsigned short *dduBuf=0;
  int length = 1;
 
  for (int event = 1; (event < maxEvents)&&(length); ++event) {
    std::cout << "---------- Event: " << event << "--------------" << std::endl;
    try {
      length= ddu.next(dduBuf);    
    } catch (std::runtime_error err ){
      std::cout <<"digi Anal:: " << err.what()<<"  end of file?" << std::endl;
      break;
    }
    unsigned short * buf = (unsigned short *) dduBuf;
        
    CSCDDUEventData dduEvent(buf);
    const std::vector<CSCEventData> & cscData = dduEvent.cscData(); 
    reportedChambers += dduEvent.header().ncsc();
    int NChambers = cscData.size();
    int repChambers = dduEvent.header().ncsc();
    std::cout << " Reported Chambers = " << repChambers <<"   "<<NChambers<< std::endl;
    if (NChambers!=repChambers) { std::cout<< "mismatched size!!!" << std::endl; mismatch++;}

    evt++;

    for (i_chamber=0; i_chamber<NChambers; i_chamber++) {
      for(ilayer = 1; ilayer <= 6; ++ilayer) {
	std::vector<CSCWireDigi> wire = cscData[i_chamber].wireDigis(ilayer) ;
	const CSCDMBHeader &thisDMBheader = cscData[i_chamber].dmbHeader();
	dmbID = cscData[i_chamber].dmbHeader().dmbID();
	crateID = cscData[i_chamber].dmbHeader().crateID();
	if(crateID == 255) continue;
	for (int i=0; i < wire.size(); i++){
	  int wire_group = wire[i].getWireGroup();
	  int wire_tbin = wire[i].getBeamCrossingTag();

	    calib_evt.event      = event;
	    calib_evt.wire       = wire_group;
	    calib_evt.tbin       = wire_tbin;
	    calib_evt.chamber    = i_chamber;
	    calib_evt.layer      = ilayer;
	    h1->Fill(wire_group,wire_tbin);
 
	    calibtree->Fill();	    
	    
	}//wire.size
      }//ilayer
    }//chamber
  }//event loop
  
  
  string test1="CSC_slice";
  string test2="pedestal";
  string test3="ped_rms";
  
  cout<<"Here is crate and DMB: "<<crateID<<"  "<<dmbID<<endl;
  map->crate_chamber(crateID,dmbID,&chamber_id,&chamber_num,&sector);
  cout<<"This is from mapping: "<< chamber_id<<"  "<<chamber_num<<"  "<<sector<<endl;

  // To send infor to DB uncomment the next line(s)  
  //cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, new_ped,2, &ret_code);
  //cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, new_rms,2, &ret_code);
  calibfile->Write();   
  calibfile->Close();
}//main

