#include "IORawData/CSCCommissioning/src/FileReaderDDU.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDDUEventData.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigi.h"
#include "EventFilter/CSCRawToDigi/interface/CSCEventData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCAnodeData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCALCTHeader.h"
#include "EventFilter/CSCRawToDigi/interface/CSCCLCTData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDDUTrailer.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDMBHeader.h"
#include "CalibMuon/CSCCalibration/interface/condbc.h"
#include "CalibMuon/CSCCalibration/interface/cscmap.h" 
#include "CalibMuon/CSCCalibration/interface/AutoCorrMat.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "string"

//root specific .h files
#include </afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TROOT.h>
#include </afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TNtuple.h>
#include </afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TFile.h>
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TH1F.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TH1.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TCanvas.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TTree.h"

int main(int argc, char **argv) {

  //for debugging purposes from Alex Tumanov's code
  //set to true if you wish to debug data
  CSCAnodeData::setDebug(false);
  CSCALCTHeader::setDebug(false);
  CSCCLCTData::setDebug(false); 
  CSCEventData::setDebug(false);  
  CSCTMBData::setDebug(false);
  CSCDDUEventData::setDebug(false);

  //data file: nr.of events,chambers read
  int maxEvents = 50000;
  int misMatch = 0;
  std::string datafile=argv[1];
  std::string chamber_id;  
  int dmbID=-999,crateID=-999,chamber_num,sector;
  int i_chamber=0,i_layer=0;
  int reportedChambers =0;
  int fff,ret_code;  

  //needed for database and mapping
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();
  
  if (argv[2]) maxEvents = (int) atof(argv[2]);
  
  FileReaderDDU ddu;
  ddu.open(datafile.c_str());

  //some variable declaration
  int evt=0;

  std::vector<int> adc;
  Chamber_AutoCorrMat cam[5];
  for(int k=0;k<5;k++)cam[k].zero();
  //opening data files and checking data (copied from Alex Tumanov's code)
  const unsigned short *dduBuf=0;
  int length = 1;


  for (int event = 0; (event < maxEvents)&&(length); ++event) {//loop over all events in file
    std::cout << "---------- Event: " << event << "--------------" << std::endl;
    
    try {
      length= ddu.next(dduBuf);    
    } catch (std::runtime_error err ){
      std::cout <<"Calibration:: " << err.what()<<"  end of file?" << std::endl;
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
       
    evt++;

    for (i_chamber=0; i_chamber<NChambers; i_chamber++) {//loop over all DMBs  
     
      for(i_layer = 1; i_layer <= 6; ++i_layer) {//loop over all layers in chambers

	std::vector<CSCStripDigi> digis = cscData[i_chamber].stripDigis(i_layer) ;
       	for (unsigned int i=0; i<digis.size(); i++){//loop over digis
	  int strip = digis[i].getStrip();
	  adc = digis[i].getADCCounts();
          int tadc[8];
          for(unsigned int j=0;j<adc.size();j++)tadc[j]=adc[j];
          cam[i_chamber].add(i_layer-1,strip-1,tadc);
        }
      }
    }
  }
  // correlation matrices calculated so dump some
  float corrmat[12];
  float *tmp;
  tmp=corrmat;
  int lay=4;int strip=10;
  tmp=cam[0].autocorrmat(lay,strip);
  for(int i=0;i<12;i++)printf(" %5.2f ",tmp[i]);printf("\n");
  lay=4;strip=11;
  tmp=cam[0].autocorrmat(lay,strip);
  for(int i=0;i<12;i++)printf(" %5.2f ",tmp[i]);printf("\n");

}
