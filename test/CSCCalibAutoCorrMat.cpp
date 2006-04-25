/*

  Author: Stan Durkin

*/

#include "IORawData/CSCCommissioning/src/FileReaderDDU.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDDUEventData.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigi.h"
#include "EventFilter/CSCRawToDigi/interface/CSCEventData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCCLCTData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDDUTrailer.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDMBHeader.h"
#include "CalibMuon/CSCCalibration/interface/condbc.h"
#include "CalibMuon/CSCCalibration/interface/cscmap.h" 
#include "CalibMuon/CSCCalibration/interface/AutoCorrMat.h"
#include "FWCore/MessageService/interface/MessageServicePresence.h"
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

//constants declaration
#define CHAMBERS 5
#define LAYERS 6
#define STRIPS 80

int main(int argc, char **argv) {
  edm::service::MessageServicePresence my_message_service;
  //for debugging purpose from Alex Tumanov's code
  //set to true if you wish to debug data
  CSCAnodeData::setDebug(false);
  CSCALCTHeader::setDebug(false);
  CSCCLCTData::setDebug(false); 
  CSCEventData::setDebug(false);  
  CSCTMBData::setDebug(false);
  CSCDDUEventData::setDebug(false);

  //data file: nr.of events,chambers read
  int maxEvents = 50000;
  int misMatch = 0,fff;
  std::string datafile=argv[1];
  std::string chamber_id;  
  int dmbID[CHAMBERS],crateID[CHAMBERS],chamber_num,sector;
  int reportedChambers =0;
  int ret_code;  
  float newMatrix1[480];
  float newMatrix2[480];
  float newMatrix3[480];
  float newMatrix4[480];
  float newMatrix5[480];
  float newMatrix6[480];
  float newMatrix7[480];
  float newMatrix8[480];
  float newMatrix9[480];
  float newMatrix10[480];
  float newMatrix11[480];
  float newMatrix12[480];

  for (int i=0;i<480;i++){
    newMatrix1[i]=0.0;
    newMatrix2[i]=0.0;
    newMatrix3[i]=0.0;
    newMatrix4[i]=0.0;
    newMatrix5[i]=0.0;
    newMatrix6[i]=0.0;
    newMatrix7[i]=0.0;
    newMatrix8[i]=0.0;
    newMatrix9[i]=0.0;
    newMatrix10[i]=0.0;
    newMatrix11[i]=0.0;
    newMatrix12[i]=0.0;
  }

  //needed for database and mapping
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();
  
  if (argv[2]) maxEvents = (int) atof(argv[2]);
  
  FileReaderDDU ddu;
  ddu.open(datafile.c_str());

  //variable declaration
  std::vector<int> adc;
  Chamber_AutoCorrMat cam[5];
  for(int k=0;k<5;k++) cam[k].zero();
  
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
       
    for (int i_chamber=0; i_chamber<NChambers; i_chamber++) {//loop over all DMBs  
     
      for(int i_layer = 1; i_layer <= LAYERS; ++i_layer) {//loop over all layers in chambers
	std::vector<CSCStripDigi> digis = cscData[i_chamber].stripDigis(i_layer) ;
	const CSCDMBHeader &thisDMBheader = cscData[i_chamber].dmbHeader();

	if (thisDMBheader.cfebAvailable()){//check that CFEB data exists
	  dmbID[i_chamber]   = cscData[i_chamber].dmbHeader().dmbID(); //get DMB ID
	  crateID[i_chamber] = cscData[i_chamber].dmbHeader().crateID(); //get crate ID
	  if(crateID[i_chamber] == 255) continue; //255 is reserved for old crate, present only 0 and 1

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
  }

  float corrmat[12];
  float *tmp;
  tmp=corrmat; 

  std::string test1="CSC_slice";
  std::string test2="elem1";
  std::string test3="elem2";
  std::string test4="elem3";
  std::string test5="elem4";
  std::string test6="elem5";
  std::string test7="elem6";
  std::string test8="elem7";
  std::string test9="elem8";
  std::string test10="elem9";
  std::string test11="elem10";
  std::string test12="elem11";
  std::string test13="elem12";
  std::string answer;

  for (int i=0; i<CHAMBERS; i++){
    for (int j=0; j<LAYERS; j++){
      for (int k=0; k<STRIPS; k++){
	for (int max=0; max<12;max++){
	  fff = (j*80)+k;
	  tmp=cam[i].autocorrmat(j,k);
	  newMatrix1[fff]=tmp[0];
	  newMatrix2[fff]=tmp[1];
	  newMatrix3[fff]=tmp[2];
	  newMatrix4[fff]=tmp[3];
	  newMatrix5[fff]=tmp[4];
	  newMatrix6[fff]=tmp[5];
	  newMatrix7[fff]=tmp[6];
	  newMatrix8[fff]=tmp[7];
	  newMatrix9[fff]=tmp[8];
	  newMatrix10[fff]=tmp[9];
	  newMatrix11[fff]=tmp[10];
	  newMatrix12[fff]=tmp[11];

	  std::cout<<"Chamber "<<i<<" Layer "<<j<<" strip "<<k<<" Matrix elements "<<tmp[max]<<std::endl;
	}
      }
    }
    //get chamber ID from Igor's mapping

    int new_crateID = crateID[i];
    int new_dmbID   = dmbID[i];
    std::cout<<"Here is crate: "<<new_crateID<<" and DMB:  "<<new_dmbID<<std::endl;
    map->crate_chamber(new_crateID,new_dmbID,&chamber_id,&chamber_num,&sector);
    std::cout<<" Above data is for chamber: "<< chamber_id<<" from sector "<<sector<<std::endl;
     
    std::cout<<" DO you want to send constants to DB? "<<std::endl;
    std::cout<<" Please answer y or n for EACH chamber present! "<<std::endl;
    
    std::cin>>answer;
    if(answer=="y"){
      //SEND CONSTANTS TO DB
      
      cdb->cdb_write(test1,chamber_id,chamber_num,test2, 480, newMatrix1, 2, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test3, 480, newMatrix2, 2, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test4, 480, newMatrix3, 2, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test5, 480, newMatrix4, 2, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test6, 480, newMatrix5, 2, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test7, 480, newMatrix6, 2, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test8, 480, newMatrix7, 2, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test9, 480, newMatrix8, 2, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test10,480, newMatrix9, 2, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test11,480, newMatrix10,2, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test12,480, newMatrix11,2, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test13,480, newMatrix12,2, &ret_code);
    }else{
      std::cout<<" NO data was sent!!! "<<std::endl;
    }
  }
  
}//main
