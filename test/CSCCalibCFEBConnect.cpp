/*

  authors: Oana Boeriu, Nicole Ippolito
  This will create arrays (480=whole ME2/3 chamber) of differences between max and min 
  ADC values. Creates root ntuple for debugging.

*/

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

class TCalibEvt { public:
  Int_t adc[8];
  Int_t strip;
  Int_t chamber;
  Int_t event;
  Int_t layer;
  Int_t adc_max[8];
  Int_t adc_min[8];

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

  calibfile = new TFile("calibstripcon.root","RECREATE","Calibration Ntuple");
  calibtree = new TTree("Calibration","Connectivity");
  calibevt = calibtree->Branch("EVENT", &calib_evt, "adc[8]/I:strip/I:chamber/I:event/I:layer/I:adc_max[8]/I:adc_min[8]/I");

  //data file: nr.of events,chambers read
  int maxEvents = 500000;
  int misMatch = 0;
  std::string datafile=argv[1];
  int dmbID[CHAMBERS],crateID[CHAMBERS],chamber_num,sector;
  int reportedChambers =0;
  int fff,ret_code,strip=-999;  
  std::string chamber_id;
  
  //needed for database and mapping
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();

  if (argv[2]) maxEvents = (int) atof(argv[2]);

  FileReaderDDU ddu;
  ddu.open(datafile.c_str());
  
  //some variable declaration
  std::vector<int> newadc;
  int evt=0;
      
  //definition of arrays
  double new_diff[480];
  int adcMin[CHAMBERS][LAYERS][STRIPS], adcMax[CHAMBERS][LAYERS][STRIPS];
  double diff[CHAMBERS][LAYERS][STRIPS];


  //initialize arrays
  for (int i=0; i<CHAMBERS; i++){
    for (int j=0; j<LAYERS; j++){
      for (int k=0; k<STRIPS;k++){
	adcMin[i][j][k]    = 9999999;
	adcMax[i][j][k]    =-9999999;
	diff[i][j][k]      = 0.;
      }
    }
  }

  for (int i=0;i<480;i++){
    new_diff[i]= 0;
  }

  //opening data files and checking data (copied from Alex Tumanov's code)
  const unsigned short *dduBuf=0;
  int length = 1;
  
  for (int event = 0; (event < maxEvents)&&(length); ++event) {//loop over all events in file
       
    std::cout << "---------- Event: " << event << "--------------" << std::endl;
    try {
      length= ddu.next(dduBuf);    
    } catch (std::runtime_error err ){
      std::cout <<"digi Anal:: " << err.what()<<"  end of file?" << std::endl;
      break;
    }

    evt++;

    unsigned short * buf = (unsigned short *) dduBuf; 
    CSCDDUEventData dduEvent(buf);
    
    const std::vector<CSCEventData> & cscData = dduEvent.cscData(); 
    reportedChambers += dduEvent.header().ncsc();
    int NChambers = cscData.size();
    int repChambers = dduEvent.header().ncsc();
    std::cout << " Reported Chambers = " << repChambers <<"   "<<NChambers<< std::endl;
    if (NChambers!=repChambers) { std::cout<< "misMatched size!!!" << std::endl; misMatch++;}

    //
    for (int i_chamber=0; i_chamber<NChambers; i_chamber++) {//loop over all DMBs

      for(unsigned i_layer = 1; i_layer <=6; ++i_layer) {//loop over all layers in chambers

	std::vector<CSCStripDigi> digis = cscData[i_chamber].stripDigis(i_layer) ;
	const CSCDMBHeader &thisDMBheader = cscData[i_chamber].dmbHeader();

	if (thisDMBheader.cfebAvailable()){//check that CFEB data exists
	  dmbID[i_chamber] = cscData[i_chamber].dmbHeader().dmbID();//get DMB ID
	  crateID[i_chamber] = cscData[i_chamber].dmbHeader().crateID();//get crate ID
	  if(crateID[i_chamber] == 255) continue; //255 is reserved for old crate, present only 0 and 1
	  
	  for (unsigned int i=0; i<digis.size(); i++){//loop over digis
	    std::vector<int> adc = digis[i].getADCCounts();
	    strip = digis[i].getStrip();

	    calib_evt.event    = event;
            calib_evt.strip    = strip;
            calib_evt.chamber  = i_chamber;
            calib_evt.layer    = i_layer;

	    for (unsigned int k=0;k<adc.size();k++){//loop over timeBins
	      
	      //find max ADC value
	      if(adc[k] > adcMax[i_chamber][i_layer-1][strip-1]) {
		adcMax[i_chamber][i_layer-1][strip-1]= adc[k];
	      }

	      //find min ADC value
	      if(adc[k] < adcMin[i_chamber][i_layer-1][strip-1]){
		adcMin[i_chamber][i_layer-1][strip-1]=adc[k];
	      }
	      
		calib_evt.adc_max[k] = adcMax[i_chamber][i_layer-1][strip-1];
		calib_evt.adc_min[k] = adcMin[i_chamber][i_layer-1][strip-1];
	    }//end timeBin loop

	    calibtree->Fill();
	  }//end digis loop
	}//end cfeb.available
      }//end layer loop
    }//end chamber loop
  }//end events loop

  for(int mychamber=0; mychamber<CHAMBERS;mychamber++){
    for(int mylayer=0; mylayer<LAYERS;mylayer++){
      for(int mystrip=0;mystrip<STRIPS;mystrip++){
	diff[mychamber][mylayer][mystrip]=adcMax[mychamber][mylayer][mystrip]-adcMin[mychamber][mylayer][mystrip];
	fff=(mylayer*80)+mystrip;
	new_diff[fff] =diff[mychamber][mylayer][mystrip];
	//std::cout<<"Chamber "<<mychamber<<" Layer "<<mylayer<<" Strip "<<mystrip<<" diff "<<new_diff[fff]<<std::endl; 
      }
    }
    //get chamber ID from Igor's mapping
    int new_dmbID = dmbID[mychamber];
    int new_crateID = crateID[mychamber];
    cout<<"Here is crate and DMB: "<<new_crateID<<"  "<<new_dmbID<<endl;
    map->crate_chamber(new_crateID,new_dmbID,&chamber_id,&chamber_num,&sector);
    cout<<"This is from mapping: "<< chamber_id<<"  "<<chamber_num<<"  "<<sector<<endl;
    
    //info needed for database
    string test1="CSC_slice";
    string test2="cfeb_con_diff";
    string test3="con_con_rms";
    
    //*******to send this array to DB uncomment the next two lines*************
    //cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, new_diff,1, &ret_code);
    //cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, new_rms, 1, &ret_code);
  }
    
  //root ntuple end
  calibfile->Write();
  calibfile->Close();

}// main




