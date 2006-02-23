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
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TH1.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TCanvas.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TTree.h"


using namespace std;

class TCalibEvt { public:
  Int_t adc[8];
    Float_t ped_mean;
    Int_t strip;
    Float_t time[8];
    Int_t chamber;
    Int_t event;
    Int_t layer;
  
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

  calibfile = new TFile("calibntuple676.root","RECREATE","Calibration Ntuple");
  calibtree = new TTree("Calibration","Pedestal");
  calibevt = calibtree->Branch("EVENT", &calib_evt, "adc[8]/I:ped_mean/F:strip/I:time[8]/F:chamber/I:event/I:layer/I");

  int maxEvents = 900000;
  int misMatch = 0;
  std::string datafile=argv[1];
  std::string chamber_id;  
  int dmbID=-999,crateID=-999,chamber_num,sector;
  int i_chamber=0,i_layer=0;
  int reportedChambers =0;
  int fff,ret_code;  
  
  
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();

  if (argv[2]) maxEvents = (int) atof(argv[2]);

  FileReaderDDU ddu;
  ddu.open(datafile.c_str());
  
  int evt=0;
  std::vector<int> newadc;
  std::vector<int> adc;
  float pedMean=0.0,t=0.0;
  int pedSum = 0, strip =-999;

  //arrays
  double arrayOfPed[6][80];
  double arrayOfPedSquare[6][80];
  double arrayPed[6][80];
  double newPed[480];
  double newRMS[480];
  
  for(int i=0;i<6;i++){
    for (int j=0; j<80; j++){
      arrayOfPed[i][j] = 0.;
      arrayOfPedSquare[i][j] = 0.;
      arrayPed[i][j]=0.;
    }
  }
  
  for (int i=0;i<480;i++){
    newPed[i]=0;
    newRMS[i]=0;
  }
  
  const unsigned short *dduBuf=0;
  int length = 1;


  for (int event = 0; (event < maxEvents)&&(length); ++event) {
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

    for (i_chamber=0; i_chamber<NChambers; i_chamber++) {  
      for(i_layer = 1; i_layer <= 6; ++i_layer) {
	std::vector<CSCStripDigi> digis = cscData[i_chamber].stripDigis(i_layer) ;
	const CSCDMBHeader &thisDMBheader = cscData[i_chamber].dmbHeader();
	if (thisDMBheader.cfebAvailable()){
	  dmbID = cscData[i_chamber].dmbHeader().dmbID();
	  crateID = cscData[i_chamber].dmbHeader().crateID();
	  if(crateID == 255) continue;
	  
	  for (int i=0; i<digis.size(); i++){
	    strip = digis[i].getStrip();
	    adc = digis[i].getADCCounts();
	    
	    pedSum =adc[0]+adc[1];
	    pedMean = (float)pedSum/2.0;
	    
	    calib_evt.event    = event;
	    calib_evt.pedMean = pedMean;
	    calib_evt.strip    = strip;
	    calib_evt.chamber  = i_chamber;
	    calib_evt.layer    = i_layer;

	    int offset = event / 20;
	    

	    for(int k=0;k<adc.size();k++){
	      
	      t = (50*k)-((event)*6.25)+194+(200*offset);
	          	      
	      calib_evt.adc[k]      = adc[k];
	      calib_evt.time[k]     = t;

	    }//adc.size

	      calibtree->Fill();

	    arrayPed[i_layer-1][strip-1] = pedMean;
	    arrayOfPed[i_layer - 1][strip - 1] += pedMean;
	    arrayOfPedSquare[i_layer - 1][strip - 1] += pedMean*pedMean ;
	  }//digi.size()
	}//cfeb available
      }//layer
    }//chamber
  }//events
  
  for (int i=0; i<6; i++){
    for (int j=0; j<80; j++){
      double mean_ped = 0.;
      fff = (i*80)+j;
      double mean_ped_sq = 0.;
      double the_rms = 0.;
      double the_ped =0.;
      double the_rsquare = 0.;
      
      the_ped = arrayPed[i][j];
      mean_ped = arrayOfPed[i][j] / evt;
      newPed[fff] = mean_ped;
      mean_ped_sq = arrayOfPedSquare[i][j] / evt;
      the_rms = sqrt(abs(mean_ped_sq - mean_ped*mean_ped));
      newRMS[fff] = the_rms;
      the_rsquare = (the_ped-mean_ped)*(the_ped-mean_ped)/(the_rms*the_rms*the_rms*the_rms); 
    }
  }
  
  string test1="CSC_slice";
  string test2="pedestal";
  string test3="ped_rms";
  
  cout<<"Here is crate and DMB: "<<crateID<<"  "<<dmbID<<endl;
  map->crate_chamber(crateID,dmbID,&chamber_id,&chamber_num,&sector);
  cout<<"This is from mapping: "<< chamber_id<<"  "<<chamber_num<<"  "<<sector<<endl;
  
  // to send this array to DB uncomment thenext two lines!

  //cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, newPed,2, &ret_code);
  //cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, newRMS,2, &ret_code);
  

  calibfile->Write();   
  calibfile->Close();
}//main


