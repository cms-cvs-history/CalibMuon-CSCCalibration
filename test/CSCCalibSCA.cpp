/*
  authors: Oana Boeriu, Nicole Ippolito
  
  This will create arrays (480=whole ME2/3 chamber) of pedestal and RMS values
  which can be sent to online database.
  It creates a root ntuple for debugging purpose.  
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
#include "EventFilter/CSCRawToDigi/interface/CSCCFEBData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCCFEBTimeSlice.h"


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
  Float_t pedMean;
  Int_t strip;
  Float_t time[8];
  Int_t chamber;
  Int_t event;
  Int_t layer;
  
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

  calibfile = new TFile("ntuples/calibsca.root","RECREATE","Calibration Ntuple");
  calibtree = new TTree("Calibration","SCA Pedestal");
  calibevt = calibtree->Branch("EVENT", &calib_evt, "adc[8]/I:pedMean/F:strip/I:time[8]/F:chamber/I:event/I:layer/I");
  
  //data file: nr.of events,chambers read
  int maxEvents = 900000;
  int misMatch = 0;
  std::string datafile=argv[1];
  //std::string chamber_id;  
  int dmbID[CHAMBERS],crateID[CHAMBERS];//,chamber_num,sector;
  int i_chamber=0,i_layer=0;
  int reportedChambers =0;
  //int fff,ret_code;  
  
  //needed for database and mapping
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();
  
  if (argv[2]) maxEvents = (int) atof(argv[2]);
  
  FileReaderDDU ddu;
  ddu.open(datafile.c_str());
  
  //some variable declaration
  int evt=0;
  // std::vector<int> newadc;
  std::vector<int> adc;
  float pedMean=0.0,time=0.0,pedSum=0.0;
  int strip =-999;
  unsigned short scablock,trigtime,lctphase;
  int cap, sca_number;

  //definition of arrays
  double arrayOfPed[CHAMBERS][LAYERS][STRIPS];
  double arrayOfPedSquare[CHAMBERS][LAYERS][STRIPS];
  double arrayPed[CHAMBERS][LAYERS][STRIPS];

  double newPed[480];
  double newRMS[480];
  
  //initialize arrays
  for(int i=0;i<CHAMBERS;i++){
    for(int j=0; j<LAYERS; j++){
      for(int k=0; k<STRIPS; k++){
	arrayOfPed[i][j][k] = 0.;
	arrayOfPedSquare[i][j][k] = 0.;
	arrayPed[i][j][k]=0.;
      }
    }
  }

  for (int i=0;i<480;i++){
    newPed[i]=0;
    newRMS[i]=0;
  }
  
  //opening data files and checking data (copied from Alex Tumanov's code)
  const unsigned short *dduBuf=0;
  int length = 1;
  int power;

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
      const CSCDMBHeader &thisDMBheader = cscData[i_chamber].dmbHeader();
      
      if (thisDMBheader.cfebAvailable()){//check that CFEB data exists
	dmbID[i_chamber]   = cscData[i_chamber].dmbHeader().dmbID(); //get DMB ID
	crateID[i_chamber] = cscData[i_chamber].dmbHeader().crateID(); //get crate ID
	if(crateID[i_chamber] == 255) continue; //255 is reserved for old crate, present only 0 and 1
	
	  for (int icfeb=0; icfeb<5;icfeb++) {//loop over cfebs in a given chamber
	    CSCCFEBData * mycfebData =  cscData[i_chamber].cfebData(icfeb);
	    if (mycfebData){

	      for (int itime=0; itime<8; itime++) {//loop over time samples
		CSCCFEBTimeSlice * mytimeSlice =  mycfebData->timeSlice(itime);
		if (mytimeSlice){
		  for(i_layer = 1; i_layer <= 6; i_layer++) {//loop over all layers in chambers 
		    std::vector<CSCStripDigi> digis = cscData[i_chamber].stripDigis(i_layer) ;
		    scablock = mytimeSlice->scaControllerWord(i_layer).sca_blk;
		    trigtime = mytimeSlice->scaControllerWord(i_layer).trig_time;
		    lctphase = mytimeSlice->scaControllerWord(i_layer).lct_phase;
		    int tmp=1;
		    for(power=0;power<8;power++){if(trigtime==tmp) lctphase=power; tmp=tmp*2;}
		    
		    cap = lctphase+itime;
		    sca_number=8*scablock+cap;
		    
		    // for (unsigned int i=0; i<digis.size(); i++){//loop over digis
		      
// 		      strip = digis[i].getStrip();
// 		      adc   = digis[i].getADCCounts();
// 		      pedSum  = adc[0]+adc[1];
// 		      pedMean = (float)pedSum/2.0;
// 		      calib_evt.strip    = strip;
// 		      calib_evt.event    = event;
// 		      calib_evt.pedMean  = pedMean;
// 		      calib_evt.chamber  = i_chamber;
// 		      calib_evt.layer    = i_layer;
		      
		    //std::cout <<"Layer "<<i_layer<<" Strip "<<strip<<" sca_block "<<scablock <<" trig_time "<<trigtime<<" lct_phase  "<<lctphase<<" sca_number "<<sca_number<<" time slice " <<itime<< std::endl;
    
		    std::cout <<"CFEB "<<icfeb<<" Layer "<<i_layer<<" sca_block "<<scablock <<" trig_time "<<trigtime<<" lct_phase "<<lctphase<<" sca_number "<<sca_number<<" time slice "<<itime<< std::endl;
		      //}
		  }
		}else std::cout<<"no time slice" <<std::endl;
	      }
	    }else std::cout<<"no cfeb data"<< std::endl;
	  }
      }
      //calibtree->Fill();
    }//chambers
    
    
    // 		     arrayPed[i_chamber][i_layer-1][strip-1] = pedMean;
    // 		     arrayOfPed[i_chamber][i_layer - 1][strip - 1] += pedMean;
    // 		     arrayOfPedSquare[i_chamber][i_layer - 1][strip - 1] += pedMean*pedMean ;
    
    
    
    
    //root ntuple end
    //calibfile->Write();   
    calibfile->Close();
  }//event
  
  //create array (480 entries) for database transfer
  //  for(int myChamber=0; myChamber<CHAMBERS; myChamber++){
  //     double meanPedestal = 0.;
  //     double meanPedestalSquare = 0.;
  //     double theRMS = 0.;
  //     double thePedestal =0.;
  //     double theRSquare = 0.;
  //     std::string test1="CSC_slice";
  //     std::string test2="pedestal";
  //     std::string test3="ped_rms";
  //     std::string answer;
  
  //     for (int i=0; i<CHAMBERS; i++){
  //       if (myChamber !=i) continue;
  
  //       for (int j=0; j<LAYERS; j++){
  // 	for (int k=0; k<STRIPS; k++){
  // 	  fff = (j*80)+k;
  // 	  thePedestal  = arrayPed[i][j][k];
  // 	  meanPedestal = arrayOfPed[i][j][k] / evt;
  // 	  newPed[fff]  = meanPedestal;
  // 	  meanPedestalSquare = arrayOfPedSquare[i][j][k] / evt;
  // 	  theRMS       = sqrt(abs(meanPedestalSquare - meanPedestal*meanPedestal));
  // 	  newRMS[fff]  = theRMS;
  // 	  theRSquare   = (thePedestal-meanPedestal)*(thePedestal-meanPedestal)/(theRMS*theRMS*theRMS*theRMS); 
  
  // 	  std::cout <<" chamber "<<i<<" layer "<<j<<" strip "<<fff<<"  ped "<<newPed[fff]<<" RMS "<<newRMS[fff]<<std::endl;
  // 	}
  //       }
  //     }
  //     //get chamber ID from Igor's mapping
  
  //     int new_crateID = crateID[myChamber];
  //     int new_dmbID   = dmbID[myChamber];
  //     std::cout<<" Crate: "<<new_crateID<<" and DMB:  "<<new_dmbID<<std::endl;
  //     map->crate_chamber(new_crateID,new_dmbID,&chamber_id,&chamber_num,&sector);
  //     std::cout<<" Above data is for chamber: "<< chamber_id<<" from sector "<<sector<<std::endl;
  
  //     std::cout<<" DO you want to send constants to DB? "<<" Please answer y or n for EACH chamber present! "<<std::endl;
  
  //     std::cin>>answer;
  //     if(answer=="y"){
  //       //SEND CONSTANTS TO DB
  //       cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, newPed,2, &ret_code);
  //       cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, newRMS,2, &ret_code);
  
  //       std::cout<<" Data SENT to DB !!! "<<std::endl;
  //     }else{
  //       std::cout<<" NO data was sent!!! "<<std::endl;
  //     }
  //   }
  
  
}//main


