/** 
 * Demo analyzer for reading digis
 * author A.Tumanov 2/22/06 
 * ripped from Jeremy's and Rick's analyzers
 *   
 */
#include <iostream>
#include <fstream>
#include <vector>
#include "string"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"
#include "CalibMuon/CSCCalibration/test/PedestalAnalyzer.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigi.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigi.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"
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
#include "DataFormats/CSCDigi/interface/CSCComparatorDigi.h"
#include "DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCRPCDigi.h"
#include "DataFormats/CSCDigi/interface/CSCRPCDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
// #include "CalibMuon/CSCCalibration/interface/condbc.h"
// #include "CalibMuon/CSCCalibration/interface/cscmap.h" 

#include "FWCore/MessageLogger/interface/MessageLogger.h"
//root specific .h files
#include </afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TROOT.h>
#include </afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TNtuple.h>
#include </afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TFile.h>
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TH1F.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TH1.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TCanvas.h"
#include "/afs/cern.ch/cms/external/lcg/external/root/5.08.00/slc3_ia32_gcc323/root/include/TTree.h"

PedestalAnalyzer::PedestalAnalyzer(edm::ParameterSet const& conf) {

  // If your module takes parameters, here is where you would define
  // their names and types, and access them to initialize internal
  // variables. Example as follows:
  //
  eventNumber = 0;
  evt = 0;
  pedMean=0.0,time=0.0,max =-9999999.,max1=-9999999.;
  pedSum = 0, strip =-999;

  //initialize arrays
  for(int i=0;i<CHAMBERS;i++){
    for(int j=0; j<LAYERS; j++){
      for(int k=0; k<STRIPS; k++){
	arrayOfPed[i][j][k]       = 0.;
	arrayOfPedSquare[i][j][k] = 0.;
	arrayPed[i][j][k]         = 0.;
	arrayPeak[i][j][k]        = 0.;
	arrayOfPeak[i][j][k]      = 0.; 
	arrayOfPeakSquare[i][j][k]= 0.;
	arraySumFive[i][j][k]     = 0.;
      }
    }
  }

  for (int i=0;i<480;i++){
    newPed[i]=0;
    newRMS[i]=0;
    newPeakRMS[i]=0.;
    newPeak[i]=0.;
    newSumFive[i]=0.;
  }
}

void PedestalAnalyzer::analyze(edm::Event const& e, edm::EventSetup const& iSetup) {
  
  // These declarations create handles to the types of records that you want
  // to retrieve from event "e".
  //
  // edm::Handle<CSCWireDigiCollection> wires;
  edm::Handle<CSCStripDigiCollection> strips;
  //edm::Handle<CSCComparatorDigiCollection> comparators;
  //edm::Handle<CSCALCTDigiCollection> alcts;
  //edm::Handle<CSCCLCTDigiCollection> clcts;
  //edm::Handle<CSCRPCDigiCollection> rpcs;
  //edm::Handle<CSCCorrelatedLCTDigiCollection> correlatedlcts;
  
  // Pass the handle to the method "getByType", which is used to retrieve
  // one and only one instance of the type in question out of event "e". If
  // zero or more than one instance exists in the event an exception is thrown.
  //
  //e.getByLabel("cscunpacker","MuonCSCWireDigi",wires);
  e.getByLabel("cscunpacker","MuonCSCStripDigi",strips);
  //e.getByLabel("cscunpacker","MuonCSCComparatorDigi",comparators);
  //e.getByLabel("cscunpacker","MuonCSCALCTDigi",alcts);
  //e.getByLabel("cscunpacker","MuonCSCCLCTDigi",clcts);
  //e.getByLabel("cscunpacker","MuonCSCRPCDigi",rpcs);
  //e.getByLabel("cscunpacker","MuonCSCCorrelatedLCTDigi",correlatedlcts);
  
  
  class TCalibEvt { public:
    Int_t adc[8];
    Float_t pedMean;
    Int_t strip;
    Float_t time[8];
    Int_t chamber;
    Int_t event;
    Int_t layer;
  };
  
  
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

  calibfile = new TFile("ntuples/calibpedestal.root","RECREATE","Calibration Ntuple");
  calibtree = new TTree("Calibration","Pedestal");
  calibevt = calibtree->Branch("EVENT", &calib_evt, "adc[8]/I:pedMean/F:strip/I:time[8]/F:chamber/I:event/I:layer/I");
  
  //data file: nr.of events,chambers read
  int misMatch = 0;
  //  std::string datafile=argv[1];
  std::string chamber_id;  
  int dmbID[CHAMBERS],crateID[CHAMBERS],chamber_num,sector;
  int i_chamber=0,i_layer=0;
  int reportedChambers =0;
  int fff,ret_code;  
  float aPeak=0.0,sumFive=0.0;
  const unsigned short *dduBuf=0;
  int length = 1;

  //needed for database and mapping
  // condbc *cdb = new condbc(); 
//   cscmap *map = new cscmap();
 
 //  for (int event = 0; (event < maxEvents)&&(length); ++event) {//loop over all events in file
//     std::cout << "---------- Event: " << event << "--------------" << std::endl;
    
   //  try {
//       length= ddu.next(dduBuf);    
//     } catch (std::runtime_error err ){
//       std::cout <<"Calibration:: " << err.what()<<"  end of file?" << std::endl;
//       break;
//     }

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
	const CSCDMBHeader &thisDMBheader = cscData[i_chamber].dmbHeader();

	if (thisDMBheader.cfebAvailable()){//check that CFEB data exists

	  dmbID[i_chamber]   = cscData[i_chamber].dmbHeader().dmbID(); //get DMB ID
	  crateID[i_chamber] = cscData[i_chamber].dmbHeader().crateID(); //get crate ID
	  if(crateID[i_chamber] == 255) continue; //255 is reserved for old crate, present only 0 and 1
	  
	  for (unsigned int i=0; i<digis.size(); i++){//loop over digis
	    strip = digis[i].getStrip();
	    adc   = digis[i].getADCCounts();
	    
	    pedSum  = adc[0]+adc[1];
	    pedMean = (float)pedSum/2.0;
	    
	    calib_evt.event    = evt;
	    calib_evt.pedMean  = pedMean;
	    calib_evt.strip    = strip;
	    calib_evt.chamber  = i_chamber;
	    calib_evt.layer    = i_layer;

	    for(unsigned int k=0;k<adc.size();k++){//loop over timeBins
	      time = (50. * k)-((evt%20)* 6.25)+116.5;	      
	      if(adc[3]>1000) aPeak = adc[3];
		
	      sumFive = adc[2]+adc[3]+adc[4];
	      
	      //if (max <aPeak){
	      //	max=aPeak;
	      // }
	      if (max1<sumFive){
		max1=sumFive;
	      }

	      calib_evt.adc[k]  = adc[k];
	      calib_evt.time[k]  = time;

	    }//adc.size
	    
	    calibtree->Fill();
	    
	    arrayPed[i_chamber][i_layer-1][strip-1] = pedMean;
	    arrayOfPed[i_chamber][i_layer - 1][strip - 1] += pedMean;
	    arrayOfPedSquare[i_chamber][i_layer - 1][strip - 1] += pedMean*pedMean ;

	    arrayPeak[i_chamber][i_layer-1][strip-1] = aPeak-pedMean;
	    arrayOfPeak[i_chamber][i_layer - 1][strip - 1] += aPeak-pedMean;
	    arrayOfPeakSquare[i_chamber][i_layer - 1][strip - 1] += (aPeak-pedMean)*(aPeak-pedMean);

	    arraySumFive[i_chamber][i_layer-1][strip-1] = (max1-pedMean)/(max-pedMean);
	    
	  }//end digis loop
	}//end if cfeb.available loop
      }//end layer loop
    }//end chamber loop
    //}//end events loop
  

  //root ntuple end
  calibfile->Write();   
  calibfile->Close();

  //create array (480 entries) for database transfer
  for(int myChamber=0; myChamber<CHAMBERS; myChamber++){
    double meanPedestal = 0.0,meanPeak=0.0,meanPeakSquare=0.;
    double meanPedestalSquare = 0.;
    double theRMS = 0.;
    double thePedestal =0.;
    double theRSquare = 0.;
    double thePeak =0.0,thePeakRMS=0.0;
    double theSumFive=0.0;

    std::string test1="CSC_slice";
    std::string test2="pedestal";
    std::string test3="ped_rms";
    std::string test4="peak_spread";
    std::string test5="pulse_shape";
    std::string answer;
    
    for (int i=0; i<CHAMBERS; i++){
      if (myChamber !=i) continue;
      
      for (int j=0; j<LAYERS; j++){
	for (int k=0; k<STRIPS; k++){
	  fff = (j*80)+k;
	  thePedestal  = arrayPed[i][j][k];
	  meanPedestal = arrayOfPed[i][j][k] / evt;
	  newPed[fff]  = meanPedestal;
	  meanPedestalSquare = arrayOfPedSquare[i][j][k] / evt;
	  theRMS       = sqrt(abs(meanPedestalSquare - meanPedestal*meanPedestal));
	  newRMS[fff]  = theRMS;
	  theRSquare   = (thePedestal-meanPedestal)*(thePedestal-meanPedestal)/(theRMS*theRMS*theRMS*theRMS); 
	  
	  thePeak = arrayPeak[i][j][k];
	  meanPeak = arrayOfPeak[i][j][k] / evt;
	  meanPeakSquare = arrayOfPeakSquare[i][j][k] / evt;
	  thePeakRMS = sqrt(abs(meanPeakSquare - meanPeak*meanPeak));
	  newPeakRMS[fff] = thePeakRMS;
	  newPeak[fff] = thePeak;

	  theSumFive = arraySumFive[i][j][k];
	  newSumFive[fff]=theSumFive;

	  std::cout <<" chamber "<<i<<" layer "<<j<<" strip "<<fff<<"  ped "<<newPed[fff]<<" RMS "<<newRMS[fff]<<" peakADC "<<newPeak[fff]<<" Peak RMS "<<newPeakRMS[fff]<<" Sum_of_four/apeak "<<newSumFive[fff]<<std::endl;
	}
      }
    }
    //get chamber ID from Igor's mapping

   //  int new_crateID = crateID[myChamber];
//     int new_dmbID   = dmbID[myChamber];
//     std::cout<<" Crate: "<<new_crateID<<" and DMB:  "<<new_dmbID<<std::endl;
//     map->crate_chamber(new_crateID,new_dmbID,&chamber_id,&chamber_num,&sector);
//     std::cout<<" Above data is for chamber: "<< chamber_id<<" from sector "<<sector<<std::endl;

//     std::cout<<" DO you want to send constants to DB? "<<std::endl;
//     std::cout<<" Please answer y or n for EACH chamber present! "<<std::endl;

//     std::cin>>answer;
//     if(answer=="y"){
      //SEND CONSTANTS TO DB
      // cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, newPed, 2, &ret_code);
//       cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, newRMS, 2, &ret_code);
//       cdb->cdb_write(test1,chamber_id,chamber_num,test4,480, newPeakRMS,2, &ret_code);
//       cdb->cdb_write(test1,chamber_id,chamber_num,test5,480, newSumFive,2, &ret_code);

    //   std::cout<<" Your results were sent to DB !!! "<<std::endl;
//     }else{
//       std::cout<<" NO data was sent!!! "<<std::endl;
//     }
  }

  eventNumber++;
  edm::LogInfo ("PedestalAnalyzer")  << "end of event number " << eventNumber;
  
}
