//This calculates gains for CFEB, makes root ntuple for debugging, array of 480 entries,sends it to DB.
//Every strip pulsed 10 times,each time increase charge amplitude by 0.2 (max amplitude charge=2.2),
//each step 0.2 corresponds to 22.4fC (correct?? must check).
//First 100 events give strip=1,second 100 events give strip=2 and so on.

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
#define NUMMODTEN 200
#define NUMLAYERS 6
#define NUMSTRIPS  80
#define NUMCHAMBERS 1
#define NUMBERPLOTTED 20


class TCalibEvt { public:
  Float_t max_ADC[6][80];
  Float_t charge[6][80];
  Float_t slope[6][80];
  Float_t intercept[6][80];
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
  
  calibfile = new TFile("calibgains.root","RECREATE","Calibration Ntuple");
  calibtree = new TTree("Calibration","Gains");
  calibevt = calibtree->Branch("EVENT", &calib_evt, "max_ADC[6][80]/F:charge[6][80]/F:slope[6][80]/F:intercept[6][80]/F");

  //data file: nr.of events,chambers read
  int maxEvents = 200;
  int misMatch = 0;
  std::string datafile=argv[1];
  int dmbID=-999,crateID=-999,chamber_num,sector;
  int reportedChambers =0;
  int fff,ret_code,strip=-999;  
  std::string chamber_id;

  //needed for database and mapping
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();

  if (argv[2]) maxEvents = (int) atof(argv[2]);

  FileReaderDDU ddu;
  ddu.open(datafile.c_str());

  //some variables declaration  
  std::vector<int> newadc;
  float min1=9999999999999.0;
  float sumOfX=0.0, sumx2=0.0,sumOfXY=0.0;
  float sumOfY=0.0;
  float chi2=0.0;
  int evt=0;
  float gainSlope=-999.0,gainIntercept=-999.0;

  //charge injected in amplitude steps of 0.2 == 22.4fC steps
 
  const float x[20] = {22.4, 44.8, 67.2, 89.6, 112, 134.4, 156.8, 179.2, 201.6, 224.0, 246.4, 268.8, 291.2, 313.6, 336.0, 358.4, 380.8, 403.2, 425.6, 448};
 
  //definition of arrays
  double adcMax[6][80];
  double adcMean_max[6][80];
  double arrayOfGain[6][80];
  double arrayOfGainSquare[6][80];
  double arrayOfIntercept[6][80];
  double arrayOfInterceptSquare[6][80];
  double newGain[480];
  double newIntercept[480];
  float maxmodten[NUMMODTEN][NUMLAYERS][NUMSTRIPS];


 //initialize arrays
  for (int i=0; i<NUMMODTEN; i++){
    for (int j=0; j<NUMLAYERS; j++){
      for (int k=0; k<NUMSTRIPS; k++){
	maxmodten[i][j][k] = -999.;
      }
    }
  }
  
  for (int i=0; i<6; i++){
    for (int j=0; j<80; j++){
      arrayOfGain[i][j] = 0.;
      arrayOfGainSquare[i][j] = 0.;
      arrayOfGain[i][j]=0.;

      arrayOfIntercept[i][j]=0.;
      arrayOfInterceptSquare[i][j]=0.;
      adcMax[i][j]=-999.;
      adcMean_max[i][j]=0.0;
    }
  }

    
  for (int i=0;i<480;i++){
    newGain[i]=0;
    newIntercept[i]=0;
  }

  //opening data files and checking data (copied from Alex Tumanov's code)
  const unsigned short *dduBuf=0;
  int length = 1;
  

  //************** START **************//

  for (int event = 0; (event < maxEvents)&&(length); ++event) { //loop over all events in file
      
    std::cout << "---------- Event: " << event << "--------------" << std::endl;
    try {
      length= ddu.next(dduBuf);    
    } catch (std::runtime_error err ){
      std::cout <<"digi Anal:: " << err.what()<<"  end of file?" << std::endl;
      break;
    }
    unsigned short * buf = (unsigned short *) dduBuf; 

    CSCDDUEventData dduEvent(buf);
    evt++;

    const std::vector<CSCEventData> & cscData = dduEvent.cscData(); 
    reportedChambers += dduEvent.header().ncsc();
    int NChambers = cscData.size();
    int repChambers = dduEvent.header().ncsc();
    std::cout << " Reported Chambers = " << repChambers <<"   "<<NChambers<< std::endl;
    if (NChambers!=repChambers) { std::cout<< "misMatched size!!!" << std::endl; misMatch++;}
    

    for (int i_chamber=0; i_chamber<1; i_chamber++) {//loop over all DMBs      

      for(int i_layer = 1; i_layer <=6; ++i_layer) {//loop over all layers in chambers

	std::vector<CSCStripDigi> digis = cscData[i_chamber].stripDigis(i_layer) ;
	const CSCDMBHeader &thisDMBheader = cscData[i_chamber].dmbHeader();

	if (thisDMBheader.cfebAvailable()){
	  dmbID = cscData[i_chamber].dmbHeader().dmbID();//get DMB ID
	  crateID = cscData[i_chamber].dmbHeader().crateID();//get crate ID
	  if(crateID == 255) continue;
	  
	  for (unsigned int i=0; i<digis.size(); i++){//loop over digis
	    std::vector<int> adc = digis[i].getADCCounts();
	    strip = digis[i].getStrip();

	    for(unsigned int k=0;k<adc.size();k++){//loop over timeBins

	      if(adc[k] > adcMax[i_layer-1][strip-1]) {
		adcMax[i_layer-1][strip-1]=adc[k];
		adcMean_max[i_layer-1][strip-1] += adcMax[i_layer-1][strip-1]/8;
		
	      }
	    }//adc.size

	    // On the 10th event
	    if (evt%10 == 0){
	      int ten = int(evt/10) - 1;
	      maxmodten[ten][i_layer-1][strip-1] = adcMean_max[i_layer-1][strip-1];	     
	    }
	  }//end digis loop
	}//end cfeb.available loop
      }//end layer loop
    }//end chamber loop

    // Fill the ntuple every 10th event
    if (evt%10 == 0){
      for(unsigned thelayer = 0; thelayer<6; thelayer++) {
	for (unsigned thestrip=0; thestrip<80; thestrip++){
	  calib_evt.max_ADC[thelayer][thestrip] = maxmodten[int(evt/10) - 1][thelayer][thestrip];
	  calib_evt.charge[thelayer][thestrip]=x[int(evt/10) - 1];
	}
      }
      
      calibtree->Fill();
    }
  }//end event loop
  

  //create array (480 entries) for database transfer
  for (int layeriter=0; layeriter<6; layeriter++){
    for (int stripiter=0; stripiter<80; stripiter++){

      for (int j=0; j<6; j++){//layer	
	if (j != layeriter) continue;
	for (int k=0; k<80; k++){//strip
	  if (k != stripiter) continue;
	  sumOfX = 0.;
	  sumOfY = 0.;
	  sumOfXY = 0.;
	  sumx2 = 0.;
	  gainSlope = 0.;
	  gainIntercept = 0.;
	  chi2 = 0.;

	  for(int ii=0;ii<20;ii++){//numbers       
	    //do straight line fit ( y = kx + m )

	    sumOfX += x[ii];
	    sumOfY += maxmodten[ii][j][k];
	    sumOfXY += (x[ii]*maxmodten[ii][j][k]);
	    sumx2 += (x[ii]*x[ii]);
	    
	    gainSlope= ((NUMBERPLOTTED*sumOfXY) - (sumOfX * sumOfY))/((NUMBERPLOTTED*sumx2) - (sumOfX*sumOfX));//k
	    gainIntercept = ((sumOfY*sumx2)-(sumOfX*sumOfXY))/((NUMBERPLOTTED*sumx2)-(sumOfX*sumOfX));//m
	    chi2  += (maxmodten[ii][j][k]-(gainIntercept+(gainSlope*x[ii])))*(maxmodten[ii][j][k]-(gainIntercept+(gainSlope*x[ii])));
	    
	    if(chi2<min1){
	      min1=chi2;
	      gainSlope=gainSlope;
	      gainIntercept=gainIntercept;
	    }
	  }
	  
	  arrayOfGain[j][k] += gainSlope;
	  arrayOfGainSquare[j][k] += gainSlope*gainSlope;
	  arrayOfIntercept[j][k] += gainIntercept;
	  arrayOfInterceptSquare[j][k] += gainIntercept*gainIntercept; 
	  
	  fff = (j*80)+k;
	  
	  double the_gain_sq = 0.;
	  double the_gain = 0.;
	  int dmb_id =0;
	  double the_intercept=0.0;
	  
	  the_gain = arrayOfGain[j][k];
	  the_gain_sq = arrayOfGainSquare[j][k];
	  the_intercept = arrayOfIntercept[j][k];
	  newIntercept[fff] = the_intercept;
	  newGain[fff] = the_gain;
	  
	}
      }
    }
  }

  for(int i=0;i<6;i++){
    for(int j=0;j<80;j++){
      fff = (i*80)+j;
      //std::cout <<"Layer:   "<<i<<" Strip:   "<<fff<<"  gainSlope:    "<<newGain[fff] <<"   gainIntercept:    "<<newIntercept[fff] <<std::endl;
    }
  }
  

  
  //get chamber ID from Igor's mapping
  cout<<"Here is crate and DMB: "<<crateID<<"  "<<dmbID<<endl;
  map->crate_chamber(crateID,dmbID,&chamber_id,&chamber_num,&sector);
  cout<<"This is from mapping: "<< chamber_id<<"  "<<chamber_num<<"  "<<sector<<endl;

  //info needed for database
  string test1="CSC_slice";
  string test2="gain_slope";
  string test3="gain_intercept";

  //*******to send this array to DB uncomment the next two lines*************
  //cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, newGain,1, &ret_code);
  //cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, newIntercept,1, &ret_code);
  
  //root ntuple end
  calibfile->Write();   
  calibfile->Close();
  
}// main



