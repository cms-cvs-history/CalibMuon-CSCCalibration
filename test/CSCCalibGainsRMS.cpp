/*
  authors: Stan Durkin, Oana Boeriu, Nicole Ippolito

 This calculates gains for CFEB, makes root ntuple for debugging, array of 480 entries,sends it to DB.
 Every strip pulsed 20 times,each time increase charge amplitude by 0.2 (max amplitude charge=2.2),
 each step 0.2 corresponds to 22.4fC.
 First 200 events give strip=1,second 200 events give strip=2 and so on.
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
#define LAYERS 6
#define STRIPS  80
#define CHAMBERS 5
#define NUMBERPLOTTED 10


class TCalibEvt { public:
  Float_t max_ADC[CHAMBERS][LAYERS][STRIPS];
  Float_t charge[CHAMBERS][LAYERS][STRIPS];
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
  calibevt = calibtree->Branch("EVENT", &calib_evt, "max_ADC[1][6][80]/F:charge[1][6][80]/F");

  //data file: nr.of events,chambers read
  int maxEvents = 500000;
  int misMatch = 0;
  std::string datafile=argv[1];
  int dmbID[CHAMBERS],crateID[CHAMBERS],chamber_num,sector;
  int reportedChambers =0;
  int fff, strip=-999;
  //int ret_code=-999;  
  std::string chamber_id;

  //needed for database and mapping
  //condbc *cdb = new condbc(); 
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
 
  const float x[10] = {22.4, 44.8, 67.2, 89.6, 112, 134.4, 156.8, 179.2, 201.6, 224.0};// 246.4, 268.8, 291.2, 313.6, 336.0, 358.4, 380.8, 403.2, 425.6, 448};
 
  //definition of arrays
  double adcMax[CHAMBERS][LAYERS][STRIPS];
  double adcMean_max[CHAMBERS][LAYERS][STRIPS];
  double arrayOfGain[CHAMBERS][LAYERS][STRIPS];
  double arrayOfGainSquare[CHAMBERS][LAYERS][STRIPS];
  double arrayOfIntercept[CHAMBERS][LAYERS][STRIPS];
  double arrayOfInterceptSquare[CHAMBERS][LAYERS][STRIPS];
  double newGain[480];
  double newIntercept[480];
  float maxmodten[NUMMODTEN][CHAMBERS][LAYERS][STRIPS];


 //initialize arrays
  for (int i=0; i<NUMMODTEN; i++){
    for (int j=0; j<CHAMBERS; j++){
      for (int k=0; k<LAYERS; k++){
	for (int l=0;l<STRIPS;l++){
	  maxmodten[i][j][k][l] = -999.;
	}
      }
    }
  }

  for (int i=0; i<CHAMBERS; i++){
    for (int j=0; j<LAYERS; j++){
      for (int k=0;k<STRIPS;k++){
	arrayOfGain[i][j][k]       = -999.;
	arrayOfGainSquare[i][j][k] = -999.;
	arrayOfGain[i][j][k]       = -999.;
	arrayOfIntercept[i][j][k]  = -999.;
	arrayOfInterceptSquare[i][j][k]=-999.;
	adcMax[i][j][k]            = -999.;
	adcMean_max[i][j][k]       = -999.;
      }
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
      
    std::cout << "---------- Event: " << event << "--------------" << std::endl;printf(" length %d \n",length);
    try {
      length= ddu.next(dduBuf);    
    } catch (std::runtime_error err ){
      std::cout <<"Calibration:: " << err.what()<<"  end of file?" << std::endl;
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
    

    for (int i_chamber=0; i_chamber<NChambers; i_chamber++) {//loop over all DMBs   
   
      for(int i_layer = 1; i_layer <=6; ++i_layer) {//loop over all layers in chambers
	
	std::vector<CSCStripDigi> digis = cscData[i_chamber].stripDigis(i_layer) ;
	const CSCDMBHeader &thisDMBheader = cscData[i_chamber].dmbHeader();
        
	if (thisDMBheader.cfebAvailable()){
	  dmbID[i_chamber] = cscData[i_chamber].dmbHeader().dmbID();//get DMB ID
	  crateID[i_chamber] = cscData[i_chamber].dmbHeader().crateID();//get crate ID
	  if(crateID[i_chamber] == 255) continue;
	  
	  for (unsigned int i=0; i<digis.size(); i++){//loop over digis
	    std::vector<int> adc = digis[i].getADCCounts();
	    strip = digis[i].getStrip();
            adcMax[i_chamber][i_layer-1][strip-1]=-99.;  
	    for(unsigned int k=0;k<adc.size();k++){//loop over timeBins
              float ped=(adc[0]+adc[1])/2.;
	      if(adc[k]-ped > adcMax[i_chamber][i_layer-1][strip-1]) {
		adcMax[i_chamber][i_layer-1][strip-1]=adc[k]-ped;
	      }
	    }//adc.size
	    adcMean_max[i_chamber][i_layer-1][strip-1] += adcMax[i_chamber][i_layer-1][strip-1]/20.;  
            
	    // // On the 10th event
	    if (evt%20 == 0&&(strip-1)%16==(evt-1)/200){
	      int ten = int((evt-1)/20)%10 ;
	      maxmodten[ten][i_chamber][i_layer-1][strip-1] = adcMean_max[i_chamber][i_layer-1][strip-1];
	    }
	  }//end digis loop
	}//end cfeb.available loop
      }//end layer loop
    }//end chamber loop
    
    if((evt-1)%20==0){
      for(int ii=0;ii<CHAMBERS;ii++){
	for(int jj=0;jj<LAYERS;jj++){
	  for(int kk=0;kk<STRIPS;kk++){
	    adcMean_max[ii][jj][kk]=0.0;
	  }
	}
      }
    }
    
    // Fill the ntuple every 20th event
    if (evt%20 == 0){
      for (int thechamber = 0; thechamber<CHAMBERS; thechamber++){
	for(int thelayer = 0; thelayer<LAYERS; thelayer++) {
	  for (int thestrip=0; thestrip<STRIPS; thestrip++){
	    calib_evt.max_ADC[thechamber][thelayer][thestrip] = maxmodten[int(evt/20) - 1][thechamber][thelayer][thestrip];
	    calib_evt.charge[thechamber][thelayer][thestrip]=x[int(evt/20) - 1];
	  }
	}
      }
      calibtree->Fill(); 
    }
  }//end event loop

  
  //create array (480 entries) for database transfer
  for(int chamberiter=0;chamberiter<1;chamberiter++){
    double the_gain_sq = 0.;
    double the_gain = 0.;
    double the_intercept=0.0; 
    
    
    for (int thischamber=0; thischamber<CHAMBERS;thischamber++){
      for (int thislayer=0; thislayer<LAYERS;thislayer++){
	for (int thisstrip=0; thisstrip<STRIPS;thisstrip++){
	  arrayOfGain[thischamber][thislayer][thisstrip]       = 0.;
	  arrayOfGainSquare[thischamber][thislayer][thisstrip] = 0.;
	  arrayOfIntercept[thischamber][thislayer][thisstrip]  = 0.;
	  arrayOfInterceptSquare[thischamber][thislayer][thisstrip] = 0.;
	}
      }
    }
        
    for (int cham=0;cham<1;cham++){
      if (cham !=chamberiter) continue;
            
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
	      
	      for(int ii=0;ii<10;ii++){//numbers       
		//do straight line fit ( y = kx + m )
		sumOfX += x[ii];
		sumOfY += maxmodten[ii][cham][j][k];
		sumOfXY += (x[ii]*maxmodten[ii][cham][j][k]);
		sumx2 += (x[ii]*x[ii]);
		
		gainSlope= ((NUMBERPLOTTED*sumOfXY) - (sumOfX * sumOfY))/((NUMBERPLOTTED*sumx2) - (sumOfX*sumOfX));//k
		gainIntercept = ((sumOfY*sumx2)-(sumOfX*sumOfXY))/((NUMBERPLOTTED*sumx2)-(sumOfX*sumOfX));//m
		chi2  += (maxmodten[ii][cham][j][k]-(gainIntercept+(gainSlope*x[ii])))*(maxmodten[ii][cham][j][k]-(gainIntercept+(gainSlope*x[ii])));
		
		if(chi2<min1){
		  min1=chi2;
		  gainSlope=gainSlope;
		  gainIntercept=gainIntercept;
		}
	      }
	      
	      arrayOfGain[cham][j][k]       += gainSlope;
	      arrayOfGainSquare[cham][j][k] += gainSlope*gainSlope;
	      arrayOfIntercept[cham][j][k]  += gainIntercept;
	      arrayOfInterceptSquare[cham][j][k] += gainIntercept*gainIntercept; 
	      
	      fff = (j*80)+k; //this is for 480 entries in the array; obsolite soon!
	    
	      the_gain          = arrayOfGain[cham][j][k];
	      the_gain_sq       = arrayOfGainSquare[cham][j][k];
	      the_intercept     = arrayOfIntercept[cham][j][k];
	      newIntercept[fff] = the_intercept;
	      newGain[fff]      = the_gain;
	      
	    }//k loop
	  }//j loop
	}//stripiter loop
      }//layiter loop
    }//cham loop
  }//chamberiter loop

  for(int i=0;i<CHAMBERS;i++){
    for(int j=0;j<LAYERS;j++){
      for (int k=0;k<STRIPS;k++){
	fff = (j*80)+k; //this is for 480 entries in the array; obsolite soon!
	//std::cout <<"Chamber: "<<i<<" Layer:   "<<i<<" Strip:   "<<fff<<"  Slope:    "<<newGain[fff] <<"  Intercept:    "<<newIntercept[fff] <<std::endl;
      }
    }
    //get chamber ID from Igor's mapping        
    int new_crateID = crateID[i];
    int new_dmbID   = dmbID[i];
    std::cout<<"Here is crate: "<<new_crateID<<" and DMB:  "<<new_dmbID<<std::endl;
    map->crate_chamber(new_crateID,new_dmbID,&chamber_id,&chamber_num,&sector);
    std::cout<<"This is from mapping: "<< chamber_id<<"  "<<chamber_num<<"  "<<sector<<std::endl;
    //info needed for database
    string test1="CSC_slice";
    string test2="gain_slope";
    string test3="gain_intercept";
    
    //*******to send this array to DB uncomment the next two lines*************
    //cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, newGain,1, &ret_code);
    //cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, newIntercept,1, &ret_code);
  }
    
  //root ntuple end
  calibfile->Write();   
  calibfile->Close();
  
}// main



