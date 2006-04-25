/*
  authors: Oana Boeriu, Nicole Ippolito
  
  This will create arrays (480=whole ME2/3 chamber) of comparator threshold and noise values
  which can be sent to online database.
  It creates a root ntuple for debugging purpose.  
*/

#include "IORawData/CSCCommissioning/src/FileReaderDDU.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDDUEventData.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigi.h"
#include "EventFilter/CSCRawToDigi/interface/CSCEventData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCCLCTData.h"
#include "DataFormats/CSCDigi/interface/CSCComparatorDigi.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDMBHeader.h"
#include "CalibMuon/CSCCalibration/interface/condbc.h"
#include "CalibMuon/CSCCalibration/interface/cscmap.h" 
#include "FWCore/MessageService/interface/MessageServicePresence.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "string"
#include <math.h>

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
#define NUMMOD 875
#define NUMBERPLOTTED 25

class TCalibEvt { public:
  Float_t mean[CHAMBERS][LAYERS][STRIPS];
  Float_t charge1[CHAMBERS][LAYERS][STRIPS];
};

//void dierfc(double y);

int main(int argc, char **argv) {
  edm::service::MessageServicePresence my_message_service;
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

  calibfile = new TFile("ntuples/calibcomparator.root","RECREATE","Calibration Ntuple");
  calibtree = new TTree("Calibration","Comparator threshold");
  calibevt = calibtree->Branch("EVENT", &calib_evt, "mean[5][6][80]/F:charge1[5][6][80]/F");
  
  //data file: nr.of events,chambers read
  int maxEvents = 900000;
  int misMatch = 0;
  std::string datafile=argv[1];
  std::string chamber_id;  
  int dmbID[CHAMBERS],crateID[CHAMBERS],chamber_num,sector;
  int i_chamber=0,i_layer=0;
  int reportedChambers =0;
  int fff; 
  
  //needed for database and mapping
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();
  
  if (argv[2]) maxEvents = (int) atof(argv[2]);
  
  FileReaderDDU ddu;
  ddu.open(datafile.c_str());
  
  //some variable declaration
  int evt=0,ret_code;
    
  //charge injected for different thresholds: increase 35 steps of 3mV each and start at 13mV 
 
  const float thresh1[35] = {13.0, 16.0, 19.0, 22.0, 25.0, 28.0, 31.0, 34.0, 37.0, 40.0, 43.0, 46.0, 49.0, 52.0, 55.0, 58.0, 61.0, 64.0, 67.0, 70.0, 73.0, 76.0, 79.0, 82.0, 85.0, 88.0, 91.0, 94.0, 97.0, 100.0, 103.0, 106.0, 109.0, 112.0, 115.0}; 
    
  //definition of arrays
  double mean[CHAMBERS][LAYERS][STRIPS];
  double meanTot[CHAMBERS][LAYERS][STRIPS];
  double theMeanThresh[CHAMBERS][LAYERS][STRIPS];
  float meanmod[NUMMOD][CHAMBERS][LAYERS][STRIPS];
  double arrayMeanThresh[CHAMBERS][LAYERS][STRIPS];
  double newThresh[480];
  //double newRMS[480];
  
  //initialize arrays
  for(int i=0;i<CHAMBERS;i++){
    for(int j=0; j<LAYERS; j++){
      for(int k=0; k<STRIPS; k++){
	theMeanThresh[i][j][k] = 0.;
	arrayMeanThresh[i][j][k] = 0.;
	mean[i][j][k]=0.;
	meanTot[i][j][k]=0.;
      }
    }
  }

  for (int i=0;i<480;i++){
    newThresh[i]=0;
    // newRMS[i]=0;
  }
  
  //opening data files and checking data (copied from Alex Tumanov's code)
  const unsigned short *dduBuf=0;
  int length = 1;
  int timebin=-999,mycompstrip=-999;

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
      if (cscData[i_chamber].nclct()) {
 	CSCCLCTData & clctData = cscData[i_chamber].clctData();
      }else {
	std::cout<<" No CLCT!" <<std::endl;
	continue;
      }
      CSCCLCTData & clctData = cscData[i_chamber].clctData();
      for(i_layer = 1; i_layer <= 6; ++i_layer) {//loop over all layers in chambers
	//std::vector<CSCStripDigi> digis = cscData[i_chamber].stripDigis(i_layer) ;
	std::vector<CSCComparatorDigi> comp = clctData.comparatorDigis(i_layer);

	const CSCDMBHeader &thisDMBheader = cscData[i_chamber].dmbHeader();
	
	if (thisDMBheader.cfebAvailable()){//check that CFEB data exists

	  dmbID[i_chamber]   = cscData[i_chamber].dmbHeader().dmbID(); //get DMB ID
	  crateID[i_chamber] = cscData[i_chamber].dmbHeader().crateID(); //get crate ID
	  if(crateID[i_chamber] == 255) continue; //255 is reserved for old crate, present only 0 and 1
	  
	  for (unsigned int i=0; i<comp.size(); i++){//loop over CFEB comparator digis
	    int comparator = comp[i].getComparator();
	    timebin = comp[i].getTimeBin() ;
	    int compstrip =  comp[i].getStrip();

	    int this_comparator[4] = {4, 5, 6, 7};

	    for (int iii=0; iii<40; iii++){
	      if ((compstrip == iii) && (comparator == this_comparator[0] || comparator == this_comparator[1])) {
		mycompstrip = 0 + iii*2;
	      } else if ((compstrip == iii) && (comparator == this_comparator[2] || comparator == this_comparator[3])) {
		mycompstrip = 1 + iii*2;
	      }
	    }

	    mean[i_chamber][i_layer-1][mycompstrip] = comparator/5;
	    
	  }//end comp loop
	  
	  meanTot[i_chamber][i_layer-1][mycompstrip] +=mean[i_chamber][i_layer-1][mycompstrip]/25.;
	  std::cout<<"layer: "<<i_layer<<"  strip  "<<mycompstrip<<"   MeanTot "<<meanTot[i_chamber][i_layer-1][mycompstrip]<<std::endl;

	  // On the 25th event
	  if (evt%25 == 0&&(mycompstrip)%16==(evt-1)/875){
	    int tmp = int((evt-1)/25)%35 ;
	    meanmod[tmp][i_chamber][i_layer-1][mycompstrip] = meanTot[i_chamber][i_layer-1][mycompstrip];
	    std::cout<<"Mean "<<meanmod[tmp][i_chamber][i_layer-1][mycompstrip]<<" chamber  "<<i_chamber<<"  layer  "<<i_layer<<"   mycompstrip "<<mycompstrip<<" charge "<<thresh1[int(evt/25) - 1]<<std::endl;
	  }
	}//end if cfeb.available loop
      }//end layer loop
    }//end chamber loop

    if((evt-1)%25==0){
      for(int ii=0;ii<CHAMBERS;ii++){
      	for(int jj=0;jj<LAYERS;jj++){
      	  for(int kk=0;kk<STRIPS;kk++){
      	    mean[ii][jj][kk]=0.0;
	    meanTot[ii][jj][kk]=0.0;
      	  }
      	}
      }
    }
    // Fill the ntuple every 25th event
    if (evt%25 == 0){
      for (int thechamber = 0; thechamber<CHAMBERS; thechamber++){
	for(int thelayer = 0; thelayer<LAYERS; thelayer++) {
	  for (int thestrip=0; thestrip<STRIPS; thestrip++){
	    calib_evt.mean[thechamber][thelayer][thestrip] = meanmod[int(evt/25) - 1][thechamber][thelayer][thestrip];
	    calib_evt.charge1[thechamber][thelayer][thestrip]=thresh1[int(evt/25) - 1];
	  }
	}
	calibtree->Fill(); 
      }
    }
    
  }//end events loop
  
  //root ntuple end
  calibfile->Write();   
  calibfile->Close();

  //create array (480 entries) for database transfer
  for(int myChamber=0; myChamber<CHAMBERS; myChamber++){
    double meanComp = 0.;
    std::string test1="CSC_slice";
    std::string test2="comparator_threshold";
    std::string answer;
    
    /////////Call erf function and do the fit here!!!/////////////////


    //print out result here
    for (int i=0; i<CHAMBERS; i++){
      if (myChamber !=i) continue;
      
      for (int j=0; j<LAYERS; j++){
	for (int k=0; k<STRIPS; k++){
	  //arrayMeanThresh[i][j][k]= meanmod[tmp][i_chamber][i_layer-1][mycompstrip];
	  fff = (j*80)+k;
	  //theMeanThresh  = arrayMeanThresh[i][j][k];
	  //newMeanThresh[fff]  = theRMS;
	  
	  //std::cout <<" chamber "<<i<<" layer "<<j<<" strip "<<fff<<"  ped "<<newPed[fff]<<" RMS "<<newRMS[fff]<<std::endl;
	}
      }
    }
    //get chamber ID from Igor's mapping

    int new_crateID = crateID[myChamber];
    int new_dmbID   = dmbID[myChamber];
    std::cout<<" Crate: "<<new_crateID<<" and DMB:  "<<new_dmbID<<std::endl;
    map->crate_chamber(new_crateID,new_dmbID,&chamber_id,&chamber_num,&sector);
    std::cout<<" Above data is for chamber: "<< chamber_id<<" from sector "<<sector<<std::endl;

   //  std::cout<<" DO you want to send constants to DB? "<<std::endl;
//     std::cout<<" Please answer y or n for EACH chamber present! "<<std::endl;

//     std::cin>>answer;
//     if(answer=="y"){
//       //SEND CONSTANTS TO DB
//       cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, newMeanThresh,2, &ret_code);
//       std::cout<<" Your results were sent to DB !!! "<<std::endl;
//     }else{
//       std::cout<<" NO data was sent!!! "<<std::endl;
//     }
  }
  
  
}//main

/* inverse of error function in double precision */
double dierfc(double y){
  double s, t, u, w, x, z;
  
  z = y;
  if (y > 1) {
    z = 2 - y;
  }
  w = 0.916461398268964 - log(z);
  u = sqrt(w);
  s = (log(u) + 0.488826640273108) / w;
  t = 1 / (u + 0.231729200323405);
  x = u * (1 - s * (s * 0.124610454613712 + 0.5)) - 
    ((((-0.0728846765585675 * t + 0.269999308670029) * t + 
       0.150689047360223) * t + 0.116065025341614) * t + 
     0.499999303439796) * t;
  t = 3.97886080735226 / (x + 3.97886080735226);
  u = t - 0.5;
  s = (((((((((0.00112648096188977922 * u + 
	       1.05739299623423047e-4) * u - 0.00351287146129100025) * u - 
	     7.71708358954120939e-4) * u + 0.00685649426074558612) * u + 
	   0.00339721910367775861) * u - 0.011274916933250487) * u - 
	 0.0118598117047771104) * u + 0.0142961988697898018) * u + 
       0.0346494207789099922) * u + 0.00220995927012179067;
  s = ((((((((((((s * u - 0.0743424357241784861) * u - 
		 0.105872177941595488) * u + 0.0147297938331485121) * u + 
	       0.316847638520135944) * u + 0.713657635868730364) * u + 
	     1.05375024970847138) * u + 1.21448730779995237) * u + 
	   1.16374581931560831) * u + 0.956464974744799006) * u + 
	 0.686265948274097816) * u + 0.434397492331430115) * u + 
       0.244044510593190935) * t - 
    z * exp(x * x - 0.120782237635245222);
  x += s * (x * s + 1);
  if (y > 1) {
    x = -x;
  }
  return x;
}
