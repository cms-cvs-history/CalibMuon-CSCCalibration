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
#define CHAMBERS 5
#define LAYERS 6
#define STRIPS 80
#define TIMEBINS 8


class TCalibEvt { public:
  Int_t adc[8];
    Float_t pedMean;
    Int_t strip;
    Float_t time[8];
    Int_t chamber;
    Int_t event;
    Int_t layer;
};


//functions
float elec(float t,float vs);
void mkbins(float vs);
void convolution(float *xleft_a, float *xleft_b, float *min_left, float* xright_a, float* xright_b, float *min_right);
float elec(float time,float vs);
float chifit_ts(float tql,float tq,float tqh);

float convd[3][120];
float nconvd[3][120];
float conve[120];
float conv[3][120];

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

  calibfile = new TFile("calibcrosstalk.root","RECREATE","Calibration Ntuple");
  calibtree = new TTree("Calibration","Pedestal");
  calibevt = calibtree->Branch("EVENT", &calib_evt, "adc[8]/I:pedMean/F:strip/I:time[8]/F:chamber/I:event/I:layer/I");

  //variables declaration
  int maxEvents = 500000;
  int misMatch = 0;
  int dmbID=-999,crateID=-999,chamber_num,sector;
  int reportedChambers =0;
  int fff,ret_code,strip=-999;  
  std::string chamber_id;
  std::string datafile=argv[1];
  float pedMean=0.0,time=0.0;
  
  //needed for database and mapping
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();
  
  //declare arrays
  float thetime[CHAMBERS][LAYERS][STRIPS][TIMEBINS];
  int theadccounts[CHAMBERS][LAYERS][STRIPS][TIMEBINS];
  float crosstalk_left[CHAMBERS][LAYERS][STRIPS];
  float crosstalk_right[CHAMBERS][LAYERS][STRIPS];
  
  //initialize arrays  
  for (int i=0; i<CHAMBERS; i++){
    for (int j=0; j<LAYERS; j++){
      for (int k=0; k<STRIPS; k++){
	for (int l=0; l<TIMEBINS; l++){
	  thetime[i][j][k][l] = 0.;
	  theadccounts[i][j][k][l] = 0;
	}
      }
    }
  }
  
  for (int i=0; i<CHAMBERS; i++){
    for (int j=0; j<LAYERS; j++){
      for (int k=0; k<STRIPS; k++){
	crosstalk_left[i][j][k] = -999.;
	crosstalk_right[i][j][k] = -999.;
      }
    }
  }
  
  // open data file
  if (argv[2]) maxEvents = (int) atof(argv[2]);
  FileReaderDDU ddu;
  ddu.open(datafile.c_str());
  
  const unsigned short *dduBuf=0;
  int length = 1;

  //////////////////////////////////////////////////////////////////////////////////////
  // First loop over all events and fill our arrays with time and ADC information
  
  for (int event = 0; (event < maxEvents) && length; ++event) {
  
    std::cout << "---------------------- Event "<< event<< " --------------------- " <<std::endl;
    try {
      length = ddu.next(dduBuf);
    } catch (std::runtime_error err){
      std::cout  << "digi Anal:: " << err.what()<< " end of file?" << std::endl;
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

    // Loop over chambers, layers and strips
    for (int chamber = 0; chamber < NChambers; chamber++){
   
      for (int layer = 1; layer <= 6; layer++){
	
	std::vector<CSCStripDigi> digis = cscData[chamber].stripDigis(layer) ;
	const CSCDMBHeader &thisDMBheader = cscData[chamber].dmbHeader();
      
	//for (int strip = 0; strip < 80; strip++){

	for (unsigned int i=0; i<digis.size(); i++){
	  int strip = digis[i].getStrip();
	  std::vector<int> adc = digis[i].getADCCounts();
	  //int which_event = event / 20;

	  //std::cout << "event / 20 = " << which_event << "  strip = " << strip << endl;
	  // if ((event / 20 != 0) && ((strip == 1) || (strip == 17) || (strip == 33) || (strip == 49) || (strip == 65))) continue;
// 	  if ((event / 20 != 1) && ((strip == 2) || (strip == 18) || (strip == 34) || (strip == 50) || (strip == 66))) continue;
// 	  if ((event / 20 != 2) && ((strip == 3) || (strip == 19) || (strip == 35) || (strip == 51) || (strip == 67))) continue;
// 	  if ((event / 20 != 3) && ((strip == 4) || (strip == 20) || (strip == 36) || (strip == 52) || (strip == 68))) continue;
// 	  if ((event / 20 != 4) && ((strip == 5) || (strip == 21) || (strip == 37) || (strip == 53) || (strip == 69))) continue;
// 	  if ((event / 20 != 5) && ((strip == 6) || (strip == 22) || (strip == 38) || (strip == 54) || (strip == 70))) continue;
// 	  if ((event / 20 != 6) && ((strip == 7) || (strip == 23) || (strip == 39) || (strip == 55) || (strip == 71))) continue;
// 	  if ((event / 20 != 7) && ((strip == 8) || (strip == 24) || (strip == 40) || (strip == 56) || (strip == 72))) continue;
// 	  if ((event / 20 != 8) && ((strip == 9) || (strip == 25) || (strip == 41) || (strip == 57) || (strip == 73))) continue;
// 	  if ((event / 20 != 9) && ((strip == 10) || (strip == 26) || (strip == 42) || (strip == 58) || (strip == 74))) continue;
// 	  if ((event / 20 != 10) && ((strip == 11) || (strip == 27) || (strip == 43) || (strip == 59) || (strip == 75))) continue;
// 	  if ((event / 20 != 11) && ((strip == 12) || (strip == 28) || (strip == 44) || (strip == 60) || (strip == 76))) continue;
// 	  if ((event / 20 != 12) && ((strip == 13) || (strip == 29) || (strip == 45) || (strip == 61) || (strip == 77))) continue;
// 	  if ((event / 20 != 13) && ((strip == 14) || (strip == 30) || (strip == 46) || (strip == 62) || (strip == 78))) continue;
// 	  if ((event / 20 != 14) && ((strip == 15) || (strip == 31) || (strip == 47) || (strip == 63) || (strip == 79))) continue;
// 	  if ((event / 20 != 15) && ((strip == 16) || (strip == 32) || (strip == 48) || (strip == 64) || (strip == 80))) continue;

	  int offset = event / 20;
	  
	  for(unsigned int k=0;k<adc.size();k++){
	    
	    time = (50. * k)-(event * 6.25)+116.5+(200*offset);
	    pedMean =(adc[0]+adc[1])/2;
	
	    calib_evt.event    = event;
	    calib_evt.pedMean  = pedMean;
	    calib_evt.strip    = strip;
	    calib_evt.chamber  = chamber;
	    calib_evt.layer    = layer;
	    calib_evt.adc[k]   = adc[k];
	    calib_evt.time[k]  = time;

	    thetime[chamber][layer-1][strip-1][k] = time;
	    theadccounts[chamber][layer-1][strip-1][k] = adc[k]; 
	  } //end loop over timebins

	  calibtree->Fill();

	}//end loop over digis
      }//end loop over layers 
    }//end loop over chambers
  }//end loop over events

  
  ////////////////////////////////////////////////////////////////////
  // Now that we have filled our array, extract convd and nconvd
  float adc_ped_sub_left = -999.;
  float adc_ped_sub = -999.;
  float adc_ped_sub_right = -999.;
  int thebin;

  for (int i=0; i<CHAMBERS; i++){
    for (int j=0; j<LAYERS; j++){
      for (int k=0; k<STRIPS; k++){

	if (k == 0 || k == (STRIPS - 1)) continue;
	if (k == 80 || k == (STRIPS + 1)) continue;

	// re-zero convd and nconvd
	for (int m=0; m<3; m++){
	  for (int n=0; n<120; n++){
	    convd[m][n] = 0.;
	    nconvd[m][n] = 0.;
	  }
	}
	
	// Need special cases for k=0, k=N
	for (int l=0; l<TIMEBINS; l++){
	  adc_ped_sub_left  = theadccounts[i][j][k-1][l] - (theadccounts[i][j][k-1][0] + theadccounts[i][j][k-1][1])/2.;	  
	  adc_ped_sub       = theadccounts[i][j][k][l]   - (theadccounts[i][j][k][0]   + theadccounts[i][j][k][1])/2.;
	  adc_ped_sub_right = theadccounts[i][j][k+1][l] - (theadccounts[i][j][k+1][0] + theadccounts[i][j][k+1][1])/2.;
	  thebin = ((thetime[i][j][k][l] - 200. + 0.001)/6.24);

	  if (thebin >= 0 && thebin < 120){
	    convd[k-1][thebin] += adc_ped_sub_left;
	    nconvd[k-1][thebin] += 1.0;
	    
	    convd[k][thebin] += adc_ped_sub;
	    nconvd[k][thebin] += 1.0;

	    convd[k+1][thebin] += adc_ped_sub_right;
	    nconvd[k+1][thebin] += 1.0;

	  }
	} //loop over timebins

	// Call our functions to calculate the cross-talk
	float xl_temp_a = 0.;
	float xl_temp_b = 0.;
	float minl_temp = 0.;
	float xr_temp_a = 0.;
	float xr_temp_b = 0.;
	float minr_temp = 0.;
	mkbins(50.);
	convolution(&xl_temp_a, &xl_temp_b, &minl_temp, &xr_temp_a, &xr_temp_b, &minr_temp);

	crosstalk_left[i][j][k] = xl_temp_a;
        crosstalk_right[i][j][k] = xr_temp_a;

      }//loop over strips
    }//loop over layers
  }//loop over chambers

  calibfile->Write();   
  calibfile->Close();

 //get chamber ID from Igor's mapping
  cout<<"Here is crate and DMB: "<<crateID<<"  "<<dmbID<<endl;
  map->crate_chamber(crateID,dmbID,&chamber_id,&chamber_num,&sector);
  cout<<"This is from mapping: "<< chamber_id<<"  "<<chamber_num<<"  "<<sector<<endl;

  //info needed for database
  string test1="CSC_slice";
  string test2="xtalk_slope_right";
  string test3="xtalk_intercept_right";
  string test4="xtalk_chi2_right";
  string test5="xtalk_slope_left";
  string test6="xtalk_intercept_left";
  string test7="xtalk_chi2_left";

  //*******to send this array to DB uncomment the next six lines*************
  //cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, new_slope_right,    1, &ret_code);
  //cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, new_intercept_right,1, &ret_code);
  //cdb->cdb_write(test1,chamber_id,chamber_num,test4,480, new_chi2_right,     1, &ret_code);
  //cdb->cdb_write(test1,chamber_id,chamber_num,test5,480, new_slope_left,     1, &ret_code);
  //cdb->cdb_write(test1,chamber_id,chamber_num,test6,480, new_intercept_left, 1, &ret_code);
  //cdb->cdb_write(test1,chamber_id,chamber_num,test7,480, new_chi2_left     , 1, &ret_code);

  //root ntuple end
  calibfile->Write();   
  calibfile->Close();

}// main


//square wave fcn convoluted with 1/(t+2.5)
float elec(float t,float vs){
  float f;
  if (t<=vs){
    f=log(t+2.5)-log(2.5);
  }
  else{
    f=log(t+2.5)-log(t-50+2.5);
  }
  return f;
}//elec

//calculate single electron distribution in 6.25 ns steps
void mkbins(float vs){
  int i,k;
  float t;
  for(i=0;i<120;i++) conve[i]=0.0;
  for(i=0;i<120;i++){
    for(k=0;k<16;k++){
      t=(6.25*i)+(k*0.625);
      conve[i]=conve[i]+elec(t,vs);
    }
  }
} //mkbins

//convolution function
void convolution(float *xleft_a, float *xleft_b, float *min_left, float *xright_a, float *xright_b, float *min_right){
  //void(convolution){  
  float max,cross0,cross2,min_l,min_r,sum_x=0.0,sumx2=0.;
  float sum_y_left=0.0,sum_y_right=0.0,sum_xy_left=0.0,sum_xy_right=0.0;
  float a_left=0.0,a_right=0.0,b_left=0.0,b_right=0.0,chi2_left=0.0,chi2_right=0.0,chi_left=0.0,chi_right=0.0;
  int i,j,k,l;

  for(l=0;l<3;l++){
    for(i=0;i<119;i++)conv[l][i]=0.0;
    for(j=0;j<119;j++){
      for(k=0;k<119;k++){
	if(j+k<119)conv[l][j+k]=conv[l][j+k]+convd[l][j]*conve[k];
      }
    }
  }
  max=0;
  min_l=9999999999999999.0;
  min_r=9999999999999999.0;
  for(i=0;i<119;i++){
    if(conv[1][i]>max) 
      max=conv[1][i];
  }
  
  for(l=0;l<3;l++){
    for(i=0;i<119;i++)conv[l][i]=conv[l][i]/max;
  }

  int nobs = 0;
  for (int j=0; j<119; j++){
    if (conv[1][j]>0.2) nobs++;
  }

  for(i=0;i<119;i++){
    cross0=0.0;
    cross2=0.0;

    if(conv[1][i]>0.2){
      cross0=conv[0][i]/(conv[0][i]+conv[1][i]+conv[2][i]);
      cross2=conv[2][i]/(conv[0][i]+conv[1][i]+conv[2][i]);
      
      sum_x += i;
      sum_y_left += cross0;
      sum_y_right += cross2;
      sumx2 += i*i;
      sum_xy_left += i*cross0;
      sum_xy_right += i*cross2;
      
      //LMS fitting straight line
            b_left=((nobs*sum_xy_left) - (sum_x * sum_y_left))/((nobs*sumx2) - (sum_x*sum_x));
            b_right=((nobs*sum_xy_right) - (sum_x * sum_y_right))/((nobs*sumx2) - (sum_x*sum_x));
            a_left=((sum_y_left*sumx2)-(sum_x*sum_xy_left))/((nobs*sumx2)-(sum_x*sum_x));
            a_right=((sum_y_right*sumx2)-(sum_x*sum_xy_right))/((nobs*sumx2)-(sum_x*sum_x));
            chi2_left += (cross0 -(a_left+(b_left*i)))*(cross0 -(a_left+(b_left*i)));
            chi_left +=sqrt(chi2_left);
            chi2_right += (cross2 -(a_right+(b_right*i)))*(cross2 -(a_right+(b_right*i)));
            chi_right =sqrt(chi2_right);	      
      if(chi_left<min_l){
	min_l=chi_left;
	b_left=b_left;
	a_left=a_left;
      }
      if(chi_right<min_r){
	min_r=chi_right;
	b_right=b_right;
	a_right=a_right;
      }
      
      *xleft_a = a_left; 
      *xleft_b = b_left;
      *min_left = min_l;
      *xright_a = a_right;
      *xright_b = b_right;
      *min_right = min_r;
    }  
   
    void close();
  }//convolution
}



