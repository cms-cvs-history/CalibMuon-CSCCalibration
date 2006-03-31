/*
  authors: Stan Durkin, Oana Boeriu, Nicole Ippolito
  
  This will create arrays (480=whole ME2/3 chamber) of fit parameters for cross-talk
  which can be sent to online database.
  It creates a root ntuple for debugging purposes.   
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
#define CHAMBERS 4
#define LAYERS 6
#define STRIPS 80
#define TIMEBINS 8
using namespace std;

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
  
  calibfile = new TFile("ntuples/calibcrosstalk.root","RECREATE","Calibration Ntuple");
  calibtree = new TTree("Calibration","Pedestal");
  calibevt = calibtree->Branch("EVENT", &calib_evt, "adc[8]/I:pedMean/F:strip/I:time[8]/F:chamber/I:event/I:layer/I");
  
    
  //variables declaration
  int maxEvents = 500000;
  int misMatch = 0;
  int dmbID[CHAMBERS],crateID[CHAMBERS],chamber_num,sector;
  int reportedChambers =0;
  int fff,ret_code;  
  std::string chamber_id;
  std::string datafile=argv[1];
  float pedMean=0.0,time=0.0;
  
  //needed for database and mapping
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();

  float new_rxtalk[480];
  float new_lxtalk[480];
 
  //declare arrays
  float thetime[CHAMBERS][LAYERS][STRIPS][TIMEBINS*20];
  int thebins[CHAMBERS][LAYERS][STRIPS][TIMEBINS*20];
  int theadccountsc[CHAMBERS][LAYERS][STRIPS][TIMEBINS*20];
  int theadccountsl[CHAMBERS][LAYERS][STRIPS][TIMEBINS*20];
  int theadccountsr[CHAMBERS][LAYERS][STRIPS][TIMEBINS*20];
  float xtalk_intercept_left[CHAMBERS][LAYERS][STRIPS];
  float xtalk_intercept_right[CHAMBERS][LAYERS][STRIPS];
  float xtalk_slope_left[CHAMBERS][LAYERS][STRIPS];
  float xtalk_slope_right[CHAMBERS][LAYERS][STRIPS];
  float xtalk_chi2_left[CHAMBERS][LAYERS][STRIPS];
  float xtalk_chi2_right[CHAMBERS][LAYERS][STRIPS];
  double new_xtalk_intercept_right[480];
  double new_xtalk_intercept_left[480];
  double new_xtalk_slope_right[480];
  double new_xtalk_slope_left[480];
  double new_rchi2[480];
  double new_lchi2[480];


   //initialize arrays  

  for (int i=0;i<480;i++){
    new_lxtalk[i]=0.;
    new_rxtalk[i]=0.;
    new_xtalk_intercept_right[i] = -999.;
    new_xtalk_intercept_left[i]  = -999.;
    new_xtalk_slope_right[i]     = -999.;
    new_xtalk_slope_left[i]      = -999.;
    new_rchi2[i]                 = -999.;
    new_lchi2[i]                 = -999.;
  }

  for (int i=0; i<CHAMBERS; i++){
    for (int j=0; j<LAYERS; j++){
      for (int k=0; k<STRIPS; k++){
	for (int l=0; l<TIMEBINS*20; l++){
	  thetime[i][j][k][l] = 0.;
          thebins[i][j][k][l] = 0;
	  theadccountsc[i][j][k][l] = 0;
	  theadccountsl[i][j][k][l] = 0;
	  theadccountsr[i][j][k][l] = 0;
	}
      }
    }
  }
  
  for (int i=0; i<CHAMBERS; i++){
    for (int j=0; j<LAYERS; j++){
      for (int k=0; k<STRIPS; k++){
	xtalk_intercept_left[i][j][k] = -999.;
	xtalk_intercept_right[i][j][k] = -999.;
	xtalk_slope_left[i][j][k] = -999.;
	xtalk_slope_right[i][j][k] = -999.;
	xtalk_chi2_left[i][j][k] = -999.;
	xtalk_chi2_right[i][j][k] = -999.;
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
      length = ddu.next(dduBuf);printf(" length %d \n",length);
    } catch (std::runtime_error err){
      std::cout  << "Calibration:: " << err.what()<< " end of file?" << std::endl;
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
    for (int chamber = 0; chamber < NChambers; chamber++){//loop over all DMBs
   
       for (int layer = 1; layer <= 6; layer++){
	
	std::vector<CSCStripDigi> digis = cscData[chamber].stripDigis(layer) ;
	const CSCDMBHeader &thisDMBheader = cscData[chamber].dmbHeader();
	
	if (thisDMBheader.cfebAvailable()){
          dmbID[chamber] = cscData[chamber].dmbHeader().dmbID();//get DMB ID
          crateID[chamber] = cscData[chamber].dmbHeader().crateID();//get crate ID
	  if(crateID[chamber] == 255) continue;
	  
	  for (unsigned int i=0; i<digis.size(); i++){
	    int strip = digis[i].getStrip();
	    std::vector<int> adc = digis[i].getADCCounts();
	    
	    int offset = event / 20;
	    int smain[5],splus[5],sminus[5]; //5 for CFEBs
	    for(int s=0;s<5;s++) smain[s]  = s*16+offset;
	    for(int s=0;s<5;s++) splus[s]  = s*16+offset+1;
	    for(int s=0;s<5;s++) sminus[s] = s*16+offset-1;
	    int iuse=-99;
	    for(int s=0; s<5; s++) {if(strip-1==smain[s])  iuse=smain[s];}
	    for(int s=0; s<5; s++) {if(strip-1==splus[s])  iuse=smain[s];}
	    for(int s=0; s<5; s++) {if(strip-1==sminus[s]) iuse=smain[s];}
	    
	    if(iuse!=-99){
	      
	      for(unsigned int k=0;k<adc.size();k++){
		
		time = (50. * k)-((event%20)* 6.25)+116.5;
		pedMean =(adc[0]+adc[1])/2;
		
		calib_evt.event    = event;
		calib_evt.pedMean  = pedMean;
		calib_evt.strip    = strip;
		calib_evt.chamber  = chamber;
		calib_evt.layer    = layer;
		calib_evt.adc[k]   = adc[k];
		calib_evt.time[k]  = time;
		
		int kk=8*k-event%20+19;//19 to zero everything, for binning 120
		
		thebins[chamber][layer-1][strip-1][kk] = 8*k-event%20+19;
		thetime[chamber][layer-1][strip-1][kk] = time;
		
		if(iuse==strip-1)  theadccountsc[chamber][layer-1][iuse][kk] = adc[k]; 
		if(iuse==strip)    theadccountsr[chamber][layer-1][iuse][kk] = adc[k];
		if(iuse==strip-2)  theadccountsl[chamber][layer-1][iuse][kk] = adc[k]; 
		
	      } //end loop over timebins
	      
	      calibtree->Fill();
	    }
	    
	  }//end loop over digis
	}//end cfeb.available loop
       }//end loop over layers 
    }//end loop over chambers
  }//end loop over events
 

   //root ntuple end
  calibfile->Write();   
  calibfile->Close(); 
  
  ////////////////////////////////////////////////////////////////////iuse==strip-1
  // Now that we have filled our array, extract convd and nconvd
  float adc_ped_sub_left = -999.;
  float adc_ped_sub = -999.;
  float adc_ped_sub_right = -999.;
  int thebin;

   for (int i=0; i<CHAMBERS; i++){
    for (int j=0; j<LAYERS; j++){
      for (int k=0; k<STRIPS; k++){
	// re-zero convd and nconvd
	for (int m=0; m<3; m++){
	  for (int n=0; n<120; n++){
	    convd[m][n]  = 0.;
	    nconvd[m][n] = 0.;
	  }
	}
     
	for (int l=0; l < TIMEBINS*20; l++){
	  adc_ped_sub_left  = theadccountsl[i][j][k][l] - (theadccountsl[i][j][k][0] + theadccountsl[i][j][k][1])/2.;	  
	  adc_ped_sub       = theadccountsc[i][j][k][l] - (theadccountsc[i][j][k][0] + theadccountsc[i][j][k][1])/2.;
	  adc_ped_sub_right = theadccountsr[i][j][k][l] - (theadccountsr[i][j][k][0] + theadccountsr[i][j][k][1])/2.;

          thebin=thebins[i][j][k][l];
	  
	  if (thebin >= 0 && thebin < 120){
	    convd[0][thebin]  += adc_ped_sub_left;
	    nconvd[0][thebin] += 1.0;
	    
	    convd[1][thebin]  += adc_ped_sub;
	    nconvd[1][thebin] += 1.0;
	    
	    convd[2][thebin]  += adc_ped_sub_right;
	    nconvd[2][thebin] += 1.0;
	    
	  }
	} //loop over timebins
      
	for (int m=0; m<3; m++){
	  for (int n=0; n<120; n++){
	    if(nconvd[m][n]>1.0 && nconvd[m][n] !=0.){
	      convd[m][n] = convd[m][n]/nconvd[m][n];
	    }
	  }
	}
	
	// Call our functions to calculate the cross-talk
	float xl_temp_a = 0.;
	float xl_temp_b = 0.;
	float minl_temp = 0.;
	float xr_temp_a = 0.;
	float xr_temp_b = 0.;
	float minr_temp = 0.;
	mkbins(50.);
	convolution(&xl_temp_a, &xl_temp_b, &minl_temp, &xr_temp_a, &xr_temp_b, &minr_temp);
	
	if (k==0){
	  xtalk_intercept_left[i][j][k]  = 0.0;
	  xtalk_slope_left[i][j][k]      = 0.0;
	  xtalk_chi2_left[i][j][k]       = 0.0;
	  //right side is calculated
	  xtalk_slope_right[i][j][k]     = xl_temp_b;
	  xtalk_intercept_right[i][j][k] = xl_temp_a;
	  xtalk_chi2_right[i][j][k]      = minl_temp;

	}else if(k==79){
	  xtalk_intercept_right[i][j][k]  = 0.0;
	  xtalk_slope_right[i][j][k]      = 0.0;
	  xtalk_chi2_right[i][j][k]       = 0.0;
	  //left side is calculated
	  xtalk_intercept_left[i][j][k]   = xr_temp_a;
	  xtalk_slope_left[i][j][k]       = xr_temp_b;
	  xtalk_chi2_left[i][j][k]        = minr_temp;

	}else{
	  xtalk_intercept_left[i][j][k]  = xl_temp_a;
	  xtalk_intercept_right[i][j][k] = xr_temp_a;
	  xtalk_slope_left[i][j][k]      = xl_temp_b;
	  xtalk_slope_right[i][j][k]     = xr_temp_b;
	  xtalk_chi2_left[i][j][k]       = minl_temp;
	  xtalk_chi2_right[i][j][k]      = minr_temp;
	}
	fff = (j*80)+k;
	double the_xtalk_left_a  = xtalk_intercept_left[i][j][k];
	double the_xtalk_right_a = xtalk_intercept_right[i][j][k];
	double the_xtalk_left_b  = xtalk_slope_left[i][j][k];
	double the_xtalk_right_b = xtalk_slope_right[i][j][k];
	double the_chi2_right    = xtalk_chi2_right[i][j][k];
	double the_chi2_left     = xtalk_chi2_left[i][j][k];
	
	new_xtalk_intercept_right[fff] = the_xtalk_right_a ;
	new_xtalk_intercept_left[fff]  = the_xtalk_left_a ;
	new_xtalk_slope_right[fff]     = the_xtalk_right_b ;
	new_xtalk_slope_left[fff]      = the_xtalk_left_b ;
	new_rchi2[fff]                 = the_chi2_right;
	new_lchi2[fff]                 = the_chi2_left;
	
	std::cout<<"Chamber "<<i<<" Layer "<<j<<" strip "<<k<<" Intercept left "<<new_xtalk_intercept_left[fff]<<"   Slope left "<<new_xtalk_slope_left[fff]<<"   Intercept right "<<new_xtalk_intercept_right[fff]<<"       Slope right "<<new_xtalk_slope_right[fff]<<endl;
      }//loop over strips
    }//loop over layers
    
    //get chamber ID from Igor's mapping    
    int new_crateID = crateID[i];
    int new_dmbID   = dmbID[i];
    std::cout<<" Crate: "<<new_crateID<<" and DMB:  "<<new_dmbID<<std::endl;
    map->crate_chamber(new_crateID,new_dmbID,&chamber_id,&chamber_num,&sector);
    std::cout<<" Above data is for chamber: "<< chamber_id<<" and sector "<<sector<<std::endl;
    
    string test1="CSC_slice";
    string test2="xtalk_slope_left";
    string test3="xtalk_intercept_left";
    string test4="xtalk_chi2_left";
    string test5="xtalk_slope_right";
    string test6="xtalk_intercept_right";
    string test7="xtalk_chi2_right";
    string answer;
    
    std::cout<<" DO you want to send constants to DB? "<<" Please answer y or n for EACH chamber present! "<<std::endl;
    
    std::cin>>answer;
    if(answer=="y"){
      cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, new_xtalk_slope_left,      4, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, new_xtalk_intercept_left,  4, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test4,480, new_lchi2,                 4, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test5,480, new_xtalk_slope_right,     4, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test6,480, new_xtalk_intercept_right, 4, &ret_code);
      cdb->cdb_write(test1,chamber_id,chamber_num,test7,480, new_rchi2,                 4, &ret_code);
      
      std::cout<<" Data SENT to DB! "<<std::endl;
    }else{
      std::cout<<" NO data was sent!!! "<<std::endl;
    }
   }//loop over chambers 
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

  float max, cross0,cross2,min_l,min_r,sum_x=0.0,sumx2=0.;
  float sum_y_left=0.0,sum_y_right=0.0,sum_xy_left=0.0,sum_xy_right=0.0;
  float a_left=0.0,a_right=0.0,b_left=0.0,b_right=0.0,chi2_left=0.0,chi2_right=0.0,chi_left=0.0,chi_right=0.0;
  float aleft=0.0,aright=0.0,bleft=0.0,bright=0.0;
  int i,j,k,l,imax=0;

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
    if(conv[1][i]>max){ 
      max=conv[1][i];
      imax=i;
    }
  }

  //find the max peak time from 3 timebins when line intersects x axis a+b*x=0 -> x=-a/b 
  float time1=-999.0, time2=-999.0, time3=-999.0;
  float data1=-999.0, data2=-999.0, data3=-999.0;
  float peakTime=-999.0;

  time1=imax-1;
  time2=imax;
  time3=imax+1;

  data1=conv[1][imax-1];
  data2=conv[1][imax];
  data3=conv[1][imax+1];

  peakTime=(0.5)*((time1*time1*(data3-data2)+time2*time2*(data1-data3)+time3*time3*(data2-data1))/(time1*(data3-data2)+time2*(data1-data3)+time3*(data2-data1)))*6.25;

  for(l=0;l<3;l++){
    for(i=0;i<119;i++)conv[l][i]=conv[l][i]/max;
  }

  int nobs = 0;
  for (int j=0; j<119; j++){
    if (conv[1][j]>0.6) nobs++;
  }

  for(i=0;i<119;i++){
    cross0=0.0;
    cross2=0.0;
    
    if(conv[1][i]>0.6){
      cross0=conv[0][i]/(conv[0][i]+conv[1][i]+conv[2][i]);
      cross2=conv[2][i]/(conv[0][i]+conv[1][i]+conv[2][i]);
  
      sum_x += i;
      sum_y_left += cross0;
      sum_y_right += cross2;
      sumx2 += i*i;
      sum_xy_left += i*cross0;
      sum_xy_right += i*cross2;
    }
  }  

  //LMS fitting straight line y=a+b*x

  bleft  = ((nobs*sum_xy_left) - (sum_x * sum_y_left))/((nobs*sumx2) - (sum_x*sum_x));
  bright = ((nobs*sum_xy_right) - (sum_x * sum_y_right))/((nobs*sumx2) - (sum_x*sum_x));

  aleft  = ((sum_y_left*sumx2)-(sum_x*sum_xy_left))/((nobs*sumx2)-(sum_x*sum_x));
  aright = ((sum_y_right*sumx2)-(sum_x*sum_xy_right))/((nobs*sumx2)-(sum_x*sum_x));

  for(i=0;i<119;i++ && conv[0][1]>0.6){
    chi2_left  += (cross0 -(aleft+(bleft*i)))*(cross0 -(aleft+(bleft*i)));
    chi2_right += (cross2 -(aright+(bright*i)))*(cross2 -(aright+(bright*i)));
  }	
  
  if(chi_left<min_l){
    min_l=chi_left;
    bleft=bleft;
    aleft=aleft;
  }
  if(chi_right<min_r){
    min_r=chi_right;
    bright=bright;
    aright=aright;
  }
  

  //Now calculating parameters in ns to compensate for drift in peak time  
  b_left  = bleft/6.25;
  b_right = bright/6.25;

  a_left  = aleft  + (bleft*peakTime/6.25);
  a_right = aright + (bright*peakTime/6.25);
  
  *xleft_a   = a_left; 
  *xleft_b   = b_left;
  *min_left  = chi2_left;
  *xright_a  = a_right;
  *xright_b  = b_right;
  *min_right = chi2_right;

} //CONVOLUTION  

