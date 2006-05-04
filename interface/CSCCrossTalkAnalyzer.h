/** 
 * Analyzer for calculating CFEB cross-talk.
 * author S.Durkin, O.Boeriu 25/04/06 
 *   
 */

#include <iostream>
#include "CalibMuon/CSCCalibration/interface/condbc.h"
#include "CalibMuon/CSCCalibration/interface/cscmap.h"
#include "CalibMuon/CSCCalibration/interface/CSCxTalk.h"

class CSCCrossTalkAnalyzer : public edm::EDAnalyzer {
 public:
  explicit CSCCrossTalkAnalyzer(edm::ParameterSet const& conf);
  virtual void analyze(edm::Event const& e, edm::EventSetup const& iSetup);
  
#define CHAMBERS 9
#define LAYERS 6
#define STRIPS 80
#define TIMEBINS 8

  ~CSCCrossTalkAnalyzer(){
    condbc *cdb = new condbc();
    cscmap *map = new cscmap();
    
    
    ////////////////////////////////////////////////////////////////////iuse==strip-1
    // Now that we have filled our array, extract convd and nconvd
    float adc_ped_sub_left = -999.;
    float adc_ped_sub = -999.;
    float adc_ped_sub_right = -999.;
    int thebin;
    float sum=0.0;
    float mean=0;
    
    for (int i=0; i<NChambers; i++){
      for (int j=0; j<LAYERS; j++){
	mean=0.,sum=0.;
	for (int s=0; s<STRIPS; s++) {
	  //re-zero convd and nconvd
	  for (int m=0; m<3; m++){
	    for (int n=0; n<120; n++){
	      binsConv.convd[m][n]  = 0.;
	      binsConv.nconvd[m][n] = 0.;
	    }
	  }
	  
	  for (int l=0; l < TIMEBINS*20; l++){
	    adc_ped_sub_left  = theadccountsl[i][j][s][l] - (theadccountsl[i][j][s][0] + theadccountsl[i][j][s][1])/2.;
	    adc_ped_sub       = theadccountsc[i][j][s][l] - (theadccountsc[i][j][s][0] + theadccountsc[i][j][s][1])/2.;
	    adc_ped_sub_right = theadccountsr[i][j][s][l] - (theadccountsr[i][j][s][0] + theadccountsr[i][j][s][1])/2.;
	    
	    thebin=thebins[i][j][s][l];
	    
	    if (thebin >= 0 && thebin < 120){
	      binsConv.convd[0][thebin]  += adc_ped_sub_left;
	      binsConv.nconvd[0][thebin] += 1.0;
	      
	      binsConv.convd[1][thebin]  += adc_ped_sub;
	      binsConv.nconvd[1][thebin] += 1.0;
	      
	      binsConv.convd[2][thebin]  += adc_ped_sub_right;
	      binsConv.nconvd[2][thebin] += 1.0;
	      
	    }
	  } //loop over timebins
	  
	  for (int m=0; m<3; m++){
	    for (int n=0; n<120; n++){
	      if(binsConv.nconvd[m][n]>1.0 && binsConv.nconvd[m][n] !=0.){
		binsConv.convd[m][n] = binsConv.convd[m][n]/binsConv.nconvd[m][n];
	      }
	    }
	  }
	  
	  // Call our functions first time to calculate mean pf peak time over a layer
	  float xl_temp_a = 0.0;
	  float xl_temp_b = 0.0;
	  float minl_temp = 0.0;
	  float xr_temp_a = 0.0;
	  float xr_temp_b = 0.0;
	  float minr_temp = 0.0;
	  float pTime     = 0.0;
	  
	  binsConv.mkbins(50.);
	  binsConv.convolution(&xl_temp_a, &xl_temp_b, &minl_temp, &xr_temp_a, &xr_temp_b, &minr_temp, &pTime);
	  myPeakTime[i][j][s] = pTime;
	  sum=sum+myPeakTime[i][j][s];
	  mean = sum/80.;
	}
	
	for (int k=0; k<STRIPS; k++){
	  // re-zero convd and nconvd 
	  for (int m=0; m<3; m++){
	    for (int n=0; n<120; n++){
	      binsConv.convd[m][n]  = 0.;
	      binsConv.nconvd[m][n] = 0.;
	    }
	  }
	  
	  for (int l=0; l < TIMEBINS*20; l++){
	    adc_ped_sub_left  = theadccountsl[i][j][k][l] - (theadccountsl[i][j][k][0] + theadccountsl[i][j][k][1])/2.;	  
	    adc_ped_sub       = theadccountsc[i][j][k][l] - (theadccountsc[i][j][k][0] + theadccountsc[i][j][k][1])/2.;
	    adc_ped_sub_right = theadccountsr[i][j][k][l] - (theadccountsr[i][j][k][0] + theadccountsr[i][j][k][1])/2.;
	    
	    thebin=thebins[i][j][k][l];
	    
	    if (thebin >= 0 && thebin < 120){
	      binsConv.convd[0][thebin]  += adc_ped_sub_left;
	      binsConv.nconvd[0][thebin] += 1.0;
	      
	      binsConv.convd[1][thebin]  += adc_ped_sub;
	      binsConv.nconvd[1][thebin] += 1.0;
	      
	      binsConv.convd[2][thebin]  += adc_ped_sub_right;
	      binsConv.nconvd[2][thebin] += 1.0;
	      
	    }
	  } //loop over timebins
	  
	  for (int m=0; m<3; m++){
	    for (int n=0; n<120; n++){
	      if(binsConv.nconvd[m][n]>1.0 && binsConv.nconvd[m][n] !=0.){
		binsConv.convd[m][n] = binsConv.convd[m][n]/binsConv.nconvd[m][n];
	      }
	    }
	  }
	  //////////////////////////////////////////////////////////////////
	  // Call our functions a second time to calculate the cross-talk //
	  //////////////////////////////////////////////////////////////////
	  float xl_temp_a = 0.;
	  float xl_temp_b = 0.;
	  float minl_temp = 0.;
	  float xr_temp_a = 0.;
	  float xr_temp_b = 0.;
	  float minr_temp = 0.;
	  float pTime     = 0.;
	  
	  binsConv.mkbins(50.);
	  binsConv.convolution(&xl_temp_a, &xl_temp_b, &minl_temp, &xr_temp_a, &xr_temp_b, &minr_temp, &pTime);
	  
	  if (k==0){
	    xtalk_intercept_left[i][j][k]  = 0.0;
	    xtalk_slope_left[i][j][k]      = 0.0;
	    xtalk_chi2_left[i][j][k]       = 0.0;
	    //right side is calculated
	    xtalk_slope_right[i][j][k]     = xl_temp_b;
	    xtalk_intercept_right[i][j][k] = xl_temp_a;
	    xtalk_chi2_right[i][j][k]      = minl_temp;
	    myPeakTime[i][j][k]            = pTime;
	  }else if(k==79){
	    xtalk_intercept_right[i][j][k]  = 0.0;
	    xtalk_slope_right[i][j][k]      = 0.0;
	    xtalk_chi2_right[i][j][k]       = 0.0;
	    //left side is calculated
	    xtalk_intercept_left[i][j][k]   = xr_temp_a;
	    xtalk_slope_left[i][j][k]       = xr_temp_b;
	    xtalk_chi2_left[i][j][k]        = minr_temp;
	    myPeakTime[i][j][k]             = pTime;
	  }else{
	    xtalk_intercept_left[i][j][k]  = xl_temp_a;
	    xtalk_intercept_right[i][j][k] = xr_temp_a;
	    xtalk_slope_left[i][j][k]      = xl_temp_b;
	    xtalk_slope_right[i][j][k]     = xr_temp_b;
	    xtalk_chi2_left[i][j][k]       = minl_temp;
	    xtalk_chi2_right[i][j][k]      = minr_temp;
	    myPeakTime[i][j][k]            = pTime;
	  }
	  
	  fff = (j*80)+k;
	  float the_xtalk_left_a  = xtalk_intercept_left[i][j][k];
	  float the_xtalk_right_a = xtalk_intercept_right[i][j][k];
	  float the_xtalk_left_b  = xtalk_slope_left[i][j][k];
	  float the_xtalk_right_b = xtalk_slope_right[i][j][k];
	  float the_chi2_right    = xtalk_chi2_right[i][j][k];
	  float the_chi2_left     = xtalk_chi2_left[i][j][k];
	  float the_peakTime      = myPeakTime[i][j][k]; 
	  
	  new_xtalk_intercept_right[fff] = the_xtalk_right_a ;
	  new_xtalk_intercept_left[fff]  = the_xtalk_left_a ;
	  new_xtalk_slope_right[fff]     = the_xtalk_right_b ;
	  new_xtalk_slope_left[fff]      = the_xtalk_left_b ;
	  new_rchi2[fff]                 = the_chi2_right;
	  new_lchi2[fff]                 = the_chi2_left;
	  newPeakTime[fff]               = the_peakTime;
	  newMeanPeakTime[fff]           = the_peakTime-mean;
	  
	  std::cout<<"Cham "<<i<<" Layer "<<j<<" strip "<<k<<" IntL "<<new_xtalk_intercept_left[fff]<<"   SlopeL "<<new_xtalk_slope_left[fff]<<"   IntR "<<new_xtalk_intercept_right[fff]<<"   SlopeR "<<new_xtalk_slope_right[fff]<<"   diff "<<newMeanPeakTime[fff]<<endl;
	}//loop over strips
      }//loop over layers

      int new_crateID = crateID[i];
      int new_dmbID   = dmbID[i];
      std::cout<<" Crate: "<<new_crateID<<" and DMB:  "<<new_dmbID<<std::endl;
      map->crate_chamber(new_crateID,new_dmbID,&chamber_id,&chamber_num,&sector);
      std::cout<<" Above data is for chamber:: "<< chamber_id<<" and sector:  "<<sector<<std::endl;
      
      std::string test1 = "CSC_slice";
      std::string test2 = "xtalk_slope_left";
      std::string test3 = "xtalk_intercept_left";
      std::string test4 = "xtalk_chi2_left";
      std::string test5 = "xtalk_slope_right";
      std::string test6 = "xtalk_intercept_right";
      std::string test7 = "xtalk_chi2_right";
      std::string test8 = "time_spread";
      std::string answer;
      std::string bad_number = "nan";
      
      std::cout<<" DO you want to send constants to DB? "<<" Please answer y or n for EACH chamber present! "<<std::endl;
      std::cin>>answer;
      if(answer=="y"){
	if(new_xtalk_slope_left[fff] != new_xtalk_slope_left[fff]) {
	  std::cout<<"this is my xtalk left "<< new_xtalk_slope_left[fff]<<std::endl;
	}
	cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, new_xtalk_slope_left,      4, &ret_code);
	cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, new_xtalk_intercept_left,  4, &ret_code);
	cdb->cdb_write(test1,chamber_id,chamber_num,test4,480, new_lchi2,                 4, &ret_code);
	cdb->cdb_write(test1,chamber_id,chamber_num,test5,480, new_xtalk_slope_right,     4, &ret_code);
	cdb->cdb_write(test1,chamber_id,chamber_num,test6,480, new_xtalk_intercept_right, 4, &ret_code);
	cdb->cdb_write(test1,chamber_id,chamber_num,test7,480, new_rchi2,                 4, &ret_code);
	cdb->cdb_write(test1,chamber_id,chamber_num,test8,480, newMeanPeakTime,           4, &ret_code);
	
	std::cout<<" Data SENT to DB! "<<std::endl;
      }else{
	std::cout<<" NO data was sent!!! "<<std::endl;
      }
    }//loop over chambers
  }
  
 private:
  // variables persistent across events should be declared here.
  int eventNumber,evt,strip,misMatch,fff,ret_code,length;
  int i_chamber,i_layer,reportedChambers,chamber_num,sector, NChambers ;
  int dmbID[CHAMBERS],crateID[CHAMBERS] ;
  std::vector<int> adc;
  std::string chamber_id;
  int thebins[CHAMBERS][LAYERS][STRIPS][TIMEBINS*20];
  int theadccountsc[CHAMBERS][LAYERS][STRIPS][TIMEBINS*20];
  int theadccountsl[CHAMBERS][LAYERS][STRIPS][TIMEBINS*20];
  int theadccountsr[CHAMBERS][LAYERS][STRIPS][TIMEBINS*20];
  float pedMean,time;
  float thetime[CHAMBERS][LAYERS][STRIPS][TIMEBINS*20];
  float xtalk_intercept_left[CHAMBERS][LAYERS][STRIPS];
  float xtalk_intercept_right[CHAMBERS][LAYERS][STRIPS];
  float xtalk_slope_left[CHAMBERS][LAYERS][STRIPS];
  float xtalk_slope_right[CHAMBERS][LAYERS][STRIPS];
  float xtalk_chi2_left[CHAMBERS][LAYERS][STRIPS];
  float xtalk_chi2_right[CHAMBERS][LAYERS][STRIPS];
  float myPeakTime[CHAMBERS][LAYERS][STRIPS];
  float myMeanPeakTime[CHAMBERS][LAYERS][STRIPS];
  float new_xtalk_intercept_right[480];
  float new_xtalk_intercept_left[480];
  float new_xtalk_slope_right[480];
  float new_xtalk_slope_left[480];
  float array_meanPeakTime[CHAMBERS][LAYERS][STRIPS];
  float new_rchi2[480];
  float new_lchi2[480];
  float newPeakTime[480];
  float newMeanPeakTime[480];
};