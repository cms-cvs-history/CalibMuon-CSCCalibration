//This will calculate left and right cross-talk, create array 480 and send it to DB.
//******BUT old ORCA code******

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "EMuCalib/interface/MuEndCalibration.h"
#include "TriDAS/emu/emuDAQ/DDUReadout/include/FileReaderDDU.h"
#include "Utilities/Configuration/interface/Architecture.h"
#include "Muon/MERawFormat/interface/MuEndDDUEventData.h"
#include "Muon/MERawFormat/interface/MuEndEventData.h"
#include "Muon/MERawFormat/interface/MuEndAnodeData.h"
#include "Muon/MERawFormat/interface/MuEndALCTHeader.h"
#include "Muon/MERawFormat/interface/MuEndCLCTData.h"
#include "Muon/MERawFormat/interface/MuEndTMBData.h"
#include "Muon/MERawFormat/interface/MuEndDDUTrailer.h"
#include "Muon/MERawFormat/interface/MuEndDMBHeader.h"
#include "Muon/MERawFormat/interface/MuEndCFEBData.h"

#include "condbc.h"
#include "cscmap.h" 
#include "string"

//functions
float elec(float t,float vs);
void mkbins(float vs);
void convolution(float *xleft_a, float *xleft_b, float *min_left, float* xright_a, float* xright_b, float *min_right);
void fit();

float conv[3][120];
float conve[120];
float convd[3][120];
float nconvd[3][120];
double adc_val[6][80];
double left_xtalk_a[7][81]; 
double left_xtalk_b[7][81];
double min_left[7][81];
double right_xtalk_a[7][81];
double right_xtalk_b[7][81];
double min_right[7][81];


float norm,t;
int maxadc=0, maxchan, maxtimeBin, adcmax;

float elec(float time,float vs);
float chifit_ts(float tql,float tq,float tqh);


int main(int argc, char **argv) {
 
  //initialize arrays
  for (int aa = 0; aa < 7; aa++){
    for (int bb = 0; bb < 81; bb++){
      adc_val[aa][bb] = -999.;
      left_xtalk_a[aa][bb] = -999.;
      left_xtalk_b[aa][bb] = -999.;
      min_left[aa][bb] = -999.;
      right_xtalk_a[aa][bb] = -999.;
      right_xtalk_b[aa][bb] = -999.;
      min_right[aa][bb] = -999.;
    }
  }
  
  int maxEvents = 500000;
  string datafile=argv[1];
  int dmbID,crateID;
  int reportedChambers =0,l=0;
  float ped0=0.0, ped1=0.0;
  float adc_ped_sub=0.0,adc=0.0;
  int count=-1, middle_strip=-999, offset=0;
  int ret_code;
  int chamber_num,sector;
  string chamber_id;
  int fff;
  condbc *cdb = new condbc(); 
  cscmap *map = new cscmap();

  double new_rxtalk[480];
  double new_lxtalk[480];
  double new_rintercept[480];
  double new_lintercept[480];

  
  for (int i=0;i<480;i++){
    new_lxtalk[i]=0.;
    new_rxtalk[i]=0.;
    new_rintercept[i]=0.;
    new_lintercept[i]=0.;
  }

  if (argv[2]) maxEvents = (int) atof(argv[2]);
  FileReaderDDU ddu;
  ddu.openFile(datafile);
  int const maxCham = 3;

  for (int evt = 0; ddu.readNextEvent() && (evt < maxEvents); ++evt) {
    std::cout << "---------- Event: " << evt << "--------------" << endl;
    count++;
    
    MuEndDDUEventData dduEvent((short unsigned int *)ddu.data());
    const vector<MuEndEventData> & cscData = dduEvent.cscData();
    reportedChambers += dduEvent.header().ncsc();
    int NChambers = cscData.size();
    int repChambers = dduEvent.header().ncsc();
    std::cout << " Reported Chambers = " << repChambers <<"   "<<NChambers<< endl;

    for (int i_chamber=0; i_chamber<maxCham; i_chamber++) {
      
      MuEndDMBHeader &thisDMBheader = cscData[i_chamber].dmbHeader();
      dmbID = cscData[i_chamber].dmbHeader().dmbID();  //DMB ID 
      crateID = cscData[i_chamber].dmbHeader().crateID(); //Crate ID
      
      for(int ilayer = 1; ilayer <= 6; ++ilayer) {
	
	for (int iii=0; iii<3; iii++){
	  for (int jjj=0; jjj<120; jjj++){
	    convd[iii][jjj] = 0.;
	    nconvd[iii][jjj] = 0.;
	  }
	}
	
	for(int icfeb =0; icfeb<5 ;icfeb++){
	  //To get the strip number (events 0->31)  do: 1 + icfeb*16
	  //To get the strip number (events 32->64) do: 2 + icfeb*16
	  //To get the strip number (events 65->96) do: 3 + icfeb*16
	  //...and so on...
	  
	  int true_strip = -999;
	  
	  if (thisDMBheader.cfebAvailable(icfeb)){
	    const MuEndCFEBData * data = cscData[i_chamber].cfebData(icfeb);
	    if (evt < 32)                   {middle_strip = 1; true_strip = 1 + icfeb* 16; offset=0;}
	    if ((evt >=32) && (evt < 64))   {middle_strip = 2; true_strip = 2 + icfeb* 16; offset=1;}
	    if ((evt >=64) && (evt < 96))   {middle_strip = 3; true_strip = 3 + icfeb* 16; offset=2;}
	    if ((evt >=96) && (evt < 128))  {middle_strip = 4; true_strip = 4 + icfeb* 16; offset=3;}
	    if ((evt >=128) && (evt < 160)) {middle_strip = 5; true_strip = 5 + icfeb* 16; offset=4;}
	    if ((evt >=160) && (evt < 192)) {middle_strip = 6; true_strip = 6 + icfeb* 16; offset=5;}
	    if ((evt >=192) && (evt < 224)) {middle_strip = 7; true_strip = 7 + icfeb* 16; offset=6;}
	    if ((evt >=224) && (evt < 256)) {middle_strip = 8; true_strip = 8 + icfeb* 16; offset=7;}
	    if ((evt >=256) && (evt < 288)) {middle_strip = 9; true_strip = 9 + icfeb* 16; offset=8;}
	    if ((evt >=288) && (evt < 320)) {middle_strip = 10; true_strip = 10 + icfeb* 16; offset=9;}
	    if ((evt >=320) && (evt < 352)) {middle_strip = 11; true_strip = 11 + icfeb* 16; offset=10;}
	    if ((evt >=352) && (evt < 384)) {middle_strip = 12; true_strip = 12 + icfeb* 16; offset=11;}
	    if ((evt >=384) && (evt < 416)) {middle_strip = 13; true_strip = 13 + icfeb* 16; offset=12;}
	    if ((evt >=416) && (evt < 448)) {middle_strip = 14; true_strip = 14 + icfeb* 16; offset=13;}
	    if ((evt >=448) && (evt < 480)) {middle_strip = 15; true_strip = 15 + icfeb* 16; offset=14;}
	    if ((evt >=480) && (evt < 512)) {middle_strip = 16; true_strip = 16 + icfeb* 16; offset=15;}

// 	    if (evt < 20)                   {middle_strip = 1; true_strip = 1 + icfeb* 16; offset=0;}
// 	    if ((evt >=20) && (evt < 40))   {middle_strip = 2; true_strip = 2 + icfeb* 16; offset=1;}
// 	    if ((evt >=40) && (evt < 60))   {middle_strip = 3; true_strip = 3 + icfeb* 16; offset=2;}
// 	    if ((evt >=60) && (evt < 80))  {middle_strip = 4; true_strip = 4 + icfeb* 16; offset=3;}
// 	    if ((evt >=80) && (evt < 100)) {middle_strip = 5; true_strip = 5 + icfeb* 16; offset=4;}
// 	    if ((evt >=100) && (evt < 120)) {middle_strip = 6; true_strip = 6 + icfeb* 16; offset=5;}
// 	    if ((evt >=120) && (evt < 140)) {middle_strip = 7; true_strip = 7 + icfeb* 16; offset=6;}
// 	    if ((evt >=140) && (evt < 160)) {middle_strip = 8; true_strip = 8 + icfeb* 16; offset=7;}
// 	    if ((evt >=160) && (evt < 180)) {middle_strip = 9; true_strip = 9 + icfeb* 16; offset=8;}
// 	    if ((evt >=180) && (evt < 200)) {middle_strip = 10; true_strip = 10 + icfeb* 16; offset=9;}
// 	    if ((evt >=200) && (evt < 220)) {middle_strip = 11; true_strip = 11 + icfeb* 16; offset=10;}
// 	    if ((evt >=220) && (evt < 240)) {middle_strip = 12; true_strip = 12 + icfeb* 16; offset=11;}
// 	    if ((evt >=240) && (evt < 260)) {middle_strip = 13; true_strip = 13 + icfeb* 16; offset=12;}
// 	    if ((evt >=260) && (evt < 280)) {middle_strip = 14; true_strip = 14 + icfeb* 16; offset=13;}
// 	    if ((evt >=280) && (evt < 300)) {middle_strip = 15; true_strip = 15 + icfeb* 16; offset=14;}
// 	    if ((evt >=300) && (evt < 320)) {middle_strip = 16; true_strip = 16 + icfeb* 16; offset=15;}
	    
	    for(int cfebStrip = middle_strip - 1; cfebStrip <= middle_strip + 1; ++cfebStrip) {

	      for(unsigned int timeBin = 0; timeBin < data->nTimeSamples(); ++timeBin) {

		if (cfebStrip == 0 &&  true_strip == 1) {
		  adc = 0.;
		  ped0 = 0.;
		  ped1 = 0.;
		  
		} else if (cfebStrip == 17 && true_strip == 80) {
		  continue;
		  
		} else if (cfebStrip == 0 && true_strip != 1){
		  int prev_cfeb = -9;
		  if ((true_strip - 1) <= 16) {
		    prev_cfeb = 0;
		  } else if (((true_strip -1) > 16) && ((true_strip -1) <= 32)){
		    prev_cfeb = 1;
		  } else if (((true_strip -1) > 32) && ((true_strip -1) <= 48)){
		    prev_cfeb = 2;
		  } else if (((true_strip -1) > 48) && ((true_strip -1) <= 64)){
		    prev_cfeb = 3;
		  } else if (((true_strip -1) > 64) && ((true_strip -1) <= 80)){
		    prev_cfeb = 4;

		  }

		  const MuEndCFEBData *tempdata = cscData[i_chamber].cfebData(prev_cfeb);
		  adc = tempdata->adcCounts(ilayer, 16, timeBin);
		  ped0 = tempdata->adcCounts(ilayer, 16, 0);
		  ped1 = tempdata->adcCounts(ilayer, 16, 1);
		  
		  
		} else if (cfebStrip == 17 && true_strip != 80){
		  int next_cfeb = -9;
		  if ((true_strip + 1) <= 16) {
		    next_cfeb = 0;
		  } else if (((true_strip + 1) > 16) && ((true_strip + 1) <= 32)){
		    next_cfeb = 1;
		  } else if (((true_strip + 1) > 32) && ((true_strip + 1) <= 48)){
		    next_cfeb = 2;
		  } else if (((true_strip + 1) > 48) && ((true_strip + 1) <= 64)){
		    next_cfeb = 3;
		  } else if (((true_strip + 1) > 64) && ((true_strip + 1) <= 80)){
		    next_cfeb = 4;
		  }
		  // Now we have to load the data one cfeb behind the current cfeb
		  const MuEndCFEBData *tempdata = cscData[i_chamber].cfebData(next_cfeb);
		  adc = tempdata->adcCounts(ilayer, 1, timeBin);
		  ped0 = tempdata->adcCounts(ilayer, 1, 0);
		  ped1 = tempdata->adcCounts(ilayer, 1, 1);
		  
		} else {
		  adc = data->adcCounts(ilayer, cfebStrip, timeBin);
		  ped0 = data->adcCounts(ilayer, cfebStrip, 0);
		  ped1 = data->adcCounts(ilayer, cfebStrip, 1);
		}
		adc_ped_sub = adc -((ped0 + ped1)/2) ;
		adc_val[ilayer][true_strip]= adc_ped_sub; 
		t = (50. * timeBin)-(evt * 6.25)+194.0+(200*offset);
		l=(t-200.+ 0.001)/6.24;
		if(l>=0 && l<120 && true_strip !=80 ){
		  convd[cfebStrip - middle_strip + 1][l]=convd[cfebStrip - middle_strip + 1][l]+adc_val[ilayer][true_strip]; 
		  nconvd[cfebStrip - middle_strip + 1][l]=nconvd[cfebStrip - middle_strip + 1][l]+1.0;
		} else if (l>=0 && l<120 && true_strip ==80){
		  convd[cfebStrip - middle_strip + 1][l]=0;
		  nconvd[cfebStrip - middle_strip + 1][l]=0;
		}
	      }// timeBin loop
	    } //cfebStrip loop
	  } //if cfebAvailable 

	  for(int iii=0;iii<119;iii++){
	    for( int kkk=0;kkk<3;kkk++){
	      convd[kkk][iii]=convd[kkk][iii]/nconvd[kkk][iii];
	      if (nconvd[kkk][iii] == 0) convd[kkk][iii] = 0;
	    }
	  }
	  
	  if (count ==32) count =0;	  
	  if (count == 31){

	    mkbins(50.);
	    float xl_temp_a = 0.;
	    float xl_temp_b = 0.;
	    float minl_temp = 0.;
	    float xr_temp_a = 0.;
	    float xr_temp_b = 0.;
	    float minr_temp = 0.;
	    convolution(&xl_temp_a, &xl_temp_b, &minl_temp, &xr_temp_a, &xr_temp_b, &minr_temp);
	    left_xtalk_a[ilayer][true_strip] = xl_temp_a;
	    left_xtalk_b[ilayer][true_strip] = xl_temp_b;
	    min_left[ilayer][true_strip] = minl_temp;
	    right_xtalk_a[ilayer][true_strip] = xr_temp_a;
	    right_xtalk_b[ilayer][true_strip] = xr_temp_b;
	    min_right[ilayer][true_strip] = minr_temp;
	    if(true_strip == 80){
	      right_xtalk_a[ilayer][true_strip] = 0;
	      right_xtalk_b[ilayer][true_strip] = 0;
	      min_right[ilayer][true_strip] = 0;
	    }
	  }
	}// icfeb loop
      }//ilayer
    } // i_chamber loop
  } //event loop
  
  for (int ij=0; ij<6; ij++){
    for (int jk=0; jk<80; jk++){
      fff = (ij*80)+jk;
      double the_xtalk_left_a = left_xtalk_a[ij+1][jk+1];
      double the_xtalk_right_a = right_xtalk_a[ij+1][jk+1];
      double the_xtalk_left_b = left_xtalk_b[ij+1][jk+1];
      double the_xtalk_right_b = right_xtalk_b[ij+1][jk+1];
      
      new_rxtalk[fff] =the_xtalk_right_a ;
      new_lxtalk[fff] =the_xtalk_left_a ;
      new_rintercept[fff] =the_xtalk_right_b ;
      new_lintercept[fff] =the_xtalk_left_b ;
      
    }
  }
  
  //info for database
  string test1="CSC_slice";
  string test2="xtalk_slope_left";
  string test3="xtalk_intercept_right";
  string test4="xtalk_slope_right";
  string test5="xtalk_intercept_right";
  
  //from mapping DB
  std::cout<<"Here is crate and DMB: "<<crateID<<"  "<<dmbID<<endl;
  map->crate_chamber(crateID,dmbID,&chamber_id,&chamber_num,&sector);
  std::cout<<"This is from mapping: "<< chamber_id<<"  "<<chamber_num<<"  "<<sector<<endl;
  
  //*******to send this array to DB uncomment the next four lines*************
  // cdb->cdb_write(test1,chamber_id,chamber_num,test2,480, new_lxtalk,2, &ret_code);
  // cdb->cdb_write(test1,chamber_id,chamber_num,test3,480, new_lintercept,2, &ret_code);
  // cdb->cdb_write(test1,chamber_id,chamber_num,test4,480, new_rxtalk,2, &ret_code);
  // cdb->cdb_write(test1,chamber_id,chamber_num,test5,480, new_rintercept,2, &ret_code);
  
}// main



//square wave fcn convoluted with 1/(t+2.5)
float elec(float t,float vs)
{
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

//do convolution
void convolution(float *xleft_a, float *xleft_b, float *min_left, float *xright_a, float *xright_b, float *min_right){
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
	//c_left=c_left;
      }
      if(chi_right<min_r){
	min_r=chi_right;
	b_right=b_right;
	a_right=a_right;
	//c_right=c_right;
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
 
  std::cout << "Left: "<<min_l <<"  "<<a_left <<"  "<< b_left<<endl;
  std::cout << "Right: "<<min_r <<" "<<a_right <<" "<<b_right<< endl;
}

//fit function
void fit()
{
  float xchi,t,xchimin,f,tc,t0,mtc,mt0;
  int i,j,k;
  xchimin=10000000.0;
  for(j=0;j<2000;j++){
    tc=10.+j*0.05;
    for(k=0;k<2000;k++){
      t0=k*0.1;  
      xchi=0.0;
      for(i=12;i<60;i++){
	t=i*6.25;
	f=0.0;
	if(t-t0>0){
	  f=pow(t-t0,4)*exp(-(t-t0)/tc)/256./pow(tc,4)/exp(-4.0);
	}
	  xchi=xchi+(conv[1][i]-f)*(conv[1][i]-f);
      }

      if(xchi<xchimin){

	mt0=t0;
	mtc=tc;
	xchimin=xchi;
      }
    } 
  }
  
  for(i=0;i<750;i++){
    t=i+mt0;
    f=pow(t-mt0,4)*exp(-(t-mt0)/mtc)/256./pow(mtc,4)/exp(-4.0);
  }
} // fit

