#include <memory>
#include "boost/shared_ptr.hpp"

// user include files
#include "FWCore/Framework/interface/SourceFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"

//FW include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

//CSCObjects
#include "CondFormats/CSCObjects/interface/CSCobject.h"
#include "CondFormats/CSCObjects/interface/CSCNoiseMatrix.h"
#include "CalibMuon/CSCCalibration/interface/CSCFakeNoiseMatrixConditions.h"
#include "CondFormats/DataRecord/interface/CSCNoiseMatrixRcd.h"

void CSCFakeNoiseMatrixMap::prefillNoiseMatrixMap(){
  
  const CSCDetId& detId = CSCDetId();
  cnmatrix = new CSCNoiseMatrix();
  
  int max_istrip,id_layer,max_ring,max_cham;
  //endcap=1 to 2,station=1 to 4, ring=1 to 4,chamber=1 to 36,layer=1 to 6 
  
  for(int iendcap=detId.minEndcapId(); iendcap<=detId.maxEndcapId(); iendcap++){
    for(int istation=detId.minStationId() ; istation<=detId.maxStationId(); istation++){
      max_ring=detId.maxRingId();
      //station 4 ring 4 not there(36 chambers*2 missing)
      //3 rings max this way of counting (ME1a & b)
      if(istation==1)    max_ring=3;
      if(istation==2)    max_ring=2;
      if(istation==3)    max_ring=2;
      if(istation==4)    max_ring=1;
	
      for(int iring=detId.minRingId(); iring<=max_ring; iring++){
	max_istrip=80;
	max_cham=detId.maxChamberId(); 
	if(istation==1 && iring==1)    max_cham=36;
	if(istation==1 && iring==2)    max_cham=36;
	if(istation==1 && iring==3)    max_cham=36;
	if(istation==2 && iring==1)    max_cham=18;
	if(istation==2 && iring==2)    max_cham=36;
	if(istation==3 && iring==1)    max_cham=18;
	if(istation==3 && iring==2)    max_cham=36;
	if(istation==4 && iring==1)    max_cham=18;
	
	for(int ichamber=detId.minChamberId(); ichamber<=max_cham; ichamber++){
	  for(int ilayer=detId.minLayerId(); ilayer<=detId.maxLayerId(); ilayer++){
	    //station 1 ring 3 has 64 strips per layer instead of 80 
	    if(istation==1 && iring==3)   max_istrip=64;

	    std::vector<CSCNoiseMatrix::Item> itemvector;
	    itemvector.resize(max_istrip);
	    id_layer = 100000*iendcap + 10000*istation + 1000*iring + 10*ichamber + ilayer;
	    
	    for(int istrip=0;istrip<max_istrip;istrip++){
	      
	      itemvector[istrip].elem33 = 10.00;
	      itemvector[istrip].elem34 = 4.00;
	      itemvector[istrip].elem35 = 10.0;
	      itemvector[istrip].elem44 = 2.50;
	      itemvector[istrip].elem45 = 3.00;
	      itemvector[istrip].elem46 = 10.0;
	      itemvector[istrip].elem55 = 2.00;
	      itemvector[istrip].elem56 = 5.00;
	      itemvector[istrip].elem57 = 10.0;
	      itemvector[istrip].elem66 = 3.00;
	      itemvector[istrip].elem67 = 4.00;
	      itemvector[istrip].elem77 = 10.0;
	      cnmatrix->matrix[id_layer]=itemvector;
	    }
	  }
	}
      }
    }
  }
}

CSCFakeNoiseMatrixConditions::CSCFakeNoiseMatrixConditions(const edm::ParameterSet& iConfig)
{
  //the following line is needed to tell the framework what
  // data is being produced
  matrix.prefillNoiseMatrixMap();
  // added by Zhen (changed since 1_2_0)
  setWhatProduced(this,&CSCFakeNoiseMatrixConditions::produceNoiseMatrix);
  findingRecord<CSCNoiseMatrixRcd>();
  //now do what ever other initialization is needed
}


CSCFakeNoiseMatrixConditions::~CSCFakeNoiseMatrixConditions()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
CSCFakeNoiseMatrixConditions::ReturnType
CSCFakeNoiseMatrixConditions::produceNoiseMatrix(const CSCNoiseMatrixRcd& iRecord)
{
    matrix.prefillNoiseMatrixMap();
    // Added by Zhen, need a new object so to not be deleted at exit
    CSCNoiseMatrix* mydata=new CSCNoiseMatrix(matrix.get());
    
    return mydata;

}

 void CSCFakeNoiseMatrixConditions::setIntervalFor(const edm::eventsetup::EventSetupRecordKey &, const edm::IOVSyncValue&,
 edm::ValidityInterval & oValidity)
 {
 oValidity = edm::ValidityInterval(edm::IOVSyncValue::beginOfTime(),edm::IOVSyncValue::endOfTime());
 
 }