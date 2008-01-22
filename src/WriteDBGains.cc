// -*- C++ -*-
//
// Package:    WriteEcalMiscalibConstants
// Class:      WriteEcalMiscalibConstants
// 
/**\class WriteEcalMiscalibConstants WriteEcalMiscalibConstants.cc CalibCalorimetry/WriteEcalMiscalibConstants/src/WriteEcalMiscalibConstants.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Stephanie BEAUCERON
//         Created:  Tue May 15 16:23:21 CEST 2007
// $Id: WriteDBGains.cc,v 1.1 2007/08/15 16:56:39 boeriu Exp $
//
//


// system include files

// user include files




#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

// DB includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

// user include files
#include "CondFormats/CSCObjects/interface/CSCDBGains.h"
#include "CondFormats/DataRecord/interface/CSCDBGainsRcd.h"
//For Checks

//this one
#include "CalibMuon/CSCCalibration/interface/WriteDBGains.h"

//
// static data member definitions
//

//
// constructors and destructor
//
WriteDBGains::WriteDBGains(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  newTagRequest_ = iConfig.getParameter< std::string > ("NewTagRequest");
}


WriteDBGains::~WriteDBGains()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void WriteDBGains::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  // Intercalib constants
  edm::ESHandle<CSCDBGains> gains;
  iSetup.get<CSCDBGainsRcd>().get(gains);
  const CSCDBGains* mygains = gains.product();
  
  edm::Service<cond::service::PoolDBOutputService> poolDbService;
  if(poolDbService.isAvailable() ){
    if (poolDbService->isNewTagRequest(newTagRequest_) ){
      std::cout<<" Creating a  new one "<<std::endl;
      poolDbService->createNewIOV<const CSCDBGains>(mygains, poolDbService->endOfTime(),newTagRequest_);
      std::cout<<"Done" << std::endl;
    }else{
      std::cout<<"Old One "<<std::endl;
      poolDbService->appendSinceTime<const CSCDBGains>(mygains, poolDbService->currentTime(),newTagRequest_);
    }
  }  
}


// ------------ method called once each job just before starting event loop  ------------
void WriteDBGains::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void WriteDBGains::endJob() {
  std::cout << "Here is the end" << std::endl; 
}
