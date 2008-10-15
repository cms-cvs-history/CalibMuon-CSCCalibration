#ifndef _CSCL1TPPARAMETERSCONDITIONS_H
#define _CSCL1TPPARAMETERSCONDITIONS_H

#include <memory>
#include <cmath>
#include "FWCore/Framework/interface/SourceFactory.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include "CondFormats/CSCObjects/interface/CSCL1TPParameters.h"
#include "CondFormats/DataRecord/interface/CSCL1TPParametersRcd.h"

class CSCL1TPParametersConditions: public edm::ESProducer, public edm::EventSetupRecordIntervalFinder  {
 public:
  CSCL1TPParametersConditions(const edm::ParameterSet&);
  ~CSCL1TPParametersConditions();
  

  inline static CSCL1TPParameters *  prefillCSCL1TPParameters();

  typedef const  CSCL1TPParameters * ReturnType;
  
  ReturnType produceCSCL1TPParameters(const CSCL1TPParametersRcd&);
  
 private:
  // ----------member data ---------------------------
  void setIntervalFor(const edm::eventsetup::EventSetupRecordKey &, const edm::IOVSyncValue&, edm::ValidityInterval & );
  CSCL1TPParameters *CSCl1TPParameters ;

};

#include<fstream>
#include<vector>
#include<iostream>

// to workaround plugin library
inline CSCL1TPParameters *  CSCL1TPParametersConditions::prefillCSCL1TPParameters()
{

  CSCL1TPParameters * cnl1tp = new CSCL1TPParameters();
    
  ////to be filled

 return cnl1tp;
}


#endif