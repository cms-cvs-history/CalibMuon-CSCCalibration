#ifndef _CSCFRONTIERCROSSTALKCONDITIONS_H
#define _CSCFRONTIERCROSSTALKCONDITIONS_H

#include <memory>
#include "FWCore/Framework/interface/SourceFactory.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "CondFormats/CSCObjects/interface/CSCobject.h"
#include "CondFormats/DataRecord/interface/CSCcrosstalkRcd.h"
#include "CalibMuon/CSCCalibration/interface/CSCFrontierCrosstalkMap.h"

class CSCFrontierCrosstalkConditions: public edm::ESProducer, public edm::EventSetupRecordIntervalFinder  {
   public:
      CSCFrontierCrosstalkConditions(const edm::ParameterSet&);
      ~CSCFrontierCrosstalkConditions();

      typedef const  CSCcrosstalk * ReturnType;

      ReturnType produceCrosstalk(const CSCcrosstalkRcd&);

   private:
      // ----------member data ---------------------------
    void setIntervalFor(const edm::eventsetup::EventSetupRecordKey &, const edm::IOVSyncValue&, edm::ValidityInterval & );
    
    CSCFrontierCrosstalkMap crosstalk;
};

#endif