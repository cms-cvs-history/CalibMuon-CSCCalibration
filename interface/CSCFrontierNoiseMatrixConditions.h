#ifndef _CSCFRONTIERNOISEMATRIXCONDITIONS_H
#define _CSCFRONTIERNOISEMATRIXCONDITIONS_H

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
#include "CondFormats/DataRecord/interface/CSCNoiseMatrixRcd.h"
#include "CalibMuon/CSCCalibration/interface/CSCFrontierNoiseMatrixMap.h"

class CSCFrontierNoiseMatrixConditions: public edm::ESProducer, public edm::EventSetupRecordIntervalFinder  {
   public:
      CSCFrontierNoiseMatrixConditions(const edm::ParameterSet&);
      ~CSCFrontierNoiseMatrixConditions();

      typedef const  CSCNoiseMatrix * ReturnType;

      ReturnType produceNoiseMatrix(const CSCNoiseMatrixRcd&);

   private:
      // ----------member data ---------------------------
    void setIntervalFor(const edm::eventsetup::EventSetupRecordKey &, const edm::IOVSyncValue&, edm::ValidityInterval & );
    
    CSCFrontierNoiseMatrixMap matrix;
};

#endif
