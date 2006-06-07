#include <iostream>
#include <vector>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CalibMuon/CSCCalibration/interface/CSCAFEBThrAnalysis.h"


class CSCAFEBAnalyzer : public edm::EDAnalyzer {
public:
  explicit CSCAFEBAnalyzer(edm::ParameterSet const& conf);
  virtual void analyze(edm::Event const& e, edm::EventSetup const& iSetup);
  virtual void endJob();
private:
  /// variables persistent across events should be declared here.

  CSCAFEBThrAnalysis analysisthr_;
};

CSCAFEBAnalyzer::CSCAFEBAnalyzer(edm::ParameterSet const& conf) {

  /// If your module takes parameters, here is where you would define
  /// their names and types, and access them to initialize internal
  /// variables. Example as follows:

  analysisthr_.setup(conf.getParameter<std::string>("HistogramFile"));
}

void CSCAFEBAnalyzer::analyze(edm::Event const& e,edm::EventSetup const& iSetup) {
   edm::Handle<CSCWireDigiCollection> wire_digis;

   /// For CSC unpacker
   const char* modtag="cscunpacker";
   e.getByLabel(modtag,"MuonCSCWireDigi",wire_digis);

   /// To get information from the event setup, you must request the "Record"
   /// which contains it and then extract the object you need (HCAL example)

   ///  edm::ESHandle<CaloGeometry> geometry;
   ///  iSetup.get<IdealGeometryRecord>().get(geometry);

   analysisthr_.analyze(*wire_digis);
}

void CSCAFEBAnalyzer::endJob() {
  analysisthr_.done();
}

/// Here are the necessary incantations to declare your module to the
/// framework, so it can be referenced in a cmsRun file.

#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(CSCAFEBAnalyzer)