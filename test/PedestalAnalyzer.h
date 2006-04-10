/** 
 * Demo analyzer for reading digis
 * author A.Tumanov 2/22/06 
 *   
 */

#include <iostream>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"


class PedestalAnalyzer : public edm::EDAnalyzer {
public:
  explicit PedestalAnalyzer(edm::ParameterSet const& conf);
  virtual void analyze(edm::Event const& e, edm::EventSetup const& iSetup);
 
  //constants declaration
#define CHAMBERS 5
#define LAYERS 6
#define STRIPS 80
  
  //some variable declaration
  int evt;
  std::vector<int> newadc;
  std::vector<int> adc;
  float pedMean,time,max,max1;
  int pedSum, strip;
  
  //definition of arrays
  double arrayOfPed[CHAMBERS][LAYERS][STRIPS];
  double arrayOfPedSquare[CHAMBERS][LAYERS][STRIPS];
  double arrayPed[CHAMBERS][LAYERS][STRIPS];
  double arrayPeak[CHAMBERS][LAYERS][STRIPS];
  double arrayOfPeak[CHAMBERS][LAYERS][STRIPS];
  double arrayOfPeakSquare[CHAMBERS][LAYERS][STRIPS];
  double arraySumFive[CHAMBERS][LAYERS][STRIPS];
  
  double newPed[480];
  double newRMS[480];
  double newPeakRMS[480];
  double newPeak[480];
  double newSumFive[480];



  //virtual void endJob();
private:
  // variables persistent across events should be declared here.
  //
  int eventNumber;
};
