#ifndef CSC_DBNOISEMATRIX_SRC_IMPL_H
#define CSC_DBNOISEMATRIX_SRC_IMPL_H

#include <vector>
#include <string>
#include <iostream>
#include <typeinfo>

#include "CondCore/PopCon/interface/PopConSourceHandler.h"
#include "CondCore/PopCon/interface/LogReader.h"
#include "CondFormats/CSCObjects/interface/CSCobject.h"
#include "CondFormats/DataRecord/interface/CSCDBNoiseMatrixRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "CSCNoiseMatrixDBConditions.h"

namespace popcon
{
  class CSCDBNoiseMatrixImpl : public popcon::PopConSourceHandler<CSCDBNoiseMatrix>
    {
      
    public:
      void getNewObjects();
      std::string id() const { return m_name;}
      ~CSCDBNoiseMatrixImpl(); 
      CSCDBNoiseMatrixImpl(const edm::ParameterSet& pset);
       
    private:
      std::string m_pop_connect; //connect string to popcon metaschema
      std::string m_name;
      std::string m_cs;
      const CSCDBNoiseMatrix * mymatrix;
      LogReader* lgrdr;
    };
}
#endif
