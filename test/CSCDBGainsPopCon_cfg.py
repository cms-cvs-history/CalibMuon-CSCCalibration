# The following comments couldn't be translated into the new config version:

# eg to write payload to the oracle database 
#   replace CondDBCommon.connect = "oracle://cms_orcoff_int2r/CMS_COND_CSC"
# Database output service

import FWCore.ParameterSet.Config as cms

process = cms.Process("ProcessOne")
#PopCon config
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.CondDBCommon.connect = cms.string("sqlite_file:DBGains.db")

process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        )
    ),
    destinations = cms.untracked.vstring('cout')
)

process.source = cms.Source("EmptyIOVSource",
    lastRun = cms.untracked.uint32(1),
    timetype = cms.string('runnumber'),
    #change the firstRun if you want a different IOV
    firstRun = cms.untracked.uint32(1),
    interval = cms.uint32(1)
)

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDBCommon,
    logconnect = cms.untracked.string('sqlite_file:gainslog.db'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('CSCDBGainsRcd'),
        tag = cms.string('CSCDBGains_new_popcon')
    ))
)

process.WriteGainsWithPopCon = cms.EDAnalyzer("CSCGainsPopConAnalyzer",
    SinceAppendMode = cms.bool(True),
    record = cms.string('CSCDBGainsRcd'),
    loggingOn = cms.untracked.bool(True),
    Source = cms.PSet(

    )
)

process.p = cms.Path(process.WriteGainsWithPopCon)

