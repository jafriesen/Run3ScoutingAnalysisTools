import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkSummary.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/320/570/00000/F2AC3A45-9494-E811-9CD1-FA163E390D83.root'
        #'/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/325/101/00000/1EDC0B34-D167-DB4E-81DA-EDDF7C3F0CF6.root',
        #'/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/325/101/00000/A9352733-3357-2D4F-8294-24D40CC3B427.root',
        #'/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/325/101/00000/922DD2E9-5732-084D-B342-3689EBC66AF7.root',
        #'/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/325/022/00000/3EBD7023-C3D9-EF4C-8662-B54DC1755181.root',
        #'/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/325/101/00000/3CACB629-CB62-1A45-BBB5-17E59814F509.root',
        #'/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/325/101/00000/518F962E-2EDE-0E42-A53C-0BE96F18D648.root',
        #'/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/325/101/00000/F4FA257A-30F2-F544-BDE9-5567BA47CA53.root',
        #'/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/325/101/00000/52D96D18-1703-2246-A56D-660A039BF2DB.root',
        #'/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/325/101/00000/D278F2B4-73B9-5348-8571-173169B30909.root',
        '/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/320/570/00000/F2AC3A45-9494-E811-9CD1-FA163E390D83.root'
    )
)

process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("scoutRun2.root")
)

#process.ScoutingFilterPath = cms.Path(process.scoutingFilter)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v37', '')

#L1Info = ["L1_DoubleMu_12_5","L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2"]
L1Info = ["L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"]

process.scoutingTree = cms.EDAnalyzer('ScoutingTreeMakerRun2',
                                      triggerresults   = cms.InputTag("TriggerResults", "", "HLT"),
                                      ReadPrescalesFromFile = cms.bool( False ),
                                      AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                      l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      doL1 = cms.bool( True ),
                                      l1Seeds           = cms.vstring(L1Info),
                                      muons            = cms.InputTag("hltScoutingMuonPackerCalo"),
                                      pfjets           = cms.InputTag("hltScoutingCaloPacker"),
                                      tracks           = cms.InputTag("hltScoutingTrackPacker"),
                                      primaryVertices  = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx"),
                                      displacedVertices  = cms.InputTag("hltScoutingMuonPackerCalo","displacedVtx"),
                                      rho         = cms.InputTag("hltScoutingCaloPacker","rho"),
                                  )

process.p = cms.Path(process.gtStage2Digis+process.scoutingTree)