import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.source = cms.Source("PoolSource",
                    fileNames = cms.untracked.vstring(
                        'root://cms-xrd-global.cern.ch//store/data/Run2022D/Muon/MINIAOD/PromptReco-v1/000/357/542/00000/050a5ab6-44bb-46b4-b624-123b55834e60.root',
                    ),
                    secondaryFileNames = cms.untracked.vstring(),
#                     lumisToProcess = cms.untracked.VLuminosityBlockRange('258158:1-258158:1786'),

)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.GlobalTag.globaltag = 'auto:phase1_2021_realistic'
process.GlobalTag.globaltag = 'auto:run3_hlt'

from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'auto:run3_hlt')

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')


process.tauNtuples =cms.EDAnalyzer("TauNtuples",
#                        offlineVtx               = cms.InputTag("offlinePrimaryVertices"),
#                        offlineMuons             = cms.InputTag("muons"),
#                        
                       triggerResult            = cms.untracked.InputTag("TriggerResults::MYHLT"),
                       triggerSummary           = cms.untracked.InputTag("hltTriggerSummaryAOD::MYHLT"),
                       triggerResultTag         = cms.untracked.InputTag("TriggerResults::HLT"),
                       triggerSummaryTag        = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
                       triggerObjMiniTag        = cms.untracked.InputTag("slimmedPatTrigger::RECO"),
                       

                       hltMuCandidates          = cms.untracked.InputTag("hltIterL3DisplacedMuonCandidates"), 

                       hltTauCandidates         = cms.untracked.InputTag("hltHpsPFTauProducerDispl", "", "MYHLT"), 
                       tauIP                    = cms.untracked.InputTag("hltHpsPFTauTransverseImpactParameters"), 
                       tauIso                   = cms.untracked.InputTag("hltHpsDisplPFTauMediumAbsOrRelChargedIsolationDiscriminator"), 
                       tauIsoValue              = cms.untracked.InputTag("hltHpsDisplPFTauMediumAbsoluteChargedIsolationValue"), 

                       hltTauCandidatesNoDispl  = cms.untracked.InputTag("hltHpsPFTauProducer", "", "MYHLT"), 
                       tauIPNoDispl             = cms.untracked.InputTag("hltHpsPFTauTransverseImpactParametersNoDispl"), 
                       tauIsoNoDispl            = cms.untracked.InputTag("hltHpsPFTauMediumAbsOrRelChargedIsolationDiscriminator"), 
                       tauIsoValueNoDispl       = cms.untracked.InputTag("hltHpsPFTauMediumAbsoluteChargedIsolationValue"), 

                       L2isoJetTag              = cms.untracked.InputTag("hltL2TauPixelIsoTagProducerL1TauSeededGlob"), 
                       onlineVertices           = cms.untracked.InputTag("hltTrimmedPixelVertices"), 
                       tauPVertices             = cms.untracked.InputTag("hltHpsPFTauPrimaryVertexProducer", "PFTauPrimaryVertices", "MYHLT"), 
                       pixelTracks              = cms.untracked.InputTag("hltPixelTracks"), 
                       allTracks                = cms.untracked.InputTag("hltPFMuonMergingForDisplTau"), 

                       L1MuonCandidates         = cms.untracked.InputTag("gmtStage2Digis", "Muon"), 
                       L1TauCandidates          = cms.untracked.InputTag("caloStage2Digis", "Tau"), 

                       lumiScalerTag            = cms.untracked.InputTag("scalersRawToDigi"),
                       puInfoTag                = cms.untracked.InputTag("addPileupInfo"),
                       beamspot                 = cms.untracked.InputTag("hltOnlineBeamSpot"),
                       )   
                       
process.mypath  = cms.Path(process.tauNtuples)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("tauNtuple.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))   


process.MessageLogger = cms.Service("MessageLogger",
   destinations   = cms.untracked.vstring('cerr'),
   cerr           = cms.untracked.PSet(
       threshold      = cms.untracked.string('ERROR'),
   ),
#    debugModules  = cms.untracked.vstring('*')
)

