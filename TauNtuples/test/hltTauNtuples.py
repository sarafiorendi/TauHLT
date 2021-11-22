import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.source = cms.Source("PoolSource",
                    fileNames = cms.untracked.vstring(
#                       '/store/group/phys_bphys/fiorendi/p5prime/displTaus/Staus_M_500_100mm_14TeV_Run3MC/crab_hlt_gmsb_100mm_v26_patatrack_iter4_onL1CandsPt24_remove130/211118_151855/0000/outputHLT_1.root',
#                       'file:/eos/cms/store/group/phys_bphys/fiorendi/p5prime/displTaus/1cf292de-a744-4511-818f-ff9f09961cd9.root',
                      'file:/afs/cern.ch/work/f/fiorendi/private/displacedTaus/hlt/CMSSW_12_1_0/src/outputHLT_2.root',
                    ),
                    secondaryFileNames = cms.untracked.vstring(),
#                     lumisToProcess = cms.untracked.VLuminosityBlockRange('258158:1-258158:1786'),

)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'auto:phase1_2021_realistic'

# process.GlobalTag.globaltag = cms.string('94X_dataRun2_v11') ## for 2017
from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'auto:phase1_2021_realistic')

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')


process.tauNtuples =cms.EDAnalyzer("TauNtuples",
#                        offlineVtx               = cms.InputTag("offlinePrimaryVertices"),
#                        offlineMuons             = cms.InputTag("muons"),
#                        
                       triggerResult            = cms.untracked.InputTag("TriggerResults::MYHLT"),
                       triggerSummary           = cms.untracked.InputTag("hltTriggerSummaryAOD::MYHLT"),
#                        tagTriggerResult         = cms.untracked.InputTag("TriggerResults::HLT"),
#                        tagTriggerSummary        = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
#                        
                       hltMuCandidates             = cms.untracked.InputTag("hltL3MuonCandidates"), 
                       hltTauCandidates            = cms.untracked.InputTag("hltHpsPFTauProducerDispl", "", "MYHLT"), 
#                        hltTauCandidates            = cms.untracked.InputTag("hltHpsSelectedPFTausTrackPt1NoIsolationGlobDispl", "", "MYHLT"), 
                       tauIP                       = cms.untracked.InputTag("hltHpsPFTauTransverseImpactParameters"), 
                       tauIso                      = cms.untracked.InputTag("hltHpsDisplPFTauMediumAbsOrRelChargedIsolationDiscriminator"), 
                       tauIsoValue                 = cms.untracked.InputTag("hltHpsDisplPFTauMediumAbsOrRelChargedIsolationValue"), 
#                        L2Candidates             = cms.untracked.InputTag("hltL2MuonCandidates"), 
                       L1MuonCandidates            = cms.untracked.InputTag("gmtStage2Digis", "Muon"), 
                       L1TauCandidates             = cms.untracked.InputTag("caloStage2Digis", "Tau"), 

#                        RhoCorrectionOnline      = cms.untracked.InputTag("hltFixedGridRhoFastjetAllCaloForMuons"), # for now, same for tag and probe muons

#                        lumiScalerTag            = cms.untracked.InputTag("scalersRawToDigi"),
#                        puInfoTag                = cms.untracked.InputTag("addPileupInfo"),
#                        genParticlesTag          = cms.untracked.InputTag("genParticles"),
#                        doOffline                = cms.untracked.bool(True)
                       )   
                       
process.mypath  = cms.Path(process.tauNtuples)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("tauNtuple.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))   


process.MessageLogger = cms.Service("MessageLogger",
   destinations   = cms.untracked.vstring('cerr'),
   cerr           = cms.untracked.PSet(
       threshold      = cms.untracked.string('ERROR'),
   ),
#    debugModules  = cms.untracked.vstring('*')
)

