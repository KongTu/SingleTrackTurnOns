import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.Timing = cms.Service("Timing")

process.ana_PbPb = cms.EDAnalyzer('singleTrackAnalyzer',
                      
                		  vertexSrc = cms.string('hiSelectedVertex'),
                		  trackSrc = cms.InputTag('hiGeneralTracks'),             
				  pfCandSrc = cms.untracked.InputTag('particleFlowTmp'),
 
				  doCaloMatched = cms.untracked.bool(True),
				  
				  reso = cms.untracked.double(2.0),#2.0	
				  offlineDCA = cms.untracked.double(3.0),#3.0
				  offlineChi2 = cms.untracked.double(0.15),#0.15
				  offlineptErr = cms.untracked.double(0.1),#0.05
				  offlinenhits = cms.untracked.double(11)#10
					
)

### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

### conditions
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '75X_mcRun2_HeavyIon_v1','')

process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 1000 ) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'root://cmsxrootd.fnal.gov//store/user/qwang/HIHardProbes/HIHardProbes_FullTrackSkim2015_v3/151216_192437/0000/FullTrack_10.root'
))

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
)
process.hltMB = process.hltHM.clone()
process.hltMB.HLTPaths = ['HLT_HIL1MinimumBiasHF1AND_v1']
process.hltMB.andOr = cms.bool(True)#Flase OR
process.hltMB.throw = cms.bool(False)

process.hltHM1 = process.hltHM.clone()
process.hltHM2 = process.hltHM.clone()
process.hltHM3 = process.hltHM.clone()
process.hltHM4 = process.hltHM.clone()
process.hltHM5 = process.hltHM.clone()

##### PbPb HLTSingleTrack ######
process.hltHM1.HLTPaths = ['HLT_HIFullTrack12_L1MinimumBiasHF1_AND_v1']
process.hltHM2.HLTPaths = ['HLT_HIFullTrack18_L1MinimumBiasHF1_AND_v1']
process.hltHM3.HLTPaths = ['HLT_HIFullTrack24_v1']
process.hltHM4.HLTPaths = ['HLT_HIFullTrack34_v1']
process.hltHM5.HLTPaths = ['HLT_HIFullTrack45_v1']

process.hltHM1.andOr = cms.bool(True)
process.hltHM1.throw = cms.bool(False)

process.hltHM2.andOr = cms.bool(True)
process.hltHM2.throw = cms.bool(False)

process.hltHM3.andOr = cms.bool(True)
process.hltHM3.throw = cms.bool(False)

process.hltHM4.andOr = cms.bool(True)
process.hltHM4.throw = cms.bool(False)

process.hltHM5.andOr = cms.bool(True)
process.hltHM5.throw = cms.bool(False)

process.ana_PbPb0 = process.ana_PbPb.clone()
process.ana_PbPb1 = process.ana_PbPb.clone()
process.ana_PbPb2 = process.ana_PbPb.clone()
process.ana_PbPb3 = process.ana_PbPb.clone()
process.ana_PbPb4 = process.ana_PbPb.clone()
process.ana_PbPb5 = process.ana_PbPb.clone()
process.ana_PbPb6 = process.ana_PbPb.clone()
process.ana_PbPb7 = process.ana_PbPb.clone()
process.ana_PbPb8 = process.ana_PbPb.clone()
process.ana_PbPb9 = process.ana_PbPb.clone()

process.d0 = cms.Path( process.hltMB*process.ana_PbPb0 )
process.d1 = cms.Path( process.hltHM1*process.ana_PbPb1 )
process.d2 = cms.Path( process.hltHM2*process.ana_PbPb2 )
process.d3 = cms.Path( process.hltHM3*process.ana_PbPb3 )
process.d4 = cms.Path( process.hltHM4*process.ana_PbPb4 )

process.p0 = cms.Path( process.hltMB*process.hltHM1*process.ana_PbPb5 )
process.p1 = cms.Path( process.hltHM1*process.hltHM2*process.ana_PbPb6 )
process.p2 = cms.Path( process.hltHM2*process.hltHM3*process.ana_PbPb7 )
process.p3 = cms.Path( process.hltHM3*process.hltHM4*process.ana_PbPb8 )
process.p4 = cms.Path( process.hltHM4*process.hltHM5*process.ana_PbPb9 )

process.schedule = cms.Schedule(process.d0,process.d1,process.d2,process.d3,process.d4,process.p0,process.p1,process.p2,process.p3,process.p4)

process.TFileService = cms.Service("TFileService",fileName = cms.string("singletrack.root"))

