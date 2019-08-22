#this version is for out-of-time photon analysis.
#edited by De-Lin Macive Xiong, Florida State University, 10.22.2018


import FWCore.ParameterSet.Config as cms
from HiggsAnalysis.HiggsTo2photons.hggPhotonIDCuts_cfi import *
#from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
#from CMGTools.External.pujetidproducer_cfi import stdalgos, chsalgos

ggNtuplizer = cms.EDAnalyzer("ggNtuplizer",
                             hggPhotonIDConfiguration = hggPhotonIDCuts,
                             doGenParticles       = cms.bool(True),
                             runOnParticleGun     = cms.bool(False),
                             runOnSherpa          = cms.bool(False),
                             dumpPFPhotons        = cms.bool(True), 
                             dumpJets             = cms.bool(False),
                             dumpAK8Jets          = cms.bool(False),
                             dumpTaus             = cms.bool(False),
                             dumpPDFSystWeight    = cms.bool(False),
                             dumpHFElectrons      = cms.bool(True),
                             development          = cms.bool(False),
                             addFilterInfoAOD     = cms.bool(False),
                             addFilterInfoMINIAOD = cms.bool(False),
                             doNoHFMET            = cms.bool(False),

                             year                 = cms.int32(2017), 

                             trgFilterDeltaPtCut  = cms.double(0.5),
                             trgFilterDeltaRCut   = cms.double(0.3),

                             triggerEvent         = cms.InputTag("slimmedPatTrigger", "", ""),
                             triggerResults       = cms.InputTag("TriggerResults", "", "HLT"),
                             patTriggerResults    = cms.InputTag("TriggerResults", "", "PAT"),
                             #patTriggerResults    = cms.InputTag("TriggerResults", "", "RECO"),
                             genParticleSrc       = cms.InputTag("prunedGenParticles"),
                             generatorLabel       = cms.InputTag("generator"),
                             LHEEventLabel        = cms.InputTag("externalLHEProducer"),
                             newParticles         = cms.vint32(1000006, 1000021, 1000022, 1000024, 1000025, 1000039, 3000001, 3000002, 35),
                             pileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
                             VtxLabel             = cms.InputTag("offlineSlimmedPrimaryVertices"),
                             VtxBSLabel           = cms.InputTag("offlinePrimaryVerticesWithBS"),
                             rhoLabel             = cms.InputTag("fixedGridRhoFastjetAll"),
                             rhoCentralLabel      = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                             pfMETLabel           = cms.InputTag("slimmedMETs"),
                             electronSrc          = cms.InputTag("slimmedElectrons"),
                             #calibelectronSrc     = cms.InputTag("calibratedPatElectrons"),
                             calibelectronSrc     = cms.InputTag("slimmedElectrons"),
                             ####oot####
                             ootPhotons	   	       = cms.InputTag("slimmedOOTPhotons"),
							 #ootPhotons	       = cms.InputTag("patOOTPhotons"),
							 #recoOOTPhotonSrc     = cms.InputTag("reducedEgamma", "reducedOOTPhotons"),
							 #calibOOTphotonSrc    = cms.InputTag("ootPhotons"),
		    			     #calibOOTphotonSrc    = cms.InputTag("patOOTPhotons"),
                             ###########
                             photonSrc            = cms.InputTag("slimmedPhotons"),
                             #calibphotonSrc       = cms.InputTag("calibratedPatPhotons"),
                             calibphotonSrc       = cms.InputTag("slimmedPhotons"),
                             recoPhotonSrc             = cms.InputTag("reducedEgamma", "reducedGedPhotonCores"),
                             muonSrc              = cms.InputTag("slimmedMuons"),
                             gsfTrackSrc          = cms.InputTag("reducedEgamma", "reducedGsfTracks"),
                             ebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
                             eeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits"),
                             esReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedESRecHits"),
                             
                             TrackLabel                = cms.InputTag("generalTracks"),
                             gsfElectronLabel          = cms.InputTag("gsfElectrons"),
                             PFAllCandidates           = cms.InputTag("particleFlow"),
                             #ak4JetSrc                 = cms.InputTag("updatedJets"),
                             ak4JetSrc                 = cms.InputTag("slimmedJets"),
                             #ak4JetSrc                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
                             ak8JetSrc                 = cms.InputTag("slimmedJetsAK8"),
                             #ak8JetSrc                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJECAK8"),
                             #boostedDoubleSVLabel      = cms.InputTag("pfBoostedDoubleSecondaryVertexAK8BJetTags"),
                             tauSrc                    = cms.InputTag("slimmedTaus"),
                             #pfLooseId                 = pfJetIDSelector.clone(),

                             packedPFCands             = cms.InputTag("packedPFCandidates"),
                             elePFClusEcalIsoProducer  = cms.InputTag("electronEcalPFClusterIsolationProducer"),
                             elePFClusHcalIsoProducer  = cms.InputTag("electronHcalPFClusterIsolationProducer"),
                             BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter"),
                             BadPFMuonFilter           = cms.InputTag("BadPFMuonFilter")
                            #####ootphotons#####
                            #ootPhotons                = cms.InputTag("slimmedOOTPhotons")
                            ###################
)

