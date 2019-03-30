import FWCore.ParameterSet.Config as cms

mtdTPGenerator = cms.EDProducer("mtdTPGenerator",#"l1pfProducer"              "PF"              "REPR"
                                L1PFObjects      = cms.InputTag("l1pfProducer","PF","REPR"),
                                mtdClusterBarrel = cms.InputTag("mtdClusters", "FTLBarrel", "RECO"),
                                mtdClusterEndcap = cms.InputTag("mtdClusters", "FTLEndcap", "RECO"),
                                FTLBarrel        = cms.InputTag("mix","FTLBarrel","HLT"),
                                FTLEndcap        = cms.InputTag("mix","FTLEndcap","HLT"),
                                recHitBarrel     = cms.InputTag("mtdRecHits","FTLBarrel","RECO"),
                                recHitEndcap     = cms.InputTag("mtdRecHits","FTLEndcap","RECO"),
                                muons            = cms.InputTag('simGmtStage2Digis'),
                                emClusters       = cms.VInputTag(cms.InputTag('pfClustersFromHGC3DClustersEM'), cms.InputTag('pfClustersFromL1EGClusters')),
                                hadClusters      = cms.VInputTag(cms.InputTag('pfClustersFromCombinedCalo:calibrated')),
                                emPtCut          = cms.double(0.5),
                                hadPtCut         = cms.double(1.0),
### Pick only one of these: UseMenu OR UseROI                                
                                useMenu          = cms.bool(True),
                                menuMuPtCut      = cms.double(5),
                                menuEmPtCut      = cms.double(10),
                                menuHcPtCut      = cms.double(15),
                                useROI           = cms.bool(False),
                                roiMuPtCut       = cms.double(2),
                                roiEmPtCut       = cms.double(2.5),
                                roiHcPtCut       = cms.double(5)
                                )
