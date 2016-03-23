import FWCore.ParameterSet.Config as cms
import sys
process = cms.Process("TEST")
runOnMC = True
#use TreeMaker Options
from AMCMSSWInterface.AMTrackProducer.CommandLineParams import CommandLineParams

## MessageLogger
parameters = CommandLineParams()
inputFiles = parameters.value("inputFiles",
    'file:/home/demattia/fdata_demattia/Seb/Test/CMSSW_6_2_0_SLHC27/src/L1Trigger/TrackFindingAM/test/AMFIT_output.root'
)

outputFile = parameters.value("outputFile","test_ntuple.root")
mode = parameters.value("mode","AM")
maxEvents = parameters.value("maxEvents", -1)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

## Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFiles)
)
## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(maxEvents) )


## Write the TTree
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(outputFile)
)


process.ana = cms.EDAnalyzer('TkTriggerParticleAnalzer' ,
    L1TkMuonsInputTag = cms.InputTag("L1TkMuonsMerge",""),
    GenPartInputTag=cms.InputTag("genParticles"),
    TrackPartTag=cms.InputTag("mix","MergedTrackTruth"),
    TTTracksInputTag=cms.InputTag("MergeFITOutput", "AML1Tracks"),
    # TTTracksInputTag=cms.InputTag("AMTrackProducer", "AML1Tracks"),
    inputTagMC = cms.InputTag('TTTrackAssociatorFromPixelDigis', 'AML1Tracks'),
    # inputTagMC = cms.InputTag('TTTrackAssociatorForAM', 'Level1TTTracks'),
    L1TkElectronsInputTag = cms.InputTag("L1TkElectrons","EG"),
    L1TkPhotonsInputTag= cms.InputTag("L1TkPhotons", "EG"),
    ParticleType=cms.int32(13),
)

process.pAna = cms.Path( process.ana )

AMTrackInputTag = cms.InputTag("AMTrackProducer", "Level1TTTracks")
if mode != "AM":AMTrackInputTag = cms.InputTag("TTTracksFromPixelDigis", "Level1TTTracks")

process.schedule = cms.Schedule(process.pAna)
