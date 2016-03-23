import FWCore.ParameterSet.Config as cms

AMTrackProducer = cms.EDProducer('AMTrackProducer',
    inputTagStub=cms.InputTag('TTStubsFromPixelDigis', 'StubAccepted'),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_ntuple*_*_*",
    )
)
