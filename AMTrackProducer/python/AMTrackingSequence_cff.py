import FWCore.ParameterSet.Config as cms

from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleGen_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleGenExtra_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleSim_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleDigi_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleL1TrackTrigger_cfi import *
from SLHCL1TrackTriggerSimulations.NTupleTools.ntupleMaker_cfi import *
from AMCMSSWInterface.AMTrackProducer.AMtracks_cfi import *

ntupleGen += ntupleGenExtra
AMTracks=cms.Sequence(ntupleGen * ntupleSim * ntupleDigi * ntupleL1TrackTrigger * AMTrackProducer)
