# Analyzer for CMSSW AM simulation

This package contains an analyzer that accepts as input a file produced by running the AM+LTF in CMSSW. The output is a root tree with reconstructed track parameters and chi2/ndf as well as parametrs for the tracking particles with $p_{T}$ > 3 GeV/c.

The analyzer can be run with the cfg under AMTrackProducer/test.
