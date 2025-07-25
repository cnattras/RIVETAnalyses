#!/bin/bash
export RIVET_ANALYSIS_PATH=$PWD


rivet-build RivetPHENIX_2025_I2820229.so PHENIX_2025_I2820229.cc


mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &

pythia8-main144 -c pp200GeV.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o pp200GeV

pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc
