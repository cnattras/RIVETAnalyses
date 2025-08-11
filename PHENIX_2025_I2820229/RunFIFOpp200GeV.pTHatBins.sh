#!/bin/bash
export RIVET_ANALYSIS_PATH=$PWD


rivet-build RivetPHENIX_2025_I2820229.so PHENIX_2025_I2820229.cc


mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &
pythia8-main144 -c pp200GeV.1.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o pp200GeV.1
pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc


mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &
pythia8-main144 -c pp200GeV.2.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o pp200GeV.2
pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc

mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &
pythia8-main144 -c pp200GeV.3.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o pp200GeV.3
pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc

mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &
pythia8-main144 -c pp200GeV.4.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o pp200GeV.4
pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc

mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &
pythia8-main144 -c pp200GeV.5.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o pp200GeV.5
pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc

mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &
pythia8-main144 -c pp200GeV.6.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o pp200GeV.6
pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc

mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &
pythia8-main144 -c pp200GeV.7.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o pp200GeV.7
pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc

mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &
pythia8-main144 -c pp200GeV.8.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o pp200GeV.8
pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc