#!/bin/bash


rivet-build RivetSPHENIX_2025_I2907537.so SPHENIX_2025_I2907537.cc
export RIVET_ANALYSIS_PATH=$PWD

mkfifo AuAu200GeV.hepmc 
cat AuAu200GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/AuAu200GeV.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o AuAu200GeV

pkill -f "cat AuAu200GeV.hepmc"
rm -f AuAu200GeV.hepmc
