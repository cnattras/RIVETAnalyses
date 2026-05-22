#!/bin/bash
rm *.so
rivet-build RivetPHENIX_2008_I778168.so PHENIX_2008_I778168.cc


date
mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/pp200GeV.cmnd -c main144HepMC.cmnd -c main144Rivet_pp.cmnd -o pp200GeV

pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc
date
