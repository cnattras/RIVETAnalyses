#!/bin/bash
rm *.so
rivet-build RivetPHENIX_2007_I731133.so PHENIX_2007_I731133.cc


date
mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/pp200GeV.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o pp200GeV

pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc
date
