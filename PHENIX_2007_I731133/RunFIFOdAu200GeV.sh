#!/bin/bash
rm *.so
rivet-build RivetPHENIX_2007_I731133.so PHENIX_2007_I731133.cc


date
mkfifo dAu200GeV.hepmc 
cat dAu200GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/dAu200GeV.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o dAu200GeV

pkill -f "cat dAu200GeV.hepmc"
rm -f dAu200GeV.hepmc
date
