#!/bin/bash
rm *.so
rivet-build RivetPHENIX_2008_I778168.so PHENIX_2008_I778168.cc


date
mkfifo AuAu200GeV.hepmc 
cat AuAu200GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/AuAu200GeV.cmnd -c main144HepMC.cmnd -c main144Rivet_AuAu.cmnd -o AuAu200GeV

pkill -f "cat AuAu200GeV.hepmc"
rm -f AuAu200GeV.hepmc
date
