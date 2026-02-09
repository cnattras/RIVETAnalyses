#!/bin/bash
rm *.so
rivet-build RivetPHENIX_2007_I731133.so PHENIX_2007_I731133.cc


date
mkfifo AuAu200GeV.hepmc 

rivet -a PHENIX_2007_I731133 \
  -p pp200GeV.yoda \
  -o AuAu200GeV.yoda \
  AuAu200GeV.hepmc &

pythia8-main144 \
  -c ../RunPYTHIAInDocker/AuAu200GeV.cmnd \
  -c main144HepMC.cmnd \
  -c main144Rivet.cmnd \
  -o AuAu200GeV

pkill -f "rivet -a PHENIX_2007_I731133"
rm -f AuAu200GeV.hepmc
date
