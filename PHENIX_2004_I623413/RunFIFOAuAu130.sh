#!/bin/bash
rm *.so
rivet-build RivetPHENIX_2004_I623413.so PHENIX_2004_I623413.cc


date
mkfifo AuAu130GeV.hepmc 
cat AuAu130GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/AuAu130GeV.cmnd -c main144HepMC.cmnd -c main144RivetAuAu130.cmnd -o AuAu130GeV

pkill -f "cat AuAu130GeV.hepmc"
rm -f AuAu130GeV.hepmc
date
