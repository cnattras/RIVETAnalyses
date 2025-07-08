#!/bin/bash

mkfifo UU193GeV.hepmc 
cat UU193GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/UU193GeV.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o UU193GeV

pkill -f "cat UU193GeV.hepmc"
rm -f UU193GeV.hepmc
