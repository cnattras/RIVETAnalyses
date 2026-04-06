#!/bin/bash
set -e

rivet-build RivetPHENIX_2002_I585561.so PHENIX_2002_I585561.cc

export RIVET_ANALYSIS_PATH=$PWD
export PYTHONPATH="$(rivet-config --pythonpath):${PYTHONPATH}"
export LD_LIBRARY_PATH="$(rivet-config --libdir):${LD_LIBRARY_PATH}"

rm -f AuAu130GeV.hepmc
mkfifo AuAu130GeV.hepmc
cat AuAu130GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/AuAu130GeV.cmnd -c ../RunPYTHIAInDocker/main144HepMC.cmnd -c main144Rivet.cmnd -o AuAu130GeV

pkill -f "cat AuAu130GeV.hepmc" || true
rm -f AuAu130GeV.hepmc
