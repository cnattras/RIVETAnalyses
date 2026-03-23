#!/bin/bash
set -e

rivet-build RivetPHENIX_2002_I585561.so PHENIX_2002_I585561.cc

export RIVET_ANALYSIS_PATH=$PWD
export PYTHONPATH="$(rivet-config --pythonpath):${PYTHONPATH}"
export LD_LIBRARY_PATH="$(rivet-config --libdir):${LD_LIBRARY_PATH}"

pythia8-main144 -n 1000 \
  -c ../RunPYTHIAInDocker/AuAu130GeV.cmnd \
  -c main144HepMC.cmnd \
  -c main144Rivet.cmnd \
  -o AuAu130GeV
