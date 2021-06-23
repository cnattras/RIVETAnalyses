#!/bin/bash
rivet-build RivetALICE_2019_I1693308.so ALICE_2019_I1693308.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -a ALICE_2019_I1693308 -o Rivet.yoda ../../../../../../../../eos/user/j/jpiel/pythia8.pp.inel.70000.hepmc
