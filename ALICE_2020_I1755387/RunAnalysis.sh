#!/bin/bash
rivet-build RivetALICE_2020_I1755387.so ALICE_2020_I1755387.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -a ALICE_2020_I1755387 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
