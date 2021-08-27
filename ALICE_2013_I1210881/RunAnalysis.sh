#!/bin/bash
rivet-build RivetALICE_2013_I1210881.so ALICE_2013_I1210881.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -a ALICE_2013_I1210881 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
