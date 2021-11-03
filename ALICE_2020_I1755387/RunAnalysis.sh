#!/bin/bash
rivet-build RivetALICE_2020_I1755387.so ALICE_2020_I1755387.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_ALICE_PbPb2760GeV.yoda -a ALICE_2020_I1755387:cent=GEN -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
