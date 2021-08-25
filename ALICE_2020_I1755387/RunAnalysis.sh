#!/bin/bash
rivet-build RivetALICE_2020_I1755387.so ALICE_2020_I1755387.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a ALICE_2020_I1755387:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
