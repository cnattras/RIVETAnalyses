#!/bin/bash
rivet-build RivetALICE_2019_I1693308.so ALICE_2019_I1693308.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a ALICE_2019_I1693308:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
