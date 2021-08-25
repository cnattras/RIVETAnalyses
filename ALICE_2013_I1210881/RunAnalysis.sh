#!/bin/bash
rivet-build RivetALICE_2013_I1210881.so ALICE_2013_I1210881.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a ALICE_2013_I1210881:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
