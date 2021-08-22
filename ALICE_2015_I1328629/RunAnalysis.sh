#!/bin/bash
rivet-build RivetALICE_2015_I1328629.so ALICE_2015_I1328629.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a ALICE_2015_I1328629 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat

### ALICE_2015_I1328629:cent=GEN:beam=dAU200
