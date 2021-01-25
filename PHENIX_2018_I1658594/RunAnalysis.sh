#!/bin/bash
rivet-build RivetPHENIX_2018_I1658594.so PHENIX_2018_I1658594.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2018_I1658594:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
