#!/bin/bash
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-build RivetPHENIX_2013_I1207323.so PHENIX_2013_I1207323.cc
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2013_I1207323:cent=GEN -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
