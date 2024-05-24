#!/bin/bash
rm PHENIX_2018_I1672859.so
rivet-build RivetPHENIX_2018_I1672859.so PHENIX_2018_I1672859.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2018_I1672859:cent=GEN:beam=CUAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
