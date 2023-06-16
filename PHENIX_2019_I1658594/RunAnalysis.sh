#!/bin/bash
rivet-build RivetPHENIX_2019_I1658594.so PHENIX_2019_I1658594.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2019_I1658594:cent=GEN:beam=AUAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
