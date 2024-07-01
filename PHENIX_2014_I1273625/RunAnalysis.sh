#!/bin/bash
rivet-build RivetPHENIX_2014_I1273625.so PHENIX_2014_I1273625.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2014_I1273625:cent=GEN:beam=AUAU200 -a PHENIX_2014_I1273625:cent=GEN:beam=AUAU130 -a PHENIX_2014_I1273625:cent=GEN:beam=AUAU62 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
