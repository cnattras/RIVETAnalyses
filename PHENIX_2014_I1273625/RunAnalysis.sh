#!/bin/bash
rivet-build RivetPHENIX_2014_I1273625.so PHENIX_2014_I1273625.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/AuAu200GeV.yoda -a PHENIX_2014_I1273625:cent=GEN:beam=AUAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat

# rivet --pwd -p ../Centralities/Calibration/AuAu200GeV.yoda -a PHENIX_2014_I1273625:cent=IMP:beam=AUAU200 -o Rivet.IMP.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat

# rivet --pwd -p ../Centralities/Calibration/AuAu200GeV.yoda -a PHENIX_2014_I1273625:cent=GEN:beam=AUAU200:fixedcentrality=5 -o Rivet.fixed.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat