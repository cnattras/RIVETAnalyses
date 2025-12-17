#!/bin/bash
rivet-build RivetPHENIX_2007_I731133.so PHENIX_2007_I731133.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/AuAu200GeV.yoda -a PHENIX_2007_I731133:cent=GEN:beam=AUAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
