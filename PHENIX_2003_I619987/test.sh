#!/bin/bash
rm RivetPHENIX_2003_I619987.so
rivet-build RivetPHENIX_2003_I619987.so PHENIX_2003_I619987.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2003_I619987:cent=GEN:beam=AUAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat

rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
