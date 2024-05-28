#!/bin/bash
rm RivetPHENIX_2010_I857187.so
rivet-build RivetPHENIX_2010_I857187.so PHENIX_2010_I857187.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2010_I857187:cent=GEN -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat

