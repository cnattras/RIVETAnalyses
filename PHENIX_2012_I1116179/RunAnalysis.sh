#!/bin/bash
rivet-build RivetPHENIX_2012_I1116179.so PHENIX_2012_I1116179.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p  ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2012_I1116179:cent=GEN -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
