#!/bin/bash
rm RivetPHENIX_2011_I900703.so
rivet-build RivetPHENIX_2011_I900703.so PHENIX_2011_I900703.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2011_I900703:cent=GEN -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
