#!/bin/bash
rm RivetPHENIX_2018_I1672476.so
rivet-build RivetPHENIX_2018_I1672476.so PHENIX_2018_I1672476.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2018_I1672476:cent=GEN -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
