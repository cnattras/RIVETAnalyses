#!/bin/bash
rm RivetPHENIX_2016_I1394433.so
rivet-build RivetPHENIX_2016_I1394433.so PHENIX_2016_I1394433.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2016_I1394433:cent=GEN -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
