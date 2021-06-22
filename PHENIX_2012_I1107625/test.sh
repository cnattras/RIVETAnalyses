#!/bin/bash
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-build RivetPHENIX_2012_I1107625.so PHENIX_2012_I1107625.cc
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2012_I1107625:cent=GEN:beam=AUAU -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
