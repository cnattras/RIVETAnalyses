#!/bin/bash
rm RivetPHENIX_2009_I815217.so
rivet-build RivetPHENIX_2009_I815217.so PHENIX_2009_I815217.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2009_I815217:cent=GEN -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
