#!/bin/bash
rivet-build RivetPHENIX_2008_I776624.so PHENIX_2008_I776624.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2008_I776624:cent=GEN: -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat