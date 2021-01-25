#!/bin/bash
rivet-build RivetBRAHMS_2007_I742956.so BRAHMS_2007_I742956.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a BRAHMS_2007_I742956:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
