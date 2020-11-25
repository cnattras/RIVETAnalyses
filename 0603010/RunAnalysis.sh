#!/bin/bash
rivet-build RivetPHENIX_2006_I711951.so PHENIX_2006_I711951.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2006_I711951:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
