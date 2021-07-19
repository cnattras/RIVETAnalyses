#!/bin/bash
rivet-build RivetPHENIX_2004_I623413.so PHENIX_2004_I623413.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2004_I623413:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
