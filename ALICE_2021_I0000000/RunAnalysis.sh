#!/bin/bash
rivet-build RivetEXPERIMENT_YEAR_ICODE.so EXPERIMENT_YEAR_ICODE.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a EXPERIMENT_YEAR_ICODE:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
