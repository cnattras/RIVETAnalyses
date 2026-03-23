#!/bin/bash
rivet-build RivetPHENIX_2002_I585561.so PHENIX_2002_I585561.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu130GeV.yoda --ignore-beams -a PHENIX_2002_I585561:cent=GEN:beam=AUAU130 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
