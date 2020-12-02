#!/bin/bash
rivet-build RivetSTAR_2009_I793126.so STAR_2009_I793126.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_AuAu_200GeV_STAR.yoda -a STAR_2009_I793126:cent=GEN:beam=AUAU130 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
