#!/bin/bash
rivet-build RivetSTAR_2016_I1420183.so STAR_2016_I1420183.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a STAR_2016_I1420183:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
