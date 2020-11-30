#!/bin/bash
rivet-build RivetSTAR_2020_I1771348.so STAR_2020_I1771348.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_AuAu_200GeV_STAR.yoda -a STAR_2020_I1771348:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
