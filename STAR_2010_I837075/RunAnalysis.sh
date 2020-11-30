#!/bin/bash
rm STAR_2010_I837075.so
rivet-build RivetSTAR_2010_I837075.so STAR_2010_I837075.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_AuAu_200GeV_STAR.yoda -a STAR_2010_I837075:cent=GEN:beam=AuAu200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat

