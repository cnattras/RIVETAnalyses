#!/bin/bash
ANALYSIS_DIR=$PWD
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-build RivetSTAR_2003_I619063.so STAR_2003_I619063.cc
rivet --pwd -p ../Centralities/Calibration/calibration_STAR_AuAu130GeV.yoda -a STAR_2003_I619063:cent=GEN:beam=AUAU130 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
rivet-merge --pwd -p ../Centralities/Calibration/calibration_AuAu_130GeV_STAR.yoda -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
