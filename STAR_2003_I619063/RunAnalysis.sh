#!/bin/bash
ANALYSIS_DIR=$PWD
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetSTAR_2003_I619063.so STAR_2003_I619063.cc
rivet --pwd -p $ANALYSIS_DIR/Centrality/calibration_AuAu_130GeV_STAR.yoda -a STAR_2003_I619063:cent=GEN:beam=AUAU130 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
rivet-merge --pwd -p $ANALYSIS_DIR/Centrality/calibration_AuAu_130GeV_STAR.yoda -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
