#!/bin/bash
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetPHENIX_2013_I1207323.so PHENIX_2013_I1207323.cc
rivet --pwd -p calibration.yoda -a PHENIX_2013_I1207323:cent=GEN -o Rivet.yoda $PWD/PYTHIAAuAuFileSMALLTEST.dat
