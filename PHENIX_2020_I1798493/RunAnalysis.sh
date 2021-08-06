#!/bin/bash
rivet-build RivetPHENIX_2020_I1798493.so PHENIX_2020_I1798493.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2020_I1798493:cent=GEN:beam=AUAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
