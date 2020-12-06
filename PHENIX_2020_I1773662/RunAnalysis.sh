#!/bin/bash
rivet-build RivetPHENIX_2020_I1773662.so PHENIX_2020_I1773662.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2020_I1773662:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
