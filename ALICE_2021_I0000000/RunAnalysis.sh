#!/bin/bash
rivet-build RivetALICE_2021_I0000000.so ALICE_2021_I0000000.cc
export RIVET_ANALYSIS_PATH=$PWD
#rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a

#rivet --pwd -a ALICE_2021_I0000000 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
rivet --pwd -a ALICE_2021_I0000000 -o ALICE_2021_I0000000.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#EXPERIMENT_YEAR_ICODE:cent=GEN:beam=dAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
