#!/bin/bash
rivet-build RivetALICE_2024_IPAT.so ALICE_2024_IPAT.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a ALICE_2024_IPAT -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
