#!/bin/bash
rivet-build RivetSPHENIX_2025_I2907573.so SPHENIX_2025_I2907573.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/AuAu200GeV.yoda -a SPHENIX_2025_I2907573:cent=GEN:beam=AUAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#rivet --pwd -p ../Centralities/Calibration/AuAu200GeV.yoda -a SPHENIX_2025_I2907573:cent=IMP:beam=AUAU130 -o Rivet.IMP.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#rivet --pwd -p ../Centralities/Calibration/AuAu200GeV.yoda -a SPHENIX_2025_I2907573:cent=GEN:beam=AUAU130:fixedcentrality=5 -o Rivet.FIXED.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
