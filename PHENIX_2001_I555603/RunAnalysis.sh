#!/bin/bash
rivet-build RivetPHENIX_2001_I555603.so PHENIX_2001_I555603.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/AuAu130GeV.yoda -a PHENIX_2001_I555603:cent=GEN:beam=AUAU130 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#rivet --pwd -p ../Centralities/Calibration/AuAu130GeV.yoda -a PHENIX_2001_I555603:cent=IMP:beam=AUAU130 -o Rivet.IMP.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#rivet --pwd -p ../Centralities/Calibration/AuAu130GeV.yoda -a PHENIX_2001_I555603:cent=GEN:beam=AUAU130:fixedcentrality=5 -o Rivet.FIXED.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
