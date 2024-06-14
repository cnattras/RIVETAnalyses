#!/bin/bash
rm RivetPHENIX_2016_I1394433.so

rivet-build RivetPHENIX_2016_I1394433.so PHENIX_2016_I1394433.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2016_I1394433:cent=GEN -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#rivet --pwd -p Centralities/Calibration/calibration_PHENIX_AuAu62GeV.yoda -a PHENIX_2016_I1394433:cent=GEN -o Rivet_AUAU_62.yoda ../testfiles/AuAu_62.hepmc
#rivet --pwd -p Centralities/Calibration/calibration_PHENIX_AuAu39GeV.yoda -a PHENIX_2016_I1394433:cent=GEN -o Rivet_AUAU_39.yoda ../testfiles/AuAu_39.hepmc







