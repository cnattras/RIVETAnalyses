#!/bin/bash
rm RivetPHENIX_2012_I1116179.so
rivet-build RivetPHENIX_2012_I1116179.so PHENIX_2012_I1116179.cc
export RIVET_ANALYSIS_PATH=$PWD
#rivet --pwd -p  ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2012_I1116179:cent=GEN -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
rivet --pwd -p  ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2012_I1116179:cent=GEN -o Rivet.yoda /phenix/scratch/cen/sampleHepData/AuAu200/hepmcAu.Au.200GeV.1000Events.seed0.pthard0to-1.hepmc
