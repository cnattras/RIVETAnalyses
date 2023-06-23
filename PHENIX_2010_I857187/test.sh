#!/bin/bash
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetPHENIX_2010_I857187.so PHENIX_2010_I857187.cc
#rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2010_I857187:cent=GEN -o Rivet.yoda $PWD/pp/hepmc_pp_200GeV_1.hepmc $PWD/pp/hepmc_pp_200GeV_2.hepmc
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2010_I857187:cent=GEN:beam=AUAU -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
