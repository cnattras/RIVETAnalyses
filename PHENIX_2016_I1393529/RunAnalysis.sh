#!/bin/bash
rm RivetPHENIX_2016_I1393529.so
rivet-build RivetPHENIX_2016_I1393529.so PHENIX_2016_I1393529.cc
export RIVET_ANALYSIS_PATH=$PWD
#rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2016_I1393529:cent=GEN:beam=AUAU -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2016_I1393529:cent=GEN:beam=AUAU -o Rivet_AuAu.yoda /phenix/scratch/cen/sampleHepData/AuAu200/hepmcAu.Au.200GeV.1000Events.seed7.pthard50to-1.hepmc.gz
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2016_I1393529:cent=GEN:beam=PP -o Rivet_pp.yoda /phenix/scratch/cen/sampleHepData/pp200/hepmcd.Au.200GeV.30000Events.seed10.pthard50to-1.hepmc.gz
rivet-merge -O beam -o Rivet.yoda Rivet_AuAu.yoda Rivet_pp.yoda
