#!/bin/bash
rm PHENIX_2018_I1672859.so
rivet-build RivetPHENIX_2018_I1672859.so PHENIX_2018_I1672859.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2018_I1672859:cent=GEN:beam=CUAU200 -o Rivet_CuAu.yoda /phenix/scratch/cen/sampleHepData/CuAu200/hepmcCu.Au.200GeV.1000Events.seed0.pthard30to40.hepmc.gz
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2018_I1672859:cent=GEN:beam=PP200 -o Rivet_pp.yoda /phenix/scratch/cen/sampleHepData/pp200/hepmcd.Au.200GeV.30000Events.seed11.pthard30to40.hepmc.gz
