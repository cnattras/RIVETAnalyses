#!/bin/bash
rm RivetPHENIX_2019_I1672133.so
rivet-build RivetPHENIX_2019_I1672133.so PHENIX_2019_I1672133.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2019_I1672133:cent=GEN -o Rivet.yoda hepmc_file.hepmc
