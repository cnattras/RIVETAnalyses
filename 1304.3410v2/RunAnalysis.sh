#!/bin/bash
rivet-build RivetPHENIX_2013_I1227971.so PHENIX_2013_I1227971.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2013_I1227971:cent=GEN:beam=AUAU200 -o Rivet.yoda hepmc_AuAu_200GeV_1.hepmc
