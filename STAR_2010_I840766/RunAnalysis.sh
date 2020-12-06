#!/bin/bash
rm RivetSTAR_2010_I840766.so
rivet-build RivetSTAR_2010_I840766.so STAR_2010_I840766.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_AuAu_200GeV_STAR.yoda -a STAR_2010_I840766:cent=GEN:beam=pp -o Rivet_pp.yoda hepmc_pp_200GeV_1.hepmc
rivet --pwd -p ../Centralities/Calibration/calibration_AuAu_200GeV_STAR.yoda -a STAR_2010_I840766:cent=GEN:beam=dAu -o Rivet_dAu.yoda PYTHIAAuAuFileSMALLTEST.dat 
rivet-merge -O beam -O cent -o Rivet_final.yoda Rivet_dAu.yoda Rivet_pp.yoda
#rivet-mkhtml --pwd Rivet_pp.yoda