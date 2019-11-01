#!/bin/bash
ANALYSIS_DIR=PATH_TO_ANALYSIS/1110.5800
PYTHIA_EX_DIR=PATH_TO_PYTHIA/pythia8243/examples
cd $ANALYSIS_DIR/CentCalibration
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetRHIC_2019_CentralityCalibration.so RHIC_2019_CentralityCalibration.cc
rivet --pwd -a RHIC_2019_CentralityCalibration:exp=STAR -o calibration.yoda --ignore-beams $PYTHIA_EX_DIR/hepmc_AuAu_MB_1.hepmc
rivet-mkhtml --pwd calibration.yoda 
cd ..
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetSTAR_2012_I943192.so STAR_2012_I943192.cc
rivet --pwd -p $ANALYSIS_DIR/CentCalibration/calibration.yoda -a STAR_2012_I943192:cent=GEN $PYTHIA_EX_DIR/hepmc_AuAu_MB_1.hepmc
rivet-mkhtml --pwd Rivet.yoda 
