#!/bin/bash
ANALYSIS_DIR=/home/antoniosilva/RivetWorkDir/BackgroundRho/ResponseMatrix
PYTHIA_EX_DIR=/media/antoniosilva/Backup/HEPMC_PbPb_MB
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetResponseMtrix_2020_Example.so ResponseMtrix_2020_Example.cc
rivet --pwd -p $ANALYSIS_DIR/calibration.yoda -a ResponseMtrix_2020_Example:cent=GEN -o RM_output.yoda --ignore-beams $PYTHIA_EX_DIR/hepMC_files.hepmc
