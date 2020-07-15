#!/bin/bash
SIMULATION_DIR="/home/antoniosilva/RivetWorkdir/GenSim/AuAu_200GeV"
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetPHENIX_2008_I778168.so PHENIX_2008_I778168.cc
rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2008_I778168:cent=GEN:beam=PP -o Rivet_pp.yoda $SIMULATION_DIR/hepmc_AuAu_200GeV_1.hepmc
#rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2008_I778168:cent=GEN:beam=AUAU -o Rivet.yoda $SIMULATION_DIR/hepmc_AuAu_200GeV_1.hepmc
