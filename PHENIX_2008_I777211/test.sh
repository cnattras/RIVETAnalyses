#!/bin/bash
#ANALYSIS_DIR=$PWD
#SIMULATION_DIR_200=/home/antoniosilva/RivetWorkdir/GenSim/AuAu_200GeV
#SIMULATION_DIR_PP=/home/antoniosilva/RivetWorkdir/GenSim/pp_200GeV
#cd $ANALYSIS_DIR/Centrality
rm *.so
export RIVET_ANALYSIS_PATH=$PWD


rivet-build RivetPHENIX_2008_I777211.so PHENIX_2008_I777211.cc
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2008_I777211:cent=GEN:beam=AUAU -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat

rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
