#!/bin/bash
ANALYSIS_DIR=$PWD
#SIMULATION_DIR_200=/home/antoniosilva/RivetWorkdir/GenSim/AuAu_200GeV
SIMULATION_DIR_PP=/home/antoniosilva/RivetWorkdir/GenSim/pp_200GeV
#cd $ANALYSIS_DIR/Centrality
#rm *.so
#export RIVET_ANALYSIS_PATH=$PWD
#rivet-buildplugin RivetRHIC_2019_CentralityCalibration.so RHIC_2019_CentralityCalibration.cc
#rivet --pwd -a RHIC_2019_CentralityCalibration:exp=STAR -o calibration_AuAu_130GeV_STAR.yoda --ignore-beams $SIMULATION_DIR/hepmc_AuAu_130GeV_1.hepmc $SIMULATION_DIR/hepmc_AuAu_130GeV_2.hepmc $SIMULATION_DIR/hepmc_AuAu_130GeV_3.hepmc $SIMULATION_DIR/hepmc_AuAu_130GeV_4.hepmc $SIMULATION_DIR/hepmc_AuAu_130GeV_5.hepmc
#cd ..
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetPHENIX_2008_I777211.so PHENIX_2008_I777211.cc
#rivet --pwd -p $ANALYSIS_DIR/Centrality/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2008_I777211:cent=GEN:beam=AUAU200 -o Angantyr_AuAu_200GeV.yoda $SIMULATION_DIR_200/hepmc_AuAu_200GeV_1.hepmc $SIMULATION_DIR_200/hepmc_AuAu_200GeV_2.hepmc $SIMULATION_DIR_200/hepmc_AuAu_200GeV_3.hepmc $SIMULATION_DIR_200/hepmc_AuAu_200GeV_4.hepmc $SIMULATION_DIR_200/hepmc_AuAu_200GeV_5.hepmc
rivet --pwd -a PHENIX_2008_I777211:cent=GEN:beam=PP -o Angantyr_pp_200GeV.yoda $SIMULATION_DIR_PP/hepmc_pp_200GeV_1.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_2.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_3.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_4.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_5.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_6.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_7.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_8.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_9.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_10.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_11.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_12.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_13.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_14.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_15.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_16.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_17.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_18.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_19.hepmc $SIMULATION_DIR_PP/hepmc_pp_200GeV_20.hepmc
#rivet-mkhtml --pwd Rivet.yoda 
