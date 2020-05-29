#!/bin/bash
SIMULATION_DIR=/media/antoniosilva/Backup/AuAu_62GeV
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetRHIC_2019_CentralityCalibration.so RHIC_2019_CentralityCalibration.cc
FILES=""
for i in {1..10}
    do
        FILES="$FILES $SIMULATION_DIR/hepmc_AuAu_62GeV_$i.hepmc"
    done
rivet --pwd -a RHIC_2019_CentralityCalibration:exp=PHENIX -o calibration_AuAu_62GeV_PHENIX.yoda --ignore-beams $FILES 
