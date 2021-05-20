#!/bin/bash
ANALYSIS_DIR=$PWD
SIMULATION_DIR_39=/media/antoniosilva/Backup/AuAu_39GeV
SIMULATION_DIR_62=/media/antoniosilva/Backup/AuAu_62GeV
SIMULATION_DIR_PP39=/media/antoniosilva/B848-CB321/pp_39GeV
SIMULATION_DIR_PP62=/media/antoniosilva/B848-CB321/pp_62GeV
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-build RivetPHENIX_2012_I1107625.so PHENIX_2012_I1107625.cc
#for i in {1..40}
 #   do
  #      FILES39="$FILES39 $SIMULATION_DIR_39/hepmc_AuAu_39GeV_$i.hepmc"
   #     FILES62="$FILES62 $SIMULATION_DIR_62/hepmc_AuAu_62GeV_$i.hepmc"
    #    FILESPP39="$FILESPP39 $SIMULATION_DIR_PP39/hepmc_pp_39GeV_$i.hepmc"
     #   FILESPP62="$FILESPP62 $SIMULATION_DIR_PP62/hepmc_pp_62GeV_$i.hepmc"
  #  done
#rivet --pwd -p $ANALYSIS_DIR/Centrality/calibration_AuAu_39GeV_PHENIX.yoda -a PHENIX_2012_I1107625:cent=GEN:beam=AUAU39 -o Angantyr_AuAu_39GeV.yoda $FILES39
#rivet --pwd -p $ANALYSIS_DIR/Centrality/calibration_AuAu_62GeV_PHENIX.yoda -a PHENIX_2012_I1107625:cent=GEN:beam=AUAU62 -o Angantyr_AuAu_62GeV.yoda $FILES62
#rivet --pwd -a PHENIX_2012_I1107625:cent=GEN:beam=PP39 -o Angantyr_pp_39GeV.yoda $FILESPP39
#rivet --pwd -a PHENIX_2012_I1107625:cent=GEN:beam=PP62 -o Angantyr_pp_62GeV.yoda $FILESPP62

rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2012_I1107625:cent=GEN:beam=AUAU -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat