#!/bin/bash
rivet-buildplugin RivetPHENIX_2006_I0611006.so PHENIX_2006_I0611006.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2006_I0611006:cent=GEN:beam=AUAU200 -o Rivet_Output.yoda hepmc_AuAu_200GeV_1.hepmc.gz
