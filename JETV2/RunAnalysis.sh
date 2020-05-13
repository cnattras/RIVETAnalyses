#!/bin/bash
rm *.so
rivet-buildplugin RivetALICE_2019_JETV2.so ALICE_2019_JETV2.cc
export RIVET_ANALYSIS_PATH=$PWD
#rivet --pwd -a ALICE_2019_JETV2 --ignore-beams /home/william/JEWEL/jewel-2.2.0/eventfiles/Merge_Test/*/*/merged.hepmc
#rivet --pwd -a ALICE_2019_JETV2 --ignore-beams /media/william/Seagate\ Expansion\ Drive/Jewel_Data/5020.10-30.20M/*/merged.hepmc
rivet --pwd -a ALICE_2019_JETV2 --ignore-beams ~/ANGATYR/*.hepmc
rivet-mkhtml --pwd --errs Rivet.yoda
