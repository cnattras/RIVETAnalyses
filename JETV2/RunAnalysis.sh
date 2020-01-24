#!/bin/bash
rm *.so
rivet-buildplugin RivetALICE_2019_JETV2.so ALICE_2019_JETV2.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -a ALICE_2019_JETV2 --ignore-beams /home/william/JEWEL/jewel-2.2.0/eventfiles/jewel.Pb-Pb.2760.30-50cent.*.hepmc
#rivet --pwd -a ALICE_2019_JETV2 --ignore-beams /home/william/JEWEL/jewel-2.2.0/eventfiles/jewel.Pb-Pb.2760.0-5cent.01.hepmc
#rivet --pwd -a ALICE_2019_JETV2 --ignore-beams /home/william/JEWEL/jewel-2.2.0/eventfiles/2357154.apollo-acf/*/*.hepmc
rivet-mkhtml --pwd --errs Rivet.yoda
