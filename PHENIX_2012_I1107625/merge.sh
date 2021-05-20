#!/bin/bash
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetPHENIX_2012_I1107625.so PHENIX_2012_I1107625.cc
rivet-merge -o Angantyr.yoda -O beam Angantyr_pp_39GeV.yoda Angantyr_pp_62GeV.yoda Angantyr_AuAu_39GeV.yoda Angantyr_AuAu_62GeV.yoda  
