#!/bin/bash
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetPHENIX_2008_I777211.so PHENIX_2008_I777211.cc
rivet-merge -o Angantyr.yoda -O beam Angantyr_AuAu_200GeV.yoda Angantyr_pp_200GeV.yoda
 
