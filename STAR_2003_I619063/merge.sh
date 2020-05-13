#!/bin/bash
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetSTAR_2003_I619063.so STAR_2003_I619063.cc
rivet-merge -o Angantyr.yoda -O beam Angantyr_AuAu_200GeV.yoda Angantyr_pp_200GeV.yoda Angantyr_AuAu_130GeV.yoda
 
