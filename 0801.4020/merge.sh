#!/bin/bash
rm *.so
export RIVET_ANALYSIS_PATH=$PWD
rivet-buildplugin RivetPHENIX_2008_I778168.so PHENIX_2008_I778168.cc
rivet-merge -o Rivet_merged.yoda -O beam Rivet.yoda Rivet_pp.yoda
