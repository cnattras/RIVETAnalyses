#!/bin/bash
rivet-build RivetSTAR_2016_I1420183.so STAR_2016_I1420183.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_AuAu_200GeV_STAR.yoda -a STAR_2016_I1420183:cent=GEN:beam=pp -o Rivet_pp.yoda PYTHIAAuAuFileSMALLTEST.dat
rivet --pwd -p calibration_AuAu_200GeV_STAR.yoda -a STAR_2016_I1420183:cent=GEN:beam=dAu -o Rivet_dAu.yoda PYTHIAAuAuFileSMALLTEST.dat
rivet-merge -O beam -o Rivet.yoda Rivet_dAu.yoda Rivet_pp.yoda

#rivet --pwd -p calibration_AuAu_200GeV_STAR.yoda -a STAR_2016_I1420183:beam=pp -o Rivet.yoda PYTHIAAuAuFileSMALLTEST.dat