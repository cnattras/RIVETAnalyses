#!/bin/bash
rivet-buildplugin RivetPHENIX_2006_I0611006.so PHENIX_2006_I0611006.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p calibration_PHENIX_AuAu39GeV.yoda -a PHENIX_2006_I0611006:cent=GEN -o Rivet_Output.yoda PYTHIAAuAuFileSMALLTEST.dat
