#!/bin/bash
SAMPLE_FILE=PYTHIAAuAuFileSMALLTEST.dat
#This compiles the code
rivet-build RivetPHENIX_2010_I857187.so PHENIX_2010_I857187.cc
#This runs it
rivet --pwd -p calibration.yoda -a PHENIX_2010_I857187:cent=GEN $SAMPLE_FILE
