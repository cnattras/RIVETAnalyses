#!/bin/bash
rivet-build RivetBRAHMS_2007_I742956.so BRAHMS_2007_I742956.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -a BRAHMS_2007_I742956 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
