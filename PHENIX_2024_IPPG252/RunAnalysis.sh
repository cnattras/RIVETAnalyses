#!/bin/bash
rivet-build RivetPHENIX_2024_IPPG252.so PHENIX_2024_IPPG252.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -a PHENIX_2024_IPPG252 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
