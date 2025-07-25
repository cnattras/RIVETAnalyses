#!/bin/bash
rivet-build RivetPHENIX_2025_I2820229.so PHENIX_2025_I2820229.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -a PHENIX_2025_I2820229 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
