#!/bin/bash
rivet-build RivetPHENIX_2011_I886590.so PHENIX_2011_I886590.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -a PHENIX_2011_I886590:beam=pp200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat

#../../sims/pp_200GeV/hepmc_pp_200GeV_28.hepmc
