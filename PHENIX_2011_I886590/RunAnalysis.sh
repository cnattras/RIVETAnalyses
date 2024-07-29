#!/bin/bash
rivet-build RivetPHENIX_2011_I886590.so PHENIX_2011_I886590.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -a PHENIX_2011_I886590:beam=AUAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#../../sims/pp_200GeV/hepmc_pp_200GeV_*.hepmc
#rivet --pwd -a PHENIX_2011_I886590:beam=pp62 -o Rivet_pp62.yoda ../../sims/pp_200GeV/hepmc_pp_200GeV_*.hepmc
#../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#../../sims/pp_200GeV/hepmc_pp_200GeV_28.hepmc

#rivet-mkhtml --pwd Rivet_pp200.yoda
#cp -r rivet-plots/ rivet-plots-pp200/

#rivet-mkhtml --pwd Rivet_pp62.yoda
#cp -r rivet-plots/ rivet-plots-pp62/


