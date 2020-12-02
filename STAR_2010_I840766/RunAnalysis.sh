#!/bin/bash
rm RivetSTAR_2010_I840766.so
docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet-build RivetSTAR_2010_I840766.so STAR_2010_I840766.cc
export RIVET_ANALYSIS_PATH=$PWD
docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet --pwd -p Centralities/Calibration/calibration_AuAu_200GeV_STAR.yoda -a STAR_2010_I840766:cent=GEN:beam=pp -o Rivet_AuAu.yoda hepmc_pp_200GeV_1.hepmc
docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet --pwd -p Centralities/Calibration/calibration_AuAu_200GeV_STAR.yoda -a STAR_2010_I840766:cent=GEN:beam=AuAu -o Rivet_pp.yoda PYTHIAAuAuFileSMALLTEST.dat 
docker run -i --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet-merge -O beam -o Rivet_final.yoda Rivet_AuAu.yoda Rivet_pp.yoda
