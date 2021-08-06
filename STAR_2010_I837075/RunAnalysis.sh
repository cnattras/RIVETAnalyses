#!/bin/bash
#if test -f Rivet.yoda; then
#    rm Rivet.yoda
#fi
rivet-build RivetSTAR_2010_I837075.so STAR_2010_I837075.cc
#make
 # docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet-build RivetSTAR_2010_I837075.so STAR_2010_I837075.cc
# export RIVET_ANALYSIS_PATH=$PWD
# docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet  --pwd -p calibration_AuAu_200GeV_STAR.yoda -a STAR_2010_I837075:cent=GEN:beam=AuAu200 -o Rivet.yoda PYTHIAAuAuFileSMALLTEST.dat
#docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet  --pwd -p calibration_AuAu_200GeV_STAR.yoda -a STAR_2010_I837075:cent=GEN:beam=AuAu200 -o Rivet_CuCu.yoda hepmcCu.Cu.200GeV.1000Events.seed7.pthard0to-1.hepmc -n 300
rivet  --pwd -p ../Centralities/Calibration/calibration_STAR_PHENIX_CuCu200GeV.yoda -a STAR_2010_I837075:cent=GEN:beam=CUCU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
# docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet  --pwd -p calibration_AuAu_200GeV_STAR.yoda -a STAR_2010_I837075:cent=GEN:beam=AuAu200 -o Rivet_pp.yoda hepmcp.p.200GeV.10000Events.seed100.pthard0to-1.hepmc

