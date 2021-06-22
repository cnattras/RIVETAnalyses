#!/bin/bash
rivet-build RivetSTAR_2009_I793126.so STAR_2009_I793126.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -p ../Centralities/Calibration/calibration_AuAu_200GeV_STAR.yoda -a STAR_2009_I793126:cent=GEN:beam=AUAU130 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
#rivet --pwd -p ../Centralities/Calibration/calibration_AuAu_200GeV_STAR.yoda -a STAR_2009_I793126:cent=GEN:beam=AUAU62 -o Rivet_AUAU62.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#rivet --pwd -p ../Centralities/Calibration/calibration_AuAu_200GeV_STAR.yoda -a STAR_2009_I793126:cent=GEN:beam=DAU200 -o Rivet_DAU200.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat


#rivet-mkhtml --pwd Rivet_AUAU130.yoda
#cp -r rivet-plots/ rivet-plots-AUAU130/
#rm -rf rivet-plots-AUAU130/rivet-plots/

#rivet-mkhtml --pwd Rivet_AUAU62.yoda
#cp -r rivet-plots/ rivet-plots-AUAU62/
#rm -rf rivet-plots-AUAU62/rivet-plots/

#rivet-mkhtml --pwd Rivet_DAU200.yoda
#cp -r rivet-plots/ rivet-plots-DAU200/
#rm -rf rivet-plots-DAU200/rivet-plots/

#rm -rf rivet-plots/
