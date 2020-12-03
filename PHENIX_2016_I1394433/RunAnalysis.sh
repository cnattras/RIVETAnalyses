#!/bin/bash
rm RivetPHENIX_2016_I1394433.so

docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet-build RivetPHENIX_2016_I1394433.so PHENIX_2016_I1394433.cc

export RIVET_ANALYSIS_PATH=$PWD

docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet --pwd -p Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2016_I1394433:cent=GEN -o Rivet_AUAU_200.yoda testfiles/PYTHIAAuAuFileSMALLTEST.dat

docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet --pwd -p Centralities/Calibration/calibration_PHENIX_AuAu62GeV.yoda -a PHENIX_2016_I1394433:cent=GEN -o Rivet_AUAU_62.yoda testfiles/AuAu_62.hepmc

docker run -i  --rm  -u `id -u $USER`:`id -g`  -v $PWD:$PWD -w $PWD  hepstore/rivet:3.1.2 rivet --pwd -p Centralities/Calibration/calibration_PHENIX_AuAu39GeV.yoda -a PHENIX_2016_I1394433:cent=GEN -o Rivet_AUAU_39.yoda testfiles/AuAu_39.hepmc







