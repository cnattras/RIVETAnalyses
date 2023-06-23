SAMPLE_FILE=../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#This compiles the code
rivet-build RivetSTAR_2016_I1442357.so STAR_2016_I1442357.cc
#This runs it
rivet --pwd -p ../Centralities/Calibration/calibration_STAR_AuAu200.yoda -a STAR_2016_I1442357:cent=GEN:beam=AUAU -o Rivet.yoda $SAMPLE_FILE
rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
