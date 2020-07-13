SAMPLE_FILE=../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#This compiles the code
rivet-build RivetSTAR_2016_I1442357.so STAR_2016_I1442357.cc
#This runs it
rivet --pwd -p calibration.yoda -a STAR_2016_I1442357:cent=GEN $SAMPLE_FILE
