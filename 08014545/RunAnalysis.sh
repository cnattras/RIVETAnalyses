SAMPLE_FILE=../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#This compiles the code
rivet-build RivetPHENIX_2008_I778396.so PHENIX_2008_I778396.cc
#This runs it
rivet --pwd -p calibration.yoda -a PHENIX_2008_I778396:cent=GEN:beam=AUAU -o Rivet.yoda $SAMPLE_FILE
