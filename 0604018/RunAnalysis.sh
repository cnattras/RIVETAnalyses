SAMPLE_FILE=../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#This compiles the code
rivet-build RivetSTAR_2006_I715470.so STAR_2006_I715470.cc
#This runs it
rivet --pwd -p calibration.yoda -a STAR_2006_I715470:cent=GEN:beam=AUAU -o Rivet.yoda $SAMPLE_FILE
