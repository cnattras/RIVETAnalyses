SAMPLE_FILE=../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#This compiles the code
rivet-build RivetSTAR_2010_I851937.so STAR_2010_I851937.cc
#This runs it
rivet --pwd -p calibration.yoda -a STAR_2010_I851937:cent=GEN:beam=AUAU -o Rivet.yoda $SAMPLE_FILE
