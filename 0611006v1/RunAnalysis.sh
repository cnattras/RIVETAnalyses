SAMPLE_FILE=../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#This compiles the code
rivet-build RivetPHENIX_2006_I0611006.so PHENIX_2006_I0611006.cc
#This runs it
rivet --pwd -p calibration.yoda -a PHENIX_2006_I0611006:cent=GEN:beam=AUAU -o Rivet.yoda $SAMPLE_FILE
