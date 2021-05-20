
#This compiles the code
rivet-build RivetALICE_2016_I1427723.so ALICE_2016_I1427723.cc
#This runs it
rivet --pwd -p ../Centralities/Calibration/calibration_ALICE.yoda -a ALICE_2016_I1427723:cent=GEN:beam=AUAU -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#rivet-mkhtml --pwd Rivet.yoda
