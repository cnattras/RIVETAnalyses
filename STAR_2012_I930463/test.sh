
#This compiles the code
rivet-build RivetSTAR_2012_I930463.so STAR_2012_I930463.cc
#This runs it
rivet --pwd -p ../Centralities/Calibration/calibration_STAR_AuAu200GeV.yoda -a STAR_2012_I930463:cent=GEN:beam=AUAU -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#rivet-mkhtml --pwd Rivet.yoda
rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
