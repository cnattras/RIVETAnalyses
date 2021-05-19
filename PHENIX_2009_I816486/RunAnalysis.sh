
#This compiles the code
rivet-build RivetPHENIX_2009_I816486.so PHENIX_2009_I816486.cc
#This runs it
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2009_I816486:cent=GEN:beam=AUAU -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#rivet-mkhtml --pwd Rivet.yoda
