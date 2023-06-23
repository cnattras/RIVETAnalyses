
#This compiles the code
rivet-build RivetPHENIX_2011_I872172.so PHENIX_2011_I872172.cc
#This runs it
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2011_I872172:cent=GEN:beam=AUAU -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat
#rivet-mkhtml --pwd Rivet.yoda
rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
