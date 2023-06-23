
#This compiles the code
rivet-build RivetSTAR_2006_I715470.so STAR_2006_I715470.cc
#This runs it

rivet --pwd -p ../Centralities/Calibration/calibration_STAR_AuAu200GeV.yoda -a STAR_2006_I715470:cent=GEN:beam=AUAU200 -o Rivet.yoda ../testfiles/PYTHIAAuAuFileSMALLTEST.dat

rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda
