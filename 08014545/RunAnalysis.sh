SAMPLE_FILE=../testfiles/PYTHIAAuAuFileSMALLTEST.dat
NEW_FILE=/lustre/haven/proj/UTK0019/Pythia8/pythia8243/examples/hepmc_AuAu_MB_1.hepmc
#This compiles the code
rivet-build RivetPHENIX_2008_I778396.so PHENIX_2008_I778396.cc
#This runs it
rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a PHENIX_2008_I778396:cent=GEN:beam=AUAU -o Rivet.yoda $NEW_FILE
#rivet-mkhtml --pwd Rivet.yoda
