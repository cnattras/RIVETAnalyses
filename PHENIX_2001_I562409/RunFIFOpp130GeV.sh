#!/bin/bash
#mkfifo /tmp/fifo.hepmc 
#export RIVET_ANALYSIS_PATH=/tmp
#/usr/local/share/Pythia8/examples/./main144 -c /work/AuAu200.cmnd -n 100 -o /tmp/fifo.hepmc & \
#rivet --pwd -a RHIC_2019_CentralityCalibration /tmp/fifo.hepmc
#!/bin/bash
#mkfifo /tmp/fifo.hepmc
#export RIVET_ANALYSIS_PATH=/tmp
#/usr/local/share/Pythia8/examples/./main144 -c /work/AuAu200.cmnd -n 100 -o /tmp/fifo.hepmc &
#pythia8-main144 -c /work/AuAu200.cmnd -o /tmp/fifo.hepmc
#MAIN_PID=$!
#rivet --pwd -a RHIC_2019_CentralityCalibration /tmp/fifo.hepmc
#wait $MAIN_PID
#rm /tmp/fifo.hepmc
rivet-build RivetPHENIX_2001_I562409.so PHENIX_2001_I562409.cc


mkfifo pp130GeV.hepmc 
cat pp130GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/pp130GeV.cmnd -c main144HepMCpp.cmnd -c main144Rivetpp130.cmnd -o pp130GeV

pkill -f "cat pp130GeV.hepmc"
rm -f pp130GeV.hepmc
