#!/bin/bash
#mkfifo /tmp/fifo.hepmc 
#export RIVET_ANALYSIS_PATH=/tmp
#/usr/local/share/Pythia8/examples/./main144 -c ../RunPYTHIAInDocker/pp200.cmnd -n 100 -o /tmp/fifo.hepmc & \
#rivet --pwd -a RHIC_2019_CentralityCalibration /tmp/fifo.hepmc
#!/bin/bash
#mkfifo /tmp/fifo.hepmc
#export RIVET_ANALYSIS_PATH=/tmp
#/usr/local/share/Pythia8/examples/./main144 -c ../RunPYTHIAInDocker/pp200.cmnd -n 100 -o /tmp/fifo.hepmc &
#pythia8-main144 -c ../RunPYTHIAInDocker/pp200.cmnd -o /tmp/fifo.hepmc
#MAIN_PID=$!
#rivet --pwd -a RHIC_2019_CentralityCalibration /tmp/fifo.hepmc
#wait $MAIN_PID
#rm /tmp/fifo.hepmc

mkfifo pp200GeV.hepmc 
cat pp200GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/pp200GeV.cmnd -c main144HepMCpp.cmnd -c main144Rivet.cmnd -o pp200GeV

pkill -f "cat pp200GeV.hepmc"
rm -f pp200GeV.hepmc
