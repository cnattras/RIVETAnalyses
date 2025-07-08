#!/bin/bash
#mkfifo /tmp/fifo.hepmc 
#export RIVET_ANALYSIS_PATH=/tmp
#/usr/local/share/Pythia8/examples/./main144 -c ../RunPYTHIAInDocker/CuAu200.cmnd -n 100 -o /tmp/fifo.hepmc & \
#rivet --pwd -a RHIC_2019_CentralityCalibration /tmp/fifo.hepmc
#!/bin/bash
#mkfifo /tmp/fifo.hepmc
#export RIVET_ANALYSIS_PATH=/tmp
#/usr/local/share/Pythia8/examples/./main144 -c ../RunPYTHIAInDocker/CuAu200.cmnd -n 100 -o /tmp/fifo.hepmc &
#pythia8-main144 -c ../RunPYTHIAInDocker/CuAu200.cmnd -o /tmp/fifo.hepmc
#MAIN_PID=$!
#rivet --pwd -a RHIC_2019_CentralityCalibration /tmp/fifo.hepmc
#wait $MAIN_PID
#rm /tmp/fifo.hepmc

mkfifo CuAu200GeV.hepmc 
cat CuAu200GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/CuAu200GeV.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o CuAu200GeV

pkill -f "cat CuAu200GeV.hepmc"
rm -f CuAu200GeV.hepmc
