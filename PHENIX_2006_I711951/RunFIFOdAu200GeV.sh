#!/bin/bash
#mkfifo /tmp/fifo.hepmc 
#export RIVET_ANALYSIS_PATH=/tmp
#/usr/local/share/Pythia8/examples/./main144 -c /work/dAu200.cmnd -n 100 -o /tmp/fifo.hepmc & \
#rivet --pwd -a RHIC_2019_CentralityCalibration /tmp/fifo.hepmc
#!/bin/bash
#mkfifo /tmp/fifo.hepmc
#export RIVET_ANALYSIS_PATH=/tmp
#/usr/local/share/Pythia8/examples/./main144 -c /work/dAu200.cmnd -n 100 -o /tmp/fifo.hepmc &
#pythia8-main144 -c /work/dAu200.cmnd -o /tmp/fifo.hepmc
#MAIN_PID=$!
#rivet --pwd -a RHIC_2019_CentralityCalibration /tmp/fifo.hepmc
#wait $MAIN_PID
#rm /tmp/fifo.hepmc
rivet-build RivetPHENIX_2006_I711951.so PHENIX_2006_I711951.cc


mkfifo dAu200GeV.hepmc 
cat dAu200GeV.hepmc > /dev/null &

pythia8-main144 -c ../RunPYTHIAInDocker/dAu200GeV.cmnd -c main144HepMC.cmnd -c main144Rivet.cmnd -o dAu200GeV

pkill -f "cat dAu200GeV.hepmc"
rm -f dAu200GeV.hepmc
