#!/bin/bash
#ATTENTION! FIFO requires a lot of virtual memory. When entering in debug more, use: qsub -I -l nodes=1:ppn=4 -q debug
#This will request 4 cores and sufficient virtual memory
#Copy the main113 file to the same directory where you have your Rivet Analysis
#On ACF you can find main113 here: /lustre/haven/proj/UTK0019/FIFO
#Change the name "PHENIX_2008_I778168" to your Rivet Analysis
#Have your centrality calibration file in the same folder or point $CALIBRATION to the right path/file
ANALYSIS="PHENIX_2004_I624474"
ANALYSIS_DIR=$PWD
source /lustre/haven/proj/UTK0019/Rivet3.1.3/local/rivetenv.sh
rm RivetPHENIX_2004_I624474.so
rivet-build RivetPHENIX_2004_I624474.so PHENIX_2004_I624474.cc
FIFOFILE="fifo$ANALYSIS.hepmc"
#Number of events
NEVENTS="1000"
#Beams
BEAM1="Au"
BEAM2="Au"
#Energy in the center of mass system
CMS_ENERGY="200"
#Pythia seed. If 0 pythia will generate a random number
GENERATOR_SEED="0"
#Min and max pT-hard. if PTHARDMIN > PTHARDMAX then pythia switches to minimum bias (No pt-hard)
PTHARDMIN="0"
PTHARDMAX="-1"
#Centrality Calibration file
CALIBRATION="calibration_PHENIX_AuAu200GeV.yoda"
#Flags of your analysis (Ex. centrality: cent=GEN)
RIVET_FLAGS=":cent=GEN"
rm /tmp/$FIFOFILE
mkfifo /tmp/$FIFOFILE
./main113 /tmp $FIFOFILE $NEVENTS $BEAM1 $BEAM2 $CMS_ENERGY $GENERATOR_SEED $PTHARDMIN $PTHARDMAX & rivet --pwd -p ../Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda -a $ANALYSIS$RIVET_FLAGS -o Rivet.yoda /tmp/$FIFOFILE
