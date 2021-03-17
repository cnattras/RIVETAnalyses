#!/bin/bash
ANALYSIS="PHENIX_2008_I778168"
ANALYSIS_DIR=$PWD
source /lustre/haven/proj/UTK0019/Rivet3.1.3/local/rivetenv.sh
rm RivetPHENIX_2008_I778168.so
rivet-build RivetPHENIX_2008_I778168.so PHENIX_2008_I778168.cc
FIFOFILE="fifo$ANALYSIS.hepmc"
#Number of events
NEVENTS="100"
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
#Flags of your analysis (Ex. centrality: cent=GEN)
RIVET_FLAGS=":cent=GEN"
rm /tmp/$FIFOFILE
mkfifo /tmp/$FIFOFILE
./main113 /tmp $FIFOFILE $NEVENTS $BEAM1 $BEAM2 $CMS_ENERGY $GENERATOR_SEED $PTHARDMIN $PTHARDMAX & rivet --pwd -p calibration_PHENIX_AuAu62GeV.yoda -a $ANALYSIS$RIVET_FLAGS -o Rivet.yoda /tmp/$FIFOFILE
