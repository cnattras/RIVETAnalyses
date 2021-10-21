#!/bin/bash
# ATTENTION! FIFO requires a lot of virtual memory. When entering in debug more, use: qsub -I -l nodes=1:ppn=4 -q debug
# This will request 4 cores and sufficient virtual memory
# Copy the main113 file to the same directory where you have your Rivet Analysis
# On ACF you can find runPythia here: /lustre/haven/proj/UTK0019/Rivet
# Change the name "PHENIX_2008_I778168" to your Rivet Analysis
# Have your centrality calibration file in the same folder or point $CALIBRATION to the right path/file
ANALYSIS="ALICE_2015_I1328629"
ANALYSIS_DIR=$PWD
source /lustre/haven/proj/UTK0019/Rivet/rivet.sh
rm RivetALICE_2015_I1328629.so
rivet-build RivetALICE_2015_I1328629.so ALICE_2015_I1328629.cc
FIFOFILE="fifoRia.hepmc"
# Number of events
NEVENTS="100000"
# Beams
BEAM1="p"
BEAM2="p"
# Energy in the center of mass system
CMS_ENERGY="7000"
# Pythia seed. If 0 pythia will generate a random number
GENERATOR_SEED="0"
# Min and max pT-hard. if PTHARDMIN > PTHARDMAX then pythia switches to minimum bias (No pt-hard)
PTHARDMIN="70"
PTHARDMAX="80"
# Centrality Calibration file
CALIBRATION="/lustre/haven/proj/UTK0019/rpatel77/RIVETAnalyses/Centralities/Calibration/calibration_PHENIX_AuAu200GeV.yoda"
# Flags of your analysis (Ex. centrality: cent=GEN)
RIVET_FLAGS=""
rm /tmp/$FIFOFILE
mkfifo /tmp/$FIFOFILE
./runPythia /tmp $FIFOFILE $NEVENTS $BEAM1 $BEAM2 $CMS_ENERGY $GENERATOR_SEED $PTHARDMIN $PTHARDMAX & rivet --pwd -a $ANALYSIS$RIVET_FLAGS -o Rivet.yoda /tmp/$FIFOFILE
