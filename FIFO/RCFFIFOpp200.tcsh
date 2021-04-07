#!/usr/local/bin/tcsh
#    1. Path where the output file will be saved
#    () Name of the output file
#    #    2. Number of events
#    #    3. Beam 1
#    #    4. Beam 2
#    #    5. Energy (GeV)
#    #    6. Seed
#    #    7. Minimum pt hard (GeV/c)
#    #    8 Maximum pt hard (GeV/c) 
#  This script takes the same arguments but also compresses the output file

source /cvmfs/sft.cern.ch/lcg/contrib/gcc/8.2.0/x86_64-centos7-gcc8-opt/setup.csh
source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc8-opt/setup.csh
source /gpfs01/star/pwg/elayavalli/RIVET/Rivet-3.1.3/rivetenv.csh
set path = ($path /cvmfs/sft.cern.ch/lcg/external/texlive/2020/bin/x86_64-linux/)


#echo $1

setenv OUTFILE "Rivet$3.$4.$5GeV.$2Events.seed$6.pthard$7to$8.yoda"
setenv FIFOFILE "fifo$6.hepmc"
rm /tmp/$FIFOFILE
mkfifo /tmp/$FIFOFILE
#These analyses need to be in the working directory and the string with analyses needs to include any flags for them
setenv WORKDIR "/phenix/scratch/cen/sampleHepData/LocalAnalyses/"
setenv ANALYSES "BRAHMS_2007_I742956"
setenv CENTCALIBRATION "calibration_PHENIX_AuAu200GeV.yoda"


date
#cd $1
cd $WORKDIR
echo "Current directory: "
pwd
echo Output file: $OUTFILE
echo Working file: $WORKDIR
echo Analyses: $ANALYSES
#Antonio's line
#./main113 /tmp $FIFOFILE 500 Au Au 200 $JOB_NUMBER 0 -1 & rivet --pwd -p calibration_PHENIX_AuAu200GeV.yoda $ANALYSES -o $OUTFILE /tmp/$FIFOFILE
#Old command
#/direct/phenix+u/cen/pythia8303/examples/main117 $1 hepmc$3.$4.$5GeV.$2Events.seed$6.pthard$7to$8.hepmc $2 $3 $4 $5 $6 $7 $8
#New command
#/direct/phenix+u/cen/pythia8244/examples/main117 /tmp fifotest.hepmc 5 p p 200 1 0 -1 > & junk1.txt & rivet --pwd -a BRAHMS_2007_I742956 -o test.yoda /tmp/fifotest.hepmc
setenv COMMAND "/direct/phenix+u/cen/pythia8244/examples/./main117 /tmp $FIFOFILE $2 $3 $4 $5 $6 $7 $8 & rivet --pwd -a $ANALYSES -o $OUTFILE /tmp/$FIFOFILE"
echo $COMMAND
/direct/phenix+u/cen/pythia8244/examples/./main117 /tmp $FIFOFILE $2 $3 $4 $5 $6 $7 $8 & rivet --pwd -a $ANALYSES -o $OUTFILE /tmp/$FIFOFILE
mv $OUTFILE $1/.
date

