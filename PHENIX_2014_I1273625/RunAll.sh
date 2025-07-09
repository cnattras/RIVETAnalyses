# cp main144RivetTEMPLATE.cmnd main144RivetAuAu7.cmnd
# sed -i 's/TEMPLATE/AuAu7/' main144RivetAuAu7.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetAuAu14.cmnd
# sed -i 's/TEMPLATE/AuAu14/' main144RivetAuAu14.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetAuAu19.cmnd
# sed -i 's/TEMPLATE/AuAu19/' main144RivetAuAu19.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetAuAu27.cmnd
# sed -i 's/TEMPLATE/AuAu27/' main144RivetAuAu27.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetAuAu39.cmnd
# sed -i 's/TEMPLATE/AuAu39/' main144RivetAuAu39.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetAuAu62.cmnd
# sed -i 's/TEMPLATE/AuAu62/' main144RivetAuAu62.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetAuAu130.cmnd
# sed -i 's/TEMPLATE/AuAu130/' main144RivetAuAu130.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetAuAu200.cmnd
# sed -i 's/TEMPLATE/AuAu200/' main144RivetAuAu200.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetCuCu62.cmnd
# sed -i 's/TEMPLATE/CuCu62/' main144RivetCuCu62.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetCuCu200.cmnd
# sed -i 's/TEMPLATE/CuCu200/' main144RivetCuCu200.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetCuAu200.cmnd
# sed -i 's/TEMPLATE/CuAu200/' main144RivetCuAu200.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetUU193.cmnd
# sed -i 's/TEMPLATE/UU193/' main144RivetUU193.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetdAu200.cmnd
# sed -i 's/TEMPLATE/dAu200/' main144RivetdAu200.cmnd
# cp main144RivetTEMPLATE.cmnd main144RivetHeAu200.cmnd
# sed -i 's/TEMPLATE/HeAu200/' main144RivetHeAu200.cmnd
# cp main144RivetTEMPLATE.cmnd main144Rivetpp200.cmnd
# sed -i 's/TEMPLATE/pp200/' main144Rivetpp200.cmnd
#!/bin/bash

rm RivetPHENIX_2014_I1273625.so

rivet-build RivetPHENIX_2014_I1273625.so PHENIX_2014_I1273625.cc
export RIVET_ANALYSIS_PATH=$PWD


mkfifo AuAu62GeV.hepmc
cat AuAu62GeV.hepmc > /dev/null &
pythia8-main144 -c ../RunPYTHIAInDocker/AuAu62GeV.cmnd -c main144HepMC.cmnd -c main144RivetAuAu62.cmnd -o AuAu62GeV
pkill -f "cat AuAu62GeV.hepmc"
rm -f AuAu62GeV.hepmc
mkfifo AuAu130GeV.hepmc
cat AuAu130GeV.hepmc > /dev/null &
pythia8-main144 -c ../RunPYTHIAInDocker/AuAu130GeV.cmnd -c main144HepMC.cmnd -c main144RivetAuAu130.cmnd -o AuAu130GeV
pkill -f "cat AuAu130GeV.hepmc"
rm -f AuAu130GeV.hepmc
mkfifo AuAu200GeV.hepmc
cat AuAu200GeV.hepmc > /dev/null &
pythia8-main144 -c ../RunPYTHIAInDocker/AuAu200GeV.cmnd -c main144HepMC.cmnd -c main144RivetAuAu200.cmnd -o AuAu200GeV
pkill -f "cat AuAu200GeV.hepmc"
rm -f AuAu200GeV.hepmc
