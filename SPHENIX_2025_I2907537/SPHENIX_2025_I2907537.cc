// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class SPHENIX_2025_I2907537 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SPHENIX_2025_I2907537);


    /// @name Analysis methods
    /// @{

//Class for finding the centrality bin 
    int findCentBin(float c) {
    // Number of bins is (number of edges - 1)

    // Handle out-of-range values explicitly (optional)
    if (c <= centBinEdges[0]) return -1;      // below range (no bin)
    if (c > centBinEdges[nCB]) return -1;  // above range (no bin)

    // Binary search for efficiency
    int low = 0, high = nCB - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (c > centBinEdges[mid] && c <= centBinEdges[mid + 1]) {
            return mid;
        } else if (c <= centBinEdges[mid]) {
            high = mid - 1;
        } else { // c > centBinEdges[mid+1]
            low = mid + 1;
        }
    }

    // If no bin found (should not happen), return -1
    return -1;
}

int findEtaBin(double c) {
    //cout<<"starting search "<<etaBinEdges[0]<<" - "<<etaBinEdges[nEB]<<endl;
    // Out of range checks:
    if (c <= etaBinEdges[0]) return -1;        // < lowest edge — no bin
    if (c > etaBinEdges[nEB]) return -1;    // > highest edge — no bin

    // Binary search to find bin
    int low = 0;
    int high = nEB; // bin indices [0..10]
    //cout<<"low "<<low<<" high "<<high<<endl;

    while (low <= high) {
        int mid = (low + high) / 2;
        if (c > etaBinEdges[mid] && c <= etaBinEdges[mid + 1]) {
    //cout<<"low "<<low<< " "<<etaBinEdges[low]<<" high "<<high<< " "<<etaBinEdges[high]<<" returning "<<mid<<endl; 
            return mid; // found bin
        } else if (c <= etaBinEdges[mid]) {
    //cout<<"low "<<low<< " "<<etaBinEdges[low]<<" high "<<high<< " "<<etaBinEdges[high]<<endl; 
            high = mid - 1; // search left
        } else {
    //cout<<"low "<<low<< " "<<etaBinEdges[low]<<" high "<<high<< " "<<etaBinEdges[high]<<endl; 
            low = mid + 1;  // search right
        }
    }
    //cout<<"low "<<low<< " "<<etaBinEdges[low]<<" high "<<high<< " "<<etaBinEdges[high]<<" FAIL "<<endl;

    return -1; // should never happen if bounds are respected
}

    /// Book histograms and initialise projections before the run
    void init() {


      // Implement beam options
      beamOpt = getOption<string>("beam", "NONE");
      fixedcentralityOpt = getOption<string>("fixedcentrality", "NONE");

      const ParticlePair& beam = beams();
      int NN = 0.;
      if (beamOpt == "NONE") {
        
            if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
          NN = 197.;
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1e-3)){ collSys = AuAu200;        }
      }
      }
    if (beamOpt == "AUAU200") collSys = AuAu200;
      // Declare Centrality
      declareCentrality(RHICCentrality("sPHENIX"), "RHIC_2019_CentralityCalibration:exp=sPHENIX", "CMULT", "CMULT");
      if(fixedcentralityOpt!= "NONE"){
        manualCentrality = std::stof(fixedcentralityOpt);
      }

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 5.0);
      declare(fs, "fs");
     declare(ALICE::PrimaryParticles(Cuts::abseta < 1.1 && Cuts::pT > 0.0*MeV && Cuts::abscharge > 0), "APRIM");


        book(_hist_NchCent, "d01-x01-y01", refData(1, 1, 1));

        for(int i=0;i<nCB;i++){
          string histoname = "_hist_NchEta_Cent_"+std::to_string(centBinEdges[i])+"_"+std::to_string(centBinEdges[i+1]);

          std::ostringstream oss;
          oss << std::setw(2) << std::setfill('0') << i+1;
          string rivethistoname = "d02-x01-y"+oss.str();

          book(_p[histoname], rivethistoname, refData(2, 1, i+1)); 
          //cout<<rivethistoname<<" "<<histoname<<endl;
        }
        //book(_hist_NchEta, "d02-x01-y01", refData(2, 1, 1));
      // book(_TEST, "NchCent", 10, 0.0, 100.0); 
      // book(_hist_NchCent, 1, 1, 1);
      // book(_hist_NchEta, 2, 1, 1);//Note this actually needs to be fleshed out so that it has the eta distributions in all of the centralities

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (collSys == AuAu200){

        double deltaeta = 0.6;//2.2; 
        int nchcounter = 0;//chargedParticles.size();
        int nchcounters[11] = {0,0,0,0,0, 0,0,0,0,0, 0};


        Particles fsParticles = apply<FinalState>(event,"fs").particles();
        Particles chargedParticles = apply<ALICE::PrimaryParticles>(event,"APRIM").particles();

        const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
         double c = cent();
      if(fixedcentralityOpt!= "NONE"){
           c = manualCentrality; // Use the float value set in init()
      }
        if (c > 70) vetoEvent;


        for(const Particle& p : chargedParticles) // loop over all final state particles
        {   
            int myetabin = findEtaBin(p.eta());
            if(myetabin!=-1) nchcounters[myetabin]++;
            //cout<<"eta bin "<<myetabin<<" eta "<<p.eta()<<endl;
            if(abs(p.eta())<0.3) nchcounter++;
        }
        int mycentbin = findCentBin(c);
       string histoname = "_hist_NchEta_Cent_"+std::to_string(centBinEdges[mycentbin])+"_"+std::to_string(centBinEdges[mycentbin+1]);
       cout<<"CB :"<<mycentbin<<" centrality "<<c<<endl;
        for(int i=0;i<nEB;i++){
            if(mycentbin!=-1) _p[histoname]->fill((etaBinEdges[i]+etaBinEdges[i+1])/2,((double) nchcounters[i])/(etaBinEdges[i+1]-etaBinEdges[i]));
          //cout<<"Filling eta "<<(etaBinEdges[i]+etaBinEdges[i+1])/2<<" dNch/deta "<<((double) nchcounters[i])/(etaBinEdges[i+1]-etaBinEdges[i])<<endl;
        }

        //get charged particles;
          _hist_NchCent->fill(c,((double)nchcounter)/deltaeta);

//cout<<"Centrality MEOW "<<c<<" dNch/deta "<<nchcounter/deltaeta<<endl;
      
    }
    }


    /// Normalise histograms etc., after the run
    void finalize() {


    }

    /// @}


    /// @name Histograms
    /// @{
    const double etaBinEdges[12] =  {-1.100000e+00, -9.000000e-01, -7.000000e-01, -5.000000e-01, -3.000000e-01, -1.000000e-01, 1.000000e-01, 3.000000e-01, 5.000000e-01, 7.000000e-01, 9.000000e-01, 1.100000e+00};
    const int nEB = 11;
    const int centBinEdges[16] = {0,3,6,10,15,  20,25,30,35,40,  45,50,55,60,65,  70};
    const int nCB = 15;
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    Profile1DPtr _hist_NchCent,_hist_NchEta;
    string beamOpt;
    string fixedcentralityOpt;
    double manualCentrality = -1.0;
    enum CollisionSystem {NONE, AuAu200};
    CollisionSystem collSys;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SPHENIX_2025_I2907537);

}
