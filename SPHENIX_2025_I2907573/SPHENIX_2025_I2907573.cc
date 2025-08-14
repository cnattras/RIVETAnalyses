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
  class SPHENIX_2025_I2907573 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SPHENIX_2025_I2907573);


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
      const FinalState fs(Cuts::abseta < etaMax);
      declare(fs, "fs");

      for(int i=0;i<nCB;i++){
          string histoname = "histETCaloCent_"+std::to_string(centBinEdges[i])+"_"+std::to_string(centBinEdges[i+1]);
          string histoname2 = "histETEmcalCent_"+std::to_string(centBinEdges[i])+"_"+std::to_string(centBinEdges[i+1]);
          string histoname3 = "histETHcalCent_"+std::to_string(centBinEdges[i])+"_"+std::to_string(centBinEdges[i+1]);

          std::ostringstream oss;
          oss << std::setw(2) << std::setfill('0') << i+1;
          string rivethistoname = "d05-x01-y"+oss.str();
          std::ostringstream oss2;
          oss2 << std::setw(2) << std::setfill('0') << i+1+nCB;
          string rivethistoname2 = "d05-x01-y"+oss2.str();
          std::ostringstream oss3;
          oss3 << std::setw(2) << std::setfill('0') << i+1+2*nCB;
          string rivethistoname3 = "d05-x01-y"+oss3.str();

          book(_p[histoname], rivethistoname, refData(5, 1, i+1)); 
          book(_p[histoname2], rivethistoname2, refData(5, 1, i+1+nCB)); 
          book(_p[histoname3], rivethistoname3, refData(5, 1, i+1+2*nCB)); 
        }


        // book(_p["histETEmcalCent"], "d06-x01-y01", refData(6, 1, 1));
        // book(_p["histETHcalCent"], "d06-x01-y02", refData(6, 1, 2));
        book(_p["histETCentDependence"], "d07-x01-y02", refData(7, 1, 2));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        Particles fsParticles = apply<FinalState>(event,"fs").particles();

        const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
         double c = cent();
      if(fixedcentralityOpt!= "NONE"){
           c = manualCentrality; // Use the float value set in init()
      }
        if (c > 70) vetoEvent;


    // double totalEt = 0;
    double totalEtMidRap = 0;
    double eTEta[7] = {0,0,0,0,0, 0,0};
    for(const Particle& p : fsParticles)
        {
            // if(abs(p.eta())<etaMax)totalEt += p.Et()/GeV;
            if(abs(p.eta())<etaMax) totalEtMidRap += p.Et()/GeV;
            int myetabin = findEtaBin(p.eta());
            if(myetabin!=-1) eTEta[myetabin]+= p.Et()/GeV;
        }

        int mycentbin = findCentBin(c);
       string histoname = "histETCaloCent_"+std::to_string(centBinEdges[mycentbin])+"_"+std::to_string(centBinEdges[mycentbin+1]);
       string histoname2 = "histETEmcalCent_"+std::to_string(centBinEdges[mycentbin])+"_"+std::to_string(centBinEdges[mycentbin+1]);
       string histoname3 = "histETHcalCent_"+std::to_string(centBinEdges[mycentbin])+"_"+std::to_string(centBinEdges[mycentbin+1]);
       //cout<<"CB :"<<mycentbin<<" centrality "<<c<<endl;
        if(mycentbin!=-1){
            for(int i=0;i<nEB;i++){
                _p[histoname]->fill((etaBinEdges[i]+etaBinEdges[i+1])/2,((double) eTEta[i])/(etaBinEdges[i+1]-etaBinEdges[i]));
                _p[histoname2]->fill((etaBinEdges[i]+etaBinEdges[i+1])/2,((double) eTEta[i])/(etaBinEdges[i+1]-etaBinEdges[i]));
                _p[histoname3]->fill((etaBinEdges[i]+etaBinEdges[i+1])/2,((double) eTEta[i])/(etaBinEdges[i+1]-etaBinEdges[i]));
          //cout<<"Filling eta cent "<<mycentbin<<" "<<(etaBinEdges[i]+etaBinEdges[i+1])/2<<" dNch/deta "<<((double) eTEta[i])/(etaBinEdges[i+1]-etaBinEdges[i])<<endl;
            } 
            //cout<<"totalEtMidRap "<<totalEtMidRap<<" total ET "<<totalEt<<endl;
        }
        _p["histETCentDependence"]->fill(c,totalEtMidRap*2.0/npart[mycentbin]/2.0/etaMax);
        cout<<" cb "<<mycentbin<<" c "<<c<<" ET/npart "<<totalEtMidRap*2.0/npart[mycentbin]/2.0/etaMax<<endl;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // _p["histETHcalCent"]->scaleW(2.0/etaMax);
      // _p["histETEmcalCent"]->scaleW(2.0/etaMax);
      // int i=0;
      // for (auto& bin : _p["histETCentDependence"]->bins()) {
      //   bin.scale(2.0/npart[i]);
      //     i++;
      // }
      
    }

    /// @}
    double const etaMax = 1.15;

    const double etaBinEdges[7] =  {-1.101400e+00, -7.351000e-01, -3.689000e-01, -2.700000e-03, 3.639000e-01, 7.309000e-01, 1.098000e+00};
    const int nEB = 6;
    const int centBinEdges[9] = {0,5,10,20,30,40,50,60,70};
    const int nCB = 8;
    const double npart[9] = {349.96, 301.76, 238.14, 170.83, 118.55, 78.26, 47.8, 26.15};

    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    string beamOpt;
    string fixedcentralityOpt;
    double manualCentrality = -1.0;
    enum CollisionSystem {NONE, AuAu200};
    CollisionSystem collSys;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SPHENIX_2025_I2907573);

}
