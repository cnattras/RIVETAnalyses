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
      cout<<"Hello! system "<<collSys<<endl;

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 5.0);
      declare(fs, "fs");
     declare(ALICE::PrimaryParticles(Cuts::abseta < 1.1 && Cuts::pT > 0.0*MeV && Cuts::abscharge > 0), "APRIM");


      book(_hist_NchCent, 1, 1, 1);
      book(_hist_NchEta, 2, 1, 1);//Note this actually needs to be fleshed out so that it has the eta distributions in all of the centralities

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (collSys == AuAu200){

        double totalEt = 0;
        double deltaeta = 2.2; 


        Particles fsParticles = apply<FinalState>(event,"fs").particles();

        const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
         double c = cent();
      if(fixedcentralityOpt!= "NONE"){
           c = manualCentrality; // Use the float value set in init()
      }
        if (c > 70) vetoEvent;

        for(const Particle& p : fsParticles) // loop over all final state particles
        {
            totalEt += p.Et()/GeV;
        }

        //get charged particles
        Particles chargedParticles = apply<ALICE::PrimaryParticles>(event,"APRIM").particles();
        int nchcounter = chargedParticles.size();
          _hist_NchCent->fill(c,nchcounter/deltaeta);

cout<<"Centrality CATS "<<c<<" dNch/deta "<<nchcounter/deltaeta<<endl;
      
    }
    }


    /// Normalise histograms etc., after the run
    void finalize() {


    }

    /// @}


    /// @name Histograms
    /// @{
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
