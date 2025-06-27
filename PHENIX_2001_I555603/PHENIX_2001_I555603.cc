// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2001_I555603 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2001_I555603);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Implement beam options
      beamOpt = getOption<string>("beam", "NONE");

      const ParticlePair& beam = beams();
      int NN = 0.;
      if (beamOpt == "NONE") {
        
            if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
          NN = 197.;
          if (fuzzyEquals(sqrtS()/GeV, 130*NN, 5)) collSys = AuAu130;
      }
      }

      if (beamOpt == "AUAU130") collSys = AuAu130;

      // Declare Centrality
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      // Final State projection
      const FinalState fs(Cuts::abseta < 0.5);
      declare(fs, "fs");

      // Book histograms
      book(_hist_E, 1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (collSys == AuAu130){

        double totalEt = 0;
        double deltaeta = 1; 

        Particles fsParticles = apply<FinalState>(event,"fs").particles();

        const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
        const double c = cent();
        if (c > 50) vetoEvent;

        for(const Particle& p : fsParticles) // loop over all final state particles
        {
            totalEt += p.Et()/GeV;
        }


        if (collSys == AuAu130){
          _hist_E->fill(c,totalEt/deltaeta);
        }
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
    Profile1DPtr _hist_E;
    string beamOpt;
    enum CollisionSystem {NONE, AuAu130};
    CollisionSystem collSys;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PHENIX_2001_I555603);

}
