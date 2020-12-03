// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2016_I1393529 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2016_I1393529);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      const FinalState fs(Cuts::abseta < 0.35 && Cuts::pT > 0.0*GeV && Cuts::pT < 20.0*GeV);
      declare(fs, "fs");

      book(_h["InvYield_charm"], 1, 1, 1);
      book(_h["InvYield_bottom"], 1, 1, 2);
      book(_h["bfrac"], 2, 1, 1);
      book(_h["RAA_c2e"], 3, 1, 1);
      book(_h["RAA_b2e"], 3, 1, 2);
      book(_h["RAA_ratio"], 4, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

    Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();

      for(const Particle& p : fsParticles) 
      {
         if(p.fromCharm()) _h["InvYield_charm"]->fill(p.pT()/GeV);
         if(p.fromBottom()) _h["InvYield_bottom"]->fill(p.pT()/GeV);
      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {


    }

    /// @name Histograms
    //@}
    map<string, Histo1DPtr> _h; 
    map<string, Profile1DPtr> _p; 
    map<string, CounterPtr> _c; 
    //@}

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2016_I1393529);

}
