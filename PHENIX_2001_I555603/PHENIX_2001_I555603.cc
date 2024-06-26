// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/HadronicFinalState.hh"

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

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      //const FinalState fs(Cuts::abseta < 4.9);
      beamOpt = getOption<string>("beam", "NONE");

      const ParticlePair& beam = beams();

      if (beamOpt == "NONE"){
      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu130;
      }

      if (beamOpt == "AUAU130") collSys = AuAu130;

      //declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");


      // Book histograms
      book(_h["avgenergydensity"], 1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      
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
    string beamOpt;
    enum CollisionSystem {AuAu130};
    CollisionSystem collSys;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PHENIX_2001_I555603);

}
