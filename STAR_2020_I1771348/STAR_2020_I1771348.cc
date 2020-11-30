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
  class STAR_2020_I1771348 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2020_I1771348);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

	declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve dressed leptons, sorted by pT

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut

      // Remove all jets within dR < 0.2 of a dressed lepton

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection

      // Veto event if there are no b-jets

      // Apply a missing-momentum cut

      // Fill histogram with leading b-jet pT

    }


    /// Normalise histograms etc., after the run
    void finalize() {


    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}


  };


  DECLARE_RIVET_PLUGIN(STAR_2020_I1771348);

}
