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

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

    }


    /// Normalise histograms etc., after the run
    void finalize() {


    }

    //@}



  };


  DECLARE_RIVET_PLUGIN(PHENIX_2016_I1393529);

}
