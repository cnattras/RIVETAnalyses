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
	const FinalState fs(Cuts::abseta < 4.9);
      //fig2
      book(_h["Fig2_Toward"], 1, 1, 1);
      book(_h["Fig2_Away"], 1, 1, 2);
      book(_h["Fig2_Transverse"], 1, 1, 3);
      //fig3
      //book(_h["BBBB"], 3, 1, 1);
      book(_h["Fig3_Transverse_pt02"], 2, 1, 1);
      book(_h["Fig3_Transverse_pt05"], 2, 1, 2);
      //book(_h["Fig3_Transverse"],2,1,3);

      book(_h["Fig4_Toward"], 3, 1, 1);
      book(_h["Fig4_Away"], 3, 1, 2);
      book(_h["Fig4_Transverse"], 3, 1, 3);

      book(_h["Fig5_Transverse_pt02"], 4, 1, 1);
      book(_h["Fig5_Transverse_pt05"], 4, 1, 2);
      //book(_h["Fig5_Transverse"], 4, 1, 3);

      book(_h["Fig6_Transverse_pt02"], 5, 1, 1);
      book(_h["Fig6_Transverse_pt05"], 5, 1, 2);
      //book(_h["Fig6_Transverse"], 5, 1, 3);

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
