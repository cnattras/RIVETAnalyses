// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BRAHMS_2007_I742956 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BRAHMS_2007_I742956);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);

      // Book histograms
      // specify custom binning
      book(_h["XXXX"], "myh1", 20, 0.0, 100.0);	
      book(_h["CrsSecPI+"], 1, 1, 1);
      book(_h["CrsSecPI-"], 2, 1, 1);
      book(_h["CrsSecK+"], 3, 1, 1);
      book(_h["CrsSecP"], 4, 1, 1);
      book(_h["CrsSecAntP"], 5, 1, 1);
      book(_h["CrsSecK+"], 6, 1, 1);
//Rapidity Changes at this point
      book(_h["CrsSecPI+R2"], 7, 1, 1);
      book(_h["CrsSecPI-R2"], 8, 1, 1);
      book(_h["CrsSecK+R2"], 9, 1, 1);
      book(_h["CrsSecPR2"], 10, 1, 1);
      book(_h["CrsSecAntPR2"], 11, 1, 1);
      book(_h["CrsSecK+R2"], 12, 1, 1);
 }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // Fill histogram with leading b-jet pT
      //_h["XXXX"]->fill(bjets[0].pT()/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h["XXXX"]); // normalize to unity

    }

    //@}


    //@name Histograms
    // @{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}


  };


  DECLARE_RIVET_PLUGIN(BRAHMS_2007_I742956);

}
