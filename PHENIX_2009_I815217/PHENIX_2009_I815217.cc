// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {

  /// @brief Add a short analysis description here
  class PHENIX_2009_I815217 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2009_I815217);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 0.35);

      // Book histograms
      book(_s["Placeholder1"], 1, 1, 1);
      book(_s["Placeholder2"], 2, 1, 1);
      book(_s["Placeholder3"], 2, 1, 2);
      book(_s["Placeholder4"], 2, 1, 3);
      book(_s["Placeholder5"], 3, 1, 1);
      book(_s["Placeholder6"], 4, 1, 1);
      book(_s["Placeholder7"], 4, 1, 2);
      book(_s["Placeholder8"], 4, 1, 3);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Fill histogram with leading b-jet pT
     // _h["XXXX"]->fill(bjets[0].pT()/GeV);

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
    map<string, Scatter2DPtr> _s;
    //@}


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2009_I815217);

}
