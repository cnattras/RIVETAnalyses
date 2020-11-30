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
  class STAR_2010_I837075 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2010_I837075);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
        declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /* // Retrieve dressed leptons, sorted by pT */
      /* vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons(); */

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

      /* // Remove all jets within dR < 0.2 of a dressed lepton */
      /* idiscardIfAnyDeltaRLess(jets, leptons, 0.2); */

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
      Jets bjets = filter_select(jets, [](const Jet& jet) {
        return  jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
      });

      // Veto event if there are no b-jets
      if (bjets.empty())  vetoEvent;

      // Apply a missing-momentum cut
      if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent;

      // Fill histogram with leading b-jet pT
      _h["XXXX"]->fill(bjets[0].pT()/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h["XXXX"]); // normalize to unity
      normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
      scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}


  };


  DECLARE_RIVET_PLUGIN(STAR_2010_I837075);

}
