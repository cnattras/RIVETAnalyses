// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2020_I1755387 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2020_I1755387);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
	const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::abscharge > 0);
	declare(aprim, "aprim");

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");


      // Book histograms
      // specify custom binning
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["hist1"], 1, 1, 1);
      book(_h["hist2"], 2, 1, 1);
      book(_h["hist3"], 3, 1, 1);
      book(_h["hist4"], 4, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
        const ALICE::PrimaryParticles aprim = apply<ALICE::PrimaryParticles>(event, "aprim");
	const Particles ALICEparticles = aprim.particles();
	FastJets FJjets = apply<FastJets>(event, "jets");
	FJjets.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
	Jets jets = FJjets.jetsByPt(); //get jets (ordered by pT)


    }


    /// Normalise histograms etc., after the run
    void finalize() {

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2020_I1755387);

}
