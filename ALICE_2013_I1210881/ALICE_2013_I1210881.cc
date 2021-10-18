// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2013_I1210881 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2013_I1210881);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // 
      // muons and neutrinos are excluded from the clustering

      FastJets jetfs04(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs04, "jets04");

	FastJets jetfs02(fs, FastJets::ANTIKT, 0.2, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs02, "jets02");

      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["SpectraR0.2"], 1, 1, 1);
      book(_h["SpectraR0.4"], 1, 1, 2);
     
	string refname = mkAxisCode(2, 1, 1);
	book(_s["Ratio"], refname);
	book(_c["sow"], "sow");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

	_c["sow"]->fill();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets04 = apply<FastJets>(event, "jets04").jetsByPt(Cuts::pT > 20*GeV);
	for(auto jet : jets04)
	{
	
		_h["SpectraR0.4"]->fill(jet.pT()/GeV);
		
	}
	 Jets jets02 = apply<FastJets>(event, "jets02").jetsByPt(Cuts::pT > 20*GeV);
        for(auto jet : jets02)
        {

                _h["SpectraR0.2"]->fill(jet.pT()/GeV);

        }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
	_h["SpectraR0.2"]->scaleW(crossSection()/_c["sow"]->sumW());
	_h["SpectraR0.4"]->scaleW(crossSection()/_c["sow"]->sumW());
	divide(_h["SpectraR0.2"], _h["SpectraR0.4"], _s["Ratio"]);

    
}

    ///@}

    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2013_I1210881);

}
