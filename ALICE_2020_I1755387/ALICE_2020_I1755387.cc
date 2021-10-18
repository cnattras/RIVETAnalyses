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
     
      FastJets jet01(fs, FastJets::ANTIKT, 0.1, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet01, "jets01");
      FastJets jet02(fs, FastJets::ANTIKT, 0.2, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet02, "jets02");
      FastJets jet03(fs, FastJets::ANTIKT, 0.3, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet03, "jets03");
      FastJets jet04(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet04, "jets04");
      FastJets jet05(fs, FastJets::ANTIKT, 0.5, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet05, "jets05");
      FastJets jet06(fs, FastJets::ANTIKT, 0.6, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet06, "jets06");

      // Book histograms
      // specify custom binning
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
     
      book(_h["ppspectraR0.1"], 1, 1, 1);
      book(_h["ppspectraR0.2"], 2, 1, 1);
      book(_h["ppspectraR0.3"], 3, 1, 1);
      book(_h["ppspectraR0.4"], 4, 1, 1);
      book(_h["ppspectraR0.5"], 5, 1, 1);
      book(_h["ppspectraR0.6"], 6, 1, 1);
//      book(_s["ppratioR0.1divR0.2"], 13, 1, 1);
      book(_s["ppratioR0.1divR0.3"], 14, 1, 1);
      book(_s["ppratioR0.1divR0.4"], 15, 1, 1);
      book(_s["ppratioR0.1divR0.5"], 16, 1, 1);
      book(_s["ppratioR0.1divR0.6"], 17, 1, 1);
      book(_s["ppratioR0.2divR0.3"], 18, 1, 1);
      book(_s["ppratioR0.2divR0.4"], 19, 1, 1);
      book(_s["ppratioR0.2divR0.5"], 20,1, 1);
      book(_s["ppratioR0.2divR0.6"], 21, 1, 1);
      book(_h["ppspectra5GeVleadtrackR0.2"], 22, 1, 1);
      book(_h["pbspectra5GeVleadtrackR0.2"], 23, 1, 1);
      book(_h["ppspectra7GeVleadtrackR0.4"], 24, 1, 1);
      book(_h["pbspectra7GeVleadtrackR0.4"], 25, 1, 1);
      book(_h["ppcrossleadtrackbias5div0"], 26, 1, 1);
      book(_h["ppcrossleadtrackbias7div0"], 27, 1, 1);
      book(_h["ppcrossleadtrackbias7div5"], 28, 1, 1);
      book(_h["pbcross7to5leadchargeR0.2"], 29, 1, 1);
      book(_h["jetRaa5GeVleadtrackR0.2"], 30, 1, 1);
      book(_h["jetRaa7GeVleadtrackR0.4"], 31, 1, 1);
      book(_h["jetRaa5GeVpbleadtrackR0.2"], 32, 1, 1);
      book(_h["jetRaa7GeVpbleadtrackR0.4"], 33, 1, 1);

string refname = mkAxisCode(13, 1, 1);
const Scatter2D& refdata = refData(refname);
book(_h["ppratioR0.1divR0.2_R0.1"], refname + "_R0.1", refdata);
book(_h["ppratioR0.1divR0.2_R0.2"], refname + "_R0.2", refdata);
book(_s["ppratioR0.1divR0.2"], refname);

      book(_c["sow"], "sow");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
       
        const ALICE::PrimaryParticles aprim = apply<ALICE::PrimaryParticles>(event, "aprim");
	const Particles ALICEparticles = aprim.particles();

	FastJets FJjets01 = apply<FastJets>(event, "jets01");
        FastJets FJjets02 = apply<FastJets>(event, "jets02");
        FastJets FJjets03 = apply<FastJets>(event, "jets03");
        FastJets FJjets04 = apply<FastJets>(event, "jets04");
        FastJets FJjets05 = apply<FastJets>(event, "jets05");
        FastJets FJjets06 = apply<FastJets>(event, "jets06");
 
	FJjets01.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
        FJjets02.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
        FJjets03.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
        FJjets04.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
        FJjets05.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
        FJjets06.calc(ALICEparticles); //give ALICE primary particles to FastJet projection

	Jets jets01 = FJjets01.jetsByPt(Cuts::pT >= 20.*GeV); //get jets (ordered by pT)
        Jets jets02 = FJjets02.jetsByPt(Cuts::pT >= 20.*GeV); //get jets (ordered by pT)
        Jets jets03 = FJjets03.jetsByPt(Cuts::pT >= 20.*GeV); //get jets (ordered by pT)
        Jets jets04 = FJjets04.jetsByPt(Cuts::pT >= 20.*GeV); //get jets (ordered by pT)
        Jets jets05 = FJjets05.jetsByPt(Cuts::pT >= 20.*GeV); //get jets (ordered by pT)
        Jets jets06 = FJjets06.jetsByPt(Cuts::pT >= 20.*GeV); //get jets (ordered by pT)
	
	_c["sow"]->fill();

	for(auto jet : jets01)
	{
		_h["ppspectraR0.1"]->fill(jet.pT()/GeV);
		_h["ppratioR0.1divR0.2_R0.1"]->fill(jet.pT()/GeV);
   	}

        for(auto jet : jets02)
        {
                _h["ppspectraR0.2"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.2_R0.2"]->fill(jet.pT()/GeV);
    }

        for(auto jet : jets03)
        {
                _h["ppspectraR0.3"]->fill(jet.pT()/GeV);
        }

        for(auto jet : jets04)
        {
                _h["ppspectraR0.4"]->fill(jet.pT()/GeV);
        }

        for(auto jet : jets05)
        {
                _h["ppspectraR0.5"]->fill(jet.pT()/GeV);
        }

        for(auto jet : jets06)
        {
                _h["ppspectraR0.6"]->fill(jet.pT()/GeV);
        }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

	_h["ppspectraR0.1"]->scaleW((crossSection()/millibarn)/_c["sow"]->sumW());
        _h["ppspectraR0.2"]->scaleW((crossSection()/millibarn)/_c["sow"]->sumW());
        _h["ppspectraR0.3"]->scaleW((crossSection()/millibarn)/_c["sow"]->sumW());
        _h["ppspectraR0.4"]->scaleW((crossSection()/millibarn)/_c["sow"]->sumW());
        _h["ppspectraR0.5"]->scaleW((crossSection()/millibarn)/_c["sow"]->sumW());
        _h["ppspectraR0.6"]->scaleW((crossSection()/millibarn)/_c["sow"]->sumW());
	divide(_h["ppratioR0.1divR0.2_R0.1"], _h["ppratioR0.1divR0.2_R0.2"], _s["ppratioR0.1divR0.2"]);
/*        divide(_h[], _h["ppspectraR0.3"], _s["ppratioR0.1divR0.3"]);
        divide(_h["ppspectraR0.1"], _h["ppspectraR0.4"], _s["ppratioR0.1divR0.4"]);
        divide(_h["ppspectraR0.1"], _h["ppspectraR0.5"], _s["ppratioR0.1divR0.5"]);
        divide(_h["ppspectraR0.1"], _h["ppspectraR0.6"], _s["ppratioR0.1divR0.6"]);
        divide(_h["ppspectraR0.2"], _h["ppspectraR0.3"], _s["ppratioR0.2divR0.3"]);
        divide(_h["ppspectraR0.2"], _h["ppspectraR0.4"], _s["ppratioR0.2divR0.4"]);
        divide(_h["ppspectraR0.2"], _h["ppspectraR0.5"], _s["ppratioR0.2divR0.5"]);
        divide(_h["ppspectraR0.2"], _h["ppspectraR0.6"], _s["ppratioR0.2divR0.6"]);
*/
	
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


  RIVET_DECLARE_PLUGIN(ALICE_2020_I1755387);

}
