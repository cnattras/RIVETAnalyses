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
     
      // Centralities
	declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");

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
//      book(_s["ppratioR0.1divR0.3"], 14, 1, 1);
//      book(_s["ppratioR0.1divR0.4"], 15, 1, 1);
//      book(_s["ppratioR0.1divR0.5"], 16, 1, 1);
//      book(_s["ppratioR0.1divR0.6"], 17, 1, 1);
//      book(_s["ppratioR0.2divR0.3"], 18, 1, 1);
//      book(_s["ppratioR0.2divR0.4"], 19, 1, 1);
//      book(_s["ppratioR0.2divR0.5"], 20,1, 1);
//      book(_s["ppratioR0.2divR0.6"], 21, 1, 1);
      book(_h["ppspectra5GeVleadtrackR0.2"], 22, 1, 1);
      book(_h["pbspectra5GeVleadtrackR0.2"], 23, 1, 1);
      book(_h["ppspectra7GeVleadtrackR0.4"], 24, 1, 1);
      book(_h["pbspectra7GeVleadtrackR0.4"], 25, 1, 1);
//      book(_s["ppcrossleadtrackbias5div0"], 26, 1, 1);
//      book(_s["ppcrossleadtrackbias7div0"], 27, 1, 1);
      book(_s["ppcrossleadtrackbias7div5"], 28, 1, 1);
      book(_s["pbcross7to5leadchargeR0.2"], 29, 1, 1);
      book(_s["jetRaa5GeVleadtrackR0.2"], 30, 1, 1);
      book(_s["jetRaa7GeVleadtrackR0.4"], 31, 1, 1);
      book(_s["jetRaa5GeVpbleadtrackR0.2"], 32, 1, 1);
      book(_s["jetRaa7GeVpbleadtrackR0.4"], 33, 1, 1);

string refname13 = mkAxisCode(13, 1, 1);
const Scatter2D& refdata13 = refData(refname13);
book(_h["ppratioR0.1divR0.2_R0.1"], refname13 + "_R0.1", refdata13);
book(_h["ppratioR0.1divR0.2_R0.2"], refname13 + "_R0.2", refdata13);
book(_s["ppratioR0.1divR0.2"], refname13);

string refname14 = mkAxisCode(14, 1, 1);
const Scatter2D& refdata14 = refData(refname14);
book(_h["ppratioR0.1divR0.3_R0.1"], refname14 + "_R0.1", refdata14);
book(_h["ppratioR0.1divR0.3_R0.3"], refname14 + "_R0.3", refdata14);
book(_s["ppratioR0.1divR0.3"], refname14);

string refname15 = mkAxisCode(15, 1, 1);
const Scatter2D& refdata15 = refData(refname15);
book(_h["ppratioR0.1divR0.4_R0.1"], refname15 + "_R0.1", refdata15);
book(_h["ppratioR0.1divR0.4_R0.4"], refname15 + "_R0.4", refdata15);
book(_s["ppratioR0.1divR0.4"], refname15);

string refname16 = mkAxisCode(16, 1, 1);
const Scatter2D& refdata16 = refData(refname16);
book(_h["ppratioR0.1divR0.5_R0.1"], refname16 + "_R0.1", refdata16);
book(_h["ppratioR0.1divR0.5_R0.5"], refname16 + "_R0.5", refdata16);
book(_s["ppratioR0.1divR0.5"], refname16);

string refname17 = mkAxisCode(17, 1, 1);
const Scatter2D& refdata17 = refData(refname17);
book(_h["ppratioR0.1divR0.6_R0.1"], refname17 + "_R0.1", refdata17);
book(_h["ppratioR0.1divR0.6_R0.6"], refname17 + "_R0.6", refdata17);
book(_s["ppratioR0.1divR0.6"], refname17);

string refname18 = mkAxisCode(18, 1, 1);
const Scatter2D& refdata18 = refData(refname18);
book(_h["ppratioR0.2divR0.3_R0.2"], refname18 + "_R0.2", refdata18);
book(_h["ppratioR0.2divR0.3_R0.3"], refname18 + "_R0.3", refdata18);
book(_s["ppratioR0.2divR0.3"], refname18);

string refname19 = mkAxisCode(19, 1, 1);
const Scatter2D& refdata19 = refData(refname19);
book(_h["ppratioR0.2divR0.4_R0.2"], refname19 + "_R0.2", refdata19);
book(_h["ppratioR0.2divR0.4_R0.4"], refname19 + "_R0.4", refdata19);
book(_s["ppratioR0.2divR0.4"], refname19);

string refname20 = mkAxisCode(20, 1, 1);
const Scatter2D& refdata20 = refData(refname20);
book(_h["ppratioR0.2divR0.5_R0.2"], refname20 + "_R0.2", refdata20);
book(_h["ppratioR0.2divR0.5_R0.5"], refname20 + "_R0.5", refdata20);
book(_s["ppratioR0.2divR0.5"], refname20);

string refname21 = mkAxisCode(21, 1, 1);
const Scatter2D& refdata21 = refData(refname21);
book(_h["ppratioR0.2divR0.6_R0.2"], refname21 + "_R0.2", refdata21);
book(_h["ppratioR0.2divR0.6_R0.6"], refname21 + "_R0.6", refdata21);
book(_s["ppratioR0.2divR0.6"], refname21);

string refname26 = mkAxisCode(26, 1, 1);
const Scatter2D& refdata26 = refData(refname26);
book(_h["ppcrossleadtrackbias5"], refname26 + "_5GeV_R0.2", refdata26);
book(_s["ppcrossleadtrackbias5div0"], refname26);

string refname27 = mkAxisCode(27, 1, 1);
const Scatter2D& refdata27 = refData(refname27);
book(_h["ppcrossleadtrackbias7"], refname27 + "_7GeV_R0.2", refdata27);
book(_s["ppcrossleadtrackbias7div0"], refname27);


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
                _h["ppratioR0.1divR0.3_R0.1"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.4_R0.1"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.5_R0.1"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.6_R0.1"]->fill(jet.pT()/GeV);

	}

        for(auto jet : jets02)
        {
		if(jet.particles(Cuts::pT > 5.*GeV).size() > 0){
			_h["ppcrossleadtrackbias5"]->fill(jet.pT()/GeV);
 			_h["ppspectra5GeVleadtrackR0.2"]->fill(jet.pT()/GeV);
		}
                if(jet.particles(Cuts::pT > 7.*GeV).size() > 0){
                        _h["ppcrossleadtrackbias7"]->fill(jet.pT()/GeV);
                        _h["ppspectra7GeVleadtrackR0.2"]->fill(jet.pT()/GeV);
               	}
	        _h["ppspectraR0.2"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.2_R0.2"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.3_R0.2"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.4_R0.2"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.5_R0.2"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.6_R0.2"]->fill(jet.pT()/GeV);

    }

        for(auto jet : jets03)
        {
                _h["ppspectraR0.3"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.3_R0.3"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.3_R0.3"]->fill(jet.pT()/GeV);

	 }

        for(auto jet : jets04)
        {
                if(jet.particles(Cuts::pT > 5.*GeV).size() > 0)
                        _h["ppspectra5GeVleadtrackR0.4"]->fill(jet.pT()/GeV);
                if(jet.particles(Cuts::pT > 7.*GeV).size() > 0)
                        _h["ppspectra7GeVleadtrackR0.4"]->fill(jet.pT()/GeV);

		 _h["ppspectraR0.4"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.4_R0.4"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.4_R0.4"]->fill(jet.pT()/GeV);

        }

        for(auto jet : jets05)
        {
                _h["ppspectraR0.5"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.5_R0.5"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.5_R0.5"]->fill(jet.pT()/GeV);

	 }

        for(auto jet : jets06)
        {
                _h["ppspectraR0.6"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.6_R0.6"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.6_R0.6"]->fill(jet.pT()/GeV);
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
        divide(_h["ppratioR0.1divR0.3_R0.1"], _h["ppratioR0.1divR0.3_R0.3"], _s["ppratioR0.1divR0.3"]);
        divide(_h["ppratioR0.1divR0.4_R0.1"], _h["ppratioR0.1divR0.4_R0.4"], _s["ppratioR0.1divR0.4"]);
        divide(_h["ppratioR0.1divR0.5_R0.1"], _h["ppratioR0.1divR0.5_R0.5"], _s["ppratioR0.1divR0.5"]);
        divide(_h["ppratioR0.1divR0.6_R0.1"], _h["ppratioR0.1divR0.6_R0.6"], _s["ppratioR0.1divR0.6"]);
        divide(_h["ppratioR0.2divR0.3_R0.2"], _h["ppratioR0.2divR0.3_R0.3"], _s["ppratioR0.2divR0.3"]);
        divide(_h["ppratioR0.2divR0.4_R0.2"], _h["ppratioR0.2divR0.4_R0.4"], _s["ppratioR0.2divR0.4"]);
        divide(_h["ppratioR0.2divR0.5_R0.2"], _h["ppratioR0.2divR0.5_R0.5"], _s["ppratioR0.2divR0.5"]);
        divide(_h["ppratioR0.2divR0.6_R0.2"], _h["ppratioR0.2divR0.6_R0.6"], _s["ppratioR0.2divR0.6"]);

        _h["ppcrossleadtrackbias5"]->scaleW((crossSection()/millibarn)/_c["sow"]->sumW());
        divide(_h["ppcrossleadtrackbias5"], _h["ppspectraR0.2"], _s["ppcrossleadtrackbias5div0"]);
        _h["ppcrossleadtrackbias7"]->scaleW((crossSection()/millibarn)/_c["sow"]->sumW());
        divide(_h["ppcrossleadtrackbias7"], _h["ppspectraR0.2"], _s["ppcrossleadtrackbias7div0"]);
        divide(_h["ppcrossleadtrackbias7"], _h["ppcrossleadtrackbias5"], _s["ppcrossleadtrackbias7div5"]);
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
