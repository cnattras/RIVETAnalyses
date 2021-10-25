// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
//#include "Rivet/Projections/DressedLeptons.hh"
//#include "Rivet/Projections/MissingMomentum.hh"
//#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2021_I0000000 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2021_I0000000);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      
      
      // Initialise and register projections
      // Declaring ALICE primary particles
      const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::abscharge > 0);
      declare(aprim, "aprim");
      
      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 0.9);
      

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.2, 0.3, 0.4, 0.5, 0.6, 0.7
      // muons and neutrinos are excluded from the clustering
      
      FastJets JetPtR02(fs, FastJets::ANTIKT, 0.2, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(JetPtR02, "JetPtR02");
      FastJets JetPtR03(fs, FastJets::ANTIKT, 0.3, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(JetPtR03, "JetPtR03");
      FastJets JetPtR04(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(JetPtR04, "JetPtR04");
      FastJets JetPtR05(fs, FastJets::ANTIKT, 0.5, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(JetPtR05, "JetPtR05");
      FastJets JetPtR06(fs, FastJets::ANTIKT, 0.6, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(JetPtR06, "JetPtR06");
      FastJets JetPtR07(fs, FastJets::ANTIKT, 0.7, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(JetPtR07, "JetPtR07");
      
      // Book histograms
      // specify custom binning
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["JetPtR02"], 1, 1, 1); //Figure 3 R = 0.2
      book(_h["JetPtR03"], 1, 1, 2); //Figure 3 R = 0.3
      book(_h["JetPtR04"], 1, 1, 3); //Figure 3 R = 0.4
      book(_h["JetPtR05"], 1, 1, 4); //Figure 3 R = 0.5
      book(_h["JetPtR06"], 1, 1, 5); //Figure 3 R = 0.6
      book(_h["JetPtR07"], 1, 1, 6); //Figure 3 R = 0.7
      
      
      
      // Calculating ratios
      /*string refname1 = mkAxisCode(20., 0., 100.); //(table,1,1)
      book(_s["Ratio1"], refname1);
      
      //const Scatter2D& refdata = refData(refname); //am I using refData?
      //book(_h["Numerator"], refname + "Numerator", refdata);
      //book(_h["Denominator"], refname + "Denominator", refdata);
      */
      
      // Booking counter (number of events)
      book(_c["sow"], "sow");
      

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Getting the jets
      const ALICE::PrimaryParticles aprim = apply<ALICE::PrimaryParticles>(event, "aprim");
      const Particles ALICEparticles = aprim.particles();
      
      // R = 0.2
      FastJets FJjetsR02 = apply<FastJets>(event, "JetPtR02");
      FJjetsR02.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
      
      // R = 0.3
      FastJets FJjetsR03 = apply<FastJets>(event, "JetPtR02");
      FJjetsR03.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
      
      // R = 0.4
      FastJets FJjetsR04 = apply<FastJets>(event, "JetPtR02");
      FJjetsR04.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
      
      // R = 0.5
      FastJets FJjetsR05 = apply<FastJets>(event, "JetPtR02");
      FJjetsR05.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
      
      // R = 0.6
      FastJets FJjetsR06 = apply<FastJets>(event, "JetPtR02");
      FJjetsR06.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
      
      // R = 0.7
      FastJets FJjetsR07 = apply<FastJets>(event, "JetPtR02");
      FJjetsR07.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
      
      
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      // R = 0.2
      Jets JetPtR02 = FJjetsR02.jetsByPt(Cuts::pT > 20.*GeV && Cuts::abseta < 0.5); //get jets (ordered by pT), only above 20GeV
      for(auto jet : JetPtR02){
      	// Normalizing histogram by number of events
      	_h["JetPtR02"]->fill(jet.pT()/GeV);
      }
      // R = 0.3
      Jets JetPtR03 = FJjetsR03.jetsByPt(Cuts::pT > 20.*GeV && Cuts::abseta < 0.5);
      for(auto jet : JetPtR03){
      	// Normalizing histogram by number of events
      	_h["JetPtR03"]->fill(jet.pT()/GeV);
      }
      // R = 0.4
      Jets JetPtR04 = FJjetsR04.jetsByPt(Cuts::pT > 20.*GeV && Cuts::abseta < 0.5);
      for(auto jet : JetPtR04){
      	// Normalizing histogram by number of events
      	_h["JetPtR04"]->fill(jet.pT()/GeV);
      }
      // R = 0.5
      Jets JetPtR05 = FJjetsR05.jetsByPt(Cuts::pT > 20.*GeV && Cuts::abseta < 0.5);
      for(auto jet : JetPtR05){
      	// Normalizing histogram by number of events
      	_h["JetPtR05"]->fill(jet.pT()/GeV);
      }
      // R = 0.6
      Jets JetPtR06 = FJjetsR06.jetsByPt(Cuts::pT > 20.*GeV && Cuts::abseta < 0.5);
      for(auto jet : JetPtR06){
      	// Normalizing histogram by number of events
      	_h["JetPtR06"]->fill(jet.pT()/GeV);
      }
      // R = 0.7
      Jets JetPtR07 = FJjetsR07.jetsByPt(Cuts::pT > 20.*GeV && Cuts::abseta < 0.5);
      for(auto jet : JetPtR07){
      	// Normalizing histogram by number of events
      	_h["JetPtR07"]->fill(jet.pT()/GeV);
      }
      
      // Sum of weights counter
      _c["sow"]->fill();
      
      
      
      

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      
      //_h["JetPtR02"]->scaleW(1./_c["sow"]->sumW());
      _h["JetPtR02"]->scaleW((crossSection()/millibarn)/(_c["sow"]->sumW()));
      _h["JetPtR03"]->scaleW((crossSection()/millibarn)/(_c["sow"]->sumW()));
      _h["JetPtR04"]->scaleW((crossSection()/millibarn)/(_c["sow"]->sumW()));
      _h["JetPtR05"]->scaleW((crossSection()/millibarn)/(_c["sow"]->sumW()));
      _h["JetPtR06"]->scaleW((crossSection()/millibarn)/(_c["sow"]->sumW()));
      _h["JetPtR07"]->scaleW((crossSection()/millibarn)/(_c["sow"]->sumW()));
      
      
      
      // Ratio
      //divide(_h["Numerator"], _h["Denominator"], _s["Ratio"]);

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    // Defining map for ratio to be calculated
    //map<string, Scatter2DPtr> _s;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2021_I0000000);

}
