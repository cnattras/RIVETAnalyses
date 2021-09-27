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
      
      // Declaring ALICE primary particles
      const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::abscharge > 0);
      declare(aprim, "aprim");
      // Booking counter (number of events)
      book(_c["sow"], "sow");
      
      // Initialise and register projections
      
      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 0.9);
      
      // Defining map for ratio to be calculated
      map<string, Scatter2DPtr> _s;

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");
      
      
      
      // Book histograms
      // specify custom binning
      book(_h["Jet spectra X"], "myh1", 20, 0., 100.);
      book(_h["Jet spectra Y"], "myh2", logspace(20, 1e-2, 1e3));
      book(_h["JetSpectrum"], "JetSpectrum", 10, 0, 100);
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      //book(_h["AAAA"], 1, 1, 1);
      
      // Calculating ratios
      string refname = mkAxisCode(20., 0., 100.);
      const Scatter2D& refdata = refData(refname);
      book(_h["Numerator"], refname + "Numerator", refdata);
      book(_h["Denominator"], refname + "Denominator", refdata);
      book(_s["Ratio"], refname);
      
      

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Getting the jets
      const ALICE::PrimaryParticles aprim = apply<ALICE::PrimaryParticles>(event, "aprim");
      const Particles ALICEparticles = aprim.particles();
      FastJets FJjets = apply<FastJets>(event, "jets");
      FJjets.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
      Jets jets = FJjets.jetsByPt(Cuts::pT > 20.*GeV && Cuts::abseta < 0.5); //get jets (ordered by pT), only above 20GeV
      
      _c["sow"]->fill();
      
      for(auto jet : jets){
      	// Normalizing histogram by number of events
      	_h["JetSpectrum"]->fill(jet.pT()/GeV);
      	
      	
      }
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      //Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 0*GeV && Cuts::abseta < 0.5);
      

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      
      //_h["JetSpectrum"]->scaleW(crossSection()/_c["sow"]->sumW());
      _h["JetSpectrum"]->scaleW(1./_c["sow"]->sumW());
      
      // Defining map for ratio to be calculated
      map<string, Scatter2DPtr> _s;
      // Ratio
      divide(_h["Numerator"], _h["Denominator"], _s["Ratio"]);

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2021_I0000000);

}
