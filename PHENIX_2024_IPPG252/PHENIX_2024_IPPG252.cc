// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2024_IPPG252 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2024_IPPG252);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

// Initialise and register projections
      const FinalState fs(Cuts::abseta < 0.35);
        //In case "NONE" is given as option
      const ParticlePair& beam = beams();

      beamOpt = getOption<string>("beam","NONE");
        

      if (beamOpt == "NONE") {
      if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = pp;
      }

      if (beamOpt =="pp") collSys = pp;
      
      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.3, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");


      // Book histograms
      // specify custom binning
     
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(h["jetpt"], 1, 1, 1);
      book(h["jetcross"], 1, 1, 2);

      book(h["zg910"], 2, 1, 1);
      book(h["zg1012"], 2, 1, 2);
      book(h["zg1214"], 2, 1, 3);
      book(h["zg1417"], 2, 1, 4);
      book(h["zg1720"], 2, 1, 5);
      book(h["zg2024"], 2, 1, 6);
      book(h["zg2429"], 2, 1, 7);
           
      book(h["ETA910"], 3, 1, 1);

      book(h["ETA1012"], 4, 1, 1);
      book(h["ETA1214"], 4, 1, 2);
      book(h["ETA1217"], 4, 1, 3);
      book(h["ETA1720"], 4, 1, 4);
      book(h["ETA2024"], 4, 1, 5);
      book(h["ETA2429"], 4, 1, 6);

      book(h["R910"], 5, 1, 1);
      book(h["R1012"], 5, 1, 2);
      book(h["R1214"], 5, 1, 3);
      book(h["R1417"], 5, 1, 4);
      book(h["R1720"], 5, 1, 5);
      book(h["R2024"], 5, 1, 6);
      book(h["R2429"], 5, 1, 7);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 5*GeV);


    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(h["XXXX"]); // normalize to unity
      //normalize(h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in pb (no cuts)
      //scale(h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> h;
    map<string, Profile1DPtr> p;
    map<string, CounterPtr> c;
    string beamOpt;
    enum CollisionSystem {pp};
    CollisionSystem collSys;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PHENIX_2024_IPPG252);

}
