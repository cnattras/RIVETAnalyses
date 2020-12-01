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
  class PHENIX_2016_I1394433 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2016_I1394433);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      declareCentrality(RHICCentrality("PHENIX"),"RHIC_2019_CentralityCalibration:exp=PHENIX","CMULT","CMULT");
    

    // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);
      book(_h["hist_E_200"],7,1,3);
      book(_h["hist_Ch_200"],8,1,3);
      book(_h["hist_E_130"],9,1,3);
      book(_h["hist_Ch_130"],10,1,3);
      book(_h["hist_E_62.4"],11,1,3);
      book(_h["hist_Ch_62.4"],12,1,3);
      book(_h["hist_E_39"],13,1,3);
      book(_h["hist_Ch_39"],14,1,3);
      book(_h["hist_E_27"],15,1,3);
      book(_h["hist_Ch_27"],16,1,3);
      book(_h["hist_E_19.6"],17,1,3);
      book(_h["hist_Ch_19.6"],18,1,3);
      book(_h["hist_E_14.5"],19,1,3);
      book(_h["hist_Ch_14.5"],20,1,3);
      book(_h["hist_E_7.7"],21,1,3);
      book(_h["hist_Ch_7.7"],22,1,3);
      
      book(_h["hist_E_200_Cu"],23,1,3);
      book(_h["hist_Ch_200_Cu"],24,1,3);
      book(_h["hist_E_62.4_Cu"],25,1,3);
      book(_h["hist_Ch_62.4_Cu"],26,1,3);
      book(_h["hist_E_200_Cu_Au"],27,1,3);
      book(_h["hist_Ch_200_Cu_Au"],28,1,3);
//      book(_h["hist_E_193_UU"],29,1,3);
//      book(_h["hist_Ch_193_UU"],30,1,3);
//      book(_h["hist_E_200_dAu"],31,1,3);
//      book(_h["hist_Ch_200_dAu"],32,1,3);
//      book(_h["hist_E_200_He_A"],33,1,3);
//      book(_h["hist_Ch_200_He_A"],34,1,3);

/*      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning
      book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
      book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["AAAA"], 1, 1, 1);
      book(_p["BBBB"], 2, 1, 1);
      book(_c["CCCC"], 3, 1, 1);
*/
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
/*
      // Retrieve dressed leptons, sorted by pT
      vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

      // Remove all jets within dR < 0.2 of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

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
*/
    }


    /// Normalise histograms etc., after the run
    void finalize() {
/*
      normalize(_h["XXXX"]); // normalize to unity
      normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
      scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)
*/
    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2016_I1394433);

}
