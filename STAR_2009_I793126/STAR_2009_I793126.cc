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
  class STAR_2009_I793126 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2009_I793126);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);

      // The final-state particles declared above are clustered using FastJet with
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
      //book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
      //book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      //book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      //book(_h["AAAA"], 1, 1, 1);
      //book(_p["BBBB"], 2, 1, 1);
      //book(_c["CCCC"], 3, 1, 1);

      //Figures from the paper
      //Figure 1 AuAu200
       book(_h["Figure_1_AuAu200"], 1, 1, 1);
      //Figure 1 AuAu62
       book(_h["Figure_1_AuAu62"], 2, 1, 1);
      //Figure 2a
       book(_h["Figure_2a"], 3, 1, 1);
      //Figure 2b
       book(_h["Figure_2b_1"], 4, 1, 1);
       book(_h["Figure_2b_2"], 4, 1, 2);
       book(_h["Figure_2b_3"], 4, 1, 3);
       book(_h["Figure_2b_4"], 4, 1, 4);
      //Figure 3
       book(_h["Figure_3"], 5, 1, 1);
      //Figure 4a AuAu200
       book(_h["Figure_4a_AuAu200_1"], 6, 1, 1);
       book(_h["Figure_4a_AuAu200_2"], 6, 1, 2);
      //Figure 4a AuAu62
       book(_h["Figure_4a_AuAu62_1"], 7, 1, 1);
       book(_h["Figure_4a_AuAu62_2"], 7, 1, 2);
      //Figure 4a pp
       book(_h["Figure_4a_pp"], 8, 1, 1);
      //Figure 4b AuAu200
       book(_h["AAAA"], 9, 1, 1);
       book(_h["AAAA"], 9, 1, 2);
      //Figure 4b AuAu62
       book(_h["AAAA"], 10, 1, 1);
       book(_h["AAAA"], 10, 1, 2);
      //Figure 4b pp
       book(_h["AAAA"], 11, 1, 1);
      //Figure 4b AuAu200
      //Figure 4b AuAu62
      //Figure 10
       book(_h["AAAA"], 12, 1, 1);
       book(_h["AAAA"], 12, 1, 2);
       book(_h["AAAA"], 12, 1, 3);
       book(_h["AAAA"], 12, 1, 4);
       book(_h["AAAA"], 12, 1, 5);
       book(_h["AAAA"], 12, 1, 6);
      //Figure 11a
       book(_h["AAAA"], 13, 1, 1);
       book(_h["AAAA"], 13, 1, 2);
      //Figure 11b
     //book(_h["AAAA"], 14, 1, 1);
     //book(_h["AAAA"], 14, 1, 2);
      //Figure 12 dAu
      //book(_h["AAAA"], 15, 1, 1);
      //Figure 12 pp
      //book(_h["AAAA"], 16, 1, 1);
      //Figure 13a
      //book(_h["AAAA"], 17, 1, 1);
     //book(_h["AAAA"], 17, 1, 2);
     //book(_h["AAAA"], 17, 1, 3);
      //Figure 13b
     //book(_h["AAAA"], 18, 1, 1);
     //book(_h["AAAA"], 18, 1, 2);
     //book(_h["AAAA"], 18, 1, 3);
      //Figure 14a pion and kaon
     //book(_h["AAAA"], 19, 1, 1);
     //book(_h["AAAA"], 19, 1, 2);
      //Figure 14a proton
     //book(_h["AAAA"], 20, 1, 1);
      //Figure 14b pion and kaon
     //book(_h["AAAA"], 21, 1, 1);
     //book(_h["AAAA"], 21, 1, 2);
      //Figure 14b proton
     //book(_h["AAAA"], 22, 1, 1);
      //Figure 15a
     //book(_h["AAAA"], 23, 1, 1);
     //book(_h["AAAA"], 23, 1, 2);
     //book(_h["AAAA"], 23, 1, 3);
     //book(_h["AAAA"], 23, 1, 4);
      //Figure 15a Bkgd
     //book(_h["AAAA"], 24, 1, 1);
      //Figure 15b
     //book(_h["AAAA"], 25, 1, 1);
     //book(_h["AAAA"], 25, 1, 2);
     //book(_h["AAAA"], 25, 1, 3);
     //book(_h["AAAA"], 25, 1, 4);
      //Figure 15b Bkgd
     //book(_h["AAAA"], 26, 1, 1);
      //Figure 15c
     //book(_h["AAAA"], 27, 1, 1);
     //book(_h["AAAA"], 27, 1, 2);
     //book(_h["AAAA"], 27, 1, 3);
     //book(_h["AAAA"], 27, 1, 4);
      //Figure 15c Bkgd
     //book(_h["AAAA"], 28, 1, 1);
      //Figure 15d
     //book(_h["AAAA"], 29, 1, 1);
     //book(_h["AAAA"], 29, 1, 2);
     //book(_h["AAAA"], 29, 1, 3);
     //book(_h["AAAA"], 29, 1, 4);
      //Figure 15d Bkgd
     //book(_h["AAAA"], 30, 1, 1);
      //Figure 16 all and weak-decay bkgd
     //book(_h["AAAA"], 31, 1, 1);
     //book(_h["AAAA"], 31, 1, 2);
      //Figure 16 muon contamination
     //book(_h["AAAA"], 32, 1, 1);
      //Figure 17 AuAu dE/dx
     //book(_h["AAAA"], 33, 1, 1);
      //Figure 17 AuAu Blast-wave fit
     //book(_h["AAAA"], 34, 1, 1);
      //Figure 17 AuAu p_T-Gaussian fit
     //book(_h["AAAA"], 35, 1, 1);
      //Figure 17 AuAu p_T-exponential fit
     //book(_h["AAAA"], 36, 1, 1);
      //Figure 17 AuAu TOF data
     //book(_h["AAAA"], 37, 1, 1);
      //Figure 17 dAu dE/dx
     //book(_h["AAAA"], 38, 1, 1);
      //Figure 17 dAu Blast-wave fit
     //book(_h["AAAA"], 39, 1, 1);
      //Figure 17 dAu p_T-Gaussian fit
     //book(_h["AAAA"], 40, 1, 1);
      //Figure 17 dAu p_T-exponential fit
     //book(_h["AAAA"], 41, 1, 1);
      //Figure 17 dAu TOF data
     //book(_h["AAAA"], 42, 1, 1);
      //Figure 18 kaon
     //book(_h["AAAA"], 43, 1, 1);
     //book(_h["AAAA"], 43, 1, 2);
     //book(_h["AAAA"], 43, 1, 3);
     //book(_h["AAAA"], 43, 1, 4);
     //book(_h["AAAA"], 43, 1, 5);
     //book(_h["AAAA"], 43, 1, 6);
     //book(_h["AAAA"], 43, 1, 7);
     //book(_h["AAAA"], 43, 1, 8);
      //Figure 18 pion
     //book(_h["AAAA"], 44, 1, 1);
     //book(_h["AAAA"], 44, 1, 2);
     //book(_h["AAAA"], 44, 1, 3);
     //book(_h["AAAA"], 44, 1, 4);
     //book(_h["AAAA"], 44, 1, 5);
     //book(_h["AAAA"], 44, 1, 6);
     //book(_h["AAAA"], 44, 1, 7);
     //book(_h["AAAA"], 44, 1, 8);
      //Figure 18 proton
     //book(_h["AAAA"], 45, 1, 1);
     //book(_h["AAAA"], 45, 1, 2);
     //book(_h["AAAA"], 45, 1, 3);
     //book(_h["AAAA"], 45, 1, 4);
     //book(_h["AAAA"], 45, 1, 5);
     //book(_h["AAAA"], 45, 1, 6);
     //book(_h["AAAA"], 45, 1, 7);
     //book(_h["AAAA"], 45, 1, 8);
      //Figure 19 kaon
     //book(_h["AAAA"], 46, 1, 1);
     //book(_h["AAAA"], 46, 1, 2);
     //book(_h["AAAA"], 46, 1, 3);
     //book(_h["AAAA"], 46, 1, 4);
     //book(_h["AAAA"], 46, 1, 5);
     //book(_h["AAAA"], 46, 1, 6);
     //book(_h["AAAA"], 46, 1, 7);
     //book(_h["AAAA"], 46, 1, 8);
     //book(_h["AAAA"], 46, 1, 9);
     //book(_h["AAAA"], 46, 1, 10);
     //book(_h["AAAA"], 46, 1, 11);
     //book(_h["AAAA"], 46, 1, 12);
     //book(_h["AAAA"], 46, 1, 13);
     //book(_h["AAAA"], 46, 1, 14);
     //book(_h["AAAA"], 46, 1, 15);
     //book(_h["AAAA"], 46, 1, 16);
     //book(_h["AAAA"], 46, 1, 17);
     //book(_h["AAAA"], 46, 1, 18);
      //Figure 19 pion
     //book(_h["AAAA"], 47, 1, 1);
     //book(_h["AAAA"], 47, 1, 2);
     //book(_h["AAAA"], 47, 1, 3);
     //book(_h["AAAA"], 47, 1, 4);
     //book(_h["AAAA"], 47, 1, 5);
     //book(_h["AAAA"], 47, 1, 6);
     //book(_h["AAAA"], 47, 1, 7);
     //book(_h["AAAA"], 47, 1, 8);
     //book(_h["AAAA"], 47, 1, 9);
     //book(_h["AAAA"], 47, 1, 10);
     //book(_h["AAAA"], 47, 1, 11);
     //book(_h["AAAA"], 47, 1, 12);
     //book(_h["AAAA"], 47, 1, 13);
     //book(_h["AAAA"], 47, 1, 14);
     //book(_h["AAAA"], 47, 1, 15);
     //book(_h["AAAA"], 47, 1, 16);
     //book(_h["AAAA"], 47, 1, 17);
     //book(_h["AAAA"], 47, 1, 18);
      //Figure 19 proton
     //book(_h["AAAA"], 48, 1, 1);
     //book(_h["AAAA"], 48, 1, 2);
     //book(_h["AAAA"], 48, 1, 3);
     //book(_h["AAAA"], 48, 1, 4);
     //book(_h["AAAA"], 48, 1, 5);
     //book(_h["AAAA"], 48, 1, 6);
     //book(_h["AAAA"], 48, 1, 7);
     //book(_h["AAAA"], 48, 1, 8);
     //book(_h["AAAA"], 48, 1, 9);
     //book(_h["AAAA"], 48, 1, 10);
     //book(_h["AAAA"], 48, 1, 11);
     //book(_h["AAAA"], 48, 1, 12);
     //book(_h["AAAA"], 48, 1, 13);
     //book(_h["AAAA"], 48, 1, 14);
     //book(_h["AAAA"], 48, 1, 15);
     //book(_h["AAAA"], 48, 1, 16);
     //book(_h["AAAA"], 48, 1, 17);
     //book(_h["AAAA"], 48, 1, 18);
      //Figure 20
     //book(_h["AAAA"], 49, 1, 1);
     //book(_h["AAAA"], 49, 1, 2);
     //book(_h["AAAA"], 49, 1, 3);
     //book(_h["AAAA"], 49, 1, 4);
     //book(_h["AAAA"], 49, 1, 5);
     //book(_h["AAAA"], 49, 1, 6);
     //book(_h["AAAA"], 49, 1, 7);
     //book(_h["AAAA"], 49, 1, 8);
     //book(_h["AAAA"], 49, 1, 9);
     //book(_h["AAAA"], 49, 1, 10);
     //book(_h["AAAA"], 49, 1, 11);
     //book(_h["AAAA"], 49, 1, 12);
     //book(_h["AAAA"], 49, 1, 13);
     //book(_h["AAAA"], 49, 1, 14);
     //book(_h["AAAA"], 49, 1, 15);
     //book(_h["AAAA"], 49, 1, 16);
      //Figure 24 Au+Au 62.4 GeV
     //book(_h["AAAA"], 50, 1, 1);
     //book(_h["AAAA"], 50, 1, 2);
      //Figure 24 Au+Au 130 GeV
     //book(_h["AAAA"], 51, 1, 1);
     //book(_h["AAAA"], 51, 1, 2);
      //Figure 24 Au+Au 200 GeV
     //book(_h["AAAA"], 52, 1, 1);
     //book(_h["AAAA"], 52, 1, 2);
      //Figure 30a
     //book(_h["AAAA"], 53, 1, 1);
      //Figure 30b
     //book(_h["AAAA"], 54, 1, 1);
      //Figure 32 AGS E859 Si+Al 5.4 GeV
     //book(_h["AAAA"], 55, 1, 1);
      //Figure 32 AGS E866 Au+Au 4.7 GeV
     //book(_h["AAAA"], 56, 1, 1);
      //Figure 32 SPS NA49 Pb+Pb 17.3 GeV
     //book(_h["AAAA"], 57, 1, 1);
      //Figure 32 SPS NA49 Pb+Pb energy scan
     //book(_h["AAAA"], 58, 1, 1);
      //Figure 32 SPS NA49 S+S 20 GeV
     //book(_h["AAAA"], 59, 1, 1);
      //Figure 32 SPS NA49 C+C/Si+Si 17.3 GeV
     //book(_h["AAAA"], 60, 1, 1);
      //Figure 32 STAR Au+Au 62.4 GeV
     //book(_h["AAAA"], 61, 1, 1);
      //Figure 32 STAR Au+Au 130 GeV
     //book(_h["AAAA"], 62, 1, 1);
      //Figure 32 STAR Au-Au 200 GeV
     //book(_h["AAAA"], 63, 1, 1);
      //Figure 33 AGS E859 Si+Al 5.4 GeV
     //book(_h["AAAA"], 64, 1, 1);
      //Figure 33 AGS E866 Au+Au 4.7 GeV
     //book(_h["AAAA"], 65, 1, 1);
      //Figure 33 SPS NA49 Pb+Pb 17.3 GeV
     //book(_h["AAAA"], 66, 1, 1);
      //Figure 33 SPS NA49 Pb+Pb energy scan
     //book(_h["AAAA"], 67, 1, 1);
      //Figure 33 SPS NA49 S+S 20 GeV
     //book(_h["AAAA"], 68, 1, 1);
      //Figure 33 SPS NA49 C+C/Si+Si 17.3 GeV
     //book(_h["AAAA"], 69, 1, 1);
      //Figure 33 STAR Au+Au 62.4 GeV
     //book(_h["AAAA"], 70, 1, 1);
      //Figure 33 STAR Au+Au 130 GeV
     //book(_h["AAAA"], 71, 1, 1);
      //Figure 33 STAR Au+Au 200 GeV
     //book(_h["AAAA"], 72, 1, 1);
      //Figure 38 Becattini et al.
     //book(_h["AAAA"], 73, 1, 1);
      //Figure 38 Andronic et al.
     //book(_h["AAAA"], 74, 1, 1);
      //Figure 38 SIS
     //book(_h["AAAA"], 75, 1, 1);
      //Figure 38 AGS
     //book(_h["AAAA"], 76, 1, 1);
      //Figure 38 SPS
     //book(_h["AAAA"], 77, 1, 1);
      //Figure 38 STAR
     //book(_h["AAAA"], 78, 1, 1);
      //Figure 39 Elementary collisions Becattini et al.
     //book(_h["AAAA"], 79, 1, 1);
      //Figure 39 Becattini et al.
     //book(_h["AAAA"], 80, 1, 1);
      //Figure 39 Andronic et al.
     //book(_h["AAAA"], 81, 1, 1);
      //Figure 39 SIS
     //book(_h["AAAA"], 82, 1, 1);
      //Figure 39 AGS
     //book(_h["AAAA"], 83, 1, 1);
      //Figure 39 SPS
     //book(_h["AAAA"], 84, 1, 1);
      //Figure 39 STAR pp
     //book(_h["AAAA"], 85, 1, 1);
      //Figure 39 STAR Tchem
     //book(_h["AAAA"], 86, 1, 1);
      //Figure 39 EOS
     //book(_h["AAAA"], 87, 1, 1);
      //Figure 39 FOPI
     //book(_h["AAAA"], 88, 1, 1);
      //Figure 39 E866
     //book(_h["AAAA"], 89, 1, 1);
      //Figure 39 NA49
     //book(_h["AAAA"], 90, 1, 1);
      //Figure 39 STAR Tkin
     //book(_h["AAAA"], 91, 1, 1);
      //Figure 40 FOPI
     //book(_h["AAAA"], 92, 1, 1);
      //Figure 40 EOS
     //book(_h["AAAA"], 93, 1, 1);
      //Figure 40 E866
     //book(_h["AAAA"], 94, 1, 1);
      //Figure 40 NA49
     //book(_h["AAAA"], 95, 1, 1);
      //Figure 40 STAR
     //book(_h["AAAA"], 96, 1, 1);
      //Figure 41 mu Andronic et al.
     //book(_h["AAAA"], 97, 1, 1);
      //Figure 41 mu SIS
     //book(_h["AAAA"], 98, 1, 1);
      //Figure 41 mu AGS 4.8 GeV
     //book(_h["AAAA"], 99, 1, 1);
      //Figure 41 mu SPS
     //book(_h["AAAA"], 100, 1, 1);
      //Figure 41 mu STAR Au+Au
     //book(_h["AAAA"], 101, 1, 1);
      //Figure 41 mu STAR pp 200 GeV
     //book(_h["AAAA"], 102, 1, 1);
      //Figure 41 Tch Andronic et al.
     //book(_h["AAAA"], 103, 1, 1);
      //Figure 41 Tch SIS
     //book(_h["AAAA"], 104, 1, 1);
      //Figure 41 Tch AGS 4.8 GeV
     //book(_h["AAAA"], 105, 1, 1);
      //Figure 41 Tch SPS
     //book(_h["AAAA"], 106, 1, 1);
      //Figure 41 Tch STAR Au+Au
     //book(_h["AAAA"], 107, 1, 1);
      //Figure 41 Tch STAR pp 200 GeV
     //book(_h["AAAA"], 108, 1, 1);
      //Figure 42 Optical Glauber
     //book(_h["AAAA"], 109, 1, 1);
      //Figure 42 MC Glauber
     //book(_h["AAAA"], 110, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

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


  DECLARE_RIVET_PLUGIN(STAR_2009_I793126);

}
