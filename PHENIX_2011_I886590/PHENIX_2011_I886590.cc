// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2011_I886590 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2011_I886590);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
cout << "Made it into initialize" << endl;
      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9 && Cuts::pT > 0.15*GeV);
      declare(fs, "fs");

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      //FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      //declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      //PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      //PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      //Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      //DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      //declare(dressed_leps, "leptons");

      // Missing momentum
      //declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning
      //book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
      //book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      //book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      
      // Histos from HEPdata at 200GeV
      book(_h["xsec_piplus_200"], 1, 1, 1);
      book(_h["xsec_piminus_200"], 1, 1, 2);

      book(_h["xsec_kplus_200"], 2, 1, 1);
      book(_h["xsec_kminus_200"], 2, 1, 2);

      book(_h["xsec_p_noFD_200_1"], 3, 1, 1);
      book(_h["xsec_p_noFD_200_2"], 9, 1, 1);
      book(_h["xsec_pbar_noFD_200_1"], 3, 1, 2);
      book(_h["xsec_pbar_noFD_200_2"], 9, 1, 2);

      book(_h["xsec_p_withFD_200_1"], 4, 1, 1);
      book(_h["xsec_p_withFD_200_2"], 10, 1, 1);
      book(_h["xsec_pbar_withFD_200_1"], 4, 1, 2);
      book(_h["xsec_pbar_withFD_200_2"], 10, 1, 2);

      // Histos from HEPdata at 62.4GeV
      book(_h["xsec_piplus_624"], 5, 1, 1);
      book(_h["xsec_piminus_624"], 5, 1, 2);

      book(_h["xsec_kplus_624"], 6, 1, 1);
      book(_h["xsec_kminus_624"], 6, 1, 2);

      book(_h["xsec_p_noFD_624_1"], 7, 1, 1);
      book(_h["xsec_pbar_noFD_624_1"], 7, 1, 2);

      book(_h["xsec_p_withFD_624_1"], 8, 1, 1);
      book(_h["xsec_pbar_withFD_624_1"], 8, 1, 2);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      cout << "Made it to analyze" << endl;

      Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();
      for( const Particle& p : fsParticles)
      {

		cout << "Made it into particle loop" << endl;
		
		// Histos 200GeV
      		if(p.pid() == 211) _h["xsec_piplus_200"]->fill(p.pT()/GeV, 1.0);
		if(p.pid() == -211) _h["xsec_piminus_200"]->fill(p.pT()/GeV, 1.0);

		if(p.pid() == 321) _h["xsec_kplus_200"]->fill(p.pT()/GeV, 1.0);
		if(p.pid() == -321) _h["xsec_kminus_200"]->fill(p.pT()/GeV, 1.0);

                cout << "Made it past pions and kaons" << endl;

		if(p.pid() == 2212) _h["xsec_p_noFD_200_1"]->fill(p.pT()/GeV, 1.0);
		if(p.pid() == 2212) _h["xsec_p_noFD_200_2"]->fill(p.pT()/GeV, 1.0);
		if(p.pid() == -2212) _h["xsec_pbar_noFD_200_1"]->fill(p.pT()/GeV, 1.0);
		if(p.pid() == -2212) _h["xsec_pbar_noFD_200_2"]->fill(p.pT()/GeV, 1.0);

		if(p.pid() == 2212) _h["xsec_p_withFD_200_1"]->fill(p.pT()/GeV, 1.0);
		if(p.pid() == 2212) _h["xsec_p_withFD_200_2"]->fill(p.pT()/GeV, 1.0);
		if(p.pid() == -2212) _h["xsec_pbar_withFD_200_1"]->fill(p.pT()/GeV, 1.0);
		if(p.pid() == -2212) _h["xsec_pbar_withFD_200_2"]->fill(p.pT()/GeV, 1.0);

                cout << "Made it past first four histos" << endl;

		// Histos 62.4GeV
                if(p.pid() == 211) _h["xsec_piplus_624"]->fill(p.pT()/GeV, 1.0);
                if(p.pid() == -211) _h["xsec_piminus_624"]->fill(p.pT()/GeV, 1.0);

                if(p.pid() == 321) _h["xsec_kplus_624"]->fill(p.pT()/GeV, 1.0);
                if(p.pid() == -321) _h["xsec_kminus_624"]->fill(p.pT()/GeV, 1.0);

                if(p.pid() == 2212) _h["xsec_p_noFD_624_1"]->fill(p.pT()/GeV, 1.0);
                if(p.pid() == -2212) _h["xsec_pbar_noFD_624_1"]->fill(p.pT()/GeV, 1.0);

                if(p.pid() == 2212) _h["xsec_p_withFD_624_1"]->fill(p.pT()/GeV, 1.0);
                if(p.pid() == -2212) _h["xsec_pbar_withFD_624_1"]->fill(p.pT()/GeV, 1.0);

                cout << "Made it past all histos" << endl;

      }

      // Retrieve dressed leptons, sorted by pT
      //vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      //Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

      // Remove all jets within dR < 0.2 of a dressed lepton
      //idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
      //Jets bjets = filter_select(jets, [](const Jet& jet) {
        //return  jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
      //});

      // Veto event if there are no b-jets
      //if (bjets.empty())  vetoEvent;

      // Apply a missing-momentum cut
      //if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent;

      // Fill histogram with leading b-jet pT
      //_h["XXXX"]->fill(bjets[0].pT()/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

                cout << "Made it into finalize" << endl;

      //normalize(_h["XXXX"]); // normalize to unity
      //normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
      //scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    //@}

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;
    //@}

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2011_I886590);

}
