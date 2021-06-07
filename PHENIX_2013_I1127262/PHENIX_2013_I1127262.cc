// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2013_I1127262 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2013_I1127262);

    //adding a function to get the acc of the event plant dectector
    //so that reaction plant dependency can be calculated
    double GetEventPlaneDetectorAcc(int n, const FinalState pf, vector<Cut> etaRxP, int nPhiSections)
    {
            double QIn = 0.;
            double QRn = 0.;

            for(auto eta : etaRxP)
            {
                    for(int iphi = 0; iphi < nPhiSections; iphi++)
                    {
                            Cut aCut = eta && Cuts::phi > iphi*(2.*M_PI/nPhiSections) && Cuts::phi < (iphi+1)*(2.*M_PI/nPhiSections);
                            Particles particles = pf.particles(aCut);
                            int weight = particles.size();
                            for(const Particle& p : particles)
                            {
                                    QIn += weight*sin(n*p.phi());
                                    QRn += weight*cos(n*p.phi());
                            }
                    }
            }

            double eventPlane = mapAngle0To2Pi((1./n)*atan2(QIn,QRn));

            return eventPlane;
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

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
      book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
      book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["AAAA"], 1, 1, 1);
      book(_p["BBBB"], 2, 1, 1);
      book(_c["CCCC"], 3, 1, 1);

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

      //Calcualte Reaction Plane Dependency:
      
      //get reaction plane positive and negative final state values
      const FinalState& RxP = apply<FinalState>(event, "RxP");
      const FinalState& RxPPos = apply<FinalState>(event, "RxPPos");
      const FinalState& RxPNeg = apply<FinalState>(event, "RxPNeg");

      //Inner and Outer rings of North and South sections of the RxP dectector
      vector<Cut> etaRxP = {Cuts::eta > 1. && Cuts::eta < 1.5, Cuts::eta > 1.5 && Cuts::eta < 2.8, Cuts::eta < -1. && Cuts::eta > -1.5, Cuts::eta < -1.5 && Cuts::eta > -2.8};
      int nPhiSections = 12;

      double evPPosNeg = GetEventPlaneDetectorAcc(2, RxP, etaRxP, nPhiSections);

      vector<Cut> etaRxPPos = {Cuts::eta > 1. && Cuts::eta < 1.5, Cuts::eta > 1.5 && Cuts::eta < 2.8};
      vector<Cut> etaRxPNeg = {Cuts::eta < -1. && Cuts::eta > -1.5, Cuts::eta < -1.5 && Cuts::eta > -2.8};

      double evPPos = GetEventPlaneDetectorAcc(2, RxPPos, etaRxPPos, nPhiSections);
      double evPNeg = GetEventPlaneDetectorAcc(2, RxPNeg, etaRxPNeg, nPhiSections);

      _p["RxPcosPos"]->fill(int(floor(c/10))+0.5, cos(2*(evPPos-evPNeg)));



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
    std::vector<double> v2centBins = {0., 10., 20., 30., 40., 50., 60.};
    //@}


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2013_I1127262);

}
