// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"
//PAT: Add header for centrality. 

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2024_IPAT : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2024_IPAT);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
//PAT: You will need a line like this but it should specify 
  // declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      const FinalState fs(Cuts::abseta < 0.9 && Cuts::pT > 0.15*GeV && Cuts::abscharge > 0);
      declare(fs, "fs");
      // associated particles
      const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::pT > 0.15*GeV && Cuts::abscharge > 0);
      declare(aprim, "aprim");
      //The second primary particles is so that the jet specctra goes over all particles instead of cutting out pT < 0.15GeV
      //const ALICE::PrimaryParticles aprimall(Cuts::abseta < 0.9 && Cuts::abscharge > 0 && (Cuts::abspid == Rivet::PID::PIPLUS || Cuts::abspid == Rivet::PID::KPLUS || Cuts::abspid == Rivet::PID::PROTON || Cuts::abspid == Rivet::PID::ELECTRON || Cuts::abspid == Rivet::PID::MUON));


      //Pat, change this so that it accepts charged primary particles OR photons
      //You also want to 
      //This is the jet constituents
      const ALICE::PrimaryParticles aprimall(Cuts::abseta < 0.7 && Cuts::abscharge > 0);
      declare(aprimall, "aprimall");


      // jets à la ALICE - Jet area will be available using the pseudojet
      fastjet::AreaType fjAreaType = fastjet::active_area_explicit_ghosts;
      fastjet::GhostedAreaSpec fjGhostAreaSpec = fastjet::GhostedAreaSpec(1., 1, 0.005, 1., 0.1, 1e-100);
      fjAreaDef = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
      FastJets jetfs(aprimall, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.2, fjAreaDef);
      declare(jetfs, "jetsfs");

      // Book histograms

      //particle types: (pi+ + pi-), (p), (K+ + K-)
      //particle momentum: (1, 1.5), (1.5, 2), (2,3), (3,4), (4,5), (5,6), (6,10)
      // jet momentum: 20-40 GeV

      // Desired plots
      // Deta and Dphi for each particle type and particle momentum
      // Deta: 60 bins from -1.5 to 1.5
      // Dphi: 48 bins from -pi/2 to 3pi/2

      // Integrated yield for each particle type as a function of momentum
      // Yield: [(1, 1.5), (1.5, 2), (2,3), (3,4), (4,5), (5,6), (6,10)]
      // Ratio of K to pi and p to pi as a function of momentum
      // Ratio: [(1, 1.5), (1.5, 2), (2,3), (3,4), (4,5), (5,6), (6,10)]

      // specify custom binning
      book(_histos["dphi_pi"], "dphi_pi", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["dphi_p"], "dphi_p", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["dphi_k"], "dphi_k", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_scatters["dphi_ktopi"], "dphi_ktopi", 48 , - 3.14 / 2, 3 * 3.14 / 2);

      book(_histos["deta_pi"], "deta_pi", 60 , - 1.5, 1.5);
      book(_histos["deta_p"], "deta_p", 60 , - 1.5, 1.5);
      book(_histos["deta_k"], "deta_k", 60 , - 1.5, 1.5);


      book(_counters["number_of_events"], "number_of_events");
      book(_counters["number_of_jets"], "number_of_jets");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      std::cout << "Event number: " << event.genEvent()->event_number() << std::endl;
      _counters["number_of_events"]->fill();

      //Get final state particles, e.g. all particles "detected"
      const FinalState fs = apply<FinalState>(event, "fs");
      //Get jets
      FastJets jetsfs = apply<FastJets>(event, "jetsfs");
      //For spectra
      const ALICE::PrimaryParticles aprimall = apply<ALICE::PrimaryParticles>(event, "aprimall");
      const Particles ALICEparticlesall = aprimall.particles();
      //For fragmentation functions
      const ALICE::PrimaryParticles aprim = apply<ALICE::PrimaryParticles>(event, "aprim");
      const Particles ALICEparticles = aprim.particles();


      //jetsfs.calc(ALICEparticles);
      jetsfs.calc(ALICEparticlesall);

      Jets jets = jetsfs.jetsByPt(Cuts::abseta < 0.5 && Cuts::pT >= 5.*GeV);
      std::cout << "Number of jets: " << jets.size() << std::endl;
      std::cout << "Number of particles: " << ALICEparticles.size() << std::endl;
      for(auto jet : jets)
      {
        _counters["number_of_jets"]->fill();
        for(auto particle : ALICEparticles)
        {
            auto dphi = jet.phi()-particle.phi();
            auto deta = jet.eta()-particle.eta();
            // make sure dphi is in [-pi/2,3pi/2] range
            if (dphi > 3*3.14159/2) dphi -= 2*3.14159;
            if (dphi < -3*3.14159/2) dphi += 2*3.14159;
            if (particle.pid() == PID::PIPLUS)
            {
                _histos["dphi_pi"]->fill(dphi);
                _histos["deta_pi"]->fill(deta);
            }
            if (particle.pid() == PID::PIMINUS)
            {
                _histos["dphi_pi"]->fill(dphi);
                _histos["deta_pi"]->fill(deta);
            }
            if (particle.pid() == PID::PROTON)
            {
                _histos["dphi_p"]->fill(dphi);
                _histos["deta_p"]->fill(deta);
            }
            if (particle.pid() == PID::KPLUS)
            {
                _histos["dphi_k"]->fill(dphi);
                _histos["deta_k"]->fill(deta);
            }
            if (particle.pid() == PID::KMINUS)
            {
                _histos["dphi_k"]->fill(dphi);
                _histos["deta_k"]->fill(deta);
            }
        }
        
      }
      if(jets.size() == 0) vetoEvent;

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      float numJets = _counters["number_of_jets"]->sumW();
      _histos["dphi_pi"]->scaleW(1/(nanobarn*numJets));
      _histos["deta_pi"]->scaleW(1/(nanobarn*numJets));
      _histos["dphi_p"]->scaleW(1/(nanobarn*numJets));
      _histos["deta_p"]->scaleW(1/(nanobarn*numJets));
      _histos["dphi_k"]->scaleW(1/(nanobarn*numJets));
      _histos["deta_k"]->scaleW(1/(nanobarn*numJets));
      divide(_histos["dphi_k"], _histos["dphi_pi"], _scatters["dphi_ktopi"]);
      // a change
      
      // normalize(_histos["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in pb (no cuts)
      // scale(_histos["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _histos;
    map<string, Profile1DPtr> _profiles;
    map<string, CounterPtr> _counters;
    map<string, Scatter2DPtr> _scatters;
    /// @}
    fastjet::AreaDefinition *fjAreaDef;


  };


  RIVET_DECLARE_PLUGIN(ALICE_2024_IPAT);

}
