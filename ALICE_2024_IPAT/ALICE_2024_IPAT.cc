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
      const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::pT > 0.15*GeV && Cuts::abscharge > 0);
      declare(aprim, "aprim");
      //The second primary particles is so that the jet specctra goes over all particles instead of cutting out pT < 0.15GeV
      //const ALICE::PrimaryParticles aprimall(Cuts::abseta < 0.9 && Cuts::abscharge > 0 && (Cuts::abspid == Rivet::PID::PIPLUS || Cuts::abspid == Rivet::PID::KPLUS || Cuts::abspid == Rivet::PID::PROTON || Cuts::abspid == Rivet::PID::ELECTRON || Cuts::abspid == Rivet::PID::MUON));
      //This is the jet constituents
      //Pat, change this so that it accepts charged primary particles OR photons
      //You also want to 
      const ALICE::PrimaryParticles aprimall(Cuts::abseta < 0.7 && Cuts::abscharge > 0);
      declare(aprimall, "aprimall");

      // jets à la ALICE - Jet area will be available using the pseudojet
      fastjet::AreaType fjAreaType = fastjet::active_area_explicit_ghosts;
      fastjet::GhostedAreaSpec fjGhostAreaSpec = fastjet::GhostedAreaSpec(1., 1, 0.005, 1., 0.1, 1e-100);
      fjAreaDef = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
      FastJets jetfs(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.4, fjAreaDef);
      declare(jetfs, "jetsfs");

      // Book histograms
      // specify custom binning
      book(_histos["dphi_pi"], "dphi_pi", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["dphi_p"], "dphi_p", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["dphi_k"], "dphi_k", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["dphi_ktopi"], "dphi_ktopi", 48 , - 3.14 / 2, 3 * 3.14 / 2);

      book(_histos["deta_pi"], "deta_pi", 60 , - 1.5, 1.5);
      book(_histos["deta_p"], "deta_p", 60 , - 1.5, 1.5);
      book(_histos["deta_k"], "deta_k", 60 , - 1.5, 1.5);
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      // book(_histos["AAAA"], 1, 1, 1);
      // book(_profilesrofiles["BBBB"], 2, 1, 1);
      // book(_counters["CCCC"], 3, 1, 1);
      // book(_counters["sow"], "sow"); // what in tarnation? This probably stands for sum of weights and I hate that a lot. Commented and renamed.
      book(_counters["number_of_events"], "number_of_events");
      book(_counters["number_of_jets"], "number_of_jets");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      std::cout << "Event number: " << event.genEvent()->event_number() << std::endl;
      _counters["number_of_events"]->fill();

      const FinalState fs = apply<FinalState>(event, "fs");
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
      // divide(_histos["dphi_k"], _histos["dphi_pi"], _scatters["dphi_ktopi"]);
      for (size_t i = 0; i < _histos["dphi_ktopi"]->numBins(); ++i)
      {
        const double binContent1 = _histos["dphi_k"]->bin(i).height();
        const double binContent2 = _histos["dphi_pi"]->bin(i).height();

        if (binContent2 != 0.0)
        {
          _histos["dphi_ktopi"]->bin(i).set(binContent1 / binContent2);
        }
        else
        {
          _histos["dhpi_ktopi"]->bin(i).set(0.0); // Handle division by zero
        }
      }
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
