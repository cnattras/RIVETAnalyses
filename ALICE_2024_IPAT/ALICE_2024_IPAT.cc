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
      book(_histos["dphi"], "dphi", -3.14/2, 3*3.14/2, 48);
      book(_histos["deta"], "deta", -1.5, 1.5, 60);
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

      for(auto jet : jets)
      {
        _counters['number_of_jets'].fill();
        for(auto particle : ALICEparticles)
        {
            auto dphi = jet.phi()-particle.phi();
            auto deta = jet.eta()-particle.eta();
            _histos["dphi"]->fill(dphi);
            _histos["deta"]->fill(deta);
        }
      }
      if(jets.size() == 0) vetoEvent;

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_histos["dphi"], 1/_counters["number_of_jets"]); // normalize by number of jets
      scale(_histos["eta"], 1/_counters["number_of_jets"]); // normalize by number of jets
      // normalize(_histos["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in pb (no cuts)
      // scale(_histos["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _histos;
    map<string, Profile1DPtr> _profiles;
    map<string, CounterPtr> _counters;
    /// @}
    fastjet::AreaDefinition *fjAreaDef;


  };


  RIVET_DECLARE_PLUGIN(ALICE_2024_IPAT);

}
