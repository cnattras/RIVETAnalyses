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
      FastJets jetfs(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.4, fjAreaDef, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jetsfs");

      // Book histograms
      // specify custom binning
      book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
      book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      // book(_h["AAAA"], 1, 1, 1);
      // book(_p["BBBB"], 2, 1, 1);
      // book(_c["CCCC"], 3, 1, 1);
      book(_c["sow"], "sow");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      _c["sow"]->fill();

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
            _h["XXXX"]->fill(jet.pT());
            _h["YYYY"]->fill(10.0*jet.pT());
            _h["ZZZZ"]->fill(jet.pT());
      }
      if(jets.size() == 0) vetoEvent;

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // normalize(_h["XXXX"]); // normalize to unity
      // normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in pb (no cuts)
      // scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    /// @}
    fastjet::AreaDefinition *fjAreaDef;


  };


  RIVET_DECLARE_PLUGIN(ALICE_2024_IPAT);

}
