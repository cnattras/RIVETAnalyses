// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2019_I1693308 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2019_I1693308);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState fs(Cuts::abseta < 0.9 && Cuts::pT > 0.15*GeV && Cuts::abscharge > 0);
      declare(fs, "fs");
      const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::pT > 0.15*GeV && Cuts::abscharge > 0);
      declare(aprim, "aprim");
      //The second primary particles is so that the jet specctra goes over all particles instead of cutting out pT < 0.15GeV
      const ALICE::PrimaryParticles aprimall(Cuts::abseta < 0.9 && Cuts::abscharge > 0 && (Cuts::abspid == Rivet::PID::PIPLUS || Cuts::abspid == Rivet::PID::KPLUS || Cuts::abspid == Rivet::PID::PROTON || Cuts::abspid == Rivet::PID::ELECTRON || Cuts::abspid == Rivet::PID::MUON));
      declare(aprimall, "aprimall");

      // jets à la ALICE - Jet area will be available using the pseudojet
      fastjet::AreaType fjAreaType = fastjet::active_area_explicit_ghosts;
      fastjet::GhostedAreaSpec fjGhostAreaSpec = fastjet::GhostedAreaSpec(1., 1, 0.005, 1., 0.1, 1e-100);
      fjAreaDef = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
      FastJets jetfs(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.4, fjAreaDef, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jetsfs");

      book(_h["PPS7"],1,1,1);
      book(_h["pTJ5-10"],2,1,1);
      book(_h["pTJ10-15"],3,1,1);
      book(_h["pTJ15-20"],4,1,1);

      book(_c["sow"], "sow");
      book(_c["sow2"],"sow2");
      book(_c["sow3"],"sow3");
      book(_c["sow4"],"sow4");

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

      if(jets.size() == 0) vetoEvent;

      double UEphiPlus = mapAngle0To2Pi(jets[0].phi() + M_PI/2.);
      double UEphiMinus = mapAngle0To2Pi(jets[0].phi() - M_PI/2.);
      double UEeta = jets[0].eta();

      double UEpT = 0.;

      for(const Particle& p : ALICEparticlesall)
      {
              double DeltaRPlus = sqrt(pow(p.eta()-UEeta, 2) + pow(p.phi()-UEphiPlus, 2));
              double DeltaRMinus = sqrt(pow(p.eta()-UEeta, 2) + pow(p.phi()-UEphiMinus, 2));
              if(DeltaRPlus > 0.4 && DeltaRMinus > 0.4) continue;

              UEpT += p.pT()/GeV;
      }

      //Dividing by 2 because two cones were used and by cone area piR^2
      //pT-density
      UEpT /= 2.*M_PI*0.16;


      for(auto jet : jets)
      {
        _h["PPS7"]->fill(jet.pT()/GeV - (UEpT*jet.pseudojet().area())); //pT-density*JetArea
      }

      jetsfs.calc(ALICEparticles);
      Jets jets2 = jetsfs.jetsByPt(Cuts::abseta < 0.5 && Cuts::pT >= 5.*GeV);
      for(auto jet : jets2)
      {

        if(jet.pT() >= 5*GeV && jet.pT() < 10*GeV){

          _c["sow2"]->fill();

        }else if(jet.pT() >= 10*GeV && jet.pT() < 15*GeV){

          _c["sow3"]->fill();

        } else if(jet.pT() >= 15*GeV && jet.pT() < 20*GeV){

          _c["sow4"]->fill();

        }

        for(const Particle&p:jet.particles()){
        if(jet.pT() >= 5*GeV && jet.pT() < 10*GeV){

          _h["pTJ5-10"]->fill((p.pT()/GeV)/(jet.pT()/GeV));

        }else if(jet.pT() >= 10*GeV && jet.pT() < 15*GeV){

          _h["pTJ10-15"]->fill((p.pT()/GeV)/(jet.pT()/GeV));

        } else if(jet.pT() >= 15*GeV && jet.pT() < 20*GeV){

          _h["pTJ15-20"]->fill((p.pT()/GeV)/(jet.pT()/GeV));

        }
      }

    }
  }


    /// Normalise histograms etc., after the run
    void finalize() {
      //double norm = (crossSection()*1.E-9)/sow->sumW();
      //Above has "sow" which isn't anything that is declared so I changed it to the ptr _c["sow"]
      //double norm = (crossSection()*1.E-9)/(_c["sow"]->sumW());
      double norm = (crossSection()/millibarn)/(_c["sow"]->sumW());


      _h["PPS7"]->scaleW(norm);

      _h["pTJ5-10"]->scaleW(1./_c["sow"]->sumW()*(_c["sow"]->effNumEntries()/_c["sow2"]->effNumEntries()));
      _h["pTJ10-15"]->scaleW(1./_c["sow"]->sumW()*(_c["sow"]->effNumEntries()/_c["sow3"]->effNumEntries()));
      _h["pTJ15-20"]->scaleW(1./_c["sow"]->sumW()*(_c["sow"]->effNumEntries()/_c["sow4"]->effNumEntries()));

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    fastjet::AreaDefinition *fjAreaDef;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(ALICE_2019_I1693308);

}
