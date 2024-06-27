// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2001_I555603 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2001_I555603);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Implement beam options
      beamOpt = getOption<string>("beam", "NONE");

      const ParticlePair& beam = beams();

      if (beamOpt == "NONE"){
      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu130;
      }

      if (beamOpt == "AUAU130") collSys = AuAu130;

      // Declare Centrality
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      // ALICE Projection
      declare(ALICE::PrimaryParticles(Cuts::abseta < 0.38 && Cuts::pT > 0.0*MeV && Cuts::abscharge > 0), "APRIM");

      // Book histograms
      book(_h["dEtdEta"], 1, 1, 1);
      book(_hist_E, "d01-x01-y01-test", refData(1, 1, 1));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

        double totalEt = 0;
        double deltaeta = .76; //does this relate to the abseta in the ALICE declaration?

        Particles chargedParticles = applyProjection<ALICE::PrimaryParticles>(event,"APRIM").particles();

        const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
        const double c = cent();
        if (c > 50) vetoEvent;

        for(const Particle& p : chargedParticles)
        {
            totalEt += p.Et()/GeV;
        }


        if (collSys == AuAu130){
          _hist_E->fill(c,totalEt/deltaeta);
        }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    Profile1DPtr _hist_E;
    string beamOpt;
    enum CollisionSystem {AuAu130};
    CollisionSystem collSys;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PHENIX_2001_I555603);

}
