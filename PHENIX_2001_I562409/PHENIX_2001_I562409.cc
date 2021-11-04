// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/HadronicFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2001_I562409 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2001_I562409);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const HadronicFinalState hfs(Cuts::abseta < 0.5 && Cuts::abscharge > 0);
      declare(hfs, "hfs");

      const UnstableParticles pi0UP(Cuts::abseta < 0.5 && Cuts::abspid == 111);
      declare(pi0UP,"pi0UP");

      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["ChHadronsCent0_10"], 5, 1, 1);
      book(_h["ChHadronsCent60_80"], 6, 1, 1);
      book(_c["Cent0_10"], "Cent0_10");
      book(_c["Cent60_80"], "Cent60_80");
      book(_h["Pi0PbScCent0_10"], 1, 1, 1);
      book(_h["Pi0PbScCent60_80"], 2, 1, 1);
      book(_h["Pi0PbGlCent0_10"], 3, 1, 1);
      book(_h["Pi0PbGlCent60_80"], 4, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
      const double cent = centProj();
      
      const HadronicFinalState hfs = apply<HadronicFinalState>(event, "hfs");
      const Particles hfsParticles = hfs.particles();

      const UnstableParticles pi0UP = apply<UnstableParticles>(event, "pi0UP");
      const Particles pi0UPParticles = pi0UP.particles();

      double inv2PI = 1./(2.*M_PI);

            if(cent < 10.) //Check centrality of the event
      {
              _c["Cent0_10"]->fill(); //fill counter for 0-10% most central events
              for(auto p : hfsParticles) //loop over charged hadrons
              {
                      _h["ChHadronsCent0_10"]->fill(p.pT()/GeV, inv2PI/(2.*p.pT()/GeV)); //additional 1/2 factor to take into account h^(+)+h^(-)/2
              }
              
              for(auto p : pi0UPParticles) //loop over pi0s
              {
                      _h["Pi0PbScCent0_10"]->fill(p.pT()/GeV, inv2PI/(p.pT()/GeV));
		      _h["Pi0PbGlCent0_10"]->fill(p.pT()/GeV, inv2PI/(p.pT()/GeV));
              }
      }
      else if(cent >= 60. && cent < 80.) //Check centrality of the event
      {
              _c["Cent60_80"]->fill(); //fill counter for 60-80% most central events
              for(auto p : hfsParticles) //loop over charged hadrons
              {
                      _h["ChHadronsCent60_80"]->fill(p.pT()/GeV, inv2PI/(2.*p.pT()/GeV)); //additional 1/2 factor to take into account h^(+)+h^(-)/2
              }

              for(auto p : pi0UPParticles) //loop over pi0s
              {
                      _h["Pi0PbScCent60_80"]->fill(p.pT()/GeV, inv2PI/(p.pT()/GeV));
		      _h["Pi0PbGlCent60_80"]->fill(p.pT()/GeV, inv2PI/(p.pT()/GeV)); 
              }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //Divide histograms by number of events
      _h["ChHadronsCent0_10"]->scaleW(1./_c["Cent0_10"]->sumW());
      _h["ChHadronsCent60_80"]->scaleW(1./_c["Cent60_80"]->sumW());
      _h["Pi0PbScCent0_10"]->scaleW(1./_c["Cent0_10"]->sumW());
      _h["Pi0PbScCent60_80"]->scaleW(1./_c["Cent60_80"]->sumW());
      _h["Pi0PbGlCent0_10"]->scaleW(1./_c["Cent0_10"]->sumW());
      _h["Pi0PbGlCent60_80"]->scaleW(1./_c["Cent0_10"]->sumW());
    }

      ///@}
      //@name Histograms
      ///@{
      map<string, Histo1DPtr> _h;
      map<string, Profile1DPtr> _p;
      map<string, CounterPtr> _c;
      ///@}
  };


  RIVET_DECLARE_PLUGIN(PHENIX_2001_I562409);

}
