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
  class PHENIX_2012_I1116179 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2012_I1116179);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

	  declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");


	  //const FinalState fs(Cuts::abseta < 0.5 && Cuts::pT > 0.15*GeV);
	  //declare(fs, "fs");

	  const PromptFinalState pfs(Cuts::abseta < 0.5 && Cuts::pT > 4.0 && Cuts::pT < 22.0);
	  declare(pfs, "pfs");

//	  const UnstableParticles up(Cuts::abseta < 0.5 && Cuts::abspid == 421);
//	  declare(up, "up");

	  book(_h["dir_photon_AuAu0005"], 1, 1, 1);
	  book(_c["sow_AuAu0005"], "sow_AuAu0005");
//	  book(_c["sow_AuAu3050"], sow_AuAu3050);
//	  book(_c["sow_pp"], "sow_pp");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();

      if(c < 5.) 
      {
        _c["sow_AuAu0005"]->fill();
      }

      Particles pfsParticles = applyProjection<FinalState>(event,"pfs").particles();

      for(const Particle& p : pfsParticles) 
      {
        if(c < 5. && p.pid() == 22) _h["dir_photon_AuAu0005"]->fill(p.pT()/GeV);
        
      }






    }


    /// Normalise histograms etc., after the run
    void finalize() {
      _h["dir_photon_AuAu0005"]->scaleW(1./_c["sow_AuAu0005"]->sumW());


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


  DECLARE_RIVET_PLUGIN(PHENIX_2012_I1116179);

}
