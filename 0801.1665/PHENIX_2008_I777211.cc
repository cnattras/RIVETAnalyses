// -*- C++ -*-

//need to remove unused packages
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <math.h>
#define _USE_MATH_DEFINES


namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2008_I777211 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2008_I777211);



    void init() {


      //needs centrality for PHENIX
      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M","V0M");

      //not sure about abscharge. does this relate to the neutral pions?
      const FinalState fs(Cuts::abseta < 0.35 && Cuts::abscharge == 0);
      declare(fs, "fs");

      book(hPion0Pt, 1, 1, 1);
      book(sow,"sow");

    }


    void analyze(const Event& event) {

      const CentralityProjection& centProj = apply<CentralityProjection>(event,"V0M");
      const double cent = centProj();

      //centrality for 0-5%
      if (cent > 5.) vetoEvent;
      sow->fill();

      //renamed for neutral particle (it was chargedParticlesY05. what is the Y05?)
      Particles neutralParticlesY05 = applyProjection<FinalState>(event,"fs").particles();

      for(const Particle& p : neutralParticlesY05)
      {
        if(p.pid() == 111)
        {
    	     double partPt = p.pT()/GeV;
           //need to talk about the paper's weight. not sure if it is the same as here.
           double pt_weight = 1./(partPt*2.*M_PI);
    	     hPion0Pt->fill(partPt, pt_weight);
         }
       }
    }



    void finalize() {
      //need to talk about this too.
      hPion0Pt->scaleW(1./sow->sumW());

    }


    Histo1DPtr hPion0Pt;
    CounterPtr sow;


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2008_I777211);

}
