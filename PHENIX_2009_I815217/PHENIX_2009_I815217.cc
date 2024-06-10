// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
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
#include "../Centralities/RHICCentrality.hh"
#include <math.h>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES

namespace Rivet {

  /// @brief Add a short analysis description here
  class PHENIX_2009_I815217 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2009_I815217);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      


      // Initialise and register projections
      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 0.35);
      declare(fs, "fs")

      const ParticlePair& beam = beams();

      
      beamOpt = getOption<string>("beam","NONE");
      if (beamOpt == "NONE") {
      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu200;
      }

      // Book histograms
      book(_s["JPsiyield"], 1, 1, 1);
      book(_s["EEyield2023"], 2, 1, 1);
      book(_s["EEyield2328"], 2, 1, 2);
      book(_s["EEyield2028"], 2, 1, 3);
      book(_s["JPsicross"], 3, 1, 1);
      book(_s["EEcross2023"], 4, 1, 1);
      book(_s["EEcross2328"], 4, 1, 2);
      book(_s["EEcross2028"], 4, 1, 3);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
    
      // Fill histogram with leading b-jet pT
     // _h["XXXX"]->fill(bjets[0].pT()/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      
     
   

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;

    string beamOpt;
    enum CollisionSystem {AuAu200};
    CollisionSystem collSys;
    //@}


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2009_I815217);

}
