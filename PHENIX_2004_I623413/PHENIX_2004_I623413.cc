// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "../Centralities/RHICCentrality.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2004_I623413 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2004_I623413);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 0.35);
      declare(fs, "fs");
      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
      book(_h["KaonPlusMinBias"], 1, 1, 1);
      book(_h["KaonMinusMinBias"], 1, 1, 2);
      book(_h["PionPlusMinBias"], 2, 1, 1);
      book(_h["PionMinusMinBias"], 2, 1, 2);
      book(_h["ProtonMinBias"], 3, 1, 1);
      book(_h["AntiProtonMinBias"], 3, 1, 2);
      book(_h["PionPlusCent0_5"], 4, 1, 1);
      book(_h["PionPlusCent5_15"], 4, 1, 2);
      book(_h["PionPlusCent15_30"], 4, 1, 3);
      book(_h["PionPlusCent30_60"], 4, 1, 4);
      book(_h["PionMinusCent0_5"], 4, 1, 5);
      book(_h["PionMinusCent5_15"], 4, 1, 6);
      book(_h["PionMinusCent15_30"], 4, 1, 7);
      book(_h["PionMinusCent30_60"], 4, 1, 8);
      book(_h["PionPlusCent60_92"], 5, 1, 1);
      book(_h["PionMinusCent60_92"], 5, 1, 2);
      book(_h["KaonPlusCent0_5"], 6, 1, 1);
      book(_h["KaonPlusCent5_15"], 6, 1, 2);
      book(_h["KaonPlusCent15_30"], 6, 1, 3);
      book(_h["KaonPlusCent30_60"], 6, 1, 4);
      book(_h["KaonMinusCent0_5"], 6, 1, 5);
      book(_h["KaonMinusCent5_15"], 6, 1, 6);
      book(_h["KaonMinusCent15_30"], 6, 1, 7);
      book(_h["KaonMinusCent30_60"], 6, 1, 8);
      book(_h["KaonPlusCent60_92"], 7, 1, 1);
      book(_h["KaonMinusCent60_92"], 7, 1, 2);
      book(_h["ProtonCent0_5"], 8, 1, 1);
      book(_h["ProtonCent0_5"], 8, 1, 2);
      book(_h["ProtonCent5_15"], 9, 1, 1);
      book(_h["AntiProtonCent5_15"], 9, 1, 2);
      book(_h["ProtonCent15_30"], 10, 1, 1);
      book(_h["AntiProtonCent15_30"], 11, 1, 1);
      book(_h["ProtonCent30_60"], 12, 1, 1);
      book(_h["AntiProtonCent30_60"], 13, 1, 1);
      book(_h["ProtonCent60_92"], 14, 1, 1);
      book(_h["AntiProtonCent60_92"], 14, 1, 2);
      book(_h["PionPlusMinBias"], 15, 1, 1);
      book(_h["PionMinusMinBias"], 15, 1, 2);
      book(_h["KaonPlusMinBias"], 16, 1, 1);
      book(_h["KaonMinusMinBias"], 16, 1, 2);
      book(_h["ProtonMinBias"], 17, 1, 1);
      book(_h["AntiProtonMinBias"], 17, 1, 2);
      book(_h["PionPlusCent0_5"], 18, 1, 1);
      book(_h["PionPlusCent5_15"], 18, 1, 2);
      book(_h["PionPlusCent15_30"], 18, 1, 3);
      book(_h["PionPlusCent30_60"], 18, 1, 4);
      book(_h["PionMinusCent0_15"], 18, 1, 5);
      book(_h["PionMinusCent5_15"], 18, 1, 6);
      book(_h["PionMinusCent15_30"], 18, 1, 7);
      book(_h["PionMinusCent30_60"], 18, 1, 8);
      book(_h["PionPlusCent60_92"], 19, 1, 1);
      book(_h["PionMinusCent60_92"], 19, 1, 2);
      book(_h["KaonPlusCent0_5"], 20, 1, 1);
      book(_h["KaonPlusCent5_15"], 20, 1, 2);
      book(_h["KaonPlusCent15_30"], 20, 1, 3);
      book(_h["KaonPlusCent30_60"], 20, 1, 4);
      book(_h["KaonPlusCent60_92"], 20, 1, 5);
      book(_h["KaonMinusCent0_5"], 20, 1, 6);
      book(_h["KaonMinusCent5_15"], 20, 1, 7);
      book(_h["KaonMinusCent15_30"], 20, 1, 8);
      book(_h["KaonMinusCent30_60"], 20, 1, 9);
      book(_h["KaonMinusCent60_92"], 20, 1, 10);
      book(_h["ProtonCent0_5"], 21, 1, 1);
      book(_h["ProtonCent5_15"], 21, 1, 2);
      book(_h["ProtonCent15_30"], 21, 1, 3);
      book(_h["ProtonCent30_60"], 21, 1, 4);
      book(_h["AntiProtonCent0_5"], 21, 1, 5);
      book(_h["AntiProtonCent5_15"], 21, 1, 6);
      book(_h["AntiProtonCent15_30"], 21, 1, 7);
      book(_h["AntiProtonCent30_60"], 21, 1, 8);
      book(_h["ProtonCent60_92"], 22, 1, 1);
      book(_h["AntiProtonCent60_92"], 22, 1, 2);
      book(_h["Table7FitParameters_Apart"], 23, 1, 1);
      book(_h["Table7FitParameters_Acoll"], 23, 1, 2);
      book(_h["Table8FitParameters_Apart"], 24, 1, 1);
      book(_h["Table8FitParameters_Acoll"], 24, 1, 2);
      book(_h["PionPlusParticle"], 25, 1, 1);
      book(_h["PionMinusParticle"], 25, 1, 2);
      book(_h["KaonPlusParticle"], 25, 1, 3);
      book(_h["KaonMinusParticle"], 25, 1, 4);
      book(_h["ProtonParticle"], 25, 1, 5);
      book(_h["AntiProtonParticle"], 25, 1, 6);
      book(_h["N_collCentrality"], 26, 1, 1);
      book(_h["N_partCentrality"], 27, 1, 1);
      book(_h["PionPlusCentrality"], 28, 1, 1);
      book(_h["PionMinusCentrality"], 28, 1, 2);
      book(_h["KaonPlusCentrality"], 28, 1, 3);
      book(_h["KaonMinusCentrality"], 28, 1, 4);
      book(_h["ProtonCentrality"], 28, 1, 5);
      book(_h["AntiProtonCentrality"], 28, 1, 6);
      book(_h["PionPlusCentrality"], 29, 1, 1);
      book(_h["PionMinusCentrality"], 29, 1, 2);
      book(_h["KaonPlusCentrality"], 29, 1, 3);
      book(_h["KaonMinusCentrality"], 29, 1, 4);
      book(_h["ProtonCentrality"], 29, 1, 5);
      book(_h["AntiProtonCentrality"], 29, 1, 6);
      book(_h["BestFitCentrality"], 30, 1, 1);
      book(_h["BestFitCentrality"], 31, 1, 1);
      book(_h["BestFitCentrality"], 31, 1, 2);
      book(_h["BestFitCentrality"], 31, 1, 3);
      book(_h["BestFitCentrality"], 32, 1, 1);
      book(_h["BestFitCentrality"], 32, 1, 2);
      book(_h["BestFitCentrality"], 32, 1, 3);
      book(_c["MinBias"], "MinBias");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
      const double cent = centProj();
      Particles particles = applyProjection<FinalState>(event,"fs").particles();
      if(cent > 92.) vetoEvent;
        _c["MinBias"]->fill();



      for(Particle& p : particles){
        if(p.pid() == 321)
        {
          _h["KaonPlusMinBias"]->fill(p.pT()/GeV);

        }
        if(p.pid() == -321)
        {
          _h["KaonMinusMinBias"]->fill(p.pT()/GeV);

        }
        if(p.pid() == 2212)
        {
          _h["ProtonMinBias"]->fill(p.pT()/GeV);

        }
        if(p.pid() == -2212)
        {
          _h["AntiProtonMinBias"]->fill(p.pT()/GeV);

        }
        if(p.pid() == 211)
        {
          _h["PionPlusMinBias"]->fill(p.pT()/GeV);

        }
        if(p.pid() == -211)
        {
          _h["PionMinusMinBias"]->fill(p.pT()/GeV);

        }
       }


    }


    /// Normalise histograms etc., after the run
    void finalize() {


    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(PHENIX_2004_I623413);

}
