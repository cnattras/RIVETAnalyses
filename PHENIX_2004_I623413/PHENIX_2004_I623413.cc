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
      book(_h["Fig10PionPlusCent0_5"], 4, 1, 1);
      book(_h["FIg10PionPlusCent5_15"], 4, 1, 2);
      book(_h["Fig10PionPlusCent15_30"], 4, 1, 3);
      book(_h["Fig10PionPlusCent30_60"], 4, 1, 4);
      book(_h["Fig10PionMinusCent0_5"], 4, 1, 5);
      book(_h["Fig10PionMinusCent5_15"], 4, 1, 6);
      book(_h["Fig10PionMinusCent15_30"], 4, 1, 7);
      book(_h["Fig10PionMinusCent30_60"], 4, 1, 8);
      book(_h["Fig10PionPlusCent60_92"], 5, 1, 1);
      book(_h["Fig10PionMinusCent60_92"], 5, 1, 2);
      book(_h["Fig10KaonPlusCent0_5"], 6, 1, 1);
      book(_h["Fig10KaonPlusCent5_15"], 6, 1, 2);
      book(_h["Fig10KaonPlusCent15_30"], 6, 1, 3);
      book(_h["Fig10KaonPlusCent30_60"], 6, 1, 4);
      book(_h["Fig10KaonMinusCent0_5"], 6, 1, 5);
      book(_h["Fig10KaonMinusCent5_15"], 6, 1, 6);
      book(_h["Fig10KaonMinusCent15_30"], 6, 1, 7);
      book(_h["Fig10KaonMinusCent30_60"], 6, 1, 8);
      book(_h["Fig10KaonPlusCent60_92"], 7, 1, 1);
      book(_h["Fig10KaonMinusCent60_92"], 7, 1, 2);
      book(_h["Fig10ProtonCent0_5"], 8, 1, 1);
      book(_h["Fig10AntiProtonCent0_5"], 8, 1, 2);
      book(_h["Fig10ProtonCent5_15"], 9, 1, 1);
      book(_h["Fig10AntiProtonCent5_15"], 9, 1, 2);
      book(_h["Fig10ProtonCent15_30"], 10, 1, 1);
      book(_h["Fig10AntiProtonCent15_30"], 11, 1, 1);
      book(_h["Fig10ProtonCent30_60"], 12, 1, 1);
      book(_h["Fig10AntiProtonCent30_60"], 13, 1, 1);
      book(_h["Fig10ProtonCent60_92"], 14, 1, 1);
      book(_h["Fig10AntiProtonCent60_92"], 14, 1, 2);
      book(_h["Fig11PionPlusMinBias"], 15, 1, 1);
      book(_h["Fig11PionMinusMinBias"], 15, 1, 2);
      book(_h["Fig11KaonPlusMinBias"], 16, 1, 1);
      book(_h["Fig11KaonMinusMinBias"], 16, 1, 2);
      book(_h["Fig11ProtonMinBias"], 17, 1, 1);
      book(_h["Fig11AntiProtonMinBias"], 17, 1, 2);
      book(_h["Fig12PionPlusCent0_5"], 18, 1, 1);
      book(_h["Fig12PionPlusCent5_15"], 18, 1, 2);
      book(_h["Fig12PionPlusCent15_30"], 18, 1, 3);
      book(_h["Fig12PionPlusCent30_60"], 18, 1, 4);
      book(_h["Fig12PionMinusCent0_5"], 18, 1, 5);
      book(_h["Fig12PionMinusCent5_15"], 18, 1, 6);
      book(_h["Fig12PionMinusCent15_30"], 18, 1, 7);
      book(_h["Fig12PionMinusCent30_60"], 18, 1, 8);
      book(_h["Fig12PionPlusCent60_92"], 19, 1, 1);
      book(_h["Fig12PionMinusCent60_92"], 19, 1, 2);
      book(_h["Fig12KaonPlusCent0_5"], 20, 1, 1);
      book(_h["Fig12KaonPlusCent5_15"], 20, 1, 2);
      book(_h["Fig12KaonPlusCent15_30"], 20, 1, 3);
      book(_h["Fig12KaonPlusCent30_60"], 20, 1, 4);
      book(_h["Fig12KaonPlusCent60_92"], 20, 1, 5);
      book(_h["Fig12KaonMinusCent0_5"], 20, 1, 6);
      book(_h["Fig12KaonMinusCent5_15"], 20, 1, 7);
      book(_h["Fig12KaonMinusCent15_30"], 20, 1, 8);
      book(_h["Fig12KaonMinusCent30_60"], 20, 1, 9);
      book(_h["Fig12KaonMinusCent60_92"], 20, 1, 10);
      book(_h["Fig12ProtonCent0_5"], 21, 1, 1);
      book(_h["Fig12ProtonCent5_15"], 21, 1, 2);
      book(_h["Fig12ProtonCent15_30"], 21, 1, 3);
      book(_h["Fig12ProtonCent30_60"], 21, 1, 4);
      book(_h["Fig12AntiProtonCent0_5"], 21, 1, 5);
      book(_h["Fig12AntiProtonCent5_15"], 21, 1, 6);
      book(_h["Fig12AntiProtonCent15_30"], 21, 1, 7);
      book(_h["Fig12AntiProtonCent30_60"], 21, 1, 8);
      book(_h["Fig12ProtonCent60_92"], 22, 1, 1);
      book(_h["Fig12AntiProtonCent60_92"], 22, 1, 2);
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
      book(_h["Table4PionPlusCentrality"], 29, 1, 1);
      book(_h["Table4PionMinusCentrality"], 29, 1, 2);
      book(_h["Table4KaonPlusCentrality"], 29, 1, 3);
      book(_h["Table4KaonMinusCentrality"], 29, 1, 4);
      book(_h["Table4ProtonCentrality"], 29, 1, 5);
      book(_h["Table4AntiProtonCentrality"], 29, 1, 6);
      book(_c["MinBias"], "MinBias");
      book(_c["Cent0_5"], "Cent0_5");
      book(_c["Cent5_15"], "Cent5_15");
      book(_c["Cent15_30"], "Cent15_30");
      book(_c["Cent30_60"], "Cent30_60");
      book(_c["Cent60_92"], "Cent60_92");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
      const double cent = centProj();
      Particles particles = applyProjection<FinalState>(event,"fs").particles();
      if(cent > 92.) vetoEvent;
        _c["MinBias"]->fill();
        if(cent >= 0. && cent < 5.) _c["Cent0_5"]->fill();
        else if(cent >= 5. && cent < 15.) _c["Cent5_15"]->fill();
        else if(cent >= 15. && cent < 30.) _c["Cent15_30"]->fill();
        else if(cent >= 30. && cent < 60.) _c["Cent30_60"]->fill();
        else if(cent >= 60. && cent < 92.) _c["Cent60_92"]->fill();

      for(Particle& p : particles){
        if(p.pid() == 321)
        {
          _h["KaonPlusMinBias"]->fill(p.pT()/GeV);
          _h["Fig11KaonPlusMinBias"]->fill(p.pT()/GeV);

          if(cent >= 0. && cent < 5.) {
            _h["KaonPlusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig10KaonPlusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig11KaonPlusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12KaonPlusCent0_5"]->fill(p.pT()/GeV);
            _h["KaonPlusParticle0_5"]->fill(p.pT()/GeV);
            _h["KaonPlusCentrality0_5"]->fill(p.pT()/GeV);
            _h["Table4KaonPlusCentrality0_5"]->fill(p.pT()/GeV);
            }
          else if(cent >= 5. && cent < 15.)  {
            _h["KaonPlusCent5_15"]->fill(p.pT()/GeV);
            _h["Fig10KaonPlusCent5_15"]->fill(p.pT()/GeV);
            _h["KaonPlusParticle5_15"]->fill(p.pT()/GeV);
            _h["Fig11KaonPlusCent5_15"]->fill(p.pT()/GeV);
            _h["Fig12KaonPlusCent5_15"]->fill(p.pT()/GeV);
            _h["KaonPlusCentrality5_15"]->fill(p.pT()/GeV);
            _h["Table4KaonPlusCentrality5_15"]->fill(p.pT()/GeV);
            }
          else if(cent >= 15. && cent < 30.)  {
            _h["KaonPlusCent15_30"]->fill(p.pT()/GeV);
            _h["Fig11KaonPlusCent15_30"]->fill(p.pT()/GeV);
            _h["Fig12KaonPlusCent15_30"]->fill(p.pT()/GeV);
            _h["KaonPlusParticle15_30"]->fill(p.pT()/GeV);
            _h["KaonPlusCentrality15_30"]->fill(p.pT()/GeV);
            _h["Table4KaonPlusCentrality15_30"]->fill(p.pT()/GeV);
            }
          else if(cent >= 30. && cent < 60.)  {
            _h["KaonPlusCent30_60"]->fill(p.pT()/GeV);
            _h["Fig10KaonPlusCent30_60"]->fill(p.pT()/GeV);
            _h["Fig11KaonPlusCent30_60"]->fill(p.pT()/GeV);
            _h["Fig12KaonPlusCent30_60"]->fill(p.pT()/GeV);
            _h["KaonPlusParticle30_60"]->fill(p.pT()/GeV);
            _h["Table4KaonPlusCentrality30_60"]->fill(p.pT()/GeV);
            }
          else if(cent >= 60. && cent < 92.)  {
            _h["KaonPlusCent60_92"]->fill(p.pT()/GeV);
            _h["Fig10KaonPlusCent60_92"]->fill(p.pT()/GeV);
            _h["Fig11KaonPlusCent60_92"]->fill(p.pT()/GeV);
            _h["Fig12KaonPlusCent60_92"]->fill(p.pT()/GeV);
            _h["KaonPlusParticle60_92"]->fill(p.pT()/GeV);
            _h["KaonPlusCentrality60_92"]->fill(p.pT()/GeV);
            _h["Table4KaonPlusCentrality60_92"]->fill(p.pT()/GeV);
            }


        }
        if(p.pid() == -321)
        {
          _h["KaonMinusMinBias"]->fill(p.pT()/GeV);
          _h["Fig11KaonMinusMinBias"]->fill(p.pT()/GeV);
          _h["KaonMinusParticle"]->fill(cent);

          if(cent >= 0. && cent < 5.) {
            _h["KaonMinusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig10KaonMinusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12KaonMinus0_5"]->fill(p.pT()/GeV);
            _h["KaonMinusCentrality0_5"]->fill(p.pT()/GeV);
            _h["Table4KaonMinusCentrality0_5"]->fill(p.pT()/GeV);
            }
          else if(cent >= 5. && cent < 15.)  {
            _h["KaonMinusCent5_15"]->fill(p.pT()/GeV);
            _h["Fig10KaonMinusCent5_15"]->fill(p.pT()/GeV);
            _h["Fig12KaonMinus5_15"]->fill(p.pT()/GeV);
            _h["KaonMinusCentrality5_15"]->fill(p.pT()/GeV);
            _h["Table4KaonMinusCentrality5_15"]->fill(p.pT()/GeV);
            }
          else if(cent >= 15. && cent < 30.)  {
            _h["KaonMinusCent15_30"]->fill(p.pT()/GeV);
            _h["Fig10KaonMinusCent15_30"]->fill(p.pT()/GeV);
            _h["Fig12KaonMinus15_30"]->fill(p.pT()/GeV);
            _h["KaonMinusCentrality15_30"]->fill(p.pT()/GeV);
            _h["Table4KaonMinusCentrality15_30"]->fill(p.pT()/GeV);
            }
            else if(cent >= 30. && cent < 60.)  {
              _h["KaonMinusCent30_60"]->fill(p.pT()/GeV);
              _h["Fig10KaonMinusCent30_60"]->fill(p.pT()/GeV);
              _h["Fig12KaonMinus30_60"]->fill(p.pT()/GeV);
              _h["KaonMinusCentrality30_60"]->fill(p.pT()/GeV);
              _h["Table4KaonMinusCentrality30_60"]->fill(p.pT()/GeV);
            }
            else if(cent >= 60. && cent < 92.)  {
              _h["KaonMinusCent60_92"]->fill(p.pT()/GeV);
              _h["Fig10KaonMinusCent60_92"]->fill(p.pT()/GeV);
              _h["Fig12KaonMinus60_92"]->fill(p.pT()/GeV);
              _h["KaonMinusCentrality60_92"]->fill(p.pT()/GeV);
              _h["Table4KaonMinusCentrality60_92"]->fill(p.pT()/GeV);
            }
        }
        if(p.pid() == 2212)
        {
          _h["ProtonMinBias"]->fill(p.pT()/GeV);
          _h["Fig11ProtonMinBias"]->fill(p.pT()/GeV);
          _h["ProtonParticle"]->fill(cent);

          if(cent >= 0. && cent < 5.) {
            _h["Fig10ProtonCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12ProtonCent0_5"]->fill(p.pT()/GeV);
          }
          else if(cent >= 5. && cent < 15.)  {
            _h["Fig10ProtonCent5_15"]->fill(p.pT()/GeV);
            _h["Fig12ProtonCent5_ 5"]->fill(p.pT()/GeV);
          }
          else if(cent >= 15. && cent < 30.)  {
            _h["Fig10ProtonCent15_30"]->fill(p.pT()/GeV);
            _h["Fig12ProtonCent15_30"]->fill(p.pT()/GeV);
          }
          else if(cent >= 30. && cent < 60.)  {
            _h["Fig10ProtonCent30_60"]->fill(p.pT()/GeV);
            _h["Fig12ProtonCent30_60"]->fill(p.pT()/GeV);
          }
          else if(cent >= 60. && cent < 92.)  {
            _h["Fig10ProtonCent60_92"]->fill(p.pT()/GeV);
            _h["Fig12ProtonCent60_92"]->fill(p.pT()/GeV);
          }
        }
        if(p.pid() == -2212)
        {
          _h["AntiProtonMinBias"]->fill(p.pT()/GeV);
          _h["Fig11AntiProtonMinBias"]->fill(p.pT()/GeV);
          _h["AntiProtonParticle"]->fill(cent);

          if(cent >= 0. && cent < 5.) {
            _h["Fig10AntiProtonCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12AntiProtonCent0_5"]->fill(p.pT()/GeV);
            }
            else if(cent >= 5. && cent < 15.)  {
              _h["Fig10AntiProtonCent5_15"]->fill(p.pT()/GeV);
              _h["Fig12AntiProtonCent5_15"]->fill(p.pT()/GeV);
            }
            else if(cent >= 15. && cent < 30.)  {
              _h["Fig10AntiProtonCent15_30"]->fill(p.pT()/GeV);
              _h["Fig12AntiProtonCent15_30"]->fill(p.pT()/GeV);
            }
            else if(cent >= 30. && cent < 60.)  {
              _h["Fig10AntiProtonCent30_60"]->fill(p.pT()/GeV);
              _h["Fig12AntiProtonCent30_60"]->fill(p.pT()/GeV);
            }
            else if(cent >= 60. && cent < 92.)  {
              _h["Fig10AntiProtonCent60_92"]->fill(p.pT()/GeV);
              _h["Fig12AntiProtonCent60_92"]->fill(p.pT()/GeV);
            }
           }

        if(p.pid() == 211)
        {
          _h["PionPlusMinBias"]->fill(p.pT()/GeV);
          _h["Fig11PionPlusMinBias"]->fill(p.pT()/GeV);
          _h["PionPlusParticle"]->fill(cent);

          if(cent >= 0. && cent < 5.) {
            _h["Fig10PionPlusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent0_5"]->fill(p.pT()/GeV);
            _h["PionPlusCentrality0_5"]->fill(p.pT()/GeV);
            _h["Table4PionPlusCentrality0_5"]->fill(p.pT()/GeV);
          }
          else if(cent >= 5. && cent < 15.)  {
            _h["Fig10PionPlusCent5_15"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent5_15"]->fill(p.pT()/GeV);
            _h["PionPlusCentrality5_15"]->fill(p.pT()/GeV);
            _h["Table4PionPlusCentrality5_15"]->fill(p.pT()/GeV);
          }
          else if(cent >= 15. && cent < 30.)  {
            _h["Fig10PionPlusCent15_30"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent15_30"]->fill(p.pT()/GeV);
            _h["PionPlusCentrality15_30"]->fill(p.pT()/GeV);
            _h["Table4PionPlusCentrality15_30"]->fill(p.pT()/GeV);
          }
          else if(cent >= 30. && cent < 60.)  {
            _h["Fig10PionPlusCent30_60"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent30_60"]->fill(p.pT()/GeV);
            _h["PionPlusCentrality30_60"]->fill(p.pT()/GeV);
            _h["Table4PionPlusCentrality30_60"]->fill(p.pT()/GeV);
          }
          else if(cent >= 60. && cent < 92.)  {
            _h["Fig10PionPlusCent60_92"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent60_92"]->fill(p.pT()/GeV);
            _h["PionPlusCentrality60_92"]->fill(p.pT()/GeV);
            _h["Table4PionPlusCentrality60_92"]->fill(p.pT()/GeV);
          }
        }
        if(p.pid() == -211)
        {
          _h["PionMinusMinBias"]->fill(p.pT()/GeV);
          _h["Fig11PionMinusMinBias"]->fill(p.pT()/GeV);
          _h["PionMinusParticle"]->fill(cent);

          if(cent >= 0. && cent < 5.) {
            _h["Fig10PionMinusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent0_5"]->fill(p.pT()/GeV);
            _h["PionMinusCentrality0_5"]->fill(p.pT()/GeV);
            _h["Table4PionMinusCentrality0_5"]->fill(p.pT()/GeV);
          }
          else if(cent >= 5. && cent < 15.)  {
            _h["Fig10PionMinusCent5_15"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent5_15"]->fill(p.pT()/GeV);
            _h["PionMinusCentrality5_15"]->fill(p.pT()/GeV);
            _h["Table4PionMinusCentrality5_15"]->fill(p.pT()/GeV);
          }
          else if(cent >= 15. && cent < 30.)  {
            _h["Fig10PionMinusCent15_30"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent15_30"]->fill(p.pT()/GeV);
            _h["PionMinusCentrality15_30"]->fill(p.pT()/GeV);
            _h["Table4PionMinusCentrality15_30"]->fill(p.pT()/GeV);
          }
          else if(cent >= 30. && cent < 60.)  {
            _h["Fig10PionMinusCent30_60"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent30_60"]->fill(p.pT()/GeV);
            _h["PionMinusCentrality30_60"]->fill(p.pT()/GeV);
            _h["Table4PionMinusCentrality30_60"]->fill(p.pT()/GeV);
          }
          else if(cent >= 60. && cent < 92.)  {
            _h["Fig10PionMinusCent60_92"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent60_92"]->fill(p.pT()/GeV);
            _h["PionMinusCentrality60_92"]->fill(p.pT()/GeV);
            _h["Table4PionMinusCentrality60_92"]->fill(p.pT()/GeV);
          }
        }
        }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      _h["KaonPlusMinBias"]->scaleW(1./_c["MinBias"]->sumW());
      _h["KaonMinusMinBias"]->scaleW(1./_c["MinBias"]->sumW());
      _h["ProtonMinBias"]->scaleW(1./_c["MinBias"]->sumW());
      _h["AntiProtonMinBias"]->scaleW(1./_c["MinBias"]->sumW());
      _h["PionPlusMinBias"]->scaleW(1./_c["MinBias"]->sumW());
      _h["PionMinusMinBias"]->scaleW(1./_c["MinBias"]->sumW());
      _h["Fig10KaonPlusCent0_5"] ->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig10KaonPlusCent5_15"] ->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig10KaonPlusCent15_30"] ->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig10KaonPlusCent30_60"] ->scaleW(1./_c["Cent30_60"]->sumW());
      _h["Fig10KaonPlusCent60_92"] ->scaleW(1./_c["Cent60_92"]->sumW());
      _h["Fig10KaonMinusCent0_5"] ->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig10KaonMinusCent5_15"] ->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig10KaonMinusCent15_30"] ->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig10KaonMinusCent30_60"] ->scaleW(1./_c["Cent30_60"]->sumW());
      _h["Fig10KaonMinusCent60_92"] ->scaleW(1./_c["Cent60_92"]->sumW());
      _h["Fig10ProtonCent0_5"] ->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig10ProtonCent0_5"] ->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig10ProtonCent5_15"]->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig10AntiProtonCent5_15"]->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig10ProtonCent15_30"]->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig10AntiProtonCent15_30"]->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig10ProtonCent30_60"] ->scaleW(1./_c["Cent30_60"]->sumW());
      _h["Fig10ProtonCent60_92"]->scaleW(1./_c["Cent60_92"]->sumW());
      _h["Fig10AntiProtonCent60_92"]->scaleW(1./_c["Cent60_92"]->sumW());
      _h["Fig11PionPlusMinBias"]->scaleW(1./_c["MinBias"]->sumW());
      _h["Fig11PionMinusMinBias"]->scaleW(1./_c["MinBias"]->sumW());
      _h["Fig11KaonPlusMinBias"] ->scaleW(1./_c["MinBias"]->sumW());
      _h["Fig11KaonMinusMinBias"]->scaleW(1./_c["MinBias"]->sumW());
      _h["Fig11ProtonMinBias"]->scaleW(1./_c["MinBias"]->sumW());
      _h["Fig11AntiProtonMinBias"]->scaleW(1./_c["MinBias"]->sumW());
      _h["Fig12PionPlusCent0_5"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig12PionPlusCent5_15"]->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig12PionPlusCent15_30"]->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig12PionPlusCent30_60"]->scaleW(1./_c["Cent30_60"]->sumW());
      _h["Fig12PionMinusCent0_5"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig12PionMinusCent5_15"]->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig12PionMinusCent15_30"]->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig12PionMinusCent30_60"]->scaleW(1./_c["Cent30_60"]->sumW());
      _h["Fig12PionPlusCent60_92"]->scaleW(1./_c["Cent60_92"]->sumW());
      _h["Fig12PionMinusCent60_92"]->scaleW(1./_c["Cent60_92"]->sumW());
      _h["Fig12KaonPlusCent0_5"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig12KaonPlusCent5_15"]->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig12KaonPlusCent15_30"]->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig12KaonPlusCent30_60"]->scaleW(1./_c["Cent30_60"]->sumW());
      _h["Fig12KaonPlusCent60_92"]->scaleW(1./_c["Cent60_92"]->sumW());
      _h["Fig12KaonMinusCent0_5"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig12KaonMinusCent5_15"]->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig12KaonMinusCent15_30"]->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig12KaonMinusCent30_60"]->scaleW(1./_c["Cent30_60"]->sumW());
      _h["Fig12KaonMinusCent60_92"]->scaleW(1./_c["Cent60_92"]->sumW());
      _h["Fig12ProtonCent0_5"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig12ProtonCent5_15"]->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig12ProtonCent15_30"]->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig12ProtonCent30_60"]->scaleW(1./_c["Cent30_60"]->sumW());
      _h["Fig12AntiProtonCent0_5"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig12AntiProtonCent5_15"]->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig12AntiProtonCent15_30"]->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig12AntiProtonCent30_60"]->scaleW(1./_c["Cent30_60"]->sumW());
      _h["Fig12ProtonCent60_92"]->scaleW(1./_c["Cent60_92"]->sumW());
      _h["Fig12AntiProtonCent60_92"]->scaleW(1./_c["Cent60_92"]->sumW());
/*
      vector<int> centBins = {0, 5, 15, 30, 60, 92};

      for(int i = 1; i <= _h["PionPlusParticle"]->numBins(); i++)
      {
              string scounter = "Cent" + to_string(centBins[i-1]) + "_" + to_string(centBins[i]);
              double normCent = 1./_c[scounter]->sumW();
              _h["PionPlusParticle"]->bin(i).scaleW(normCent);
              _h["PionMinusParticle"]->bin(i).scaleW(normCent);
              _h["KaonPlusParticle"]->bin(i).scaleW(normCent);
              _h["KaonMinusParticle"]->bin(i).scaleW(normCent);
              _h["ProtonParticle"]->bin(i).scaleW(normCent);
              _h["AntiProtonParticle"]->bin(i).scaleW(normCent);
      }

*/



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
