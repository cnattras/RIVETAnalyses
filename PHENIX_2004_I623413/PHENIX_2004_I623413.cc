// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2004_I623413 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2004_I623413);


    /// @name Analysis methods
    ///@{


//create binShift function
void binShift(YODA::Histo1D& histogram) {
    std::vector<YODA::HistoBin1D> binlist = histogram.bins();
    int n = 0;
    for (YODA::HistoBin1D bins : binlist) {
        double p_high = bins.xMax();
        double p_low = bins.xMin();
        //Now calculate f_corr
        if (bins.xMin() == binlist[0].xMin()) { //Check if we are working with first bin
            float b = 1 / (p_high - p_low) * log(binlist[0].height()/binlist[1].height());
            float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
            histogram.bin(n).scaleW(f_corr);
            n += 1;
        } else if (bins.xMin() == binlist.back().xMin()){ //Check if we are working with last bin
            float b = 1 / (p_high - p_low) * log(binlist[binlist.size()-2].height() / binlist.back().height());
            float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
            histogram.bin(n).scaleW(f_corr);
        } else { //Check if we are working with any middle bin
            float b = 1 / (p_high - p_low) * log(binlist[n-1].height() / binlist[n+1].height());
            float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
            histogram.bin(n).scaleW(f_corr);
            n += 1;
        }
    }
}
    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const ALICE::PrimaryParticles cp(Cuts::absrap < 0.5 && Cuts::pT > 0.5*GeV && Cuts::pT < 9*GeV);
      declare(cp,"cp");  
      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      beamOpt = getOption<string>("beam", "NONE");
      const ParticlePair& beam = beams();

     
      if (beamOpt == "NONE"){
      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu130;
      }
      
      if (beamOpt == "AUAU130") collSys = AuAu130; 

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
      book(_h["KaonPlusMinBias"], 1, 1, 1);
      book(_h["KaonMinusMinBias"], 1, 1, 2);
      book(_h["PionPlusMinBias"], 2, 1, 1);
      book(_h["PionMinusMinBias"], 2, 1, 2);
      book(_h["ProtonMinBias"], 3, 1, 1);
      book(_h["AntiProtonMinBias"], 3, 1, 2);
      book(_h["Fig10PionPlusCent0_5"], 4, 1, 1);
      book(_h["Fig10PionPlusCent5_15"], 4, 1, 2);
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

      string refnameRcpPionPlus = mkAxisCode(15,1,1);
      const Scatter2D& refdataRcpPionPlus =refData(refnameRcpPionPlus);
      book(_h["Fig11PionPlusCentral"], refnameRcpPionPlus + "_Central", refdataRcpPionPlus);
      book(_h["Fig11PionPlusPeripheral"], refnameRcpPionPlus + "_Peripheral", refdataRcpPionPlus);
      book(_s["Fig11PionPlus"], refnameRcpPionPlus);

      string refnameRcpPionMinus = mkAxisCode(15,1,2);
      const Scatter2D& refdataRcpPionMinus =refData(refnameRcpPionMinus);
      book(_h["Fig11PionMinusCentral"], refnameRcpPionMinus + "_Central", refdataRcpPionMinus);
      book(_h["Fig11PionMinusPeripheral"], refnameRcpPionMinus + "_Peripheral", refdataRcpPionMinus);
      book(_s["Fig11PionMinus"], refnameRcpPionMinus);

      string refnameRcpKaonPlus = mkAxisCode(16,1,1);
      const Scatter2D& refdataRcpKaonPlus =refData(refnameRcpKaonPlus);
      book(_h["Fig11KaonPlusCentral"], refnameRcpKaonPlus + "_Central", refdataRcpKaonPlus);
      book(_h["Fig11KaonPlusPeripheral"], refnameRcpKaonPlus + "_Peripheral", refdataRcpKaonPlus);
      book(_s["Fig11KaonPlus"], refnameRcpKaonPlus);

      string refnameRcpKaonMinus = mkAxisCode(16,1,2);
      const Scatter2D& refdataRcpKaonMinus =refData(refnameRcpKaonMinus);
      book(_h["Fig11KaonMinusCentral"], refnameRcpKaonMinus + "_Central", refdataRcpKaonMinus);
      book(_h["Fig11KaonMinusPeripheral"], refnameRcpKaonMinus + "_Peripheral", refdataRcpKaonMinus);
      book(_s["Fig11KaonMinus"], refnameRcpKaonMinus);

      string refnameRcpProton = mkAxisCode(17,1,1);
      const Scatter2D& refdataRcpProton =refData(refnameRcpProton);
      book(_h["Fig11ProtonCentral"], refnameRcpProton + "_Central", refdataRcpProton);
      book(_h["Fig11ProtonPeripheral"], refnameRcpProton + "_Peripheral", refdataRcpProton);
      book(_s["Fig11Proton"], refnameRcpProton);

      string refnameRcpAntiProton = mkAxisCode(17,1,2);
      const Scatter2D& refdataRcpAntiProton =refData(refnameRcpAntiProton);
      book(_h["Fig11AntiProtonCentral"], refnameRcpAntiProton + "_Central", refdataRcpAntiProton);
      book(_h["Fig11AntiProtonPeripheral"], refnameRcpAntiProton + "_Peripheral", refdataRcpAntiProton);
      book(_s["Fig11AntiProton"], refnameRcpAntiProton);

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
      book(_p["PionPlusCentrality"], 28, 1, 1);
      book(_p["PionMinusCentrality"], 28, 1, 2);
      book(_p["KaonPlusCentrality"], 28, 1, 3);
      book(_p["KaonMinusCentrality"], 28, 1, 4);
      book(_p["ProtonCentrality"], 28, 1, 5);
      book(_p["AntiProtonCentrality"], 28, 1, 6);
      
      // book(_h["Table4PionPlusCentrality"], 29, 1, 1);
      // book(_h["Table4PionMinusCentrality"], 29, 1, 2);
      // book(_h["Table4KaonPlusCentrality"], 29, 1, 3);
      // book(_h["Table4KaonMinusCentrality"], 29, 1, 4);
      // book(_h["Table4ProtonCentrality"], 29, 1, 5);
      // book(_h["Table4AntiProtonCentrality"], 29, 1, 6);
      
      book(_c["MinBias"], "_MinBias");
      book(_c["Cent0_5"], "_Cent0_5");
      book(_c["Cent5_15"], "_Cent5_15");
      book(_c["Cent15_30"], "_Cent15_30");
      book(_c["Cent30_60"], "_Cent30_60");
      book(_c["Cent60_92"], "_Cent60_92");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
      const double cent = centProj();
      Particles particles = applyProjection<PrimaryParticles>(event,"cp").particles();
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
          _h["KaonPlusParticle"]->fill(cent);
          _p["KaonPlusCentrality"]->fill(cent, p.pT()/GeV);

          if(cent >= 0. && cent < 5.) {
            _h["Fig10KaonPlusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12KaonPlusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig11KaonPlusCentral"]->fill(p.pT()/GeV);

            }
          else if(cent >= 5. && cent < 15.)  {
            _h["Fig10KaonPlusCent5_15"]->fill(p.pT()/GeV);
            _h["Fig12KaonPlusCent5_15"]->fill(p.pT()/GeV);
            }
          else if(cent >= 15. && cent < 30.)  {
            _h["Fig12KaonPlusCent15_30"]->fill(p.pT()/GeV);
            }
          else if(cent >= 30. && cent < 60.)  {
            _h["Fig10KaonPlusCent30_60"]->fill(p.pT()/GeV);
            _h["Fig12KaonPlusCent30_60"]->fill(p.pT()/GeV);
            }
          else if(cent >= 60. && cent < 92.)  {
            _h["Fig10KaonPlusCent60_92"]->fill(p.pT()/GeV);
            _h["Fig12KaonPlusCent60_92"]->fill(p.pT()/GeV);
            _h["Fig11KaonPlusPeripheral"]->fill(p.pT()/GeV);
            }


        }
        if(p.pid() == -321)
        {
          _h["KaonMinusMinBias"]->fill(p.pT()/GeV);
          _h["KaonMinusParticle"]->fill(cent);
          _p["KaonMinusCentrality"]->fill(cent, p.pT()/GeV);

          if(cent >= 0. && cent < 5.) {
            _h["Fig10KaonMinusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12KaonMinusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig11KaonMinusCentral"]->fill(p.pT()/GeV);

            //_h["Table4KaonMinusCentrality0_5"]->fill(p.pT()/GeV);
            }
          else if(cent >= 5. && cent < 15.)  {
            _h["Fig10KaonMinusCent5_15"]->fill(p.pT()/GeV);
            _h["Fig12KaonMinusCent5_15"]->fill(p.pT()/GeV);
            //_h["Table4KaonMinusCentrality5_15"]->fill(p.pT()/GeV);
            }
          else if(cent >= 15. && cent < 30.)  {
            _h["Fig10KaonMinusCent15_30"]->fill(p.pT()/GeV);
            _h["Fig12KaonMinusCent15_30"]->fill(p.pT()/GeV);
            //_h["Table4KaonMinusCentrality15_30"]->fill(p.pT()/GeV);
            }
            else if(cent >= 30. && cent < 60.)  {
              _h["Fig10KaonMinusCent30_60"]->fill(p.pT()/GeV);
              _h["Fig12KaonMinusCent30_60"]->fill(p.pT()/GeV);
             // _h["Table4KaonMinusCentrality30_60"]->fill(p.pT()/GeV);
            }
            else if(cent >= 60. && cent < 92.)  {
              _h["Fig10KaonMinusCent60_92"]->fill(p.pT()/GeV);
              _h["Fig12KaonMinusCent60_92"]->fill(p.pT()/GeV);
              _h["Fig11KaonMinusPeripheral"]->fill(p.pT()/GeV);
             // _h["Table4KaonMinusCentrality60_92"]->fill(p.pT()/GeV);
            }
        }
        if(p.pid() == 2212)
        {
          _h["ProtonMinBias"]->fill(p.pT()/GeV);
          _h["ProtonParticle"]->fill(cent);
          _p["ProtonCentrality"]->fill(cent, p.pT()/GeV);

          if(cent >= 0. && cent < 5.) {
            _h["Fig10ProtonCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12ProtonCent0_5"]->fill(p.pT()/GeV);
            _h["Fig11ProtonCentral"]->fill(p.pT()/GeV);
          }
          else if(cent >= 5. && cent < 15.)  {
            _h["Fig10ProtonCent5_15"]->fill(p.pT()/GeV);
            _h["Fig12ProtonCent5_15"]->fill(p.pT()/GeV);
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
            _h["Fig11ProtonPeripheral"]->fill(p.pT()/GeV);
          }
        }
        if(p.pid() == -2212)
        {
          _h["AntiProtonMinBias"]->fill(p.pT()/GeV);
          _h["AntiProtonParticle"]->fill(cent);
          _p["AntiProtonCentrality"]->fill(cent, p.pT()/GeV);

          if(cent >= 0. && cent < 5.) {
            _h["Fig10AntiProtonCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12AntiProtonCent0_5"]->fill(p.pT()/GeV);
            _h["Fig11AntiProtonCentral"]->fill(p.pT()/GeV);
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
              _h["Fig11AntiProtonPeripheral"]->fill(p.pT()/GeV);
            }
           }

        if(p.pid() == 211)
        {

          _h["PionPlusMinBias"]->fill(p.pT()/GeV);
          _h["PionPlusParticle"]->fill(cent);
          _p["PionPlusCentrality"]->fill(cent, p.pT()/GeV);

          if(cent >= 0. && cent < 5.) {
            _h["Fig10PionPlusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig11PionPlusCentral"]->fill(p.pT()/GeV);

          }
          else if(cent >= 5. && cent < 15.)  {
            _h["Fig10PionPlusCent5_15"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent5_15"]->fill(p.pT()/GeV);

          }
          else if(cent >= 15. && cent < 30.)  {
            _h["Fig10PionPlusCent15_30"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent15_30"]->fill(p.pT()/GeV);
          }
          else if(cent >= 30. && cent < 60.)  {
            _h["Fig10PionPlusCent30_60"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent30_60"]->fill(p.pT()/GeV);
          }
          else if(cent >= 60. && cent < 92.)  {
            _h["Fig10PionPlusCent60_92"]->fill(p.pT()/GeV);
            _h["Fig12PionPlusCent60_92"]->fill(p.pT()/GeV);
            _h["Fig11PionPlusPeripheral"]->fill(p.pT()/GeV);
          }
        }
        if(p.pid() == -211)
        {
          _h["PionMinusMinBias"]->fill(p.pT()/GeV);
          _h["PionMinusParticle"]->fill(cent);
          _p["PionMinusCentrality"]->fill(cent, p.pT()/GeV);

          if(cent >= 0. && cent < 5.) {
            _h["Fig10PionMinusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig12PionMinusCent0_5"]->fill(p.pT()/GeV);
            _h["Fig11PionMinusCentral"]->fill(p.pT()/GeV);

          }
          else if(cent >= 5. && cent < 15.)  {
            _h["Fig10PionMinusCent5_15"]->fill(p.pT()/GeV);
            _h["Fig12PionMinusCent5_15"]->fill(p.pT()/GeV);
          }
          else if(cent >= 15. && cent < 30.)  {
            _h["Fig10PionMinusCent15_30"]->fill(p.pT()/GeV);
            _h["Fig12PionMinusCent15_30"]->fill(p.pT()/GeV);
          }
          else if(cent >= 30. && cent < 60.)  {
            _h["Fig10PionMinusCent30_60"]->fill(p.pT()/GeV);
            _h["Fig12PionMinusCent30_60"]->fill(p.pT()/GeV);
          }
          else if(cent >= 60. && cent < 92.)  {
            _h["Fig10PionMinusCent60_92"]->fill(p.pT()/GeV);
            _h["Fig12PionMinusCent60_92"]->fill(p.pT()/GeV);
            _h["Fig11PionMinusPeripheral"]->fill(p.pT()/GeV);
          }
        }
        }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

    binShift(*_h["KaonPlusMinBias"]);
    binShift(*_h["KaonMinusMinBias"]);
    binShift(*_h["ProtonMinBias"]);
    binShift(*_h["AntiProtonMinBias"]);
    binShift(*_h["PionPlusMinBias"]);
    binShift(*_h["PionMinusMinBias"]);
    binShift(*_h["Fig10KaonPlusCent0_5"]);
    binShift(*_h["Fig10KaonPlusCent5_15"]);
    binShift(*_h["Fig10KaonPlusCent15_30"]);
    binShift(*_h["Fig10KaonPlusCent30_60"]);
    binShift(*_h["Fig10KaonPlusCent60_92"]);
    binShift(*_h["Fig10KaonMinusCent0_5"]);
    binShift(*_h["Fig10KaonMinusCent5_15"]);
    binShift(*_h["Fig10KaonMinusCent15_30"]);
    binShift(*_h["Fig10KaonMinusCent30_60"]);
    binShift(*_h["Fig10KaonMinusCent60_92"]);
    binShift(*_h["Fig10ProtonCent0_5"]);
    binShift(*_h["Fig10AntiProtonCent0_5"]);
    binShift(*_h["Fig10ProtonCent5_15"]);
    binShift(*_h["Fig10AntiProtonCent5_15"]);
    binShift(*_h["Fig10ProtonCent15_30"]);
    binShift(*_h["Fig10AntiProtonCent15_30"]);
    binShift(*_h["Fig10ProtonCent30_60"]);
    binShift(*_h["Fig10AntiProtonCent30_60"]);
    binShift(*_h["Fig10ProtonCent60_92"]);
    binShift(*_h["Fig10AntiProtonCent60_92"]);
    binShift(*_h["Fig11PionPlusCentral"]);
    binShift(*_h["Fig11PionPlusPeripheral"]);
    binShift(*_h["Fig11PionMinusCentral"]);
    binShift(*_h["Fig11PionMinusPeripheral"]);
    binShift(*_h["Fig11KaonPlusCentral"]);
    binShift(*_h["Fig11KaonPlusPeripheral"]);
    binShift(*_h["Fig11KaonMinusCentral"]);
    binShift(*_h["Fig11KaonMinusPeripheral"]);
    binShift(*_h["Fig11ProtonCentral"]);
    binShift(*_h["Fig11ProtonPeripheral"]);
    binShift(*_h["Fig11AntiProtonCentral"]);
    binShift(*_h["Fig11AntiProtonPeripheral"]);
    binShift(*_h["Fig12PionPlusCent0_5"]);
    binShift(*_h["Fig12PionPlusCent5_15"]);
    binShift(*_h["Fig12PionPlusCent15_30"]);
    binShift(*_h["Fig12PionPlusCent30_60"]);
    binShift(*_h["Fig12PionMinusCent0_5"]);
    binShift(*_h["Fig12PionMinusCent5_15"]);
    binShift(*_h["Fig12PionMinusCent15_30"]);
    binShift(*_h["Fig12PionMinusCent30_60"]);
    binShift(*_h["Fig12PionPlusCent60_92"]);
    binShift(*_h["Fig12PionMinusCent60_92"]);
    binShift(*_h["Fig12KaonPlusCent0_5"]);
    binShift(*_h["Fig12KaonPlusCent5_15"]);
    binShift(*_h["Fig12KaonPlusCent15_30"]);
    binShift(*_h["Fig12KaonPlusCent30_60"]);
    binShift(*_h["Fig12KaonPlusCent60_92"]);
    binShift(*_h["Fig12KaonMinusCent0_5"]);
    binShift(*_h["Fig12KaonMinusCent5_15"]);
    binShift(*_h["Fig12KaonMinusCent15_30"]);
    binShift(*_h["Fig12KaonMinusCent30_60"]);
    binShift(*_h["Fig12KaonMinusCent60_92"]);
    binShift(*_h["Fig12ProtonCent0_5"]);
    binShift(*_h["Fig12ProtonCent5_15"]);
    binShift(*_h["Fig12ProtonCent15_30"]);
    binShift(*_h["Fig12ProtonCent30_60"]);
    binShift(*_h["Fig12AntiProtonCent0_5"]);
    binShift(*_h["Fig12AntiProtonCent5_15"]);
    binShift(*_h["Fig12AntiProtonCent15_30"]);
    binShift(*_h["Fig12AntiProtonCent30_60"]);
    binShift(*_h["Fig12ProtonCent60_92"]);
    binShift(*_h["Fig12AntiProtonCent60_92"]);

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
      _h["Fig10AntiProtonCent0_5"] ->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig10ProtonCent5_15"]->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig10AntiProtonCent5_15"]->scaleW(1./_c["Cent5_15"]->sumW());
      _h["Fig10ProtonCent15_30"]->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig10AntiProtonCent15_30"]->scaleW(1./_c["Cent15_30"]->sumW());
      _h["Fig10ProtonCent30_60"] ->scaleW(1./_c["Cent30_60"]->sumW());
      _h["Fig10AntiProtonCent30_60"] ->scaleW(1./_c["Cent30_60"]->sumW());
      _h["Fig10ProtonCent60_92"]->scaleW(1./_c["Cent60_92"]->sumW());
      _h["Fig10AntiProtonCent60_92"]->scaleW(1./_c["Cent60_92"]->sumW());

      _h["Fig11PionPlusCentral"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig11PionPlusPeripheral"]->scaleW(1./_c["Cent60_92"]->sumW());
      divide(_h["Fig11PionPlusCentral"], _h["Fig11PionPlusPeripheral"], _s["Fig11PionPlus"]);

      _h["Fig11PionMinusCentral"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig11PionMinusPeripheral"]->scaleW(1./_c["Cent60_92"]->sumW());
      divide(_h["Fig11PionMinusCentral"], _h["Fig11PionMinusPeripheral"], _s["Fig11PionMinus"]);

      _h["Fig11KaonPlusCentral"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig11KaonPlusPeripheral"]->scaleW(1./_c["Cent60_92"]->sumW());
      divide(_h["Fig11KaonPlusCentral"], _h["Fig11KaonPlusPeripheral"], _s["Fig11KaonPlus"]);

      _h["Fig11KaonMinusCentral"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig11KaonMinusPeripheral"]->scaleW(1./_c["Cent60_92"]->sumW());
      divide(_h["Fig11KaonMinusCentral"], _h["Fig11KaonMinusPeripheral"], _s["Fig11KaonMinus"]);

      _h["Fig11ProtonCentral"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig11ProtonPeripheral"]->scaleW(1./_c["Cent60_92"]->sumW());
      divide(_h["Fig11ProtonCentral"], _h["Fig11ProtonPeripheral"], _s["Fig11Proton"]);

      _h["Fig11AntiProtonCentral"]->scaleW(1./_c["Cent0_5"]->sumW());
      _h["Fig11AntiProtonPeripheral"]->scaleW(1./_c["Cent60_92"]->sumW());
      divide(_h["Fig11AntiProtonCentral"], _h["Fig11AntiProtonPeripheral"], _s["Fig11AntiProton"]);


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

      vector<int> centBins = {0, 5, 15, 30, 60, 92};

      for(int i = 0; i < _h["PionPlusParticle"]->numBins(); i++)
      {
              string scounter = "Cent" + to_string(centBins[i]) + "_" + to_string(centBins[i+1]);
              double normCent = 1./_c[scounter]->sumW();
              _h["PionPlusParticle"]->bin(i).scaleW(normCent);
              _h["PionMinusParticle"]->bin(i).scaleW(normCent);
              _h["KaonPlusParticle"]->bin(i).scaleW(normCent);
              _h["KaonMinusParticle"]->bin(i).scaleW(normCent);
              _h["ProtonParticle"]->bin(i).scaleW(normCent);
              _h["AntiProtonParticle"]->bin(i).scaleW(normCent);
      }





    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;
    string beamOpt;
    enum CollisionSystem {AuAu130};
    CollisionSystem collSys;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(PHENIX_2004_I623413);

}
