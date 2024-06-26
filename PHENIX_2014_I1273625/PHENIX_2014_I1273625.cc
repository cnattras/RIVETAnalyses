// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "../Centralities/RHICCentrality.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2014_I1273625 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2014_I1273625);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projection
      // The basic final-state projection:
    
      const FinalState fs(Cuts::abseta < 4.9);
      declare(fs,"fs");

      const ParticlePair& beam = beams();

      
      beamOpt = getOption<string>("beam","NONE");

      if (beamOpt == "NONE") {
      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AUAU;
      else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = PP;
      else if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000791970) collSys = dAU200;
      else if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000010020) collSys = dAU200;
      }


      if (beamOpt =="PP") collSys = PP;
      else if (beamOpt == "dAU200") collSys = dAU200;
      else if (beamOpt == "AUAU") collSys = AUAU;
      
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
    
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(h["ETEMC62AUAU"], 1, 1, 1);
      book(h["ETEMC62AUAU05"], 2, 1, 1);
      book(h["ETEMC62AUAU510"], 3, 1, 1);
      book(h["ETEMC62AUAU1015"], 4, 1, 1);
      book(h["ETEMC62AUAU1520"], 5, 1, 1);
      book(h["ETEMC62AUAU2025"], 6, 1, 1);
      book(h["ETEMC62AUAU2530"], 7, 1, 1);
      book(h["ETEMC62AUAU3035"], 8, 1, 1);
      book(h["ETEMC62AUAU3540"], 9, 1, 1);
      book(h["ETEMC62AUAU4045"], 10, 1, 1);
      book(h["ETEMC62AUAU4550"], 11, 1, 1);
      book(h["ETEMC62AUAU5055"], 12, 1, 1);
      book(h["DETDETAAUAU"], 13, 1, 1);
      book(h["DETDETAPP"], 14, 1, 1);
      book(h["DETDETADAU"], 15, 1, 1);
      book(h["QUARK200AUAU"], 16, 1, 1);
      book(h["QUARK130AUAU"], 17, 1, 1);
      book(h["QUARK62AUAU"], 18, 1, 1);
      book(h["RATIO200AUAU"], 19, 1, 1);
      book(h["RATIO130AUAU"], 20, 1, 1);
      book(h["RATIO62AUAU"], 21, 1, 1);
      book(h["NORMALDETDETA200AUAU"], 22, 1, 1);
      book(h["NORMALDETDETA130AUAU"], 23, 1, 1);
      book(h["NORMALDETDETA62AUAU"], 24, 1, 1);
      book(h["NORMALQUARK200AUAU"], 25, 1, 1);
      book(h["NORMALQUARK130AUAU"], 26, 1, 1);
      book(h["NORMALQUARK62AUAU"], 27, 1, 1);
      book(h["DETDETAQUARK200AUAU"], 28, 1, 1);
      book(h["DETDETAQUARK130AUAU"], 29, 1, 1);
      book(h["DETDETAQUARK62AUAU"], 30, 1, 1);
      book(h["PP200GAMMA"], 31, 1, 1);
      book(h["QUARKPART200AUAU"], 32, 1, 1);
      book(h["DECONPP200"], 33, 1, 1);
      book(h["ETAUAU200"], 34, 1, 1);
      book(h["DAUAQMNQP"], 35, 1, 1);
      book(h["ETDAUNQP"], 36, 1, 1);
      book(h["CHECKAUAU200"], 37, 1, 1);
      book(h["CHECKDAU200"], 38, 1, 1);
      book(h["ETAUAU200NQP"], 39, 1, 1);
      book(h["DETDETAAUAU200NQP"], 40, 1, 1);
      book(h["DETDETADAU200NQP"], 41, 1, 1);
      book(h["NPARTAUAU200"], 42, 1, 1);
      book(h["NAPRTPP200"], 43, 1, 1);

    }

    void analyze(const Event& event) {

      
    }
    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h["XXXX"]); // normalize to unity
      //normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in pb (no cuts)
      scale(h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> h;
    map<string, Profile1DPtr> p;
    map<string, CounterPtr> _h;
    string beamOpt;
    enum CollisionSystem { PP, AUAU, dAU200 };
    CollisionSystem collSys;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PHENIX_2014_I1273625);

}
