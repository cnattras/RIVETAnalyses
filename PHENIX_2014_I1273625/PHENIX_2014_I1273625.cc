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
    
      const FinalState fs(Cuts::pT > 1*GeV && Cuts::abseta < 0.35);
      declare(fs, "fs");
  
      const ParticlePair& beam = beams();
      int NN = 0;

      
      beamOpt = getOption<string>("beam","NONE");

      if (beamOpt == "NONE") {
        if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
          NN = 197.;
          if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 1E-3)) collSys = AUAU62;
          if (fuzzyEquals(sqrtS()/GeV, 130*NN, 1E-3)) collSys = AUAU130;
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSys = AUAU200;
      }
      else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = PP200;
      else if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000791970) collSys = dAU200;
      else if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000010020) collSys = dAU200;
      }


      if (beamOpt =="PP200") collSys = PP200;
      else if (beamOpt == "dAU200") collSys = dAU200;
      else if (beamOpt == "AUAU62") collSys = AUAU62;
      else if (beamOpt == "AUAU130") collSys = AUAU130;
      else if (beamOpt == "AUAU200") collSys = AUAU200;

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
    
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)

      //Figure 2
      /*book(h["ETEMC62AUAU"], 1, 1, 1);
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

      //figure 3a
      book(h["DETDETAAUAU"], 13, 1, 1);

      //figure 3b
      book(h["DETDETAPP"], 14, 1, 1);
      book(h["DETDETADAU"], 15, 1, 1);

      //figure 4a
      book(h["QUARK200AUAU"], 16, 1, 1);
      book(h["QUARK130AUAU"], 17, 1, 1);
      book(h["QUARK62AUAU"], 18, 1, 1);

      //figure 4b
      book(h["RATIO200AUAU"], 19, 1, 1);
      book(h["RATIO130AUAU"], 20, 1, 1);
      book(h["RATIO62AUAU"], 21, 1, 1);

      //figure 5
      book(h["NORMALDETDETA200AUAU"], 22, 1, 1);
      book(h["NORMALDETDETA130AUAU"], 23, 1, 1);
      book(h["NORMALDETDETA62AUAU"], 24, 1, 1);

      //figure 6
      book(h["NORMALQUARK200AUAU"], 25, 1, 1);
      book(h["NORMALQUARK130AUAU"], 26, 1, 1);
      book(h["NORMALQUARK62AUAU"], 27, 1, 1);*/

      //figure 7 
      book(_hist_200, "d28-x01-y01", refData(28, 1, 1));
      book(_hist_130, "d29-x01-y01", refData(29, 1, 1));
      book(_hist_62, "d30-x01-y01", refData(30, 1, 1));

      //figure 8
     /* book(h["PP200GAMMA"], 31, 1, 1);

      //figure 9 
      book(h["QUARKPART200AUAU"], 32, 1, 1);

      //figure 10 
      book(h["DECONPP200"], 33, 1, 1);

      //figure 11
      book(h["ETAUAU200"], 34, 1, 1);

      //figure 12
      book(h["DAUAQMNQP"], 35, 1, 1);

      //figure 13
      book(h["ETDAUNQP"], 36, 1, 1);

      //figure 14
      book(h["CHECKAUAU200"], 37, 1, 1);
      book(h["CHECKDAU200"], 38, 1, 1);

      //figure 15
      book(h["ETAUAU200NQP"], 39, 1, 1);

      //figure 16
      book(h["DETDETAAUAU200NQP"], 40, 1, 1);
      book(h["DETDETADAU200NQP"], 41, 1, 1);

      //figure 17 
      book(h["NPARTAUAU200"], 42, 1, 1);
      book(h["NAPRTPP200"], 43, 1, 1);*/

    }

    void analyze(const Event& event) {


    double totalEt = 0;
    double deltaeta = 1; 

    Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();

    const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
        const double c = cent();
        if (c > 65) vetoEvent;
    
    for(const Particle& p : fsParticles)
        {
            totalEt += p.Et()/GeV;
        }


        if(collSys == AUAU62)
        {
          _hist_62->fill(c,totalEt/deltaeta);
        }
        else if(collSys == AUAU130)
        {
           _hist_130->fill(c,totalEt/deltaeta);
        }
        else if(collSys == AUAU200)
        {
           _hist_200->fill(c,totalEt/deltaeta); 

        }

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
    map<string, Scatter2DPtr> s;
    Profile1DPtr _hist_130;
    Profile1DPtr _hist_200;
    Profile1DPtr _hist_62;

    string beamOpt;
    enum CollisionSystem {PP200, AUAU62, AUAU130, AUAU200, dAU200};
    CollisionSystem collSys;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PHENIX_2014_I1273625);

}
