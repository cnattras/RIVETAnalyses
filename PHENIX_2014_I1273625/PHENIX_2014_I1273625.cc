// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
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
    
      const FinalState fs(Cuts::abseta < 0.5);
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
    
      }


      if (beamOpt == "AUAU62") collSys = AUAU62;
      else if (beamOpt == "AUAU130") collSys = AUAU130;
      else if (beamOpt == "AUAU200") collSys = AUAU200;

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
    
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)


      //figure 7 
      book(_hist_200, "d28-x01-y01", refData(28, 1, 1));
      book(_hist_130, "d29-x01-y01", refData(29, 1, 1));
      book(_hist_62, "d30-x01-y01", refData(30, 1, 1));

      
    }

    void analyze(const Event& event) {


    double totalEt = 0;
    double deltaeta = 1; 

    Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();

    const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
        const double c = cent();
        if (collSys == AUAU62){
          if (c > 55) vetoEvent;
        }
        else{
        if (c > 65) vetoEvent;
        }

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
