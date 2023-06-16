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
#include "../Centralities/RHICCentrality.hh" //external header for Centrality calculation
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES


namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2007_I731133 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2007_I731133);


    /// Book histograms and initialise projections before the run
      void init() {
          //Particles: eta, pi^+, pi^-, pi^0, gamma (respectively)
          std::initializer_list<int> pdgIds = { 221, 211, -211, 111, 22};
          
          //declare cuts; most of these are found in section II of the paper
          //For charged particles:
          const PrimaryParticles cp(pdgIds, Cuts::abseta < 2.2 && Cuts::abseta > 1.2 && Cuts::abscharge > 0 && Cuts::absrap < 0.35);
          declare(cp, "cp");
          
          //Uncharged particles
          //const UnstableParticles np(pdgIds, Cuts::abseta < 2.2 && Cuts::abseta > 1.2 && Cuts::abscharge == 0 && Cuts::absrap < 0.35);
          //declare(np, "np");
          
          beamOpt = getOption<string>("beam", "NONE");

          //check the collision system
          if (beamOpt == "PP") collSys = pp;
          else if (beamOpt == "AUAU200") collSys = AuAu200;
          else if (beamOpt == "dAU200") collSys = dAu200;
	
          //declaration for collision systems that are not p+p
          if (!(collSys == pp)) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

          //Counters, for N_{evt}
          book(sow["sow_pp"], "_sow_pp");
          
          //figure 13
          //d01-x01-y01
          string refname1 = mkAxisCode(1, 1, 1);
          //mkAxisCode gives us the internal histogram name for a given d, x, and y
          const Scatter2D& refdata1 = refData(refname1);
          //here we have to define a Scatter2D& for the next part, we can use the refData function on refname
          book(hCrossSec["ppMesontoGamma"], refname1 + "_pp_GammaGamma", refdata1);
          //here, we are using book() as: book(1DPtr&, const string &name, const Scatter2D), and this books a histograms with binning using d01-x01-y01 from our yoda as a reference.
          
    }



    /// Perform the per-event analysis
      void analyze(const Event& event) {
          Particles chargedParticles = applyProjection<PrimaryParticles>(event,"cp").particles();
          
          
          if (collSys == pp)
          {
              //Fill our counter; this will be useful at the end when we want to run ->sum() to find N_{evt}
              sow["sow_pp"]->fill();
              
              //a conditional for the charged particles: eta, pi+, pi-, gamma
              for (Particle p : chargedParticles) {
                  
                  
                  //define what we will fill our histograms with
                  //we technically don't have to do this, but it's convenient
                  double partPt = p.pT() / GeV;
                  double pt_weight = 1. / (partPt * 2. * M_PI);
                  
                  switch (p.pid()) {
                          //particle ids for convenience: { 221, 211, -211, 22};
                      case 221: {  //eta
                          hCrossSec["ppMesontoGamma"]->fill(partPt, pt_weight);
                          break;
                      }
                      case 211: { //pi+
                          break;
                      }
                      case -211: { //pi-
                          break;
                      }
                      case 22: { //gamma
                          break;
                      }
                  }
              }
          }
          
          //Nik, SEE HERE: this is where your first figure will go 
          if (collSys == dAu200){
              
              for (Particle p : chargedParticles) {
                  double partPt = p.pT() / GeV;
                  double pt_weight = 1. / (partPt * 2. * M_PI);
                  
                  switch (p.pid()) {
                      case 221: {  //eta
                          break;
                      }
                      case 211: { //pi+
                          break;
                      }
                      case -211: { //pi-
                          break;
                      }
                      case 22: { //gamma
                          break;
                      }
                  }
              }
          }
      }


    /// Normalise histograms etc., after the run
      void finalize() {
          double cross = crossSection() / picobarn;
          
          hCrossSec["ppMesontoGamma"]->scaleW(cross);
      }

      //histograms
      map<string, Histo1DPtr> hCrossSec;
      //map<string, Histo1DPtr> hdAuYields;
      //map<string, Histo1DPtr> hAuAuYields;
      
      //ratios
      //map<string, Scatter2DPtr> RatioAuAu;
      //map<string, Scatter2DPtr> RatiodA;
      
      //RdA, Raa
      //map<string, Scatter2DPtr> hRda;
      //map<string, Scatter2DPtr> hRaa;
      
      //Counters
      map<string, CounterPtr> sow;
      
      
      string beamOpt;
      enum CollisionSystem { pp, AuAu200, dAu200 };
      CollisionSystem collSys;
      
      //these two vectors are only initialized if we would like to for loop over centrality bins
      //vector<int> AuAuCentralityBins{ 20, 60, 92 };
      //vector<int> dAuCentralityBins{ 20, 40, 60, 88 };

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2007_I731133);

}
