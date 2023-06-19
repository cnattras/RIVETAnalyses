// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/UnstableParticles.hh"
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
          const UnstableParticles np(Cuts::abspid == 111 && Cuts::abseta < 2.2 && Cuts::abseta > 1.2 && Cuts::abscharge == 0 && Cuts::absrap < 0.35);
          declare(np, "np");
          
          beamOpt = getOption<string>("beam", "NONE");

          //check the collision system
          if (beamOpt == "PP") collSys = pp;
          else if (beamOpt == "AUAU200") collSys = AuAu200;
          else if (beamOpt == "dAU200") collSys = dAu200;
	
          //declaration for collision systems that are not p+p
          if (!(collSys == pp)) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

          //Counters
          book(sow["sow_pp"], "_sow_pp");
          book(sow["sow_dAu"], "_sow_dAu");
          book(sow["sow_dAuc0020"], "_sow_dAuc0020");
          
          
          //figure 13 (12.1 in hepdata)
          //d01-x01-y01
          string refname1 = mkAxisCode(1, 1, 1);
          //mkAxisCode gives us the internal histogram name for a given d, x, and y
          const Scatter2D& refdata1 = refData(refname1);
          //here we have to define a Scatter2D& for the next part, we can use the refData function on refname
          book(hCrossSec["ppEtatoGamma"], refname1 + "_pp_GammaGamma", refdata1);
          //here, we are using book() as: book(1DPtr&, const string &name, const Scatter2D), and this books a histograms with binning using d01-x01-y01 from our yoda as a reference.
          
          //d02-x01-y01 (12.2 in hepdata)
          string refname2 = mkAxisCode(2, 1, 1);
          const Scatter2D& refdata2 = refData(refname2);
          book(hCrossSec["ppEtatoPion"], refname2 + "_pp_Pion", refdata2);
          
          //d03-x01-y01 FOR NIK
          
          //d04-x01-y01 FOR NIK
          
          
          //Figure 14
          //d05-x01-y01
          string refname5 = mkAxisCode(5, 1, 1);
          const Scatter2D& refdata5 = refData(refname5);
          book(hdAuYields["ptyieldsdAuc0020"], refname5 + "_dAuc0020_Eta", refdata5);
          
          //d05-x01-y02
          
          //d05-x01-y03
          
          //d05-x01-y04
          
          
    }



    /// Perform the per-event analysis
      void analyze(const Event& event) {
          Particles chargedParticles = applyProjection<PrimaryParticles>(event,"cp").particles();
          Particles neutralParticles = applyProjection<UnstableParticles>(event,"np").particles();
          
          if (collSys == pp)
          {
              //Fill our counter; this will be useful at the end when we want to run ->sumW()
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
                          hCrossSec["ppEtatoGamma"]->fill(partPt, pt_weight);
                          break;
                      }
                      case 211: { //pi+
                          hCrossSec["ppEtatoPion"]->fill(partPt,pt_weight);
                          break;
                      }
                      case -211: { //pi-
                          hCrossSec["ppEtatoPion"]->fill(partPt,pt_weight);
                          break;
                      }
                      case 22: { //gamma
                          break;
                      }
                  }
              }
              
              //a conditional for the neutral particles: pi^0
              for (Particle p : neutralParticles){
                  
                  
                  //define what we will fill our histograms with
                  //we technically don't have to do this, but it's convenient
                  double partPt = p.pT() / GeV;
                  double pt_weight = 1. / (partPt * 2. * M_PI);
                  
                  hCrossSec["ppEtatoPion"]->fill(partPt, pt_weight);
              }
          }
          
          const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
          const double c = cent();
          
          
          //Nik, SEE HERE: down below is where you will do your edit for figure 13.1 and 13.2
          if (collSys == dAu200){
              
              if ((c < 0.) && (c > 100.)) vetoEvent; //veto all nonexisting events
              
              sow["sow_dAu"]->fill();
              
              for (Particle p : chargedParticles) {
                  //comment these next two lines are only commented to avoid an unused variable warning
                  //double partPt = p.pT() / GeV;
                  //double pt_weight = 1. / (partPt * 2. * M_PI);
                  
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
              
              //a conditional for the neutral particles: pi^0
              for (Particle p : neutralParticles){
                  //define what we will fill our histograms with below
                  //we technically don't have to do this, but it's convenient
                  //double partPt = p.pT() / GeV;
                  //double pt_weight = 1. / (partPt * 2. * M_PI);
                  break;
              }
              
              
              //Centrality inclusion begins here:
              if ((c >= 0.) && (c < 20.)){
                  
                  //fill our counter for this centrality
                  sow["sow_dAuc0020"]->fill();
                  
                  
                  for (Particle p : chargedParticles) {
                      //comment these next two lines are only commented to avoid an unused variable warning
                      double partPt = p.pT() / GeV;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      sow["sow_dAuc0020"]->fill();
                      
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hdAuYields["ptyieldsdAuc0020"]->fill(partPt, pt_weight);
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
          
          
      }


    /// Normalise histograms etc., after the run
      void finalize() {
          double cross = crossSection() / picobarn;
          
          //Figure 13: d01-x01-y01
          hCrossSec["ppEtatoGamma"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtatoGamma"]->scaleW(cross);
          
          //Figure 13: d02-x01-y01
          hCrossSec["ppEtatoPion"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtatoPion"]->scaleW(cross);
          
          //Figure 13: d03-x01-y01
          
          //Figure 13: d04-x01-y01
          
          //Figure 14: d05-x01-y01
          hdAuYields["ptyieldsdAuc0020"]->scaleW(1. / sow["sow_dAuc0020"]->sumW());
      }

      //histograms
      map<string, Histo1DPtr> hCrossSec;
      map<string, Histo1DPtr> hdAuYields;
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
