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
          book(sow["sow_dAuc2040"], "_sow_dAuc2040");
          book(sow["sow_dAuc4060"], "_sow_dAuc4060");
          book(sow["sow_dAuc6088"], "_sow_dAuc6088");
          
          //figure 13 (12.1 in hepdata)
          //d01-x01-y01
          string refname1 = mkAxisCode(1, 1, 1);
          //mkAxisCode gives us the internal histogram name for a given d, x, and y
          const Scatter2D& refdata1 = refData(refname1);
          //here we have to define a Scatter2D& for the next part, we can use the refData function on refname
          book(hCrossSec["ppEtatoGammaa"], refname1 + "_pp_GammaGamma", refdata1);
          //here, we are using book() as: book(1DPtr&, const string &name, const Scatter2D), and this books a histograms with binning using d01-x01-y01 from our yoda as a reference.
          
          //d02-x01-y01 (12.2 in hepdata)
          string refname2 = mkAxisCode(2, 1, 1);
          const Scatter2D& refdata2 = refData(refname2);
          book(hCrossSec["ppEtatoPiona"], refname2 + "_pp_Pion", refdata2);
          
          //d03-x01-y01 FOR NIK
          
          //d04-x01-y01 FOR NIK
          
          
          //Figure 14
          //d05-x01-y01
          string refname5 = mkAxisCode(5, 1, 1);
          const Scatter2D& refdata5 = refData(refname5);
          book(hdAuYields["ptyieldsdAuc0020a"], refname5 + "_dAuc0020_Eta", refdata5);
          
          //d05-x01-y02
          string refname6 = mkAxisCode(5, 1, 2);
          const Scatter2D& refdata6 = refData(refname6);
          book(hdAuYields["ptyieldsdAuc2040a"], refname6 + "_dAuc2040_Eta", refdata6);
          
          //d05-x01-y03
          string refname7 = mkAxisCode(5, 1, 3);
          const Scatter2D& refdata7 = refData(refname7);
          book(hdAuYields["ptyieldsdAuc4060a"], refname7 + "_dAuc4060_Eta", refdata7);
          
          //d05-x01-y04
          string refname8 = mkAxisCode(5, 1, 4);
          const Scatter2D& refdata8 = refData(refname8);
          book(hdAuYields["ptyieldsdAuc6088a"], refname8 + "_dAuc6088_Eta", refdata8);
          
          //Figure 15
          //d06-x01-y01
          
          //d06-x01-y02
          
          //d06-x01-y03
          
          //d06-x01-y04
          
          
          //Figure 16
          //d07-x01-y01
          string refname13 = mkAxisCode(7, 1, 1);
          const Scatter2D& refdata13 = refData(refname13);
          book(hdAuYields["ptyieldsdAuc0088b"], refname13 + "_dAuc0088_Eta", refdata13);
          book(hCrossSec["ppEtadAuc0088"], refname13 + "_pp_Eta", refdata13);
          book(hRda["EtadAuc0088"], refname13);
          
          //d07-x01-y02
          string refname14 = mkAxisCode(7, 1, 2);
          const Scatter2D& refdata14 = refData(refname14);
          book(hdAuYields["ptyieldsdAuc0020b"], refname14 + "_dAuc0020_Eta", refdata14);
          book(hCrossSec["ppEtadAuc0020"], refname14 + "_pp_Eta", refdata14);
          book(hRda["EtadAuc0020"], refname14);
          
          //d07-x01-y03
          string refname15 = mkAxisCode(7, 1, 3);
          const Scatter2D& refdata15 = refData(refname15);
          book(hdAuYields["ptyieldsdAuc2040b"], refname15 + "_dAuc2040_Eta", refdata15);
          book(hCrossSec["ppEtadAuc2040"], refname15 + "_pp_Eta", refdata15);
          book(hRda["EtadAuc2040"], refname15);
          
          //d07-x01-y04
          string refname16 = mkAxisCode(7, 1, 4);
          const Scatter2D& refdata16 = refData(refname16);
          book(hdAuYields["ptyieldsdAuc4060b"], refname16 + "_dAuc4060_Eta", refdata16);
          book(hCrossSec["ppEtadAuc4060"], refname16 + "_pp_Eta", refdata16);
          book(hRda["EtadAuc4060"], refname16);
          
          //d07-x01-y05
          string refname17 = mkAxisCode(7, 1, 5);
          const Scatter2D& refdata17 = refData(refname17);
          book(hdAuYields["ptyieldsdAuc6088b"], refname17 + "_dAuc6088_Eta", refdata17);
          book(hCrossSec["ppEtadAuc6088"], refname17 + "_pp_Eta", refdata17);
          book(hRda["EtadAuc6088"], refname17);
          
        
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
                          hCrossSec["ppEtatoGammaa"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc0088"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc0020"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc2040"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc4060"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc6088"]->fill(partPt, pt_weight);
                          break;
                      }
                      case 211: { //pi+
                          hCrossSec["ppEtatoPiona"]->fill(partPt,pt_weight);
                          break;
                      }
                      case -211: { //pi-
                          hCrossSec["ppEtatoPiona"]->fill(partPt,pt_weight);
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
                  
                  hCrossSec["ppEtatoPiona"]->fill(partPt, pt_weight);
              }
              
          }
          
          
          
          //Nik, SEE HERE: down below is where you will do your edit for figure 13.1 and 13.2
          if (collSys == dAu200){
            
              const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
              const double c = cent();
              
              if ((c < 0.) || (c > 88.)) vetoEvent;
              
              sow["sow_dAu"]->fill();
              
              for (Particle p : chargedParticles) {
                  //comment these next two lines are only commented to avoid an unused variable warning
                  double partPt = p.pT() / GeV;
                  double pt_weight = 1. / (partPt * 2. * M_PI);
                  
                  switch (p.pid()) {
                      case 221: {  //eta
                          hdAuYields["ptyieldsdAuc0088b"]->fill(partPt, pt_weight);
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
                      
                      double partPt = p.pT() / GeV;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      sow["sow_dAuc0020"]->fill();
                      
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hdAuYields["ptyieldsdAuc0020a"]->fill(partPt, pt_weight);
                              hdAuYields["ptyieldsdAuc0020b"]->fill(partPt, pt_weight);
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
              else if ((c >= 20.) && (c < 40.)){
                  for (Particle p : chargedParticles) {
                      
                      double partPt = p.pT() / GeV;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      sow["sow_dAuc2040"]->fill();
                      
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hdAuYields["ptyieldsdAuc2040a"]->fill(partPt, pt_weight);
                              hdAuYields["ptyieldsdAuc2040b"]->fill(partPt, pt_weight);
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
              else if ((c >= 40.) && (c < 60.)){
                  for (Particle p : chargedParticles) {
                      
                      double partPt = p.pT() / GeV;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      sow["sow_dAuc4060"]->fill();
                      
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hdAuYields["ptyieldsdAuc4060a"]->fill(partPt, pt_weight);
                              hdAuYields["ptyieldsdAuc4060b"]->fill(partPt, pt_weight);
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
              else if ((c >= 60.) && (c < 88.)){
                  for (Particle p : chargedParticles) {
                      
                      double partPt = p.pT() / GeV;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      sow["sow_dAuc6088"]->fill();
                      
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hdAuYields["ptyieldsdAuc6088a"]->fill(partPt, pt_weight);
                              hdAuYields["ptyieldsdAuc6088b"]->fill(partPt, pt_weight);
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
          
          //Figure 13:
          //d01-x01-y01
          hCrossSec["ppEtatoGammaa"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtatoGammaa"]->scaleW(cross);
          
          //d02-x01-y01
          hCrossSec["ppEtatoPiona"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtatoPiona"]->scaleW(cross);
          
          //d03-x01-y01
          
          //d04-x01-y01
          
          //Figure 14:
          //d05-x01-y01
          hdAuYields["ptyieldsdAuc0020a"]->scaleW(1. / sow["sow_dAuc0020"]->sumW());
          hdAuYields["ptyieldsdAuc2040a"]->scaleW(1. / sow["sow_dAuc2040"]->sumW());
          hdAuYields["ptyieldsdAuc4060a"]->scaleW(1. / sow["sow_dAuc4060"]->sumW());
          hdAuYields["ptyieldsdAuc6088a"]->scaleW(1. / sow["sow_dAuc6088"]->sumW());
          
          //Figure 16:
          //d07-x01-y01
          //denominator
          hCrossSec["ppEtadAuc0088"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc0088"]->scaleW(cross);
          hCrossSec["ppEtadAuc0088"]->scaleW(0.2); //scaling by <TdA> as shown in Table II
          //numerator
          hdAuYields["ptyieldsdAuc0088b"]->scaleW(1. / sow["sow_dAu"]->sumW());
          //Rda
          divide(hdAuYields["ptyieldsdAuc0088b"], hCrossSec["ppEtadAuc0088"], hRda["EtadAuc0088"]);
          
          //d07-x01-y02
          hCrossSec["ppEtadAuc0020"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc0020"]->scaleW(cross);
          hCrossSec["ppEtadAuc0020"]->scaleW(0.36); //scaling by <TdA> as shown in Table II
          //numerator
          hdAuYields["ptyieldsdAuc0020b"]->scaleW(1. / sow["sow_dAuc0020"]->sumW());
          //Rda
          divide(hdAuYields["ptyieldsdAuc0020b"], hCrossSec["ppEtadAuc0020"], hRda["EtadAuc0020"]);
          
          //d07-x01-y03
          hCrossSec["ppEtadAuc2040"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc2040"]->scaleW(cross);
          hCrossSec["ppEtadAuc2040"]->scaleW(0.25); //scaling by <TdA> as shown in Table II
          //numerator
          hdAuYields["ptyieldsdAuc2040b"]->scaleW(1. / sow["sow_dAuc2040"]->sumW());
          //Rda
          divide(hdAuYields["ptyieldsdAuc2040b"], hCrossSec["ppEtadAuc2040"], hRda["EtadAuc2040"]);
          
          //d07-x01-y04
          hCrossSec["ppEtadAuc4060"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc4060"]->scaleW(cross);
          hCrossSec["ppEtadAuc4060"]->scaleW(0.17); //scaling by <TdA> as shown in Table II
          //numerator
          hdAuYields["ptyieldsdAuc4060b"]->scaleW(1. / sow["sow_dAuc4060"]->sumW());
          //Rda
          divide(hdAuYields["ptyieldsdAuc4060b"], hCrossSec["ppEtadAuc4060"], hRda["EtadAuc4060"]);
          
          //d07-x01-y05
          hCrossSec["ppEtadAuc6088"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc6088"]->scaleW(cross);
          hCrossSec["ppEtadAuc6088"]->scaleW(0.073); //scaling by <TdA> as shown in Table II
          //numerator
          hdAuYields["ptyieldsdAuc6088b"]->scaleW(1. / sow["sow_dAuc6088"]->sumW());
          //Rda
          divide(hdAuYields["ptyieldsdAuc6088b"], hCrossSec["ppEtadAuc6088"],hRda["EtadAuc6088"]);
      }

      //histograms
      map<string, Histo1DPtr> hCrossSec;
      map<string, Histo1DPtr> hdAuYields;
      //map<string, Histo1DPtr> hAuAuYields;
      
      //ratios
      //map<string, Scatter2DPtr> RatioAuAu;
      //map<string, Scatter2DPtr> RatiodA;
      
      //RdA, Raa
      map<string, Scatter2DPtr> hRda;
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
