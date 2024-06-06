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
#include "../Centralities/RHICCentrality.hh"
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
          //Particles: eta (respectively)
          std::initializer_list<int> pdgIds = { 221 };
          
          //declare cuts; most of these are found in section II of the paper
          //For charged particles:
          //consider adding: && Cuts::pT < 12*GeV, && Cuts::abscharge == 0, && Cuts::abscharge > 0 Cuts::abseta < 2.2 && Cuts::abseta > 1.2 &&
          const ALICE::PrimaryParticles cp(Cuts::absrap < 0.5 && Cuts::pT > 1*GeV);
          declare(cp, "cp");
          
          //Uncharged particles
          const UnstableParticles np(Cuts::abspid == 111 && Cuts::absrap < 0.5 && Cuts::pT > 1*GeV);
          declare(np, "np");
          
          beamOpt = getOption<string>("beam", "NONE");
          
          if (beamOpt == "NONE") {
          if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu200;
          else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = pp;
          else if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000791970) collSys = DAu200;
          else if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000010020) collSys = DAu200;
          }

          //check the collision system
          if (beamOpt == "PP200") collSys = pp;
          else if (beamOpt == "AUAU200") collSys = AuAu200;
          else if (beamOpt == "DAU200") collSys = DAu200;
          
          //declaration for collision systems that are not p+p
          if (!(collSys == pp)) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

          //Counters
          book(sow["sow_pp"], "_sow_pp");
          book(sow["sow_dAu"], "_sow_dAu");
          book(sow["sow_dAuc0020"], "_sow_dAuc0020");
          book(sow["sow_dAuc2040"], "_sow_dAuc2040");
          book(sow["sow_dAuc4060"], "_sow_dAuc4060");
          book(sow["sow_dAuc6088"], "_sow_dAuc6088");
          book(sow["sow_AuAuc0092"], "_sow_AuAuc0092");
          book(sow["sow_AuAuc0020"], "_sow_AuAuc0020");
          book(sow["sow_AuAuc2060"], "_sow_AuAuc2060");
          book(sow["sow_AuAuc6092"], "_sow_AuAuc6092");
          
          //figure 13 (12.1 in hepdata)
          //d01-x01-y01
          string refname1 = mkAxisCode(1, 1, 1);
          //mkAxisCode gives us the internal histogram name for a given d, x, and y
          //const Scatter2D& refdata1 = refData(refname1);
          //here we have to define a Scatter2D& for the next part, we can use the refData function on refname
          book(hCrossSec["ppEtaa"], 1, 1, 1);
          //here, we are using book() as: book(1DPtr&, const string &name, const Scatter2D), and this books a histograms with binning using d01-x01-y01 from our yoda as a reference.
          
          //d02-x01-y01 We don't worry about the decay channels, so we won't need this one since it is a duplicate
          //string refname2 = mkAxisCode(2, 1, 1);
          //const Scatter2D& refdata2 = refData(refname2);
          //book(hCrossSec["ppEtatoPiona"], refname2 + "_pp_Pion", refdata2);
          
          //d03-x01-y01 FOR NIK
          string refname2 = mkAxisCode(3, 1, 1);
          //mkAxisCode gives us the internal histogram name for a given d, x, and y
          //const Scatter2D& refdata2 = refData(refname2);
          //here we have to define a Scatter2D& for the next part, we can use the refData function on refname
          book(hCrossSec["dAuEta"], 3, 1, 1);
          
          //Figure 14: Invariant yields of eta in d+Au collisions
          //d05-x01-y01: 00-20% centrality
          string refname5 = mkAxisCode(5, 1, 1);
          //const Scatter2D& refdata5 = refData(refname5);
          book(hEtaPt["ptyieldsdAuc0020a"], 5, 1, 1);
          
          //d05-x01-y02: 20-40% centrality
          string refname6 = mkAxisCode(5, 1, 2);
          //const Scatter2D& refdata6 = refData(refname6);
          book(hEtaPt["ptyieldsdAuc2040a"], 5, 1, 2);
          
          //d05-x01-y03: 40-60% centrality
          string refname7 = mkAxisCode(5, 1, 3);
          //const Scatter2D& refdata7 = refData(refname7);
          book(hEtaPt["ptyieldsdAuc4060a"], 5, 1, 3);
          
          //d05-x01-y04: 60-88% centrality
          string refname8 = mkAxisCode(5, 1, 4);
          //const Scatter2D& refdata8 = refData(refname8);
          book(hEtaPt["ptyieldsdAuc6088a"], 5, 1, 4);
          
          //Figure 15
          //d06-x01-y01 min. bias
          string refname9 = mkAxisCode(6, 1, 1);
          //const Scatter2D& refdata9 = refDa6, 1, 1
          book(hEtaPt["ptyieldsAuAuc0092a"], 6, 1, 1);         
          //d06-x01-y02 0-20%
          string refname10 = mkAxisCode(6, 1, 2);
          //const Scatter2D& refdata10 = refData(refname10);
          book(hEtaPt["ptyieldsAuAuc0020a"], 6, 1, 2);         
          //d06-x01-y03 20-40%
          string refname11 = mkAxisCode(6, 1, 3);
          //const Scatter2D& refdata11 = refData(refname11);
          book(hEtaPt["ptyieldsAuAuc2060a"], 6, 1, 3);          
          //d06-x01-y04 60-92%
          string refname12 = mkAxisCode(6, 1, 4);
          //const Scatter2D& refdata12 = refData(refname12);
          book(hEtaPt["ptyieldsAuAuc6092a"], 6, 1, 4);         
          
          //Figure 16: Rda for measured etas
          //d07-x01-y01: 00-88% centrality (minimum bias)
          string refname13 = mkAxisCode(7, 1, 1);
          const Scatter2D& refdata13 = refData(refname13);
          book(hEtaPt["ptyieldsdAuc0088b"], refname13 + "_dAuc0088_Eta", refdata13);
          book(hCrossSec["ppEtadAuc0088"], refname13 + "_pp_Eta", refdata13);
          book(hRda["EtadAuc0088"], refname13);
          
          //d07-x01-y02: 00-20% centrality
          string refname14 = mkAxisCode(7, 1, 2);
          const Scatter2D& refdata14 = refData(refname14);
          book(hEtaPt["ptyieldsdAuc0020b"], refname14 + "_dAuc0020_Eta", refdata14);
          book(hCrossSec["ppEtadAuc0020"], refname14 + "_pp_Eta", refdata14);
          book(hRda["EtadAuc0020"], refname14);
          
          //d07-x01-y03: 20-40% centrality
          string refname15 = mkAxisCode(7, 1, 3);
          const Scatter2D& refdata15 = refData(refname15);
          book(hEtaPt["ptyieldsdAuc2040b"], refname15 + "_dAuc2040_Eta", refdata15);
          book(hCrossSec["ppEtadAuc2040"], refname15 + "_pp_Eta", refdata15);
          book(hRda["EtadAuc2040"], refname15);
          
          //d07-x01-y04: 40-60% centrality
          string refname16 = mkAxisCode(7, 1, 4);
          const Scatter2D& refdata16 = refData(refname16);
          book(hEtaPt["ptyieldsdAuc4060b"], refname16 + "_dAuc4060_Eta", refdata16);
          book(hCrossSec["ppEtadAuc4060"], refname16 + "_pp_Eta", refdata16);
          book(hRda["EtadAuc4060"], refname16);
          
          //d07-x01-y05: 60-88% centrality
          string refname17 = mkAxisCode(7, 1, 5);
          const Scatter2D& refdata17 = refData(refname17);
          book(hEtaPt["ptyieldsdAuc6088b"], refname17 + "_dAuc6088_Eta", refdata17);
          book(hCrossSec["ppEtadAuc6088"], refname17 + "_pp_Eta", refdata17);
          book(hRda["EtadAuc6088"], refname17);
          
          //Figure 17
          //d08-x01-y01: 00-20% centrality
          string refname18 = mkAxisCode(8, 1, 1);
          const Scatter2D& refdata18 = refData(refname18);
          book(hEtaPt["ptyieldsAuAuc0020b"], refname18 + "_AuAuc0020_Eta", refdata18);
          book(hCrossSec["ppEtaAuAuc0020"], refname18 + "_pp_Eta", refdata18);
          book(hRaa["EtaAuAuc0020"], refname18);
          
          //d08-x01-y02: 20-60% centrality
          string refname19 = mkAxisCode(8, 1, 2);
          const Scatter2D& refdata19 = refData(refname19);
          book(hEtaPt["ptyieldsAuAuc2060b"], refname19 + "_AuAuc2060_Eta", refdata19);
          book(hCrossSec["ppEtaAuAuc2060"], refname19 + "_pp_Eta", refdata19);
          book(hRaa["EtaAuAuc2060"], refname19);
          
          //d08-x01-y03: 60-92% centrality
          string refname20 = mkAxisCode(8, 1, 3);
          const Scatter2D& refdata20 = refData(refname20);
          book(hEtaPt["ptyieldsAuAuc6092b"], refname20 + "_AuAuc6092_Eta", refdata20);
          book(hCrossSec["ppEtaAuAuc6092"], refname20 + "_pp_Eta", refdata20);
          book(hRaa["EtaAuAuc6092"], refname20);    
          
          //Figure 18: eta/pion^0 ratio in pp collisions
          //d09-x01-y01
          string refname21 = mkAxisCode(9, 1, 1);
          const Scatter2D& refdata21 = refData(refname21);
          book(hEtaPt["ptyieldsEta"], refname21 + "_pp_Eta", refdata21);
          book(hPionPt["ptyieldsPion"], refname21 + "_pp_Pion", refdata21);
          book(Ratiopp["EtaToPion"], refname21);
          
          
          //Figure 19: eta/pion^0 ratios in dAu collisions
          //d10-x01-y01: 00-88% centrality (minimum bias)
          string refname22 = mkAxisCode(10, 1, 1);
          const Scatter2D& refdata22 = refData(refname22);
          book(hEtaPt["ptyieldsdAuc0088c"], refname22 + "_pp_Eta", refdata22);
          book(hPionPt["ptyieldsdAuc0088"], refname22 + "_pp_Pion", refdata22);
          book(RatiodAu["EtaToPion0088"], refname22);
          
          //d10-x01-y02: 00-20% centrality
          string refname23 = mkAxisCode(10, 1, 2);
          const Scatter2D& refdata23 = refData(refname23);
          book(hEtaPt["ptyieldsdAuc0020c"], refname23 + "_pp_Eta", refdata23);
          book(hPionPt["ptyieldsdAuc0020"], refname23 + "_pp_Pion", refdata23);
          book(RatiodAu["EtaToPion0020"], refname23);
          
          //d10-x01-y03: 20-40% centrality
          string refname24 = mkAxisCode(10, 1, 3);
          const Scatter2D& refdata24 = refData(refname24);
          book(hEtaPt["ptyieldsdAuc2040c"], refname24 + "_pp_Eta", refdata24);
          book(hPionPt["ptyieldsdAuc2040"], refname24 + "_pp_Pion", refdata24);
          book(RatiodAu["EtaToPion2040"], refname24);
          
          //d10-x01-y04: 40-60% centrality
          string refname25 = mkAxisCode(10, 1, 4);
          const Scatter2D& refdata25 = refData(refname25);
          book(hEtaPt["ptyieldsdAuc4060c"], refname25 + "_pp_Eta", refdata25);
          book(hPionPt["ptyieldsdAuc4060"], refname25 + "_pp_Pion", refdata25);
          book(RatiodAu["EtaToPion4060"], refname25);
          
          //d10-x01-y05: 60-88% centrality
          string refname26 = mkAxisCode(10, 1, 5);
          const Scatter2D& refdata26 = refData(refname26);
          book(hEtaPt["ptyieldsdAuc6088c"], refname26 + "_pp_Eta", refdata26);
          book(hPionPt["ptyieldsdAuc6088"], refname26 + "_pp_Pion", refdata26);
          book(RatiodAu["EtaToPion6088"], refname26);

          //Figure 20: eta/pion^0 ratios in AuAu collisions
          //d11-x01-y01: 00-92% centrality (minimum bias)
          string refname27 = mkAxisCode(11, 1, 1);
          const Scatter2D& refdata27 = refData(refname27);
          book(hEtaPt["ptyieldsAuAuc0092c"], refname27 + "_pp_Eta", refdata27);
          book(hPionPt["ptyieldsAuAuc0092"], refname27 + "_pp_Pion", refdata27);
          book(RatioAuAu["EtaToPion0092"], refname27);
          
          //d11-x01-y02: 00-20% centrality
          string refname28 = mkAxisCode(11, 1, 2);
          const Scatter2D& refdata28 = refData(refname28);
          book(hEtaPt["ptyieldsAuAuc0020c"], refname28 + "_pp_Eta", refdata28);
          book(hPionPt["ptyieldsAuAuc0020"], refname28 + "_pp_Pion", refdata28);
          book(RatioAuAu["EtaToPion0020"], refname28);
          
          //d11-x01-y03: 20-60% centrality
          string refname29 = mkAxisCode(11, 1, 3);
          const Scatter2D& refdata29 = refData(refname29);
          book(hEtaPt["ptyieldsAuAuc2060c"], refname29 + "_pp_Eta", refdata29);
          book(hPionPt["ptyieldsAuAuc2060"], refname29 + "_pp_Pion", refdata29);
          book(RatioAuAu["EtaToPion2060"], refname29);
          
          //d11-x01-y04: 60-92% centrality
          string refname30 = mkAxisCode(11, 1, 4);
          const Scatter2D& refdata30 = refData(refname30);
          book(hEtaPt["ptyieldsAuAuc6092c"], refname30 + "_pp_Eta", refdata30);
          book(hPionPt["ptyieldsAuAuc6092"], refname30 + "_pp_Pion", refdata30);
          book(RatioAuAu["EtaToPion6092"], refname30); 
    }



    /// Perform the per-event analysis
      void analyze(const Event& event) {
          Particles chargedParticles = applyProjection<PrimaryParticles>(event,"cp").particles();
          Particles neutralParticles = applyProjection<UnstableParticles>(event,"np").particles();
          
          const ParticlePair& beam = beams();

          /*if (beamOpt == "NONE") {
          if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu200;
          else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = pp;
          else if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000791970) collSys = DAu200;
          else if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000010020) collSys = DAu200;
          }
          //check the collision system
          if (beamOpt == "PP200") collSys = pp;
          else if (beamOpt == "AUAU200") collSys = AuAu200;
          else if (beamOpt == "DAU200") collSys = DAu200;*/
          
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
                          hCrossSec["ppEtaa"]->fill(partPt, pt_weight);
                          hEtaPt["ptyieldsEta"]->fill(partPt);
                          hCrossSec["ppEtadAuc0088"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc0020"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc2040"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc4060"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc6088"]->fill(partPt, pt_weight);
                          break;
                      }
                  }
              }
              
              //a conditional for the neutral particles: pi^0
              for (Particle p : neutralParticles){
                  
                  
                  //define what we will fill our histograms with
                  //we technically don't have to do this, but it's convenient
                  double partPt = p.pT() / GeV;
                  //double pt_weight = 1. / (partPt * 2. * M_PI);
                
                  switch(p.pid()) {
                      case 111: { //pi0
                          hPionPt["ptyieldsPion"]->fill(partPt);
                          break;
                      }
                  }
                  break;
              }
              
          }
          
          
          
          //Nik, SEE HERE: down below is where you will do your edit for figure 13.1 and 13.2
          if (collSys == DAu200){
            
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
                          hCrossSec["dAuEta"]->fill(partPt, pt_weight);
                          hEtaPt["ptyieldsdAuc0088b"]->fill(partPt, pt_weight);
                          hEtaPt["ptyieldsdAuc0088c"]->fill(partPt);
                          break;
                      }
                  }
              }
              
              //a conditional for the neutral particles: pi^0
              for (Particle p : neutralParticles){
                  //define what we will fill our histograms with below
                  //we technically don't have to do this, but it's convenient
                  double partPt = p.pT() / GeV;
                  //double pt_weight = 1. / (partPt * 2. * M_PI);
                  
                  switch(p.pid()) {
                      case 111: { //pi0
                          hPionPt["ptyieldsdAuc0088"]->fill(partPt);
                          break;
                      }
                          
                  }
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
                              hEtaPt["ptyieldsdAuc0020a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc0020b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc0020c"]->fill(partPt);
                              break;
                          }
                      }
                  }
                  
                  for (Particle p : neutralParticles){
                      //define what we will fill our histograms with below
                      //we technically don't have to do this, but it's convenient
                      double partPt = p.pT() / GeV;
                      //double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      switch(p.pid()) {
                          case 111: { //pi0
                              hPionPt["ptyieldsdAuc0020"]->fill(partPt);
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
                              hEtaPt["ptyieldsdAuc2040a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc2040b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc2040c"]->fill(partPt);
                              break;
                          }
                      }
                  }
                  for (Particle p : neutralParticles){
                      //define what we will fill our histograms with below
                      //we technically don't have to do this, but it's convenient
                      double partPt = p.pT() / GeV;
                      //double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      switch(p.pid()) {
                          case 111: { //pi0
                              hPionPt["ptyieldsdAuc2040"]->fill(partPt);
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
                              hEtaPt["ptyieldsdAuc4060a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc4060b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc4060c"]->fill(partPt);
                              break;
                          }
                      }
                  }
                  for (Particle p : neutralParticles){
                      //define what we will fill our histograms with below
                      //we technically don't have to do this, but it's convenient
                      double partPt = p.pT() / GeV;
                      //double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      switch(p.pid()) {
                          case 111: { //pi0
                              hPionPt["ptyieldsdAuc4060"]->fill(partPt);
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
                              hEtaPt["ptyieldsdAuc6088a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc6088b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc6088c"]->fill(partPt);
                              break;
                          }
                      }
                  }
                  for (Particle p : neutralParticles){
                      //define what we will fill our histograms with below
                      //we technically don't have to do this, but it's convenient
                      double partPt = p.pT() / GeV;
                      //double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      switch(p.pid()) {
                          case 111: { //pi0
                              hPionPt["ptyieldsdAuc6088"]->fill(partPt);
                              break;
                          }
                              
                      }
                  }
              }

          } 

          if (collSys == AuAu200){
            
              const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
              const double c = cent();
              
              if ((c < 0.) || (c > 92.)) vetoEvent;
              
              sow["sow_AuAuc0092"]->fill();
              
              for (Particle p : chargedParticles) {
                  //comment these next two lines are only commented to avoid an unused variable warning
                  double partPt = p.pT() / GeV;
                  double pt_weight = 1. / (partPt * 2. * M_PI);
                  
                  switch (p.pid()) {
                      case 221: {  //eta
                          hEtaPt["ptyieldsAuAuc0092a"]->fill(partPt, pt_weight);
                          hEtaPt["ptyieldsAuAuc0092c"]->fill(partPt);
                          break;
                      }
                  }
              }
              
              //a conditional for the neutral particles: pi^0
              for (Particle p : neutralParticles){
                  //define what we will fill our histograms with below
                  //we technically don't have to do this, but it's convenient
                  double partPt = p.pT() / GeV;
                  //double pt_weight = 1. / (partPt * 2. * M_PI);
                  
                  switch(p.pid()) {
                      case 111: { //pi0
                          hPionPt["ptyieldsAuAuc0092"]->fill(partPt);
                          break;
                      }
                          
                  }
              }
              
              
              //Centrality inclusion begins here:
              if ((c >= 0.) && (c < 20.)){
                  for (Particle p : chargedParticles) {
                      
                      double partPt = p.pT() / GeV;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      sow["sow_AuAuc0020"]->fill();
                      
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hEtaPt["ptyieldsAuAuc0020a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc0020b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc0020c"]->fill(partPt);
                              break;
                          }
                      }
                  }
                  
                  for (Particle p : neutralParticles){
                      //define what we will fill our histograms with below
                      //we technically don't have to do this, but it's convenient
                      double partPt = p.pT() / GeV;
                      //double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      switch(p.pid()) {
                          case 111: { //pi0
                              hPionPt["ptyieldsAuAuc0020"]->fill(partPt);
                              break;
                          }
                              
                      }
                  }
              }
              else if ((c >= 20.) && (c < 60.)){
                  for (Particle p : chargedParticles) {
                      
                      double partPt = p.pT() / GeV;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      sow["sow_AuAuc2060"]->fill();
                      
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hEtaPt["ptyieldsAuAuc2060a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc2060b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc2060c"]->fill(partPt);
                              break;
                          }
                      }
                  }
                  for (Particle p : neutralParticles){
                      //define what we will fill our histograms with below
                      //we technically don't have to do this, but it's convenient
                      double partPt = p.pT() / GeV;
                      //double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      switch(p.pid()) {
                          case 111: { //pi0
                              hPionPt["ptyieldsAuAuc2060"]->fill(partPt);
                              break;
                          }
                              
                      }
                  }
              }
              else if ((c >= 60.) && (c < 92.)){
                  for (Particle p : chargedParticles) {
                      
                      double partPt = p.pT() / GeV;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      sow["sow_AuAuc6092"]->fill();
                      
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hEtaPt["ptyieldsAuAuc6092a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc6092b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc6092c"]->fill(partPt);
                              break;
                          }
                      }
                  }
                  for (Particle p : neutralParticles){
                      //define what we will fill our histograms with below
                      //we technically don't have to do this, but it's convenient
                      double partPt = p.pT() / GeV;
                      //double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      switch(p.pid()) {
                          case 111: { //pi0
                              hPionPt["ptyieldsAuAuc6092"]->fill(partPt);
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
          binShift(*hCrossSec["ppEtaa"]);
          hCrossSec["ppEtaa"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtaa"]->scaleW(cross);
          
          //d03-x01-y01
          binShift(*hCrossSec["dAuEta"]);
          hCrossSec["dAuEta"]->scaleW(1. / sow["sow_dAu"]->sumW());
          hCrossSec["dAuEta"]->scaleW(cross);          
          
          //Figure 14:
          //d05-x01-y01
          binShift(*hEtaPt["ptyieldsdAuc0020a"]);
          binShift(*hEtaPt["ptyieldsdAuc2040a"]);
          binShift(*hEtaPt["ptyieldsdAuc4060a"]);
          binShift(*hEtaPt["ptyieldsdAuc6088a"]);
          hEtaPt["ptyieldsdAuc0020a"]->scaleW(1. / sow["sow_dAuc0020"]->sumW());
          hEtaPt["ptyieldsdAuc2040a"]->scaleW(1. / sow["sow_dAuc2040"]->sumW());
          hEtaPt["ptyieldsdAuc4060a"]->scaleW(1. / sow["sow_dAuc4060"]->sumW());
          hEtaPt["ptyieldsdAuc6088a"]->scaleW(1. / sow["sow_dAuc6088"]->sumW());
          
          //Figure 15:
          //d05-x01-y01
          binShift(*hEtaPt["ptyieldsAuAuc0092a"]);
          binShift(*hEtaPt["ptyieldsAuAuc0020a"]);
          binShift(*hEtaPt["ptyieldsAuAuc2060a"]);
          binShift(*hEtaPt["ptyieldsAuAuc6092a"]);
          hEtaPt["ptyieldsAuAuc0092a"]->scaleW(1. / sow["sow_AuAuc0092"]->sumW());
          hEtaPt["ptyieldsAuAuc0020a"]->scaleW(1. / sow["sow_AuAuc0020"]->sumW());
          hEtaPt["ptyieldsAuAuc2060a"]->scaleW(1. / sow["sow_AuAuc2060"]->sumW());
          hEtaPt["ptyieldsAuAuc6092a"]->scaleW(1. / sow["sow_AuAuc6092"]->sumW());

          //Figure 16: Rda
          //d07-x01-y01:
          //denominator: must do our process to our cross section as done in Figure 13 plus multiplying it by <Tda>
          binShift(*hCrossSec["ppEtadAuc0088"]);
          binShift(*hEtaPt["ptyieldsdAuc0088b"]);
          hCrossSec["ppEtadAuc0088"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc0088"]->scaleW(cross);
          hCrossSec["ppEtadAuc0088"]->scaleW(0.2); //scaling by <TdA> as shown in Table II
          //numerator: must do our process to our invariant yield like in Figure 14
          hEtaPt["ptyieldsdAuc0088b"]->scaleW(1. / sow["sow_dAu"]->sumW());
          //Rda
          divide(hEtaPt["ptyieldsdAuc0088b"], hCrossSec["ppEtadAuc0088"], hRda["EtadAuc0088"]);
          
          //d07-x01-y02
          //denominator
          binShift(*hCrossSec["ppEtadAuc0020"]);
          hCrossSec["ppEtadAuc0020"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc0020"]->scaleW(cross);
          hCrossSec["ppEtadAuc0020"]->scaleW(0.36); //scaling by <TdA> as shown in Table II
          //numerator
          binShift(*hEtaPt["ptyieldsdAuc0020b"]);
          hEtaPt["ptyieldsdAuc0020b"]->scaleW(1. / sow["sow_dAuc0020"]->sumW());
          //Rda
          divide(hEtaPt["ptyieldsdAuc0020b"], hCrossSec["ppEtadAuc0020"], hRda["EtadAuc0020"]);
          
          //d07-x01-y03
          binShift(*hCrossSec["ppEtadAuc2040"]);
          hCrossSec["ppEtadAuc2040"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc2040"]->scaleW(cross);
          hCrossSec["ppEtadAuc2040"]->scaleW(0.25); //scaling by <TdA> as shown in Table II
          //numerator
          binShift(*hEtaPt["ptyieldsdAuc2040b"]);
          hEtaPt["ptyieldsdAuc2040b"]->scaleW(1. / sow["sow_dAuc2040"]->sumW());
          //Rda
          divide(hEtaPt["ptyieldsdAuc2040b"], hCrossSec["ppEtadAuc2040"], hRda["EtadAuc2040"]);
          
          //d07-x01-y04
          binShift(*hCrossSec["ppEtadAuc4060"]);
          hCrossSec["ppEtadAuc4060"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc4060"]->scaleW(cross);
          hCrossSec["ppEtadAuc4060"]->scaleW(0.17); //scaling by <TdA> as shown in Table II
          //numerator
          binShift(*hEtaPt["ptyieldsdAuc4060b"]);
          hEtaPt["ptyieldsdAuc4060b"]->scaleW(1. / sow["sow_dAuc4060"]->sumW());
          //Rda
          divide(hEtaPt["ptyieldsdAuc4060b"], hCrossSec["ppEtadAuc4060"], hRda["EtadAuc4060"]);
          
          //d07-x01-y05
          binShift(*hCrossSec["ppEtadAuc6088"]);
          hCrossSec["ppEtadAuc6088"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc6088"]->scaleW(cross);
          hCrossSec["ppEtadAuc6088"]->scaleW(0.073); //scaling by <TdA> as shown in Table II
          //numerator
          binShift(*hEtaPt["ptyieldsdAuc6088b"]);
          hEtaPt["ptyieldsdAuc6088b"]->scaleW(1. / sow["sow_dAuc6088"]->sumW());
          //Rda
          divide(hEtaPt["ptyieldsdAuc6088b"], hCrossSec["ppEtadAuc6088"],hRda["EtadAuc6088"]);

          //Figure 17: RAA  
          //d08-x01-y01:       
          binShift(*hCrossSec["ppEtaAuAuc0020"]);
          hCrossSec["ppEtaAuAuc0020"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtaAuAuc0020"]->scaleW(cross);
          hCrossSec["ppEtaAuAuc0020"]->scaleW(18.5); //scaling by <TdA> as shown in Table II
          //numerator: must do our process to our invariant yield like in Figure 14
          binShift(*hEtaPt["ptyieldsAuAuc0020b"]);
          hEtaPt["ptyieldsAuAuc0020b"]->scaleW(1. / sow["sow_AuAuc0020"]->sumW());
          //RAA
          divide(hEtaPt["ptyieldsAuAuc0020b"], hCrossSec["ppEtaAuAuc0020"], hRaa["EtaAuAuc0020"]);
          
          //d08-x01-y02
          //denominator
          binShift(*hCrossSec["ppEtaAuAuc2060"]);
          hCrossSec["ppEtaAuAuc2060"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtaAuAuc2060"]->scaleW(cross);
          hCrossSec["ppEtaAuAuc2060"]->scaleW(4.6); //scaling by <TdA> as shown in Table II
          //numerator
          binShift(*hEtaPt["ptyieldsAuAuc2060b"]);
          hEtaPt["ptyieldsAuAuc2060b"]->scaleW(1. / sow["sow_AuAuc2060"]->sumW());
          //RAA
          divide(hEtaPt["ptyieldsAuAuc2060b"], hCrossSec["ppEtaAuAuc2060"], hRaa["EtaAuAuc2060"]);
          
          //d08-x01-y03
          binShift(*hCrossSec["ppEtaAuAuc6092"]);
          hCrossSec["ppEtaAuAuc6092"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtaAuAuc6092"]->scaleW(cross);
          hCrossSec["ppEtaAuAuc6092"]->scaleW(0.3); //scaling by <TdA> as shown in Table II
          //numerator
          binShift(*hEtaPt["ptyieldsAuAuc6092b"]);
          hEtaPt["ptyieldsAuAuc6092b"]->scaleW(1. / sow["sow_AuAuc6092"]->sumW());
          //RAA
          divide(hEtaPt["ptyieldsAuAuc6092b"], hCrossSec["ppEtaAuAuc6092"], hRaa["EtaAuAuc6092"]);

          //Figure 18: ratio
          //d09-x01-y01
          binShift(*hEtaPt["ptyieldsEta"]);
          binShift(*hPionPt["ptyieldsPion"]);
          divide(hEtaPt["ptyieldsEta"], hPionPt["ptyieldsPion"], Ratiopp["EtaToPion"]);
          
          
          //Figure 19: ratios
          //d10-x01-y01
          binShift(*hEtaPt["ptyieldsdAuc0088c"]);
          binShift(*hPionPt["ptyieldsdAuc0088"]);
          divide(hEtaPt["ptyieldsdAuc0088c"], hPionPt["ptyieldsdAuc0088"], RatiodAu["EtaToPion0088"]);
          
          //d10-x01-y02
          binShift(*hEtaPt["ptyieldsdAuc0020c"]);
          binShift(*hPionPt["ptyieldsdAuc0020"]);
          divide(hEtaPt["ptyieldsdAuc0020c"], hPionPt["ptyieldsdAuc0020"], RatiodAu["EtaToPion0020"]);
    
          //d10-x01-y03
          binShift(*hEtaPt["ptyieldsdAuc2040c"]);
          binShift(*hPionPt["ptyieldsdAuc2040"]);
          divide(hEtaPt["ptyieldsdAuc2040c"], hPionPt["ptyieldsdAuc2040"], RatiodAu["EtaToPion2040"]);
          
          //d10-x01-y04
          binShift(*hEtaPt["ptyieldsdAuc4060c"]);
          binShift(*hPionPt["ptyieldsdAuc4060"]);
          divide(hEtaPt["ptyieldsdAuc4060c"], hPionPt["ptyieldsdAuc4060"], RatiodAu["EtaToPion4060"]);
          
          //d10-x01-y05
          binShift(*hEtaPt["ptyieldsdAuc6088c"]);
          binShift(*hPionPt["ptyieldsdAuc6088"]);
          divide(hEtaPt["ptyieldsdAuc6088c"], hPionPt["ptyieldsdAuc6088"], RatiodAu["EtaToPion6088"]);

          //Figure 20: ratios
          //d11-x01-y01
          binShift(*hEtaPt["ptyieldsAuAuc0092c"]);
          binShift(*hPionPt["ptyieldsAuAuc0092"]);
          divide(hEtaPt["ptyieldsAuAuc0092c"], hPionPt["ptyieldsAuAuc0092"], RatioAuAu["EtaToPion0092"]);
          
          //d11-x01-y02
          binShift(*hEtaPt["ptyieldsAuAuc0020c"]);
          binShift(*hPionPt["ptyieldsAuAuc0020"]);
          divide(hEtaPt["ptyieldsAuAuc0020c"], hPionPt["ptyieldsAuAuc0020"], RatioAuAu["EtaToPion0020"]);
    
          //d11-x01-y03
          binShift(*hEtaPt["ptyieldsAuAuc2060c"]);
          binShift(*hPionPt["ptyieldsAuAuc2060"]);
          divide(hEtaPt["ptyieldsAuAuc2060c"], hPionPt["ptyieldsAuAuc2060"], RatioAuAu["EtaToPion2060"]);

          //d11-x01-y04
          binShift(*hEtaPt["ptyieldsAuAuc6092c"]);
          binShift(*hPionPt["ptyieldsAuAuc6092"]);
          divide(hEtaPt["ptyieldsAuAuc6092c"], hPionPt["ptyieldsAuAuc6092"], RatioAuAu["EtaToPion6092"]);
          
      }

      //histograms
      map<string, Histo1DPtr> hCrossSec;
      map<string, Histo1DPtr> hPionPt;
      map<string, Histo1DPtr> hEtaPt;
      
      //ratios
      map<string, Scatter2DPtr> RatioAuAu;
      map<string, Scatter2DPtr> RatiodAu;
      map<string, Scatter2DPtr> Ratiopp;
      
      //RdA, Raa
      map<string, Scatter2DPtr> hRda;
      map<string, Scatter2DPtr> hRaa;
      
      //Counters
      map<string, CounterPtr> sow;
      
      
      string beamOpt;
      enum CollisionSystem { pp, AuAu200, DAu200 };
      CollisionSystem collSys;
      
      //these two vectors are only initialized if we would like to for loop over centrality bins
      //vector<int> AuAuCentralityBins{ 20, 60, 92 };
      //vector<int> dAuCentralityBins{ 20, 40, 60, 88 };

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2007_I731133);

}
