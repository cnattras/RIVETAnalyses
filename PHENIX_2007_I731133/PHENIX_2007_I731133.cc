// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/UnstableParticles.hh"
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
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2007_I731133);

    //create binShift function
void binShift(YODA::Histo1D& histogram) {

    const auto& binlist = histogram.bins();

    for (size_t n = 1; n < binlist.size()-1; ++n) {

        const auto& bins = binlist[n];

        double p_high = bins.xMax();
        double p_low  = bins.xMin();

        double left;
        double right;

        // Define neighboring bins
        if (n == 1) {
            left  = binlist[1].sumW();
            right = binlist[2].sumW();

        } else if (n == binlist.size()-2) {
            left  = binlist[binlist.size()-3].sumW();
            right = binlist[binlist.size()-2].sumW();

        } else {
            left  = binlist[n-1].sumW();
            right = binlist[n+1].sumW();
        }

        // Protect logarithm
        if (left <= 0 || right <= 0) {
            continue;
        }

        double b = 1.0 / (p_high - p_low) * log(left / right);

        double denom =
            pow(M_E, -b * p_high) -
            pow(M_E, -b * p_low);

        // Protect division by zero
        if (std::abs(denom) < 1e-12) {
            continue;
        }

        double f_corr =
            -b * (p_high - p_low)
            * pow(M_E, -b * (p_high+p_low) / 2)
            / denom;

        histogram.bin(n).scaleW(f_corr);
    }
}

    /// Book histograms and initialise projections before the run
      void init() {
          //Particles: eta (respectively)
          
          //declare cuts; most of these are found in section II of the paper
          //For charged particles:
          //consider adding: && Cuts::pT < 12*GeV, && Cuts::abscharge == 0, && Cuts::abscharge > 0 Cuts::abseta < 2.2 && Cuts::abseta > 1.2 &&
          const ALICE::PrimaryParticles cp(Cuts::absrap < 0.5 && Cuts::pT > 1*GeV);
          declare(cp, "cp");
          
          //Uncharged particles
          const UnstableParticles np((Cuts::abspid == 111 || Cuts::abspid == 221) && Cuts::absrap < 0.5 && Cuts::pT > 1*GeV);
          declare(np, "np");
          
          beamOpt = getOption<string>("beam", "NONE");
          fixedcentralityOpt = getOption<string>("fixedcentrality", "NONE");

          const ParticlePair& beam = beams();
          if (beamOpt == "NONE") {
        if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
          {
            float NN = 197.;
            if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1e-3)) collSys = AuAu200;
          }  
          else if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
          {
            float NN = 1.;
            if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1e-3)){ collSys = pp;}
          }
          else if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000791970)
          {
            // checking energy using form sqrt{s_{NN}}/(num_nucleons1*num_nucleons2)
            if (fuzzyEquals(sqrtS()/GeV, 3974., 1e-3)) collSys = DAu200;
          }
          else if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000010020)
          {
            if (fuzzyEquals(sqrtS()/GeV, 3974. , 1e-3)) collSys = DAu200;
          }
          }
          //check the collision system
          if (beamOpt == "PP200") collSys = pp;
          else if (beamOpt == "AUAU200") collSys = AuAu200;
          else if (beamOpt == "DAU200") collSys = DAu200;

/*          auto printC = [&]() {
          cerr << "sqrtS = " << sqrtS()/GeV
           << " GeV, collSys = "
          << (collSys == pp ? "pp" :
           collSys == AuAu200 ? "AuAu200" :
           collSys == DAu200 ? "DAu200" :
           "UNKNOWN")
          << ", fixed centrality = "
          << fixedcentralityOpt
          << endl;
          };
            printC();
*/          
          //declaration for collision systems that are not p+p
          declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
          if(fixedcentralityOpt!= "NONE"){
	        manualCentrality = std::stof(fixedcentralityOpt);}

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
         
          //figure 13 (12.1 in hepdata) pp-eta_pt_bins
          //d01-x01-y01
          string refname1 = mkAxisCode(1, 1, 1);
          //mkAxisCode gives us the internal histogram name for a given d, x, and y
          //const Estimate1D& refdata1 = refData(refname1);
          //here we have to define a Estimate1D& for the next part, we can use the refData function on refname
          book(hCrossSec["ppEtaa"], 1, 1, 1);
          //here, we are using book() as: book(1DPtr&, const string &name, const Estimate1D), and this books a histograms with binning using d01-x01-y01 from our yoda as a reference.

          //d02-x01-y01 We don't worry about the decay channels, so we won't need this one since it is a duplicate
          //string refname2 = mkAxisCode(2, 1, 1);
          //const Estimate1D& refdata2 = refData(refname2);
          //book(hCrossSec["ppEtatoPiona"], refname2 + "_pp_Pion", refdata2);

          //d03-x01-y01 (hepdata fig 13.1) dAu-to-eta_pt_bins
          string refname2 = mkAxisCode(3, 1, 1);
          //mkAxisCode gives us the internal histogram name for a given d, x, and y
          //const Estimate1D& refdata2 = refData(refname2);
          //here we have to define a Estimate1D& for the next part, we can use the refData function on refname
          book(hCrossSec["dAuEta"], 3, 1, 1);

          //Figure 14.1: Invariant yields of eta in d+Au collisions
          //d05-x01-y01: 00-20% centrality
          string refname5 = mkAxisCode(5, 1, 1);
          //const Estimate1D& refdata5 = refData(refname5);
          book(hEtaPt["ptyieldsdAuc0020a"], 5, 1, 1);

          //d05-x01-y02: 20-40% centrality
          string refname6 = mkAxisCode(5, 1, 2);
          //const Estimate1D& refdata6 = refData(refname6);
          book(hEtaPt["ptyieldsdAuc2040a"], 5, 1, 2);

          //d05-x01-y03: 40-60% centrality
          string refname7 = mkAxisCode(5, 1, 3);
          //const Estimate1D& refdata7 = refData(refname7);
          book(hEtaPt["ptyieldsdAuc4060a"], 5, 1, 3);

          //d06-x01-y01: (hepdata fig 14.2) 60-88% centrality
          string refname8 = mkAxisCode(6, 1, 1);
          //const Estimate1D& refdata8 = refData(refname8);
          book(hEtaPt["ptyieldsdAuc6088a"], 6, 1, 1);

          //Figure 15 AuAu
          //d07-x01-y01 min. bias
          string refname9 = mkAxisCode(7, 1, 1);
          //const Estimate1D& refdata9 = refDa6, 1, 1
          book(hEtaPt["ptyieldsAuAuc0092a"], 7, 1, 1);         

          //d07-x01-y02 0-20%
          string refname10 = mkAxisCode(7, 1, 2);
          //const Estimate1D& refdata10 = refData(refname10);
          book(hEtaPt["ptyieldsAuAuc0020a"], 7, 1, 2);         

          //d07-x01-y03 20-40%
          string refname11 = mkAxisCode(7, 1, 3);
          //const Estimate1D& refdata11 = refData(refname11);
          book(hEtaPt["ptyieldsAuAuc2060a"], 7, 1, 3);          

          //d08-x01-y01 60-92%
          string refname12 = mkAxisCode(8, 1, 1);
          //const Estimate1D& refdata12 = refData(refname12);
          book(hEtaPt["ptyieldsAuAuc6092a"], 8, 1, 1);         

          //Figure 16.1: Rda for measured etas
          //d09-x01-y01: 00-88% centrality (minimum bias)
          string refname13 = mkAxisCode(9, 1, 1);
          const Estimate1D& refdata13 = refData(refname13);
          book(hEtaPt["ptyieldsdAuc0088b"],  /* "_" + */ refname13 + "_dAuc0088_Eta", refdata13);
          book(hCrossSec["ppEtadAuc0088"],  /* "_" + */ refname13 + "_pp_Eta", refdata13);
          book(hRda["EtadAuc0088"], 9, 1, 1);

          //d09-x01-y02: 00-20% centrality
          string refname14 = mkAxisCode(9, 1, 2);
          const Estimate1D& refdata14 = refData(refname14);
          book(hEtaPt["ptyieldsdAuc0020b"],  /* "_" + */ refname14 + "_dAuc0020_Eta", refdata14);
          book(hCrossSec["ppEtadAuc0020"],  /* "_" + */ refname14 + "_pp_Eta", refdata14);
          book(hRda["EtadAuc0020"], refname14);

          //d09-x01-y03: 19-40% centrality
          string refname15 = mkAxisCode(9, 1, 3);
          const Estimate1D& refdata15 = refData(refname15);
          book(hEtaPt["ptyieldsdAuc2040b"],  /* "_" + */ refname15 + "_dAuc2040_Eta", refdata15);
          book(hCrossSec["ppEtadAuc2040"],  /* "_" + */ refname15 + "_pp_Eta", refdata15);
          book(hRda["EtadAuc2040"], refname15);

          //d09-x01-y04: 40-60% centrality
          string refname16 = mkAxisCode(9, 1, 4);
          const Estimate1D& refdata16 = refData(refname16);
          book(hEtaPt["ptyieldsdAuc4060b"],  /* "_" + */ refname16 + "_dAuc4060_Eta", refdata16);
          book(hCrossSec["ppEtadAuc4060"],  /* "_" + */ refname16 + "_pp_Eta", refdata16);
          book(hRda["EtadAuc4060"], refname16);

          //d10-x01-y01: Figure 16.2, 60-88% centrality
          string refname17 = mkAxisCode(10, 1, 1);
          const Estimate1D& refdata17 = refData(refname17);
          book(hEtaPt["ptyieldsdAuc6088b"],  /* "_" + */ refname17 + "_dAuc6088_Eta", refdata17);
          book(hCrossSec["ppEtadAuc6088"],  /* "_" + */ refname17 + "_pp_Eta", refdata17);
          book(hRda["EtadAuc6088"], refname17);

          //Figure 17.1
          //d11-x01-y01: 00-20% centrality
          string refname18 = mkAxisCode(11, 1, 1);
          const Estimate1D& refdata18 = refData(refname18);
          book(hEtaPt["ptyieldsAuAuc0020b"],  /* "_" + */ refname18 + "_AuAuc0020_Eta", refdata18);
          book(hCrossSec["ppEtaAuAuc0020"],  /* "_" + */ refname18 + "_pp_Eta", refdata18);
          book(hRaa["EtaAuAuc0020"], refname18);

          //d11-x01-y02: 20-60% centrality
          string refname19 = mkAxisCode(11, 1, 2);
          const Estimate1D& refdata19 = refData(refname19);
          book(hEtaPt["ptyieldsAuAuc2060b"],  /* "_" + */ refname19 + "_AuAuc2060_Eta", refdata19);
          book(hCrossSec["ppEtaAuAuc2060"],  /* "_" + */ refname19 + "_pp_Eta", refdata19);
          book(hRaa["EtaAuAuc2060"], refname19);

          //d12-x01-y01: Figure 17.2, 60-92% centrality
          string refname20 = mkAxisCode(12, 1, 1);
          const Estimate1D& refdata20 = refData(refname20);
          book(hEtaPt["ptyieldsAuAuc6092b"],  /* "_" + */ refname20 + "_AuAuc6092_Eta", refdata20);
          book(hCrossSec["ppEtaAuAuc6092"],  /* "_" + */ refname20 + "_pp_Eta", refdata20);
          book(hRaa["EtaAuAuc6092"], refname20);    

          //Figure 18: eta/pion^0 ratio in pp collisions
          //d13-x01-y01
          string refname21 = mkAxisCode(13, 1, 1);
          const Estimate1D& refdata21 = refData(refname21);
          book(hEtaPt["ptyieldsEta"],  /* "_" + */ refname21 + "_pp_Eta", refdata21);
          book(hPionPt["ptyieldsPion"],  /* "_" + */ refname21 + "_pp_Pion", refdata21);
          book(Ratiopp["EtaToPion"], 13, 1, 1);


          //Figure 19.1: eta/pion^0 ratios in dAu collisions
          //d14-x01-y01: 00-88% centrality (minimum bias)
          string refname22 = mkAxisCode(14, 1, 1);
          const Estimate1D& refdata22 = refData(refname22);
          book(hEtaPt["ptyieldsdAuc0088c"], /* "_" + */  refname22 + "_pp_Eta", refdata22);
          book(hPionPt["ptyieldsdAuc0088"],  /* "_" + */ refname22 + "_pp_Pion", refdata22);
          book(RatiodAu["EtaToPion0088"], refname22);

          //d14-x01-y02: 00-20% centrality
          string refname23 = mkAxisCode(14, 1, 2);
          const Estimate1D& refdata23 = refData(refname23);
          book(hEtaPt["ptyieldsdAuc0020c"],  /* "_" + */ refname23 + "_pp_Eta", refdata23);
          book(hPionPt["ptyieldsdAuc0020"],  /* "_" + */ refname23 + "_pp_Pion", refdata23);
          book(RatiodAu["EtaToPion0020"], refname23);

          //d14-x01-y03: 20-40% centrality
          string refname24 = mkAxisCode(14, 1, 3);
          const Estimate1D& refdata24 = refData(refname24);
          book(hEtaPt["ptyieldsdAuc2040c"],  /* "_" + */ refname24 + "_pp_Eta", refdata24);
          book(hPionPt["ptyieldsdAuc2040"],  /* "_" + */ refname24 + "_pp_Pion", refdata24);
          book(RatiodAu["EtaToPion2040"], refname24);

          //d14-x01-y04: 40-60% centrality
          string refname25 = mkAxisCode(14, 1, 4);
          const Estimate1D& refdata25 = refData(refname25);
          book(hEtaPt["ptyieldsdAuc4060c"],  /* "_" + */ refname25 + "_pp_Eta", refdata25);
          book(hPionPt["ptyieldsdAuc4060"],  /* "_" + */ refname25 + "_pp_Pion", refdata25);
          book(RatiodAu["EtaToPion4060"], refname25);

          //d15-x01-y01: Figure 19.2 60-88% centrality
          string refname26 = mkAxisCode(15, 1, 1);
          const Estimate1D& refdata26 = refData(refname26);
          book(hEtaPt["ptyieldsdAuc6088c"],  /* "_" + */ refname26 + "_pp_Eta", refdata26);
          book(hPionPt["ptyieldsdAuc6088"],  /* "_" + */ refname26 + "_pp_Pion", refdata26);
          book(RatiodAu["EtaToPion6088"], refname26);

          //Figure 20.1: eta/pion^0 ratios in AuAu collisions
          //d16-x01-y01: 00-92% centrality (minimum bias)
          string refname27 = mkAxisCode(16, 1, 1);
          const Estimate1D& refdata27 = refData(refname27);
          book(hEtaPt["ptyieldsAuAuc0092c"],  /* "_" + */ refname27 + "_pp_Eta", refdata27);
          book(hPionPt["ptyieldsAuAuc0092"],  /* "_" + */ refname27 + "_pp_Pion", refdata27);
          book(RatioAuAu["EtaToPion0092"], refname27);

          //d16-x01-y02: 00-20% centrality
          string refname28 = mkAxisCode(16, 1, 2);
          const Estimate1D& refdata28 = refData(refname28);
          book(hEtaPt["ptyieldsAuAuc0020c"],  /* "_" + */ refname28 + "_pp_Eta", refdata28);
          book(hPionPt["ptyieldsAuAuc0020"],  /* "_" + */ refname28 + "_pp_Pion", refdata28);
          book(RatioAuAu["EtaToPion0020"], refname28);

          //d16-x01-y03: 20-60% centrality
          string refname29 = mkAxisCode(16, 1, 3);
          const Estimate1D& refdata29 = refData(refname29);
          book(hEtaPt["ptyieldsAuAuc2060c"],  /* "_" + */ refname29 + "_pp_Eta", refdata29);
          book(hPionPt["ptyieldsAuAuc2060"],  /* "_" + */ refname29 + "_pp_Pion", refdata29);
          book(RatioAuAu["EtaToPion2060"], refname29);

          //d17-x01-y01: (Figure 20.2) 60-92% centrality
          string refname30 = mkAxisCode(17, 1, 1);
          const Estimate1D& refdata30 = refData(refname30);
          book(hEtaPt["ptyieldsAuAuc6092c"],  /* "_" + */ refname30 + "_pp_Eta", refdata30);
          book(hPionPt["ptyieldsAuAuc6092"],  /* "_" + */ refname30 + "_pp_Pion", refdata30);
          book(RatioAuAu["EtaToPion6092"], refname30); 
        }

    /// Perform the per-event analysis
      void analyze(const Event& event) {
        //Particles chargedParticles = apply<PrimaryParticles>(event,"cp").particles();
          Particles neutralParticles = apply<UnstableParticles>(event,"np").particles();
          
          /*const ParticlePair& beam = beams();

          if (beamOpt == "NONE") {
          if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu200;
          else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = pp;
          else if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000791970) collSys = DAu200;
          else if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000010020) collSys = DAu200;
          }
          //check the collision system
          if (beamOpt == "PP200") collSys = pp;
          else if (beamOpt == "AUAU200") collSys = AuAu200;
          else if (beamOpt == "DAU200") collSys = DAu200;*/
          
          if (collSys == pp){    
              //a conditional for the charged particles: eta, pi+, pi-, gamma
              for (Particle p : neutralParticles) {
                  
                  
                  //define what we will fill our histograms with
                  //we technically don't have to do this, but it's convenient
                  double partPt = p.pT() / GeV;
                  if (!std::isfinite(partPt) || partPt <= 0.) continue;
                  double pt_weight = 1. / (partPt * 2. * M_PI);
                  
                  switch (p.pid()) {
                          //particle ids for convenience: {221 is eta, 111 is pi0};
                      case 221: {  //eta
                          hCrossSec["ppEtaa"]->fill(partPt, pt_weight);
                          hEtaPt["ptyieldsEta"]->fill(partPt);
                          hCrossSec["ppEtadAuc0088"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc0020"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc2040"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc4060"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtadAuc6088"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtaAuAuc0020"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtaAuAuc2060"]->fill(partPt, pt_weight);
                          hCrossSec["ppEtaAuAuc6092"]->fill(partPt, pt_weight);
                          break;
                      }
                      case 111: { //pi0
                          hPionPt["ptyieldsPion"]->fill(partPt);
                          break;
                      }
                  }
                  }
              }                     
          
          
          //Nik, SEE HERE: down below is where you will do your edit for figure 13.1 and 13.2
          if (collSys == DAu200){
            
              const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
              const double c = cent();
              
              if ((c < 0.) || (c > 88.)) vetoEvent;
              
              sow["sow_dAu"]->fill();
              
              for (Particle p : neutralParticles) {
                  //comment these next two lines are only commented to avoid an unused variable warning
                  double partPt = p.pT() / GeV;
                  if (!std::isfinite(partPt) || partPt <= 0.) continue;
                  double pt_weight = 1. / (partPt * 2. * M_PI);
                  
                  switch (p.pid()) {
                      case 221: {  //eta
                          hCrossSec["dAuEta"]->fill(partPt, pt_weight);
                          hEtaPt["ptyieldsdAuc0088b"]->fill(partPt, pt_weight);
                          hEtaPt["ptyieldsdAuc0088c"]->fill(partPt);
                          break;
                      }
                      case 111: { //pi0
                          hPionPt["ptyieldsdAuc0088"]->fill(partPt);
                          break;
                      }
                  }
              }              
              
              //Centrality inclusion begins here:
              if ((c >= 0.) && (c < 20.)){
                  
                  //fill our counter for this centrality
                  sow["sow_dAuc0020"]->fill(); //This was double filling the sow counter as it is below as well
                  
                  
                  for (Particle p : neutralParticles) {
                      
                      double partPt = p.pT() / GeV;
                      if (!std::isfinite(partPt) || partPt <= 0.) continue;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hEtaPt["ptyieldsdAuc0020a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc0020b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc0020c"]->fill(partPt);
                              break;
                          }
                          case 111: { //pi0
                              hPionPt["ptyieldsdAuc0020"]->fill(partPt);
                              break;
                          }
                      }
                  }
              }
              else if ((c >= 20.) && (c < 40.)){
                
                sow["sow_dAuc2040"]->fill();
                
                for (Particle p : neutralParticles) {
                      
                      double partPt = p.pT() / GeV;
                      if (!std::isfinite(partPt) || partPt <= 0.) continue;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hEtaPt["ptyieldsdAuc2040a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc2040b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc2040c"]->fill(partPt);
                              break;
                          }
                          case 111: { //pi0
                              hPionPt["ptyieldsdAuc2040"]->fill(partPt);
                              break;
                          }
                      }
                  }
              }
              else if ((c >= 40.) && (c < 60.)){
                  
                sow["sow_dAuc4060"]->fill();

                for (Particle p : neutralParticles) {
                      
                      double partPt = p.pT() / GeV;
                      if (!std::isfinite(partPt) || partPt <= 0.) continue;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                                            
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hEtaPt["ptyieldsdAuc4060a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc4060b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc4060c"]->fill(partPt);
                              break;
                          }
                          case 111: { //pi0
                              hPionPt["ptyieldsdAuc4060"]->fill(partPt);
                              break;
                          }
                      }
                  }
              }
              else if ((c >= 60.) && (c < 88.)){
                  
                sow["sow_dAuc6088"]->fill();
                
                for (Particle p : neutralParticles) {
                      
                      double partPt = p.pT() / GeV;
                      if (!std::isfinite(partPt) || partPt <= 0.) continue;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hEtaPt["ptyieldsdAuc6088a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc6088b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsdAuc6088c"]->fill(partPt);
                              break;
                          }
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
              
              if (c > 92.) {vetoEvent;}
              
              sow["sow_AuAuc0092"]->fill();
              
              for (Particle p : neutralParticles) {
                  //comment these next two lines are only commented to avoid an unused variable warning
                  double partPt = p.pT() / GeV;
                  if (!std::isfinite(partPt) || partPt <= 0.) continue;
                  double pt_weight = 1. / (partPt * 2. * M_PI);
                  
                  switch (p.pid()) {
                      case 221: {  //eta
                          hEtaPt["ptyieldsAuAuc0092a"]->fill(partPt, pt_weight);
                          hEtaPt["ptyieldsAuAuc0092c"]->fill(partPt);
                          break;
                      }
                      case 111: { //pi0
                          hPionPt["ptyieldsAuAuc0092"]->fill(partPt);
                          break;
                      }
                  }
              }   
              
              //Centrality inclusion begins here:
              if ((c >= 0.) && (c < 20.)){
                  
                sow["sow_AuAuc0020"]->fill();
                
                for (Particle p : neutralParticles) {
                      
                      double partPt = p.pT() / GeV;
                      if (!std::isfinite(partPt) || partPt <= 0.) continue;
                      double pt_weight = 1. / (partPt * 2. * M_PI);

                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hEtaPt["ptyieldsAuAuc0020a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc0020b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc0020c"]->fill(partPt);
                              break;
                          }
                          case 111: { //pi0
                              hPionPt["ptyieldsAuAuc0020"]->fill(partPt);
                              break;
                          }
                      }
                  }
              }
              else if ((c >= 20.) && (c < 60.)){
                  
                sow["sow_AuAuc2060"]->fill();
                
                for (Particle p : neutralParticles) {
                      
                      double partPt = p.pT() / GeV;
                      if (!std::isfinite(partPt) || partPt <= 0.) continue;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hEtaPt["ptyieldsAuAuc2060a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc2060b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc2060c"]->fill(partPt);
                              break;
                          }
                          case 111: { //pi0
                              hPionPt["ptyieldsAuAuc2060"]->fill(partPt);
                              break;
                          }
                      }
                  }
              }
              else if ((c >= 60.) && (c < 92.)){
                  
                sow["sow_AuAuc6092"]->fill();
                
                for (Particle p : neutralParticles) {
                      
                      double partPt = p.pT() / GeV;
                      if (!std::isfinite(partPt) || partPt <= 0.) continue;
                      double pt_weight = 1. / (partPt * 2. * M_PI);
                      //switch statement for each charged particle
                      switch (p.pid()) {
                          case 221: {  //eta
                              hEtaPt["ptyieldsAuAuc6092a"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc6092b"]->fill(partPt, pt_weight);
                              hEtaPt["ptyieldsAuAuc6092c"]->fill(partPt);
                              break;
                          }
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
      double cross = crossSection() / millibarn;


          //Figure 13:
          if(sow["sow_pp"]->sumW()>0){
          //d01-x01-y01
          binShift(*hCrossSec["ppEtaa"]);
          hCrossSec["ppEtaa"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtaa"]->scaleW(cross);
        }
          
          if(sow["sow_dAu"]->sumW()>0){
          //d03-x01-y01
          binShift(*hCrossSec["dAuEta"]);
          hCrossSec["dAuEta"]->scaleW(1. / sow["sow_dAu"]->sumW());
          hCrossSec["dAuEta"]->scaleW(cross);  
          }        
          
          //Figure 14.1:
          //d05-x01-y01
          binShift(*hEtaPt["ptyieldsdAuc0020a"]);
          if(sow["sow_dAuc0020"]->sumW()>0){
          hEtaPt["ptyieldsdAuc0020a"]->scaleW(1. / sow["sow_dAuc0020"]->sumW());}
          //d05-x01-y02
          binShift(*hEtaPt["ptyieldsdAuc2040a"]);
          if(sow["sow_dAuc2040"]->sumW()>0){
          hEtaPt["ptyieldsdAuc2040a"]->scaleW(1. / sow["sow_dAuc2040"]->sumW());}
          //d05-x01-y03
          binShift(*hEtaPt["ptyieldsdAuc4060a"]);
          if(sow["sow_dAuc4060"]->sumW()>0){
          hEtaPt["ptyieldsdAuc4060a"]->scaleW(1. / sow["sow_dAuc4060"]->sumW());}
          //d06-x01-y01
          binShift(*hEtaPt["ptyieldsdAuc6088a"]);
          if(sow["sow_dAuc6088"]->sumW()>0){
          hEtaPt["ptyieldsdAuc6088a"]->scaleW(1. / sow["sow_dAuc6088"]->sumW());}
          
          //Figure 15:
          //d07-x01-y01
          binShift(*hEtaPt["ptyieldsAuAuc0092a"]);
          if(sow["sow_AuAuc0092"]->sumW()>0){
          hEtaPt["ptyieldsAuAuc0092a"]->scaleW(1. / sow["sow_AuAuc0092"]->sumW());}
          //d07-x01-y02
          binShift(*hEtaPt["ptyieldsAuAuc0020a"]);
          if(sow["sow_AuAuc0020"]->sumW()>0){
          hEtaPt["ptyieldsAuAuc0020a"]->scaleW(1. / sow["sow_AuAuc0020"]->sumW());}
          //d07-x01-y03
          binShift(*hEtaPt["ptyieldsAuAuc2060a"]);
          if(sow["sow_AuAuc2060"]->sumW()>0){
          hEtaPt["ptyieldsAuAuc2060a"]->scaleW(1. / sow["sow_AuAuc2060"]->sumW());}
          //d08-x01-y01
          binShift(*hEtaPt["ptyieldsAuAuc6092a"]);
          if(sow["sow_AuAuc6092"]->sumW()>0){
          hEtaPt["ptyieldsAuAuc6092a"]->scaleW(1. / sow["sow_AuAuc6092"]->sumW());}

          //Figure 16: R_dA
        //d09-x01-y01:
          //denominator: must do our process to our cross section as done in Figure 13 plus multiplying it by <Tda>
          binShift(*hCrossSec["ppEtadAuc0088"]);
          binShift(*hEtaPt["ptyieldsdAuc0088b"]);
          if(sow["sow_pp"]->sumW()>0){
          hCrossSec["ppEtadAuc0088"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc0088"]->scaleW(cross);
          hCrossSec["ppEtadAuc0088"]->scaleW(0.2); //scaling by <TdA> as shown in Table II
        }
         if(sow["sow_dAu"]->sumW()>0){
          //numerator: must do our process to our invariant yield like in Figure 14
          hEtaPt["ptyieldsdAuc0088b"]->scaleW(1. / sow["sow_dAu"]->sumW());
        }
        if(sow["sow_dAu"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          //Rda
          divide(hEtaPt["ptyieldsdAuc0088b"], hCrossSec["ppEtadAuc0088"], hRda["EtadAuc0088"]);
        }
        //d09-x01-y02
          //denominator
          binShift(*hCrossSec["ppEtadAuc0020"]);
        if(sow["sow_pp"]->sumW()>0){
          hCrossSec["ppEtadAuc0020"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc0020"]->scaleW(cross);
          hCrossSec["ppEtadAuc0020"]->scaleW(0.36); //scaling by <TdA> as shown in Table II
        }  
          //numerator
        if(sow["sow_dAuc0020"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsdAuc0020b"]);
          hEtaPt["ptyieldsdAuc0020b"]->scaleW(1. / sow["sow_dAuc0020"]->sumW());
        }
          //Rda
        if(sow["sow_dAuc0020"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc0020b"], hCrossSec["ppEtadAuc0020"], hRda["EtadAuc0020"]);
        } 
        //d09-x01-y03
          //denominator
        if(sow["sow_pp"]->sumW()>0){
          binShift(*hCrossSec["ppEtadAuc2040"]);
          hCrossSec["ppEtadAuc2040"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc2040"]->scaleW(cross);
          hCrossSec["ppEtadAuc2040"]->scaleW(0.25); //scaling by <TdA> as shown in Table II
        }  
          //numerator
        if(sow["sow_dAuc2040"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsdAuc2040b"]);
          hEtaPt["ptyieldsdAuc2040b"]->scaleW(1. / sow["sow_dAuc2040"]->sumW());
        }  
          //Rda
        if(sow["sow_dAuc2040"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc2040b"], hCrossSec["ppEtadAuc2040"], hRda["EtadAuc2040"]);
        }  
        //d09-x01-y04
          //denominator
        if(sow["sow_pp"]->sumW()>0){
          binShift(*hCrossSec["ppEtadAuc4060"]);
          hCrossSec["ppEtadAuc4060"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc4060"]->scaleW(cross);
          hCrossSec["ppEtadAuc4060"]->scaleW(0.17); //scaling by <TdA> as shown in Table II
        }
          //numerator
        if(sow["sow_dAuc4060"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsdAuc4060b"]);
          hEtaPt["ptyieldsdAuc4060b"]->scaleW(1. / sow["sow_dAuc4060"]->sumW());
        }
          //Rda
        if(sow["sow_dAuc4060"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc4060b"], hCrossSec["ppEtadAuc4060"], hRda["EtadAuc4060"]);
        }        

        //d10-x01-y01
          //denominator
        if(sow["sow_pp"]->sumW()>0){
          binShift(*hCrossSec["ppEtadAuc6088"]);
          hCrossSec["ppEtadAuc6088"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc6088"]->scaleW(cross);
          hCrossSec["ppEtadAuc6088"]->scaleW(0.073); //scaling by <TdA> as shown in Table II
        }
          //numerator
        if(sow["sow_dAuc6088"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsdAuc6088b"]);
          hEtaPt["ptyieldsdAuc6088b"]->scaleW(1. / sow["sow_dAuc6088"]->sumW());
        }
          //Rda
        if(sow["sow_dAuc6088"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc6088b"], hCrossSec["ppEtadAuc6088"],hRda["EtadAuc6088"]);
        }

        //Figure 17.1: RAA  
          //d11-x01-y01:
          //denominator
        if (sow["sow_pp"]->sumW() > 0) {       
          binShift(*hCrossSec["ppEtaAuAuc0020"]);
          hCrossSec["ppEtaAuAuc0020"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtaAuAuc0020"]->scaleW(cross);
          hCrossSec["ppEtaAuAuc0020"]->scaleW(18.5); //scaling by <TdA> as shown in Table II
        }
          //numerator: must do our process to our invariant yield like in Figure 14
        if(sow["sow_AuAuc0020"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsAuAuc0020b"]);
          hEtaPt["ptyieldsAuAuc0020b"]->scaleW(1. / sow["sow_AuAuc0020"]->sumW());
        }
          //RAA
        if(sow["sow_AuAuc0020"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsAuAuc0020b"], hCrossSec["ppEtaAuAuc0020"], hRaa["EtaAuAuc0020"]);
        }

          //d11-x01-y02
          //denominator
        if (sow["sow_pp"]->sumW() > 0) {
          binShift(*hCrossSec["ppEtaAuAuc2060"]);
          hCrossSec["ppEtaAuAuc2060"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtaAuAuc2060"]->scaleW(cross);
          hCrossSec["ppEtaAuAuc2060"]->scaleW(4.6); //scaling by <TdA> as shown in Table II
        }
          //numerator
        if(sow["sow_AuAuc2060"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsAuAuc2060b"]);
          hEtaPt["ptyieldsAuAuc2060b"]->scaleW(1. / sow["sow_AuAuc2060"]->sumW());
        }
          //RAA
        if(sow["sow_AuAuc2060"]->sumW()>0 && sow["sow_pp"]->sumW()>0 && hCrossSec["ppEtaAuAuc2060"]->sumW() > 0){
          divide(hEtaPt["ptyieldsAuAuc2060b"], hCrossSec["ppEtaAuAuc2060"], hRaa["EtaAuAuc2060"]);
        }

        //d12-x01-y01
          //denominator
        if (sow["sow_pp"]->sumW() > 0) {
          binShift(*hCrossSec["ppEtaAuAuc6092"]);
          hCrossSec["ppEtaAuAuc6092"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtaAuAuc6092"]->scaleW(cross);
          hCrossSec["ppEtaAuAuc6092"]->scaleW(0.3); //scaling by <TdA> as shown in Table II
        }
          //numerator
        if(sow["sow_AuAuc6092"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsAuAuc6092b"]);
          hEtaPt["ptyieldsAuAuc6092b"]->scaleW(1. / sow["sow_AuAuc6092"]->sumW());
        }
          //RAA
        if(sow["sow_AuAuc6092"]->sumW()>0 && sow["sow_pp"]->sumW()>0 && hCrossSec["ppEtaAuAuc6092"]->sumW() > 0){
          divide(hEtaPt["ptyieldsAuAuc6092b"], hCrossSec["ppEtaAuAuc6092"], hRaa["EtaAuAuc6092"]);
        }

          //Figure 18: eta/pion^0 ratio in pp collisions
          //d13-x01-y01
          binShift(*hEtaPt["ptyieldsEta"]);
          binShift(*hPionPt["ptyieldsPion"]);
          if(hEtaPt["ptyieldsEta"]->sumW()>0 && hPionPt["ptyieldsPion"]->sumW()>0){
          divide(hEtaPt["ptyieldsEta"], hPionPt["ptyieldsPion"], Ratiopp["EtaToPion"]);}
        
          //Figure 19.1: eta/pion^0 ratios in dAu collisions
          //d14-x01-y01
          binShift(*hEtaPt["ptyieldsdAuc0088c"]);
          binShift(*hPionPt["ptyieldsdAuc0088"]);
          if(hEtaPt["ptyieldsdAuc0088c"]->sumW()>0 && hPionPt["ptyieldsdAuc0088"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc0088c"], hPionPt["ptyieldsdAuc0088"], RatiodAu["EtaToPion0088"]);}
          
          //d14-x01-y02
          binShift(*hEtaPt["ptyieldsdAuc0020c"]);
          binShift(*hPionPt["ptyieldsdAuc0020"]);
          if(hEtaPt["ptyieldsdAuc0020c"]->sumW()>0 && hPionPt["ptyieldsdAuc0020"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc0020c"], hPionPt["ptyieldsdAuc0020"], RatiodAu["EtaToPion0020"]);}
    
          //d14-x01-y03
          binShift(*hEtaPt["ptyieldsdAuc2040c"]);
          binShift(*hPionPt["ptyieldsdAuc2040"]);
          if(hEtaPt["ptyieldsdAuc2040c"]->sumW()>0 && hPionPt["ptyieldsdAuc2040"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc2040c"], hPionPt["ptyieldsdAuc2040"], RatiodAu["EtaToPion2040"]);}
          
          //d14-x01-y04
          binShift(*hEtaPt["ptyieldsdAuc4060c"]);
          binShift(*hPionPt["ptyieldsdAuc4060"]);
          if(hEtaPt["ptyieldsdAuc4060c"]->sumW()>0 && hPionPt["ptyieldsdAuc4060"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc4060c"], hPionPt["ptyieldsdAuc4060"], RatiodAu["EtaToPion4060"]);}
          
          //d15-x01-y01 (Figure 19.2)
          binShift(*hEtaPt["ptyieldsdAuc6088c"]);
          binShift(*hPionPt["ptyieldsdAuc6088"]);
          if(hEtaPt["ptyieldsdAuc6088c"]->sumW()>0 && hPionPt["ptyieldsdAuc6088"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc6088c"], hPionPt["ptyieldsdAuc6088"], RatiodAu["EtaToPion6088"]);}

          ///Figure 20.1: eta/pion^0 ratios in AuAu collisions
          //d16-x01-y01
          binShift(*hEtaPt["ptyieldsAuAuc0092c"]);
          binShift(*hPionPt["ptyieldsAuAuc0092"]);
          if(hEtaPt["ptyieldsAuAuc0092c"]->sumW()>0 && hPionPt["ptyieldsAuAuc0092"]->sumW()>0){
          divide(hEtaPt["ptyieldsAuAuc0092c"], hPionPt["ptyieldsAuAuc0092"], RatioAuAu["EtaToPion0092"]);}
          
          //d16-x01-y02
          binShift(*hEtaPt["ptyieldsAuAuc0020c"]);
          binShift(*hPionPt["ptyieldsAuAuc0020"]);
          if(hEtaPt["ptyieldsAuAuc0020c"]->sumW()>0 && hPionPt["ptyieldsAuAuc0020"]->sumW()>0){
          divide(hEtaPt["ptyieldsAuAuc0020c"], hPionPt["ptyieldsAuAuc0020"], RatioAuAu["EtaToPion0020"]);}
    
          //d16-x01-y03
          binShift(*hEtaPt["ptyieldsAuAuc2060c"]);
          binShift(*hPionPt["ptyieldsAuAuc2060"]);
          if(hEtaPt["ptyieldsAuAuc2060c"]->sumW()>0 && hPionPt["ptyieldsAuAuc2060"]->sumW()>0){
          divide(hEtaPt["ptyieldsAuAuc2060c"], hPionPt["ptyieldsAuAuc2060"], RatioAuAu["EtaToPion2060"]);}

          //d17-x01-y01
          binShift(*hEtaPt["ptyieldsAuAuc6092c"]);
          binShift(*hPionPt["ptyieldsAuAuc6092"]);
          if(hEtaPt["ptyieldsAuAuc6092c"]->sumW()>0 && hPionPt["ptyieldsAuAuc6092"]->sumW()>0){
          divide(hEtaPt["ptyieldsAuAuc6092c"], hPionPt["ptyieldsAuAuc6092"], RatioAuAu["EtaToPion6092"]);}
          
      }

      //histograms
      map<string, Histo1DPtr> hCrossSec;
      map<string, Histo1DPtr> hPionPt;
      map<string, Histo1DPtr> hEtaPt;
      map<string, Histo1DPtr> sigma_pp;
      
      //ratios
      map<string, Estimate1DPtr> RatioAuAu;
      map<string, Estimate1DPtr> RatiodAu;
      map<string, Estimate1DPtr> Ratiopp;
      
      //RdA, Raa
      map<string, Estimate1DPtr> hRda;
      map<string, Estimate1DPtr> hRaa;
      
      //Counters
      map<string, CounterPtr> sow;
      
      string beamOpt;
      string fixedcentralityOpt;
      double manualCentrality = -1.0;
      enum CollisionSystem { Unidentified, pp, AuAu200, DAu200 };
      CollisionSystem collSys;
      
      //these two vectors are only initialized if we would like to for loop over centrality bins
      //vector<int> AuAuCentralityBins{ 20, 60, 92 };
      //vector<int> dAuCentralityBins{ 20, 40, 60, 88 };

  };


  RIVET_DECLARE_PLUGIN(PHENIX_2007_I731133);

}
