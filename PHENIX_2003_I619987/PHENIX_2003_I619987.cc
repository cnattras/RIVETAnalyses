// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include <math.h>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2003_I619987 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2003_I619987);
    
    // I figure that this is establishing binning for our pT
    bool getDeltaPt(YODA::Histo1D hist, double pT, double &deltaPt)
    {
        if(pT > xMin() && hist.xMax())
        {
            deltaPt = hist.bin(hist.binIndexAt(pT)).xMid();
            return true;
        }
        else return false;
    }
      
      
    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // Particles: pi^+, pi^-, pi^0, p, p_bar
      //pids (respectively): 211, -211, 111, 2212, -2212
      // all final-state particles within 
      // the given eta acceptance
      /// Found the cuts on page 3, paragraph 3
      const FinalState fs(Cuts::absrap < 0.35 && Cuts::phi == 0.392);
      declare(fs,"fs");

      // Declare centrality projection for centrality estimation
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX","CMULT","CMULT");
        
      //*****Counters******
      book(sow["sow_AuAu10"], "sow_AuAu10");
      book(sow["sow_AuAu20"], "sow_AuAu30");
      book(sow["sow_AuAu50"], "sow_AuAu50");
      book(sow["sow_AuAu50"], "sow_AuAu60");
      book(sow["sow_AuAu92"], "sow_AuAu92");
        
      //Ratio of protons/pions
      //d01-x01-y01
      book(hTemp_ratio_AuAu["PC10"], "PC10_AuAu", refData(1,1,1));
      book(hTemp_ratio_AuAu["PC30"], "PC30_AuAu", refData(1,1,2));
      book(hTemp_ratio_AuAu["PC92"], "PC92_AuAu", refData(1,1,3));
      book(hTemp_ratio_AuAu["P_barC10"], "P_barC10_AuAu", refData(1,1,4));
      book(hTemp_ratio_AuAu["P_barC30"], "P_barC30_AuAu", refData(1,1,5));
      book(hTemp_ratio_AuAu["P_barC92"], "P_barC92_AuAu", refData(1,1,6));
      book(hTemp_ratio_AuAu["PiplusC10"], "PiplusC10_AuAu", refData(1,1,1));
      book(hTemp_ratio_AuAu["PiplusC30"], "PiplusC30_AuAu", refData(1,1,2));
      book(hTemp_ratio_AuAu["PiplusC92"], "PiplusC92_AuAu", refData(1,1,3));
      book(hTemp_ratio_AuAu["PiminusC10"], "Piminus10_AuAu", refData(1,1,4));
      book(hTemp_ratio_AuAu["PiminusC30"], "PiminusC30_AuAu", refData(1,1,5));
      book(hTemp_ratio_AuAu["PiminusC92"], "PiminusC92_AuAu", refData(1,1,6));
      book(RatioAuAu["PtoPiplusC10"], 1,1,1);
      book(RatioAuAu["PtoPiplusC30"], 1,1,2);
      book(RatioAuAu["PtoPiplusC92"], 1,1,3);
      book(RatioAuAu["P_bartoPiminusC10"], 1,1,4);
      book(RatioAuAu["P_bartoPiminusC30"], 1,1,5);
      book(RatioAuAu["P_bartoPiminusC92"], 1,1,6);
        
       
    
        
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      //A reminder of pids as writing switch statements:
        // Particles: pi^+, pi^-, pi^0, p, p_bar
        //pids (respectively): 211, -211, 111, 2212, -2212
        
      Particles chargedP = applyProjection<FinalState>(event,"fs").partices();

      // All figures are for S_NN = 200 GeV collisions, so no if statement required
        
      /// We will need to write for centrality
      const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
      const double cent = centProj();
        
      // events with centrality < 0 or > 92 are invalid. We use vetoEvent.
      //*** @Christal, is a . needed after 92? ***
      if (cent < 0. || cent > 92) vetoEvent;
        
      // 0-10% centrality
      if (cent ) 0. && cent < 10.) {
          //fill our counter
          sow["sow_AUAU10"]->fill();
          for (const Particle& p : chargedParticles)
          {
              double partPt = p.pT() / GeV;
              double pt_weight = 1. / (partPt * 2. * M_PI);
              
              switch(p.pid()) {
                  case 211: //pi^+
                      if (getDeltaPt(*hTemp_ratio_AuAu["PiplusC10"], partPt, deltaPt){
                          pt_weight /= deltaPt;
                          hTemp_ratio_AuAu["PiplusC10"]->fill(partPt, pt_weight);
                      }
                  case -211: //pi^-
                  case 111: //pi^0
                  case 2212: //p
                  case -2212: //p_bar
                
              }
          }
          
            }
      else if ((c >= 10.) && (c < 20.)){
          sow["sow_AUAU20"]->fill();
            }
        
      else if ((c >= 20.) && (c < 30.)){
          sow["sow_AUAU30"]->fill();
            }
          
      else if ((c >= 30.) && (c < 40.)){
          sow["sow_AUAU40"]->fill();
            }
        
      else if ((c >= 40.) && (c < 50.)){
          sow["sow_AUAU50"]->fill();
            }
        
      else if ((c >= 50.) && (c < 60.)){
          sow["sow_AUAU60"]->fill();
            }
        
      else{
          sow["sow_AUAU92"]->fill();
            }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h["YYYY"]); // normalize to unity
      scale(_h["ZZZZ"], crossSection()/picobarn/sumOfWeights()); // norm to cross section

    }



    map<string, Histo1DPtr> hAuAu_Yields;
    
    map<string, Histo1DPtr> hRcp;
    map<string, Histo1DPtr> hRAA;
    
    map<string, CounterPtr> sow;
    CollisionSystem collSys;
    vector<int> AUAUCentralityBins{ 10, 20, 30, 40, 50, 60, 92};

      
  };


  DECLARE_RIVET_PLUGIN(PHENIX_2003_I619987);


}
