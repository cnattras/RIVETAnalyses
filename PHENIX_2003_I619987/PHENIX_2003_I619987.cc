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
    
    // I figure that this is establishing binning for out pT
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
      // all final-state particles within 
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);

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

      /// @todo Do the event by event analysis here

      // retrieve dressed leptons, sorted by pT
      vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();

      // retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

      // remove all jets within dR < 0.2 of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

      // select jets ghost-associated to B-hadrons with a certain fiducial selection
      Jets bjets = filter_select(jets, [](const Jet& jet) {
        return  jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
      });

      // veto event if there are no b-jets
      if (bjets.empty())  vetoEvent;

      // apply a missing-momentum cut
      if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent;

      // fill histogram with leading b-jet pT
      _h["XXXX"]->fill(bjets[0].pT()/GeV);

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


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2003_I619987);


}
