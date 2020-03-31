// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <math.h>
#include "Rivet/Projections/BackgroundRho.hh"
#include <Rivet/Projections/HepMCHeavyIon.hh>
#include<bits/stdc++.h> 
#define _USE_MATH_DEFINES

namespace Rivet {
    
    void BuildResponseMatrix(YODA::Histo1D& DeltaPt, YODA::Histo2D& RM)
    {
        
        YODA::Histo1D axisX(RM.numBinsX()-1, RM.xMin(), RM.xMax());
        YODA::Histo1D axisY(RM.numBinsY()-1, RM.yMin(), RM.yMax());
                
        for(int iBinY = 0; iBinY < RM.numBinsY(); iBinY++)
        {
            double binYCenter = (axisY.bin(iBinY).xMax() + axisY.bin(iBinY).xMin())/2.;
            
            for(int iBinX = 0; iBinX < RM.numBinsX(); iBinX++)
            {
                double binXCenter = (axisX.bin(iBinX).xMax() + axisX.bin(iBinX).xMin())/2.;
                double binValue = 0.;
                int entries = 0;
                                
                for(int i = 0; i < DeltaPt.numBins(); i++)
                {
                    if(DeltaPt.bin(i).xMin() >= (axisX.bin(iBinX).xMin()-binYCenter) && DeltaPt.bin(i).xMax() <= (axisX.bin(iBinX).xMax()-binYCenter))
                    {
                        binValue += DeltaPt.bin(i).sumW();
                        entries += DeltaPt.bin(i).numEntries();
                    }
                    
                }
                
                RM.fill(binXCenter, binYCenter, binValue, entries);
                
            }
        }
    }
    
    double Distance(Particle& p, double eta, double phi)
    {
        double deltaEta = abs(p.eta() - eta);
        double deltaPhi = abs(p.phi() - phi);
        if(deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;
        
        double distance = sqrt((deltaEta*deltaEta) + (deltaPhi*deltaPhi));
        
        return distance;
    }
    
  /// @brief Add a short analysis description here
  class ResponseMtrix_2020_Example : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ResponseMtrix_2020_Example);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      ALICE::V0AndTrigger v0and;
      declare<ALICE::V0AndTrigger>(v0and,"V0-AND");
      
      // Centrality projection.
      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M","V0M");
      
      //HepMC
      declare(HepMCHeavyIon(), "HepMC");
        
      // the basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::pT > 150*MeV && Cuts::abseta < 0.9);
      declare(fs,"fs");

      // jets à la ALICE - Jet area will be available using the pseudojet
      fastjet::AreaType fjAreaType = fastjet::active_area_explicit_ghosts;
      //fastjet::GhostedAreaSpec(MaxRap, NGhostRepeats, GhostArea, GridScatter, KtScatter, MeanGhostKt);
      fastjet::GhostedAreaSpec fjGhostAreaSpec = fastjet::GhostedAreaSpec(1., 1, 0.005, 1., 0.1, 1e-100);
      fjAreaDef = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
      fjAreaDefKT = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
      
      //Anti-kt algorithm: For signal
      const FastJets jetsFJ(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, jetR, fjAreaDef, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetsFJ, "jets");
    
      //kt algorithm: For background
      const FastJets jetsKT(fs, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::pt_scheme, jetR, fjAreaDefKT, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetsKT, "jetsKT");
      
      Cut c = Cuts::pT > 0.15*GeV && Cuts::abseta < (0.9-jetR);
            
      BackgroundRho projRho(jetsKT, removeNLeadJets, jetAreaCut, c);
      declare(projRho,"projRho");
      
      book(hJetPt_0_10, "hJetPt_0_10", 100, 0., 200.);
      book(hJetPt_30_50, "hJetPt_30_50", 100, 0., 200.);
      
      book(hJetPt_0_10_BkgSubtracted, "hJetPt_0_10_BkgSubtracted", 15, -50., 100.);
      book(hJetPt_30_50_BkgSubtracted, "hJetPt_30_50_BkgSubtracted", 15, -50., 100.);
      
      book(hDeltaPt_0_10, "hDeltaPt_0_10", 150, -50., 100.);
      book(hDeltaPt_30_50, "hDeltaPt_30_50", 150, -50., 100.);
      
      book(ResponseMatrix_0_10, "ResponseMatrix_0_10", 15, -50., 100, 9, 10, 100);
      book(ResponseMatrix_30_50, "ResponseMatrix_30_50", 15, -50., 100, 9, 10, 100);
      
      book(sow10, "sow10");
      book(sow30, "sow30");
      
}


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        
      // Event trigger.
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")()) vetoEvent;

      // The centrality projection.
      const CentralityProjection& centProj = apply<CentralityProjection>(event,"V0M");

      // The centrality.
      const double cent = centProj();
      
      //Event selection
      if(cent > 50.)
      {
          vetoEvent;          
      }
      
      //Get jets
      FastJets fj = applyProjection<FastJets>(event, "jets");
      Jets jets = fj.jetsByPt(Cuts::abseta < 0.5 && Cuts::pT > 0.15*GeV);
      
      //Get rho
      BackgroundRho projRho = applyProjection<BackgroundRho>(event, "projRho");
      double rho = projRho.getRho();
      
      //Jet loop
      if(cent < 10.)
      {
          sow10->fill();
          for(Jet jet : jets)
          {
              if(jet.pseudojet().area() < jetAreaCut) continue;
              hJetPt_0_10->fill(jet.pT()/GeV);
              double jetPt_corr = jet.pT()/GeV - rho*jet.pseudojet().area();
              hJetPt_0_10_BkgSubtracted->fill(jetPt_corr);
          }
          
      }
      else if(cent > 30. && cent < 50.)
      {
          sow30->fill();
          for(Jet jet : jets)
          {
              if(jet.pseudojet().area() < jetAreaCut) continue;
              hJetPt_30_50->fill(jet.pT()/GeV);
              double jetPt_corr = jet.pT()/GeV - rho*jet.pseudojet().area();
              hJetPt_30_50_BkgSubtracted->fill(jetPt_corr);
          }
      }
      
      FinalState fs = applyProjection<FinalState>(event, "fs");
      Particles fsparticles = fs.particles();
      
      
      
      //Random Numbers Generator
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::default_random_engine generator(seed);
      std::uniform_real_distribution<double> etaDistrib(-0.5,0.5);
      std::uniform_real_distribution<double> phiDistrib(0.,2*M_PI);
      
      int Ncones = 1;
      
      for(int icone = 0; icone < Ncones; icone++)
      {
          double eta_roll = etaDistrib(generator);
          double phi_roll = phiDistrib(generator);
          
          double conePt = 0.;
          
          for(Particle& p : fsparticles)
          {
              if(Distance(p,eta_roll,phi_roll) < jetR) conePt += p.pT()/GeV;
          }
              
          if(conePt < 0.15)
          {
              continue;
              cout << "Cone zero" << endl;
          }
              
          double deltaPt = conePt - (rho*M_PI*jetR*jetR);
          
          if(cent < 10.)
          {
              hDeltaPt_0_10->fill(deltaPt);
          }
          else if(cent > 30. && cent < 50.)
          {              
              hDeltaPt_30_50->fill(deltaPt);
          }
      }
      
      
      
            
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      
      hJetPt_0_10->scaleW(1./sow10->sumW());
      hJetPt_0_10_BkgSubtracted->scaleW(1./sow10->sumW());
      
      hDeltaPt_0_10->normalize();
      
      hJetPt_30_50->scaleW(1./sow30->sumW());
      hJetPt_30_50_BkgSubtracted->scaleW(1./sow30->sumW());
      
      hDeltaPt_30_50->normalize();
      
      BuildResponseMatrix(*hDeltaPt_0_10, *ResponseMatrix_0_10);
      BuildResponseMatrix(*hDeltaPt_30_50, *ResponseMatrix_30_50);
      
    }
    
    

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr hJetPt_0_10;
    Histo1DPtr hJetPt_0_10_BkgSubtracted;
    Histo1DPtr hDeltaPt_0_10;
    Histo2DPtr ResponseMatrix_0_10;
    
    Histo1DPtr hJetPt_30_50;
    Histo1DPtr hJetPt_30_50_BkgSubtracted;
    Histo1DPtr hDeltaPt_30_50;
    Histo2DPtr ResponseMatrix_30_50;
    
    CounterPtr sow10;
    CounterPtr sow30;
        
    const int removeNLeadJets = 2;

    double jetR = 0.4;
    double jetAreaCut = 0.557*M_PI*jetR*jetR;
    fastjet::AreaDefinition *fjAreaDef;
    fastjet::AreaDefinition *fjAreaDefKT;
  };


  DECLARE_RIVET_PLUGIN(ResponseMtrix_2020_Example);

}
