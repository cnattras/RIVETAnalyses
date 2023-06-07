// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
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
class PHENIX_2003_I619987 : public Analysis {
public:
    
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2003_I619987);
    
    
    
    
    /// Book histograms and initialise projections before the run
    void init() {
        
        // Initialise and register projections
        
        // Particles: pi^+, pi^-, pi^0, p, p_bar
        //pids (respectively): 211, -211, 111, 2212, -2212
        // all final-state particles within
        // the given eta acceptance
        /// Found the cuts on page 3, paragraph 3
        std::initializer_list<int> pdgIds = { 211, -211, 111, 2212, -2212 };
        const PrimaryParticles fs(pdgIds, Cuts::absrap < 0.35 && Cuts::phi == 0.392);
        declare(fs,"fs");
        
        //Collision system
        //if (beamOpt == "PP") collSys = pp;
        //else if (beamOpt == "AUAU200") collSys = AuAu200;
        
        // Declare centrality projection for centrality estimation
        //if (!(collSys == pp))
        declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX","CMULT","CMULT");
        
        //*****Counters******
        //book(sow["sow_AuAu10"], "sow_AuAu10");
        
        
        //Ratio of protons/pions
        //d01-x01-y01
        string refname1 = mkAxisCode(1, 1, 1);
        const Scatter2D& refdata1 = refData(refname1);
        book(hProtonPt["AuAuc0010"], refname1 + "_AuAuc0010_Proton",refdata1);
        book(hPionPosPt["AuAuc0010"], refname1 + "_AuAuc0010_PiPos",refdata1);
        book(RatioPtoPiPos["Proton/PiPos"], refname1);
        
        
        
        
    }
    
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
        
        //A reminder of pids as writing switch statements:
        // Particles: pi^+, pi^-, pi^0, p, p_bar
        //pids (respectively): 211, -211, 111, 2212, -2212
        
        Particles chargedParticles = applyProjection<PrimaryParticles>(event,"fs").particles();
        
        // All figures are for S_NN = 200 GeV collisions, so no if statement required for collSys
        for(Particle p : chargedParticles)
        {
            /// We will need to write for centrality
            const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
            const double c = cent();
            
            // events with centrality < 0 or > 92 are invalid. We use vetoEvent.
            if (c < 0. || c > 92.) vetoEvent;
            
            // 0-10% centrality
            if ((c > 0.) && (c < 10.))
            {
                for (const Particle& p : chargedParticles)
                {
                    double partPt = p.pT() / GeV;
                    //double pt_weight = 1. / (partPt * 2. * M_PI);  //Commented to avoid warning
                    
                    switch(p.pid()) {
                        case 211: //pi^+
                        {
                            hPionPosPt["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case -211: //pi^-
                        {
                            break;
                        }
                        case 111: //pi^0
                        {
                            break;
                        }
                        case 2212: //p
                        {
                            hProtonPt["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            break;
                        }
                    }
                    
                }
            }
            else if ((c >= 10.) && (c < 20.)){
                break;
            }
            
            else if ((c >= 20.) && (c < 30.)){
                break;
            }
            
            else if ((c >= 30.) && (c < 40.)){
                break;
            }
            
            else if ((c >= 40.) && (c < 50.)){
                break;
            }
            
            else if ((c >= 50.) && (c < 60.)){
                break;
            }
            
            else{
                break;
            }
        }
        
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
        
        divide(hProtonPt["AuAuc0010"], hPionPosPt["AuAuc0010"], RatioPtoPiPos["AuAuc0010"]);
        
    }
    
    //Particles
    map<string, Histo1DPtr> hProtonPt;
    map<string, Histo1DPtr> hPionPosPt;
    
    //Ratios
    map<string, Scatter2DPtr> RatioPtoPiPos;
    
    //Rcp, Raa **TBD**
    //map<string, Scatter2DPtr> hRcp;
    //map<string, Scatter2DPtr> hRaa;
    
    //Counter
    //map<string, CounterPtr> sow;
    
    //Initialize collision system and AUAU centrality bins
    string beamOpt;
    enum CollisionSystem {AuAu200};
    CollisionSystem collSys;
    vector<int> AUAUCentralityBins{ 10, 20, 30, 40, 50, 60, 92};
    
    
    
  };
  DECLARE_RIVET_PLUGIN(PHENIX_2003_I619987);
}
