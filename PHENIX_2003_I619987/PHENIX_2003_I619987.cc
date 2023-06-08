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
        //d01-x01-y02
        string refname2 = mkAxisCode(1, 1, 2);
        const Scatter2D& refdata2 = refData(refname2);
        book(hProtonPt["AuAuc2030"], refname2 + "_AuAuc2030_Proton",refdata2);
        book(hPionPosPt["AuAuc2030"], refname2 + "_AuAuc2030_PiPos",refdata2);
        book(RatioPtoPiPos["Proton/PiPos"], refname2);
        //d01-x01-y03
        string refname3 = mkAxisCode(1, 1, 3);
        const Scatter2D& refdata3 = refData(refname3);
        book(hProtonPt["AuAuc6092"], refname3 + "_AuAuc6092_Proton",refdata3);
        book(hPionPosPt["AuAuc6092"], refname3 + "_AuAuc6092_PiPos",refdata3);
        book(RatioPtoPiPos["Proton/PiPos"], refname3);
        //d01-x01-y04
        string refname4 = mkAxisCode(1, 1, 4);
        const Scatter2D& refdata4 = refData(refname4);
        book(hProBarPt["AuAuc0010"], refname4 + "_AuAuc0010_ProBar",refdata4);
        book(hPionNegPt["AuAuc0010"], refname4 + "_AuAuc0010_PiNeg",refdata4);
        book(RatioPBartoPiNeg["ProBar/PiNeg"], refname4);
        //d01-x01-y05
        string refname5 = mkAxisCode(1, 1, 5);
        const Scatter2D& refdata5 = refData(refname5);
        book(hProBarPt["AuAuc2030"], refname5 + "_AuAuc2030_ProBar",refdata5);
        book(hPionNegPt["AuAuc2030"], refname5 + "_AuAuc2030_PiNeg",refdata5);
        book(RatioPBartoPiNeg["ProBar/PiNeg"], refname5);
        //d01-x01-y06
        string refname6 = mkAxisCode(1, 1, 6);
        const Scatter2D& refdata6 = refData(refname6);
        book(hProBarPt["AuAuc6092"], refname6 + "_AuAuc6092_ProBar",refdata6);
        book(hPionNegPt["AuAuc6092"], refname6 + "_AuAuc6092_PiNeg",refdata6);
        book(RatioPBartoPiNeg["ProBar/PiNeg"], refname6);
        //d02-x01-y01
        string refname7 = mkAxisCode(2, 1, 1);
        const Scatter2D& refdata7 = refData(refname7);
        book(hProtonPt["AuAuc0010"], refname7 + "_AuAuc0010_Proton",refdata7);
        book(hPionNegPt["AuAuc0010"], refname7 + "_AuAuc0010_PiNeg",refdata7);
        book(RatioPtoPiNeg["Proton/PiNeg"], refname7);
        //d02-x01-y02
        string refname8 = mkAxisCode(2, 1, 2);
        const Scatter2D& refdata8 = refData(refname8);
        book(hProtonPt["AuAuc2030"], refname8 + "_AuAuc2030_Proton",refdata8);
        book(hPionNegPt["AuAuc2030"], refname8 + "_AuAuc2030_PiNeg",refdata8);
        book(RatioPtoPiNeg["Proton/PiNeg"], refname8);
        //d02-x01-y03
        string refname9 = mkAxisCode(2, 1, 3);
        const Scatter2D& refdata9 = refData(refname9);
        book(hProtonPt["AuAuc0010"], refname9 + "_AuAuc0010_Proton",refdata9);
        book(hPionNegPt["AuAuc0010"], refname9 + "_AuAuc0010_PiNeg",refdata9);
        book(RatioPtoPiNeg["Proton/PiNeg"], refname9);
        //d02-x01-y04
        string refname10 = mkAxisCode(2, 1, 4);
        const Scatter2D& refdata10 = refData(refname10);
        book(hProBarPt["AuAuc0010"], refname10 + "_AuAuc0010_ProBar",refdata10);
        book(hPionNegPt["AuAuc0010"], refname10 + "_AuAuc0010_PiNeg",refdata10);
        book(RatioPBartoPion["ProBar/Pion"], refname10);
        //d02-x01-y05
        string refname11 = mkAxisCode(2, 1, 5);
        const Scatter2D& refdata11 = refData(refname11);
        book(hProBarPt["AuAuc2030"], refname11 + "_AuAuc2030_ProBar",refdata11);
        book(hPionNegPt["AuAuc2030"], refname11 + "_AuAuc2030_PiNeg",refdata11);
        book(RatioPBartoPion["ProBar/Pion"], refname11);
        //d02-x01-y06
        string refname12 = mkAxisCode(2, 1, 6);
        const Scatter2D& refdata12 = refData(refname12);
        book(hProBarPt["AuAuc6092"], refname12 + "_AuAuc6092_ProBar",refdata12);
        book(hPionNegPt["AuAuc6092"], refname12 + "_AuAuc6092_PiNeg",refdata12);
        book(RatioPBartoPion["ProBar/Pion"], refname12);
        
        
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
                            hPionNegPt["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case 111: //pi^0
                        {
                            hPionPt["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case 2212: //p
                        {
                            hProtonPt["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            hProBarPt["AuAuc0010"]->fill(partPt);
                            break;
                        }
                    }
                    
                }
            }
            else if ((c >= 10.) && (c < 20.)){
                break;
            }
            
            else if ((c >= 20.) && (c < 30.))
            {
                for (const Particle& p : chargedParticles)
                {
                    double partPt = p.pT() / GeV;
                    //double pt_weight = 1. / (partPt * 2. * M_PI);  //Commented to avoid warning
                    
                    switch(p.pid()) {
                        case 211: //pi^+
                        {
                            hPionPosPt["AuAuc2030"]->fill(partPt);
                            break;
                        }
                        case -211: //pi^-
                        {
                            hPionNegPt["AuAuc2030"]->fill(partPt);
                            break;
                        }
                        case 111: //pi^0
                        {
                            hPionPt["AuAuc2030"]->fill(partPt);
                            break;
                        }
                        case 2212: //p
                        {
                            hProtonPt["AuAuc2030"]->fill(partPt);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            hProBarPt["AuAuc2030"]->fill(partPt);
                            break;
                        }
                    }
                    
                }
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
            else if ((c >= 60.) && (c < 92.)){
                for (const Particle& p : chargedParticles)
                {
                    double partPt = p.pT() / GeV;
                    //double pt_weight = 1. / (partPt * 2. * M_PI);  //Commented to avoid warning
                    
                    switch(p.pid()) {
                        case 211: //pi^+
                        {
                            hPionPosPt["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case -211: //pi^-
                        {
                            hPionNegPt["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case 111: //pi^0
                        {
                            hPionPt["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case 2212: //p
                        {
                            hProtonPt["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            hProBarPt["AuAuc6092"]->fill(partPt);
                            break;
                        }
                    }
                    
                }
            }
            
            else{
                break;
            }
        }
        
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
        
        //d01: ratios
        divide(hProtonPt["AuAuc0010"], hPionPosPt["AuAuc0010"], RatioPtoPiPos["AuAuc0010"]);
        divide(hProtonPt["AuAuc2030"], hPionPosPt["AuAuc2030"], RatioPtoPiPos["AuAuc2030"]);
        divide(hProtonPt["AuAuc6092"], hPionPosPt["AuAuc6092"], RatioPtoPiPos["AuAuc6092"]);
        divide(hProBarPt["AuAuc0010"], hPionNegPt["AuAuc0010"], RatioPBartoPiNeg["AuAuc0010"]);
        divide(hProBarPt["AuAuc2030"], hPionNegPt["AuAuc2030"], RatioPBartoPiNeg["AuAuc2030"]);
        divide(hProBarPt["AuAuc6092"], hPionNegPt["AuAuc6092"], RatioPBartoPiNeg["AuAuc6092"]);
        
        //d02
        divide(hProtonPt["AuAuc0010"], hPionNegPt["AuAuc0010"], RatioPtoPiNeg["AuAuc0010"]);
        divide(hProBarPt["AuAuc0010"], hPionPt["AuAuc0010"], RatioPBartoPion["AuAuc0010"]);
        divide(hProtonPt["AuAuc2030"], hPionNegPt["AuAuc2030"], RatioPtoPiNeg["AuAuc2030"]);
        divide(hProBarPt["AuAuc2030"], hPionPt["AuAuc2030"], RatioPBartoPion["AuAuc2030"]);
        divide(hProtonPt["AuAuc6092"], hPionNegPt["AuAuc6092"], RatioPtoPiNeg["AuAuc6092"]);
        divide(hProBarPt["AuAuc6092"], hPionPt["AuAuc6092"], RatioPBartoPion["AuAuc6092"]);
    }
    
    //Particles
    map<string, Histo1DPtr> hProtonPt; //p
    map<string, Histo1DPtr> hPionPosPt; //pi^+
    map<string, Histo1DPtr> hProBarPt; //p_bar
    map<string, Histo1DPtr> hPionNegPt; //pi^-
    map<string, Histo1DPtr> hPionPt; //pi^0
    
    //Ratios
    map<string, Scatter2DPtr> RatioPtoPiPos;
    map<string, Scatter2DPtr> RatioPBartoPiNeg;
    map<string, Scatter2DPtr> RatioPtoPiNeg;
    map<string, Scatter2DPtr> RatioPBartoPion;
    
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
