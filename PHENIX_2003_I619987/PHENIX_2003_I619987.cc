// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/UnstableParticles.hh"
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
        std::initializer_list<int> pdgIds = { 211, -211, 2212, -2212 };

        
        //Charged particles
        //consider adding back && Cuts::phi == 0.392
        const PrimaryParticles cp(pdgIds, Cuts::absrap < 0.35 && Cuts::pT > 0.5*GeV);
        declare(cp,"cp");
        //Neutral particles
        const UnstableParticles np(Cuts::absrap < 0.35 && Cuts::abspid == 111 && Cuts::pT > 0.5*GeV);
        declare(np,"np");
        
        
        //Collision system
        //if (beamOpt == "PP200") collSys = pp;
        //else if (beamOpt == "AUAU200") collSys = AuAu200;
        
        // Declare centrality projection for centrality estimation
        //if (!(collSys == pp))
        declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX","CMULT","CMULT");
        
        
        //Ratio of protons/pions
        //note: the "a" added to hProtonPt; since there are different binnings and both d01 and d02 use hProtonPt
        //Figure 1
        //d01-x01-y01
        string refname1 = mkAxisCode(1, 1, 1);
        const Scatter2D& refdata1 = refData(refname1);
        book(hProtonPt["AuAuc0010a"], refname1 + "_AuAuc0010_Proton",refdata1);
        book(hPionPosPt["AuAuc0010"], refname1 + "_AuAuc0010_PiPos",refdata1);
        book(RatioPtoPiPos["AuAuc0010"], refname1);
        //d01-x01-y02
        string refname2 = mkAxisCode(1, 1, 2);
        const Scatter2D& refdata2 = refData(refname2);
        book(hProtonPt["AuAuc2030a"], refname2 + "_AuAuc2030_Proton",refdata2);
        book(hPionPosPt["AuAuc2030"], refname2 + "_AuAuc2030_PiPos",refdata2);
        book(RatioPtoPiPos["AuAuc2030"], refname2);
        //d01-x01-y03
        string refname3 = mkAxisCode(1, 1, 3);
        const Scatter2D& refdata3 = refData(refname3);
        book(hProtonPt["AuAuc6092a"], refname3 + "_AuAuc6092_Proton",refdata3);
        book(hPionPosPt["AuAuc6092"], refname3 + "_AuAuc6092_PiPos",refdata3);
        book(RatioPtoPiPos["AuAuc6092"], refname3);
        //d01-x01-y04
        string refname4 = mkAxisCode(1, 1, 4);
        const Scatter2D& refdata4 = refData(refname4);
        book(hProBarPt["AuAuc0010a"], refname4 + "_AuAuc0010_ProBar",refdata4);
        book(hPionNegPt["AuAuc0010a"], refname4 + "_AuAuc0010_PiNeg",refdata4);
        book(RatioPBartoPiNeg["AuAuc0010"], refname4);
        //d01-x01-y05
        string refname5 = mkAxisCode(1, 1, 5);
        const Scatter2D& refdata5 = refData(refname5);
        book(hProBarPt["AuAuc2030a"], refname5 + "_AuAuc2030_ProBar",refdata5);
        book(hPionNegPt["AuAuc2030a"], refname5 + "_AuAuc2030_PiNeg",refdata5);
        book(RatioPBartoPiNeg["AuAuc2030"], refname5);
        //d01-x01-y06
        string refname6 = mkAxisCode(1, 1, 6);
        const Scatter2D& refdata6 = refData(refname6);
        book(hProBarPt["AuAuc6092a"], refname6 + "_AuAuc6092_ProBar",refdata6);
        book(hPionNegPt["AuAuc6092a"], refname6 + "_AuAuc6092_PiNeg",refdata6);
        book(RatioPBartoPiNeg["AuAuc6092"], refname6);
        //d02-x01-y01
        string refname7 = mkAxisCode(2, 1, 1);
        const Scatter2D& refdata7 = refData(refname7);
        book(hProtonPt["AuAuc0010b"], refname7 + "_AuAuc0010_Proton",refdata7);
        book(hPionNegPt["AuAuc0010b"], refname7 + "_AuAuc0010_PiNeg",refdata7);
        book(RatioPtoPiNeg["AuAuc0010"], refname7);
        //d02-x01-y02
        string refname8 = mkAxisCode(2, 1, 2);
        const Scatter2D& refdata8 = refData(refname8);
        book(hProtonPt["AuAuc2030b"], refname8 + "_AuAuc2030_Proton",refdata8);
        book(hPionNegPt["AuAuc2030b"], refname8 + "_AuAuc2030_PiNeg",refdata8);
        book(RatioPtoPiNeg["AuAuc2030"], refname8);
        //d02-x01-y03
        string refname9 = mkAxisCode(2, 1, 3);
        const Scatter2D& refdata9 = refData(refname9);
        book(hProtonPt["AuAuc6092b"], refname9 + "_AuAuc6092_Proton",refdata9);
        book(hPionNegPt["AuAuc6092b"], refname9 + "_AuAuc6092_PiNeg",refdata9);
        book(RatioPtoPiNeg["AuAuc6092"], refname9);
        //d02-x01-y04
        string refname10 = mkAxisCode(2, 1, 4);
        const Scatter2D& refdata10 = refData(refname10);
        book(hProBarPt["AuAuc0010b"], refname10 + "_AuAuc0010_ProBar",refdata10);
        book(hPionPt["AuAuc0010a"], refname10 + "_AuAuc0010_Pion",refdata10);
        book(RatioPBartoPion["AuAuc0010"], refname10);
        //d02-x01-y05
        string refname11 = mkAxisCode(2, 1, 5);
        const Scatter2D& refdata11 = refData(refname11);
        book(hProBarPt["AuAuc2030b"], refname11 + "_AuAuc2030_ProBar",refdata11);
        book(hPionPt["AuAuc2030a"], refname11 + "_AuAuc2030_Pion",refdata11);
        book(RatioPBartoPion["AuAuc2030"], refname11);
        //d02-x01-y06
        string refname12 = mkAxisCode(2, 1, 6);
        const Scatter2D& refdata12 = refData(refname12);
        book(hProBarPt["AuAuc6092b"], refname12 + "_AuAuc6092_ProBar",refdata12);
        book(hPionPt["AuAuc6092a"], refname12 + "_AuAuc6092_Pion",refdata12);
        book(RatioPBartoPion["AuAuc6092"], refname12);
        
        //yields
        //figure 2
        //counters per each centrality
        book(sow["AuAuc0010"], "_sow_AuAuc0010");
        book(sow["AuAuc2030"], "_sow_AuAuc2030");
        book(sow["AuAuc4050"], "_sow_AuAuc4050");
        book(sow["AuAuc6092"], "_sow_AuAuc6092");
        
        //d03-x01-y01
        string refname13 = mkAxisCode(3, 1, 1);
        //const Scatter2D& refdata13 = refData(refname13);
        book(hProtonPt["ptyieldsAuAuc0010"], 3, 1, 1);
        //d03-x01-y02
        string refname14 = mkAxisCode(3, 1, 2);
        //const Scatter2D& refdata14 = refData(refname14);
        book(hProtonPt["ptyieldsAuAuc2030"], 3, 1, 2);
        //d03-x01-y03
        string refname15 = mkAxisCode(3, 1, 3);
        //const Scatter2D& refdata15 = refData(refname15);
        book(hProtonPt["ptyieldsAuAuc4050"], 3, 1, 3);
        //d03-x01-y04
        string refname16 = mkAxisCode(3, 1, 4);
        //const Scatter2D& refdata16 = refData(refname16);
        book(hProtonPt["ptyieldsAuAuc6092"], 3, 1, 4);
        //d03-x01-y05
        string refname17 = mkAxisCode(3, 1, 5);
        //const Scatter2D& refdata17 = refData(refname17);
        book(hProBarPt["ptyieldsAuAuc0010"], 3, 1, 5);
        //d03-x01-y06
        string refname18 = mkAxisCode(3, 1, 6);
        //const Scatter2D& refdata18 = refData(refname18);
        book(hProBarPt["ptyieldsAuAuc2030"], 3, 1, 6);
        //d03-x01-y07
        string refname19 = mkAxisCode(3, 1, 7);
        //const Scatter2D& refdata19 = refData(refname19);
        book(hProBarPt["ptyieldsAuAuc4050"], 3, 1, 7);
        //d03-x01-y08
        string refname20 = mkAxisCode(3, 1, 8);
        //const Scatter2D& refdata20 = refData(refname20);
        book(hProBarPt["ptyieldsAuAuc6092"], 3, 1, 8);
        
        //Rcp
        //Figure 3a
        //d04-x01-y01
        string refname21 = mkAxisCode(4, 1, 1);
        const Scatter2D& refdata21 = refData(refname21);
        book(hPPlusPBarPt["ppluspbarAuAuc0010"], refname21 + "_AuAuc0010",refdata21);
        book(hPPlusPBarPt["ppluspbarAuAuc6092"], refname21 + "_AuAuc6092",refdata21);
        book(hRcp["ppluspbar"], refname21);
        
        //Figure 3b
        //d05-x01-y01
        string refname22 = mkAxisCode(5, 1, 1);
        const Scatter2D& refdata22 = refData(refname22);
        //next two lines: not needed since hPionPt is already booked and filled for //d02
        book(hPionPt["AuAuc0010b"], refname22 + "_AuAuc0010",refdata22);
        book(hPionPt["AuAuc6092b"], refname22 + "_AuAuc6092",refdata22);
        book(hRcp["pion"], refname22);
        
        //Figure 4a
        //d06-x01-y01
        string refname23 = mkAxisCode(6, 1, 1);
        const Scatter2D& refdata23 = refData(refname23);
        book(hChHadrons["AuAuc0010"], refname23 + "_AuAuc0010",refdata23);
        book(hPionPt["AuAuc0010c"], refname23 + "_AuAuc0010_Pion",refdata23);
        book(RatioHadtoPion["AuAuc0010"], refname23);
        
        //Figure 4b
        //d07-x01-y01
        string refname24 = mkAxisCode(7, 1, 1);
        const Scatter2D& refdata24 = refData(refname24);
        book(hChHadrons["AuAuc6092"], refname24 + "_AuAuc6092",refdata24);
        book(hPionPt["AuAuc6092c"], refname24 + "_AuAuc6092_Pion",refdata24);
        book(RatioHadtoPion["AuAuc6092"], refname24);
    }
    
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
        
        //A reminder of pids as writing switch statements:
        // Particles: pi^+, pi^-, pi^0, p, p_bar
        //pids (respectively): 211, -211, 111, 2212, -2212
        
        Particles chargedParticles = applyProjection<PrimaryParticles>(event,"cp").particles();
        Particles neutralParticles = applyProjection<UnstableParticles>(event,"np").particles();
        // All figures are for S_NN = 200 GeV collisions, so no if statement required for collSys
        
        /// Case for charged particles
        for(Particle p : chargedParticles)
        {
            // We will need to write for centrality
            const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
            const double c = cent();
            
            // events with centrality < 0 or > 92 are invalid. We use vetoEvent.
            if (c < 0. || c > 92.) vetoEvent;
            
            // 0-10% centrality
            if ((c > 0.) && (c < 10.))
            {
                //sow["sow_AuAuc0010"]->fill();
                for (const Particle& p : chargedParticles)
                {
                    double partPt = p.pT() / GeV;
                    double pt_weight = 1. / (partPt * 2. * M_PI);
                    
                    switch(p.pid()) {
                        case 211: //pi^+
                        {
                            hPionPosPt["AuAuc0010"]->fill(partPt);
                            hChHadrons["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case -211: //pi^-
                        {
                            hPionNegPt["AuAuc0010a"]->fill(partPt);
                            hPionNegPt["AuAuc0010b"]->fill(partPt);
                            hChHadrons["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case 2212: //p
                        {
                            hProtonPt["AuAuc0010a"]->fill(partPt);
                            hProtonPt["AuAuc0010b"]->fill(partPt);
                            hProtonPt["ptyieldsAuAuc0010"]->fill(partPt, pt_weight);
                            hPPlusPBarPt["ppluspbarAuAuc0010"]->fill(partPt);
                            hChHadrons["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            hProBarPt["AuAuc0010a"]->fill(partPt);
                            hProBarPt["AuAuc0010b"]->fill(partPt);
                            hProBarPt["ptyieldsAuAuc0010"]->fill(partPt, pt_weight);
                            hPPlusPBarPt["ppluspbarAuAuc0010"]->fill(partPt);
                            hChHadrons["AuAuc0010"]->fill(partPt);
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
                    double pt_weight = 1. / (partPt * 2. * M_PI);  //Commented to avoid warning
                    
                    switch(p.pid()) {
                        case 211: //pi^+
                        {
                            hPionPosPt["AuAuc2030"]->fill(partPt);
                            break;
                        }
                        case -211: //pi^-
                        {
                            hPionNegPt["AuAuc2030a"]->fill(partPt);
                            hPionNegPt["AuAuc2030b"]->fill(partPt);
                            break;
                        }
                        case 2212: //p
                        {
                            hProtonPt["AuAuc2030a"]->fill(partPt);
                            hProtonPt["AuAuc2030b"]->fill(partPt);
                            hProtonPt["ptyieldsAuAuc2030"]->fill(partPt, pt_weight);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            hProBarPt["AuAuc2030a"]->fill(partPt);
                            hProBarPt["AuAuc2030b"]->fill(partPt);
                            hProBarPt["ptyieldsAuAuc2030"]->fill(partPt, pt_weight);
                            break;
                        }
                    }
                    
                }
                
            }
            
            else if ((c >= 30.) && (c < 40.))
            {
                break;
            }
            
            else if ((c >= 40.) && (c < 50.))
            {
                for (const Particle& p : chargedParticles)
                {
                    double partPt = p.pT() / GeV;
                    double pt_weight = 1. / (partPt * 2. * M_PI);  //Commented to avoid warning
                    
                    switch(p.pid()) {
                        case 211: //pi^+
                        {
                            break;
                        }
                        case -211: //pi^-
                        {
                            break;
                        }
                        case 2212: //p
                        {
                            hProtonPt["ptyieldsAuAuc4050"]->fill(partPt, pt_weight);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            hProBarPt["ptyieldsAuAuc4050"]->fill(partPt, pt_weight);
                            break;
                        }
                    }
                    
                }
            }
            
            else if ((c >= 50.) && (c < 60.)){
                break;
            }
            
            else if ((c >= 60.) && (c < 92.)){
                for (const Particle& p : chargedParticles)
                {
                    double partPt = p.pT() / GeV;
                    double pt_weight = 1. / (partPt * 2. * M_PI);
                    
                    switch(p.pid()) {
                        case 211: //pi^+
                        {
                            hPionPosPt["AuAuc6092"]->fill(partPt);
                            hChHadrons["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case -211: //pi^-
                        {
                            hPionNegPt["AuAuc6092a"]->fill(partPt);
                            hPionNegPt["AuAuc6092b"]->fill(partPt);
                            hChHadrons["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case 2212: //p
                        {
                            hProtonPt["AuAuc6092a"]->fill(partPt);
                            hProtonPt["AuAuc6092b"]->fill(partPt);
                            hProtonPt["ptyieldsAuAuc6092"]->fill(partPt, pt_weight);
                            hPPlusPBarPt["ppluspbarAuAuc6092"]->fill(partPt);
                            hChHadrons["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            hProBarPt["AuAuc6092a"]->fill(partPt);
                            hProBarPt["AuAuc6092b"]->fill(partPt);
                            hProBarPt["ptyieldsAuAuc6092"]->fill(partPt, pt_weight);
                            hPPlusPBarPt["ppluspbarAuAuc6092"]->fill(partPt);
                            hChHadrons["AuAuc6092"]->fill(partPt);
                            break;
                        }
                    }
                    
                }
            }
            
            else{
                break;
            }
        }
        
        /// Case for neutral particles
        for(Particle p : neutralParticles)
        {
            /// We will need to write for centrality
            const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
            const double c = cent();
            
            // events with centrality < 0 or > 92 are invalid. We use vetoEvent.
            if (c < 0. || c > 92.) vetoEvent;
            
            // 0-10% centrality
            if ((c > 0.) && (c < 10.))
            {
                for (const Particle& p : neutralParticles)
                {
                    double partPt = p.pT() / GeV;
                    //double pt_weight = 1. / (partPt * 2. * M_PI);  //Commented to avoid warning
                    
                    switch(p.pid()) {
                        case 111: //pi^0
                        {
                            hPionPt["AuAuc0010a"]->fill(partPt);
                            hPionPt["AuAuc0010b"]->fill(partPt);
                            hPionPt["AuAuc0010c"]->fill(partPt);
                            break;
                        }
                    }
                    
                }
            }
            else if ((c >= 10.) && (c < 20.))
            {
                break;
            }
            
            else if ((c >= 20.) && (c < 30.))
            {
                for (const Particle& p : neutralParticles)
                {
                    double partPt = p.pT() / GeV;
                    //double pt_weight = 1. / (partPt * 2. * M_PI);  //Commented to avoid warning
                    
                    switch(p.pid()) {
                        case 111: //pi^0
                        {
                            hPionPt["AuAuc2030a"]->fill(partPt);
                            break;
                        }
                    }
                    
                }
                
            }
            
            else if ((c >= 30.) && (c < 40.))
            {
                break;
            }
            
            else if ((c >= 40.) && (c < 50.))
            {
                break;
            }
            
            else if ((c >= 50.) && (c < 60.))
            {
                break;
            }
            else if ((c >= 60.) && (c < 92.))
            {
                for (const Particle& p : neutralParticles)
                {
                    double partPt = p.pT() / GeV;
                    //double pt_weight = 1. / (partPt * 2. * M_PI);  //Commented to avoid warning
                    
                    switch(p.pid()) {
                        case 111: //pi^0
                        {
                            hPionPt["AuAuc6092a"]->fill(partPt);
                            hPionPt["AuAuc6092b"]->fill(partPt);
                            hPionPt["AuAuc6092c"]->fill(partPt);
                            break;
                        }
                    }
                    
                }
            }
            
            else
            {
                break;
            }
        }
        
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
        
        //d01: p and p_bar ratios
        divide(hProtonPt["AuAuc0010a"], hPionPosPt["AuAuc0010"], RatioPtoPiPos["AuAuc0010"]);
        divide(hProtonPt["AuAuc2030a"], hPionPosPt["AuAuc2030"], RatioPtoPiPos["AuAuc2030"]);
        divide(hProtonPt["AuAuc6092a"], hPionPosPt["AuAuc6092"], RatioPtoPiPos["AuAuc6092"]);
        divide(hProBarPt["AuAuc0010a"], hPionNegPt["AuAuc0010a"], RatioPBartoPiNeg["AuAuc0010"]);
        divide(hProBarPt["AuAuc2030a"], hPionNegPt["AuAuc2030a"], RatioPBartoPiNeg["AuAuc2030"]);
        divide(hProBarPt["AuAuc6092a"], hPionNegPt["AuAuc6092a"], RatioPBartoPiNeg["AuAuc6092"]);
        
        //d02 p and p_bar ratios
        divide(hProtonPt["AuAuc0010b"], hPionNegPt["AuAuc0010b"], RatioPtoPiNeg["AuAuc0010"]);
        divide(hProtonPt["AuAuc2030b"], hPionNegPt["AuAuc2030b"], RatioPtoPiNeg["AuAuc2030"]);
        divide(hProtonPt["AuAuc6092b"], hPionNegPt["AuAuc6092b"], RatioPtoPiNeg["AuAuc6092"]);
        divide(hProBarPt["AuAuc0010b"], hPionPt["AuAuc0010a"], RatioPBartoPion["AuAuc0010"]);
        divide(hProBarPt["AuAuc2030b"], hPionPt["AuAuc2030a"], RatioPBartoPion["AuAuc2030"]);
        divide(hProBarPt["AuAuc6092b"], hPionPt["AuAuc6092a"], RatioPBartoPion["AuAuc6092"]);
        
        //d03 p and p_bar yields
        hProtonPt["ptyieldsAuAuc0010"]->scaleW(1. / 955.4);
        hProtonPt["ptyieldsAuAuc0010"]->scaleW(1. / sow["AuAuc0010"]->sumW());
        
        hProtonPt["ptyieldsAuAuc2030"]->scaleW(1. / 373.8);
        hProtonPt["ptyieldsAuAuc2030"]->scaleW(1. / sow["AuAuc2030"]->sumW());
        
        hProtonPt["ptyieldsAuAuc4050"]->scaleW(1. / 120.3);
        hProtonPt["ptyieldsAuAuc4050"]->scaleW(1. / sow["AuAuc4050"]->sumW());
        
        hProtonPt["ptyieldsAuAuc6092"]->scaleW(1. / 14.5);
        hProtonPt["ptyieldsAuAuc6092"]->scaleW(1. / sow["AuAuc6092"]->sumW());
        
        hProBarPt["ptyieldsAuAuc0010"]->scaleW(1. / 955.4);
        hProtonPt["ptyieldsAuAuc0010"]->scaleW(1. / sow["AuAuc0010"]->sumW());
        
        hProBarPt["ptyieldsAuAuc2030"]->scaleW(1. / 373.8);
        hProtonPt["ptyieldsAuAuc2030"]->scaleW(1. / sow["AuAuc2030"]->sumW());
        
        hProBarPt["ptyieldsAuAuc4050"]->scaleW(1. / 120.3);
        hProtonPt["ptyieldsAuAuc4050"]->scaleW(1. / sow["AuAuc4050"]->sumW());
        
        hProBarPt["ptyieldsAuAuc6092"]->scaleW(1. / 14.5);
        hProtonPt["ptyieldsAuAuc6092"]->scaleW(1. / sow["AuAuc6092"]->sumW());
        
        //d04
        //0-10%
        //hPPlusPBarPt["ppluspbarAuAuc0010"]->scaleW(1. / 2.); //Scale by two to account for sum
        hPPlusPBarPt["ppluspbarAuAuc0010"]->scaleW(1. / (2. * 955.4)); //Scaling by N_coll and 2 for sum
        //60-92%
        //hPPlusPBarPt["ppluspbarAuAuc6092"]->scaleW(1. / 2.); //Scale by two to account for sum
        hPPlusPBarPt["ppluspbarAuAuc6092"]->scaleW(1. / (2. * 14.5)); //Scaling by N_coll and 2 for sum
        
        //Rcp
        divide(hPPlusPBarPt["ppluspbarAuAuc0010"], hPPlusPBarPt["ppluspbarAuAuc6092"], hRcp["ppluspbar"]);
        //divide(hPPlusPBarPt["ppluspbarAuAuc0010"]->scaleW(1. / 955.4), hPPlusPBarPt["ppluspbarAuAuc6092"]->scaleW(1. / 14.5), hRcp["ppluspbar"]);    >>>not successful when ran
        
        //d05
        //0-10%
        hPionPt["AuAuc0010b"]->scaleW(1. / 955.4);
        //60-92%
        hPionPt["AuAuc6092b"]->scaleW(1. / 14.5);
        
        //Rcp
        divide(hPionPt["AuAuc0010b"], hPionPt["AuAuc6092b"], hRcp["pion"]);
        
        //d06
        hChHadrons["AuAuc0010"]->scaleW(1. / 2.); //scale by two before divide
        divide(hChHadrons["AuAuc0010"], hPionPt["AuAuc0010c"], RatioHadtoPion["AuAuc0010"]);
        
        //d07
        hChHadrons["AuAuc6092"]->scaleW(1. / 2.); //scale by two before divide
        divide(hChHadrons["AuAuc6092"], hPionPt["AuAuc6092c"], RatioHadtoPion["AuAuc6092"]);
        
    }
    
    //Particles
    map<string, Histo1DPtr> hProtonPt; //p
    map<string, Histo1DPtr> hPionPosPt; //pi^+
    map<string, Histo1DPtr> hProBarPt; //p_bar
    map<string, Histo1DPtr> hPionNegPt; //pi^-
    map<string, Histo1DPtr> hPionPt; //pi^0
    map<string, Histo1DPtr> hPPlusPBarPt; //p + p_bar
    map<string, Histo1DPtr> hChHadrons; //p + p_bar + pi^+ + pi^-
    
    //Ratios
    map<string, Scatter2DPtr> RatioPtoPiPos;
    map<string, Scatter2DPtr> RatioPBartoPiNeg;
    map<string, Scatter2DPtr> RatioPtoPiNeg;
    map<string, Scatter2DPtr> RatioPBartoPion;
    map<string, Scatter2DPtr> RatioHadtoPion;
    
    //Rcp, Raa **TBD**
    map<string, Scatter2DPtr> hRcp;
    //map<string, Scatter2DPtr> hRaa;
    
    //Counter
    map<string, CounterPtr> sow;
    
    //Initialize collision system and AUAU centrality bins
    string beamOpt;
    enum CollisionSystem {AuAu200};
    CollisionSystem collSys;
    vector<int> AUAUCentralityBins{ 10, 20, 30, 40, 50, 60, 92};
    
    
    
  };
  DECLARE_RIVET_PLUGIN(PHENIX_2003_I619987);
}
