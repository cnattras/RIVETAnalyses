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
        

        // Initialise and register projections
        
        // Particles: pi^+, pi^-, pi^0, p, p_bar
        //pids (respectively): 211, -211, 111, 2212, -2212
        // all final-state particles within
        // the given eta acceptance
        /// Found the cuts on page 3, paragraph 3
        std::initializer_list<int> pdgIds = { 211, -211, 2212, -2212 };
        
        //Charged hadrons: protons, antiprotons, pions+, pions-, kaons+, kaons-
        std::initializer_list<int> chHIds = { 211, -211, 2212, -2212 , 321, -321};

        
        //Charged particles
        //consider adding back && Cuts::phi == 0.392
        const ALICE::PrimaryParticles cp(Cuts::absrap < 0.5 && Cuts::pT > 0.5*GeV && Cuts::pT < 9*GeV);
        declare(cp,"cp");
        //Inclusive charged hadrons
        const ALICE::PrimaryParticles ich(Cuts::absrap < 0.5 && Cuts::pT > 0.5*GeV && Cuts::pT < 9*GeV);
        declare(ich,"ich");
        //Neutral particles
        const UnstableParticles np(Cuts::absrap < 0.5 && Cuts::abspid == 111 && Cuts::pT > 0.5*GeV && Cuts::pT < 9*GeV);
        declare(np,"np");
        
        //Reading in the beam option
        beamOpt = getOption<string>("beam","NONE");
        if (beamOpt == "AUAU200") collSys = AuAu200;
        //In case "NONE" is given as option
        const ParticlePair& beam = beams();
        if (beamOpt == "NONE"){
            if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu200;
        }
        
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
        Particles inclusiveParticles = applyProjection<PrimaryParticles>(event,"ich").particles();
        Particles neutralParticles = applyProjection<UnstableParticles>(event,"np").particles();
        // All figures are for S_NN = 200 GeV collisions, so no if statement required for collSys
        if (collSys == AuAu200){
        /// Case for identified charged particles
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
                sow["AuAuc0010"]->fill();
                for (const Particle& p : chargedParticles)
                {
                    double partPt = p.pT() / GeV;
                    double pt_weight = 1. / (partPt * 2. * M_PI);
                    
                    switch(p.pid()) {
                        case 211: //pi^+
                        {
                            hPionPosPt["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case -211: //pi^-
                        {
                            hPionNegPt["AuAuc0010a"]->fill(partPt);
                            hPionNegPt["AuAuc0010b"]->fill(partPt);
                            break;
                        }
                        case 2212: //p
                        {
                            hProtonPt["AuAuc0010a"]->fill(partPt);
                            hProtonPt["AuAuc0010b"]->fill(partPt);
                            hProtonPt["ptyieldsAuAuc0010"]->fill(partPt, pt_weight);
                            hPPlusPBarPt["ppluspbarAuAuc0010"]->fill(partPt);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            hProBarPt["AuAuc0010a"]->fill(partPt);
                            hProBarPt["AuAuc0010b"]->fill(partPt);
                            hProBarPt["ptyieldsAuAuc0010"]->fill(partPt, pt_weight);
                            hPPlusPBarPt["ppluspbarAuAuc0010"]->fill(partPt);
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
                sow["AuAuc2030"]->fill();
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
                sow["AuAuc4050"]->fill();
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
                sow["AuAuc6092"]->fill();
                for (const Particle& p : chargedParticles)
                {
                    double partPt = p.pT() / GeV;
                    double pt_weight = 1. / (partPt * 2. * M_PI);
                    
                    switch(p.pid()) {
                        case 211: //pi^+
                        {
                            hPionPosPt["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case -211: //pi^-
                        {
                            hPionNegPt["AuAuc6092a"]->fill(partPt);
                            hPionNegPt["AuAuc6092b"]->fill(partPt);
                            break;
                        }
                        case 2212: //p
                        {
                            hProtonPt["AuAuc6092a"]->fill(partPt);
                            hProtonPt["AuAuc6092b"]->fill(partPt);
                            hProtonPt["ptyieldsAuAuc6092"]->fill(partPt, pt_weight);
                            hPPlusPBarPt["ppluspbarAuAuc6092"]->fill(partPt);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            hProBarPt["AuAuc6092a"]->fill(partPt);
                            hProBarPt["AuAuc6092b"]->fill(partPt);
                            hProBarPt["ptyieldsAuAuc6092"]->fill(partPt, pt_weight);
                            hPPlusPBarPt["ppluspbarAuAuc6092"]->fill(partPt);
                            break;
                        }
                    }
                    
                }
            }
            
            else{
                break;
            }
        }
        

        for(Particle p : inclusiveParticles)
        {
            // We will need to write for centrality
            const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
            const double c = cent();

            // events with centrality < 0 or > 92 are invalid. We use vetoEvent.
            if (c < 0. || c > 92.) vetoEvent;

            // 0-10% centrality
            if ((c > 0.) && (c < 10.))
            {
                for (const Particle& p : inclusiveParticles)
                {
                    double partPt = p.pT() / GeV;
                    double pt_weight = 1. / (partPt * 2. * M_PI);

                    switch(p.pid()) {
                        case 211: //pi^+
                        {
                            hChHadrons["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case -211: //pi^-
                        {
                            hChHadrons["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case 2212: //p
                        {
                            hChHadrons["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            hChHadrons["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case 321: //kaon^+
                        {
                            hChHadrons["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case -321: //kaon^-
                        {
                            hChHadrons["AuAuc0010"]->fill(partPt);
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

            else if ((c >= 30.) && (c < 40.))
            {
                break;
            }

            else if ((c >= 40.) && (c < 50.))
            {
                break;
            }

            else if ((c >= 50.) && (c < 60.)){
                break;
            }

            else if ((c >= 60.) && (c < 92.)){
                for (const Particle& p : inclusiveParticles)
                {
                    double partPt = p.pT() / GeV;
                    double pt_weight = 1. / (partPt * 2. * M_PI);

                    switch(p.pid()) {
                        case 211: //pi^+
                        {
                            hChHadrons["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case -211: //pi^-
                        {
                            hChHadrons["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case 2212: //p
                        {
                            hChHadrons["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case -2212: //p_bar
                        {
                            hChHadrons["AuAuc6092"]->fill(partPt);
                            break;
                        }
                        case 321: //kaon^+
                        {
                            hChHadrons["AuAuc0010"]->fill(partPt);
                            break;
                        }
                        case -321: //kaon^-
                        {
                            hChHadrons["AuAuc0010"]->fill(partPt);
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
        
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
        
        //d01: p and p_bar ratios
        //divide(hProtonPt["AuAuc0010a"], hPionPosPt["AuAuc0010"], RatioPtoPiPos["AuAuc0010"]);
        
        /*binShift(*hProtonPt["AuAuc0010a"]);
        binShift(*hPionPosPt["AuAuc0010"]);
        divide(hProtonPt["AuAuc0010a"], hPionPosPt["AuAuc0010"], RatioPtoPiPos["AuAuc0010"]);

        binShift(*hProtonPt["AuAuc2030a"]);
        binShift(*hPionPosPt["AuAuc2030"]);
        divide(hProtonPt["AuAuc2030a"], hPionPosPt["AuAuc2030"], RatioPtoPiPos["AuAuc2030"]);

        binShift(*hProtonPt["AuAuc6092a"]);
        binShift(*hPionPosPt["AuAuc6092"]);
        divide(hProtonPt["AuAuc6092a"], hPionPosPt["AuAuc6092"], RatioPtoPiPos["AuAuc6092"]);

        binShift(*hProBarPt["AuAuc0010a"]);
        binShift(*hPionNegPt["AuAuc0010a"]);
        divide(hProBarPt["AuAuc0010a"], hPionNegPt["AuAuc0010a"], RatioPBartoPiNeg["AuAuc0010"]);

        binShift(*hProBarPt["AuAuc2030a"]);
        binShift(*hPionNegPt["AuAuc2030a"]);
        divide(hProBarPt["AuAuc2030a"], hPionNegPt["AuAuc2030a"], RatioPBartoPiNeg["AuAuc2030"]);

        binShift(*hProBarPt["AuAuc6092a"]);
        binShift(*hPionNegPt["AuAuc6092a"]);
        divide(hProBarPt["AuAuc6092a"], hPionNegPt["AuAuc6092a"], RatioPBartoPiNeg["AuAuc6092"]);

        //d02 p and p_bar ratios
        binShift(*hProtonPt["AuAuc0010b"]);
        binShift(*hPionNegPt["AuAuc0010b"]);
        divide(hProtonPt["AuAuc0010b"], hPionNegPt["AuAuc0010b"], RatioPtoPiNeg["AuAuc0010"]);

        binShift(*hProtonPt["AuAuc2030b"]);
        binShift(*hPionNegPt["AuAuc2030b"]);
        divide(hProtonPt["AuAuc2030b"], hPionNegPt["AuAuc2030b"], RatioPtoPiNeg["AuAuc2030"]);

        binShift(*hProtonPt["AuAuc6092b"]);
        binShift(*hPionNegPt["AuAuc6092b"]);
        divide(hProtonPt["AuAuc6092b"], hPionNegPt["AuAuc6092b"], RatioPtoPiNeg["AuAuc6092"]);

        binShift(*hProBarPt["AuAuc0010b"]);
        binShift(*hPionPt["AuAuc0010a"]);
        divide(hProBarPt["AuAuc0010b"], hPionPt["AuAuc0010a"], RatioPBartoPion["AuAuc0010"]);

        binShift(*hProBarPt["AuAuc2030b"]);
        binShift(*hPionPt["AuAuc2030a"]);
        divide(hProBarPt["AuAuc2030b"], hPionPt["AuAuc2030a"], RatioPBartoPion["AuAuc2030"]);

        binShift(*hProBarPt["AuAuc6092b"]);
        binShift(*hPionPt["AuAuc6092a"]);
        divide(hProBarPt["AuAuc6092b"], hPionPt["AuAuc6092a"], RatioPBartoPion["AuAuc6092"]);*/
        
        //d03 p and p_bar yields
        //first scalling is commented out because this is actually scaling by the N_coll provided in the paper
        //sow is our N_coll in this simulation
        
        if(sow["AuAuc0010"]->sumW() != 0){
            hProtonPt["ptyieldsAuAuc0010"]->scaleW(1. / 955.4);
            hProtonPt["ptyieldsAuAuc0010"]->scaleW(1. / sow["AuAuc0010"]->sumW());
            binShift(*hProtonPt["ptyieldsAuAuc0010"]);
        }else{
            std::cerr << "Error: Divide by zero encountered. Unable to scale histogram." << std::endl;
        }

        if(sow["AuAuc2030"]->sumW() != 0){
            hProtonPt["ptyieldsAuAuc2030"]->scaleW(1. / 373.8);
            hProtonPt["ptyieldsAuAuc2030"]->scaleW(1. / sow["AuAuc2030"]->sumW());
        binShift(*hProtonPt["ptyieldsAuAuc2030"]);
        }else{
            std::cerr << "Error: Divide by zero encountered. Unable to scale histogram." << std::endl;
        }

        if(sow["AuAuc4050"]->sumW() != 0){
            hProtonPt["ptyieldsAuAuc4050"]->scaleW(1. / 120.3);
            hProtonPt["ptyieldsAuAuc4050"]->scaleW(1. / sow["AuAuc4050"]->sumW());
            binShift(*hProtonPt["ptyieldsAuAuc4050"]);
        }else{
            std::cerr << "Error: Divide by zero encountered. Unable to scale histogram." << std::endl;
        }
        
        if(sow["AuAuc6092"]->sumW() != 0){
            hProtonPt["ptyieldsAuAuc6092"]->scaleW(1. / 14.5);
            hProtonPt["ptyieldsAuAuc6092"]->scaleW(1. / sow["AuAuc6092"]->sumW());
            binShift(*hProtonPt["ptyieldsAuAuc6092"]);
        }else{
            std::cerr << "Error: Divide by zero encountered. Unable to scale histogram." << std::endl;
        }
        
        if(sow["AuAuc0010"]->sumW() != 0){
            hProBarPt["ptyieldsAuAuc0010"]->scaleW(1. / 955.4);
            hProBarPt["ptyieldsAuAuc0010"]->scaleW(1. / sow["AuAuc0010"]->sumW());
            binShift(*hProBarPt["ptyieldsAuAuc0010"]);
        }else{
            std::cerr << "Error: Divide by zero encountered. Unable to scale histogram." << std::endl;
        }
        
        if(sow["AuAuc2030"]->sumW() != 0){
            hProBarPt["ptyieldsAuAuc2030"]->scaleW(1. / 373.8);
            hProBarPt["ptyieldsAuAuc2030"]->scaleW(1. / sow["AuAuc2030"]->sumW());
            binShift(*hProBarPt["ptyieldsAuAuc2030"]);
        }else{
            std::cerr << "Error: Divide by zero encountered. Unable to scale histogram." << std::endl;
        }
        
        if(sow["AuAuc4050"]->sumW() != 0){
            hProBarPt["ptyieldsAuAuc4050"]->scaleW(1. / 120.3);
            hProBarPt["ptyieldsAuAuc4050"]->scaleW(1. / sow["AuAuc4050"]->sumW());
            binShift(*hProBarPt["ptyieldsAuAuc4050"]);
        }else{
            std::cerr << "Error: Divide by zero encountered. Unable to scale histogram." << std::endl;
        }
        
        if(sow["AuAuc6092"]->sumW() != 0){
            hProBarPt["ptyieldsAuAuc6092"]->scaleW(1. / 14.5);
            hProBarPt["ptyieldsAuAuc6092"]->scaleW(1. / sow["AuAuc6092"]->sumW());
            binShift(*hProBarPt["ptyieldsAuAuc6092"]);
        }else{
            std::cerr << "Error: Divide by zero encountered. Unable to scale histogram." << std::endl;
        }
        
        //d04
        //0-10%
        //hPPlusPBarPt["ppluspbarAuAuc0010"]->scaleW(1. / 2.); //Scale by two to account for sum
        //hPPlusPBarPt["ppluspbarAuAuc0010"]->scaleW(1. / 955.4); //Scaling by N_coll and 2 for sum
        //hPPlusPBarPt["ppluspbarAuAuc0010"]->scaleW(1. / 2.);
        if(sow["AuAuc0010"]->sumW() != 0){
            hPPlusPBarPt["ppluspbarAuAuc0010"]->scaleW(1. / (2. * 955.4));
            hPPlusPBarPt["ppluspbarAuAuc0010"]->scaleW(1. / sow["AuAuc0010"]->sumW()); 
            binShift(*hPPlusPBarPt["ppluspbarAuAuc0010"]);
        }else{
            std::cerr << "Error: Divide by zero encountered. Unable to scale histogram." << std::endl;
        }
        //60-92%
        //hPPlusPBarPt["ppluspbarAuAuc6092"]->scaleW(1. / 2.); //Scale by two to account for sum
        //hPPlusPBarPt["ppluspbarAuAuc6092"]->scaleW(1. / 14.5); //Scaling by N_coll and 2 for sum
        //hPPlusPBarPt["ppluspbarAuAuc6092"]->scaleW(1. / 2.);
        if(sow["AuAuc6092"]->sumW() != 0){
            hPPlusPBarPt["ppluspbarAuAuc6092"]->scaleW(1. / (2. * 14.5));
            hPPlusPBarPt["ppluspbarAuAuc6092"]->scaleW(1. / sow["AuAuc6092"]->sumW());
            binShift(*hPPlusPBarPt["ppluspbarAuAuc6092"]);
        }else{
            std::cerr << "Error: Divide by zero encountered. Unable to scale histogram." << std::endl;
        }

        //Rcp
        binShift(*hPPlusPBarPt["ppluspbarAuAuc0010"]);
        binShift(*hPPlusPBarPt["ppluspbarAuAuc6092"]);
        divide(hPPlusPBarPt["ppluspbarAuAuc0010"], hPPlusPBarPt["ppluspbarAuAuc6092"], hRcp["ppluspbar"]);
        //divide(hPPlusPBarPt["ppluspbarAuAuc0010"]->scaleW(1. / 955.4), hPPlusPBarPt["ppluspbarAuAuc6092"]->scaleW(1. / 14.5), hRcp["ppluspbar"]);    >>>not successful when ran
        
        //d05
        //0-10%
        hPionPt["AuAuc0010b"]->scaleW(1. / 955.4);
        hPionPt["AuAuc0010b"]->scaleW(1. / sow["AuAuc0010"]->sumW());
        binShift(*hPionPt["AuAuc0010b"]);
        //60-92%
        hPionPt["AuAuc6092b"]->scaleW(1. / 14.5);
        hPionPt["AuAuc6092b"]->scaleW(1. / sow["AuAuc6092"]->sumW());
        binShift(*hPionPt["AuAuc6092b"]);
        
        //Rcp
        binShift(*hPionPt["AuAuc0010b"]);
        binShift(*hPionPt["AuAuc6092b"]);
        divide(hPionPt["AuAuc0010b"], hPionPt["AuAuc6092b"], hRcp["pion"]);
        
        //d06
        hChHadrons["AuAuc0010"]->scaleW(1. / 2.); //scale by two before divide
        binShift(*hChHadrons["AuAuc0010"]);
        binShift(*hPionPt["AuAuc0010c"]);
        divide(hChHadrons["AuAuc0010"], hPionPt["AuAuc0010c"], RatioHadtoPion["AuAuc0010"]);
        
        //d07
        hChHadrons["AuAuc6092"]->scaleW(1. / 2.); //scale by two before divide
        binShift(*hChHadrons["AuAuc6092"]);
        binShift(*hPionPt["AuAuc6092c"]);
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
    enum CollisionSystem {AuAu200};
    CollisionSystem collSys;
    string beamOpt = "NONE";
    vector<int> AUAUCentralityBins{ 10, 20, 30, 40, 50, 60, 92};
    
    
    
  };
  DECLARE_RIVET_PLUGIN(PHENIX_2003_I619987);
}
