// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
//#include "RHICCentrality.hh" //external header for Centrality calculation
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES

namespace Rivet {

  class CMS_2020_I064906 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2020_I064906);
    
    /// Function to get bin center to weight entries with
        bool getDeltaPt(YODA::Histo1D hist, double pT, double &deltaPt) {
        if(pT > hist.xMin() && pT < hist.xMax())
        {
            deltaPt = hist.bin(hist.binIndexAt(pT)).xMid();

            return true;
        }
        else return false;
    }

    /// Book histograms and initialise projections before the run
    void init() {
    
//      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M","V0M");

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within the given eta acceptance
      //Particles: K0S, Lambda, Xi-, Anti-Xi+, Omega-, Anti-Omega+, anti-Lambda(?)
    	std::initializer_list<int> pdgIds = {310, 3122, 3312,-3312, 3334, -3334};
    	const PrimaryParticles fs(pdgIds, Cuts::abscharge > 0 && Cuts::absrap < 1.8);
    	declare(fs, "fs");
    	const PrimaryParticles ns(pdgIds, Cuts::abscharge == 0 && Cuts::absrap < 1.8);
    	declare(ns, "ns");
    	
    	beamOpt = getOption<string>("beam", "NONE");

	if (beamOpt == "PP") collSys = pp;
	else if (beamOpt == "pPB") collSys = pPB;
	
	//Create various counters
	book(sow["sow_pp"], "sow_pp");
	book(sow["sow_pPB"], "sow_pPB");

        //--------Begin booking histograms-----------
    	
    	//Inv. pT of K0s for various y_CM (fig 2.1)
	book(hInvariantPTK0S["pT_K0S_pp_-1.8<yCM<0"], 1, 1, 1);
	book(hInvariantPTK0S["pT_K0S_pp_-1.8<yCM<1.8"], 1, 1, 2);
	book(hInvariantPTK0S["pT_K0S_pp_0<yCM<1.8"], 1, 1, 3);
	book(hInvariantPTK0S["pT_K0S_pPB_-1.8<yCM<0"], 1, 1, 4);
	book(hInvariantPTK0S["pT_K0S_pPB_-1.8<yCM<1.8"], 1, 1, 5);
	book(hInvariantPTK0S["pT_K0S_pPB_0<yCM<1.8"], 1, 1, 6);
	
	//Inv. pT of Lambda for various y_CM (fig 2.2)
	book(hInvariantPTLambda["pT_Lambda_pp_-1.8<yCM<0"], 2, 1, 1);
	book(hInvariantPTLambda["pT_Lambda_pp_-1.8<yCM<1.8"], 2, 1, 2);
	book(hInvariantPTLambda["pT_Lambda_pp_0<yCM<1.8"], 2, 1, 3);
	book(hInvariantPTLambda["pT_Lambda_pPB_-1.8<yCM<0"], 2, 1, 4);
	book(hInvariantPTLambda["pT_Lambda_pPB_-1.8<yCM<1.8"], 2, 1, 5);
	book(hInvariantPTLambda["pT_Lambda_pPB_0<yCM<1.8"], 2, 1, 6);
	
	//Inv. pT of Xi for various y_CM (fig 2.3)
	book(hInvariantPTXi["pT_Xi_pp_-1.8<yCM<0"], 3, 1, 1);
	book(hInvariantPTXi["pT_Xi_pp_-1.8<yCM<1.8"], 3, 1, 2);
	book(hInvariantPTXi["pT_Xi_pp_0<yCM<1.8"], 3, 1, 3);
	book(hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<0"], 3, 1, 4);
	book(hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<1.8"], 3, 1, 5);
	book(hInvariantPTXi["pT_Xi_pPB_0<yCM<1.8"], 3, 1, 6);
	
	//Inv. pT of Omega for various y_CM (fig 2.4)
	book(hInvariantPTOmega["pT_Omega_pp_-1.8<yCM<1.8"], 4, 1, 1);
	book(hInvariantPTOmega["pT_Omega_pPB_-1.8<yCM<1.8"], 4, 1, 2);
	
	//R_pPB for -1.8<y_CM<1.8 in pPB (fig 3)
	//book(hRpPBFullyCM["RpPB_K0S_-1.8<yCM<1.8"], 5, 1, 1);
	//book(hRpPBFullyCM["RpPB_Lambda_-1.8<yCM<1.8"], 5, 1, 2);
	//book(hRpPBFullyCM["RpPB_Xi_-1.8<yCM<1.8"], 5, 1, 3);
	//book(hRpPBFullyCM["RpPB_Omega_-1.8<yCM<1.8"], 5, 1, 4);
	
	//R_pPB for -1.8<y_CM<0 in pPB (fig 4.1)
	//book(hRpPBLowyCM["RpPB_K0S_-1.8<yCM<0"], 6, 1, 1);
	//book(hRpPBLowyCM["RpPB_Lambda_-1.8<yCM<0"], 6, 1, 2);
	//book(hRpPBLowyCM["RpPB_Xi_-1.8<yCM<0"], 6, 1, 3);
	
	//R_pPB for 0<y_CM<1.8 in pPB (fig 4.2)
	//book(hRpPBHighyCM["RpPB_K0S_0<yCM<1.8"], 7, 1, 1);
	//book(hRpPBHighyCM["RpPB_Lambda_0<yCM<1.8"], 7, 1, 2);
	//book(hRpPBHighyCM["RpPB_Xi_0<yCM<1.8"], 7, 1, 3);
	
	//Inv. pT of K0S for various y_CM in pPB (fig 5.1)
	book(hInvariantPTK0SpPB["pT_K0S_pPB_-1.8<yCM<-1.3"], 8, 1, 1);
	book(hInvariantPTK0SpPB["pT_K0S_pPB_-1.3<yCM<-0.8"], 8, 1, 2);
	book(hInvariantPTK0SpPB["pT_K0S_pPB_-0.8<yCM<-0.3"], 8, 1, 3);
	book(hInvariantPTK0SpPB["pT_K0S_pPB_0.3<yCM<0.8"], 8, 1, 4);
	book(hInvariantPTK0SpPB["pT_K0S_pPB_0.8<yCM<1.3"], 8, 1, 5);
	book(hInvariantPTK0SpPB["pT_K0S_pPB_1.3<yCM<1.8"], 8, 1, 6);
	
	//Inv. pT of Lambda for various y_CM in pPB (fig 5.2)
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_-1.8<yCM<-1.3"], 9, 1, 1);
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_-1.3<yCM<-0.8"], 9, 1, 2);
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_-0.8<yCM<-0.3"], 9, 1, 3);
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_0.3<yCM<0.8"], 9, 1, 4);
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_0.8<yCM<1.3"], 9, 1, 5);
	book(hInvariantPTLambdapPB["pT_Lambda_pPB_1.3<yCM<1.8"], 9, 1, 6);
	
//	book(h["K0S_yCM_low"], 10, 1, 1);
//	book(h["Lambda_yCM_low"], 11, 1, 1);
	
//	book(h["charged_yCM_low"], 10, 1, 3);
//	book(h["charged_yCM_mid"], 11, 1, 3);
//	book(h["charged_yCM_high"], 12, 1, 3);
	
	//-----Book Scatter Plots from division-------
	
	//R_pPB for -1.8<y_CM<1.8 in pPB (fig 3)
	string refname1 = mkAxisCode(5, 1, 1);
	const Scatter2D& refdata1 = refData(refname1);
	book(hInvariantPTK0S["pT_K0S_pPB_full"], refname1 + "_-1.8<yCM<1.8_pPB", refdata1);
	book(hInvariantPTK0S["pT_K0S_pp_full"], refname1 + "_-1.8<yCM<1.8_pp", refdata1);
	book(RpPBFullyCM["pPB_pT_K0S/pp_pT_K0S_full"], refname1);
	
	string refname2 = mkAxisCode(5, 1, 2);
	const Scatter2D& refdata2 = refData(refname2);
	book(hInvariantPTLambda["pT_Lambda_pPB_full"], refname2 + "_-1.8<yCM<1.8_pPB", refdata2);
	book(hInvariantPTLambda["pT_Lambda_pp_full"], refname2 + "_-1.8<yCM<1.8_pp", refdata2);
	book(RpPBFullyCM["pPB_pT_Lambda/pp_pT_Lambda_full"], refname2);
	
	string refname3 = mkAxisCode(5, 1, 3);
	const Scatter2D& refdata3 = refData(refname3);
	book(hInvariantPTXi["pT_Xi_pPB_full"], refname3 + "_-1.8<yCM<1.8_pPB", refdata3);
	book(hInvariantPTXi["pT_Xi_pp_full"], refname3 + "_-1.8<yCM<1.8_pp", refdata3);
	book(RpPBFullyCM["pPB_pT_Xi/pp_pT_Xi_full"], refname3);
	
	string refname4 = mkAxisCode(5, 1, 4);
	const Scatter2D& refdata4 = refData(refname4);
	book(hInvariantPTOmega["pT_Omega_pPB_full"], refname4 + "_-1.8<yCM<1.8_pPB", refdata4);
	book(hInvariantPTOmega["pT_Omega_pp_full"], refname4 + "_-1.8<yCM<1.8_pp", refdata4);
	book(RpPBFullyCM["pPB_pT_Omega/pp_pT_Omega_full"], refname4);
	
	//R_pPB for -1.8<y_CM<0 in pPB (fig 4.1)
	string refname5 = mkAxisCode(6, 1, 1);
	const Scatter2D& refdata5 = refData(refname5);
	book(hInvariantPTK0S["pT_K0S_pPB_low"], refname5 + "_-1.8<yCM<0_pPB", refdata5);
	book(hInvariantPTK0S["pT_K0S_pp_low"], refname5 + "_-1.8<yCM<0_pp", refdata5);
	book(RpPBLowyCM["pPB_pT_K0S/pp_pT_K0S_low"], refname5);
	
	string refname6 = mkAxisCode(6, 1, 2);
	const Scatter2D& refdata6 = refData(refname6);
	book(hInvariantPTLambda["pT_Lambda_pPB_low"], refname6 + "_-1.8<yCM<0_pPB", refdata6);
	book(hInvariantPTLambda["pT_Lambda_pp_low"], refname6 + "_-1.8<yCM<0_pp", refdata6);
	book(RpPBLowyCM["pPB_pT_Lambda/pp_pT_Lambda_low"], refname6);
	
	string refname7 = mkAxisCode(6, 1, 3);
	const Scatter2D& refdata7 = refData(refname7);
	book(hInvariantPTXi["pT_Xi_pPB_low"], refname7 + "_-1.8<yCM<0_pPB", refdata7);
	book(hInvariantPTXi["pT_Xi_pp_low"], refname7 + "_-1.8<yCM<0_pp", refdata7);
	book(RpPBLowyCM["pPB_pT_Xi/pp_pT_Xi_low"], refname7);
	
	//R_pPB for 0<y_CM<1.8 in pPB (fig 4.2)
	string refname8 = mkAxisCode(7, 1, 1);
	const Scatter2D& refdata8 = refData(refname8);
	book(hInvariantPTK0S["pT_K0S_pPB_high"], refname8 + "_0<yCM<1.8_pPB", refdata8);
	book(hInvariantPTK0S["pT_K0S_pp_high"], refname8 + "_0<yCM<1.8_pp", refdata8);
	book(RpPBHighyCM["pPB_pT_K0S/pp_pT_K0S_high"], refname8);
	
	string refname9 = mkAxisCode(7, 1, 2);
	const Scatter2D& refdata9 = refData(refname9);
	book(hInvariantPTLambda["pT_Lambda_pPB_high"], refname9 + "_0<yCM<1.8_pPB", refdata9);
	book(hInvariantPTLambda["pT_Lambda_pp_high"], refname9 + "_0<yCM<1.8_pp", refdata9);
	book(RpPBHighyCM["pPB_pT_Lambda/pp_pT_Lambda_high"], refname9);
	
	string refname10 = mkAxisCode(7, 1, 3);
	const Scatter2D& refdata10 = refData(refname10);
	book(hInvariantPTXi["pT_Xi_pPB_high"], refname10 + "_0<yCM<1.8_pPB", refdata10);
	book(hInvariantPTXi["pT_Xi_pp_high"], refname10 + "_0<yCM<1.8_pp", refdata10);
	book(RpPBHighyCM["pPB_pT_Xi/pp_pT_Xi_high"], refname10);
	
	//Y_asym for 0.3 < |yCM| < 0.8
	string refname11 = mkAxisCode(10, 1, 1);
	const Scatter2D& refdata11 = refData(refname11);
	book(hInvariantPTK0SpPB["K0S_pT_pPB_LowNeg"], refname11 + "_-0.8<yCM<-0.3", refdata11);
	book(hInvariantPTK0SpPB["K0S_pT_pPB_LowPos"], refname11 + "_0.3<yCM<0.8", refdata11);
	book(YasymLow["K0S_-0.8<yCM<-0.3/0.3<yCM<0.8"], refname11);
	
	string refname12 = mkAxisCode(10, 1, 2);
	const Scatter2D& refdata12 = refData(refname12);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_LowNeg"], refname12 + "_-0.8<yCM<-0.3", refdata12);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_LowPos"], refname12 + "_0.3<yCM<0.8", refdata12);
	book(YasymLow["Lambda_-0.8<yCM<-0.3/0.3<yCM<0.8"], refname12);
	
	//Y_asym Low for h+/-
	string refname17 = mkAxisCode(10, 1, 3);
	const Scatter2D& refdata17 = refData(refname17);
	book(h["negative_charged_yCM_low"], refname17 + "_-0.8<yCM<-0.3", refdata17);
	book(h["positive_charged_yCM_low"], refname17 + "_0.3<yCM<0.8", refdata17);
	book(YasymLow["h+/-_-0.8<yCM<-0.3/0.3<yCM<0.8"], refname17);
	
	//Y_asym for 0.8 < |yCM| < 1.3
	string refname13 = mkAxisCode(11, 1, 1);
	const Scatter2D& refdata13 = refData(refname13);
	book(hInvariantPTK0SpPB["K0S_pT_pPB_MidNeg"], refname13 + "_-1.3<yCM<-0.8", refdata13);
	book(hInvariantPTK0SpPB["K0S_pT_pPB_MidPos"], refname13 + "_0.8<yCM<1.3", refdata13);
	book(YasymMid["K0S_-1.3<yCM<-0.8/0.8<yCM<1.3"], refname13);
	
	string refname14 = mkAxisCode(11, 1, 2);
	const Scatter2D& refdata14 = refData(refname14);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_MidNeg"], refname14 + "_-1.3<yCM<-0.8", refdata14);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_MidPos"], refname14 + "_0.8<yCM<1.3", refdata14);
	book(YasymMid["Lambda_-1.3<yCM<-0.8/0.8<yCM<1.3"], refname14);
	
	//Y_asym Mid for h+/-
	string refname18 = mkAxisCode(11, 1, 3);
	const Scatter2D& refdata18 = refData(refname18);
	book(h["negative_charged_yCM_mid"], refname18 + "_-1.3<yCM<-0.8", refdata18);
	book(h["positive_charged_yCM_mid"], refname18 + "_0.8<yCM<1.3", refdata18);
	book(YasymMid["h+/-_-1.3<yCM<-0.8/0.8<yCM<1.3"], refname18);
	
	//Y_asym for 1.3 < |yCM| < 1.8
	string refname15 = mkAxisCode(12, 1, 1);
	const Scatter2D& refdata15 = refData(refname15);
	book(hInvariantPTK0SpPB["K0S_pT_pPB_HighNeg"], refname15 + "_-1.8<yCM<-1.3", refdata15);
	book(hInvariantPTK0SpPB["K0S_pT_pPB_HighPos"], refname15 + "_1.3<yCM<1.8", refdata15);
	book(YasymHigh["K0S_-1.8<yCM<-1.3/1.3<yCM<1.8"], refname15);
	
	string refname16 = mkAxisCode(12, 1, 2);
	const Scatter2D& refdata16 = refData(refname16);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_HighNeg"], refname16 + "_-1.8<yCM<-1.3", refdata16);
	book(hInvariantPTLambdapPB["Lambda_pT_pPB_HighPos"], refname16 + "_1.3<yCM<1.8", refdata16);
	book(YasymHigh["Lambda_-1.8<yCM<-1.3/1.3<yCM<1.8"], refname16);
	
	//Y_asym High for h+/-
	string refname19 = mkAxisCode(12, 1, 3);
	const Scatter2D& refdata19 = refData(refname19);
	book(h["negative_charged_yCM_high"], refname19 + "_-1.8<yCM<-1.3", refdata19);
	book(h["positive_charged_yCM_high"], refname19 + "_1.3<yCM<1.8", refdata19);
	book(YasymHigh["h+/-_-1.8<yCM<-1.3/1.3<yCM<1.8"], refname19);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles chargedParticles = applyProjection<PrimaryParticles>(event, "fs").particles();
      Particles neutralParticles = applyProjection<PrimaryParticles>(event, "ns").particles();
      
//      // The centrality projection.
//      const CentralityProjection& centProj = apply<CentralityProjection>(event,"V0M");
//      // The centrality.
//      const double cent = centProj();
//      // Veto event for too large centralities since those are not used in the analysis at all
//      if ((cent < 0.) || (cent > 90.)) vetoEvent;
    
      if (collSys == pp)
	{
		sow["sow_pp"]->fill();
		for (Particle p : neutralParticles)
		{
			double partPt = p.pT() / GeV;
			double pt_weight = 1. / (partPt * 2. * M_PI);
			double deltaPt = 0.;
			switch (p.pid()) {
			case 310: // K0S
			{
				if(p.rap() < 0) hInvariantPTK0S["pT_K0S_pp_-1.8<yCM<0"]->fill(partPt, pt_weight);
				
				//if(getDeltaPt(*hInvariantPTK0S->second, partPt, deltaPt)) {
                        	//pt_weight /= deltaPt;
                        	//hInvariantPTK0S->second->fill(partPt, pt_weight);
                    		//}
                    			
				hInvariantPTK0S["pT_K0S_pp_-1.8<yCM<1.8"]->fill(partPt, pt_weight);
				if(p.rap() > 0) hInvariantPTK0S["pT_K0S_pp_0<yCM<1.8"]->fill(partPt, pt_weight);
				
				hInvariantPTK0S["pT_K0S_pp_full"]->fill(partPt, pt_weight);
				if(p.rap() < 0) hInvariantPTK0S["pT_K0S_pp_low"]->fill(partPt, pt_weight);
				if(p.rap() > 0) hInvariantPTK0S["pT_K0S_pp_high"]->fill(partPt, pt_weight);
				
				break;
			}
			case 3122: // Lambda
			{
				if(p.rap() < 0) hInvariantPTLambda["pT_Lambda_pp_-1.8<yCM<0"]->fill(partPt, pt_weight);
				hInvariantPTLambda["pT_Lambda_pp_-1.8<yCM<1.8"]->fill(partPt, pt_weight);
				if(p.rap() > 0) hInvariantPTLambda["pT_Lambda_pp_0<yCM<1.8"]->fill(partPt, pt_weight);
				
				hInvariantPTLambda["pT_Lambda_pp_full"]->fill(partPt, pt_weight);
				if(p.rap() < 0) hInvariantPTLambda["pT_Lambda_pp_low"]->fill(partPt, pt_weight);
				if(p.rap() > 0) hInvariantPTLambda["pT_Lambda_pp_high"]->fill(partPt, pt_weight);
				
				break;
			}
		}
		}
		
		for (Particle p : chargedParticles)
		{
			double partPt = p.pT() / GeV;
			double pt_weight = 1. / (partPt * 2. * M_PI);
			switch (p.pid()) {
			case 3312: // Xi-
			{
				if(p.rap() < 0) hInvariantPTXi["pT_Xi_pp_-1.8<yCM<0"]->fill(partPt, pt_weight);
				hInvariantPTXi["pT_Xi_pp_-1.8<yCM<1.8"]->fill(partPt, pt_weight);
				if(p.rap() > 0) hInvariantPTXi["pT_Xi_pp_0<yCM<1.8"]->fill(partPt, pt_weight);
				
				hInvariantPTXi["pT_Xi_pp_full"]->fill(partPt, pt_weight);
				if(p.rap() < 0) hInvariantPTXi["pT_Xi_pp_low"]->fill(partPt, pt_weight);
				if(p.rap() > 0) hInvariantPTXi["pT_Xi_pp_high"]->fill(partPt, pt_weight);
				
				break;
			}
			case -3312: // Xi+
			{
				if(p.rap() < 0) hInvariantPTXi["pT_Xi_pp_-1.8<yCM<0"]->fill(partPt, pt_weight);
				hInvariantPTXi["pT_Xi_pp_-1.8<yCM<1.8"]->fill(partPt, pt_weight);
				if(p.rap() > 0) hInvariantPTXi["pT_Xi_pp_0<yCM<1.8"]->fill(partPt, pt_weight);
				
				hInvariantPTXi["pT_Xi_pp_full"]->fill(partPt, pt_weight);
				if(p.rap() < 0) hInvariantPTXi["pT_Xi_pp_low"]->fill(partPt, pt_weight);
				if(p.rap() > 0) hInvariantPTXi["pT_Xi_pp_high"]->fill(partPt, pt_weight);
				
				break;
			}
			case 3334: // Omega-
			{
				hInvariantPTOmega["pT_Omega_pp_-1.8<yCM<1.8"]->fill(partPt, pt_weight);
				
				hInvariantPTOmega["pT_Omega_pp_full"]->fill(partPt, pt_weight);
				
				break;
			}
			case -3334: // Omega+
			{
				hInvariantPTOmega["pT_Omega_pp_-1.8<yCM<1.8"]->fill(partPt, pt_weight);
				
				hInvariantPTOmega["pT_Omega_pp_full"]->fill(partPt, pt_weight);
				
				break;
			}
		}
	}
	return;
	}
	
	if (collSys == pPB)
	{
		sow["sow_pPB"]->fill();
		for (Particle p : neutralParticles)
		{
			double partPt = p.pT() / GeV;
			double pt_weight = 1. / (partPt * 2. * M_PI);
			switch (p.pid()) {
			case 310: // K0S
			{
				if(p.rap() < 0)hInvariantPTK0S["pT_K0S_pPB_-1.8<yCM<0"]->fill(partPt, pt_weight);
				hInvariantPTK0S["pT_K0S_pPB_-1.8<yCM<1.8"]->fill(partPt, pt_weight);
				if(p.rap() > 0)hInvariantPTK0S["pT_K0S_pPB_0<yCM<1.8"]->fill(partPt, pt_weight);
	
				//hRpPBFullyCM["RpPB_K0S_-1.8<yCM<1.8"]->fill(partPt);
				//hRpPBLowyCM["RpPB_K0S_-1.8<yCM<0"]->fill(partPt);
				//hRpPBHighyCM["RpPB_K0S_0<yCM<1.8"]->fill(partPt);
				
				if(p.rap() < -1.3001)hInvariantPTK0SpPB["pT_K0S_pPB_-1.8<yCM<-1.3"]->fill(partPt, pt_weight);
				if(p.rap() < -0.8001 & p.rap() > -1.2999)hInvariantPTK0SpPB["pT_K0S_pPB_-1.3<yCM<-0.8"]->fill(partPt, pt_weight);
				if(p.rap() < -0.3001 & p.rap() > -0.7999)hInvariantPTK0SpPB["pT_K0S_pPB_-0.8<yCM<-0.3"]->fill(partPt, pt_weight);
				if(p.rap() > 0.2999 & p.rap() < 0.8001)hInvariantPTK0SpPB["pT_K0S_pPB_0.3<yCM<0.8"]->fill(partPt, pt_weight);
				if(p.rap() > 0.7999 & p.rap() < 1.3001)hInvariantPTK0SpPB["pT_K0S_pPB_0.8<yCM<1.3"]->fill(partPt, pt_weight);
				if(p.rap() > 1.2999)hInvariantPTK0SpPB["pT_K0S_pPB_1.3<yCM<1.8"]->fill(partPt, pt_weight);
	
				hInvariantPTK0S["pT_K0S_pPB_full"]->fill(partPt, pt_weight);
				if(p.rap() < 0)hInvariantPTK0S["pT_K0S_pPB_low"]->fill(partPt, pt_weight);
				if(p.rap() > 0)hInvariantPTK0S["pT_K0S_pPB_high"]->fill(partPt, pt_weight);
				
				if(p.rap() < -0.3001 & p.rap() > -0.7999)hInvariantPTK0SpPB["K0S_pT_pPB_LowNeg"]->fill(partPt, pt_weight);
				if(p.rap() > 0.2999 & p.rap() < 0.8001)hInvariantPTK0SpPB["K0S_pT_pPB_LowPos"]->fill(partPt, pt_weight);
				if(p.rap() < -0.8001 & p.rap() > -1.2999)hInvariantPTK0SpPB["K0S_pT_pPB_MidNeg"]->fill(partPt, pt_weight);
				if(p.rap() > 0.7999 & p.rap() < 1.3001)hInvariantPTK0SpPB["K0S_pT_pPB_MidPos"]->fill(partPt, pt_weight);
				if(p.rap() < -1.3001)hInvariantPTK0SpPB["K0S_pT_pPB_HighNeg"]->fill(partPt, pt_weight);
				if(p.rap() > 1.2999)hInvariantPTK0SpPB["K0S_pT_pPB_HighPos"]->fill(partPt, pt_weight);
				
				break;
			}
			case 3122: // Lambda
			{
				hInvariantPTLambda["pT_Lambda_pPB_-1.8<yCM<0"]->fill(partPt, pt_weight);
				hInvariantPTLambda["pT_Lambda_pPB_-1.8<yCM<1.8"]->fill(partPt, pt_weight);
				hInvariantPTLambda["pT_Lambda_pPB_0<yCM<1.8"]->fill(partPt, pt_weight);
	
				//hRpPBFullyCM["RpPB_Lambda_-1.8<yCM<1.8"]->fill(partPt);
				//hRpPBLowyCM["RpPB_Lambda_-1.8<yCM<0"]->fill(partPt);
				//hRpPBHighyCM["RpPB_Lambda_0<yCM<1.8"]->fill(partPt);
				
				if(p.rap() < -1.3001)hInvariantPTLambdapPB["pT_Lambda_pPB_-1.8<yCM<-1.3"]->fill(partPt, pt_weight);
				if(p.rap() < -0.8001 & p.rap() > -1.2999)hInvariantPTLambdapPB["pT_Lambda_pPB_-1.3<yCM<-0.8"]->fill(partPt, pt_weight);
				if(p.rap() < -0.3001 & p.rap() > -0.7999)hInvariantPTLambdapPB["pT_Lambda_pPB_-0.8<yCM<-0.3"]->fill(partPt, pt_weight);
				if(p.rap() > 0.2999 & p.rap() < 0.8001)hInvariantPTLambdapPB["pT_Lambda_pPB_0.3<yCM<0.8"]->fill(partPt, pt_weight);
				if(p.rap() > 0.7999 & p.rap() < 1.3001)hInvariantPTLambdapPB["pT_Lambda_pPB_0.8<yCM<1.3"]->fill(partPt, pt_weight);
				if(p.rap() > 1.2999)hInvariantPTLambdapPB["pT_Lambda_pPB_1.3<yCM<1.8"]->fill(partPt, pt_weight);
	
				hInvariantPTLambda["pT_Lambda_pPB_full"]->fill(partPt, pt_weight);
				if(p.rap() < 0)hInvariantPTLambda["pT_Lambda_pPB_low"]->fill(partPt, pt_weight);
				if(p.rap() > 0)hInvariantPTLambda["pT_Lambda_pPB_high"]->fill(partPt, pt_weight);
				
				if(p.rap() < -0.3001 & p.rap() > -0.7999)hInvariantPTLambdapPB["Lambda_pT_pPB_LowNeg"]->fill(partPt, pt_weight);
				if(p.rap() > 0.2999 & p.rap() < 0.8001)hInvariantPTLambdapPB["Lambda_pT_pPB_LowPos"]->fill(partPt, pt_weight);
				if(p.rap() < -0.8001 & p.rap() > -1.2999)hInvariantPTLambdapPB["Lambda_pT_pPB_MidNeg"]->fill(partPt, pt_weight);
				if(p.rap() > 0.7999 & p.rap() < 1.3001)hInvariantPTLambdapPB["Lambda_pT_pPB_MidPos"]->fill(partPt, pt_weight);
				if(p.rap() < -1.3001)hInvariantPTLambdapPB["Lambda_pT_pPB_HighNeg"]->fill(partPt, pt_weight);
				if(p.rap() > 1.2999)hInvariantPTLambdapPB["Lambda_pT_pPB_HighPos"]->fill(partPt, pt_weight);
				
				break;
			}
		}
		}
		
		for (Particle p : chargedParticles)
		{
			double partPt = p.pT() / GeV;
			double pt_weight = 1. / (partPt * 2. * M_PI);
			switch (p.pid()) {
			case 3312: // Xi-
			{
				if(p.rap() < 0)hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<0"]->fill(partPt, pt_weight);
				hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<1.8"]->fill(partPt, pt_weight);
				if(p.rap() > 0)hInvariantPTXi["pT_Xi_pPB_0<yCM<1.8"]->fill(partPt, pt_weight);
				
				//hRpPBFullyCM["RpPB_Xi_-1.8<yCM<1.8"]->fill(partPt);
				//hRpPBLowyCM["RpPB_Xi_-1.8<yCM<0"]->fill(partPt);
				//hRpPBHighyCM["RpPB_Xi_0<yCM<1.8"]->fill(partPt);
				
				hInvariantPTXi["pT_Xi_pPB_full"]->fill(partPt, pt_weight);
				if(p.rap() < 0)hInvariantPTXi["pT_Xi_pPB_low"]->fill(partPt, pt_weight);
				if(p.rap() > 0)hInvariantPTXi["pT_Xi_pPB_high"]->fill(partPt, pt_weight);
				
				if(p.rap() < -0.3001 & p.rap() > -0.7999)h["negative_charged_yCM_low"]->fill(partPt);
				if(p.rap() > 0.2999 & p.rap() < 0.8001)h["positive_charged_yCM_low"]->fill(partPt);
				if(p.rap() < -0.8001 & p.rap() > -1.2999)h["negative_charged_yCM_mid"]->fill(partPt);
				if(p.rap() > 0.7999 & p.rap() < 1.3001)h["positive_charged_yCM_mid"]->fill(partPt);
				if(p.rap() < -1.3001)h["negative_charged_yCM_high"]->fill(partPt);
				if(p.rap() > 1.2999)h["positive_charged_yCM_high"]->fill(partPt);
				
				break;
			}
			case -3312: // Xi+
			{
				if(p.rap() < 0)hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<0"]->fill(partPt, pt_weight);
				hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<1.8"]->fill(partPt, pt_weight);
				if(p.rap() > 0)hInvariantPTXi["pT_Xi_pPB_0<yCM<1.8"]->fill(partPt, pt_weight);
				
				//hRpPBFullyCM["RpPB_Xi_-1.8<yCM<1.8"]->fill(partPt);
				//hRpPBLowyCM["RpPB_Xi_-1.8<yCM<0"]->fill(partPt);
				//hRpPBHighyCM["RpPB_Xi_0<yCM<1.8"]->fill(partPt);
				
				hInvariantPTXi["pT_Xi_pPB_full"]->fill(partPt, pt_weight);
				if(p.rap() < 0)hInvariantPTXi["pT_Xi_pPB_low"]->fill(partPt, pt_weight);
				if(p.rap() > 0)hInvariantPTXi["pT_Xi_pPB_high"]->fill(partPt, pt_weight);
				
				if(p.rap() < -0.3001 & p.rap() > -0.7999)h["negative_charged_yCM_low"]->fill(partPt);
				if(p.rap() > 0.2999 & p.rap() < 0.8001)h["positive_charged_yCM_low"]->fill(partPt);
				if(p.rap() < -0.8001 & p.rap() > -1.2999)h["negative_charged_yCM_mid"]->fill(partPt);
				if(p.rap() > 0.7999 & p.rap() < 1.3001)h["positive_charged_yCM_mid"]->fill(partPt);
				if(p.rap() < -1.3001)h["negative_charged_yCM_high"]->fill(partPt);
				if(p.rap() > 1.2999)h["positive_charged_yCM_high"]->fill(partPt);
				
				break;
			}
			case 3334: // Omega-
			{
				
				hInvariantPTOmega["pT_Omega_pPB_-1.8<yCM<1.8"]->fill(partPt, pt_weight);

				//hRpPBFullyCM["RpPB_Omega_-1.8<yCM<1.8"]->fill(partPt);
				
				hInvariantPTOmega["pT_Omega_pPB_full"]->fill(partPt, pt_weight);
				
				if(p.rap() < -0.3001 & p.rap() > -0.7999)h["negative_charged_yCM_low"]->fill(partPt);
				if(p.rap() > 0.2999 & p.rap() < 0.8001)h["positive_charged_yCM_low"]->fill(partPt);
				if(p.rap() < -0.8001 & p.rap() > -1.2999)h["negative_charged_yCM_mid"]->fill(partPt);
				if(p.rap() > 0.7999 & p.rap() < 1.3001)h["positive_charged_yCM_mid"]->fill(partPt);
				if(p.rap() < -1.3001)h["negative_charged_yCM_high"]->fill(partPt);
				if(p.rap() > 1.2999)h["positive_charged_yCM_high"]->fill(partPt);
				
				break;
			}
			case -3334: // Omega+
			{
				
				hInvariantPTOmega["pT_Omega_pPB_-1.8<yCM<1.8"]->fill(partPt, pt_weight);

				//hRpPBFullyCM["RpPB_Omega_-1.8<yCM<1.8"]->fill(partPt);
				
				hInvariantPTOmega["pT_Omega_pPB_full"]->fill(partPt, pt_weight);
				
				if(p.rap() < -0.3001 & p.rap() > -0.7999)h["negative_charged_yCM_low"]->fill(partPt);
				if(p.rap() > 0.2999 & p.rap() < 0.8001)h["positive_charged_yCM_low"]->fill(partPt);
				if(p.rap() < -0.8001 & p.rap() > -1.2999)h["negative_charged_yCM_mid"]->fill(partPt);
				if(p.rap() > 0.7999 & p.rap() < 1.3001)h["positive_charged_yCM_mid"]->fill(partPt);
				if(p.rap() < -1.3001)h["negative_charged_yCM_high"]->fill(partPt);
				if(p.rap() > 1.2999)h["positive_charged_yCM_high"]->fill(partPt);
				
				break;
			}
		}
	}
	return;
	}
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    
//    bool pp_available = false;
//    bool pPB_available = false;
    
//    for (auto element : hInvariantPTK0S)
//	{
//		string name = element.second->name();
//		if (name.find("pp") != std::string::npos)
//		{
//			if (element.second->numEntries() > 0) pp_available = true;
//			else
//			{
//				pp_available = false;
//				break;
//			}
//		}
//		else if (name.find("pPB") != std::string::npos)
//		{
//			if (element.second->numEntries() > 0) pPB_available = true;
//			else
//			{
//				pPB_available = false;
//				break;
//			}
//		}
//	}
	
//	if (!(pp_available && pPB_available)) return;

	//*******Scale histograms*********
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pp_-1.8<yCM<0"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pp_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pp_0<yCM<1.8"]->scaleW(1. / sow["sow_pp"]->sumW());
		
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pp_full"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pp_low"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pp_high"]->scaleW(1. / sow["sow_pp"]->sumW());
				
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pp_-1.8<yCM<0"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pp_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pp_0<yCM<1.8"]->scaleW(1. / sow["sow_pp"]->sumW());	
		
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pp_full"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pp_low"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pp_high"]->scaleW(1. / sow["sow_pp"]->sumW());

	if(sow["sow_pp"]->sumW() > 0) hInvariantPTXi["pT_Xi_pp_-1.8<yCM<0"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTXi["pT_Xi_pp_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTXi["pT_Xi_pp_0<yCM<1.8"]->scaleW(1. / sow["sow_pp"]->sumW());
					
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTXi["pT_Xi_pp_full"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTXi["pT_Xi_pp_low"]->scaleW(1. / sow["sow_pp"]->sumW());
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTXi["pT_Xi_pp_high"]->scaleW(1. / sow["sow_pp"]->sumW());

	if(sow["sow_pp"]->sumW() > 0) hInvariantPTOmega["pT_Omega_pp_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pp"]->sumW());
		
	if(sow["sow_pp"]->sumW() > 0) hInvariantPTOmega["pT_Omega_pp_full"]->scaleW(1. / sow["sow_pp"]->sumW());
				
		
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pPB_-1.8<yCM<0"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pPB_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pPB_0<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
	
//	if(sow["sow_pPB"]->sumW() > 0) hRpPBFullyCM["RpPB_K0S_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
//	if(sow["sow_pPB"]->sumW() > 0) hRpPBLowyCM["RpPB_K0S_-1.8<yCM<0"]->scaleW(1. / sow["sow_pPB"]->sumW());
//	if(sow["sow_pPB"]->sumW() > 0) hRpPBHighyCM["RpPB_K0S_0<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
			
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0SpPB["pT_K0S_pPB_-1.8<yCM<-1.3"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0SpPB["pT_K0S_pPB_-1.3<yCM<-0.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0SpPB["pT_K0S_pPB_-0.8<yCM<-0.3"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0SpPB["pT_K0S_pPB_0.3<yCM<0.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0SpPB["pT_K0S_pPB_0.8<yCM<1.3"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0SpPB["pT_K0S_pPB_1.3<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
				
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pPB_full"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pPB_low"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTK0S["pT_K0S_pPB_high"]->scaleW(1. / sow["sow_pPB"]->sumW());
	
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pPB_-1.8<yCM<0"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pPB_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pPB_0<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
	
//	if(sow["sow_pPB"]->sumW() > 0) hRpPBFullyCM["RpPB_Lambda_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
//	if(sow["sow_pPB"]->sumW() > 0) hRpPBLowyCM["RpPB_Lambda_-1.8<yCM<0"]->scaleW(1. / sow["sow_pPB"]->sumW());
//	if(sow["sow_pPB"]->sumW() > 0) hRpPBHighyCM["RpPB_Lambda_0<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
				
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambdapPB["pT_Lambda_pPB_-1.8<yCM<-1.3"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambdapPB["pT_Lambda_pPB_-1.3<yCM<-0.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambdapPB["pT_Lambda_pPB_-0.8<yCM<-0.3"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambdapPB["pT_Lambda_pPB_0.3<yCM<0.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambdapPB["pT_Lambda_pPB_0.8<yCM<1.3"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambdapPB["pT_Lambda_pPB_1.3<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
	
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pPB_full"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pPB_low"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTLambda["pT_Lambda_pPB_high"]->scaleW(1. / sow["sow_pPB"]->sumW());
	
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<0"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTXi["pT_Xi_pPB_0<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
				
//	if(sow["sow_pPB"]->sumW() > 0) hRpPBFullyCM["RpPB_Xi_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
//	if(sow["sow_pPB"]->sumW() > 0) hRpPBLowyCM["RpPB_Xi_-1.8<yCM<0"]->scaleW(1. / sow["sow_pPB"]->sumW());
//	if(sow["sow_pPB"]->sumW() > 0) hRpPBHighyCM["RpPB_Xi_0<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
				
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTXi["pT_Xi_pPB_full"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTXi["pT_Xi_pPB_low"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTXi["pT_Xi_pPB_high"]->scaleW(1. / sow["sow_pPB"]->sumW());
				
	hInvariantPTOmega["pT_Omega_pPB_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());

//	hRpPBFullyCM["RpPB_Omega_-1.8<yCM<1.8"]->scaleW(1. / sow["sow_pPB"]->sumW());
				
	if(sow["sow_pPB"]->sumW() > 0) h["negative_charged_yCM_low"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) h["positive_charged_yCM_low"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) h["negative_charged_yCM_mid"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) h["positive_charged_yCM_mid"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) h["negative_charged_yCM_high"]->scaleW(1. / sow["sow_pPB"]->sumW());
	if(sow["sow_pPB"]->sumW() > 0) h["positive_charged_yCM_high"]->scaleW(1. / sow["sow_pPB"]->sumW());
	
	if(sow["sow_pPB"]->sumW() > 0) hInvariantPTOmega["pT_Omega_pPB_full"]->scaleW(1. / sow["sow_pPB"]->sumW());
	
//	//Figure 2 Y scalings
//	//hInvariantPTK0S["pT_K0S_pp_-1.8<yCM<0"]->scaleY(10.0);
//	//hInvariantPTK0S["pT_K0S_pp_0<yCM<1.8"]->scaleY(1. / 10.0);
//	//hInvariantPTLambda["pT_Lambda_pp_-1.8<yCM<0"]->scaleY(10.0);
//	//hInvariantPTLambda["pT_Lambda_pp_0<yCM<1.8"]->scaleY(1. / 10.0);
//	//hInvariantPTXi["pT_Xi_pp_-1.8<yCM<0"]->scaleY(10.0);
//	//hInvariantPTXi["pT_Xi_pp_0<yCM<1.8"]->scaleY(1. / 10.0);
//	//hInvariantPTOmega["pT_Omega_pp_-1.8<yCM<0"]->scaleY(10.0);
//	//hInvariantPTOmega["pT_Omega_pp_0<yCM<1.8"]->scaleY(1. / 10.0);
	
	//scale by # of binary collisions, N_coll
//	hInvariantPTK0S["pT_K0S_pPB_-1.8<yCM<0"]->scaleW(1. / 6.9);
//	hInvariantPTK0S["pT_K0S_pPB_-1.8<yCM<1.8"]->scaleW(1. / 6.9);
//	hInvariantPTK0S["pT_K0S_pPB_0<yCM<1.8"]->scaleW(1. / 6.9);
//	hInvariantPTLambda["pT_Lambda_pPB_-1.8<yCM<0"]->scaleW(1. / 6.9);
//	hInvariantPTLambda["pT_Lambda_pPB_-1.8<yCM<1.8"]->scaleW(1. / 6.9);
//	hInvariantPTLambda["pT_Lambda_pPB_0<yCM<1.8"]->scaleW(1. / 6.9);
//	hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<0"]->scaleW(1. / 6.9);
//	hInvariantPTXi["pT_Xi_pPB_-1.8<yCM<1.8"]->scaleW(1. / 6.9);
//	hInvariantPTXi["pT_Xi_pPB_0<yCM<1.8"]->scaleW(1. / 6.9);
//	hInvariantPTOmega["pT_Omega_pPB_-1.8<yCM<1.8"]->scaleW(1. / 6.9);
	
	//Figure 5 Y scalings
//	//hInvariantPTK0SpPB["pT_K0S_pPB_-1.8<yCM<-1.3"]->scaleY(100.0);
//	//hInvariantPTK0SpPB["pT_K0S_pPB_-1.3<yCM<-0.8"]->scaleY(10.0);
//	//hInvariantPTK0SpPB["pT_K0S_pPB_0.3<yCM<0.8"]->scaleY(1. / 10.0);
//	//hInvariantPTK0SpPB["pT_K0S_pPB_0.8<yCM<1.3"]->scaleY(1. / 100.0);
//	//hInvariantPTK0SpPB["pT_K0S_pPB_1.3<yCM<1.8"]->scaleY(1. / 1000.0);
//	//hInvariantPTK0SpPB["pT_Lambda_pPB_-1.8<yCM<-1.3"]->scaleY(100.0);
//	//hInvariantPTK0SpPB["pT_Lambda_pPB_-1.3<yCM<-0.8"]->scaleY(10.0);
//	//hInvariantPTK0SpPB["pT_Lambda_pPB_0.3<yCM<0.8"]->scaleY(1. / 10.0);
//	//hInvariantPTK0SpPB["pT_Lambda_pPB_0.8<yCM<1.3"]->scaleY(1. / 100.0);
//	//hInvariantPTK0SpPB["pT_Lambda_pPB_1.3<yCM<1.8"]->scaleY(1. / 1000.0);
	
		return;
	
	//Scatter plots from division
//	divide(hInvariantPTXi["pT_Xi_pPB_high"], hInvariantPTXi["pT_Xi_pp_high"], RpPBHighyCM["pPB_pT_Xi/pp_pT_Xi_high"]);
//	divide(hInvariantPTXi["pT_Xi_pPB_low"], hInvariantPTXi["pT_Xi_pp_low"], RpPBLowyCM["pPB_pT_Xi/pp_pT_Xi_low"]);
//	divide(hInvariantPTXi["pT_Xi_pPB_full"], hInvariantPTXi["pT_Xi_pp_full"], RpPBHighyCM["pPB_pT_Xi/pp_pT_Xi_full"]);
	
//	divide(hInvariantPTLambda["pT_Lambda_pPB_high"], hInvariantPTLambda["pT_Lambda_pp_high"], RpPBHighyCM["pPB_pT_Lambda/pp_pT_Lambda_high"]);
//	divide(hInvariantPTLambda["pT_Lambda_pPB_low"], hInvariantPTLambda["pT_Lambda_pp_low"], RpPBHighyCM["pPB_pT_Lambda/pp_pT_Lambda_low"]);
//	divide(hInvariantPTLambda["pT_Lambda_pPB_full"], hInvariantPTLambda["pT_Lambda_pp_full"], RpPBHighyCM["pPB_pT_Lambda/pp_pT_Lambda_full"]);
	
//	divide(hInvariantPTK0S["pT_K0S_pPB_high"], hInvariantPTK0S["pT_K0S_pp_high"], RpPBHighyCM["pPB_pT_K0S/pp_pT_K0S_high"]);
//	divide(hInvariantPTK0S["pT_K0S_pPB_low"], hInvariantPTK0S["pT_K0S_pp_low"], RpPBHighyCM["pPB_pT_K0S/pp_pT_K0S_low"]);
//	divide(hInvariantPTK0S["pT_K0S_pPB_full"], hInvariantPTK0S["pT_K0S_pp_full"], RpPBHighyCM["pPB_pT_K0S/pp_pT_K0S_full"]);
	
//	divide(hInvariantPTOmega["pT_Omega_pPB_full"], hInvariantPTOmega["pT_Omega_pp_full"], RpPBHighyCM["pPB_pT_Omega/pp_pT_Omega_full"]);
	
//	divide(hInvariantPTK0SpPB["K0S_pT_pPB_LowNeg"], hInvariantPTK0SpPB["K0S_pT_pPB_LowPos"], YasymLow["K0S_-0.8<yCM<-0.3/0.3<yCM<0.8"]);
//	divide(hInvariantPTLambdapPB["Lambda_pT_pPB_LowNeg"], hInvariantPTLambdapPB["Lambda_pT_pPB_LowPos"], YasymLow["Lambda_-0.8<yCM<-0.3/0.3<yCM<0.8"]);
//	divide(h["negative_charged_yCM_low"], h["positive_charged_yCM_low"], YasymLow["h+/-_-0.8<yCM<-0.3/0.3<yCM<0.8"]);

//	divide(hInvariantPTK0SpPB["K0S_pT_pPB_MidNeg"], hInvariantPTK0SpPB["K0S_pT_pPB_MidPos"], YasymMid["K0S_-1.3<yCM<-0.8/0.8<yCM<1.3"]);
//	divide(hInvariantPTLambdapPB["Lambda_pT_pPB_MidNeg"], hInvariantPTLambdapPB["Lambda_pT_pPB_MidPos"], YasymMid["Lambda_-1.3<yCM<-0.8/0.8<yCM<1.3"]);
//	divide(h["negative_charged_yCM_mid"], h["positive_charged_yCM_mid"], YasymMid["h+/-_-1.3<yCM<-0.8/0.8<yCM<1.3"]);
	
//	divide(hInvariantPTK0SpPB["K0S_pT_pPB_HighNeg"], hInvariantPTK0SpPB["K0S_pT_pPB_HighPos"], YasymHigh["K0S_-1.8<yCM<-1.3/1.3<yCM<1.8"]);
//	divide(hInvariantPTLambdapPB["Lambda_pT_pPB_HighNeg"], hInvariantPTLambdapPB["Lambda_pT_pPB_HighPos"], YasymHigh["Lambda_-1.8<yCM<-1.3/1.3<yCM<1.8"]);
//	divide(h["negative_charged_yCM_high"], h["positive_charged_yCM_high"], YasymHigh["h+/-_-1.8<yCM<-1.3/1.3<yCM<1.8"]);
	
    }

    	map<string, Histo1DPtr> hInvariantPTK0S;
    	map<string, Histo1DPtr> hInvariantPTLambda;
    	map<string, Histo1DPtr> hInvariantPTXi;
    	map<string, Histo1DPtr> hInvariantPTOmega;
    	map<string, Histo1DPtr> hRpPBFullyCM;
    	map<string, Histo1DPtr> hRpPBLowyCM;
    	map<string, Histo1DPtr> hRpPBHighyCM;
    	map<string, Histo1DPtr> hInvariantPTK0SpPB;
    	map<string, Histo1DPtr> hInvariantPTLambdapPB;
    	map<string, Histo1DPtr> h;
    	
	map<string, Scatter2DPtr> RpPBFullyCM;
	map<string, Scatter2DPtr> RpPBLowyCM;
	map<string, Scatter2DPtr> RpPBHighyCM;
	map<string, Scatter2DPtr> YasymLow;
	map<string, Scatter2DPtr> YasymMid;
	map<string, Scatter2DPtr> YasymHigh;
	
	map<string, CounterPtr> sow;
	string beamOpt;
	enum CollisionSystem {pp, pPB};
	CollisionSystem collSys;
  };


  DECLARE_RIVET_PLUGIN(CMS_2020_I064906);

}
