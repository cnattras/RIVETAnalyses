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
#include "RHICCentrality.hh" //external header for Centrality calculation
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES


namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2006_I0611006 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2006_I0611006);


    /// Book histograms and initialise projections before the run
    void init() {
    	//Particles: eta, pi^+, pi^-, pi^0, gamma
    	std::initializer_list<int> pdgIds = { 221, 211, -211, 111, 22};
    	const PrimaryParticles fs(pdgIds, Cuts::abseta < 0.35 && Cuts::abscharge > 0);
    	declare(fs, "fs");
    	
    	const PrimaryParticles ns(pdgIds, Cuts::abseta < 0.35 && Cuts::abscharge == 0);
    	declare(ns, "ns");
    	
    	beamOpt = getOption<string>("beam", "NONE");

	if (beamOpt == "PP") collSys = pp;
	else if (beamOpt == "AUAU200") collSys = AuAu200;
	else if (beamOpt == "dAU200") collSys = dAu200;
	
	if (!(collSys == pp)) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

	//Create various counters
	book(sow["sow_pp"], "sow_pp");
	book(sow["sow_AUAUc20"], "sow_AUAUc20");
	book(sow["sow_AUAUc40"], "sow_AUAUc40");
	book(sow["sow_AUAUc60"], "sow_AUAUc60");
	book(sow["sow_AUAUc92"], "sow_AUAUc92");
	book(sow["sow_dAUc20"], "sow_dAUc20");
	book(sow["sow_dAUc40"], "sow_dAUc40");
	book(sow["sow_dAUc60"], "sow_dAUc60");
	book(sow["sow_dAUc88"], "sow_dAUc88");
	
    	//--------Begin booking histograms-----------
    	
    	//Cross sections in pp (fig 12.1)
	book(hCrossSec["xSection_pp_GammaGamma"], 1, 1, 1);
	
	//Cross sections in pp (fig 12.2)
	book(hCrossSec["xSection_pp_Pions"], 2, 1, 1);
	
	//Cross sections in dAu (fig 13.1)
	book(hCrossSec["xSection_dAu_GammaGamma"], 3, 1, 1);
	
	//Cross sections in dAu (fig 13.2)
	book(hCrossSec["xSection_dAu_Pions"], 4, 1, 1);
	
	//Ratio eta/pi^0 in pp (fig 18)
	book(hRatiopp["ratioEtaPi0_pp"], 9, 1, 1);
	
	//R_dAu 0-20 (fig 16)
	book(hRdAu20["RdAu_c20"], 7, 1, 2);
	
	//R_dAu 20-40 (fig 16)
	book(hRdAu40["RdAu_c40"], 7, 1, 3);
	
	//R_dAu 40-60 (fig 16)
	book(hRdAu60["RdAu_c60"], 7, 1, 4);
	
	//R_dAu 60-88 (fig 16)
	book(hRdAu88["RdAu_cMB"], 7, 1, 5);
	
	//R_dAu 0-88(MB) (fig 16)
	book(hRdAuMB["RdAu_cMB"], 7, 1, 1);
	
	//R_AuAu 0-20 (fig 17)
	book(hRAuAu["RAuAu_c20"], 8, 1, 1);
	
	//R_AuAu 20-60 (fig 17)
	book(hRAuAu["RAuAu_c60"], 8, 1, 2);
	
	//R_AuAu 60-92 (fig 17)
	book(hRAuAu["RAuAu_c92"], 8, 1, 3);
	
	//Ratio eta/pi^0 in Au+Au 0-92 MB (fig 20) ***UNUSED IN HISTO***
	//string refname31 = mkAxisCode(11, 1, 1);
	//const Scatter2d& refdata31 = refData(refname31);
	//book(hRatioAuAu["ratioEtaPi0_cMB"], refname31 + "_AuAu", refdata31);
	
	//Ratio eta/pi^0 in Au+Au 0-20 (fig 20)
	book(hRatioAuAu["ratioEtaPi0_c20"], 11, 1, 2);
	
	//Ratio eta/pi^0 in Au+Au 20-60 (fig 20)
	book(hRatioAuAu["ratioEtaPi0_c60"], 11, 1, 3);
	
	//Ratio eta/pi^0 in Au+Au 60-92 (fig 20)
	book(hRatioAuAu["ratioEtaPi0_c92"], 11, 1, 4);
	
	//Ratio eta/pi^0 in Au+Au d+Au MB comparison (fig 20) 
	book(hRatioAuAu["ratioEtaPi0_cMB_dAu"], 10, 1, 1);
	
	//invariant yields in d+Au 0-20 (fig 14)
	book(hdAuYields["pTyields_c20"], 5, 1, 1);
	
	//invariant yields in d+Au 20-40 (fig 14)
	book(hdAuYields["pTyields_c40"], 5, 1, 2);
	
	//invariant yields in d+Au 40-60 (fig 14)
	book(hdAuYields["pTyields_c60"], 5, 1, 3);
	
	//invariant yields in d+Au 60-88 (fig 14)
	book(hdAuYields["pTyields_c88"], 5, 1, 4);
	
	//invariant yields in Au+Au 0-92 (fig 15)
	book(hAuAuYields["pTyields_cMB"], 6, 1, 1);
	
	//invariant yields in Au+Au 0-20 (fig 15)
	book(hAuAuYields["pTyields_c20"], 6, 1, 2);
	
	//invariant yields in Au+Au 20-40 (fig 15)
	book(hAuAuYields["pTyields_c40"], 6, 1, 3);
	
	//invariant yields in Au+Au 60-92 (fig 15)
	book(hAuAuYields["pTyields_c92"], 6, 1, 4);
	
	//Ratio eta/pi^0 in d+Au 0-88 (fig 19)  ***UNUSED IN HISTO***
	//string refname26 = mkAxisCode(10, 1, 1);
	//const Scatter2d& refdata26 = refData(refname26);
	//book(hRatiodAu["ratioEtaPi0_cMB"], refname26 + "_dAu", refdata26);
	
	//Ratio eta/pi^0 in d+Au 0-20 (fig 19)
	book(hRatiodAu["ratioEtaPi0_c20"], 10, 1, 2);
	
	//Ratio eta/pi^0 in d+Au 19-40 (fig 19)
	book(hRatiodAu["ratioEtaPi0_c40"], 10, 1, 3);
	
	//Ratio eta/pi^0 in d+Au 40-60 (fig 19)
	book(hRatiodAu["ratioEtaPi0_c60"], 10, 1, 4);
	
	//Ratio eta/pi^0 in d+Au 60-88 (fig 19)
	book(hRatiodAu["ratioEtaPi0_c88"], 10, 1, 5);
	
	//-----Book Scatter Plots from division-------
	
	string refname1 = mkAxisCode(11, 1, 2);
	const Scatter2D& refdata1 = refData(refname1);
	book(hRAuAu["RAuAu_c20"], refname1 + "_RAuAu", refdata1);
	book(hCrossSec["xSection_pp_GammaGamma"], refname1 + "_pp_GammaGamma", refdata1);
	book(RatioAuAu["AuAuc/pp_c20"], refname1);
	
	string refname2 = mkAxisCode(11, 1, 3);
	const Scatter2D& refdata2 = refData(refname2);
	book(hRAuAu["RAuAu_c60"], refname2 + "_RAuAu", refdata2);
	book(hCrossSec["xSection_pp_GammaGamma"], refname2 + "_pp_GammaGamma", refdata2);
	book(RatioAuAu["AuAuc/pp_c60"], refname2);
	
	string refname3 = mkAxisCode(11, 1, 4);
	const Scatter2D& refdata3 = refData(refname3);
	book(hRAuAu["RAuAu_c92"], refname3 + "_RAuAu", refdata3);
	book(hCrossSec["xSection_pp_GammaGamma"], refname3 + "_pp_GammaGamma", refdata3);
	book(RatioAuAu["AuAuc/pp_c92"], refname3);
	
	string refname4 = mkAxisCode(10, 1, 2);
	const Scatter2D& refdata4 = refData(refname4);
	book(hRdAu20["RdAu_c20"], refname4 + "_RdAu", refdata4);
	book(hCrossSec["xSection_pp_GammaGamma"], refname4 + "_pp_GammaGamma", refdata4);
	book(RatiodAu20["dAuc/pp_c20"], refname4);
	
	string refname5 = mkAxisCode(10, 1, 3);
	const Scatter2D& refdata5 = refData(refname5);
	book(hRdAu40["RdAu_c40"], refname5 + "_RdAu", refdata5);
	book(hCrossSec["xSection_pp_GammaGamma"], refname5 + "_pp_GammaGamma", refdata5);
	book(RatiodAu40["dAuc/pp_c40"], refname5);
	
	string refname6 = mkAxisCode(10, 1, 4);
	const Scatter2D& refdata6 = refData(refname6);
	book(hRdAu60["RdAu_c60"], refname6 + "_RdAu", refdata6);
	book(hCrossSec["xSection_pp_GammaGamma"], refname6 + "_pp_GammaGamma", refdata6);
	book(RatiodAu60["dAuc/pp_c60"], refname6);
	
	string refname7 = mkAxisCode(10, 1, 5);
	const Scatter2D& refdata7 = refData(refname7);
	book(hRdAu88["RdAu_c88"], refname7 + "_RdAu", refdata7);
	book(hCrossSec["xSection_pp_GammaGamma"], refname7 + "_pp_GammaGamma", refdata7);
	book(RatiodAu88["dAuc/pp_c88"], refname7);
	
	string refname8 = mkAxisCode(10, 1, 1);
	const Scatter2D& refdata8 = refData(refname8);
	book(hRdAuMB["RdAu_cMB"], refname8 + "_RdAu", refdata8);
	book(hCrossSec["xSection_pp_GammaGamma"], refname8 + "_pp_GammaGamma", refdata8);
	book(RatiodAuMB["dAuc/pp_cMB"], refname8);
	
    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {
    
      Particles chargedParticles = applyProjection<PrimaryParticles>(event, "fs").particles();
      Particles neutralParticles = applyProjection<PrimaryParticles>(event, "ns").particles();
    
      if (collSys == pp)
	{
		sow["sow_pp"]->fill();
		for (Particle p : chargedParticles)
		{
			double partPt = p.pT() / GeV;
			switch (p.pid()) {
			case 211: // pi+
			{
				//hCrossSec["xSection_pp_GammaGamma"]->fill(partPt);
				hCrossSec["xSection_pp_Pions"]->fill(partPt);
				//hRatiopp["ratioEtaPi0_pp"]->fill(partPt);
				
				break;
			}
			case -211: // pi-
			{
				//hCrossSec["xSection_pp_GammaGamma"]->fill(partPt);
				hCrossSec["xSection_pp_Pions"]->fill(partPt);
				//hRatiopp["ratioEtaPi0_pp"]->fill(partPt);
				
				break;
			}
		}
		}
		
		for (Particle p : neutralParticles)
		{
			double partPt = p.pT() / GeV;
			switch (p.pid()) {
			case 221: // eta
			{
				hCrossSec["xSection_pp_GammaGamma"]->fill(partPt);
				hCrossSec["xSection_pp_Pions"]->fill(partPt);
				hRatiopp["ratioEtaPi0_pp"]->fill(partPt);
				
				break;
			}
			case 111: // pi^0
			{
				//hCrossSec["xSection_pp_GammaGamma"]->fill(partPt);
				hCrossSec["xSection_pp_Pions"]->fill(partPt);
				hRatiopp["ratioEtaPi0_pp"]->fill(partPt);
				
				break;
			}
			case 22: // gamma
			{
				hCrossSec["xSection_pp_GammaGamma"]->fill(partPt);
				//hCrossSec["xSection_pp_Pions"]->fill(partPt);
				//hRatiopp["ratioEtaPi0_pp"]->fill(partPt);
				
				break;
			}
		}
	}
	return;
	}
      
        const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
	const double c = cent();


	if (collSys == AuAu200)
	{
		if ((c < 0.) || (c > 92.)) vetoEvent;

		if ((c >= 0.) && (c < 20.))
		{
			sow["sow_AUAUc20"]->fill();
			for (const Particle& p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 211: // pi+
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case -211: // pi-
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
			for (const Particle& p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 221: // eta
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 111: // pi^0
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 22: // gamma
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
		}
		
		else if ((c >= 20.) && (c < 40.))
		{
			sow["sow_AUAUc40"]->fill();

			for (const Particle& p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 211: // pi+
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case -211: // pi-
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
			for (const Particle& p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 221: // eta
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 111: // pi^0
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 22: // gamma
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
		}
		
		else if ((c >= 40.) && (c < 60.))
		{
			sow["sow_AUAUc60"]->fill();
			for (const Particle& p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 211: // pi+
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case -211: // pi-
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
			for (const Particle& p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 221: // eta
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 111: // pi^0
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 22: // gamma
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
		}
		
		else if ((c >= 60.) && (c < 92.))
		{
			sow["sow_AUAUc92"]->fill();
			for (const Particle& p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 211: // pi+
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case -211: // pi-
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
			for (const Particle& p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 221: // eta
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 111: // pi^0
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 22: // gamma
				{
					hAuAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB"]->fill(partPt, pt_weight);
					
					hRAuAu["RAuAu_c20"]->fill(partPt);
					hRAuAu["RAuAu_c60"]->fill(partPt);
					hRAuAu["RAuAu_c92"]->fill(partPt);
					
					//hRatioAuAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_c92"]->fill(partPt);
					//////hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
		}
		return;
	}
	
	
	if (collSys == dAu200)
	{
		if ((c < 0.) || (c > 88.)) vetoEvent;
                   
		if ((c >= 0.) && (c < 20.))
		{
			sow["sow_dAUc20"]->fill();
			for (const Particle& p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 211: // pi+
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case -211: // pi-
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
			for (const Particle& p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 221: // eta
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 111: // pi^0
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 22: // gamma
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					//hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
		}
		
		else if ((c >= 20.) && (c < 40.))
		{
			sow["sow_dAUc40"]->fill();
			for (const Particle& p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 211: // pi+
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case -211: // pi-
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
			for (const Particle& p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 221: // eta
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 111: // pi^0
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 22: // gamma
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					//hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
		}
		
		else if ((c >= 40.) && (c < 60.))
		{
			sow["sow_dAUc60"]->fill();
			for (const Particle& p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 211: // pi+
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case -211: // pi-
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
			for (const Particle& p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 221: // eta
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 111: // pi^0
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 22: // gamma
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					//hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
		}

		else if ((c >= 60.) && (c < 88.))
		{
			sow["sow_dAUc88"]->fill();
			for (const Particle& p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 211: // pi+
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case -211: // pi-
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
			
			for (const Particle& p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {
				case 221: // eta
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 111: // pi^0
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					//hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				case 22: // gamma
				{
					hdAuYields["pTyields_c20"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c40"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c60"]->fill(partPt, pt_weight);
					hdAuYields["pTyields_c88"]->fill(partPt, pt_weight);
					
					hCrossSec["xSection_dAu_GammaGamma"]->fill(partPt);
					//hCrossSec["xSection_dAu_Pions"]->fill(partPt);
					
					hRdAu20["RdAu_c20"]->fill(partPt);
					hRdAu40["RdAu_c40"]->fill(partPt);
					hRdAu60["RdAu_c60"]->fill(partPt);
					hRdAu88["RdAu_c88"]->fill(partPt);
					hRdAuMB["RdAu_cMB"]->fill(partPt);
					
					//hRatiodAu["ratioEtaPi0_c20"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c40"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c60"]->fill(partPt);
					//hRatiodAu["ratioEtaPi0_c88"]->fill(partPt);
					//hRatioAuAu["ratioEtaPi0_cMB_dAu"]->fill(partPt);
					
					break;
				}
				}
			}
		}
		return;
	}

    }


    /// Normalise histograms etc., after the run
    void finalize() {

	bool AuAu200_available = false;
	bool dAu200_available = false;
	bool pp_available = false;
			
	for (auto element : hdAuYields)
	{
		string name = element.second->name();
		if (name.find("AuAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) AuAu200_available = true;
			else
			{
				AuAu200_available = false;
				break;
			}
		}
		else if (name.find("dAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) dAu200_available = true;
			else
			{
				dAu200_available = false;
				break;
			}
		}
		else if (name.find("pp") != std::string::npos)
		{
			if (element.second->numEntries() > 0) pp_available = true;
			else
			{
				pp_available = false;
				break;
			}
		}
	}
	
	
	for (auto element : hAuAuYields)
	{
		string name = element.second->name();
		if (name.find("AuAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) AuAu200_available = true;
			else
			{
				AuAu200_available = false;
				break;
			}
		}
		else if (name.find("dAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) dAu200_available = true;
			else
			{
				dAu200_available = false;
				break;
			}
		}
		else if (name.find("pp") != std::string::npos)
		{
			if (element.second->numEntries() > 0) pp_available = true;
			else
			{
				pp_available = false;
				break;
			}
		}
	}
	
	for (auto element : hCrossSec)
	{
		string name = element.second->name();
		if (name.find("AuAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) AuAu200_available = true;
			else
			{
				AuAu200_available = false;
				break;
			}
		}
		else if (name.find("dAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) dAu200_available = true;
			else
			{
				dAu200_available = false;
				break;
			}
		}
		else if (name.find("pp") != std::string::npos)
		{
			if (element.second->numEntries() > 0) pp_available = true;
			else
			{
				pp_available = false;
				break;
			}
		}
	}
	
	for (auto element : hRatiopp)
	{
		string name = element.second->name();
		if (name.find("AuAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) AuAu200_available = true;
			else
			{
				AuAu200_available = false;
				break;
			}
		}
		else if (name.find("dAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) dAu200_available = true;
			else
			{
				dAu200_available = false;
				break;
			}
		}
		else if (name.find("pp") != std::string::npos)
		{
			if (element.second->numEntries() > 0) pp_available = true;
			else
			{
				pp_available = false;
				break;
			}
		}
	}
	
	for (auto element : hRatioAuAu)
	{
		string name = element.second->name();
		if (name.find("AuAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) AuAu200_available = true;
			else
			{
				AuAu200_available = false;
				break;
			}
		}
		else if (name.find("dAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) dAu200_available = true;
			else
			{
				dAu200_available = false;
				break;
			}
		}
		else if (name.find("pp") != std::string::npos)
		{
			if (element.second->numEntries() > 0) pp_available = true;
			else
			{
				pp_available = false;
				break;
			}
		}
	}
	
	for (auto element : hRatiodAu)
	{
		string name = element.second->name();
		if (name.find("AuAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) AuAu200_available = true;
			else
			{
				AuAu200_available = false;
				break;
			}
		}
		else if (name.find("dAu") != std::string::npos)
		{
			if (element.second->numEntries() > 0) dAu200_available = true;
			else
			{
				dAu200_available = false;
				break;
			}
		}
		else if (name.find("pp") != std::string::npos)
		{
			if (element.second->numEntries() > 0) pp_available = true;
			else
			{
				pp_available = false;
				break;
			}
		}
	}
	
	if (!(AuAu200_available && dAu200_available && pp_available)) return;
	
	cout << "Here ************" << endl;

	//Scale histograms
	
	hCrossSec["xSection_pp_GammaGamma"]->scaleW(1. / sow["sow_pp"]->sumW());
	
	hCrossSec["xSection_pp_Pions"]->scaleW(1. / sow["sow_pp"]->sumW());
	
	hCrossSec["xSection_dAu_GammaGamma"]->scaleW(1. / sow["sow_dAUc20"]->sumW());
	hCrossSec["xSection_dAu_GammaGamma"]->scaleW(1. / sow["sow_dAUc40"]->sumW());
	hCrossSec["xSection_dAu_GammaGamma"]->scaleW(1. / sow["sow_dAUc60"]->sumW());
	hCrossSec["xSection_dAu_GammaGamma"]->scaleW(1. / sow["sow_dAUc88"]->sumW());
	
	hCrossSec["xSection_dAu_Pions"]->scaleW(1. / sow["sow_dAUc20"]->sumW());
	hCrossSec["xSection_dAu_Pions"]->scaleW(1. / sow["sow_dAUc40"]->sumW());
	hCrossSec["xSection_dAu_Pions"]->scaleW(1. / sow["sow_dAUc60"]->sumW());
	hCrossSec["xSection_dAu_Pions"]->scaleW(1. / sow["sow_dAUc88"]->sumW());
	
	hRatiopp["ratioEtaPi0_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
	
	hRatioAuAu["ratioEtaPi0_c20"]->scaleW(1. / sow["sow_AUAUc20"]->sumW());
	hRatioAuAu["ratioEtaPi0_c60"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
	hRatioAuAu["ratioEtaPi0_c92"]->scaleW(1. / sow["sow_AUAUc92"]->sumW());
	hRatioAuAu["ratioEtaPi0_cMB_dAu"]->scaleW(1. / sow["sow_dAUc20"]->sumW());
	hRatioAuAu["ratioEtaPi0_cMB_dAu"]->scaleW(1. / sow["sow_dAUc40"]->sumW());
	hRatioAuAu["ratioEtaPi0_cMB_dAu"]->scaleW(1. / sow["sow_dAUc60"]->sumW());
	hRatioAuAu["ratioEtaPi0_cMB_dAu"]->scaleW(1. / sow["sow_dAUc88"]->sumW());
	
	hRatiodAu["ratioEtaPi0_c20"]->scaleW(1. / sow["sow_dAUc20"]->sumW());
	hRatiodAu["ratioEtaPi0_c40"]->scaleW(1. / sow["sow_dAUc40"]->sumW());
	hRatiodAu["ratioEtaPi0_c60"]->scaleW(1. / sow["sow_dAUc60"]->sumW());
	hRatiodAu["ratioEtaPi0_c88"]->scaleW(1. / sow["sow_dAUc88"]->sumW());
	
	hRdAu20["RdAu_c20"]->scaleW(1. / sow["sow_dAUc20"]->sumW());
	hRdAu40["RdAu_c40"]->scaleW(1. / sow["sow_dAUc40"]->sumW());
	hRdAu60["RdAu_c60"]->scaleW(1. / sow["sow_dAUc60"]->sumW());
	hRdAu88["RdAu_c88"]->scaleW(1. / sow["sow_dAUc88"]->sumW());
	hRdAuMB["RdAu_cMB"]->scaleW(1. / sow["sow_dAUc20"]->sumW());
	hRdAuMB["RdAu_cMB"]->scaleW(1. / sow["sow_dAUc40"]->sumW());
	hRdAuMB["RdAu_cMB"]->scaleW(1. / sow["sow_dAUc60"]->sumW());
	hRdAuMB["RdAu_cMB"]->scaleW(1. / sow["sow_dAUc88"]->sumW());
	divide(hRdAu20["RdAu_c20"], hCrossSec["xSection_pp_GammaGamma"], RatiodAu20["dAuc/pp_c20"]);
	divide(hRdAu40["RdAu_c40"], hCrossSec["xSection_pp_GammaGamma"], RatiodAu40["dAuc/pp_c40"]);
	divide(hRdAu60["RdAu_c60"], hCrossSec["xSection_pp_GammaGamma"], RatiodAu60["dAuc/pp_c60"]);
	divide(hRdAu88["RdAu_c88"], hCrossSec["xSection_pp_GammaGamma"], RatiodAu88["dAuc/pp_c88"]);
	divide(hRdAuMB["RdAu_cMB"], hCrossSec["xSection_pp_GammaGamma"], RatiodAuMB["dAuc/pp_cMB"]);
	
	hRAuAu["RAuAu_c20"]->scaleW(1. / sow["sow_AUAUc20"]->sumW());
	hRAuAu["RAuAu_c60"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
	hRAuAu["RAuAu_c92"]->scaleW(1. / sow["sow_AUAUc92"]->sumW());
	divide(hRAuAu["RAuAu_c20"], hCrossSec["xSection_pp_GammaGamma"], RatioAuAu["AuAuc/pp_c20"]);
	divide(hRAuAu["RAuAu_c60"], hCrossSec["xSection_pp_GammaGamma"], RatioAuAu["AuAuc/pp_c60"]);
	divide(hRAuAu["RAuAu_c92"], hCrossSec["xSection_pp_GammaGamma"], RatioAuAu["AuAuc/pp_c92"]);
	
	hdAuYields["pTyields_c20"]->scaleW(1. / sow["sow_dAUc20"]->sumW());
	hdAuYields["pTyields_c40"]->scaleW(1. / sow["sow_dAUc40"]->sumW());
	hdAuYields["pTyields_c60"]->scaleW(1. / sow["sow_dAUc60"]->sumW());
	hdAuYields["pTyields_c88"]->scaleW(1. / sow["sow_dAUc88"]->sumW());
	
	hAuAuYields["pTyields_c20"]->scaleW(1. / sow["sow_AUAUc20"]->sumW());
	hAuAuYields["pTyields_c40"]->scaleW(1. / sow["sow_AUAUc40"]->sumW());
	hAuAuYields["pTyields_c92"]->scaleW(1. / sow["sow_AUAUc92"]->sumW());
	hAuAuYields["pTyields_cMB"]->scaleW(1. / sow["sow_AUAUc20"]->sumW());
	hAuAuYields["pTyields_cMB"]->scaleW(1. / sow["sow_AUAUc40"]->sumW());
	hAuAuYields["pTyields_cMB"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
	hAuAuYields["pTyields_cMB"]->scaleW(1. / sow["sow_AUAUc92"]->sumW());
	
	RatioAuAu["AuAuc/pp_c20"]->scaleY(1. / 6.14);
	RatioAuAu["AuAuc/pp_c60"]->scaleY(1. / 4.6);
	RatioAuAu["AuAuc/pp_c92"]->scaleY(1. / 0.3);
	
	RatiodAu20["dAuc/pp_c20"]->scaleY(1. / 0.36);
	RatiodAu40["dAuc/pp_c40"]->scaleY(1. / 0.25);
	RatiodAu60["dAuc/pp_c60"]->scaleY(1. / 0.17);
	RatiodAu88["dAuc/pp_c88"]->scaleY(1. / 0.078);
	RatiodAuMB["dAuc/pp_cMB"]->scaleY(1. / 0.2);

    }

	map<string, Histo1DPtr> hCrossSec;
	map<string, Histo1DPtr> hRatiopp;
	map<string, Histo1DPtr> hRdAu20;
	map<string, Histo1DPtr> hRdAu40;
	map<string, Histo1DPtr> hRdAu60;
	map<string, Histo1DPtr> hRdAu88;
	map<string, Histo1DPtr> hRdAuMB;
	map<string, Histo1DPtr> hRAuAu;
	map<string, Histo1DPtr> hRatioAuAu;
	map<string, Scatter2DPtr> RatioAuAu;
	map<string, Histo1DPtr> hdAuYields;
	map<string, Histo1DPtr> hAuAuYields;
	map<string, Histo1DPtr> hRatiodAu;
	map<string, Scatter2DPtr> RatiodAu20;
	map<string, Scatter2DPtr> RatiodAu40;
	map<string, Scatter2DPtr> RatiodAu60;
	map<string, Scatter2DPtr> RatiodAu88;
	map<string, Scatter2DPtr> RatiodAuMB;
	
	map<string, CounterPtr> sow;
	string beamOpt;
	enum CollisionSystem { pp, AuAu200, dAu200 };
	CollisionSystem collSys;

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2006_I0611006);

}
