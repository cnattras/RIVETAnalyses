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
#include "/home/antoniosilva/RivetWorkdir/Christal/1304.3410/RHICCentrality.hh" //external header for Centrality calculation
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
    	//Particles: eta, pi^+, pi^-, pi^0, gamma (Not listed: deuteron, Au)
    	std::initializer_list<int> pdgIds = { 221, 211, -211, 111, 22};
    	const PrimaryParticles fs(Cuts::abseta < 0.35 && Cuts::abscharge > 0);
    	declare(fs, "fs");
    	
    	const PrimaryParticles ns(Cuts::abseta < 0.35 && Cuts::abscharge == 0);
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
    	string refname1 = mkAxisCode(1, 1, 1);
	const Scatter2D& refdata1 = refData(refname1);
	book(hCrossSec["xSection_pp_GammaGamma"], refname1 + "_pp", refdata1);
	
	//Cross sections in pp (fig 12.2)
	string refname2 = mkAxisCode(1, 1, 2);
	const Scatter2D& refdata2 = refData(refname2);
	book(hCrossSec["xSection_pp_Pions"], refname2 + "_pp", refdata2);
	
	//Cross sections in dAu (fig 13.1)
    	string refname3 = mkAxisCode(2, 1, 1);
	const Scatter2D& refdata3 = refData(refname3);
	book(hCrossSec["xSection_dAu_GammaGamma"], refname3 + "_dAu", refdata3);
	
	//Cross sections in dAu (fig 13.2)
	string refname4 = mkAxisCode(2, 1, 2);
	const Scatter2D& refdata4 = refData(refname4);
	book(hCrossSec["xSection_dAu_Pions"], refname4 + "_dAu", refdata4);
	
	//Ratio eta/pi^0 in pp (fig 18)
	string refname5 = mkAxisCode(9, 1, 1);
	const Scatter2D& refdata5 = refData(refname5);
	book(hRatiopp["ratioEtaPi0_pp"], refname5 + "_pp", refdata5);
	
	//R_dAu 0-20 (fig 16)
	string refname6 = mkAxisCode(7, 1, 2);
	const Scatter2d& refdata6 = refData(refname6);
	book(hRdAu20["RdAu_c20"], refname6 + "_dAu", refdata6);
	
	//R_dAu 20-40 (fig 16)
	string refname7 = mkAxisCode(7, 1, 3);
	const Scatter2d& refdata7 = refData(refname7);
	book(hRdAu40["RdAu_c40"], refname7 + "_dAu", refdata7);
	
	//R_dAu 40-60 (fig 16)
	string refname8 = mkAxisCode(7, 1, 4);
	const Scatter2d& refdata8 = refData(refname8);
	book(hRdAu60["RdAu_c60"], refname8 + "_dAu", refdata8);
	
	//R_dAu 60-88 (fig 16)
	string refname9 = mkAxisCode(7, 1, 5);
	const Scatter2d& refdata9 = refData(refname9);
	book(hRdAu88["RdAu_c88"], refname9 + "_dAu", refdata9);
	
	//R_dAu 0-88(MB) (fig 16)
	string refname10 = mkAxisCode(7, 1, 1);
	const Scatter2d& refdata10 = refData(refname10);
	book(hRdAuMB["RdAu_cMB"], refname10 + "_dAu", refdata10);
	
	//R_AuAu 0-20 (fig 17)
	string refname15 = mkAxisCode(8, 1, 1);
	const Scatter2d& refdata15 = refData(refname15);
	book(hRAuAu["RAuAu_c20"], refname15 + "_AuAu", refdata15);
	
	//R_AuAu 20-60 (fig 17)
	string refname16 = mkAxisCode(8, 1, 2);
	const Scatter2d& refdata16 = refData(refname16);
	book(hRAuAu["RAuAu_c60"], refname16 + "_AuAu", refdata16);
	
	//R_AuAu 60-92 (fig 17)
	string refname17 = mkAxisCode(8, 1, 3);
	const Scatter2d& refdata17 = refData(refname17);
	book(hRAuAu["RAuAu_c92"], refname17 + "_AuAu", refdata17);
	
	//Ratio eta/pi^0 in Au+Au 0-92 MB (fig 20) ***UNUSED IN HISTO***
	//string refname31 = mkAxisCode(11, 1, 1);
	//const Scatter2d& refdata31 = refData(refname31);
	//book(hRatioAuAu["ratioEtaPi0_cMB"], refname31 + "_AuAu", refdata31);
	
	//Ratio eta/pi^0 in Au+Au 0-20 (fig 20)
	string refname11 = mkAxisCode(11, 1, 2);
	const Scatter2d& refdata11 = refData(refname11);
	book(hRatioAuAu["ratioEtaPi0_c20"], refname11 + "_AuAu", refdata11);
	
	//Ratio eta/pi^0 in Au+Au 20-60 (fig 20)
	string refname12 = mkAxisCode(11, 1, 3);
	const Scatter2d& refdata12 = refData(refname12);
	book(hRatioAuAu["ratioEtaPi0_c60"], refname12 + "_AuAu", refdata12);
	
	//Ratio eta/pi^0 in Au+Au 60-92 (fig 20)
	string refname13 = mkAxisCode(11, 1, 4);
	const Scatter2d& refdata13 = refData(refname13);
	book(hRatioAuAu["ratioEtaPi0_c92"], refname13 + "_AuAu", refdata13);
	
	//Ratio eta/pi^0 in Au+Au d+Au MB comparison (fig 20) 
	string refname14 = mkAxisCode(10, 1, 1);
	const Scatter2d& refdata14 = refData(refname14);
	book(hRatioAuAu["ratioEtaPi0_cMB_dAu"], refname11 + "_dAu", refdata11);
	
	//invariant yields in d+Au 0-20 (fig 14)
	string refname18 = mkAxisCode(5, 1, 1);
	const Scatter2d& refdata18 = refData(refname18);
	book(hdAuYields["pTyields_c20"], refname18 + "_dAu", refdata18);
	
	//invariant yields in d+Au 20-40 (fig 14)
	string refname19 = mkAxisCode(5, 1, 2);
	const Scatter2d& refdata19 = refData(refname19);
	book(hdAuYields["pTyields_c40"], refname19 + "_dAu", refdata19);
	
	//invariant yields in d+Au 40-60 (fig 14)
	string refname20 = mkAxisCode(5, 1, 3);
	const Scatter2d& refdata20 = refData(refname20);
	book(hdAuYields["pTyields_c60"], refname20 + "_dAu", refdata20);
	
	//invariant yields in d+Au 60-88 (fig 14)
	string refname21 = mkAxisCode(5, 1, 4);
	const Scatter2d& refdata21 = refData(refname21);
	book(hdAuYields["pTyields_c88"], refname21 + "_dAu", refdata21);
	
	//invariant yields in Au+Au 0-92 (fig 15)
	string refname22 = mkAxisCode(6, 1, 1);
	const Scatter2d& refdata22 = refData(refname22);
	book(hAuAuYields["pTyields_cMB"], refname22 + "_AuAu", refdata22);
	
	//invariant yields in Au+Au 0-20 (fig 15)
	string refname23 = mkAxisCode(6, 1, 2);
	const Scatter2d& refdata23 = refData(refname23);
	book(hAuAuYields["pTyields_c20"], refname23 + "_AuAu", refdata23);
	
	//invariant yields in Au+Au 20-40 (fig 15)
	string refname24 = mkAxisCode(6, 1, 3);
	const Scatter2d& refdata24 = refData(refname24);
	book(hAuAuYields["pTyields_c40"], refname24 + "_AuAu", refdata24);
	
	//invariant yields in Au+Au 60-92 (fig 15)
	string refname25 = mkAxisCode(6, 1, 4);
	const Scatter2d& refdata25 = refData(refname25);
	book(hAuAuYields["pTyields_c92"], refname25 + "_AuAu", refdata25);
	
	//Ratio eta/pi^0 in d+Au 0-88 (fig 19)  ***UNUSED IN HISTO***
	//string refname26 = mkAxisCode(10, 1, 1);
	//const Scatter2d& refdata26 = refData(refname26);
	//book(hRatiodAu["ratioEtaPi0_cMB"], refname26 + "_dAu", refdata26);
	
	//Ratio eta/pi^0 in d+Au 0-20 (fig 19)
	string refname27 = mkAxisCode(10, 1, 2);
	const Scatter2d& refdata27 = refData(refname27);
	book(hRatiodAu["ratioEtaPi0_c20"], refname27 + "_dAu", refdata27);
	
	//Ratio eta/pi^0 in d+Au 19-40 (fig 19)
	string refname28 = mkAxisCode(10, 1, 3);
	const Scatter2d& refdata28 = refData(refname28);
	book(hRatiodAu["ratioEtaPi0_c40"], refname28 + "_dAu", refdata28);
	
	//Ratio eta/pi^0 in d+Au 40-60 (fig 19)
	string refname29 = mkAxisCode(10, 1, 4);
	const Scatter2d& refdata29 = refData(refname29);
	book(hRatiodAu["ratioEtaPi0_c60"], refname29 + "_dAu", refdata29);
	
	//Ratio eta/pi^0 in d+Au 60-88 (fig 19)
	string refname30 = mkAxisCode(10, 1, 5);
	const Scatter2d& refdata30 = refData(refname30);
	book(hRatiodAu["ratioEtaPi0_c88"], refname30 + "_dAu", refdata30);

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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
					hAuAuYields["pTyields_c20_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c40_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_c92_AuAu"]->fill(partPt, pt_weight);
					hAuAuYields["pTyields_cMB_AuAu"]->fill(partPt, pt_weight);
					
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
	divide(hRdAu20["RdAu_c20"], hCrossSec["xSection_pp_GammaGamma"], RatiodAu20["dAuc/pp"]);
	divide(hRdAu40["RdAu_c40"], hCrossSec["xSection_pp_GammaGamma"], RatiodAu40["dAuc/pp"]);
	divide(hRdAu60["RdAu_c60"], hCrossSec["xSection_pp_GammaGamma"], RatiodAu60["dAuc/pp"]);
	divide(hRdAu88["RdAu_c88"], hCrossSec["xSection_pp_GammaGamma"], RatiodAu88["dAuc/pp"]);
	divide(hRdAuMB["RdAu_cMB"], hCrossSec["xSection_pp_GammaGamma"], RatiodAuMB["dAuc/pp"]);
	
	hRAuAu["RAuAu_c20"]->scaleW(1. / sow["sow_AUAUc20"]->sumW());
	hRAuAu["RAuAu_c60"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
	hRAuAu["RAuAu_c92"]->scaleW(1. / sow["sow_AUAUc92"]->sumW());
	divide(hRAuAu["RAuAu_c20"], hCrossSec["xSection_pp_GammaGamma"], RatioAuAu["AuAuc/pp"]);
	divide(hRAuAu["RAuAu_c60"], hCrossSec["xSection_pp_GammaGamma"], RatioAuAu["AuAuc/pp"]);
	divide(hRAuAu["RAuAu_c92"], hCrossSec["xSection_pp_GammaGamma"], RatioAuAu["AuAuc/pp"]);
	
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
	
	hRAuAu["RAuAu_c20"]->scaleY(1. / 6.14);
	hRAuAu["RAuAu_c60"]->scaleY(1. / 4.6);
	hRAuAu["RAuAu_c92"]->scaleY(1. / 0.3);
	
	hRdAu20["RdAu_c20"]->scaleY(1. / 0.36);
	hRdAu20["RdAu_c40"]->scaleY(1. / 0.25);
	hRdAu20["RdAu_c60"]->scaleY(1. / 0.17);
	hRdAu20["RdAu_c88"]->scaleY(1. / 0.078);
	hRdAu20["RdAu_cMB"]->scaleY(1. / 0.2);

    }

	map<string, Scatter2DPtr> hCrossSec;
	map<string, Scatter2DPtr> hRatiopp;
	map<string, Scatter2DPtr> hRdAu20;
	map<string, Scatter2DPtr> hRdAu40;
	map<string, Scatter2DPtr> hRdAu60;
	map<string, Scatter2DPtr> hRdAu88;
	map<string, Scatter2DPtr> hRdAuMB;
	map<string, Scatter2DPtr> hRAuAu;
	map<string, Scatter2DPtr> hRatioAuAu;
	map<string, Scatter2DPtr> RatioAuAu;
	map<string, Scatter2DPtr> hdAuYields;
	map<string, Scatter2DPtr> hAuAuYields;
	map<string, Scatter2DPtr> hRatiodAu;
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
