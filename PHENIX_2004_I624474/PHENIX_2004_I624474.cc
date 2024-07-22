// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "../Centralities/RHICCentrality.hh"
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES


namespace Rivet {


	/// @brief Add a short analysis description here
	class PHENIX_2004_I624474 : public Analysis {
	public:

	/// Constructor
	DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2004_I624474);

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
	
	// Particles: pi+, pi-, K+, K-, p, pbar, pi0
	
	std::initializer_list<int> pdgIds ={211, 321, 2212, -211, -321, -2212};

	const PrimaryParticles cp(pdgIds, Cuts::abseta < .35 && Cuts::abscharge > 0);
	//const ALICE::PrimaryParticles cp(pdgIds, Cuts::abseta < .35 && Cuts::abscharge > 0);
	
	declare(cp, "cp");

	const UnstableParticles np(Cuts::abseta < .35 && Cuts::abspid == 111);
	declare(np, "np");

	declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

	// Book histograms
	

	//___PiPlus yields___

	//Invariant yield of piplus minimum bias
	book(hAUAU_Yields["Piplusmin"],1,1,1);

	//Invariant yield of piplus 0-5%
	book(hAUAU_Yields["Piplus5"],1,1,2);

	//Invariant yield of piplus 5-10%
	book(hAUAU_Yields["Piplus10"],1,1,3);

	//Invariant yield of piplus 10-15% 
	book(hAUAU_Yields["Piplus15"],1,1,4);

	//Invariant yield of piplus 15-20%
	book(hAUAU_Yields["Piplus20"],1,1,5);

	//Invariant yield of piplus 20-30%
	book(hAUAU_Yields["Piplus30"],1,1,6);

	//Invariant yield of piplus 30-40%
	book(hAUAU_Yields["Piplus40"],1,1,7);

	//Invariant yield of piplus 40-50%
	book(hAUAU_Yields["Piplus50"],1,1,8);

	//Invariant yield of piplus 50-60%
	book(hAUAU_Yields["Piplus60"],1,1,9);

	//Invariant yield of piplus 60-70%
	book(hAUAU_Yields["Piplus70"],1,1,10);

	//Invariant yield of piplus 70-80%
	book(hAUAU_Yields["Piplus80"],1,1,11);

	//Invariant yield of piplus 80-92%
	book(hAUAU_Yields["Piplus92"],1,1,12);

	//Invariant yield of piplus 60-92%
	book(hAUAU_Yields["Piplus60_92"],1,1,13);


	//___PiMinus yields___

	//Invariant yield of piminus minimum bias
	book(hAUAU_Yields["Piminusmin"],2,1,1);

	//Invariant yield of piminus 0-5%
	book(hAUAU_Yields["Piminus5"],2,1,2);

	//Invariant yield of piminus 5-10%
	book(hAUAU_Yields["Piminus10"],2,1,3);

	//Invariant yield of piminus 10-15% 
	book(hAUAU_Yields["Piminus15"],2,1,4);

	//Invariant yield of piminus 15-20%
	book(hAUAU_Yields["Piminus20"],2,1,5);

	//Invariant yield of piminus 20-30%
	book(hAUAU_Yields["Piminus30"],2,1,6);

	//Invariant yield of piminus 30-40%
	book(hAUAU_Yields["Piminus40"],2,1,7);

	//Invariant yield of piminus 40-50%
	book(hAUAU_Yields["Piminus50"],2,1,8);

	//Invariant yield of piminus 50-60%
	book(hAUAU_Yields["Piminus60"],2,1,9);

	//Invariant yield of piminus 60-70%
	book(hAUAU_Yields["Piminus70"],2,1,10);

	//Invariant yield of piminus 70-80%
	book(hAUAU_Yields["Piminus80"],2,1,11);

	//Invariant yield of piminus 80-92%
	book(hAUAU_Yields["Piminus92"],2,1,12);

	//Invariant yield of piminus 60-92%
	book(hAUAU_Yields["Piminus60_92"],2,1,13);


	//___KPlus yields___

	//Invariant yield of Kplus minimum bias
	book(hAUAU_Yields["Kplusmin"],3,1,1);

	//Invariant yield of Kplus 0-5%
	book(hAUAU_Yields["Kplus5"],3,1,2);

	//Invariant yield of Kplus 5-10%
	book(hAUAU_Yields["Kplus10"],3,1,3);

	//Invariant yield of Kplus 10-15% 
	book(hAUAU_Yields["Kplus15"],3,1,4);

	//Invariant yield of Kplus 15-20%
	book(hAUAU_Yields["Kplus20"],3,1,5);

	//Invariant yield of Kplus 20-30%
	book(hAUAU_Yields["Kplus30"],3,1,6);

	//Invariant yield of Kplus 30-40%
	book(hAUAU_Yields["Kplus40"],3,1,7);

	//Invariant yield of Kplus 40-50%
	book(hAUAU_Yields["Kplus50"],3,1,8);

	//Invariant yield of Kplus 50-60%
	book(hAUAU_Yields["Kplus60"],3,1,9);

	//Invariant yield of Kplus 60-70%
	book(hAUAU_Yields["Kplus70"],3,1,10);

	//Invariant yield of Kplus 70-80%
	book(hAUAU_Yields["Kplus80"],3,1,11);

	//Invariant yield of Kplus 80-92%
	book(hAUAU_Yields["Kplus92"],3,1,12);

	//Invariant yield of Kplus 60-92%
	book(hAUAU_Yields["Kplus60_92"],3,1,13);


	//___KMinus yields___

	//Invariant yield of Kminus minimum bias
	book(hAUAU_Yields["Kminusmin"],4,1,1);

	//Invariant yield of Kminus 0-5%
	book(hAUAU_Yields["Kminus5"],4,1,2);

	//Invariant yield of Kminus 5-10%
	book(hAUAU_Yields["Kminus10"],4,1,3);

	//Invariant yield of Kminus 10-15% 
	book(hAUAU_Yields["Kminus15"],4,1,4);

	//Invariant yield of Kminus 15-20%
	book(hAUAU_Yields["Kminus20"],4,1,5);

	//Invariant yield of Kminus 20-30%
	book(hAUAU_Yields["Kminus30"],4,1,6);

	//Invariant yield of Kminus 30-40%
	book(hAUAU_Yields["Kminus40"],4,1,7);

	//Invariant yield of Kminus 40-50%
	book(hAUAU_Yields["Kminus50"],4,1,8);

	//Invariant yield of Kminus 50-60%
	book(hAUAU_Yields["Kminus60"],4,1,9);

	//Invariant yield of Kminus 60-70%
	book(hAUAU_Yields["Kminus70"],4,1,10);

	//Invariant yield of Kminus 70-80%
	book(hAUAU_Yields["Kminus80"],4,1,11);

	//Invariant yield of Kminus 80-92%
	book(hAUAU_Yields["Kminus92"],4,1,12);

	//Invariant yield of Kminus 60-92%
	book(hAUAU_Yields["Kminus60_92"],4,1,13);


	//___Proton yields___

	//Invariant yield of protons minimum bias
	book(hAUAU_Yields["Protonsmin"],5,1,1);

	//Invariant yield of protons 0-5%
	book(hAUAU_Yields["Protons5"],5,1,2);

	//Invariant yield of protons 5-10%
	book(hAUAU_Yields["Protons10"],5,1,3);

	//Invariant yield of protons 10-15% 
	book(hAUAU_Yields["Protons15"],5,1,4);

	//Invariant yield of protons 15-20%
	book(hAUAU_Yields["Protons20"],5,1,5);

	//Invariant yield of protons 20-30%
	book(hAUAU_Yields["Protons30"],5,1,6);

	//Invariant yield of protons 30-40%
	book(hAUAU_Yields["Protons40"],5,1,7);

	//Invariant yield of protons 40-50%
	book(hAUAU_Yields["Protons50"],5,1,8);

	//Invariant yield of protons 50-60%
	book(hAUAU_Yields["Protons60"],5,1,9);

	//Invariant yield of protons 60-70%
	book(hAUAU_Yields["Protons70"],5,1,10);

	//Invariant yield of protons 70-80%
	book(hAUAU_Yields["Protons80"],5,1,11);

	//Invariant yield of protons 80-92%
	book(hAUAU_Yields["Protons92"],5,1,12);

	//Invariant yield of protons 60-92%
	book(hAUAU_Yields["Protons60_92"],5,1,13);


	//___Pbar yields___

	//Invariant yield of pbar minimum bias
	book(hAUAU_Yields["Pbarmin"],6,1,1);

	//Invariant yield of pbar 0-5%
	book(hAUAU_Yields["Pbar5"],6,1,2);

	//Invariant yield of pbar 5-10%
	book(hAUAU_Yields["Pbar10"],6,1,3);

	//Invariant yield of pbar 10-15% 
	book(hAUAU_Yields["Pbar15"],6,1,4);

	//Invariant yield of pbar 15-20%
	book(hAUAU_Yields["Pbar20"],6,1,5);

	//Invariant yield of pbar 20-30%
	book(hAUAU_Yields["Pbar30"],6,1,6);

	//Invariant yield of pbar 30-40%
	book(hAUAU_Yields["Pbar40"],6,1,7);

	//Invariant yield of pbar 40-50%
	book(hAUAU_Yields["Pbar50"],6,1,8);

	//Invariant yield of pbar 50-60%
	book(hAUAU_Yields["Pbar60"],6,1,9);

	//Invariant yield of pbar 60-70%
	book(hAUAU_Yields["Pbar70"],6,1,10);

	//Invariant yield of pbar 70-80%
	book(hAUAU_Yields["Pbar80"],6,1,11);

	//Invariant yield of pbar 80-92%
	book(hAUAU_Yields["Pbar92"],7,1,1);

	//Invariant yield of pbar 60-92%
	book(hAUAU_Yields["Pbar60_92"],6,1,12);



	//____Counters____

	book(sow["sow_AUAUmin"],"_sow_AUAUmin");
	book(sow["sow_AUAU5"],"_sow_AUAU5");
	book(sow["sow_AUAU10"],"_sow_AUAU10");
	book(sow["sow_AUAU15"],"_sow_AUAU15");
	book(sow["sow_AUAU20"],"_sow_AUAU20");
	book(sow["sow_AUAU30"],"_sow_AUAU30");
	book(sow["sow_AUAU40"],"_sow_AUAU40");
	book(sow["sow_AUAU50"],"_sow_AUAU50");
	book(sow["sow_AUAU60"],"_sow_AUAU60");
	book(sow["sow_AUAU70"],"_sow_AUAU70");
	book(sow["sow_AUAU80"],"_sow_AUAU80");
	book(sow["sow_AUAU92"],"_sow_AUAU92");
	book(sow["sow_AUAU_60_92"],"_sow_AUAU_60_92");
	book(sow["sow_AUAU0_10"],"_sow_AUAU0_10");
	book(sow["sow_AUAUall"],"_sow_AUAUall");

	//____Rcp____
	
	string refnameRcpPi = mkAxisCode(21,1,1);
	const Scatter2D& refdataRcpPi = refData(refnameRcpPi);
	book(hPi["AUAU0_10"], refnameRcpPi + "_0_10Pion", refdataRcpPi);
	book(hPi["AUAU60_92"], refnameRcpPi + "_60_92Pion", refdataRcpPi);
	book(hRcp["Pions"], refnameRcpPi);

	string refnameRcpK = mkAxisCode(22,1,1);
	const Scatter2D& refdataRcpK = refData(refnameRcpK);
	book(hK["AUAU0_10"], refnameRcpK + "_0_10Kaon", refdataRcpK);
	book(hK["AUAU60_92"], refnameRcpK + "_60_92Kaon", refdataRcpK);
	book(hRcp["Kaons"], refnameRcpK);

	string refnameRcpP = mkAxisCode(23,1,1);
	const Scatter2D& refdataRcpP = refData(refnameRcpP);
	book(hP["AUAU0_10"], refnameRcpP + "_0_10Proton", refdataRcpP);
	book(hP["AUAU60_92"], refnameRcpP + "_60_92Proton", refdataRcpP);
	book(hRcp["Pbar+P"], refnameRcpP);

	string refnameRcpPi0 = mkAxisCode(24,1,1);
	const Scatter2D& refdataRcpPi0 = refData(refnameRcpPi0);
	book(hPi0["AUAU0_10"], refnameRcpPi0 + "_0_10Pion0", refdataRcpPi0);
	book(hPi0["AUAU60_92"], refnameRcpPi0 + "_60_92Pion0", refdataRcpPi0);
	book(hRcp["Pi0"], refnameRcpPi0);
	
	//____ Rcp 1.5GeV+____
	
	string refnameRcpPi1_5 = mkAxisCode(25,1,1);
	const Scatter2D& refdataRcpPi1_5 = refData(refnameRcpPi1_5);
	book(hPi["AUAU0_10GeV1_5"], refnameRcpPi1_5 + "_0_10Pion1_5", refdataRcpPi1_5);
	book(hPi["AUAU60_92GeV1_5"], refnameRcpPi1_5 + "_60_92Pion1_5", refdataRcpPi1_5);
	book(hRcp["Pions1_5"], refnameRcpPi1_5);
	//These two need to be added to yoda and uncommented
	string refnameRcpK1_5 = mkAxisCode(25,1,2);
	const Scatter2D& refdataRcpK1_5 = refData(refnameRcpK1_5);
	book(hK["AUAU0_10GeV1_5"], refnameRcpK1_5 + "_0_10Kaon1_5", refdataRcpK1_5);
	book(hK["AUAU60_92GeV1_5"], refnameRcpK1_5 + "_60_92Kaon1_5", refdataRcpK1_5);
	book(hRcp["Kaons1_5"], refnameRcpK1_5);

	string refnameRcpP1_5 = mkAxisCode(25,1,3);
	const Scatter2D& refdataRcpP1_5 = refData(refnameRcpP1_5);
	book(hP["AUAU0_10GeV1_5"], refnameRcpP1_5 + "_0_10Proton1_5", refdataRcpP1_5);
	book(hP["AUAU60_92GeV1_5"], refnameRcpP1_5 + "_60_92Proton1_5", refdataRcpP1_5);
	book(hRcp["Pbar+P1_5"], refnameRcpP1_5);

	//____mean Pt____
	
	book(pmeanPt["Piplus"],12,1,1);
	book(pmeanPt["Piminus"],12,1,2);
	book(pmeanPt["Kplus"],12,1,3);
	book(pmeanPt["Kminus"],12,1,4);
	book(pmeanPt["Protons"],12,1,5);
	book(pmeanPt["Pbar"],12,1,6);


	//____dN/dy____
	
	book(pdN_dy["Piplus"],13,1,1);
	book(pdN_dy["Piminus"],13,1,2);
	book(pdN_dy["Kplus"],13,1,3);
	book(pdN_dy["Kminus"],13,1,4);
	book(pdN_dy["Protons"],13,1,5);
	book(pdN_dy["Pbar"],13,1,6);


	//____Ratios____
	
	string refnameratio0_5PiPi = mkAxisCode(14,1,1);	//Pi-/Pi+
	const Scatter2D& refdataratio0_5PiPi = refData(refnameratio0_5PiPi);
	book(hPiPi["AUAU0_5Piminus"], refnameratio0_5PiPi + "_0_5Piminus", refdataratio0_5PiPi);
	book(hPiPi["AUAU0_5Piplus"], refnameratio0_5PiPi + "_0_5Piplus", refdataratio0_5PiPi);
	book(hRatio["Piminus_Piplus0_5"], refnameratio0_5PiPi);

	string refnameratio60_92PiPi = mkAxisCode(14,1,2);
	const Scatter2D& refdataratio60_92PiPi = refData(refnameratio60_92PiPi);
	book(hPiPi["AUAU60_92Piminus"], refnameratio60_92PiPi + "_60_92Piminus", refdataratio60_92PiPi);
	book(hPiPi["AUAU60_92Piplus"], refnameratio60_92PiPi + "_60_92Piplus", refdataratio60_92PiPi);
	book(hRatio["Piminus_Piplus60_92"], refnameratio60_92PiPi);

	string refnameratio0_5KK = mkAxisCode(15,1,1);	//Kaon-/Kaon+
	const Scatter2D& refdataratio0_5KK = refData(refnameratio0_5KK);
	book(hKK["AUAU0_5Kminus"], refnameratio0_5KK + "_0_5Kminus", refdataratio0_5KK);
	book(hKK["AUAU0_5Kplus"], refnameratio0_5KK + "_0_5Kplus", refdataratio0_5KK);
	book(hRatio["Kiminus_Kiplus0_5"], refnameratio0_5KK);

	string refnameratio60_92KK = mkAxisCode(15,1,2);
	const Scatter2D& refdataratio60_92KK = refData(refnameratio60_92KK);
	book(hKK["AUAU60_92Kminus"], refnameratio60_92KK + "_60_92Kminus", refdataratio60_92KK);
	book(hKK["AUAU60_92Kplus"], refnameratio60_92KK + "_60_92Kplus", refdataratio60_92KK);
	book(hRatio["Kiminus_Kiplus60_92"], refnameratio60_92KK);

	string refnameratio0_5PbarP = mkAxisCode(16,1,1);	//Antiproton/Proton
	const Scatter2D& refdataratio0_5PbarP = refData(refnameratio0_5PbarP);
	book(hPP["AUAU0_5Pbar"], refnameratio0_5PbarP + "_0_5Pbar", refdataratio0_5PbarP);
	book(hPP["AUAU0_5P"], refnameratio0_5PbarP + "_0_5P", refdataratio0_5PbarP);
	book(hRatio["Pbar_P0_5"], refnameratio0_5PbarP);

	string refnameratio60_92PbarP = mkAxisCode(16,1,2);
	const Scatter2D& refdataratio60_92PbarP = refData(refnameratio60_92PbarP);
	book(hPP["AUAU60_92Pbar"], refnameratio60_92PbarP + "_60_92Pbar", refdataratio60_92PbarP);
	book(hPP["AUAU60_92P"], refnameratio60_92PbarP + "_60_92P", refdataratio60_92PbarP);
	book(hRatio["Pbar_P60_92"], refnameratio60_92PbarP);

	string refnameratiominPbarP = mkAxisCode(16,1,3);
	const Scatter2D& refdataratiominPbarP = refData(refnameratiominPbarP);
	book(hPP["AUAUminPbar"], refnameratiominPbarP + "_minPbar", refdataratiominPbarP);
	book(hPP["AUAUminP"], refnameratiominPbarP + "_minP", refdataratiominPbarP);
	book(hRatio["Pbar_Pmin"], refnameratiominPbarP);

	string refnameratio0_5KplusPiplus = mkAxisCode(17,1,1);	//Kaon+/Pi+
	const Scatter2D& refdataratio0_5KplusPiplus = refData(refnameratio0_5KplusPiplus);
	book(hKPi["AUAU0_5Kplus"], refnameratio0_5KplusPiplus + "_0_5Kplus", refdataratio0_5KplusPiplus);
	book(hKPi["AUAU0_5Piplus"], refnameratio0_5KplusPiplus + "_0_5Piplus", refdataratio0_5KplusPiplus);
	book(hRatio["Kplus_Piplus0_5"], refnameratio0_5KplusPiplus);

	string refnameratio60_92KplusPiplus = mkAxisCode(17,1,2);
	const Scatter2D& refdataratio60_92KplusPiplus = refData(refnameratio60_92KplusPiplus);
	book(hKPi["AUAU60_92Kplus"], refnameratio60_92KplusPiplus + "_60_92Kplus", refdataratio60_92KplusPiplus);
	book(hKPi["AUAU60_92Piplus"], refnameratio60_92KplusPiplus + "_60_92Piplus", refdataratio60_92KplusPiplus);
	book(hRatio["Kplus_Piplus60_92"], refnameratio60_92KplusPiplus);
	
	string refnameratio0_5KminusPiminus = mkAxisCode(17,1,3);	//Kaon-/Pi-
	const Scatter2D& refdataratio0_5KminusPiminus = refData(refnameratio0_5KminusPiminus);
	book(hKPi["AUAU0_5Kminus"], refnameratio0_5KminusPiminus + "_0_5Kminus", refdataratio0_5KminusPiminus);
	book(hKPi["AUAU0_5Piminus"], refnameratio0_5KminusPiminus + "_0_5Piminus", refdataratio0_5KminusPiminus);
	book(hRatio["Kminus_Piminus0_5"], refnameratio0_5KminusPiminus);

	string refnameratio60_92KminusPiminus = mkAxisCode(17,1,4);
	const Scatter2D& refdataratio60_92KminusPiminus = refData(refnameratio60_92KminusPiminus);
	book(hKPi["AUAU60_92Kminus"], refnameratio60_92KminusPiminus + "_60_92Kminus", refdataratio60_92KminusPiminus);
	book(hKPi["AUAU60_92Piminus"], refnameratio60_92KminusPiminus + "_60_92Piminus", refdataratio60_92KminusPiminus);
	book(hRatio["Kminus_Piminus60_92"], refnameratio60_92KminusPiminus);

	string refnameratio0_10PPiplus = mkAxisCode(18,1,1);	//Proton/Pi+
	const Scatter2D& refdataratio0_10PPiplus = refData(refnameratio0_10PPiplus);
	book(hPPi["AUAU0_10ProtonPiplus"], refnameratio0_10PPiplus + "_0_10P", refdataratio0_10PPiplus);
	book(hPPi["AUAU0_10Piplus"], refnameratio0_10PPiplus + "_0_10Piplus", refdataratio0_10PPiplus);
	book(hRatio["P_Piplus0_10"], refnameratio0_10PPiplus);

	string refnameratio20_30PPiplus = mkAxisCode(18,1,2);
	const Scatter2D& refdataratio20_30PPiplus = refData(refnameratio20_30PPiplus);
	book(hPPi["AUAU20_30ProtonPiplus"], refnameratio20_30PPiplus + "_20_30P", refdataratio20_30PPiplus);
	book(hPPi["AUAU20_30Piplus"], refnameratio20_30PPiplus + "_20_30Piplus", refdataratio20_30PPiplus);
	book(hRatio["P_Piplus20_30"], refnameratio20_30PPiplus);

	string refnameratio60_92PPiplus = mkAxisCode(18,1,3);
	const Scatter2D& refdataratio60_92PPiplus = refData(refnameratio60_92PPiplus);
	book(hPPi["AUAU60_92ProtonPiplus"], refnameratio60_92PPiplus + "_60_92P", refdataratio60_92PPiplus);
	book(hPPi["AUAU60_92Piplus"], refnameratio60_92PPiplus + "_60_92Piplus", refdataratio60_92PPiplus);
	book(hRatio["P_Piplus60_92"], refnameratio60_92PPiplus);
	
	string refnameratio0_10PbarPiminus = mkAxisCode(18,1,4);	//Antiproton/Pi-
	const Scatter2D& refdataratio0_10PbarPiminus = refData(refnameratio0_10PbarPiminus);
	book(hPPi["AUAU0_10PbarPiminus"], refnameratio0_10PbarPiminus + "_0_10Pbar", refdataratio0_10PbarPiminus);
	book(hPPi["AUAU0_10Piminus"], refnameratio0_10PbarPiminus + "_0_10Piminus", refdataratio0_10PbarPiminus);
	book(hRatio["Pbar_Piminus0_10"], refnameratio0_10PbarPiminus);

	string refnameratio20_30PbarPiminus = mkAxisCode(18,1,5);
	const Scatter2D& refdataratio20_30PbarPiminus = refData(refnameratio20_30PbarPiminus);
	book(hPPi["AUAU20_30PbarPiminus"], refnameratio20_30PbarPiminus + "_20_30Pbar", refdataratio20_30PbarPiminus);
	book(hPPi["AUAU20_30Piminus"], refnameratio20_30PbarPiminus + "_20_30Piminus", refdataratio20_30PbarPiminus);
	book(hRatio["Pbar_Piminus20_30"], refnameratio20_30PbarPiminus);

	string refnameratio60_92PbarPiminus = mkAxisCode(18,1,6);
	const Scatter2D& refdataratio60_92PbarPiminus = refData(refnameratio60_92PbarPiminus);
	book(hPPi["AUAU60_92PbarPiminus"], refnameratio60_92PbarPiminus + "_60_92Pbar", refdataratio60_92PbarPiminus);
	book(hPPi["AUAU60_92Piminus"], refnameratio60_92PbarPiminus + "_60_92Piminus", refdataratio60_92PbarPiminus);
	book(hRatio["Pbar_Piminus60_92"], refnameratio60_92PbarPiminus);

	string refnameratio0_10PPi0 = mkAxisCode(19,1,1);	//Proton/Pi0
	const Scatter2D& refdataratio0_10PPi0 = refData(refnameratio0_10PPi0);
	book(hPPi["AUAU0_10P"], refnameratio0_10PPi0 + "_0_10P", refdataratio0_10PPi0);
	book(hPPi["AUAU0_10PPi0"], refnameratio0_10PPi0 + "_0_10Pi0", refdataratio0_10PPi0);
	book(hRatio["P_Pi00_10"], refnameratio0_10PPi0);

	string refnameratio20_30PPi0 = mkAxisCode(19,1,2);
	const Scatter2D& refdataratio20_30PPi0 = refData(refnameratio20_30PPi0);
	book(hPPi["AUAU20_30P"], refnameratio20_30PPi0 + "_20_30P", refdataratio20_30PPi0);
	book(hPPi["AUAU20_30PPi0"], refnameratio20_30PPi0 + "_20_30Pi0", refdataratio20_30PPi0);
	book(hRatio["P_Pi020_30"], refnameratio20_30PPi0);

	string refnameratio60_92PPi0 = mkAxisCode(19,1,3);
	const Scatter2D& refdataratio60_92PPi0 = refData(refnameratio60_92PPi0);
	book(hPPi["AUAU60_92P"], refnameratio60_92PPi0 + "_60_92P", refdataratio60_92PPi0);
	book(hPPi["AUAU60_92PPi0"], refnameratio60_92PPi0 + "_60_92Pi0", refdataratio60_92PPi0);
	book(hRatio["P_Pi060_92"], refnameratio60_92PPi0);
	
	string refnameratio0_10PbarPi0 = mkAxisCode(19,1,4);	//Antiproton/Pi0
	const Scatter2D& refdataratio0_10PbarPi0 = refData(refnameratio0_10PbarPi0);
	book(hPPi["AUAU0_10Pbar"], refnameratio0_10PbarPi0 + "_0_10Pbar", refdataratio0_10PbarPi0);
	book(hPPi["AUAU0_10PbarPi0"], refnameratio0_10PbarPi0 + "_0_10Pi0", refdataratio0_10PbarPi0);
	book(hRatio["Pbar_Pi00_10"], refnameratio0_10PbarPi0);

	string refnameratio20_30PbarPi0 = mkAxisCode(19,1,5);
	const Scatter2D& refdataratio20_30PbarPi0 = refData(refnameratio20_30PbarPi0);
	book(hPPi["AUAU20_30Pbar"], refnameratio20_30PbarPi0 + "_20_30Pbar", refdataratio20_30PbarPi0);
	book(hPPi["AUAU20_30PbarPi0"], refnameratio20_30PbarPi0 + "_20_30Pi0", refdataratio20_30PbarPi0);
	book(hRatio["Pbar_Pi020_30"], refnameratio20_30PbarPi0);

	string refnameratio60_92PbarPi0 = mkAxisCode(19,1,6);
	const Scatter2D& refdataratio60_92PbarPi0 = refData(refnameratio60_92PbarPi0);
	book(hPPi["AUAU60_92Pbar"], refnameratio60_92PbarPi0 + "_60_92Pbar", refdataratio60_92PbarPi0);
	book(hPPi["AUAU60_92PbarPi0"], refnameratio60_92PbarPi0 + "_60_92Pi0", refdataratio60_92PbarPi0);
	book(hRatio["Pbar_Pi060_92"], refnameratio60_92PbarPi0);

	
	//____Ratios vs Centrality____
	
	string refnamePiminusPiplus = mkAxisCode(20,1,1);	//Pi-/Pi+
	const Scatter2D& refdataPiminusPiplus = refData(refnamePiminusPiplus);
	book(hPiPi["AUAUPiminus"], refnamePiminusPiplus + "Piminus", refdataPiminusPiplus);
	book(hPiPi["AUAUPiplus"], refnamePiminusPiplus + "Piplus", refdataPiminusPiplus);
	book(hRatio["PiminusPiplus"], refnamePiminusPiplus);

	string refnameKminusKplus = mkAxisCode(20,1,2);		//K-/K+
	const Scatter2D& refdataKminusKplus = refData(refnameKminusKplus);
	book(hKK["AUAUKminus"], refnameKminusKplus + "Kminus", refdataKminusKplus);
	book(hKK["AUAUKplus"], refnameKminusKplus + "Kplus", refdataKminusKplus);
	book(hRatio["KminusKplus"], refnameKminusKplus);

	string refnamePbarP = mkAxisCode(20,1,3);
	const Scatter2D& refdataPbarP = refData(refnamePbarP);	//Antiproton/Proton
	book(hPP["AUAUPbar"], refnamePbarP + "Pbar", refdataPbarP);
	book(hPP["AUAUP"], refnamePbarP + "P", refdataPbarP);
	book(hRatio["PbarP"], refnamePbarP);

	string refnameKplusPiplus = mkAxisCode(20,1,4);		//K+/Pi+
	const Scatter2D& refdataKplusPiplus = refData(refnameKplusPiplus);
	book(hKPi["AUAUKplus"], refnameKplusPiplus + "Kplus", refdataKplusPiplus);
	book(hKPi["AUAUPiplus"], refnameKplusPiplus + "Piplus", refdataKplusPiplus);
	book(hRatio["KplusPiplus"], refnameKplusPiplus);

	string refnameKminusPiminus = mkAxisCode(20,1,5);	//K-/Pi-
	const Scatter2D& refdataKminusPiminus = refData(refnameKminusPiminus);
	book(hKPi["AUAUKminus"], refnameKminusPiminus + "Kminus", refdataKminusPiminus);
	book(hKPi["AUAUPiminus"], refnameKminusPiminus + "Piminus", refdataKminusPiminus);
	book(hRatio["KminusPiminus"], refnameKminusPiminus);

	string refnamePPiplus = mkAxisCode(20,1,6);		//Proton/Pi+
	const Scatter2D& refdataPPiplus = refData(refnamePPiplus);
	book(hPPi["AUAUP"], refnamePPiplus + "P", refdataPPiplus);
	book(hPPi["AUAUPiplus"], refnamePPiplus + "Pi", refdataPPiplus);
	book(hRatio["PPiplus"], refnamePPiplus);

	string refnamePbarPiminus = mkAxisCode(20,1,7);		//Antiproton/Pi-
	const Scatter2D& refdataPbarPiminus = refData(refnamePbarPiminus);
	book(hPPi["AUAUPbar"], refnamePbarPiminus + "Pbar", refdataPbarPiminus);
	book(hPPi["AUAUPiminus"], refnamePbarPiminus + "Pi", refdataPbarPiminus);
	book(hRatio["PbarPiminus"], refnamePbarPiminus);

	}

	void analyze(const Event& event) {

		const PrimaryParticles cp = apply<PrimaryParticles>(event, "cp");
		const UnstableParticles np = apply<UnstableParticles>(event, "np");
		const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
		const double c = cent();
		const Particles chargedParticles = cp.particles();
		const Particles neutralParticles = np.particles();
		const Particles Piplus = cp.particles(Cuts::pid == 211);
		const Particles Piminus = cp.particles(Cuts::pid == -211);
		const Particles Kplus = cp.particles(Cuts::pid == 321);
		const Particles Kminus = cp.particles(Cuts::pid == -321);
		const Particles Protons = cp.particles(Cuts::pid == 2212);
		const Particles Pbar = cp.particles(Cuts::pid == -2212);

		const ParticlePair& beam = beams();
	

		beamOpt = getOption<string>("beam","NONE");


    	if (beamOpt == "NONE") {
    	if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
          float NN = 197.;
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = AUAU200;
      }
		}
		if (collSys == AUAU200){
		pdN_dy["Piplus"]->fill(c, Piplus.size());
		pdN_dy["Piminus"]->fill(c, Piminus.size());
		pdN_dy["Kplus"]->fill(c, Kplus.size());
		pdN_dy["Kminus"]->fill(c, Kminus.size());
		pdN_dy["Protons"]->fill(c, Protons.size());
		pdN_dy["Pbar"]->fill(c, Pbar.size());

		if ((c < 0.) || (c > 92.2)) vetoEvent;

		sow["sow_AUAUall"]->fill();
			
		if ((c >= 0.) && (c < 5.)) {
				
			sow["sow_AUAU5"]->fill();
			sow["sow_AUAU0_10"]->fill();
			sow["sow_AUAUmin"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double partPtM = p.pT()/MeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus5"]->fill(partPt, pt_weight);
						pmeanPt["Piplus"]->fill(c, partPtM);
						hPi["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for Pi
						hPiPi["AUAU0_5Piplus"]->fill(partPt);	//ratio denominator for Pi-/Pi+ 0-5%
						hKPi["AUAU0_5Piplus"]->fill(partPt);	//ratio denominator for K+/Pi+ 0-5%
						hPPi["AUAU0_10Piplus"]->fill(partPt);	//ratio denominator for P/Pi+ 0-10%
						hPiPi["AUAUPiplus"]->fill(c);	//ratio denominator for Pi-/Pi+ vs c				
						hKPi["AUAUPiplus"]->fill(c);	//ratio denominator for K+/Pi+ vs c
						hPPi["AUAUPiplus"]->fill(c);	//ratio denominator for P/Pi+ vs c

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus5"]->fill(partPt, pt_weight);
						pmeanPt["Piminus"]->fill(c, partPtM);
						hPi["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for Pi
						hPiPi["AUAU0_5Piminus"]->fill(partPt);	//ratio numerator for Pi-/Pi+ 0-5%
						hKPi["AUAU0_5Piminus"]->fill(partPt);	//ratio denominator for K-/Pi- 0-5%
						hPPi["AUAU0_10Piminus"]->fill(partPt);	//ratio denominator for Pbar/Pi- 0-10%
						hPiPi["AUAUPiminus"]->fill(c);	//ratio numerator for Pi-/Pi+ vs c
						hKPi["AUAUPiminus"]->fill(c);	//ratio denominator for K-/Pi- vs c
						hPPi["AUAUPiminus"]->fill(c);	//ratio denominator for Pbar/Pi- vs c

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus5"]->fill(partPt, pt_weight);			
						pmeanPt["Kplus"]->fill(c, partPtM);
						hK["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for K
						hKK["AUAU0_5Kplus"]->fill(partPt);	//ratio denominator for K-/K+ 0-5%
						hKPi["AUAU0_5Kplus"]->fill(partPt);	//ratio numerator for K+/Pi+ 0-5%
						hKK["AUAUKplus"]->fill(c);	//ratio denominator for K-/K+ vs c
						hKPi["AUAUKplus"]->fill(c);	//ratio numerator for K+/Pi+ vs c

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus5"]->fill(partPt, pt_weight);
						pmeanPt["Kminus"]->fill(c, partPtM);
						hK["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for K
						hKK["AUAU0_5Kplus"]->fill(partPt);	//ratio numerator for K-/K+ 0-5%
						hKPi["AUAU0_5Kminus"]->fill(partPt);	//ratio numerator for K-/Pi- 0-5%
						hKK["AUAUKminus"]->fill(c);	//ratio numerator for K-/K+ vs c
						hKPi["AUAUKminus"]->fill(c);	//ratio numerator for K-/Pi- vs c

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons5"]->fill(partPt, pt_weight);			
						pmeanPt["Protons"]->fill(c, partPtM);
						hP["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for P+Pbar
						hPP["AUAU0_5P"]->fill(partPt);		//ratio denominator for Pbar/P 0-5%
						hPP["AUAUminP"]->fill(partPt);		//ratio denominator for Pbar/P minimum bias	
						hPPi["AUAU0_10ProtonPiplus"]->fill(partPt);	//ratio numerator for P/Pi+ 0-10%
						hPPi["AUAU0_10P"]->fill(partPt);	//ratio numerator for P/Pi0 0-10%
						hPP["AUAUP"]->fill(c);		//ratio denominator for Pbar/P vs c
						hPPi["AUAUP"]->fill(c);		//ratio numerator for P/Pi+ vs c

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar5"]->fill(partPt, pt_weight);	
						pmeanPt["Pbar"]->fill(c, partPtM);
						hP["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for P+Pbar
						hPP["AUAU0_5Pbar"]->fill(partPt);	//ratio numerator for Pbar/P 0-5%
						hPP["AUAUminPbar"]->fill(partPt);	//ratio numerator for Pbar/P minimum bias
						hPPi["AUAU0_10PbarPiminus"]->fill(partPt);	//ratio numerator for Pbar/Pi- 0-10%
						hPPi["AUAU0_10Pbar"]->fill(partPt);	//ratio numerator for Pbar/Pi0 0-10%
						hPP["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/P vs c
						hPPi["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/Pi- vs c

						break;
				}

				if (partPt > 1.5) {
				
					switch(p.pid()) {
			
						case 211:	//pi+
						
							hPi["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Pions above 1.5GeV
						
							break;

						case -211:	//pi-
						
							hPi["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Pions above 1.5GeV
						
							break;

						case 321:	//K+
	
							//hK["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Kaons above 1.5GeV
						
							break;

						case -321:	//K-
				
							//hK["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Kaons above 1.5GeV
						
							break;

						case 2212:	//proton

							//hP["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Protons above 1.5GeV
						
							break;

						case -2212:	//anti-proton

							//hP["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Protons above 1.5GeV
						
							break;
					}
				}	
			}
	
			for (const Particle& p : neutralParticles) {
		
				double partPt = p.pT()/GeV;

				hPPi["AUAU0_10PPi0"]->fill(partPt);	//ratio denominator for P/Pi0 0-10%
				hPPi["AUAU0_10PbarPi0"]->fill(partPt);	//ratio denominator for Pbar/Pi0 0-10%
				hPi0["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for Pi0
			}
		}

		else if ((c >= 5.) && (c < 10.)) {
				
			sow["sow_AUAU10"]->fill();
			sow["sow_AUAU0_10"]->fill();
			sow["sow_AUAUmin"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double partPtM = p.pT()/MeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus10"]->fill(partPt, pt_weight);
						pmeanPt["Piplus"]->fill(c, partPtM);
						hPi["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for Pi
						hPPi["AUAU0_10Piplus"]->fill(partPt);	//ratio denominator for P/Pi+ 0-10%
						hPiPi["AUAUPiplus"]->fill(c);	//ratio denominator for Pi-/Pi+ vs c				
						hKPi["AUAUPiplus"]->fill(c);	//ratio denominator for K+/Pi+
						hPPi["AUAUPiplus"]->fill(c);	//ratio denominator for P/Pi+ vs c

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus10"]->fill(partPt, pt_weight);
						pmeanPt["Piminus"]->fill(c, partPtM);
						hPi["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for Pi
						hPPi["AUAU0_10Piminus"]->fill(partPt);	//ratio denominator for Pbar/Pi- 0-10%
						hPiPi["AUAUPiminus"]->fill(c);	//ratio numerator for Pi-/Pi+ vs c
						hKPi["AUAUPiminus"]->fill(c);	//ratio denominator for K-/Pi- vs c
						hPPi["AUAUPiminus"]->fill(c);	//ratio denominator for Pbar/Pi- vs c

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus10"]->fill(partPt, pt_weight);			
						pmeanPt["Kplus"]->fill(c, partPtM);
						hK["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for K
						hKK["AUAUKplus"]->fill(c);	//ratio denominator for K-/K+ vs c
						hKPi["AUAUKplus"]->fill(c);	//ratio numerator for K+/Pi+ vs c

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus10"]->fill(partPt, pt_weight);
						pmeanPt["Kminus"]->fill(c, partPtM);
						hK["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for K
						hKK["AUAUKminus"]->fill(c);	//ratio numerator for K-/K+ vs c
						hKPi["AUAUKminus"]->fill(c);	//ratio numerator for K-/Pi- vs c

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons10"]->fill(partPt, pt_weight);			
						pmeanPt["Protons"]->fill(c, partPtM);
						hP["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for P+Pbar
						hPP["AUAUminP"]->fill(partPt);		//ratio denominator for Pbar/P minimum bias	
						hPPi["AUAU0_10ProtonPiplus"]->fill(partPt);	//ratio numerator for P/Pi+ 0-10%
						hPPi["AUAU0_10P"]->fill(partPt);	//ratio numerator for P/Pi0 0-10%
						hPP["AUAUP"]->fill(c);		//ratio denominator for Pbar/P vs c
						hPPi["AUAUP"]->fill(c);		//ratio numerator for P/Pi+ vs c

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar10"]->fill(partPt, pt_weight);	
						pmeanPt["Pbar"]->fill(c, partPtM);
						hP["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for P+Pbar
						hPP["AUAUminPbar"]->fill(partPt);	//ratio numerator for Pbar/P minimum bias
						hPPi["AUAU0_10PbarPiminus"]->fill(partPt);	//ratio numerator for Pbar/Pi- 0-10%
						hPPi["AUAU0_10Pbar"]->fill(partPt);	//ratio numerator for Pbar/Pi0 0-10%
						hPP["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/P vs c
						hPPi["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/Pi- vs c

						break;
				}

				if (partPt > 1.5) {
				
					switch(p.pid()) {
			
						case 211:	//pi+
						
							hPi["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Pions above 1.5GeV
						
							break;

						case -211:	//pi-
						
							hPi["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Pions above 1.5GeV
						
							break;

						case 321:	//K+
	
							//hK["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Kaons above 1.5GeV
						
							break;

						case -321:	//K-
				
							//hK["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Kaons above 1.5GeV
						
							break;

						case 2212:	//proton

							//hP["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Protons above 1.5GeV
						
							break;

						case -2212:	//anti-proton

							//hP["AUAU0_10GeV1_5"]->fill(p.pT()/GeV);	//Rcp numerator for Protons above 1.5GeV
						
							break;
					}
				}
			}

			for (const Particle& p : neutralParticles) {
		
				double partPt = p.pT()/GeV;

				hPPi["AUAU0_10PPi0"]->fill(partPt);	//ratio denominator for P/Pi0 0-10%
				hPPi["AUAU0_10PbarPi0"]->fill(partPt);	//ratio denominator for Pbar/Pi0 0-10%
				hPi0["AUAU0_10"]->fill(p.pT()/GeV);	//Rcp numerator for Pi0
			}
		}

		else if ((c >= 10.) && (c < 15.)) {
				
			sow["sow_AUAU15"]->fill();
			sow["sow_AUAUmin"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double partPtM = p.pT()/MeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus15"]->fill(partPt, pt_weight);
						pmeanPt["Piplus"]->fill(c, partPtM);
						hPiPi["AUAUPiplus"]->fill(c);	//ratio denominator for Pi-/Pi+ vs c				
						hKPi["AUAUPiplus"]->fill(c);	//ratio denominator for K+/Pi+
						hPPi["AUAUPiplus"]->fill(c);	//ratio denominator for P/Pi+ vs c

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus15"]->fill(partPt, pt_weight);
						pmeanPt["Piminus"]->fill(c, partPtM);
						hPiPi["AUAUPiminus"]->fill(c);	//ratio numerator for Pi-/Pi+ vs c
						hKPi["AUAUPiminus"]->fill(c);	//ratio denominator for K-/Pi- vs c
						hPPi["AUAUPiminus"]->fill(c);	//ratio denominator for Pbar/Pi- vs c

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus15"]->fill(partPt, pt_weight);			
						pmeanPt["Kplus"]->fill(c, partPtM);
						hKK["AUAUKplus"]->fill(c);	//ratio denominator for K-/K+ vs c
						hKPi["AUAUKplus"]->fill(c);	//ratio numerator for K+/Pi+ vs c

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus15"]->fill(partPt, pt_weight);
						pmeanPt["Kminus"]->fill(c, partPtM);
						hKK["AUAUKminus"]->fill(c);	//ratio numerator for K-/K+ vs c
						hKPi["AUAUKminus"]->fill(c);	//ratio numerator for K-/Pi- vs c

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons15"]->fill(partPt, pt_weight);			
						pmeanPt["Protons"]->fill(c, partPtM);
						hPP["AUAUminP"]->fill(partPt);		//ratio denominator for Pbar/P minimum bias	
						hPP["AUAUP"]->fill(c);		//ratio denominator for Pbar/P vs c
						hPPi["AUAUP"]->fill(c);		//ratio numerator for P/Pi+ vs c

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar15"]->fill(partPt, pt_weight);	
						pmeanPt["Pbar"]->fill(c, partPtM);
						hPP["AUAUminPbar"]->fill(partPt);	//ratio numerator for Pbar/P minimum bias
						hPP["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/P vs c
						hPPi["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/Pi- vs c

						break;
				}
			}
		}

		else if ((c >= 15.) && (c < 20.)) {
				
			sow["sow_AUAU20"]->fill();
			sow["sow_AUAUmin"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double partPtM = p.pT()/MeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus20"]->fill(partPt, pt_weight);
						pmeanPt["Piplus"]->fill(c, partPtM);
						hPiPi["AUAUPiplus"]->fill(c);	//ratio denominator for Pi-/Pi+ vs c				
						hKPi["AUAUPiplus"]->fill(c);	//ratio denominator for K+/Pi+
						hPPi["AUAUPiplus"]->fill(c);	//ratio denominator for P/Pi+ vs c

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus20"]->fill(partPt, pt_weight);
						pmeanPt["Piminus"]->fill(c, partPtM);
						hPiPi["AUAUPiminus"]->fill(c);	//ratio numerator for Pi-/Pi+ vs c
						hKPi["AUAUPiminus"]->fill(c);	//ratio denominator for K-/Pi- vs c
						hPPi["AUAUPiminus"]->fill(c);	//ratio denominator for Pbar/Pi- vs c

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus20"]->fill(partPt, pt_weight);			
						pmeanPt["Kplus"]->fill(c, partPtM);
						hKK["AUAUKplus"]->fill(c);	//ratio denominator for K-/K+ vs c
						hKPi["AUAUKplus"]->fill(c);	//ratio numerator for K+/Pi+ vs c

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus20"]->fill(partPt, pt_weight);
						pmeanPt["Kminus"]->fill(c, partPtM);
						hKK["AUAUKminus"]->fill(c);	//ratio numerator for K-/K+ vs c
						hKPi["AUAUKminus"]->fill(c);	//ratio numerator for K-/Pi- vs c

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons20"]->fill(partPt, pt_weight);			
						pmeanPt["Protons"]->fill(c, partPtM);
						hPP["AUAUminP"]->fill(partPt);		//ratio denominator for Pbar/P minimum bias	
						hPP["AUAUP"]->fill(c);		//ratio denominator for Pbar/P vs c
						hPPi["AUAUP"]->fill(c);		//ratio numerator for P/Pi+ vs c

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar20"]->fill(partPt, pt_weight);	
						pmeanPt["Pbar"]->fill(c, partPtM);
						hPP["AUAUminPbar"]->fill(partPt);	//ratio numerator for Pbar/P minimum bias
						hPP["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/P vs c
						hPPi["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/Pi- vs c

						break;
				}
			}
		}

		else if ((c >= 20.) && (c < 30.)) {
				
			sow["sow_AUAU30"]->fill();
			sow["sow_AUAUmin"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double partPtM = p.pT()/MeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus30"]->fill(partPt, pt_weight);
						pmeanPt["Piplus"]->fill(c, partPtM);
						hPPi["AUAU20_30Piplus"]->fill(partPt);	//ratio denominator for P/Pi+ 20-30%
						hPiPi["AUAUPiplus"]->fill(c);	//ratio denominator for Pi-/Pi+ vs c				
						hKPi["AUAUPiplus"]->fill(c);	//ratio denominator for K+/Pi+
						hPPi["AUAUPiplus"]->fill(c);	//ratio denominator for P/Pi+ vs c

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus30"]->fill(partPt, pt_weight);
						pmeanPt["Piminus"]->fill(c, partPtM);
						hPPi["AUAU20_30Piminus"]->fill(partPt);	//ratio denominator for Pbar/Pi- 20-30%
						hPiPi["AUAUPiminus"]->fill(c);	//ratio numerator for Pi-/Pi+ vs c
						hKPi["AUAUPiminus"]->fill(c);	//ratio denominator for K-/Pi- vs c
						hPPi["AUAUPiminus"]->fill(c);	//ratio denominator for Pbar/Pi- vs c

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus30"]->fill(partPt, pt_weight);			
						pmeanPt["Kplus"]->fill(c, partPtM);
						hKK["AUAUKplus"]->fill(c);	//ratio denominator for K-/K+ vs c
						hKPi["AUAUKplus"]->fill(c);	//ratio numerator for K+/Pi+ vs c

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus30"]->fill(partPt, pt_weight);
						pmeanPt["Kminus"]->fill(c, partPtM);
						hKK["AUAUKminus"]->fill(c);	//ratio numerator for K-/K+ vs c
						hKPi["AUAUKminus"]->fill(c);	//ratio numerator for K-/Pi- vs c

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons30"]->fill(partPt, pt_weight);			
						pmeanPt["Protons"]->fill(c, partPtM);
						hPP["AUAUminP"]->fill(partPt);		//ratio denominator for Pbar/P minimum bias	
						hPPi["AUAU20_30ProtonPiplus"]->fill(partPt);	//ratio numerator for P/Pi+ 20-30%
						hPPi["AUAU20_30P"]->fill(partPt);	//ratio numerator for P/Pi0 20-30%
						hPP["AUAUP"]->fill(c);		//ratio denominator for Pbar/P vs c
						hPPi["AUAUP"]->fill(c);		//ratio numerator for P/Pi+ vs c

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar30"]->fill(partPt, pt_weight);	
						pmeanPt["Pbar"]->fill(c, partPtM);
						hPP["AUAUminPbar"]->fill(partPt);	//ratio numerator for Pbar/P minimum bias
						hPPi["AUAU20_30PbarPiminus"]->fill(partPt);	//ratio numerator for Pbar/Pi- 20-30%
						hPPi["AUAU20_30Pbar"]->fill(partPt);	//ratio numerator for Pbar/Pi0 20-30%
						hPP["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/P vs c
						hPPi["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/Pi- vs c

						break;
				}
			}

			for (const Particle& p : neutralParticles) {
		
				double partPt = p.pT()/GeV;

				hPPi["AUAU20_30PPi0"]->fill(partPt);	//ratio denominator for P/Pi0 20-30%
				hPPi["AUAU20_30PbarPi0"]->fill(partPt);	//ratio denominator for Pbar/Pi0 20-30%
			}
		}

		else if ((c >= 30.) && (c < 40.)) {
				
			sow["sow_AUAU40"]->fill();
			sow["sow_AUAUmin"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double partPtM = p.pT()/MeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus40"]->fill(partPt, pt_weight);
						pmeanPt["Piplus"]->fill(c, partPtM);
						hPiPi["AUAUPiplus"]->fill(c);	//ratio denominator for Pi-/Pi+ vs c				
						hKPi["AUAUPiplus"]->fill(c);	//ratio denominator for K+/Pi+
						hPPi["AUAUPiplus"]->fill(c);	//ratio denominator for P/Pi+ vs c

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus40"]->fill(partPt, pt_weight);
						pmeanPt["Piminus"]->fill(c, partPtM);
						hPiPi["AUAUPiminus"]->fill(c);	//ratio numerator for Pi-/Pi+ vs c
						hKPi["AUAUPiminus"]->fill(c);	//ratio denominator for K-/Pi- vs c
						hPPi["AUAUPiminus"]->fill(c);	//ratio denominator for Pbar/Pi- vs c

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus40"]->fill(partPt, pt_weight);			
						pmeanPt["Kplus"]->fill(c, partPtM);
						hKK["AUAUKplus"]->fill(c);	//ratio denominator for K-/K+ vs c
						hKPi["AUAUKplus"]->fill(c);	//ratio numerator for K+/Pi+ vs c

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus40"]->fill(partPt, pt_weight);
						pmeanPt["Kminus"]->fill(c, partPtM);
						hKK["AUAUKminus"]->fill(c);	//ratio numerator for K-/K+ vs c
						hKPi["AUAUKminus"]->fill(c);	//ratio numerator for K-/Pi- vs c

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons40"]->fill(partPt, pt_weight);			
						pmeanPt["Protons"]->fill(c, partPtM);
						hPP["AUAUminP"]->fill(partPt);		//ratio denominator for Pbar/P minimum bias	
						hPP["AUAUP"]->fill(c);		//ratio denominator for Pbar/P vs c
						hPPi["AUAUP"]->fill(c);		//ratio numerator for P/Pi+ vs c

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar40"]->fill(partPt, pt_weight);	
						pmeanPt["Pbar"]->fill(c, partPtM);
						hPP["AUAUminPbar"]->fill(partPt);	//ratio numerator for Pbar/P minimum bias
						hPP["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/P vs c
						hPPi["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/Pi- vs c

						break;
				}
			}
		}

		else if ((c >= 40.) && (c < 50.)) {
				
			sow["sow_AUAU50"]->fill();
			sow["sow_AUAUmin"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double partPtM = p.pT()/MeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus50"]->fill(partPt, pt_weight);
						pmeanPt["Piplus"]->fill(c, partPtM);
						hPiPi["AUAUPiplus"]->fill(c);	//ratio denominator for Pi-/Pi+ vs c				
						hKPi["AUAUPiplus"]->fill(c);	//ratio denominator for K+/Pi+
						hPPi["AUAUPiplus"]->fill(c);	//ratio denominator for P/Pi+ vs c

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus50"]->fill(partPt, pt_weight);
						pmeanPt["Piminus"]->fill(c, partPtM);
						hPiPi["AUAUPiminus"]->fill(c);	//ratio numerator for Pi-/Pi+ vs c
						hKPi["AUAUPiminus"]->fill(c);	//ratio denominator for K-/Pi- vs c
						hPPi["AUAUPiminus"]->fill(c);	//ratio denominator for Pbar/Pi- vs c

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus50"]->fill(partPt, pt_weight);			
						pmeanPt["Kplus"]->fill(c, partPtM);
						hKK["AUAUKplus"]->fill(c);	//ratio denominator for K-/K+ vs c
						hKPi["AUAUKplus"]->fill(c);	//ratio numerator for K+/Pi+ vs c

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus50"]->fill(partPt, pt_weight);
						pmeanPt["Kminus"]->fill(c, partPtM);
						hKK["AUAUKminus"]->fill(c);	//ratio numerator for K-/K+ vs c
						hKPi["AUAUKminus"]->fill(c);	//ratio numerator for K-/Pi- vs c

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons50"]->fill(partPt, pt_weight);			
						pmeanPt["Protons"]->fill(c, partPtM);
						hPP["AUAUminP"]->fill(partPt);		//ratio denominator for Pbar/P minimum bias	
						hPP["AUAUP"]->fill(c);		//ratio denominator for Pbar/P vs c
						hPPi["AUAUP"]->fill(c);		//ratio numerator for P/Pi+ vs c

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar50"]->fill(partPt, pt_weight);	
						pmeanPt["Pbar"]->fill(c, partPtM);
						hPP["AUAUminPbar"]->fill(partPt);	//ratio numerator for Pbar/P minimum bias
						hPP["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/P vs c
						hPPi["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/Pi- vs c

						break;
				}
			}
		}

		else if ((c >= 50.) && (c < 60.)) {
				
			sow["sow_AUAU60"]->fill();
			sow["sow_AUAUmin"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double partPtM = p.pT()/MeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus60"]->fill(partPt, pt_weight);
						pmeanPt["Piplus"]->fill(c, partPtM);
						hPiPi["AUAUPiplus"]->fill(c);	//ratio denominator for Pi-/Pi+ vs c				
						hKPi["AUAUPiplus"]->fill(c);	//ratio denominator for K+/Pi+
						hPPi["AUAUPiplus"]->fill(c);	//ratio denominator for P/Pi+ vs c

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus60"]->fill(partPt, pt_weight);
						pmeanPt["Piminus"]->fill(c, partPtM);
						hPiPi["AUAUPiminus"]->fill(c);	//ratio numerator for Pi-/Pi+ vs c
						hKPi["AUAUPiminus"]->fill(c);	//ratio denominator for K-/Pi- vs c
						hPPi["AUAUPiminus"]->fill(c);	//ratio denominator for Pbar/Pi- vs c

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus60"]->fill(partPt, pt_weight);			
						pmeanPt["Kplus"]->fill(c, partPtM);
						hKK["AUAUKplus"]->fill(c);	//ratio denominator for K-/K+ vs c
						hKPi["AUAUKplus"]->fill(c);	//ratio numerator for K+/Pi+ vs c

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus60"]->fill(partPt, pt_weight);
						pmeanPt["Kminus"]->fill(c, partPtM);
						hKK["AUAUKminus"]->fill(c);	//ratio numerator for K-/K+ vs c
						hKPi["AUAUKminus"]->fill(c);	//ratio numerator for K-/Pi- vs c

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons60"]->fill(partPt, pt_weight);			
						pmeanPt["Protons"]->fill(c, partPtM);
						hPP["AUAUminP"]->fill(partPt);		//ratio denominator for Pbar/P minimum bias	
						hPP["AUAUP"]->fill(c);		//ratio denominator for Pbar/P vs c
						hPPi["AUAUP"]->fill(c);		//ratio numerator for P/Pi+ vs c

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar60"]->fill(partPt, pt_weight);	
						pmeanPt["Pbar"]->fill(c, partPtM);
						hPP["AUAUminPbar"]->fill(partPt);	//ratio numerator for Pbar/P minimum bias
						hPP["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/P vs c
						hPPi["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/Pi- vs c

						break;
				}
			}
		}

		else if ((c >= 60.) && (c < 70.)) {
				
			sow["sow_AUAU70"]->fill();
			sow["sow_AUAUmin"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double partPtM = p.pT()/MeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus70"]->fill(partPt, pt_weight);
						pmeanPt["Piplus"]->fill(c, partPtM);
						hPiPi["AUAUPiplus"]->fill(c);	//ratio denominator for Pi-/Pi+ vs c				
						hKPi["AUAUPiplus"]->fill(c);	//ratio denominator for K+/Pi+
						hPPi["AUAUPiplus"]->fill(c);	//ratio denominator for P/Pi+ vs c

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus70"]->fill(partPt, pt_weight);
						pmeanPt["Piminus"]->fill(c, partPtM);
						hPiPi["AUAUPiminus"]->fill(c);	//ratio numerator for Pi-/Pi+ vs c
						hKPi["AUAUPiminus"]->fill(c);	//ratio denominator for K-/Pi- vs c
						hPPi["AUAUPiminus"]->fill(c);	//ratio denominator for Pbar/Pi- vs c

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus70"]->fill(partPt, pt_weight);			
						pmeanPt["Kplus"]->fill(c, partPtM);
						hKK["AUAUKplus"]->fill(c);	//ratio denominator for K-/K+ vs c
						hKPi["AUAUKplus"]->fill(c);	//ratio numerator for K+/Pi+ vs c

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus70"]->fill(partPt, pt_weight);
						pmeanPt["Kminus"]->fill(c, partPtM);
						hKK["AUAUKminus"]->fill(c);	//ratio numerator for K-/K+ vs c
						hKPi["AUAUKminus"]->fill(c);	//ratio numerator for K-/Pi- vs c

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons70"]->fill(partPt, pt_weight);			
						pmeanPt["Protons"]->fill(c, partPtM);
						hPP["AUAUminP"]->fill(partPt);		//ratio denominator for Pbar/P minimum bias	
						hPP["AUAUP"]->fill(c);		//ratio denominator for Pbar/P vs c
						hPPi["AUAUP"]->fill(c);		//ratio numerator for P/Pi+ vs c

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar70"]->fill(partPt, pt_weight);	
						pmeanPt["Pbar"]->fill(c, partPtM);
						hPP["AUAUminPbar"]->fill(partPt);	//ratio numerator for Pbar/P minimum bias
						hPP["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/P vs c
						hPPi["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/Pi- vs c

						break;
				}
			}
		}

		else if ((c >= 70.) && (c < 80.)) {
				
			sow["sow_AUAU80"]->fill();
			sow["sow_AUAUmin"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double partPtM = p.pT()/MeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus80"]->fill(partPt, pt_weight);
						pmeanPt["Piplus"]->fill(c, partPtM);
						hPiPi["AUAUPiplus"]->fill(c);	//ratio denominator for Pi-/Pi+ vs c				
						hKPi["AUAUPiplus"]->fill(c);	//ratio denominator for K+/Pi+
						hPPi["AUAUPiplus"]->fill(c);	//ratio denominator for P/Pi+ vs c

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus80"]->fill(partPt, pt_weight);
						pmeanPt["Piminus"]->fill(c, partPtM);
						hPiPi["AUAUPiminus"]->fill(c);	//ratio numerator for Pi-/Pi+ vs c
						hKPi["AUAUPiminus"]->fill(c);	//ratio denominator for K-/Pi- vs c
						hPPi["AUAUPiminus"]->fill(c);	//ratio denominator for Pbar/Pi- vs c

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus80"]->fill(partPt, pt_weight);			
						pmeanPt["Kplus"]->fill(c, partPtM);
						hKK["AUAUKplus"]->fill(c);	//ratio denominator for K-/K+ vs c
						hKPi["AUAUKplus"]->fill(c);	//ratio numerator for K+/Pi+ vs c

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus80"]->fill(partPt, pt_weight);
						pmeanPt["Kminus"]->fill(c, partPtM);
						hKK["AUAUKminus"]->fill(c);	//ratio numerator for K-/K+ vs c
						hKPi["AUAUKminus"]->fill(c);	//ratio numerator for K-/Pi- vs c

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons80"]->fill(partPt, pt_weight);			
						pmeanPt["Protons"]->fill(c, partPtM);
						hPP["AUAUminP"]->fill(partPt);		//ratio denominator for Pbar/P minimum bias	
						hPP["AUAUP"]->fill(c);		//ratio denominator for Pbar/P vs c
						hPPi["AUAUP"]->fill(c);		//ratio numerator for P/Pi+ vs c

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar80"]->fill(partPt, pt_weight);	
						pmeanPt["Pbar"]->fill(c, partPtM);
						hPP["AUAUminPbar"]->fill(partPt);	//ratio numerator for Pbar/P minimum bias
						hPP["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/P vs c
						hPPi["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/Pi- vs c

						break;
				}
			}
		}

		if ((c >= 80.) && (c <= 92.2)) {
				
			sow["sow_AUAU92"]->fill();
			sow["sow_AUAUmin"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double partPtM = p.pT()/MeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus92"]->fill(partPt, pt_weight);
						pmeanPt["Piplus"]->fill(c, partPtM);
						hPiPi["AUAUPiplus"]->fill(c);	//ratio denominator for Pi-/Pi+ vs c				
						hKPi["AUAUPiplus"]->fill(c);	//ratio denominator for K+/Pi+
						hPPi["AUAUPiplus"]->fill(c);	//ratio denominator for P/Pi+ vs c

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus92"]->fill(partPt, pt_weight);
						pmeanPt["Piminus"]->fill(c, partPtM);
						hPiPi["AUAUPiminus"]->fill(c);	//ratio numerator for Pi-/Pi+ vs c
						hKPi["AUAUPiminus"]->fill(c);	//ratio denominator for K-/Pi- vs c
						hPPi["AUAUPiminus"]->fill(c);	//ratio denominator for Pbar/Pi- vs c

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus92"]->fill(partPt, pt_weight);			
						pmeanPt["Kplus"]->fill(c, partPtM);
						hKK["AUAUKplus"]->fill(c);	//ratio denominator for K-/K+ vs c
						hKPi["AUAUKplus"]->fill(c);	//ratio numerator for K+/Pi+ vs c

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus92"]->fill(partPt, pt_weight);
						pmeanPt["Kminus"]->fill(c, partPtM);
						hKK["AUAUKminus"]->fill(c);	//ratio numerator for K-/K+ vs c
						hKPi["AUAUKminus"]->fill(c);	//ratio numerator for K-/Pi- vs c

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons92"]->fill(partPt, pt_weight);			
						pmeanPt["Protons"]->fill(c, partPtM);
						hPP["AUAUminP"]->fill(partPt);		//ratio denominator for Pbar/P minimum bias	
						hPP["AUAUP"]->fill(c);		//ratio denominator for Pbar/P vs c
						hPPi["AUAUP"]->fill(c);		//ratio numerator for P/Pi+ vs c

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar92"]->fill(partPt, pt_weight);	
						pmeanPt["Pbar"]->fill(c, partPtM);
						hPP["AUAUminP"]->fill(partPt);	//ratio numerator for Pbar/P minimum bias
						hPP["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/P vs c
						hPPi["AUAUPbar"]->fill(c);	//ratio numerator for Pbar/Pi- vs c

						break;
				}
			}
		}

		if ((c >= 60.) && (c <= 92.2)) {
				
			sow["sow_AUAU_60_92"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid()) {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplus60_92"]->fill(partPt, pt_weight);
						hPi["AUAU60_92"]->fill(p.pT()/GeV);	//Rcp denominator for Pi
						hPiPi["AUAU60_92Piplus"]->fill(partPt);	//ratio denominator for Pi-/Pi+ 60-92%
						hKPi["AUAU60_92Piplus"]->fill(partPt);	//ratio denominator for K+/Pi+ 60-92%
						hPPi["AUAU60_92Piplus"]->fill(partPt);	//ratio denominator for P/Pi+ 60-92%

						break;
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminus60_92"]->fill(partPt, pt_weight);
						hPi["AUAU60_92"]->fill(p.pT()/GeV);	//Rcp denominator for Pi
						hPiPi["AUAU60_92Piminus"]->fill(partPt);	//ratio numerator for Pi-/Pi+ 60-92%
						hKPi["AUAU60_92Piminus"]->fill(partPt);	//ratio denominator for K-/Pi- 60-92%
						hPPi["AUAU60_92Piminus"]->fill(partPt);	//ratio denominator for Pbar/Pi- 60-92%

						break;
					
					case 321:	//K+
						
						hAUAU_Yields["Kplus60_92"]->fill(partPt, pt_weight);			
						hK["AUAU60_92"]->fill(p.pT()/GeV);	//Rcp denominator for K
						hKK["AUAU60_92Kplus"]->fill(partPt);	//ratio denominator for K-/K+ 60-92%
						hKPi["AUAU60_92Kplus"]->fill(partPt);	//ratio numerator for K+/Pi+ 60-92%

						break;

					case -321:	//K-
						
						hAUAU_Yields["Kminus60_92"]->fill(partPt, pt_weight);
						hK["AUAU60_92"]->fill(p.pT()/GeV);	//Rcp denominator for K
						hKK["AUAU60_92Kplus"]->fill(partPt);	//ratio numerator for K-/K+ 60-92%
						hKPi["AUAU60_92Kminus"]->fill(partPt);	//ratio numerator for K-/Pi- 60-92%

						break;
					
					case 2212:	//proton
						
						hAUAU_Yields["Protons60_92"]->fill(partPt, pt_weight);			
						hP["AUAU60_92"]->fill(p.pT()/GeV);	//Rcp denominator for P+Pbar
						hPP["AUAU60_92P"]->fill(partPt);	//ratio denominator for Pbar/P 60-92%	
						hPPi["AUAU60_92ProtonPiplus"]->fill(partPt);	//ratio numerator for P/Pi+ 60-92%
						hPPi["AUAU60_92P"]->fill(partPt);	//ratio numerator for P/Pi0 60-92%

						break;

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbar60_92"]->fill(partPt, pt_weight);	
						hP["AUAU60_92"]->fill(p.pT()/GeV);	//Rcp denominator for P+Pbar
						hPP["AUAU60_92Pbar"]->fill(partPt);	//ratio numerator for Pbar/P 60-92%
						hPPi["AUAU60_92PbarPiminus"]->fill(partPt);	//ratio numerator for Pbar/Pi- 60-92%
						hPPi["AUAU60_92Pbar"]->fill(partPt);	//ratio numerator for Pbar/Pi0 60-92%

						break;
				}

				if (partPt > 1.5) {
				
					switch(p.pid()) {
			
						case 211:	//pi+
						
							hPi["AUAU60_92GeV1_5"]->fill(p.pT()/GeV);	//Rcp denominator for Pions above 1.5GeV
						
							break;

						case -211:	//pi-
						
							hPi["AUAU60_92GeV1_5"]->fill(p.pT()/GeV);	//Rcp denominator for Pions above 1.5GeV
						
							break;

						case 321:	//K+
	
							//hK["AUAU60_92GeV1_5"]->fill(p.pT()/GeV);	//Rcp denominator for Kaons above 1.5GeV
						
							break;

						case -321:	//K-
				
							//hK["AUAU60_92GeV1_5"]->fill(p.pT()/GeV);	//Rcp denominator for Kaons above 1.5GeV
						
							break;

						case 2212:	//proton

							//hP["AUAU60_92GeV1_5"]->fill(p.pT()/GeV);	//Rcp denominator for Protons above 1.5GeV
						
							break;

						case -2212:	//anti-proton

							//hP["AUAU60_92GeV1_5"]->fill(p.pT()/GeV);	//Rcp denominator for Protons above 1.5GeV
						
							break;
					}
				}
			}

			for (const Particle& p : neutralParticles) {
		
				double partPt = p.pT()/GeV;

				hPPi["AUAU60_92PPi0"]->fill(partPt);	//ratio denominator for P/Pi0 60-92%
				hPPi["AUAU60_92PbarPi0"]->fill(partPt);	//ratio denominator for Pbar/Pi0 60-92%
				hPi0["AUAU60_92"]->fill(p.pT()/GeV);	//Rcp denominator for Pi0
			}
		}
	}
	}


	void finalize() {
		binShift(*hAUAU_Yields["Piplusmin"]);
		binShift(*hAUAU_Yields["Piplus5"]);
		binShift(*hAUAU_Yields["Piplus10"]);
		binShift(*hAUAU_Yields["Piplus15"]);
		binShift(*hAUAU_Yields["Piplus20"]);
		binShift(*hAUAU_Yields["Piplus30"]);
		binShift(*hAUAU_Yields["Piplus40"]);
		binShift(*hAUAU_Yields["Piplus50"]);
		binShift(*hAUAU_Yields["Piplus60"]);
		binShift(*hAUAU_Yields["Piplus70"]);
		binShift(*hAUAU_Yields["Piplus80"]);
		binShift(*hAUAU_Yields["Piplus92"]);
		binShift(*hAUAU_Yields["Piplus60_92"]);

		binShift(*hAUAU_Yields["Piminusmin"]);
		binShift(*hAUAU_Yields["Piminus5"]);
		binShift(*hAUAU_Yields["Piminus10"]);
		binShift(*hAUAU_Yields["Piminus15"]);
		binShift(*hAUAU_Yields["Piminus20"]);
		binShift(*hAUAU_Yields["Piminus30"]);
		binShift(*hAUAU_Yields["Piminus40"]);
		binShift(*hAUAU_Yields["Piminus50"]);
		binShift(*hAUAU_Yields["Piminus60"]);
		binShift(*hAUAU_Yields["Piminus70"]);
		binShift(*hAUAU_Yields["Piminus80"]);
		binShift(*hAUAU_Yields["Piminus92"]);
		binShift(*hAUAU_Yields["Piminus60_92"]);

		binShift(*hAUAU_Yields["Kplusmin"]);
		binShift(*hAUAU_Yields["Kplus5"]);
		binShift(*hAUAU_Yields["Kplus10"]);
		binShift(*hAUAU_Yields["Kplus15"]);
		binShift(*hAUAU_Yields["Kplus20"]);
		binShift(*hAUAU_Yields["Kplus30"]);
		binShift(*hAUAU_Yields["Kplus40"]);
		binShift(*hAUAU_Yields["Kplus50"]);
		binShift(*hAUAU_Yields["Kplus60"]);
		binShift(*hAUAU_Yields["Kplus70"]);
		binShift(*hAUAU_Yields["Kplus80"]);
		binShift(*hAUAU_Yields["Kplus92"]);
		binShift(*hAUAU_Yields["Kplus60_92"]);

		binShift(*hAUAU_Yields["Kminusmin"]);
		binShift(*hAUAU_Yields["Kminus5"]);
		binShift(*hAUAU_Yields["Kminus10"]);
		binShift(*hAUAU_Yields["Kminus15"]);
		binShift(*hAUAU_Yields["Kminus20"]);
		binShift(*hAUAU_Yields["Kminus30"]);
		binShift(*hAUAU_Yields["Kminus40"]);
		binShift(*hAUAU_Yields["Kminus50"]);
		binShift(*hAUAU_Yields["Kminus60"]);
		binShift(*hAUAU_Yields["Kminus70"]);
		binShift(*hAUAU_Yields["Kminus80"]);
		binShift(*hAUAU_Yields["Kminus92"]);
		binShift(*hAUAU_Yields["Kminus60_92"]);

		binShift(*hAUAU_Yields["Protonsmin"]);
		binShift(*hAUAU_Yields["Protons5"]);
		binShift(*hAUAU_Yields["Protons10"]);
		binShift(*hAUAU_Yields["Protons15"]);
		binShift(*hAUAU_Yields["Protons20"]);
		binShift(*hAUAU_Yields["Protons30"]);
		binShift(*hAUAU_Yields["Protons40"]);
		binShift(*hAUAU_Yields["Protons50"]);
		binShift(*hAUAU_Yields["Protons60"]);
		binShift(*hAUAU_Yields["Protons70"]);
		binShift(*hAUAU_Yields["Protons80"]);
		binShift(*hAUAU_Yields["Protons92"]);
		binShift(*hAUAU_Yields["Protons60_92"]);

		binShift(*hAUAU_Yields["Pbarmin"]);
		binShift(*hAUAU_Yields["Pbar5"]);
		binShift(*hAUAU_Yields["Pbar10"]);
		binShift(*hAUAU_Yields["Pbar15"]);
		binShift(*hAUAU_Yields["Pbar20"]);
		binShift(*hAUAU_Yields["Pbar30"]);
		binShift(*hAUAU_Yields["Pbar40"]);
		binShift(*hAUAU_Yields["Pbar50"]);
		binShift(*hAUAU_Yields["Pbar60"]);
		binShift(*hAUAU_Yields["Pbar70"]);
		binShift(*hAUAU_Yields["Pbar80"]);
		binShift(*hAUAU_Yields["Pbar92"]);
		binShift(*hAUAU_Yields["Pbar60_92"]);

		binShift(*hPiPi["AUAU0_5Piminus"]); 
		binShift(*hPiPi["AUAU0_5Piplus"]); 

		/*binShift(*hK["AUAU0_10"]);
		binShift(*hK["AUAU0_92"]);*/

		/*binShift(*hP["AUAU0_10"]);
		binShift(*hP["AUAU60_92"]);*/

		/*binShift(*hPi0["AUAU0_10"]);
		binShift(*hPi0["AUAU60_92"]);

		binShift(*hPiPi["AUAU0_5Piminus"]);
		binShift(*hPiPi["AUAU0_5Piplus"]);

		binShift(*hPiPi["AUAU0_92Piminus"]);
		binShift(*hPiPi["AUAU0_92Piplus"]);*/

		binShift(*hPiPi["AUAU0_5Kminus"]);
		binShift(*hPiPi["AUAU0_5kplus"]);
		
		


		//____Yields____
		hAUAU_Yields["Piplusmin"]->scaleW(1./sow["sow_AUAUmin"]->sumW());	//minimum bias centrality
		hAUAU_Yields["Piminusmin"]->scaleW(1./sow["sow_AUAUmin"]->sumW());
		hAUAU_Yields["Kplusmin"]->scaleW(1./sow["sow_AUAUmin"]->sumW());
		hAUAU_Yields["Kminusmin"]->scaleW(1./sow["sow_AUAUmin"]->sumW());
		hAUAU_Yields["Protonsmin"]->scaleW(1./sow["sow_AUAUmin"]->sumW());
		hAUAU_Yields["Pbarmin"]->scaleW(1./sow["sow_AUAUmin"]->sumW());

		hAUAU_Yields["Piplus5"]->scaleW(1./sow["sow_AUAU5"]->sumW());	//0-5% centrality
		hAUAU_Yields["Piminus5"]->scaleW(1./sow["sow_AUAU5"]->sumW());
		hAUAU_Yields["Kplus5"]->scaleW(1./sow["sow_AUAU5"]->sumW());
		hAUAU_Yields["Kminus5"]->scaleW(1./sow["sow_AUAU5"]->sumW());
		hAUAU_Yields["Protons5"]->scaleW(1./sow["sow_AUAU5"]->sumW());
		hAUAU_Yields["Pbar5"]->scaleW(1./sow["sow_AUAU5"]->sumW());

		hAUAU_Yields["Piplus10"]->scaleW(1./sow["sow_AUAU10"]->sumW());	//5-10% centrality
		hAUAU_Yields["Piminus10"]->scaleW(1./sow["sow_AUAU10"]->sumW());
		hAUAU_Yields["Kplus10"]->scaleW(1./sow["sow_AUAU10"]->sumW());
		hAUAU_Yields["Kminus10"]->scaleW(1./sow["sow_AUAU10"]->sumW());
		hAUAU_Yields["Protons10"]->scaleW(1./sow["sow_AUAU10"]->sumW());
		hAUAU_Yields["Pbar10"]->scaleW(1./sow["sow_AUAU10"]->sumW());

		hAUAU_Yields["Piplus15"]->scaleW(1./sow["sow_AUAU15"]->sumW());	//10-15% centrality
		hAUAU_Yields["Piminus15"]->scaleW(1./sow["sow_AUAU15"]->sumW());
		hAUAU_Yields["Kplus15"]->scaleW(1./sow["sow_AUAU15"]->sumW());
		hAUAU_Yields["Kminus15"]->scaleW(1./sow["sow_AUAU15"]->sumW());
		hAUAU_Yields["Protons15"]->scaleW(1./sow["sow_AUAU15"]->sumW());
		hAUAU_Yields["Pbar15"]->scaleW(1./sow["sow_AUAU15"]->sumW());

		hAUAU_Yields["Piplus20"]->scaleW(1./sow["sow_AUAU20"]->sumW());	//15-20% centrality
		hAUAU_Yields["Piminus20"]->scaleW(1./sow["sow_AUAU20"]->sumW());
		hAUAU_Yields["Kplus20"]->scaleW(1./sow["sow_AUAU20"]->sumW());
		hAUAU_Yields["Kminus20"]->scaleW(1./sow["sow_AUAU20"]->sumW());
		hAUAU_Yields["Protons20"]->scaleW(1./sow["sow_AUAU20"]->sumW());
		hAUAU_Yields["Pbar20"]->scaleW(1./sow["sow_AUAU20"]->sumW());

		hAUAU_Yields["Piplus30"]->scaleW(1./sow["sow_AUAU30"]->sumW());	//20-30% centrality
		hAUAU_Yields["Piminus30"]->scaleW(1./sow["sow_AUAU30"]->sumW());
		hAUAU_Yields["Kplus30"]->scaleW(1./sow["sow_AUAU30"]->sumW());
		hAUAU_Yields["Kminus30"]->scaleW(1./sow["sow_AUAU30"]->sumW());
		hAUAU_Yields["Protons30"]->scaleW(1./sow["sow_AUAU30"]->sumW());
		hAUAU_Yields["Pbar30"]->scaleW(1./sow["sow_AUAU30"]->sumW());

		hAUAU_Yields["Piplus40"]->scaleW(1./sow["sow_AUAU40"]->sumW());	//30-40% centrality
		hAUAU_Yields["Piminus40"]->scaleW(1./sow["sow_AUAU40"]->sumW());
		hAUAU_Yields["Kplus40"]->scaleW(1./sow["sow_AUAU40"]->sumW());
		hAUAU_Yields["Kminus40"]->scaleW(1./sow["sow_AUAU40"]->sumW());
		hAUAU_Yields["Protons40"]->scaleW(1./sow["sow_AUAU40"]->sumW());
		hAUAU_Yields["Pbar40"]->scaleW(1./sow["sow_AUAU40"]->sumW());

		hAUAU_Yields["Piplus50"]->scaleW(1./sow["sow_AUAU50"]->sumW());	//40-50% centrality
		hAUAU_Yields["Piminus50"]->scaleW(1./sow["sow_AUAU50"]->sumW());
		hAUAU_Yields["Kplus50"]->scaleW(1./sow["sow_AUAU50"]->sumW());
		hAUAU_Yields["Kminus50"]->scaleW(1./sow["sow_AUAU50"]->sumW());
		hAUAU_Yields["Protons50"]->scaleW(1./sow["sow_AUAU50"]->sumW());
		hAUAU_Yields["Pbar50"]->scaleW(1./sow["sow_AUAU50"]->sumW());

		hAUAU_Yields["Piplus60"]->scaleW(1./sow["sow_AUAU60"]->sumW());	//50-60% centrality
		hAUAU_Yields["Piminus60"]->scaleW(1./sow["sow_AUAU60"]->sumW());
		hAUAU_Yields["Kplus60"]->scaleW(1./sow["sow_AUAU60"]->sumW());
		hAUAU_Yields["Kminus60"]->scaleW(1./sow["sow_AUAU60"]->sumW());
		hAUAU_Yields["Protons60"]->scaleW(1./sow["sow_AUAU60"]->sumW());
		hAUAU_Yields["Pbar60"]->scaleW(1./sow["sow_AUAU60"]->sumW());

		hAUAU_Yields["Piplus70"]->scaleW(1./sow["sow_AUAU70"]->sumW());	//60-70% centrality
		hAUAU_Yields["Piminus70"]->scaleW(1./sow["sow_AUAU70"]->sumW());
		hAUAU_Yields["Kplus70"]->scaleW(1./sow["sow_AUAU70"]->sumW());
		hAUAU_Yields["Kminus70"]->scaleW(1./sow["sow_AUAU70"]->sumW());
		hAUAU_Yields["Protons70"]->scaleW(1./sow["sow_AUAU70"]->sumW());
		hAUAU_Yields["Pbar70"]->scaleW(1./sow["sow_AUAU70"]->sumW());

		hAUAU_Yields["Piplus80"]->scaleW(1./sow["sow_AUAU80"]->sumW());	//70-80% centrality
		hAUAU_Yields["Piminus80"]->scaleW(1./sow["sow_AUAU80"]->sumW());
		hAUAU_Yields["Kplus80"]->scaleW(1./sow["sow_AUAU80"]->sumW());
		hAUAU_Yields["Kminus80"]->scaleW(1./sow["sow_AUAU80"]->sumW());
		hAUAU_Yields["Protons80"]->scaleW(1./sow["sow_AUAU80"]->sumW());
		hAUAU_Yields["Pbar80"]->scaleW(1./sow["sow_AUAU80"]->sumW());

		hAUAU_Yields["Piplus92"]->scaleW(1./sow["sow_AUAU92"]->sumW());	//80-92% centrality
		hAUAU_Yields["Piminus92"]->scaleW(1./sow["sow_AUAU92"]->sumW());
		hAUAU_Yields["Kplus92"]->scaleW(1./sow["sow_AUAU92"]->sumW());
		hAUAU_Yields["Kminus92"]->scaleW(1./sow["sow_AUAU92"]->sumW());
		hAUAU_Yields["Protons92"]->scaleW(1./sow["sow_AUAU92"]->sumW());
		hAUAU_Yields["Pbar92"]->scaleW(1./sow["sow_AUAU92"]->sumW());

		hAUAU_Yields["Piplus60_92"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());	//60-92% centrality
		hAUAU_Yields["Piminus60_92"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());
		hAUAU_Yields["Kplus60_92"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());
		hAUAU_Yields["Kminus60_92"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());
		hAUAU_Yields["Protons60_92"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());
		hAUAU_Yields["Pbar60_92"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());


		//____Rcp____

		hPi["AUAU0_10"]->scaleW(1./sow["sow_AUAU0_10"]->sumW());
		hPi["AUAU60_92"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());
		divide(hPi["AUAU0_10"], hPi["AUAU60_92"], hRcp["Pions"]);	//Rcp Pi
		
		hK["AUAU0_10"]->scaleW(1./sow["sow_AUAU0_10"]->sumW());
		hK["AUAU60_92"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());
		divide(hK["AUAU0_10"], hK["AUAU60_92"], hRcp["Kaons"]);		//Rcp K
		
		hP["AUAU0_10"]->scaleW(1./sow["sow_AUAU0_10"]->sumW());
		hP["AUAU60_92"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());
		divide(hP["AUAU0_10"], hP["AUAU60_92"], hRcp["Pbar+P"]);	//Rcp Pbar+P

		hPi0["AUAU0_10"]->scaleW(1./sow["sow_AUAU0_10"]->sumW());
		hPi0["AUAU60_92"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());
		divide(hPi0["AUAU0_10"], hPi0["AUAU60_92"], hRcp["Pi0"]);	//Rcp Pi0		


		//____Rcp 1.5GeV+____

		divide(hPi["AUAU0_10GeV1_5"], hPi["AUAU60_92GeV1_5"], hRcp["Pions1_5"]);	//Rcp Pi 1.5GeV

		//divide(hK["AUAU0_10GeV1_5"], hK["AUAU60_92GeV1_5"], hRcp["Kaons1_5"]);	//Rcp K 1.5GeV

		//divide(hP["AUAU0_10GeV1_5"], hP["AUAU60_92GeV1_5"], hRcp["Pbar+P1_5"]);	//RcP Proton 1.5GeV


		//____Particle Ratios____

//		hPiPi["AUAU0_5Piminus"]->scaleW(1./sow["sow_AUAU5"]->sumW());
//		hPiPi["AUAU0_5Piplus"]->scaleW(1./sow["sow_AUAU5"]->sumW());
		divide(hPiPi["AUAU0_5Piminus"], hPiPi["AUAU0_5Piplus"], hRatio["Piminus_Piplus0_5"]);	//Pi-/Pi+

//		hPiPi["AUAU60_92Piminus"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());
//		hPiPi["AUAU60_92Piplus"]->scaleW(1./sow["sow_AUAU_60_92"]->sumW());
		divide(hPiPi["AUAU60_92Piminus"], hPiPi["AUAU60_92Piplus"], hRatio["Piminus_Piplus60_92"]);
		
//		hPiPi["AUAU0_5Kminus"]->scaleW(1./sow["sow_AUAU5"]->sumW());
//		hPiPi["AUAU0_5Kplus"]->scaleW(1./sow["sow_AUAU5"]->sumW());
		divide(hKK["AUAU0_5Kminus"], hKK["AUAU0_5Kplus"], hRatio["Kiminus_Kiplus0_5"]);	//K-/K+

		divide(hKK["AUAU60_92Kminus"], hKK["AUAU60_92Kplus"], hRatio["Kiminus_Kiplus60_92"]);

		divide(hPP["AUAU0_5Pbar"], hPP["AUAU0_5P"], hRatio["Pbar_P0_5"]);	//Pbar/P

		divide(hPP["AUAU60_92Pbar"], hPP["AUAU60_92P"], hRatio["Pbar_P60_92"]);

		divide(hPP["AUAUminPbar"], hPP["AUAUminP"], hRatio["Pbar_Pmin"]);
	
		divide(hKPi["AUAU0_5Kplus"], hKPi["AUAU0_5Piplus"], hRatio["Kplus_Piplus0_5"]);	//K+/Pi+

		divide(hKPi["AUAU60_92Kplus"], hKPi["AUAU60_92Piplus"], hRatio["Kplus_Piplus60_92"]);

		divide(hKPi["AUAU0_5Kminus"], hKPi["AUAU0_5Piminus"], hRatio["Kminus_Piminus0_5"]);	//K-/Pi-

		divide(hKPi["AUAU60_92Kminus"], hKPi["AUAU60_92Piminus"], hRatio["Kminus_Piminus60_92"]);

		divide(hPPi["AUAU0_10ProtonPiplus"], hPPi["AUAU0_10Piplus"], hRatio["P_Piplus0_10"]);	//P/Pi+

		divide(hPPi["AUAU20_30ProtonPiplus"], hPPi["AUAU20_30Piplus"], hRatio["P_Piplus20_30"]);

		divide(hPPi["AUAU60_92ProtonPiplus"], hPPi["AUAU60_92Piplus"], hRatio["P_Piplus60_92"]);

		divide(hPPi["AUAU0_10PbarPiminus"], hPPi["AUAU0_10Piminus"], hRatio["Pbar_Piminus0_10"]);	//Pbar/Pi-

		divide(hPPi["AUAU20_30PbarPiminus"], hPPi["AUAU20_30Piminus"], hRatio["Pbar_Piminus20_30"]);

		divide(hPPi["AUAU60_92PbarPiminus"], hPPi["AUAU60_92Piminus"], hRatio["Pbar_Piminus60_92"]);

		divide(hPPi["AUAU0_10P"], hPPi["AUAU0_10PPi0"], hRatio["P_Pi00_10"]);	//P/Pi0

		divide(hPPi["AUAU20_30P"], hPPi["AUAU20_30PPi0"], hRatio["P_Pi020_30"]);

		divide(hPPi["AUAU60_92P"], hPPi["AUAU60_92PPi0"], hRatio["P_Pi060_92"]);

		divide(hPPi["AUAU0_10Pbar"], hPPi["AUAU0_10PbarPi0"], hRatio["Pbar_Pi00_10"]);	//Pbar/Pi0

		divide(hPPi["AUAU20_30Pbar"], hPPi["AUAU20_30PbarPi0"], hRatio["Pbar_Pi020_30"]);

		divide(hPPi["AUAU60_92Pbar"], hPPi["AUAU60_92PbarPi0"], hRatio["Pbar_Pi060_92"]);


		//____Ratios vs Centrality____

		divide(hPiPi["AUAUPiminus"], hPiPi["AUAUPiplus"], hRatio["PiminusPiplus"]);	//Pi-/Pi+
		divide(hKK["AUAUKminus"], hKK["AUAUKplus"], hRatio["KminusKplus"]);	//K-/K+
		divide(hPP["AUAUPbar"], hPP["AUAUP"], hRatio["PbarP"]);	//Pbar/P
		divide(hKPi["AUAUKplus"], hKPi["AUAUPiplus"], hRatio["KplusPiplus"]);	//K+/Pi+
		divide(hKPi["AUAUKminus"], hKPi["AUAUPiminus"], hRatio["KminusPiminus"]);	//K-/Pi-
		divide(hPPi["AUAUP"], hPPi["AUAUPiplus"], hRatio["PPiplus"]);	//P/Pi+
		divide(hPPi["AUAUPbar"], hPPi["AUAUPiminus"], hRatio["PbarPiminus"]);	//Pbar/Pi-

	}

	map<string, Histo1DPtr> hAUAU_Yields;
	map<string, CounterPtr> sow;
	map<string, Profile1DPtr> pmeanPt;
	map<string, Profile1DPtr> pdN_dy;
	map<string, Histo1DPtr> hPi;
	map<string, Histo1DPtr> hK;
	map<string, Histo1DPtr> hP;
	map<string, Scatter2DPtr> hRcp;
	map<string, Histo1DPtr> hPi0;
	map<string, Histo1DPtr> hPiPi;
	map<string, Scatter2DPtr> hRatio;
	map<string, Histo1DPtr> hKK;
	map<string, Histo1DPtr> hPP;
	map<string, Histo1DPtr> hKPi;
	map<string, Histo1DPtr> hPPi;

	string beamOpt;
    enum CollisionSystem {AUAU200};
    CollisionSystem collSys;
	};

	DECLARE_RIVET_PLUGIN(PHENIX_2004_I624474);
}
