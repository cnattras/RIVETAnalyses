// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "../Centralities/RHICCentrality.hh"
#include <math.h>
#include <iostream>
#define _USE_MATH_DEFINES

namespace Rivet {


  /// @brief Add a short analysis description here
	class PHENIX_2010_I856259 : public Analysis {
	public:

		/// Constructor
		RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2010_I856259);


		/// Book histograms and initialise projections before the run
		void init() {
		
			//Particles: eta meson, Pi0?
			
			std::initializer_list<int> pdgIds = {221, 111};

			const UnstableParticles np(pdgIds, Cuts::abseta < .35);
			declare(np, "np");

			declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
		
			//____Cross Section vs. pT____

			book(AUAU_cross["cross0_5"],1,1,1);
			book(AUAU_cross["cross0_10"],2,1,1);
			book(AUAU_cross["cross10_20"],3,1,1);
			book(AUAU_cross["cross0_20"],4,1,1);
			book(AUAU_cross["cross20_40"],5,1,1);
			book(AUAU_cross["cross40_60"],6,1,1);
			book(AUAU_cross["cross20_60"],7,1,1);
			book(AUAU_cross["cross60_92"],8,1,1);
			book(AUAU_cross["cross0_92"],9,1,1);
			book(AUAU_cross["crossPP"],10,1,1);


			//____Eta Raa____

			string refnameRaaEta = mkAxisCode(11,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield"], refnameRaaEta + "_0_5Eta", refdataRaaEta);
			book(hEta["AUAUcross"], refnameRaaEta + "_0_5Eta", refdataRaaEta);
			book(sRaa{"RaaEta0_5"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(12,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield"], refnameRaaEta + "_0_10Eta", refdataRaaEta);
			book(hEta["AUAUcross"], refnameRaaEta + "_0_10Eta", refdataRaaEta);
			book(sRaa{"RaaEta0_10"], refnameRaaEta);
	
			string refnameRaaEta = mkAxisCode(13,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield"], refnameRaaEta + "_10_20Eta", refdataRaaEta);
			book(hEta["AUAUcross"], refnameRaaEta + "_10_20Eta", refdataRaaEta);
			book(sRaa{"RaaEta10_20"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(14,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield"], refnameRaaEta + "_0_20Eta", refdataRaaEta);
			book(hEta["AUAUcross"], refnameRaaEta + "_0_20Eta", refdataRaaEta);
			book(sRaa{"RaaEta0_20"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(15,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield"], refnameRaaEta + "_20_40Eta", refdataRaaEta);
			book(hEta["AUAUcross"], refnameRaaEta + "_20_40Eta", refdataRaaEta);
			book(sRaa{"RaaEta20_40"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(16,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield"], refnameRaaEta + "_40_60Eta", refdataRaaEta);
			book(hEta["AUAUcross"], refnameRaaEta + "_40_60Eta", refdataRaaEta);
			book(sRaa{"RaaEta40_60"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(17,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield"], refnameRaaEta + "_20_60Eta", refdataRaaEta);
			book(hEta["AUAUcross"], refnameRaaEta + "_20_60Eta", refdataRaaEta);
			book(sRaa{"RaaEta20_60"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(18,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield"], refnameRaaEta + "_60_92Eta", refdataRaaEta);
			book(hEta["AUAUcross"], refnameRaaEta + "_60_92Eta", refdataRaaEta);
			book(sRaa{"RaaEta60_92"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(19,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield"], refnameRaaEta + "_0_92Eta", refdataRaaEta);
			book(hEta["AUAUcross"], refnameRaaEta + "_0_92Eta", refdataRaaEta);
			book(sRaa{"RaaEta0_92"], refnameRaaEta);

			string refnameRaaPi = mkAxisCode(20,1,1);
			const Scatter2D& refdataRaaPi = refData(refnameRaaPi);
			book(hPi["AUAUyield"], refnameRaaPi + "_0_92Pi", refdataRaaPi);
			book(hPi["AUAUcross"], refnameRaaPi + "_0_92Pi", refdataRaaPi);
			book(sRaa{"RaaPi0_92"], refnameRaaPi);





		}


		/// Perform the per-event analysis
		void analyze(const Event& event) {

			const UnstableParticles np = apply<UnstableParticles>(event, "np");
			const double c = cent();
			const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
			const Particles unstableParticles = np.particles();

			
		}


		/// Normalise histograms etc., after the run
		void finalize() {
		
		}

		map<string, Histo1DPtr> AUAU_cross;
		map<string, CounterPtr> sow;
		map<string, Scatter2DPtr> sRaa;
		map<string, Histo1DPtr> hPi;
		map<string, Histo1DPtr> hEta;
	};


	RIVET_DECLARE_PLUGIN(PHENIX_2010_I856259);
}
