// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "../Centralities/RHICCentrality.hh"
#include <math.h>
#include "Rivet/Tools/Cuts.hh"
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

			const UnstableParticles up(pdgIds, Cuts::abseta < .35);
			declare(up, "up");

			declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
		
			//____Counters____
			
			book(sow["sow_AUAU0_5"],"sow_AUAU0_5");
			book(sow["sow_AUAU0_10"],"sow_AUAU0_10");
			book(sow["sow_AUAU10_20"],"sow_AUAU10_20");
			book(sow["sow_AUAU0_20"],"sow_AUAU0_20");
			book(sow["sow_AUAU20_40"],"sow_AUAU20_40");
			book(sow["sow_AUAU40_60"],"sow_AUAU40_60");
			book(sow["sow_AUAU20_60"],"sow_AUAU20_60");
			book(sow["sow_AUAU60_92"],"sow_AUAU60_92");
			book(sow["sow_AUAU0_92"],"sow_AUAU0_92");
			book(sow["sow_PP"],"sow_PP");


			//____Yields vs. pT____

			book(AUAU_yield["yield0_5"],1,1,1);
			book(AUAU_yield["yield0_10"],2,1,1);
			book(AUAU_yield["yield10_20"],3,1,1);
			book(AUAU_yield["yield0_20"],4,1,1);
			book(AUAU_yield["yield20_40"],5,1,1);
			book(AUAU_yield["yield40_60"],6,1,1);
			book(AUAU_yield["yield20_60"],7,1,1);
			book(AUAU_yield["yield60_92"],8,1,1);
			book(AUAU_yield["yield0_92"],9,1,1);
			book(PP_yield["yieldPP"],10,1,1);


			//____Eta Raa____

			string refnameRaaEta = mkAxisCode(11,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield0_5"], refnameRaaEta + "_0_5Eta", refdataRaaEta);
			book(hEta["PPcross0_5"], refnameRaaEta + "_0_5Eta", refdataRaaEta);
			book(sRaa["RaaEta0_5"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(12,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield0_10"], refnameRaaEta + "_0_10Eta", refdataRaaEta);
			book(hEta["PPcross0_10"], refnameRaaEta + "_0_10Eta", refdataRaaEta);
			book(sRaa["RaaEta0_10"], refnameRaaEta);
	
			string refnameRaaEta = mkAxisCode(13,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield10_20"], refnameRaaEta + "_10_20Eta", refdataRaaEta);
			book(hEta["PPcross10_20"], refnameRaaEta + "_10_20Eta", refdataRaaEta);
			book(sRaa["RaaEta10_20"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(14,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield0_20"], refnameRaaEta + "_0_20Eta", refdataRaaEta);
			book(hEta["PPcross0_20"], refnameRaaEta + "_0_20Eta", refdataRaaEta);
			book(sRaa["RaaEta0_20"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(15,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield20_40"], refnameRaaEta + "_20_40Eta", refdataRaaEta);
			book(hEta["PPcross20_40"], refnameRaaEta + "_20_40Eta", refdataRaaEta);
			book(sRaa["RaaEta20_40"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(16,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield40_60"], refnameRaaEta + "_40_60Eta", refdataRaaEta);
			book(hEta["PPcross40_60"], refnameRaaEta + "_40_60Eta", refdataRaaEta);
			book(sRaa["RaaEta40_60"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(17,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield20_60"], refnameRaaEta + "_20_60Eta", refdataRaaEta);
			book(hEta["PPcross20_60"], refnameRaaEta + "_20_60Eta", refdataRaaEta);
			book(sRaa["RaaEta20_60"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(18,1,1);
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield60_92"], refnameRaaEta + "_60_92Eta", refdataRaaEta);
			book(hEta["PPcross60_92"], refnameRaaEta + "_60_92Eta", refdataRaaEta);
			book(sRaa["RaaEta60_92"], refnameRaaEta);

			string refnameRaaEta = mkAxisCode(19,1,1);	//Figures 3-9 and 4-2 in HEPData; should be used twice in finalize
			const Scatter2D& refdataRaaEta = refData(refnameRaaEta);
			book(hEta["AUAUyield0_92"], refnameRaaEta + "_0_92Eta", refdataRaaEta);
			book(hEta["PPcross0_92"], refnameRaaEta + "_0_92Eta", refdataRaaEta);
			book(sRaa["RaaEta0_92"], refnameRaaEta);

			string refnameRaaPi = mkAxisCode(20,1,1);
			const Scatter2D& refdataRaaPi = refData(refnameRaaPi);
			book(hPi["AUAUyield0_92Pi"], refnameRaaPi + "_0_92Pi", refdataRaaPi);
			book(hPi["PPcross0_92Pi"], refnameRaaPi + "_0_92Pi", refdataRaaPi);
			book(sRaa["RaaPi0_92"], refnameRaaPi);

			book(pcross["cross_section"]1,0,1);



		}


		/// Perform the per-event analysis
		void analyze(const Event& event) {

			const ParticlePair& beam = beams();

			if (beam.first.pid() == 1000791970 && beam.secondpid() == 1000791970) collSys = AuAu200;
			
			else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = pp;
			
			const UnstableParticles up = apply<UnstableParticles>(event, "up");
			const double c = cent();
			const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
			const Particles unstableParticles = up.particles();

			if ((c < 0.) || (c > 92.)) vetoEvent;

			if (collSys == pp) {

				for (const Particle& : unstableParticles) {
			
					double partPt = p.pT()/GeV;

					if (p.pid() == 221) {

						PP_yieldEta["yieldPP"]->fill(partPt);
					}
				}
			}

			if (collSys == AuAu200) {

				if ((c >= 0.) && (c < 5.)) {
	
					sow["AUAU0_5"]->fill();

					for (const Particle& : unstableParticles) {
			
						double partPt = p.pT()/GeV;

						if (p.pid == 221) {

							AUAU_yieldEta["yield0_5"]->fill(partPt);
						}
					}
				}

				if ((c >= 0.) && (c < 10.)) {
	
					sow["AUAU0_10"]->fill();

					for (const Particle& : unstableParticles) {
			
						double partPt = p.pT()/GeV;

						if (p.pid == 221) {

							AUAU_yieldEta["yield0_10"]->fill(partPt);
						}
					}
				}

				if ((c >= 10.) && (c < 20.)) {
	
					sow["AUAU10_20"]->fill();

					for (const Particle& : unstableParticles) {
			
						double partPt = p.pT()/GeV;

						if (p.pid == 221) {

							AUAU_yieldEta["yield10_20"]->fill(partPt);
						}
					}
				}

				if ((c >= 0.) && (c < 20.)) {
	
					sow["AUAU0_20"]->fill();

					for (const Particle& : unstableParticles) {
			
						double partPt = p.pT()/GeV;

						if (p.pid == 221) {

							AUAU_yieldEta["yield0_20"]->fill(partPt);
						}
					}
				}

				if ((c >= 20.) && (c < 40.)) {
	
					sow["AUAU20_40"]->fill();

					for (const Particle& : unstableParticles) {
			
						double partPt = p.pT()/GeV;

						if (p.pid == 221) {

							AUAU_yieldEta["yield20_40"]->fill(partPt);
						}
					}
				}

				if ((c >= 40.) && (c < 60.)) {
	
					sow["AUAU40_60"]->fill();

					for (const Particle& : unstableParticles) {
			
						double partPt = p.pT()/GeV;

						if (p.pid == 221) {

							AUAU_yieldEta["yield40_60"]->fill(partPt);
						}
					}
				}

				if ((c >= 20.) && (c < 60.)) {
	
					sow["AUAU20_60"]->fill();

					for (const Particle& : unstableParticles) {
			
						double partPt = p.pT()/GeV;

						if (p.pid == 221) {

							AUAU_yieldEta["yield20_60"]->fill(partPt);
						}
					}
				}

				if ((c >= 60.) && (c < 92.)) {
	
					sow["AUAU60_92"]->fill();

					for (const Particle& : unstableParticles) {
			
						double partPt = p.pT()/GeV;

						if (p.pid == 221) {

							AUAU_yieldEta["yield60_92"]->fill(partPt);
						}
					}
				}

				if ((c >= 0.) && (c < 92.)) {
	
					sow["AUAU0_92"]->fill();

					for (const Particle& : unstableParticles) {
			
						double partPt = p.pT()/GeV;

						if (p.pid == 221) {

							AUAU_yieldEta["yield0_92"]->fill(partPt);
						}
					}
				}
			}
		}

		/// Normalise histograms etc., after the run

		void finalize() {

			if(collSys == pp) pcross["cross_section"]->fill(0.5, crossSection());
			double xsec = pcross["cross_section"]->bin(0).mean()/millibarn;

			//____Yields vs. pT____
			AUAU_yieldEta["yield0_5"]->scaleW(1./sow["AUAU0_5"]);
			AUAU_yieldEta["yield0_10"]->scaleW(1./sow["sow_AUAU0_10"]);
			AUAU_yieldEta["yield10_20"]->scaleW(1./sow["sow_AUAU10_20"]);
			AUAU_yieldEta["yield0_20"]->scaleW(1./sow["sow_AUAU0_20"]);
			AUAU_yieldEta["yield20_40"]->scaleW(1./sow["sow_AUAU20_40"]);
			AUAU_yieldEta["yield40_60"]->scaleW(1./sow["sow_AUAU40_60"]);
			AUAU_yieldEta["yield20_60"]->scaleW(1./sow["sow_AUAU20_60"]);
			AUAU_yieldEta["yield60_92"]->scaleW(1./sow["sow_AUAU60_92"]);
			AUAU_yieldEta["yield0_92"]->scaleW(1./sow["sow_AUAU0_92"]);
			PP_yieldEta["yieldPP"]->scaleW(1./sow["sow_PP"]);


			hEta["AUAUyield0_5"]->scaleW(1./sow["sow_AUAU0_5"]->sumW());
			hEta["PPcross0_5"]->scaleW(Taa*xsec/sow["sow_PP"]->sumW());
			divide(hEta["AUAUyield0_5"], hEta["PPcross0_5"], sRaa["RaaEta0_5"]);
		}

		map<string, Histo1DPtr> AUAU_yieldEta;
		map<string, Histo1DPtr> PP_yieldEta;
		map<string, CounterPtr> sow;
		map<string, Scatter2DPtr> sRaa;
		map<string, Histo1DPtr> hPi;
		map<string, Histo1DPtr> hEta;
		map<string, Profile1DPtr> pcross;

		enum CollisionSystem {empty, pp, AuAu200};
		CollisionSystem collSys = empty;
	};


	RIVET_DECLARE_PLUGIN(PHENIX_2010_I856259);
}
