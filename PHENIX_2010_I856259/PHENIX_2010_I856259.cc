// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FinalState.hh"
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
		
			//Particles: eta meson, Pi0?
			
			beamOpt = getOption<string>("beam", "NONE");

			const UnstableParticles up((Cuts::abspid == 221 || Cuts::abspid == 111) && Cuts::abseta < .35);
			declare(up, "up");

			declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
		
			//____Counters____
			
			book(sow["sow_AUAU0_5"],"_sow_AUAU0_5");
			book(sow["sow_AUAU0_10"],"_sow_AUAU0_10");
			book(sow["sow_AUAU10_20"],"_sow_AUAU10_20");
			book(sow["sow_AUAU0_20"],"_sow_AUAU0_20");
			book(sow["sow_AUAU20_40"],"_sow_AUAU20_40");
			book(sow["sow_AUAU40_60"],"_sow_AUAU40_60");
			book(sow["sow_AUAU20_60"],"_sow_AUAU20_60");
			book(sow["sow_AUAU60_92"],"_sow_AUAU60_92");
			book(sow["sow_AUAU0_92"],"_sow_AUAU0_92");
			book(sow["sow_PP"],"_sow_PP");
			book(sow["sow_AUAU0_92Pi"], "_sow_AUAU0_92Pi");

			//____Yields vs. pT____

			book(AUAU_yieldEta["yield0_5"],1,1,1);
			book(AUAU_yieldEta["yield0_10"],2,1,1);
			book(AUAU_yieldEta["yield10_20"],3,1,1);
			book(AUAU_yieldEta["yield0_20"],4,1,1);
			book(AUAU_yieldEta["yield20_40"],5,1,1);
			book(AUAU_yieldEta["yield40_60"],6,1,1);
			book(AUAU_yieldEta["yield20_60"],7,1,1);
			book(AUAU_yieldEta["yield60_92"],8,1,1);
			book(AUAU_yieldEta["yield0_92"],9,1,1);
			book(PP_yieldEta["yieldPP"],10,1,1);


			//____Eta Raa____

			string refnameRaaEta0_5 = mkAxisCode(11,1,1);
			const Scatter2D& refdataRaaEta0_5 = refData(refnameRaaEta0_5);
			book(hEta["AUAUyield0_5"], refnameRaaEta0_5 + "_0_5Eta", refdataRaaEta0_5);
			book(hEta["PPcross0_5"], refnameRaaEta0_5 + "_0_5Pcross", refdataRaaEta0_5);
			book(sRaa["RaaEta0_5"], refnameRaaEta0_5);

			string refnameRaaEta0_10 = mkAxisCode(12,1,1);
			const Scatter2D& refdataRaaEta0_10 = refData(refnameRaaEta0_10);
			book(hEta["AUAUyield0_10"], refnameRaaEta0_10 + "_0_10Eta", refdataRaaEta0_10);
			book(hEta["PPcross0_10"], refnameRaaEta0_10 + "_0_10Pcross", refdataRaaEta0_10);
			book(sRaa["RaaEta0_10"], refnameRaaEta0_10);
	
			string refnameRaaEta10_20 = mkAxisCode(13,1,1);
			const Scatter2D& refdataRaaEta10_20 = refData(refnameRaaEta10_20);
			book(hEta["AUAUyield10_20"], refnameRaaEta10_20 + "_10_20Eta", refdataRaaEta10_20);
			book(hEta["PPcross10_20"], refnameRaaEta10_20 + "_10_20Pcross", refdataRaaEta10_20);
			book(sRaa["RaaEta10_20"], refnameRaaEta10_20);

			string refnameRaaEta0_20 = mkAxisCode(14,1,1);
			const Scatter2D& refdataRaaEta0_20 = refData(refnameRaaEta0_20);
			book(hEta["AUAUyield0_20"], refnameRaaEta0_20 + "_0_20Eta", refdataRaaEta0_20);
			book(hEta["PPcross0_20"], refnameRaaEta0_20 + "_0_20Pcross", refdataRaaEta0_20);
			book(sRaa["RaaEta0_20"], refnameRaaEta0_20);

			string refnameRaaEta20_40 = mkAxisCode(15,1,1);
			const Scatter2D& refdataRaaEta20_40 = refData(refnameRaaEta20_40);
			book(hEta["AUAUyield20_40"], refnameRaaEta20_40 + "_20_40Eta", refdataRaaEta20_40);
			book(hEta["PPcross20_40"], refnameRaaEta20_40 + "_20_40Pcross", refdataRaaEta20_40);
			book(sRaa["RaaEta20_40"], refnameRaaEta20_40);

			string refnameRaaEta40_60 = mkAxisCode(16,1,1);
			const Scatter2D& refdataRaaEta40_60 = refData(refnameRaaEta40_60);
			book(hEta["AUAUyield40_60"], refnameRaaEta40_60 + "_40_60Eta", refdataRaaEta40_60);
			book(hEta["PPcross40_60"], refnameRaaEta40_60 + "_40_60Pcross", refdataRaaEta40_60);
			book(sRaa["RaaEta40_60"], refnameRaaEta40_60);

			string refnameRaaEta20_60 = mkAxisCode(17,1,1);
			const Scatter2D& refdataRaaEta20_60 = refData(refnameRaaEta20_60);
			book(hEta["AUAUyield20_60"], refnameRaaEta20_60 + "_20_60Eta", refdataRaaEta20_60);
			book(hEta["PPcross20_60"], refnameRaaEta20_60 + "_20_60Pcross", refdataRaaEta20_60);
			book(sRaa["RaaEta20_60"], refnameRaaEta20_60);

			string refnameRaaEta60_92 = mkAxisCode(18,1,1);
			const Scatter2D& refdataRaaEta60_92 = refData(refnameRaaEta60_92);
			book(hEta["AUAUyield60_92"], refnameRaaEta60_92 + "_60_92Eta", refdataRaaEta60_92);
			book(hEta["PPcross60_92"], refnameRaaEta60_92 + "_60_92Pcross", refdataRaaEta60_92);
			book(sRaa["RaaEta60_92"], refnameRaaEta60_92);

			string refnameRaaEta0_92 = mkAxisCode(19,1,1);	//Figures 3-9 and 4-2 in HEPData; should maybe be used twice in finalize
			const Scatter2D& refdataRaaEta0_92 = refData(refnameRaaEta0_92);
			book(hEta["AUAUyield0_92"], refnameRaaEta0_92 + "_0_92Eta", refdataRaaEta0_92);
			book(hEta["PPcross0_92"], refnameRaaEta0_92 + "_0_92Pcross", refdataRaaEta0_92);
			book(sRaa["RaaEta0_92"], refnameRaaEta0_92);

			string refnameRaaPi0_92 = mkAxisCode(20,1,1);
			const Scatter2D& refdataRaaPi0_92 = refData(refnameRaaPi0_92);
			book(hPi["AUAUyield0_92Pi"], refnameRaaPi0_92 + "_0_92Pi", refdataRaaPi0_92);
			book(hPi["PPcross0_92Pi"], refnameRaaPi0_92 + "_0_92PiPcross", refdataRaaPi0_92);
			book(sRaa["RaaPi0_92"], refnameRaaPi0_92);

			//string refnameRaaEtapT20 = mkAxisCode(22,1,1);
			//const Scatter2D& refdataRaaEtapT20 = refData(refnameRaaEtapT20);
			//book(hPi["AUAUyieldEtapT20"], refnameRaaEtapT20 + "_EtapT20", refdataRaaEtapT20);
			//book(hPi["PPcrossEtapT20"], refnameRaaEtapT20 + "_PcrossEtapT20", refdataRaaEtapT20);
			//book(sRaa["RaaEtapT20"], refnameRaaEtapT20);

			//string refnameRaaEtapT5 = mkAxisCode(23,1,1);
			//const Scatter2D& refdataRaaEtapT5 = refData(refnameRaaEtapT5);
			//book(hPi["AUAUyieldEtapT5"], refnameRaaEtapT5 + "_EtapT5", refdataRaaEtapT5);
			//book(hPi["PPcrossEtapT5"], refnameRaaEtapT5 + "_PcrossEtapT5", refdataRaaEtapT5);
			//book(sRaa["RaaEtapT5"], refnameRaaEtapT5);


			book(pcross["cross_section"],"cross_section",1,0,1);	//cross section binning profile



		}


		/// Perform the per-event analysis
		void analyze(const Event& event) {

			const ParticlePair& beam = beams();

			if (beamOpt == "NONE") {

				if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu200;
			
				else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = pp;
			}

			else if (beamOpt == "PP200") collSys = pp;
			else if (beamOpt == "AUAU200") collSys = AuAu200;

			const UnstableParticles up = apply<UnstableParticles>(event, "up");
			const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
			const double c = cent();
			const Particles unstableParticles = up.particles();

			if ((c < 0.) || (c > 92.)) vetoEvent;

			if (collSys == pp) {

				sow["sow_PP"]->fill();
			
				for (const Particle& p : unstableParticles) {
			
					double partPt = p.pT()/GeV;
					double pt_weight = 1./(partPt*2.*M_PI);

					if (p.pid() == 221) {

						PP_yieldEta["yieldPP"]->fill(partPt, pt_weight);
					}
				}
			}

			if (collSys == AuAu200) {

				if ((c >= 0.) && (c < 5.)) {
	
					sow["sow_AUAU0_5"]->fill();

					for (const Particle& p : unstableParticles) {
			
						double partPt = p.pT()/GeV;
						double pt_weight = 1./(partPt*2.*M_PI);


						if (p.pid() == 221) {

							AUAU_yieldEta["yield0_5"]->fill(partPt, pt_weight);
							hEta["AUAUyield0_5"]->fill(partPt, pt_weight);
							hEta["PPcross0_5"]->fill(partPt, pt_weight);
						}
					}
				}

				if ((c >= 0.) && (c < 10.)) {
	
					sow["sow_AUAU0_10"]->fill();

					for (const Particle& p : unstableParticles) {
			
						double partPt = p.pT()/GeV;
						double pt_weight = 1./(partPt*2.*M_PI);

						if (p.pid() == 221) {

							AUAU_yieldEta["yield0_10"]->fill(partPt, pt_weight);
							hEta["AUAUyield0_10"]->fill(partPt, pt_weight);
							hEta["PPcross0_10"]->fill(partPt, pt_weight);
						}
					}
				}

				if ((c >= 10.) && (c < 20.)) {
	
					sow["sow_AUAU10_20"]->fill();

					for (const Particle& p : unstableParticles) {
			
						double partPt = p.pT()/GeV;
						double pt_weight = 1./(partPt*2.*M_PI);

						if (p.pid() == 221) {

							AUAU_yieldEta["yield10_20"]->fill(partPt, pt_weight);
							hEta["AUAUyield10_20"]->fill(partPt, pt_weight);
							hEta["PPcross10_20"]->fill(partPt, pt_weight);
						}
					}
				}

				if ((c >= 0.) && (c < 20.)) {
	
					sow["sow_AUAU0_20"]->fill();

					for (const Particle& p : unstableParticles) {
			
						double partPt = p.pT()/GeV;
						double pt_weight = 1./(partPt*2.*M_PI);

						if (p.pid() == 221) {

							AUAU_yieldEta["yield0_20"]->fill(partPt, pt_weight);
							hEta["AUAUyield0_20"]->fill(partPt, pt_weight);
							hEta["PPcross0_20"]->fill(partPt, pt_weight);
						}
					}
				}

				if ((c >= 20.) && (c < 40.)) {
	
					sow["sow_AUAU20_40"]->fill();

					for (const Particle& p : unstableParticles) {
			
						double partPt = p.pT()/GeV;
						double pt_weight = 1./(partPt*2.*M_PI);

						if (p.pid() == 221) {

							AUAU_yieldEta["yield20_40"]->fill(partPt, pt_weight);
							hEta["AUAUyield20_40"]->fill(partPt, pt_weight);
							hEta["PPcross20_40"]->fill(partPt, pt_weight);
						}
					}
				}

				if ((c >= 40.) && (c < 60.)) {
	
					sow["sow_AUAU40_60"]->fill();

					for (const Particle& p : unstableParticles) {
			
						double partPt = p.pT()/GeV;
						double pt_weight = 1./(partPt*2.*M_PI);

						if (p.pid() == 221) {

							AUAU_yieldEta["yield40_60"]->fill(partPt, pt_weight);
							hEta["AUAUyield40_60"]->fill(partPt, pt_weight);
							hEta["PPcross40_60"]->fill(partPt, pt_weight);
						}
					}
				}

				if ((c >= 20.) && (c < 60.)) {
	
					sow["sow_AUAU20_60"]->fill();

					for (const Particle& p : unstableParticles) {
			
						double partPt = p.pT()/GeV;
						double pt_weight = 1./(partPt*2.*M_PI);

						if (p.pid() == 221) {

							AUAU_yieldEta["yield20_60"]->fill(partPt, pt_weight);
							hEta["AUAUyield20_60"]->fill(partPt, pt_weight);
							hEta["PPcross20_60"]->fill(partPt, pt_weight);
						}
					}
				}

				if ((c >= 60.) && (c < 92.)) {
	
					sow["sow_AUAU60_92"]->fill();

					for (const Particle& p : unstableParticles) {
			
						double partPt = p.pT()/GeV;
						double pt_weight = 1./(partPt*2.*M_PI);

						if (p.pid() == 221) {

							AUAU_yieldEta["yield60_92"]->fill(partPt, pt_weight);
							hEta["AUAUyield60_92"]->fill(partPt, pt_weight);
							hEta["PPcross60_92"]->fill(partPt, pt_weight);
						}
					}
				}

				if ((c >= 0.) && (c <= 92.)) {
	
					sow["sow_AUAU0_92"]->fill();

					for (const Particle& p : unstableParticles) {
			
						double partPt = p.pT()/GeV;
						double pt_weight = 1./(partPt*2.*M_PI);

						if (p.pid() == 221) {

							AUAU_yieldEta["yield0_92"]->fill(partPt, pt_weight);
							hEta["AUAUyield0_92"]->fill(partPt, pt_weight);
							hEta["PPcross0_92"]->fill(partPt, pt_weight);
						}

						if (p.pid() == 111) {

							hPi["AUAUyield0_92Pi"]->fill(partPt, pt_weight);
							hPi["PPcross0_92Pi"]->fill(partPt, pt_weight);
						}
					}
				}
			}
		}

		/// Normalise histograms etc., after the run

		void finalize() {

			if (collSys == pp) pcross["cross_section"]->fill(0.5, crossSection());
			
			double xsec = 1;

			if (pcross["cross_section"]->bin(0).numEntries() > 0) {

				xsec = pcross["cross_section"]->bin(0).mean()/millibarn;
			}


			binShift(*hEta["AUAUyield0_5"]); 
			binShift(*hEta["PPcross0_5"]);
		
		
			binShift(*hEta["AUAUyield0_10"]);
			binShift(*hEta["PPcross0_10"]);
			

			binShift(*hEta["AUAUyield10_20"]);
			binShift(*hEta["PPcross10_20"]);

			binShift(*hEta["AUAUyield0_20"]);
			binShift(*hEta["PPcross0_20"]);
			

			binShift(*hEta["AUAUyield20_40"]);
			binShift(*hEta["PPcross20_40"]);

			binShift(*hEta["AUAUyield40_60"]);
			binShift(*hEta["PPcross40_60"]);

			binShift(*hEta["AUAUyield20_60"]);
			binShift(*hEta["PPcross20_60"]);

			binShift(*hEta["AUAUyield60_92"]);
			binShift(*hEta["PPcross60_92"]);

			binShift(*hEta["AUAUyield0_92"]);
			binShift(*hEta["PPcross0_92"]);
			binShift(*hPi["AUAUyield0_92Pi"]);
			binShift(*hPi["PPcross0_92Pi"]);

			//____Yields vs. pT____
			AUAU_yieldEta["yield0_5"]->scaleW(1./sow["sow_AUAU0_5"]->sumW());
			AUAU_yieldEta["yield0_10"]->scaleW(1./sow["sow_AUAU0_10"]->sumW());
			AUAU_yieldEta["yield10_20"]->scaleW(1./sow["sow_AUAU10_20"]->sumW());
			AUAU_yieldEta["yield0_20"]->scaleW(1./sow["sow_AUAU0_20"]->sumW());
			AUAU_yieldEta["yield20_40"]->scaleW(1./sow["sow_AUAU20_40"]->sumW());
			AUAU_yieldEta["yield40_60"]->scaleW(1./sow["sow_AUAU40_60"]->sumW());
			AUAU_yieldEta["yield20_60"]->scaleW(1./sow["sow_AUAU20_60"]->sumW());
			AUAU_yieldEta["yield60_92"]->scaleW(1./sow["sow_AUAU60_92"]->sumW());
			AUAU_yieldEta["yield0_92"]->scaleW(1./sow["sow_AUAU0_92"]->sumW());
			PP_yieldEta["yieldPP"]->scaleW(1./sow["sow_PP"]->sumW());

			//____Raa____
			
		if (sow["sow_AUAU0_5"]->sumW() != 0) {
    		hEta["AUAUyield0_5"]->scaleW(1. / sow["sow_AUAU0_5"]->sumW());
   		 	hEta["PPcross0_5"]->scaleW(25.37 * xsec / sow["sow_PP"]->sumW());
   		 	divide(hEta["AUAUyield0_5"], hEta["PPcross0_5"], sRaa["RaaEta0_5"]);
		} else {
    		std::cerr << "Error: Divide by zero encountered for sow_AUAU0_5. Unable to scale histogram." << std::endl;
}

		if (sow["sow_AUAU0_10"]->sumW() != 0) {
    		hEta["AUAUyield0_10"]->scaleW(1. / sow["sow_AUAU0_10"]->sumW());
    		hEta["PPcross0_10"]->scaleW(22.75 * xsec / sow["sow_PP"]->sumW());
    		divide(hEta["AUAUyield0_10"], hEta["PPcross0_10"], sRaa["RaaEta0_10"]);
		} else {
    		std::cerr << "Error: Divide by zero encountered for sow_AUAU0_10. Unable to scale histogram." << std::endl;
}

		if (sow["sow_AUAU10_20"]->sumW() != 0) {
   			hEta["AUAUyield10_20"]->scaleW(1. / sow["sow_AUAU10_20"]->sumW());
    		hEta["PPcross10_20"]->scaleW(14.35 * xsec / sow["sow_PP"]->sumW());
   			divide(hEta["AUAUyield10_20"], hEta["PPcross10_20"], sRaa["RaaEta10_20"]);
		} else {
    		std::cerr << "Error: Divide by zero encountered for sow_AUAU10_20. Unable to scale histogram." << std::endl;
}

		if (sow["sow_AUAU0_20"]->sumW() != 0) {
    		hEta["AUAUyield0_20"]->scaleW(1. / sow["sow_AUAU0_20"]->sumW());
    		hEta["PPcross0_20"]->scaleW(18.55 * xsec / sow["sow_PP"]->sumW());
    		divide(hEta["AUAUyield0_20"], hEta["PPcross0_20"], sRaa["RaaEta0_20"]);
		} else {
    		std::cerr << "Error: Divide by zero encountered for sow_AUAU0_20. Unable to scale histogram." << std::endl;
}

		if (sow["sow_AUAU20_40"]->sumW() != 0) {
    		hEta["AUAUyield20_40"]->scaleW(1. / sow["sow_AUAU20_40"]->sumW());
    		hEta["PPcross20_40"]->scaleW(7.065 * xsec / sow["sow_PP"]->sumW());
    		divide(hEta["AUAUyield20_40"], hEta["PPcross20_40"], sRaa["RaaEta20_40"]);
		} else {
    		std::cerr << "Error: Divide by zero encountered for sow_AUAU20_40. Unable to scale histogram." << std::endl;
}

		if (sow["sow_AUAU40_60"]->sumW() != 0) {
   			hEta["AUAUyield40_60"]->scaleW(1. / sow["sow_AUAU40_60"]->sumW());
    		hEta["PPcross40_60"]->scaleW(2.155 * xsec / sow["sow_PP"]->sumW());
    		divide(hEta["AUAUyield40_60"], hEta["PPcross40_60"], sRaa["RaaEta40_60"]);
		} else {
    		std::cerr << "Error: Divide by zero encountered for sow_AUAU40_60. Unable to scale histogram." << std::endl;
}		

		if (sow["sow_AUAU20_60"]->sumW() != 0) {
   			hEta["AUAUyield20_60"]->scaleW(1. / sow["sow_AUAU20_60"]->sumW());
    		hEta["PPcross20_60"]->scaleW(4.61 * xsec / sow["sow_PP"]->sumW());
    		divide(hEta["AUAUyield20_60"], hEta["PPcross20_60"], sRaa["RaaEta20_60"]);
} 		else {
    		std::cerr << "Error: Divide by zero encountered for sow_AUAU20_60. Unable to scale histogram." << std::endl;
}

		if (sow["sow_AUAU60_92"]->sumW() != 0) {
    		hEta["AUAUyield60_92"]->scaleW(1. / sow["sow_AUAU60_92"]->sumW());
    		hEta["PPcross60_92"]->scaleW(.35 * xsec / sow["sow_PP"]->sumW());
    		divide(hEta["AUAUyield60_92"], hEta["PPcross60_92"], sRaa["RaaEta60_92"]);
		} else {
    		std::cerr << "Error: Divide by zero encountered for sow_AUAU60_92. Unable to scale histogram." << std::endl;
}

		if (sow["sow_AUAU0_92"]->sumW() != 0) {
    		hEta["AUAUyield0_92"]->scaleW(1. / sow["sow_AUAU0_92"]->sumW());
    		hEta["PPcross0_92"]->scaleW(6.14 * xsec / sow["sow_PP"]->sumW());
    		divide(hEta["AUAUyield0_92"], hEta["PPcross0_92"], sRaa["RaaEta0_92"]);
		} else {
    		std::cerr << "Error: Divide by zero encountered for sow_AUAU0_92. Unable to scale histogram." << std::endl;
}

		if (sow["sow_AUAU0_92Pi"]->sumW() != 0) {
    		hPi["AUAUyield0_92Pi"]->scaleW(1. / sow["sow_AUAU0_92Pi"]->sumW());
    		hPi["PPcross0_92Pi"]->scaleW(6.14 * xsec / sow["sow_PP"]->sumW());
    		divide(hPi["AUAUyield0_92Pi"], hPi["PPcross0_92Pi"], sRaa["RaaPi0_92"]);
		} else {
    		std::cerr << "Error: Divide by zero encountered for sow_AUAU0_92Pi. Unable to scale histogram." << std::endl;
}

			hEta["AUAUyield0_5"]->scaleW(1./sow["sow_AUAU0_5"]->sumW());
			hEta["PPcross0_5"]->scaleW(25.37*xsec/sow["sow_PP"]->sumW());
			divide(hEta["AUAUyield0_5"], hEta["PPcross0_5"], sRaa["RaaEta0_5"]);
		
			hEta["AUAUyield0_10"]->scaleW(1./sow["sow_AUAU0_10"]->sumW());
			hEta["PPcross0_10"]->scaleW(22.75*xsec/sow["sow_PP"]->sumW());
			divide(hEta["AUAUyield0_10"], hEta["PPcross0_10"], sRaa["RaaEta0_10"]);

			hEta["AUAUyield10_20"]->scaleW(1./sow["sow_AUAU10_20"]->sumW());
			hEta["PPcross10_20"]->scaleW(14.35*xsec/sow["sow_PP"]->sumW());
			divide(hEta["AUAUyield10_20"], hEta["PPcross10_20"], sRaa["RaaEta10_20"]);

			hEta["AUAUyield0_20"]->scaleW(1./sow["sow_AUAU0_20"]->sumW());
			hEta["PPcross0_20"]->scaleW(18.55*xsec/sow["sow_PP"]->sumW());
			divide(hEta["AUAUyield0_20"], hEta["PPcross0_20"], sRaa["RaaEta0_20"]);

			hEta["AUAUyield20_40"]->scaleW(1./sow["sow_AUAU20_40"]->sumW());
			hEta["PPcross20_40"]->scaleW(7.065*xsec/sow["sow_PP"]->sumW());
			divide(hEta["AUAUyield20_40"], hEta["PPcross20_40"], sRaa["RaaEta20_40"]);

			hEta["AUAUyield40_60"]->scaleW(1./sow["sow_AUAU40_60"]->sumW());
			hEta["PPcross40_60"]->scaleW(2.155*xsec/sow["sow_PP"]->sumW());
			divide(hEta["AUAUyield40_60"], hEta["PPcross40_60"], sRaa["RaaEta40_60"]);

			hEta["AUAUyield20_60"]->scaleW(1./sow["sow_AUAU20_60"]->sumW());
			hEta["PPcross20_60"]->scaleW(4.61*xsec/sow["sow_PP"]->sumW());
			divide(hEta["AUAUyield20_60"], hEta["PPcross20_60"], sRaa["RaaEta20_60"]);

			hEta["AUAUyield60_92"]->scaleW(1./sow["sow_AUAU60_92"]->sumW());
			hEta["PPcross60_92"]->scaleW(.35*xsec/sow["sow_PP"]->sumW());
			divide(hEta["AUAUyield60_92"], hEta["PPcross60_92"], sRaa["RaaEta60_92"]);

			hEta["AUAUyield0_92"]->scaleW(1./sow["sow_AUAU0_92"]->sumW());
			hEta["PPcross0_92"]->scaleW(6.14*xsec/sow["sow_PP"]->sumW());
			divide(hEta["AUAUyield0_92"], hEta["PPcross0_92"], sRaa["RaaEta0_92"]);

			hPi["AUAUyield0_92Pi"]->scaleW(1./sow["sow_AUAU0_92Pi"]->sumW());
			hPi["PPcross0_92Pi"]->scaleW(6.14*xsec/sow["sow_PP"]->sumW());
			divide(hPi["AUAUyield0_92Pi"], hPi["PPcross0_92Pi"], sRaa["RaaPi0_92"]);

			//Taa values
			//0-5: 25.37
			//0-10: 22.75
			//10-20: 14.35
			//0-20: 18.55
			//20-40: 7.065
			//40-60: 2.155
			//20-60: 4.61
			//60-92: 0.35
			//0-92: 6.14
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
		string beamOpt = "NONE";
	};


	RIVET_DECLARE_PLUGIN(PHENIX_2010_I856259);
}
