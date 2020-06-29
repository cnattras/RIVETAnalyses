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
#include "Centrality/RHICCentrality.hh" //external header for Centrality calculation
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES

namespace Rivet {


	class PHENIX_2013_I1227971 : public Analysis {
	public:

		DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2013_I1227971);


		void init() {

			//Need to verify if this can be done like this
			std::initializer_list<int> pdgIds = { 321 };  // Kaon +
			std::initializer_list<int> pdgIds = { -321 };  // Kaon -
			std::initializer_list<int> pdgIds = { 211 };  // Pion +
			std::initializer_list<int> pdgIds = { -211 };  // Pion -
			std::initializer_list<int> pdgIds = { 2212 };  // Proton
			std::initializer_list<int> pdgIds = { -2212 };  // Proton Bar

			const PrimaryParticles fs(pdgIds, Cuts::abseta < 0.35 && Cuts::abscharge < 0);
			declare(fs, "fs");

			beamOpt = getOption<string>("beam", "NONE");

			if (beamOpt == "PP") collSys = pp;
			else if (beamOpt == "AUAU200") collSys = AuAu200;
			else if (beamOpt == "dAU200") collSys = dAu200;


			if (!(collSys == pp)) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

			book(sow["sow_pp"], "sow_pp");
			for (int i = 0, N = AUAUCentralityBins.size(); i < N; ++i)
			{
				//yields (fig 4)_________________
				book(hKaonNegPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"], 1, 1, 1 + i);
				book(hPionNegPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"], 3, 1, 1 + i);
				book(hProtNegPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"], 5, 1, 1 + i);
				book(hKaonPosPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"], 1, 1, 6 + i);
				book(hPionPosPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"], 3, 1, 6 + i);
				book(hProtPosPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"], 5, 1, 6 + i);
				book(hKaonNegPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"], 2, 1, 1 + i);
				book(hPionNegPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"], 4, 1, 1 + i);
				book(hProtNegPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"], 6, 1, 1 + i);
				book(hKaonPosPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"], 2, 1, 6 + i);
				book(hPionPosPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"], 4, 1, 6 + i);
				book(hProtPosPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"], 6, 1, 6 + i);

				book(sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"], "sow_AUAUc + std::to_string(AUAUCentralityBins[i])");
				book(sow["sow_dAUc + std::to_string(dAUCentralityBins[i])"], "sow_dAUc + std::to_string(dAUCentralityBins[i])");
			

				//Ratio of yields (figs 5-9)_________________
			
				//Histograms for the ratios Neg/Pos as a function of pT in different centrality classes				
				string refnameK_K = mkAxisCode(7, 1, 1 + i);
					const Scatter2D& refdataK_K = refData(refnameK_K);
				book(rKK_KaonNeg["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnameK_K + "_KaonNeg", refdataK_K);
				book(rKK_KaonPos["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnameK_K + "_KaonPos", refdataK_K);
				book(RatioKaon["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnameK_K);

				string refnameK_K = mkAxisCode(8, 1, 1 + i);
					const Scatter2D& refdataK_K = refData(refnameK_K);
				book(rKK_KaonNeg["dAuc + std::to_string(dAUCentralityBins[i])"], refnameK_K + "_KaonNeg", refdataK_K);
				book(rKK_KaonPos["dAuc + std::to_string(dAUCentralityBins[i])"], refnameK_K + "_KaonPos", refdataK_K);
				book(RatioKaon["dAuc + std::to_string(dAUCentralityBins[i])"], refnameK_K);
				
				string refnamepi_pi = mkAxisCode(9, 1, 1 + i);
					const Scatter2D& refdatapi_pi = refData(refnamepi_pi);
				book(rKK_PionNeg["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnamepi_pi + "_PionNeg", refdatapi_pi);
				book(rKK_PionPos["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnamepi_pi + "_PionPos", refdatapi_pi);
				book(RatioPion["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnamepi_pi);

				string refnamepi_pi = mkAxisCode(10, 1, 1 + i);
					const Scatter2D& refdatapi_pi = refData(refnamepi_pi);
				book(rKK_PionNeg["dAuc + std::to_string(dAUCentralityBins[i])"], refnamepi_pi + "_PionNeg", refdatapi_pi);
				book(rKK_PionPos["dAuc + std::to_string(dAUCentralityBins[i])"], refnamepi_pi + "_PionPos", refdatapi_pi);
				book(RatioPion["dAuc + std::to_string(dAUCentralityBins[i])"], refnamepi_pi);

				string refnamep_p = mkAxisCode(11, 1, 1 + i);
					const Scatter2D& refdatap_p = refData(refnamep_p);
				book(rKK_ProtNeg["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnamep_p + "_ProtNeg", refdatap_p);
				book(rKK_ProtPos["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnamep_p + "_ProtPos", refdatap_p);
				book(RatioProt["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnamep_p);

				string refnamep_p = mkAxisCode(12, 1, 1 + i);
					const Scatter2D& refdatap_p = refData(refnamep_p);
				book(rKK_ProtNeg["dAuc + std::to_string(dAUCentralityBins[i])"], refnamep_p + "_ProtNeg", refdatap_p);
				book(rKK_ProtPos["dAuc + std::to_string(dAUCentralityBins[i])"], refnamep_p + "_ProtPos", refdatap_p);
				book(RatioProt["dAuc + std::to_string(dAUCentralityBins[i])"], refnamep_p);
						
				//Histograms for the ratios Kaon/Pion as a function of pT in different centrality classes
				string refnameK_pi = mkAxisCode(13, 1, 1 + i);
					const Scatter2D& refdataK_pi = refData(refnameK_pi);
				book(rKpi_Kaon["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnameK_pi + "_kaons", refdataK_pi);
				book(rKpi_Pion["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnameK_pi + "_pions", refdataK_pi);
				book(RatioK_pi["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnameK_pi);

				string refnameK_pi = mkAxisCode(14, 1, 1 + i);
					const Scatter2D& refdataK_pi = refData(refnameK_pi);
				book(rKpi_Kaon["dAuc + std::to_string(dAUCentralityBins[i])"], refnameK_pi + "_kaons", refdataK_pi);
				book(rKpi_Pion["dAuc + std::to_string(dAUCentralityBins[i])"], refnameK_pi + "_pions", refdataK_pi);
				book(RatioK_pi["dAuc + std::to_string(dAUCentralityBins[i])"], refnameK_pi);

				//Histograms for the ratios Proton/Pion as a function of pT in different centrality classes
				string refnamep_pi = mkAxisCode(15, 1, 1 + i);
					const Scatter2D& refdatap_pi = refData(refnamep_pi);
				book(rppi_Proton["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnamep_pi + "_protons", refdatap_pi);
				book(rppi_Pion["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnamep_pi + "_pions", refdatap_pi);
				book(ratiop_pi["AuAuc + std::to_string(AUAUCentralityBins[i])"], refnamep_pi);
				
				string refnamep_pi = mkAxisCode(16, 1, 1 + i);
					const Scatter2D& refdatap_pi = refData(refnamep_pi);
				book(rppi_Proton["dAuc + std::to_string(AUAUCentralityBins[i])"], refnamep_pi + "_protons", refdatap_pi);
				book(rppi_Pion["dAuc + std::to_string(AUAUCentralityBins[i])"], refnamep_pi + "_pions", refdatap_pi);
				book(ratiop_pi["dAuc + std::to_string(AUAUCentralityBins[i])"], refnamep_pi);

				//RCP (fig 10) _________________(Assuming booking is the same as RAA)
			
				string refnameRaa = mkAxisCode(17, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hKaonPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa + "_AuAu", refdataRaa);
				book(hKaonPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRcp["K_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa);

				string refnameRaa = mkAxisCode(18, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hPionPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa + "_AuAu", refdataRaa);
				book(hPionPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRcp["pi_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa);

				string refnameRaa = mkAxisCode(19, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hProtPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa + "_AuAu", refdataRaa);
				book(hProtPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRcp["p_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa);
			

				//RAA (fig 11)_________________
				
				string refnameRaa = mkAxisCode(20, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hKaonPt["Raa_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa + "_AuAu", refdataRaa);
				book(hKaonPt["Raa_c + std::to_string(AUAUCentralityBins[i])_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRaa["K_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa);

				string refnameRaa = mkAxisCode(21, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hPionPt["Raa_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa + "_AuAu", refdataRaa);
				book(hPionPt["Raa_c + std::to_string(AUAUCentralityBins[i])_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRaa["pi_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa);

				string refnameRaa = mkAxisCode(22, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hProtPt["Raa_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa + "_AuAu", refdataRaa);
				book(hProtPt["Raa_c + std::to_string(AUAUCentralityBins[i])_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRaa["p_c + std::to_string(AUAUCentralityBins[i])_AuAu"], refnameRaa);
			

				//RdA (fig 12)_________________(Assuming booking is the same as RAA)
			
				string refnameRaa = mkAxisCode(23, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hKaonPt["Rda_c + std::to_string(dAUCentralityBins[i])_dAu"], refnameRaa + "_dAu", refdataRaa);
				book(hKaonPt["Rda_c + std::to_string(dAUCentralityBins[i])_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRda["K_c + std::to_string(dAUCentralityBins[i])_dAu"], refnameRaa);

				string refnameRaa = mkAxisCode(24, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hPionPt["Rda_c + std::to_string(dAUCentralityBins[i])_dAu"], refnameRaa + "_dAu", refdataRaa);
				book(hPionPt["Rda_c + std::to_string(dAUCentralityBins[i])_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRda["pi_c + std::to_string(dAUCentralityBins[i])_dAu"], refnameRaa);

				string refnameRaa = mkAxisCode(25, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hProtPt["Rda_c + std::to_string(dAUCentralityBins[i])_dAu"], refnameRaa + "_dAu", refdataRaa);
				book(hProtPt["Rda_c + std::to_string(dAUCentralityBins[i])_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRda["p_c + std::to_string(dAUCentralityBins[i])_dAu"], refnameRaa);


				// Ratio of Spectra(fig 15)_________________(Need to book)


			}

		}


	void analyze(const Event& event) {


		Particles chargedParticles = applyProjection<PrimaryParticles>(event, "fs").particles();

		if (collSys == pp)
		{
			sow["sow_pp"]->fill();
			for (Particle p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);
				hPion0Pt["c" + std::to_string(CentralityBins[i + 2]) + "pt_pp"]->fill(partPt, pt_weight);
			}
			return;
		}


		const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
		const double c = cent();


		if (collSys == AuAu200)
		{
			if ((c < 0.) || (c > 92.)) vetoEvent;
			
			//sow["sow_AuAu"]->fill();
			double partPt = p.pT() / GeV;
			double pt_weight = 1. / (partPt * 2. * M_PI);

			for (const Particle& p : chargedParticles)
			{

				//hPion0Pt["Pion0Pt_AuAu"]->fill(partPt, pt_weight);
			}



			switch (p.pid()) {
			case 211: // pi+
			{
				if ((c >= 0.) && (c < 10.))
				{
					sow["sow_AUAUc10"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionPosPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);	
					}
				}
				
				else if ((c >= 10.) && (c < 20.))
				{
					sow["sow_AUAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionPosPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_AUAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionPosPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
					}
				}
				
				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_AUAUc60"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionPosPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 60.) && (c < 92.))
				{
					sow["sow_AUAUc92"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionPosPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
					}
				}

				break;
			}
			case -211: // pi-
			{
				if ((c >= 0.) && (c < 10.))
				{
					sow["sow_AUAUc10"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionNegPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 10.) && (c < 20.))
				{
					sow["sow_AUAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionNegPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_AUAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionNegPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_AUAUc60"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionNegPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 60.) && (c < 92.))
				{
					sow["sow_AUAUc92"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionNegPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
					}
				}

				break;
			}
			case 321: // K+
			{
				if ((c >= 0.) && (c < 10.))
				{
					sow["sow_AUAUc10"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonPosPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 10.) && (c < 20.))
				{
					sow["sow_AUAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonPosPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_AUAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonPosPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_AUAUc60"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonPosPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 60.) && (c < 92.))
				{
					sow["sow_AUAUc92"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonPosPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
					}
				}

				break;
			}
			case -321: // K-
			{
				if ((c >= 0.) && (c < 10.))
				{
					sow["sow_AUAUc10"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonNegPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 10.) && (c < 20.))
				{
					sow["sow_AUAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonNegPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_AUAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonNegPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_AUAUc60"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonNegPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 60.) && (c < 92.))
				{
					sow["sow_AUAUc92"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonNegPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
					}
				}

				break;
			}
			case 2212: // proton
			{
				if ((c >= 0.) && (c < 10.))
				{
					sow["sow_AUAUc10"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtPosPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 10.) && (c < 20.))
				{
					sow["sow_AUAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtPosPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_AUAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtPosPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_AUAUc60"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtPosPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 60.) && (c < 92.))
				{
					sow["sow_AUAUc92"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtPosPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
					}
				}

				break;
			}
			case -2212: // anti-proton
			{
				if ((c >= 0.) && (c < 10.))
				{
					sow["sow_AUAUc10"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtNegPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 10.) && (c < 20.))
				{
					sow["sow_AUAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtNegPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_AUAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtNegPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_AUAUc60"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtNegPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 60.) && (c < 92.))
				{
					sow["sow_AUAUc92"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtNegPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
					}
				}

				break;
			}
			}
		}

		}
		if (collSys == dAu200)
		{
			if ((c < 0.) || (c > 100.)) vetoEvent;

			//sow["sow_AuAu"]->fill();
			double partPt = p.pT() / GeV;
			double pt_weight = 1. / (partPt * 2. * M_PI);

			for (const Particle& p : chargedParticles)
			{
				//hPion0Pt["Pion0Pt_AuAu"]->fill(partPt, pt_weight);
			}



			switch (p.pid()) {
			case 211: // pi+
			{
				sow["sow_dAUc100"]->fill();
				hPionPosPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);

				if ((c >= 0.) && (c < 20.))
				{
					sow["sow_dAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionPosPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_dAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionPosPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_dAUc60"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionPosPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 60.) && (c < 88.))
				{
					sow["sow_dAUc88"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionPosPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
					}
				}
				break;
			}
			case -211: // pi-
			{
				sow["sow_dAUc100"]->fill();
				hPionNegPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);

				if ((c >= 0.) && (c < 20.))
				{
					sow["sow_dAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionNegPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_dAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionNegPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_dAUc60"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionNegPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 60.) && (c < 88.))
				{
					sow["sow_dAUc88"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hPionNegPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
					}
				}
				break;
			}
			case 321: // K+
			{
				sow["sow_dAUc100"]->fill();
				hKaonPosPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);

				if ((c >= 0.) && (c < 20.))
				{

					sow["sow_dAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonPosPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_dAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonPosPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_dAUc60"]->fill();

					for (const Particle& p : chargedParticles)
					{
						hKaonPosPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 60.) && (c < 88.))
				{
					sow["sow_dAUc88"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonPosPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
					}
				}
				break;
			}
			case -321: // K-
			{
				sow["sow_dAUc100"]->fill();
				hKaonNegPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);

				if ((c >= 0.) && (c < 20.))
				{
					sow["sow_dAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonNegPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_dAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonNegPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_dAUc60"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonNegPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 60.) && (c < 88.))
				{
					sow["sow_dAUc88"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hKaonNegPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
					}
				}
				break;
			}
			case 2212: // proton
			{
				sow["sow_dAUc100"]->fill();
				hProtPosPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);

				if ((c >= 0.) && (c < 20.))
				{
					sow["sow_dAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtPosPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_dAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtPosPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_dAUc60"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtPosPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 60.) && (c < 88.))
				{
					sow["sow_dAUc88"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtPosPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
					}
				}
				break;
			}
			case -2212: // anti-proton
			{
				sow["sow_dAUc100"]->fill();
				hProtNegPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);

				if ((c >= 0.) && (c < 20.))
				{
					sow["sow_dAUc20"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtNegPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_dAUc40"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtNegPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
					}
				}
				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_dAUc60"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtNegPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
					}
				}

				else if ((c >= 60.) && (c < 88.))
				{
					sow["sow_dAUc88"]->fill();
					for (const Particle& p : chargedParticles)
					{
						hProtNegPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
					}
				}
				break;
			}
			}
			

		}


		


	}

	void finalize() {

		bool AuAu200_available = false;
		bool pp200_available = false;

		//need to correct for other particles
		for (auto element : hPion0Pt)
		{
			string name = element.second->name();
			if (name.find("AuAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) AuAu200_available = true;
			}
			else if (name.find("pp") != std::string::npos)
			{
				if (element.second->numEntries() > 0) pp200_available = true;
			}
		}

		if (!(AuAu200_available && pp200_available)) return;

		//yields (fig 4)_________________ 
		for (int i = 0, N = AUAUCentralityBins.size(); i < N; ++i)
		{
			hKaonNegPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
			hPionNegPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
			hProtNegPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
			hKaonPosPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
			hPionPosPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
			hProtPosPt["ptyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
			hKaonNegPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"]->scaleW(1. / sow["sow_dAUc + std::to_string(dAUCentralityBins[i])"]->sumW());
			hPionNegPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"]->scaleW(1. / sow["sow_dAUc + std::to_string(dAUCentralityBins[i])"]->sumW());
			hProtNegPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"]->scaleW(1. / sow["sow_dAUc + std::to_string(dAUCentralityBins[i])"]->sumW());
			hKaonPosPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"]->scaleW(1. / sow["sow_dAUc + std::to_string(dAUCentralityBins[i])"]->sumW());
			hPionPosPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"]->scaleW(1. / sow["sow_dAUc + std::to_string(dAUCentralityBins[i])"]->sumW());
			hProtPosPt["ptyieldsdAuc + std::to_string(dAUCentralityBins[i])"]->scaleW(1. / sow["sow_dAUc + std::to_string(dAUCentralityBins[i])"]->sumW());
		}

		//Ratio of yields (figs 5-9)_________________(Draft)
		for (int i = 0, N = AUAUCentralityBins.size(); i < N; ++i)
		{
			hKaonRatio["RatioyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"];
			hPionRatio["RatioyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"];
			hProtRatio["RatioyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"];
			hKaonRatio["RatioyieldsdAuc + std::to_string(dAUCentralityBins[i])"];
			hPionRatio["RatioyieldsdAuc + std::to_string(dAUCentralityBins[i])"];
			hProtRatio["RatioyieldsdAuc + std::to_string(dAUCentralityBins[i])"];

			hKpiRatio["RatioyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"];
			hppiRatio["RatioyieldsAuAuc + std::to_string(AUAUCentralityBins[i])"];
			hKPiRatio["RatioyieldsdAuc + std::to_string(dAUCentralityBins[i])"];
			hppiRatio["RatioyieldsdAuc + std::to_string(dAUCentralityBins[i])"];
		}

		//RCP (fig 10) _________________(Draft and need scaling)



		hKaonPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
		hKaonPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hKaonPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_AuAu"], hKaonPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_pp"], hRcp["K_c + std::to_string(AUAUCentralityBins[i])_AuAu"]);
		
		hRcp["K_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleY(1. / 777.2);


		hPionPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
		hPionPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hPionPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_AuAu"], hPionPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_pp"], hRcp["pi_c + std::to_string(AUAUCentralityBins[i])_AuAu"]);
		
		hRcp["pi_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleY(1. / 777.2);


		hProtPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
		hProtPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hProtPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_AuAu"], hProtPt["Rcp_c + std::to_string(AUAUCentralityBins[i])_pp"], hRcp["p_c + std::to_string(AUAUCentralityBins[i])_AuAu"]);
		
		hRcp["p_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleY(1. / 777.2);

		


		//RAA (fig 11)_________________(Draft and need scaling)

		hKaonPt["Raa_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
		hKaonPt["Raa_c + std::to_string(AUAUCentralityBins[i])_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hKaonPt["Raa_c + std::to_string(AUAUCentralityBins[i])_AuAu"], hKaonPt["Raa_c + std::to_string(AUAUCentralityBins[i])_pp"], hRaa["K_c + std::to_string(AUAUCentralityBins[i])_AuAu"]);

		hRaa["K_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleY(1. / 777.2);


		hPionPt["Raa_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
		hPionPt["Raa_c + std::to_string(AUAUCentralityBins[i])_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hPionPt["Raa_c + std::to_string(AUAUCentralityBins[i])_AuAu"], hPionPt["Raa_c + std::to_string(AUAUCentralityBins[i])_pp"], hRaa["pi_c + std::to_string(AUAUCentralityBins[i])_AuAu"]);

		hRaa["pi_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleY(1. / 777.2);


		hProtPt["Raa_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleW(1. / sow["sow_AUAUc + std::to_string(AUAUCentralityBins[i]"]->sumW());
		hProtPt["Raa_c + std::to_string(AUAUCentralityBins[i])_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hProtPt["Raa_c + std::to_string(AUAUCentralityBins[i])_AuAu"], hProtPt["Raa_c + std::to_string(AUAUCentralityBins[i])_pp"], hRaa["p_c + std::to_string(AUAUCentralityBins[i])_AuAu"]);

		hRaa["p_c + std::to_string(AUAUCentralityBins[i])_AuAu"]->scaleY(1. / 777.2);


		//RdA (fig 12)_________________(Draft and need scaling)

		hKaonPt["Rda_c + std::to_string(dAUCentralityBins[i])_dAu"]->scaleW(1. / sow["sow_dAUc + std::to_string(dAUCentralityBins[i]"]->sumW());
		hKaonPt["Rda_c + std::to_string(dAUCentralityBins[i])_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hKaonPt["Rda_c + std::to_string(dAUCentralityBins[i])_dAu"], hKaonPt["Rda_c + std::to_string(dAUCentralityBins[i])_pp"], hRda["K_c + std::to_string(dAUCentralityBins[i])_dAu"]);

		hRda["K_c + std::to_string(dAUCentralityBins[i])_dAu"]->scaleY(1. / 777.2);


		hPionPt["Rda_c + std::to_string(dAUCentralityBins[i])_dAu"]->scaleW(1. / sow["sow_dAUc + std::to_string(dAUCentralityBins[i]"]->sumW());
		hPionPt["Rda_c + std::to_string(dAUCentralityBins[i])_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hPionPt["Raa_c + std::to_string(dAUCentralityBins[i])_dAu"], hPionPt["Rda_c + std::to_string(dAUCentralityBins[i])_pp"], hRda["pi_c + std::to_string(dAUCentralityBins[i])_dAu"]);

		hRda["pi_c + std::to_string(dAUCentralityBins[i])_dAu"]->scaleY(1. / 777.2);


		hProtPt["Rda_c + std::to_string(dAUCentralityBins[i])_dAu"]->scaleW(1. / sow["sow_dAUc + std::to_string(dAUCentralityBins[i]"]->sumW());
		hProtPt["Rda_c + std::to_string(dAUCentralityBins[i])_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hProtPt["Rda_c + std::to_string(dAUCentralityBins[i])_dAu"], hProtPt["Rda_c + std::to_string(dAUCentralityBins[i])_pp"], hRda["p_c + std::to_string(dAUCentralityBins[i])_dAu"]);

		hRda["p_c + std::to_string(dAUCentralityBins[i])_dAu"]->scaleY(1. / 777.2);

		// Ratio of Spectra(fig 15)_________________

	}


	map<string, Histo1DPtr> hKaonNegPt;
	map<string, Histo1DPtr> hKaonPosPt;
	map<string, Histo1DPtr> hPionNegPt;
	map<string, Histo1DPtr> hPionPosPt;
	map<string, Histo1DPtr> hProtNegPt;
	map<string, Histo1DPtr> hProtPosPt;
	map<string, Histo1DPtr> hKaonPt;
	map<string, Histo1DPtr> hPionPt;
	map<string, Histo1DPtr> hProtPt;



	map<double, Histo1DPtr> rKK_KaonNeg;
	map<double, Histo1DPtr> rKK_KaonPos;
	map<double, Scatter2DPtr> RatioKaon;

	map<double, Histo1DPtr> rKK_PionNeg;
	map<double, Histo1DPtr> rKK_PionPos;
	map<double, Scatter2DPtr> RatioPion;

	map<double, Histo1DPtr> rKK_ProtNeg;
	map<double, Histo1DPtr> rKK_ProtPos;
	map<double, Scatter2DPtr> RatioProt;

	map<double, Histo1DPtr> rKPi_Kaon;
	map<double, Histo1DPtr> rKPi_Pion;
	map<double, Scatter2DPtr> RatioKpi;

	map<double, Histo1DPtr> rKPi_Prot;
	map<double, Histo1DPtr> rKPi_Pion;
	map<double, Scatter2DPtr> Ratioppi;

	Scatter2DPtr hRcp; 
	Scatter2DPtr hRaa;
	Scatter2DPtr hRda;
	map<string, CounterPtr> sow;
	string beamOpt;
	enum CollisionSystem {pp, AuAu200, dAu200};
	CollisionSystem collSys;
	vector<int> AUAUCentralityBins{10, 20, 40, 60, 92};
	vector<int> dAUCentralityBins{20, 100, 40, 60, 88};

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2013_I1227971);

}
