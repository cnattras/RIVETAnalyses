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

			std::initializer_list<int> pdgIds = { 321, -321, 211, -211, 2212, -2212 };

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
				book(hKaonNegPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 1, 1, 1 + i);
				book(hPionNegPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 3, 1, 1 + i);
				book(hProtNegPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 5, 1, 1 + i);
				book(hKaonPosPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 1, 1, 6 + i);
				book(hPionPosPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 3, 1, 6 + i);
				book(hProtPosPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 5, 1, 6 + i);
				book(hKaonNegPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])], 2, 1, 1 + i);
				book(hPionNegPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])], 4, 1, 1 + i);
				book(hProtNegPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])], 6, 1, 1 + i);
				book(hKaonPosPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])], 2, 1, 6 + i);
				book(hPionPosPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])], 4, 1, 6 + i);
				book(hProtPosPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])], 6, 1, 6 + i);

				book(sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i]], "sow_AUAUc" + std::to_string(AUAUCentralityBins[i]));
				book(sow["sow_dAUc" + std::to_string(dAUCentralityBins[i])], "sow_dAUc" + std::to_string(dAUCentralityBins[i]));
			

				//Ratio of yields (figs 5-9)_________________
			
				//Histograms for the ratios Neg/Pos			
				string refnameK_K = mkAxisCode(7, 1, 1 + i);
					const Scatter2D& refdataK_K = refData(refnameK_K);
				book(rKK_KaonNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnameK_K + "_KaonNeg", refdataK_K);
				book(rKK_KaonPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnameK_K + "_KaonPos", refdataK_K);
				book(RatioKaon["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnameK_K);

				string refnameK_K = mkAxisCode(8, 1, 1 + i);
					const Scatter2D& refdataK_K = refData(refnameK_K);
				book(rKK_KaonNeg["dAuc" + std::to_string(dAUCentralityBins[i])], refnameK_K + "_KaonNeg", refdataK_K);
				book(rKK_KaonPos["dAuc" + std::to_string(dAUCentralityBins[i])], refnameK_K + "_KaonPos", refdataK_K);
				book(RatioKaon["dAuc" + std::to_string(dAUCentralityBins[i])], refnameK_K);
				
				string refnamepi_pi = mkAxisCode(9, 1, 1 + i);
					const Scatter2D& refdatapi_pi = refData(refnamepi_pi);
				book(rpipi_PionNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnamepi_pi + "_PionNeg", refdatapi_pi);
				book(rpipi_PionPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnamepi_pi + "_PionPos", refdatapi_pi);
				book(RatioPion["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnamepi_pi);

				string refnamepi_pi = mkAxisCode(10, 1, 1 + i);
					const Scatter2D& refdatapi_pi = refData(refnamepi_pi);
				book(rpipi_PionNeg["dAuc" + std::to_string(dAUCentralityBins[i])], refnamepi_pi + "_PionNeg", refdatapi_pi);
				book(rpipi_PionPos["dAuc" + std::to_string(dAUCentralityBins[i])], refnamepi_pi + "_PionPos", refdatapi_pi);
				book(RatioPion["dAuc" + std::to_string(dAUCentralityBins[i])], refnamepi_pi);

				string refnamep_p = mkAxisCode(11, 1, 1 + i);
					const Scatter2D& refdatap_p = refData(refnamep_p);
				book(rpp_ProtNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnamep_p + "_ProtNeg", refdatap_p);
				book(rpp_ProtPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnamep_p + "_ProtPos", refdatap_p);
				book(RatioProt["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnamep_p);

				string refnamep_p = mkAxisCode(12, 1, 1 + i);
					const Scatter2D& refdatap_p = refData(refnamep_p);
				book(rpp_ProtNeg["dAuc" + std::to_string(dAUCentralityBins[i])], refnamep_p + "_ProtNeg", refdatap_p);
				book(rpp_ProtPos["dAuc" + std::to_string(dAUCentralityBins[i])], refnamep_p + "_ProtPos", refdatap_p);
				book(RatioProt["dAuc" + std::to_string(dAUCentralityBins[i])], refnamep_p);
						
				//Histograms for the ratios Kaon/Pion 
				string refnameK_pi = mkAxisCode(13, 1, 1 + i);
					const Scatter2D& refdataK_pi = refData(refnameK_pi);
				book(rKpi_Kaon["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnameK_pi + "_kaons", refdataK_pi);
				book(rKpi_Pion["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnameK_pi + "_pions", refdataK_pi);
				book(RatioK_pi["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnameK_pi);

				string refnameK_pi = mkAxisCode(14, 1, 1 + i);
					const Scatter2D& refdataK_pi = refData(refnameK_pi);
				book(rKpi_Kaon["dAuc" + std::to_string(dAUCentralityBins[i])], refnameK_pi + "_kaons", refdataK_pi);
				book(rKpi_Pion["dAuc" + std::to_string(dAUCentralityBins[i])], refnameK_pi + "_pions", refdataK_pi);
				book(RatioK_pi["dAuc" + std::to_string(dAUCentralityBins[i])], refnameK_pi);

				//Histograms for the ratios Proton/Pion 
				string refnamep_pi = mkAxisCode(15, 1, 1 + i);
					const Scatter2D& refdatap_pi = refData(refnamep_pi);
				book(rppi_Proton["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnamep_pi + "_protons", refdatap_pi);
				book(rppi_Pion["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnamep_pi + "_pions", refdatap_pi);
				book(Ratiop_pi["AuAuc" + std::to_string(AUAUCentralityBins[i])], refnamep_pi);
				
				string refnamep_pi = mkAxisCode(16, 1, 1 + i);
					const Scatter2D& refdatap_pi = refData(refnamep_pi);
				book(rppi_Proton["dAuc" + std::to_string(dAUCentralityBins[i])], refnamep_pi + "_protons", refdatap_pi);
				book(rppi_Pion["dAuc" + std::to_string(dAUCentralityBins[i])], refnamep_pi + "_pions", refdatap_pi);
				book(Ratiop_pi["dAuc" + std::to_string(dAUCentralityBins[i])], refnamep_pi);


				//RAA (fig 11)_________________
				
				string refnameRaa = mkAxisCode(20, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refnameRaa + "_AuAu", refdataRaa);
				book(hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRaa["K_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refnameRaa);

				string refnameRaa = mkAxisCode(21, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refnameRaa + "_AuAu", refdataRaa);
				book(hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRaa["pi_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refnameRaa);

				string refnameRaa = mkAxisCode(22, 1, 1 + i);
					const Scatter2D& refdataRaa = refData(refnameRaa);
				book(hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refnameRaa + "_AuAu", refdataRaa);
				book(hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], refnameRaa + "_pp", refdataRaa);
				book(hRaa["p_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refnameRaa);
			

				//RdA (fig 12)_________________
			
				string refnameRda = mkAxisCode(23, 1, 1 + i);
					const Scatter2D& refdataRda = refData(refnameRda);
				book(hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refnameRda + "_dAu", refdataRda);
				book(hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"], refnameRda + "_pp", refdataRda);
				book(hRda["K_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refnameRda);

				string refnameRda = mkAxisCode(24, 1, 1 + i);
					const Scatter2D& refdataRda = refData(refnameRda);
				book(hPionPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refnameRda + "_dAu", refdataRda);
				book(hPionPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"], refnameRda + "_pp", refdataRda);
				book(hRda["pi_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refnameRda);

				string refnameRda = mkAxisCode(25, 1, 1 + i);
					const Scatter2D& refdataRda = refData(refnameRda);
				book(hProtPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refnameRda + "_dAu", refdataRda);
				book(hProtPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"], refnameRda + "_pp", refdataRda);
				book(hRda["p_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refnameRda);


			}

				//RCP (fig 10) _________________ Need to check if this is correct
							   
				string refname = mkAxisCode(17, 1, 1);
					const Scatter2D& refdata = refData(refname);
				book(hKaonPosPt["AuAuc0010a"], refname + "_KaonPos", refdata);
				book(hKaonPosPt["AuAuc4060"], refname + "_KaonPos", refdata);
				book(hRcp["Kpos_c00104060_AuAu"], refname);
				
				string refname = mkAxisCode(17, 1, 2);
					const Scatter2D& refdata = refData(refname);
				book(hKaonPosPt["AuAuc0010b"], refname + "_KaonPos", refdata);
				book(hKaonPosPt["AuAuc6092"], refname + "_KaonPos", refdata);
				book(hRcp["Kpos_c00106092_AuAu"], refname);
				
				string refname = mkAxisCode(17, 1, 3);
					const Scatter2D& refdata = refData(refname);
				book(hKaonNegPt["AuAuc0010a"], refname + "_KaonNeg", refdata);
				book(hKaonNegPt["AuAuc4060"], refname + "_KaonNeg", refdata);
				book(hRcp["Kneg_c00104060_AuAu"], refname);

				string refname = mkAxisCode(17, 1, 4);
					const Scatter2D& refdata = refData(refname);
				book(hKaonNegPt["AuAuc0010b"], refname + "_KaonNeg", refdata);
				book(hKaonNegPt["AuAuc6092"], refname + "_KaonNeg", refdata);
				book(hRcp["Kneg_c00106092" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refname);

				string refname = mkAxisCode(18, 1, 1);
					const Scatter2D& refdata = refData(refname);
				book(hPionPosPt["AuAuc0010a"], refname + "_PionPos", refdata);
				book(hPionPosPt["AuAuc4060"], refname + "_PionPos", refdata);
				book(hRcp["pipos_c00104060_AuAu"], refname);

				string refname = mkAxisCode(18, 1, 2);
					const Scatter2D& refdata = refData(refname);
				book(hPionPosPt["AuAuc0010b"], refname + "_PionPos", refdata);
				book(hPionPosPt["AuAuc6092"], refname + "_PionPos", refdata);
				book(hRcp["pipos_c00106092 _AuAu"], refname);

				string refname = mkAxisCode(18, 1, 3);
					const Scatter2D& refdata = refData(refname);
				book(hPionNegPt["AuAuc0010a"], refname + "_PionNeg", refdata);
				book(hPionNegPt["AuAuc4060"], refname + "_PionNeg", refdata);
				book(hRcp["pineg_c00104060_AuAu"], refname);

				string refname = mkAxisCode(18, 1, 4);
					const Scatter2D& refdata = refData(refname);
				book(hPionNegPt["AuAuc0010b"], refname + "_PionNeg", refdata);
				book(hPionNegPt["AuAuc6092"], refname + "_PionNeg", refdata);
				book(hRcp["pineg_c00106092_AuAu"], refname);

				string refname = mkAxisCode(19, 1, 1);
					const Scatter2D& refdata = refData(refname);
				book(hProtPosPt["AuAuc0010a"], refname + "_ProtPos", refdata);
				book(hProtPosPt["AuAuc4060"], refname + "_ProtPos", refdata);
				book(hRcp["ppos_c00104060_AuAu"], refname);

				string refname = mkAxisCode(19, 1, 2);
					const Scatter2D& refdata = refData(refname);
				book(hProtPosPt["AuAuc0010b"], refname + "_ProtPos", refdata);
				book(hProtPosPt["AuAuc6092"], refname + "_ProtPos", refdata);
				book(hRcp["ppos_c00106092_AuAu"], refname);

				string refname = mkAxisCode(19, 1, 3);
					const Scatter2D& refdata = refData(refname);
				book(hProtNegPt["AuAuc0010a"], refname + "_ProtNeg", refdata);
				book(hProtNegPt["AuAuc4060"], refname + "_ProtNeg", refdata);
				book(hRcp["pneg_c00104060_AuAu"], refname);

				string refname = mkAxisCode(19, 1, 4);
					const Scatter2D& refdata = refData(refname);
				book(hProtNegPt["AuAuc0010b"], refname + "_ProtNeg", refdata);
				book(hProtNegPt["AuAuc6092"], refname + "_ProtNeg", refdata);
				book(hRcp["pneg_c00106092_AuAu"], refname);


				// Ratio of Spectra(fig 15)_________________
				
				string refname = mkAxisCode(26, 1, 1);
					const Scatter2D& refdata = refData(refname);
				book(hKaonPt["AuAuc6092"], refname + "_Kaon", refdata);
				book(hKaonPt["dAuc0020"], refname + "_Kaon", refdata);
				book(RatioK["AuAuc/dAU"], refname);

				string refname = mkAxisCode(27, 1, 1);
					const Scatter2D& refdata = refData(refname);
				book(hPionPt["AuAuc6092"], refname + "_Pion", refdata);
				book(hPionPt["dAuc0020"], refname + "_Pion", refdata);
				book(Ratiopi["AuAuc/dAU"], refname);

				string refname = mkAxisCode(28, 1, 1);
					const Scatter2D& refdata = refData(refname);
				book(hProtPt["AuAuc6092"], refname + "_Prot", refdata);
				book(hProtPt["dAuc0020"], refname + "_Prot", refdata);
				book(Ratiop["AuAuc/dAU"], refname);

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
		bool dAu200_available = false;
		bool pp_available = false;

		//check if this is correct
		for (auto element : hKaonNegPt)
		{
			string name = element.second->name();
			if (name.find("AuAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) AuAu200_available = true;
			}
			else if (name.find("dAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) dAu_available = true;
			}
			else if (name.find("pp") != std::string::npos)
			{
				if (element.second->numEntries() > 0) pp_available = true;
			}
		}
		for (auto element : hKaonPosPt)
		{
			string name = element.second->name();
			if (name.find("AuAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) AuAu200_available = true;
			}
			else if (name.find("dAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) dAu_available = true;
			}
			else if (name.find("pp") != std::string::npos)
			{
				if (element.second->numEntries() > 0) pp_available = true;
			}
		}
		for (auto element : hPionNegPt)
		{
			string name = element.second->name();
			if (name.find("AuAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) AuAu200_available = true;
			}
			else if (name.find("dAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) dAu_available = true;
			}
			else if (name.find("pp") != std::string::npos)
			{
				if (element.second->numEntries() > 0) pp_available = true;
			}
		}
		for (auto element : hPionPosPt)
		{
			string name = element.second->name();
			if (name.find("AuAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) AuAu200_available = true;
			}
			else if (name.find("dAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) dAu_available = true;
			}
			else if (name.find("pp") != std::string::npos)
			{
				if (element.second->numEntries() > 0) pp_available = true;
			}
		}
		for (auto element : hProtNegPt)
		{
			string name = element.second->name();
			if (name.find("AuAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) AuAu200_available = true;
			}
			else if (name.find("dAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) dAu_available = true;
			}
			else if (name.find("pp") != std::string::npos)
			{
				if (element.second->numEntries() > 0) pp_available = true;
			}
		}
		for (auto element : hProtPosPt)
		{
			string name = element.second->name();
			if (name.find("AuAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) AuAu200_available = true;
			}
			else if (name.find("dAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) dAu_available = true;
			}
			else if (name.find("pp") != std::string::npos)
			{
				if (element.second->numEntries() > 0) pp_available = true;
			}
		}
		for (auto element : hKaonPt)
		{
			string name = element.second->name();
			if (name.find("AuAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) AuAu200_available = true;
			}
			else if (name.find("dAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) dAu_available = true;
			}
			else if (name.find("pp") != std::string::npos)
			{
				if (element.second->numEntries() > 0) pp_available = true;
			}
		}
		for (auto element : hPionPt)
		{
			string name = element.second->name();
			if (name.find("AuAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) AuAu200_available = true;
			}
			else if (name.find("dAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) dAu_available = true;
			}
			else if (name.find("pp") != std::string::npos)
			{
				if (element.second->numEntries() > 0) pp_available = true;
			}
		}
		for (auto element : hProtPt)
		{
			string name = element.second->name();
			if (name.find("AuAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) AuAu200_available = true;
			}
			else if (name.find("dAu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) dAu_available = true;
			}
			else if (name.find("pp") != std::string::npos)
			{
				if (element.second->numEntries() > 0) pp_available = true;
			}
		}

		if (!(AuAu200_available && dAu200_available && pp_available)) return;

		
		for (int i = 0, N = AUAUCentralityBins.size(); i < N; ++i)
		{
			//yields (fig 4)_________________ 
			hKaonNegPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i]]->sumW());
			hPionNegPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i]]->sumW());
			hProtNegPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i]]->sumW());
			hKaonPosPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i]]->sumW());
			hPionPosPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i]]->sumW());
			hProtPosPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i]]->sumW());
			hKaonNegPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])]->scaleW(1. / sow["sow_dAUc" + std::to_string(dAUCentralityBins[i])]->sumW());
			hPionNegPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])]->scaleW(1. / sow["sow_dAUc" + std::to_string(dAUCentralityBins[i])]->sumW());
			hProtNegPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])]->scaleW(1. / sow["sow_dAUc" + std::to_string(dAUCentralityBins[i])]->sumW());
			hKaonPosPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])]->scaleW(1. / sow["sow_dAUc" + std::to_string(dAUCentralityBins[i])]->sumW());
			hPionPosPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])]->scaleW(1. / sow["sow_dAUc" + std::to_string(dAUCentralityBins[i])]->sumW());
			hProtPosPt["ptyieldsdAuc" + std::to_string(dAUCentralityBins[i])]->scaleW(1. / sow["sow_dAUc" + std::to_string(dAUCentralityBins[i])]->sumW());

			//Ratio of yields (figs 5-9)_________________(Do I need to scale Neg and Pos particles first (like yields section)?)
			divide(rKK_KaonNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], rKK_KaonPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], RatioKaon["AuAuc" + std::to_string(AUAUCentralityBins[i])]);
			divide(rKK_KaonNeg["dAuc" + std::to_string(dAUCentralityBins[i])], rKK_KaonPos["dAuc" + std::to_string(dAUCentralityBins[i])], RatioKaon["dAuc" + std::to_string(dAUCentralityBins[i])]);

			divide(rpipi_PionNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], rpipi_PionPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], RatioPion["AuAuc" + std::to_string(AUAUCentralityBins[i])]);
			divide(rpipi_PionNeg["dAuc" + std::to_string(dAUCentralityBins[i])], rpipi_PionPos["dAuc" + std::to_string(dAUCentralityBins[i])], RatioPion["dAuc" + std::to_string(dAUCentralityBins[i])]);

			divide(rpp_ProtNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], rpp_ProtPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], RatioProt["AuAuc" + std::to_string(AUAUCentralityBins[i])]);
			divide(rpp_ProtNeg["dAuc" + std::to_string(dAUCentralityBins[i])], rpp_ProtPos["dAuc" + std::to_string(dAUCentralityBins[i])], RatioProt["dAuc" + std::to_string(dAUCentralityBins[i])]);

			divide(rKpi_Kaon["AuAuc" + std::to_string(AUAUCentralityBins[i])], rKpi_Pion["AuAuc" + std::to_string(AUAUCentralityBins[i])], RatioK_pi["AuAuc" + std::to_string(AUAUCentralityBins[i])]);
			divide(rKpi_Kaon["dAuc" + std::to_string(dAUCentralityBins[i])], rKpi_Pion["dAuc" + std::to_string(dAUCentralityBins[i])], RatioK_pi["dAuc" + std::to_string(dAUCentralityBins[i])]);

			divide(rppi_Proton["AuAuc" + std::to_string(AUAUCentralityBins[i])], rppi_Pion["AuAuc" + std::to_string(AUAUCentralityBins[i])], Ratiop_pi["AuAuc" + std::to_string(AUAUCentralityBins[i])]);
			divide(rppi_Proton["dAuc" + std::to_string(dAUCentralityBins[i])], rppi_Pion["dAuc" + std::to_string(dAUCentralityBins[i])], Ratiop_pi["dAuc" + std::to_string(dAUCentralityBins[i])]);

			//RAA (fig 11)_________________(need scaling (will be out of for loop))		
			hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i]]->sumW());
			hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
			divide(hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], hRaa["K_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]);

			hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i]]->sumW());
			hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
			divide(hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], hRaa["pi_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]);

			hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i]]->sumW());
			hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
			divide(hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], hRaa["p_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]);

			//RdA (fig 12)_________________(need scaling (will be out of for loop))

			hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]->scaleW(1. / sow["sow_dAUc " + std::to_string(dAUCentralityBins[i]]->sumW());
			hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
			divide(hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"], hRda["K_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]);

			hPionPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]->scaleW(1. / sow["sow_dAUc " + std::to_string(dAUCentralityBins[i]]->sumW());
			hPionPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
			divide(hPionPt["Raa_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], hPionPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"], hRda["pi_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]);

			hProtPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]->scaleW(1. / sow["sow_dAUc " + std::to_string(dAUCentralityBins[i]]->sumW());
			hProtPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
			divide(hProtPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], hProtPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"], hRda["p_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]);

		}


		//RCP (fig 10) _________________(is ->scaleY needed?)

		hKaonPosPt["AuAuc0010a"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hKaonPosPt["AuAuc4060"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hKaonPosPt["AuAuc0010a"], hKaonPosPt["AuAuc4060"], hRcp["Kpos_c00104060_AuAu"]);

		hKaonPosPt["AuAuc0010b"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hKaonPosPt["AuAuc6092"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hKaonPosPt["AuAuc0010b"], hKaonPosPt["AuAuc6092"], hRcp["Kpos_c00106092_AuAu"]);

		hKaonNegPt["AuAuc0010a"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hKaonNegPt["AuAuc4060"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hKaonNegPt["AuAuc0010a"], hKaonNegPt["AuAuc4060"], hRcp["Kneg_c00104060_AuAu"]);

		hKaonNegPt["AuAuc0010b"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hKaonNegPt["AuAuc6092"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hKaonNegPt["AuAuc0010b"], hKaonNegPt["AuAuc6092"], hRcp["Kneg_c00106092_AuAu"]);



		hPionPosPt["AuAuc0010a"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hPionPosPt["AuAuc4060"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hPionPosPt["AuAuc0010a"], hPionPosPt["AuAuc4060"], hRcp["pipos_c00104060_AuAu"]);

		hPionPosPt["AuAuc0010b"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hPionPosPt["AuAuc6092"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hPionPosPt["AuAuc0010b"], hPionPosPt["AuAuc6092"], hRcp["pipos_c00106092_AuAu"]);

		hPionNegPt["AuAuc0010a"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hPionNegPt["AuAuc4060"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hPionNegPt["AuAuc0010a"], hPionNegPt["AuAuc4060"], hRcp["pineg_c00104060_AuAu"]);

		hPionNegPt["AuAuc0010b"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hPionNegPt["AuAuc6092"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hPionNegPt["AuAuc0010b"], hPionNegPt["AuAuc6092"], hRcp["pineg_c00106092_AuAu"]);



		hProtPosPt["AuAuc0010a"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hProtPosPt["AuAuc4060"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hProtPosPt["AuAuc0010a"], hProtPosPt["AuAuc4060"], hRcp["ppos_c00104060_AuAu"]);

		hProtPosPt["AuAuc0010b"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hProtPosPt["AuAuc6092"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hProtPosPt["AuAuc0010b"], hProtPosPt["AuAuc6092"], hRcp["ppos_c00106092_AuAu"]);

		hProtNegPt["AuAuc0010a"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hProtNegPt["AuAuc4060"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hProtNegPt["AuAuc0010a"], hProtNegPt["AuAuc4060"], hRcp["pneg_c00104060_AuAu"]);

		hProtNegPt["AuAuc0010b"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
		hProtNegPt["AuAuc6092"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
		divide(hProtNegPt["AuAuc0010b"], hProtNegPt["AuAuc6092"], hRcp["pneg_c00106092_AuAu"]);





		// Ratio of Spectra(fig 15)_________________(Do I treat this the same as ratio of yields section?)

		divide(hKaonPt["AuAuc6092"], hKaonPt["dAuc0020"], RatioK["AuAuc/dAU"]);
		divide(hPionPt["AuAuc6092"], hPionPt["dAuc0020"], Ratiopi["AuAuc/dAU"]);
		divide(hProtPt["AuAuc6092"], hProtPt["dAuc0020"], Ratiop["AuAuc/dAU"]);


		//scaling for Raa and Rda
		hRaa["K_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]->scaleY(1. / 777.2);
		hRaa["pi_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]->scaleY(1. / 777.2);
		hRaa["p_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]->scaleY(1. / 777.2);

		hRda["K_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]->scaleY(1. / 777.2);
		hRda["pi_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]->scaleY(1. / 777.2);
		hRda["p_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]->scaleY(1. / 777.2);


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

	map<double, Histo1DPtr> rpipi_PionNeg;
	map<double, Histo1DPtr> rpipi_PionPos;
	map<double, Scatter2DPtr> RatioPion;

	map<double, Histo1DPtr> rpp_ProtNeg;
	map<double, Histo1DPtr> rpp_ProtPos;
	map<double, Scatter2DPtr> RatioProt;

	map<double, Histo1DPtr> rKPi_Kaon;
	map<double, Histo1DPtr> rKPi_Pion;
	map<double, Scatter2DPtr> RatioKpi;

	map<double, Histo1DPtr> rKPi_Prot;
	map<double, Histo1DPtr> rKPi_Pion;
	map<double, Scatter2DPtr> Ratioppi;

	map<double, Scatter2DPtr> RatioK;
	map<double, Scatter2DPtr> Ratiopi;
	map<double, Scatter2DPtr> Ratiop;

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
