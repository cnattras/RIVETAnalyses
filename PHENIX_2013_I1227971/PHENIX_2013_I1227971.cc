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
#include "../Centralities/RHICCentrality.hh" //external header for Centrality calculation
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

				const PrimaryParticles fs(pdgIds, Cuts::abseta < 0.35 && Cuts::abscharge > 0);
				declare(fs, "fs");

				beamOpt = getOption<string>("beam", "NONE");

				if (beamOpt == "PP") collSys = pp;
				else if (beamOpt == "AUAU200") collSys = AuAu200;
				else if (beamOpt == "dAU200") collSys = dAu200;


				if (!(collSys == pp)) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

				book(sow["sow_pp"], "_sow_pp");

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

					book(sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])], "_sow_AUAUc" + std::to_string(AUAUCentralityBins[i]));
					book(sow["sow_dAUc" + std::to_string(dAUCentralityBins[i])], "_sow_dAUc" + std::to_string(dAUCentralityBins[i]));


					//Ratio of yields (figs 5-9)_________________

					//Histograms for the ratios Neg/Pos
					string refname1 = mkAxisCode(7, 1, 1 + i);
					const Scatter2D& refdata1 = refData(refname1);
					book(rKK_KaonNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname1 + "_KaonNeg", refdata1);
					book(rKK_KaonPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname1 + "_KaonPos", refdata1);
					book(RatioKaon["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname1);

					string refname2 = mkAxisCode(8, 1, 1 + i);
					const Scatter2D& refdata2 = refData(refname2);
					book(rKK_KaonNeg["dAuc" + std::to_string(dAUCentralityBins[i])], refname2 + "_KaonNeg", refdata2);
					book(rKK_KaonPos["dAuc" + std::to_string(dAUCentralityBins[i])], refname2 + "_KaonPos", refdata2);
					book(RatioKaon["dAuc" + std::to_string(dAUCentralityBins[i])], refname2);

					string refname3 = mkAxisCode(9, 1, 1 + i);
					const Scatter2D& refdata3 = refData(refname3);
					book(rpipi_PionNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname3 + "_PionNeg", refdata3);
					book(rpipi_PionPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname3 + "_PionPos", refdata3);
					book(RatioPion["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname3);

					string refname4 = mkAxisCode(10, 1, 1 + i);
					const Scatter2D& refdata4 = refData(refname4);
					book(rpipi_PionNeg["dAuc" + std::to_string(dAUCentralityBins[i])], refname4 + "_PionNeg", refdata4);
					book(rpipi_PionPos["dAuc" + std::to_string(dAUCentralityBins[i])], refname4 + "_PionPos", refdata4);
					book(RatioPion["dAuc" + std::to_string(dAUCentralityBins[i])], refname4);

					string refname5 = mkAxisCode(11, 1, 1 + i);
					const Scatter2D& refdata5 = refData(refname5);
					book(rpp_ProtNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname5 + "_ProtNeg", refdata5);
					book(rpp_ProtPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname5 + "_ProtPos", refdata5);
					book(RatioProt["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname5);

					string refname6 = mkAxisCode(12, 1, 1 + i);
					const Scatter2D& refdata6 = refData(refname6);
					book(rpp_ProtNeg["dAuc" + std::to_string(dAUCentralityBins[i])], refname6 + "_ProtNeg", refdata6);
					book(rpp_ProtPos["dAuc" + std::to_string(dAUCentralityBins[i])], refname6 + "_ProtPos", refdata6);
					book(RatioProt["dAuc" + std::to_string(dAUCentralityBins[i])], refname6);

					//Histograms for the ratios Kaon/Pion
					string refname7 = mkAxisCode(13, 1, 1 + i);
					const Scatter2D& refdata7 = refData(refname7);
					book(rKpi_KaonPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname7 + "_kaons", refdata7);
					book(rKpi_PionPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname7 + "_pions", refdata7);
					book(RatioK_pipos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname7);

					string refname8 = mkAxisCode(13, 1, 6 + i);
					const Scatter2D& refdata8 = refData(refname8);
					book(rKpi_KaonNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname8 + "_kaons", refdata8);
					book(rKpi_PionNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname8 + "_pions", refdata8);
					book(RatioK_pineg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname8);

					string refname9 = mkAxisCode(14, 1, 1 + i);
					const Scatter2D& refdata9 = refData(refname9);
					book(rKpi_KaonPos["dAuc" + std::to_string(dAUCentralityBins[i])], refname9 + "_kaons", refdata9);
					book(rKpi_PionPos["dAuc" + std::to_string(dAUCentralityBins[i])], refname9 + "_pions", refdata9);
					book(RatioK_pipos["dAuc" + std::to_string(dAUCentralityBins[i])], refname9);

					string refname10 = mkAxisCode(14, 1, 6 + i);
					const Scatter2D& refdata10 = refData(refname10);
					book(rKpi_KaonNeg["dAuc" + std::to_string(dAUCentralityBins[i])], refname10 + "_kaons", refdata10);
					book(rKpi_PionNeg["dAuc" + std::to_string(dAUCentralityBins[i])], refname10 + "_pions", refdata10);
					book(RatioK_pineg["dAuc" + std::to_string(dAUCentralityBins[i])], refname10);

					//Histograms for the ratios Proton/Pion
					string refname11 = mkAxisCode(15, 1, 1 + i);
					const Scatter2D& refdata11 = refData(refname11);
					book(rppi_ProtonPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname11 + "_protons", refdata11);
					book(rppi_PionPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname11 + "_pions", refdata11);
					book(Ratiop_pipos["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname11);

					string refname12 = mkAxisCode(15, 1, 6 + i);
					const Scatter2D& refdata12 = refData(refname12);
					book(rppi_ProtonNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname12 + "_protons", refdata12);
					book(rppi_PionNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname12 + "_pions", refdata12);
					book(Ratiop_pineg["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname12);

					string refname13 = mkAxisCode(16, 1, 1 + i);
					const Scatter2D& refdata13 = refData(refname13);
					book(rppi_ProtonPos["dAuc" + std::to_string(dAUCentralityBins[i])], refname13 + "_protons", refdata13);
					book(rppi_PionPos["dAuc" + std::to_string(dAUCentralityBins[i])], refname13 + "_pions", refdata13);
					book(Ratiop_pipos["dAuc" + std::to_string(dAUCentralityBins[i])], refname13);

					string refname14 = mkAxisCode(16, 1, 6 + i);
					const Scatter2D& refdata14 = refData(refname14);
					book(rppi_ProtonNeg["dAuc" + std::to_string(dAUCentralityBins[i])], refname14 + "_protons", refdata14);
					book(rppi_PionNeg["dAuc" + std::to_string(dAUCentralityBins[i])], refname14 + "_pions", refdata14);
					book(Ratiop_pineg["dAuc" + std::to_string(dAUCentralityBins[i])], refname14);


					//RAA (fig 11)_________________

					string refname15 = mkAxisCode(20, 1, 1 + i);
					const Scatter2D& refdata15 = refData(refname15);
					book(hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refname15 + "_AuAu", refdata15);
					book(hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], refname15 + "_pp", refdata15);
					book(hRaa["K_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refname15);

					string refname16 = mkAxisCode(21, 1, 1 + i);
					const Scatter2D& refdata16 = refData(refname16);
					book(hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refname16 + "_AuAu", refdata16);
					book(hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], refname16 + "_pp", refdata16);
					book(hRaa["pi_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refname16);

					string refname17 = mkAxisCode(22, 1, 1 + i);
					const Scatter2D& refdata17 = refData(refname17);
					book(hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refname17 + "_AuAu", refdata17);
					book(hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], refname17 + "_pp", refdata17);
					book(hRaa["p_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], refname17);


					//RdA (fig 12)_________________

					string refname18 = mkAxisCode(23, 1, 1 + i);
					const Scatter2D& refdata18 = refData(refname18);
					book(hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refname18 + "_dAu", refdata18);
					book(hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"], refname18 + "_pp", refdata18);
					book(hRda["K_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refname18);

					string refname19 = mkAxisCode(24, 1, 1 + i);
					const Scatter2D& refdata19 = refData(refname19);
					book(hPionPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refname19 + "_dAu", refdata19);
					book(hPionPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"], refname19 + "_pp", refdata19);
					book(hRda["pi_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refname19);

					string refname20 = mkAxisCode(25, 1, 1 + i);
					const Scatter2D& refdata20 = refData(refname20);
					book(hProtPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refname20 + "_dAu", refdata20);
					book(hProtPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"], refname20 + "_pp", refdata20);
					book(hRda["p_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], refname20);


				}

				//RCP (fig 10) _________________ Need to check if this is correct

				string refname21 = mkAxisCode(17, 1, 1);
				const Scatter2D& refdata21 = refData(refname21);
				book(hKaonPosPt["AuAuc0010a"], refname21 + "_AuAuc0010a_KaonPos", refdata21);
				book(hKaonPosPt["AuAuc4060"], refname21 + "_AuAuc4060_KaonPos", refdata21);
				book(hRcp["Kpos_c00104060_AuAu"], refname21);

				string refname22 = mkAxisCode(17, 1, 2);
				const Scatter2D& refdata22 = refData(refname22);
				book(hKaonPosPt["AuAuc0010b"], refname22 + "_AuAuc0010b_KaonPos", refdata22);
				book(hKaonPosPt["AuAuc6092"], refname22 + "_AuAuc6092_KaonPos", refdata22);
				book(hRcp["Kpos_c00106092_AuAu"], refname22);

				string refname23 = mkAxisCode(17, 1, 3);
				const Scatter2D& refdata23 = refData(refname23);
				book(hKaonNegPt["AuAuc0010a"], refname23 + "_AuAuc0010a_KaonNeg", refdata23);
				book(hKaonNegPt["AuAuc4060"], refname23 + "_AuAuc4060_KaonNeg", refdata23);
				book(hRcp["Kneg_c00104060_AuAu"], refname23);

				string refname24 = mkAxisCode(17, 1, 4);
				const Scatter2D& refdata24 = refData(refname24);
				book(hKaonNegPt["AuAuc0010b"], refname24 + "_AuAuc0010b_KaonNeg", refdata24);
				book(hKaonNegPt["AuAuc6092"], refname24 + "_AuAuc6092_KaonNeg", refdata24);
				book(hRcp["Kneg_c00106092_AuAu"], refname24);

				string refname25 = mkAxisCode(18, 1, 1);
				const Scatter2D& refdata25 = refData(refname25);
				book(hPionPosPt["AuAuc0010a"], refname25 + "_AuAuc0010a_PionPos", refdata25);
				book(hPionPosPt["AuAuc4060"], refname25 + "_AuAuc4060_PionPos", refdata25);
				book(hRcp["pipos_c00104060_AuAu"], refname25);

				string refname26 = mkAxisCode(18, 1, 2);
				const Scatter2D& refdata26 = refData(refname26);
				book(hPionPosPt["AuAuc0010b"], refname26 + "_AuAuc0010b_PionPos", refdata26);
				book(hPionPosPt["AuAuc6092"], refname26 + "_AuAuc6092_PionPos", refdata26);
				book(hRcp["pipos_c00106092_AuAu"], refname26);

				string refname27 = mkAxisCode(18, 1, 3);
				const Scatter2D& refdata27 = refData(refname27);
				book(hPionNegPt["AuAuc0010a"], refname27 + "_AuAuc0010a_PionNeg", refdata27);
				book(hPionNegPt["AuAuc4060"], refname27 + "_AuAuc4060_PionNeg", refdata27);
				book(hRcp["pineg_c00104060_AuAu"], refname27);

				string refname28 = mkAxisCode(18, 1, 4);
				const Scatter2D& refdata28 = refData(refname28);
				book(hPionNegPt["AuAuc0010b"], refname28 + "_AuAuc0010b_PionNeg", refdata28);
				book(hPionNegPt["AuAuc6092"], refname28 + "_AuAuc6092_PionNeg", refdata28);
				book(hRcp["pineg_c00106092_AuAu"], refname28);

				string refname29 = mkAxisCode(19, 1, 1);
				const Scatter2D& refdata29 = refData(refname29);
				book(hProtPosPt["AuAuc0010a"], refname29 + "_AuAuc0010a_ProtPos", refdata29);
				book(hProtPosPt["AuAuc4060"], refname29 + "_AuAuc4060_ProtPos", refdata29);
				book(hRcp["ppos_c00104060_AuAu"], refname29);

				string refname30 = mkAxisCode(19, 1, 2);
				const Scatter2D& refdata30 = refData(refname30);
				book(hProtPosPt["AuAuc0010b"], refname30 + "_AuAuc0010b_ProtPos", refdata30);
				book(hProtPosPt["AuAuc6092"], refname30 + "_AuAuc6092_ProtPos", refdata30);
				book(hRcp["ppos_c00106092_AuAu"], refname30);

				string refname31 = mkAxisCode(19, 1, 3);
				const Scatter2D& refdata31 = refData(refname31);
				book(hProtNegPt["AuAuc0010a"], refname31 + "_AuAuc0010a_ProtNeg", refdata31);
				book(hProtNegPt["AuAuc4060"], refname31 + "_AuAuc4060_ProtNeg", refdata31);
				book(hRcp["pneg_c00104060_AuAu"], refname31);

				string refname32 = mkAxisCode(19, 1, 4);
				const Scatter2D& refdata32 = refData(refname32);
				book(hProtNegPt["AuAuc0010b"], refname32 + "_AuAuc0010b_ProtNeg", refdata32);
				book(hProtNegPt["AuAuc6092"], refname32 + "_AuAuc6092_ProtNeg", refdata32);
				book(hRcp["pneg_c00106092_AuAu"], refname32);


				// Ratio of Spectra(fig 15)_________________

				string refname33 = mkAxisCode(26, 1, 1);
				const Scatter2D& refdata33 = refData(refname33);
				book(hKaonPt["AuAuc6092"], refname33 + "_AuAuc6092_Kaon", refdata33);
				book(hKaonPt["dAuc0020"], refname33 + "_dAuc0020_Kaon", refdata33);
				book(RatioK["AuAuc/dAU"], refname33);

				string refname34 = mkAxisCode(27, 1, 1);
				const Scatter2D& refdata34 = refData(refname34);
				book(hPionPt["AuAuc6092"], refname34 + "_AuAuc6092_Pion", refdata34);
				book(hPionPt["dAuc0020"], refname34 + "_dAuc0020_Pion", refdata34);
				book(Ratiopi["AuAuc/dAU"], refname34);

				string refname35 = mkAxisCode(28, 1, 1);
				const Scatter2D& refdata35 = refData(refname35);
				book(hProtPt["AuAuc6092"], refname35 + "_AuAuc6092_Prot", refdata35);
				book(hProtPt["dAuc0020"], refname35 + "_dAuc0020_Prot", refdata35);
				book(Ratiop["AuAuc/dAU"], refname35);

			}


			void analyze(const Event& event) {
				Particles chargedParticles = applyProjection<PrimaryParticles>(event, "fs").particles();

				if (collSys == pp)
				{
					sow["sow_pp"]->fill();
					for (Particle p : chargedParticles)
					{
						double partPt = p.pT() / GeV;

						switch (p.pid()) {
						case 211: // pi+
						{
							hPionPt["Raa_c10_pp"]->fill(partPt);
							hPionPt["Raa_c20_pp"]->fill(partPt);
							hPionPt["Raa_c40_pp"]->fill(partPt);
							hPionPt["Raa_c60_pp"]->fill(partPt);
							hPionPt["Raa_c92_pp"]->fill(partPt);

							hPionPt["Rda_c100_pp"]->fill(partPt);
							hPionPt["Rda_c20_pp"]->fill(partPt);
							hPionPt["Rda_c40_pp"]->fill(partPt);
							hPionPt["Rda_c60_pp"]->fill(partPt);
							hPionPt["Rda_c88_pp"]->fill(partPt);
							break;
						}
						case -211: // pi-
						{
							hPionPt["Raa_c10_pp"]->fill(partPt);
							hPionPt["Raa_c20_pp"]->fill(partPt);
							hPionPt["Raa_c40_pp"]->fill(partPt);
							hPionPt["Raa_c60_pp"]->fill(partPt);
							hPionPt["Raa_c92_pp"]->fill(partPt);

							hPionPt["Rda_c100_pp"]->fill(partPt);
							hPionPt["Rda_c20_pp"]->fill(partPt);
							hPionPt["Rda_c40_pp"]->fill(partPt);
							hPionPt["Rda_c60_pp"]->fill(partPt);
							hPionPt["Rda_c88_pp"]->fill(partPt);
							break;
						}
						case  321: // K+
						{
							hKaonPt["Raa_c10_pp"]->fill(partPt);
							hKaonPt["Raa_c20_pp"]->fill(partPt);
							hKaonPt["Raa_c40_pp"]->fill(partPt);
							hKaonPt["Raa_c60_pp"]->fill(partPt);
							hKaonPt["Raa_c92_pp"]->fill(partPt);

							hKaonPt["Rda_c100_pp"]->fill(partPt);
							hKaonPt["Rda_c20_pp"]->fill(partPt);
							hKaonPt["Rda_c40_pp"]->fill(partPt);
							hKaonPt["Rda_c60_pp"]->fill(partPt);
							hKaonPt["Rda_c88_pp"]->fill(partPt);

							break;
						}
						case  -321: // K-
						{
							hKaonPt["Raa_c10_pp"]->fill(partPt);
							hKaonPt["Raa_c20_pp"]->fill(partPt);
							hKaonPt["Raa_c40_pp"]->fill(partPt);
							hKaonPt["Raa_c60_pp"]->fill(partPt);
							hKaonPt["Raa_c92_pp"]->fill(partPt);

							hKaonPt["Rda_c100_pp"]->fill(partPt);
							hKaonPt["Rda_c20_pp"]->fill(partPt);
							hKaonPt["Rda_c40_pp"]->fill(partPt);
							hKaonPt["Rda_c60_pp"]->fill(partPt);
							hKaonPt["Rda_c88_pp"]->fill(partPt);

							break;
						}
						case 2212: // proton
						{
							hProtPt["Raa_c10_pp"]->fill(partPt);
							hProtPt["Raa_c20_pp"]->fill(partPt);
							hProtPt["Raa_c40_pp"]->fill(partPt);
							hProtPt["Raa_c60_pp"]->fill(partPt);
							hProtPt["Raa_c92_pp"]->fill(partPt);

							hProtPt["Rda_c100_pp"]->fill(partPt);
							hProtPt["Rda_c20_pp"]->fill(partPt);
							hProtPt["Rda_c40_pp"]->fill(partPt);
							hProtPt["Rda_c60_pp"]->fill(partPt);
							hProtPt["Rda_c88_pp"]->fill(partPt);
							break;
						}
						case -2212: // anti-proton
						{
							hProtPt["Raa_c10_pp"]->fill(partPt);
							hProtPt["Raa_c20_pp"]->fill(partPt);
							hProtPt["Raa_c40_pp"]->fill(partPt);
							hProtPt["Raa_c60_pp"]->fill(partPt);
							hProtPt["Raa_c92_pp"]->fill(partPt);

							hProtPt["Rda_c100_pp"]->fill(partPt);
							hProtPt["Rda_c20_pp"]->fill(partPt);
							hProtPt["Rda_c40_pp"]->fill(partPt);
							hProtPt["Rda_c60_pp"]->fill(partPt);
							hProtPt["Rda_c88_pp"]->fill(partPt);
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



					if ((c >= 0.) && (c < 10.))
					{
						sow["sow_AUAUc10"]->fill();
						for (const Particle& p : chargedParticles)
						{
							double partPt = p.pT() / GeV;
							double pt_weight = 1. / (partPt * 2. * M_PI);

							switch (p.pid()) {
							case 211: // pi+
							{
								hPionPosPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
								rpipi_PionPos["AuAuc10"]->fill(partPt);
								rKpi_PionPos["AuAuc10"]->fill(partPt);
								rppi_PionPos["AuAuc10"]->fill(partPt);
								hPionPt["Raa_c10_AuAu"]->fill(partPt);
								hPionPosPt["AuAuc0010a"]->fill(partPt);
								hPionPosPt["AuAuc0010b"]->fill(partPt);
								break;
							}
							case -211: // pi-
							{
								hPionNegPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
								rpipi_PionNeg["AuAuc10"]->fill(partPt);
								rKpi_PionNeg["AuAuc10"]->fill(partPt);
								rppi_PionNeg["AuAuc10"]->fill(partPt);
								hPionPt["Raa_c10_AuAu"]->fill(partPt);
								hPionNegPt["AuAuc0010a"]->fill(partPt);
								hPionNegPt["AuAuc0010b"]->fill(partPt);
								break;
							}

							case 321: // K+
							{
								hKaonPosPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
								rKK_KaonPos["AuAuc10"]->fill(partPt);
								rKpi_KaonPos["AuAuc10"]->fill(partPt);
								hKaonPt["Raa_c10_AuAu"]->fill(partPt);
								hKaonPosPt["AuAuc0010a"]->fill(partPt);
								hKaonPosPt["AuAuc0010b"]->fill(partPt);
								break;
							}

							case -321: // K-
							{
								hKaonNegPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
								rKK_KaonNeg["AuAuc10"]->fill(partPt);
								rKpi_KaonNeg["AuAuc10"]->fill(partPt);
								hKaonPt["Raa_c10_AuAu"]->fill(partPt);
								hKaonNegPt["AuAuc0010a"]->fill(partPt);
								hKaonNegPt["AuAuc0010b"]->fill(partPt);
								break;
							}

							case 2212: // proton
							{
								hProtPosPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
								rpp_ProtPos["AuAuc10"]->fill(partPt);
								rppi_ProtonPos["AuAuc10"]->fill(partPt);
								hProtPt["Raa_c10_AuAu"]->fill(partPt);
								hProtPosPt["AuAuc0010a"]->fill(partPt);
								hProtPosPt["AuAuc0010b"]->fill(partPt);
								break;
							}

							case -2212: // anti-proton
							{
								hProtNegPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
								rpp_ProtNeg["AuAuc10"]->fill(partPt);
								rppi_ProtonNeg["AuAuc10"]->fill(partPt);
								hProtPt["Raa_c10_AuAu"]->fill(partPt);
								hProtNegPt["AuAuc0010a"]->fill(partPt);
								hProtNegPt["AuAuc0010b"]->fill(partPt);
								break;
							}
							}
						}
					}
					else if ((c >= 10.) && (c < 20.))
					{
						sow["sow_AUAUc20"]->fill();
						for (const Particle& p : chargedParticles)
						{
							double partPt = p.pT() / GeV;
							double pt_weight = 1. / (partPt * 2. * M_PI);

							switch (p.pid()) {
							case 211: // pi+
							{
								hPionPosPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
								rpipi_PionPos["AuAuc20"]->fill(partPt);
								rKpi_PionPos["AuAuc20"]->fill(partPt);
								rppi_PionPos["AuAuc20"]->fill(partPt);
								hPionPt["Raa_c20_AuAu"]->fill(partPt);
								break;
							}
							case -211: // pi-
							{
								hPionNegPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
								rpipi_PionNeg["AuAuc20"]->fill(partPt);
								rKpi_PionNeg["AuAuc20"]->fill(partPt);
								rppi_PionNeg["AuAuc20"]->fill(partPt);
								hPionPt["Raa_c20_AuAu"]->fill(partPt);
								break;
							}

							case 321: // K+
							{
								hKaonPosPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
								rKK_KaonPos["AuAuc20"]->fill(partPt);
								rKpi_KaonPos["AuAuc20"]->fill(partPt);
								hKaonPt["Raa_c20_AuAu"]->fill(partPt);
								break;
							}

							case -321: // K-
							{
								hKaonNegPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
								rKK_KaonNeg["AuAuc20"]->fill(partPt);
								rKpi_KaonNeg["AuAuc20"]->fill(partPt);
								hKaonPt["Raa_c20_AuAu"]->fill(partPt);
								break;
							}

							case 2212: // proton
							{
								hProtPosPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
								rpp_ProtPos["AuAuc20"]->fill(partPt);
								rppi_ProtonPos["AuAuc20"]->fill(partPt);
								hProtPt["Raa_c20_AuAu"]->fill(partPt);
								break;
							}

							case -2212: // anti-proton
							{
								hProtNegPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
								rpp_ProtNeg["AuAuc20"]->fill(partPt);
								rppi_ProtonNeg["AuAuc20"]->fill(partPt);
								hProtPt["Raa_c20_AuAu"]->fill(partPt);
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
								hPionPosPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
								rpipi_PionPos["AuAuc40"]->fill(partPt);
								rKpi_PionPos["AuAuc40"]->fill(partPt);
								rppi_PionPos["AuAuc40"]->fill(partPt);
								hPionPt["Raa_c40_AuAu"]->fill(partPt);
								break;
							}
							case -211: // pi-
							{
								hPionNegPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
								rpipi_PionNeg["AuAuc40"]->fill(partPt);
								rKpi_PionNeg["AuAuc40"]->fill(partPt);
								rppi_PionNeg["AuAuc40"]->fill(partPt);
								hPionPt["Raa_c40_AuAu"]->fill(partPt);
								break;
							}

							case 321: // K+
							{
								hKaonPosPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
								rKK_KaonPos["AuAuc40"]->fill(partPt);
								rKpi_KaonPos["AuAuc40"]->fill(partPt);
								hKaonPt["Raa_c40_AuAu"]->fill(partPt);
								break;
							}

							case -321: // K-
							{
								hKaonNegPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
								rKK_KaonNeg["AuAuc40"]->fill(partPt);
								rKpi_KaonNeg["AuAuc40"]->fill(partPt);
								hKaonPt["Raa_c40_AuAu"]->fill(partPt);
								break;
							}

							case 2212: // proton
							{
								hProtPosPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
								rpp_ProtPos["AuAuc40"]->fill(partPt);
								rppi_ProtonPos["AuAuc40"]->fill(partPt);
								hProtPt["Raa_c40_AuAu"]->fill(partPt);
								break;
							}

							case -2212: // anti-proton
							{
								hProtNegPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
								rpp_ProtNeg["AuAuc40"]->fill(partPt);
								rppi_ProtonNeg["AuAuc40"]->fill(partPt);
								hProtPt["Raa_c40_AuAu"]->fill(partPt);
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
								hPionPosPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
								rpipi_PionPos["AuAuc60"]->fill(partPt);
								rKpi_PionPos["AuAuc60"]->fill(partPt);
								rppi_PionPos["AuAuc60"]->fill(partPt);
								hPionPt["Raa_c60_AuAu"]->fill(partPt);
								hPionPosPt["AuAuc4060"]->fill(partPt);
								break;
							}
							case -211: // pi-
							{
								hPionNegPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
								rpipi_PionNeg["AuAuc60"]->fill(partPt);
								rKpi_PionNeg["AuAuc60"]->fill(partPt);
								rppi_PionNeg["AuAuc60"]->fill(partPt);
								hPionPt["Raa_c60_AuAu"]->fill(partPt);
								hPionNegPt["AuAuc4060"]->fill(partPt);
								break;
							}

							case 321: // K+
							{
								hKaonPosPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
								rKK_KaonPos["AuAuc60"]->fill(partPt);
								rKpi_KaonPos["AuAuc60"]->fill(partPt);
								hKaonPt["Raa_c60_AuAu"]->fill(partPt);
								hKaonPosPt["AuAuc4060"]->fill(partPt);
								break;
							}

							case -321: // K-
							{
								hKaonNegPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
								rKK_KaonNeg["AuAuc60"]->fill(partPt);
								rKpi_KaonNeg["AuAuc60"]->fill(partPt);
								hKaonPt["Raa_c60_AuAu"]->fill(partPt);
								hKaonNegPt["AuAuc4060"]->fill(partPt);
								break;
							}

							case 2212: // proton
							{
								hProtPosPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
								rpp_ProtPos["AuAuc60"]->fill(partPt);
								rppi_ProtonPos["AuAuc60"]->fill(partPt);
								hProtPt["Raa_c60_AuAu"]->fill(partPt);
								hProtPosPt["AuAuc4060"]->fill(partPt);
								break;
							}

							case -2212: // anti-proton
							{
								hProtNegPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
								rpp_ProtNeg["AuAuc60"]->fill(partPt);
								rppi_ProtonNeg["AuAuc60"]->fill(partPt);
								hProtPt["Raa_c60_AuAu"]->fill(partPt);
								hProtNegPt["AuAuc4060"]->fill(partPt);
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
								hPionPosPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
								rpipi_PionPos["AuAuc92"]->fill(partPt);
								rKpi_PionPos["AuAuc92"]->fill(partPt);
								rppi_PionPos["AuAuc92"]->fill(partPt);
								hPionPt["Raa_c92_AuAu"]->fill(partPt);
								hPionPosPt["AuAuc6092"]->fill(partPt);
								hPionPt["AuAuc6092"]->fill(partPt);
								break;
							}
							case -211: // pi-
							{
								hPionNegPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
								rpipi_PionNeg["AuAuc92"]->fill(partPt);
								rKpi_PionNeg["AuAuc92"]->fill(partPt);
								rppi_PionNeg["AuAuc92"]->fill(partPt);
								hPionPt["Raa_c92_AuAu"]->fill(partPt);
								hPionNegPt["AuAuc6092"]->fill(partPt);
								hPionPt["AuAuc6092"]->fill(partPt);
								break;
							}

							case 321: // K+
							{
								hKaonPosPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
								rKK_KaonPos["AuAuc92"]->fill(partPt);
								rKpi_KaonPos["AuAuc92"]->fill(partPt);
								hKaonPt["Raa_c92_AuAu"]->fill(partPt);
								hKaonPosPt["AuAuc6092"]->fill(partPt);
								hKaonPt["AuAuc6092"]->fill(partPt);
								break;
							}

							case -321: // K-
							{
								hKaonNegPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
								rKK_KaonNeg["AuAuc92"]->fill(partPt);
								rKpi_KaonNeg["AuAuc92"]->fill(partPt);
								hKaonPt["Raa_c92_AuAu"]->fill(partPt);
								hKaonNegPt["AuAuc6092"]->fill(partPt);
								hKaonPt["AuAuc6092"]->fill(partPt);
								break;
							}

							case 2212: // proton
							{
								hProtPosPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
								rpp_ProtPos["AuAuc92"]->fill(partPt);
								rppi_ProtonPos["AuAuc92"]->fill(partPt);
								hProtPt["Raa_c92_AuAu"]->fill(partPt);
								hProtPosPt["AuAuc6092"]->fill(partPt);
								hProtPt["AuAuc6092"]->fill(partPt);
								break;
							}

							case -2212: // anti-proton
							{
								hProtNegPt["ptyieldsAuAuc92"]->fill(partPt, pt_weight);
								rpp_ProtNeg["AuAuc92"]->fill(partPt);
								rppi_ProtonNeg["AuAuc92"]->fill(partPt);
								hProtPt["Raa_c92_AuAu"]->fill(partPt);
								hProtNegPt["AuAuc6092"]->fill(partPt);
								hProtPt["AuAuc6092"]->fill(partPt);
								break;
							}
							}
						}
					}
					return;
				}

				if (collSys == dAu200)
				{
					if ((c < 0.) || (c > 100.)) vetoEvent;


                    //cout << "Centrality: " << c << endl;

					sow["sow_dAUc100"]->fill();
					for (const Particle& p : chargedParticles)
					{
						double partPt = p.pT() / GeV;
						double pt_weight = 1. / (partPt * 2. * M_PI);

                        //cout << "PID: " << p.pid() << endl;

						switch (p.pid()) {
						case 211: // pi+
						{
							hPionPosPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);
							rpipi_PionPos["dAuc100"]->fill(partPt);
							rKpi_PionPos["dAuc100"]->fill(partPt);
							rppi_PionPos["dAuc100"]->fill(partPt);
							hPionPt["Rda_c100_dAu"]->fill(partPt);
							break;
						}
						case -211: // pi-
						{
							hPionNegPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);
							rpipi_PionNeg["dAuc100"]->fill(partPt);
							rKpi_PionNeg["dAuc100"]->fill(partPt);
							rppi_PionNeg["dAuc100"]->fill(partPt);
							hPionPt["Rda_c100_dAu"]->fill(partPt);
							break;
						}

						case 321: // K+
						{
							hKaonPosPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);
							rKK_KaonPos["dAuc100"]->fill(partPt);
							rKpi_KaonPos["dAuc100"]->fill(partPt);
							hKaonPt["Rda_c100_dAu"]->fill(partPt);
							break;
						}

						case -321: // K-
						{
							hKaonNegPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);
							rKK_KaonNeg["dAuc100"]->fill(partPt);
							rKpi_KaonNeg["dAuc100"]->fill(partPt);
							hKaonPt["Rda_c100_dAu"]->fill(partPt);
							break;
						}

						case 2212: // proton
						{
							hProtPosPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);
							rpp_ProtPos["dAuc100"]->fill(partPt);
							rppi_ProtonPos["dAuc100"]->fill(partPt);
							hProtPt["Rda_c100_dAu"]->fill(partPt);
							break;
						}

						case -2212: // anti-proton
						{
							hProtNegPt["ptyieldsdAuc100"]->fill(partPt, pt_weight);
							rpp_ProtNeg["dAuc100"]->fill(partPt);
							rppi_ProtonNeg["dAuc100"]->fill(partPt);
							hProtPt["Rda_c100_dAu"]->fill(partPt);
							break;
						}
						}
					}

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
								hPionPosPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
								rpipi_PionPos["dAuc20"]->fill(partPt);
								rKpi_PionPos["dAuc20"]->fill(partPt);
								rppi_PionPos["dAuc20"]->fill(partPt);
								hPionPt["Rda_c20_dAu"]->fill(partPt);
								hPionPt["dAuc0020"]->fill(partPt);
								break;
							}
							case -211: // pi-
							{
								hPionNegPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
								rpipi_PionNeg["dAuc20"]->fill(partPt);
								rKpi_PionNeg["dAuc20"]->fill(partPt);
								rppi_PionNeg["dAuc20"]->fill(partPt);
								hPionPt["Rda_c20_dAu"]->fill(partPt);
								hPionPt["dAuc0020"]->fill(partPt);
								break;
							}

							case 321: // K+
							{
								hKaonPosPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
								rKK_KaonPos["dAuc20"]->fill(partPt);
								rKpi_KaonPos["dAuc20"]->fill(partPt);
								hKaonPt["Rda_c20_dAu"]->fill(partPt);
								hKaonPt["dAuc0020"]->fill(partPt);
								break;
							}

							case -321: // K-
							{
								hKaonNegPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
								rKK_KaonNeg["dAuc20"]->fill(partPt);
								rKpi_KaonNeg["dAuc20"]->fill(partPt);
								hKaonPt["Rda_c20_dAu"]->fill(partPt);
								hKaonPt["dAuc0020"]->fill(partPt);
								break;
							}

							case 2212: // proton
							{
								hProtPosPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
								rpp_ProtPos["dAuc20"]->fill(partPt);
								rppi_ProtonPos["dAuc20"]->fill(partPt);
								hProtPt["Rda_c20_dAu"]->fill(partPt);
								hProtPt["dAuc0020"]->fill(partPt);
								break;
							}

							case -2212: // anti-proton
							{
								hProtNegPt["ptyieldsdAuc20"]->fill(partPt, pt_weight);
								rpp_ProtNeg["dAuc20"]->fill(partPt);
								rppi_ProtonNeg["dAuc20"]->fill(partPt);
								hProtPt["Rda_c20_dAu"]->fill(partPt);
								hProtPt["dAuc0020"]->fill(partPt);
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
								hPionPosPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
								rpipi_PionPos["dAuc40"]->fill(partPt);
								rKpi_PionPos["dAuc40"]->fill(partPt);
								rppi_PionPos["dAuc40"]->fill(partPt);
								hPionPt["Rda_c40_dAu"]->fill(partPt);
								break;
							}
							case -211: // pi-
							{
								hPionNegPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
								rpipi_PionNeg["dAuc40"]->fill(partPt);
								rKpi_PionNeg["dAuc40"]->fill(partPt);
								rppi_PionNeg["dAuc40"]->fill(partPt);
								hPionPt["Rda_c40_dAu"]->fill(partPt);
								break;
							}

							case 321: // K+
							{
								hKaonPosPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
								rKK_KaonPos["dAuc40"]->fill(partPt);
								rKpi_KaonPos["dAuc40"]->fill(partPt);
								hKaonPt["Rda_c40_dAu"]->fill(partPt);
								break;
							}

							case -321: // K-
							{
								hKaonNegPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
								rKK_KaonNeg["dAuc40"]->fill(partPt);
								rKpi_KaonNeg["dAuc40"]->fill(partPt);
								hKaonPt["Rda_c40_dAu"]->fill(partPt);
								break;
							}

							case 2212: // proton
							{
								hProtPosPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
								rpp_ProtPos["dAuc40"]->fill(partPt);
								rppi_ProtonPos["dAuc40"]->fill(partPt);
								hProtPt["Rda_c40_dAu"]->fill(partPt);
								break;
							}

							case -2212: // anti-proton
							{
								hProtNegPt["ptyieldsdAuc40"]->fill(partPt, pt_weight);
								rpp_ProtNeg["dAuc40"]->fill(partPt);
								rppi_ProtonNeg["dAuc40"]->fill(partPt);
								hProtPt["Rda_c40_dAu"]->fill(partPt);
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
								hPionPosPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
								rpipi_PionPos["dAuc60"]->fill(partPt);
								rKpi_PionPos["dAuc60"]->fill(partPt);
								rppi_PionPos["dAuc60"]->fill(partPt);
								hPionPt["Rda_c60_dAu"]->fill(partPt);
								break;
							}
							case -211: // pi-
							{
								hPionNegPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
								rpipi_PionNeg["dAuc60"]->fill(partPt);
								rKpi_PionNeg["dAuc60"]->fill(partPt);
								rppi_PionNeg["dAuc60"]->fill(partPt);
								hPionPt["Rda_c60_dAu"]->fill(partPt);
								break;
							}

							case 321: // K+
							{
								hKaonPosPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
								rKK_KaonPos["dAuc60"]->fill(partPt);
								rKpi_KaonPos["dAuc60"]->fill(partPt);
								hKaonPt["Rda_c60_dAu"]->fill(partPt);
								break;
							}

							case -321: // K-
							{
								hKaonNegPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
								rKK_KaonNeg["dAuc60"]->fill(partPt);
								rKpi_KaonNeg["dAuc60"]->fill(partPt);
								hKaonPt["Rda_c60_dAu"]->fill(partPt);
								break;
							}

							case 2212: // proton
							{
								hProtPosPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
								rpp_ProtPos["dAuc60"]->fill(partPt);
								rppi_ProtonPos["dAuc60"]->fill(partPt);
								hProtPt["Rda_c60_dAu"]->fill(partPt);
								break;
							}

							case -2212: // anti-proton
							{
								hProtNegPt["ptyieldsdAuc60"]->fill(partPt, pt_weight);
								rpp_ProtNeg["dAuc60"]->fill(partPt);
								rppi_ProtonNeg["dAuc60"]->fill(partPt);
								hProtPt["Rda_c60_dAu"]->fill(partPt);
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
								hPionPosPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
								rpipi_PionPos["dAuc88"]->fill(partPt);
								rKpi_PionPos["dAuc88"]->fill(partPt);
								rppi_PionPos["dAuc88"]->fill(partPt);
								hPionPt["Rda_c88_dAu"]->fill(partPt);
								break;
							}
							case -211: // pi-
							{
								hPionNegPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
								rpipi_PionNeg["dAuc88"]->fill(partPt);
								rKpi_PionNeg["dAuc88"]->fill(partPt);
								rppi_PionNeg["dAuc88"]->fill(partPt);
								hPionPt["Rda_c88_dAu"]->fill(partPt);
								break;
							}

							case 321: // K+
							{
								hKaonPosPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
								rKK_KaonPos["dAuc88"]->fill(partPt);
								rKpi_KaonPos["dAuc88"]->fill(partPt);
								hKaonPt["Rda_c88_dAu"]->fill(partPt);
								break;
							}

							case -321: // K-
							{
								hKaonNegPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
								rKK_KaonNeg["dAuc88"]->fill(partPt);
								rKpi_KaonNeg["dAuc88"]->fill(partPt);
								hKaonPt["Rda_c88_dAu"]->fill(partPt);
								break;
							}

							case 2212: // proton
							{
								hProtPosPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
								rpp_ProtPos["dAuc88"]->fill(partPt);
								rppi_ProtonPos["dAuc88"]->fill(partPt);
								hProtPt["Rda_c88_dAu"]->fill(partPt);
								break;
							}

							case -2212: // anti-proton
							{
								hProtNegPt["ptyieldsdAuc88"]->fill(partPt, pt_weight);
								rpp_ProtNeg["dAuc88"]->fill(partPt);
								rppi_ProtonNeg["dAuc88"]->fill(partPt);
								hProtPt["Rda_c88_dAu"]->fill(partPt);
								break;
							}
							}
						}
					}
					return;
				}
			}

			void finalize() {

				for (int i = 0, N = AUAUCentralityBins.size(); i < N; ++i)
				{
					//yields (fig 4)_________________
					hKaonNegPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
					hPionNegPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
					hProtNegPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
					hKaonPosPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
					hPionPosPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
					hProtPosPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
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

					divide(rKpi_KaonPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], rKpi_PionPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], RatioK_pipos["AuAuc" + std::to_string(AUAUCentralityBins[i])]);
					divide(rKpi_KaonPos["dAuc" + std::to_string(dAUCentralityBins[i])], rKpi_PionPos["dAuc" + std::to_string(dAUCentralityBins[i])], RatioK_pipos["dAuc" + std::to_string(dAUCentralityBins[i])]);

					divide(rKpi_KaonNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], rKpi_PionNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], RatioK_pineg["AuAuc" + std::to_string(AUAUCentralityBins[i])]);
					divide(rKpi_KaonNeg["dAuc" + std::to_string(dAUCentralityBins[i])], rKpi_PionNeg["dAuc" + std::to_string(dAUCentralityBins[i])], RatioK_pineg["dAuc" + std::to_string(dAUCentralityBins[i])]);

					divide(rppi_ProtonPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], rppi_PionPos["AuAuc" + std::to_string(AUAUCentralityBins[i])], Ratiop_pipos["AuAuc" + std::to_string(AUAUCentralityBins[i])]);
					divide(rppi_ProtonPos["dAuc" + std::to_string(dAUCentralityBins[i])], rppi_PionPos["dAuc" + std::to_string(dAUCentralityBins[i])], Ratiop_pipos["dAuc" + std::to_string(dAUCentralityBins[i])]);

					divide(rppi_ProtonNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], rppi_PionNeg["AuAuc" + std::to_string(AUAUCentralityBins[i])], Ratiop_pineg["AuAuc" + std::to_string(AUAUCentralityBins[i])]);
					divide(rppi_ProtonNeg["dAuc" + std::to_string(dAUCentralityBins[i])], rppi_PionNeg["dAuc" + std::to_string(dAUCentralityBins[i])], Ratiop_pineg["dAuc" + std::to_string(dAUCentralityBins[i])]);

					//RAA (fig 11)_________________
					hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
					hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
					divide(hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], hKaonPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], hRaa["K_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]);

					hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
					hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
					divide(hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], hPionPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], hRaa["pi_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]);

					hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
					hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
					divide(hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"], hProtPt["Raa_c" + std::to_string(AUAUCentralityBins[i]) + "_pp"], hRaa["p_c" + std::to_string(AUAUCentralityBins[i]) + "_AuAu"]);

					//RdA (fig 12)________________

					hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]->scaleW(1. / sow["sow_dAUc" + std::to_string(dAUCentralityBins[i])]->sumW());
					hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
					divide(hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], hKaonPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"], hRda["K_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]);

					hPionPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]->scaleW(1. / sow["sow_dAUc" + std::to_string(dAUCentralityBins[i])]->sumW());
					hPionPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
					divide(hPionPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"], hPionPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_pp"], hRda["pi_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]);

					hProtPt["Rda_c" + std::to_string(dAUCentralityBins[i]) + "_dAu"]->scaleW(1. / sow["sow_dAUc" + std::to_string(dAUCentralityBins[i])]->sumW());
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
				hRaa["K_c10_AuAu"]->scaleY(1. / 960.2);
				hRaa["K_c20_AuAu"]->scaleY(1. / 609.5);
				hRaa["K_c40_AuAu"]->scaleY(1. / 300.8);
				hRaa["K_c60_AuAu"]->scaleY(1. / 94.2);
				hRaa["K_c92_AuAu"]->scaleY(1. / 14.8);

				hRaa["pi_c10_AuAu"]->scaleY(1. / 960.2);
				hRaa["pi_c20_AuAu"]->scaleY(1. / 609.5);
				hRaa["pi_c40_AuAu"]->scaleY(1. / 300.8);
				hRaa["pi_c60_AuAu"]->scaleY(1. / 94.2);
				hRaa["pi_c92_AuAu"]->scaleY(1. / 14.8);

				hRaa["p_c10_AuAu"]->scaleY(1. / 960.2);
				hRaa["p_c20_AuAu"]->scaleY(1. / 609.5);
				hRaa["p_c40_AuAu"]->scaleY(1. / 300.8);
				hRaa["p_c60_AuAu"]->scaleY(1. / 94.2);
				hRaa["p_c92_AuAu"]->scaleY(1. / 14.8);

				hRda["K_c20_dAu"]->scaleY(1. / 15.1);
				hRda["K_c40_dAu"]->scaleY(1. / 10.2);
				hRda["K_c60_dAu"]->scaleY(1. / 6.6);
				hRda["K_c88_dAu"]->scaleY(1. / 3.1);
				hRda["K_c100_dAu"]->scaleY(1. / 7.6);

				hRda["pi_c20_dAu"]->scaleY(1. / 15.1);
				hRda["pi_c40_dAu"]->scaleY(1. / 10.2);
				hRda["pi_c60_dAu"]->scaleY(1. / 6.6);
				hRda["pi_c88_dAu"]->scaleY(1. / 3.1);
				hRda["pi_c100_dAu"]->scaleY(1. / 7.6);

				hRda["p_c20_dAu"]->scaleY(1. / 15.1);
				hRda["p_c40_dAu"]->scaleY(1. / 10.2);
				hRda["p_c60_dAu"]->scaleY(1. / 6.6);
				hRda["p_c88_dAu"]->scaleY(1. / 3.1);
				hRda["p_c100_dAu"]->scaleY(1. / 7.6);


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



			map<string, Histo1DPtr> rKK_KaonNeg;
			map<string, Histo1DPtr> rKK_KaonPos;
			map<string, Scatter2DPtr> RatioKaon;

			map<string, Histo1DPtr> rpipi_PionNeg;
			map<string, Histo1DPtr> rpipi_PionPos;
			map<string, Scatter2DPtr> RatioPion;

			map<string, Histo1DPtr> rpp_ProtNeg;
			map<string, Histo1DPtr> rpp_ProtPos;
			map<string, Scatter2DPtr> RatioProt;

			map<string, Histo1DPtr> rKpi_KaonPos;
			map<string, Histo1DPtr> rKpi_PionPos;
			map<string, Scatter2DPtr> RatioK_pipos;

			map<string, Histo1DPtr> rKpi_KaonNeg;
			map<string, Histo1DPtr> rKpi_PionNeg;
			map<string, Scatter2DPtr> RatioK_pineg;

			map<string, Histo1DPtr> rKpi_ProtPos;
			map<string, Scatter2DPtr> Ratiop_pipos;

			map<string, Histo1DPtr> rKpi_ProtNeg;
			map<string, Scatter2DPtr> Ratiop_pineg;

			map<string, Histo1DPtr> rppi_ProtonPos;
			map<string, Histo1DPtr> rppi_PionPos;
			map<string, Scatter2DPtr> Ratiop_piPos;

			map<string, Histo1DPtr> rppi_ProtonNeg;
			map<string, Histo1DPtr> rppi_PionNeg;
			map<string, Scatter2DPtr> Ratiop_piNeg;


			map<string, Scatter2DPtr> RatioK;
			map<string, Scatter2DPtr> Ratiopi;
			map<string, Scatter2DPtr> Ratiop;

			map<string, Scatter2DPtr> hRcp;
			map<string, Scatter2DPtr> hRaa;
			map<string, Scatter2DPtr> hRda;

			map<string, CounterPtr> sow;
			string beamOpt;
			enum CollisionSystem { pp, AuAu200, dAu200 };
			CollisionSystem collSys;
			vector<int> AUAUCentralityBins{ 10, 20, 40, 60, 92 };
			vector<int> dAUCentralityBins{ 20, 100, 40, 60, 88 };


	};


	DECLARE_RIVET_PLUGIN(PHENIX_2013_I1227971);

}
