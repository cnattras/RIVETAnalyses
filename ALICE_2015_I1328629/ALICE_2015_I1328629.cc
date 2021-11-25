// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include <stdio.h>

namespace Rivet {

	/// @brief Add a short analysis description here
	class ALICE_2015_I1328629 : public Analysis {
		public:

			/// Constructor
			RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2015_I1328629);

			/// @name Analysis methods
			///@{

			void FillRadialDistribution(Profile1DPtr prof1D, Jet jet)
			{
				std::vector<double> sum(prof1D->numBins(), 0.);
				for(auto p : jet.particles())
				{
					double r = deltaR(p.eta(), p.phi(), jet.eta(), jet.phi());
					int bin = prof1D->binIndexAt(r);
					if(bin >= 0) sum[bin] += p.pT()/GeV;
				}
				for(unsigned int i = 0; i < sum.size(); i++) {
					prof1D->fill(prof1D->bin(i).xMid(),sum[i]/prof1D->bin(i).xWidth());
				}
			}

			void FillRadius80PercPt(Profile1DPtr prof1D, Jet jet)
			{
				std::vector<pair<double, Particle>> jetConsRPt;
				double Jet80Pt = 0.;
				double r80 = 0.;

				for(Particle p : jet.particles())
				{
					double r = deltaR(p.eta(), p.phi(), jet.eta(), jet.phi());
					jetConsRPt.push_back(make_pair(r, p));
				}
				// Sort the vector jetConsRPt by the first argument r (c++14)
				std::sort(jetConsRPt.begin(), jetConsRPt.end(),[](auto &left, auto &right) {
					return left.first < right.first;
				});

				for(auto pair : jetConsRPt)
				{
					Jet80Pt += pair.second.pT()/GeV;
					if (Jet80Pt >= 0.8 * jet.pT()/GeV){
						r80 = deltaR(pair.second.eta(), pair.second.phi(), jet.eta(), jet.phi());
						break;
					}
				}
				prof1D->fill(jet.pT()/GeV, r80);
			}

			double GetRho(Particles eventParticles, double jetR, Jet leadingJet)
			{
				double UEphiPlus = mapAngle0To2Pi(leadingJet.phi() + M_PI/2.);
				double UEphiMinus = mapAngle0To2Pi(leadingJet.phi() - M_PI/2.);
				double UEeta = leadingJet.eta();

				double rho = 0.;

				for(const Particle& p : eventParticles)
				{
					double DeltaRPlus = deltaR(p.eta(), p.phi(), UEeta, UEphiPlus);
					double DeltaRMinus = deltaR(p.eta(), p.phi(), UEeta, UEphiMinus);
					if(DeltaRPlus > jetR && DeltaRMinus > jetR) continue;

					rho += p.pT()/GeV;
				}

				//Dividing by 2 because two cones were used and by cone area piR^2
				//pT-density
				rho /= 2.*M_PI*jetR*jetR;

				return rho;
			}

                        void FillUEDist(Particles eventParticles, double jetR, Jet leadingJet, Histo1DPtr HistoUEPt, Histo1DPtr HistoUEZ, Histo1DPtr HistoUEXi)
			{
				double UEphiPlus = mapAngle0To2Pi(leadingJet.phi() + M_PI/2.);
				double UEphiMinus = mapAngle0To2Pi(leadingJet.phi() - M_PI/2.);
				double UEeta = leadingJet.eta();

				for(const Particle& p : eventParticles)
				{
					double DeltaRPlus = deltaR(p.eta(), p.phi(), UEeta, UEphiPlus);
					double DeltaRMinus = deltaR(p.eta(), p.phi(), UEeta, UEphiMinus);
                                        if(DeltaRPlus > jetR && DeltaRMinus > jetR) continue;

                                        HistoUEPt->fill(p.pT()/GeV);
                                        HistoUEZ->fill(p.pT()/leadingJet.pT());
                                        HistoUEXi->fill(log(leadingJet.pT()/p.pT()));
				}

			}

			double GetJetPtCorr(Jet jet, double rho)
			{
				return jet.pT()/GeV - (rho*jet.pseudojet().area());
			}

                        void Subtract(Histo1DPtr Histo, Histo1DPtr hSub)
                        {
                                const string path = Histo->path();
                                *Histo = YODA::subtract(*Histo, *hSub);
                                Histo->setPath(path);
                        }

			/// Book histograms and initialise projections before the run
			void init() {

				// The basic final-state projection: all final-state particles within the given eta acceptance
				const FinalState fs(Cuts::abseta < 0.75);

				/*
					The final-state particles declared above are clustered using FastJet with the anti-kT
					algorithm and a jet-radius parameter 0.4 muons & neutrinos are excluded from the clustering
				*/

				// For the area of the jets
				fastjet::AreaType fjAreaType = fastjet::active_area_explicit_ghosts;
				fastjet::GhostedAreaSpec fjGhostAreaSpec = fastjet::GhostedAreaSpec(1., 1, 0.005, 1., 0.1, 1e-100);

				fjAreaDef02 = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
				FastJets jetsAKTR02FJ(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.2, fjAreaDef02, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKTR02FJ, "jetsAKTR02FJ");

				fjAreaDef03 = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
				FastJets jetsAKTR03FJ(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.3, fjAreaDef03, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKTR03FJ, "jetsAKTR03FJ");

				fjAreaDef04 = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
				FastJets jetsAKTR04FJ(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.4, fjAreaDef04, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKTR04FJ, "jetsAKTR04FJ");

				fjAreaDef04KT = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
				FastJets jetsKTR04FJ(fs, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.4, fjAreaDef04KT, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsKTR04FJ, "jetsKTR04FJ");

				FastJets jetsCONER04FJ(fs, FastJets::SISCONE, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsCONER04FJ, "jetsCONER04FJ");

				fjAreaDef06 = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
				FastJets jetsAKTR06FJ(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.6, fjAreaDef06, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKTR06FJ, "jetsAKTR06FJ");

				// Initialize ALICE primary particles
				const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::abscharge > 0);
				declare(aprim, "aprim");

				// Create counters
				book(_c["sow"], "sow");
				book(_c["sow2030"], "sow2030");
				book(_c["sow3040"], "sow3040");
				book(_c["sow4060"], "sow4060");
				book(_c["sow6080"], "sow6080");
                                book(_c["sow2030UESub"], "sow2030UESub");
				book(_c["sow3040UESub"], "sow3040UESub");
				book(_c["sow4060UESub"], "sow4060UESub");
				book(_c["sow6080UESub"], "sow6080UESub");

				// Creates histograms
				book(_h["CrossSectionAntikT_R04"], 1, 1, 1);  // Figure 2 Anti-kT
				book(_h["CrossSectionkT_R04"], 2, 1, 1);      // Figure 2 kT
				book(_h["CrossSectionSisCone_R04"], 3, 1, 1); // Figure 2 SisCone

				// Figure 3
				book(_h["antiKTR02"], 4, 1, 1);
				book(_h["antiKTR03"], 5, 1, 1);
				book(_h["antiKTR04"], 6, 1, 1);
				book(_h["antiKTR06"], 7, 1, 1);

				// Figure 4
				book(_h["antiKT_R04_ALICE_ATLAS"], 8, 1, 1);
				book(_h["antiKT_R06_ALICE_ATLAS"], 9, 1, 1);

				// Figure 5
				book(_h["ALICEvsMC_R02_Eta07"], 10, 1, 1);
				book(_h["ALICEvsMC_R04_Eta05"], 11, 1, 1);
				book(_h["ALICEvsMC_R06_Eta03"], 12, 1, 1);

				// Figure 6 - Ratios
				string refname0204 = mkAxisCode(13, 1, 1);
				const Scatter2D& refdata0204 = refData(refname0204);
				book(_h["Ratio02_Numerator"], refname0204 + "Ratio02_Numerator", refdata0204);
				book(_h["Ratio04_Denominator"], refname0204 + "Ratio04_Denominator", refdata0204);
				book(_s["Ratio0204"], refname0204);
				string refname0206 = mkAxisCode(14, 1, 1);
				const Scatter2D& refdata0206 = refData(refname0206);
				book(_h["Ratio06_Denominator"], refname0206 + "Ratio06_Denominator", refdata0206);
				book(_s["Ratio0206"], refname0206);

				// Figure 7
				book(_p["mean_ALICEvsMC_R02_Eta_07"], 15, 1, 1);
				book(_p["mean_ALICEvsMC_R04_Eta_05"], 16, 1, 1);
				book(_p["mean_ALICEvsMC_R06_Eta_03"], 17, 1, 1);

				// Figure 8
				book(_p["pTR02_Eta07_2030"], 18, 1, 1);
				book(_p["pTR02_Eta07_3040"], 19, 1, 1);
				book(_p["pTR02_Eta07_4060"], 20, 1, 1);
				book(_p["pTR02_Eta07_6080"], 21, 1, 1);

				// Figure 9
				book(_p["pTR04_Eta05_2030"], 22, 1, 1);
				book(_p["pTR04_Eta05_3040"], 23, 1, 1);
				book(_p["pTR04_Eta05_4060"], 24, 1, 1);
				book(_p["pTR04_Eta05_6080"], 25, 1, 1);

				// Figure 10
				book(_p["pTR06_Eta03_2030"], 26, 1, 1);
				book(_p["pTR06_Eta03_3040"], 27, 1, 1);
				book(_p["pTR06_Eta03_4060"], 28, 1, 1);
				book(_p["pTR06_Eta03_6080"], 29, 1, 1);

				// Figure 11
				book(_p["avgpT_R02_Eta07"], 30, 1, 1);
				book(_p["avgpT_R04_Eta05"], 31, 1, 1);
				book(_p["avgpT_R06_Eta03"], 32, 1, 1);

				// Figure 12
                                string refnameFig12_33 = mkAxisCode(33, 1, 1);
				const Scatter2D& refdataFig12_33 = refData(refnameFig12_33);
                                book(_h["pTSpectraR04_Eta05_2030_ALICEvsMC"], refnameFig12_33, refdataFig12_33);
                                book(_h["pTSpectraR04_Eta05_2030_ALICEvsMC_UE"], refnameFig12_33 + "_UE", refdataFig12_33);

                                string refnameFig12_34 = mkAxisCode(34, 1, 1);
				const Scatter2D& refdataFig12_34 = refData(refnameFig12_34);
                                book(_h["pTSpectraR04_Eta05_3040_ALICEvsMC"], refnameFig12_34, refdataFig12_34);
                                book(_h["pTSpectraR04_Eta05_3040_ALICEvsMC_UE"], refnameFig12_34 + "_UE", refdataFig12_34);

                                string refnameFig12_35 = mkAxisCode(35, 1, 1);
				const Scatter2D& refdataFig12_35 = refData(refnameFig12_35);
                                book(_h["pTSpectraR04_Eta05_4060_ALICEvsMC"], refnameFig12_35, refdataFig12_35);
                                book(_h["pTSpectraR04_Eta05_4060_ALICEvsMC_UE"], refnameFig12_35 + "_UE", refdataFig12_35);

                                string refnameFig12_36 = mkAxisCode(36, 1, 1);
				const Scatter2D& refdataFig12_36 = refData(refnameFig12_36);
                                book(_h["pTSpectraR04_Eta05_6080_ALICEvsMC"], refnameFig12_36, refdataFig12_36);
                                book(_h["pTSpectraR04_Eta05_6080_ALICEvsMC_UE"], refnameFig12_36 + "_UE", refdataFig12_36);

				// Figure 13
                                string refnameFig13_37 = mkAxisCode(37, 1, 1);
				const Scatter2D& refdataFig13_37 = refData(refnameFig13_37);
				book(_h["pTSpectraR04_Eta05_2030_0to1"], refnameFig13_37, refdataFig13_37);
                                book(_h["pTSpectraR04_Eta05_2030_0to1_UE"], refnameFig13_37 + "_UE", refdataFig13_37);

                                string refnameFig13_38 = mkAxisCode(38, 1, 1);
				const Scatter2D& refdataFig13_38 = refData(refnameFig13_38);
				book(_h["pTSpectraR04_Eta05_3040_0to1"], refnameFig13_38, refdataFig13_38);
                                book(_h["pTSpectraR04_Eta05_3040_0to1_UE"], refnameFig13_38 + "_UE", refdataFig13_38);

                                string refnameFig13_39 = mkAxisCode(39, 1, 1);
				const Scatter2D& refdataFig13_39 = refData(refnameFig13_39);
				book(_h["pTSpectraR04_Eta05_4060_0to1"], refnameFig13_39, refdataFig13_39);
                                book(_h["pTSpectraR04_Eta05_4060_0to1_UE"], refnameFig13_39 + "_UE", refdataFig13_39);

                                string refnameFig13_40 = mkAxisCode(40, 1, 1);
				const Scatter2D& refdataFig13_40 = refData(refnameFig13_40);
				book(_h["pTSpectraR04_Eta05_6080_0to1"], refnameFig13_40, refdataFig13_40);
                                book(_h["pTSpectraR04_Eta05_6080_0to1_UE"], refnameFig13_40 + "_UE", refdataFig13_40);

				// Figure 14
                                string refnameFig14_41 = mkAxisCode(41, 1, 1);
				const Scatter2D& refdataFig14_41 = refData(refnameFig14_41);
				book(_h["pTSpectraR04_Eta05_2030_0to6"], refnameFig14_41, refdataFig14_41);
                                book(_h["pTSpectraR04_Eta05_2030_0to6_UE"], refnameFig14_41 + "_UE", refdataFig14_41);

                                string refnameFig14_42 = mkAxisCode(42, 1, 1);
				const Scatter2D& refdataFig14_42 = refData(refnameFig14_42);
				book(_h["pTSpectraR04_Eta05_3040_0to6"], refnameFig14_42, refdataFig14_42);
                                book(_h["pTSpectraR04_Eta05_3040_0to6_UE"], refnameFig14_42 + "_UE", refdataFig14_42);

                                string refnameFig14_43 = mkAxisCode(43, 1, 1);
				const Scatter2D& refdataFig14_43 = refData(refnameFig14_43);
				book(_h["pTSpectraR04_Eta05_4060_0to6"], refnameFig14_43, refdataFig14_43);
                                book(_h["pTSpectraR04_Eta05_4060_0to6_UE"], refnameFig14_43 + "_UE", refdataFig14_43);

                                string refnameFig14_44 = mkAxisCode(44, 1, 1);
				const Scatter2D& refdataFig14_44 = refData(refnameFig14_44);
				book(_h["pTSpectraR04_Eta05_6080_0to6"], refnameFig14_44, refdataFig14_44);
                                book(_h["pTSpectraR04_Eta05_6080_0to6_UE"], refnameFig14_44 + "_UE", refdataFig14_44);

				// Figure A.1 - Like Figure 5
				book(_h["ALICEvsMC_R04_Eta05_WithoutUESub"], 45, 1, 1);
				book(_h["ALICEvsMC_R06_Eta03_WithoutUESub"], 46, 1, 1);

				// Figure A.2 - Like Figure 7
				book(_p["mean_ALICEvsMC_R02_Eta_07_WithoutUESub"], 47, 1, 1);
				book(_p["mean_ALICEvsMC_R04_Eta_05_WithoutUESub"], 48, 1, 1);
				book(_p["mean_ALICEvsMC_R06_Eta_03_WithoutUESub"], 49, 1, 1);

				// Figure A.3 - Like Figure 8
				book(_p["pTR02_Eta07_2030_WithoutUESub"], 50, 1, 1);
				book(_p["pTR02_Eta07_3040_WithoutUESub"], 51, 1, 1);
				book(_p["pTR02_Eta07_4060_WithoutUESub"], 52, 1, 1);
				book(_p["pTR02_Eta07_6080_WithoutUESub"], 53, 1, 1);

				// Figure A.4 - Like Figure 9
				book(_p["pTR04_Eta05_2030_WithoutUESub"], 54, 1, 1);
				book(_p["pTR04_Eta05_3040_WithoutUESub"], 55, 1, 1);
				book(_p["pTR04_Eta05_4060_WithoutUESub"], 56, 1, 1);
				book(_p["pTR04_Eta05_6080_WithoutUESub"], 57, 1, 1);

				// Figure A.5 - Like Figure 10
				book(_p["pTR06_Eta03_2030_WithoutUESub"], 58, 1, 1);
				book(_p["pTR06_Eta03_3040_WithoutUESub"], 59, 1, 1);
				book(_p["pTR06_Eta03_4060_WithoutUESub"], 60, 1, 1);
				book(_p["pTR06_Eta03_6080_WithoutUESub"], 61, 1, 1);

				// Figure A.6 - Like Figure 12
				book(_h["pTSpectraR04_Eta05_2030_ALICEvsMC_WithoutUESub"], 62, 1, 1);
				book(_h["pTSpectraR04_Eta05_3040_ALICEvsMC_WithoutUESub"], 63, 1, 1);
				book(_h["pTSpectraR04_Eta05_4060_ALICEvsMC_WithoutUESub"], 64, 1, 1);
				book(_h["pTSpectraR04_Eta05_6080_ALICEvsMC_WithoutUESub"], 65, 1, 1);

				// Figure A.7 - Like Figure 13
				book(_h["pTSpectraR04_Eta05_2030_0to1_WithoutUESub"], 66, 1, 1);
				book(_h["pTSpectraR04_Eta05_3040_0to1_WithoutUESub"], 67, 1, 1);
				book(_h["pTSpectraR04_Eta05_4060_0to1_WithoutUESub"], 68, 1, 1);
				book(_h["pTSpectraR04_Eta05_6080_0to1_WithoutUESub"], 69, 1, 1);

				// Figure A.8 - Like Figure 14
				book(_h["pTSpectraR04_Eta05_2030_0to6_WithoutUESub"], 70, 1, 1);
				book(_h["pTSpectraR04_Eta05_3040_0to6_WithoutUESub"], 71, 1, 1);
				book(_h["pTSpectraR04_Eta05_4060_0to6_WithoutUESub"], 72, 1, 1);
				book(_h["pTSpectraR04_Eta05_6080_0to6_WithoutUESub"], 73, 1, 1);


			}


			/// Perform the per-event analysis
			void analyze(const Event& event) {
				//Increment the counter to keep track of the cross section
				_c["sow"]->fill();

				// Retrieve clustered jets, sorted by pT, with a minimum pT cut
				const ALICE::PrimaryParticles aprim = apply<ALICE::PrimaryParticles>(event, "aprim");
				const Particles ALICEparticles = aprim.particles();

				// Anti-KT alg. - Resolution = 0.2, Eta = 0.7 (From 0.9 - 0.2)
				FastJets jetsAKTR02FJ = apply<FastJets>(event, "jetsAKTR02FJ");
				jetsAKTR02FJ.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
				Jets jetsAKTR02 = jetsAKTR02FJ.jetsByPt(Cuts::pT >= 20*GeV && Cuts::abseta < 0.7);
				double rho02 = 0.;

				if(jetsAKTR02.size() != 0) {
					rho02 = GetRho(ALICEparticles, 0.2, jetsAKTR02[0]);

					_p["mean_ALICEvsMC_R02_Eta_07"]->fill(GetJetPtCorr(jetsAKTR02[0], rho02), jetsAKTR02[0].particles().size()); // Figure 7
					_p["mean_ALICEvsMC_R02_Eta_07_WithoutUESub"]->fill(jetsAKTR02[0].pT()/GeV, jetsAKTR02[0].particles().size()); // Figure A.2
					FillRadius80PercPt(_p["avgpT_R02_Eta07"], jetsAKTR02[0]); // Figure 11

					if ((GetJetPtCorr(jetsAKTR02[0], rho02) >= 20*GeV) && (GetJetPtCorr(jetsAKTR02[0], rho02) < 30*GeV)) {
						FillRadialDistribution(_p["pTR02_Eta07_2030"], jetsAKTR02[0]); // Figure 8
					}
					if ((GetJetPtCorr(jetsAKTR02[0], rho02) >= 30*GeV) && (GetJetPtCorr(jetsAKTR02[0], rho02) < 40*GeV)) {
						FillRadialDistribution(_p["pTR02_Eta07_3040"], jetsAKTR02[0]); // Figure 8
					}
					if ((GetJetPtCorr(jetsAKTR02[0], rho02) >= 40*GeV) && (GetJetPtCorr(jetsAKTR02[0], rho02) < 60*GeV)) {
						FillRadialDistribution(_p["pTR02_Eta07_4060"], jetsAKTR02[0]); // Figure 8
					}
					if ((GetJetPtCorr(jetsAKTR02[0], rho02) >= 60*GeV) && (GetJetPtCorr(jetsAKTR02[0], rho02) < 80*GeV)) {
						FillRadialDistribution(_p["pTR02_Eta07_6080"], jetsAKTR02[0]); // Figure 8
					}

					if(jetsAKTR02[0].pT() >= 20*GeV && jetsAKTR02[0].pT() < 30*GeV) {
						FillRadialDistribution(_p["pTR02_Eta07_2030_WithoutUESub"], jetsAKTR02[0]); // Figure A.3
					}
					if(jetsAKTR02[0].pT() >= 30*GeV && jetsAKTR02[0].pT() < 40*GeV) {
						FillRadialDistribution(_p["pTR02_Eta07_3040_WithoutUESub"], jetsAKTR02[0]); // Figure A.3
					}
					if(jetsAKTR02[0].pT() >= 40*GeV && jetsAKTR02[0].pT() < 60*GeV) {
						FillRadialDistribution(_p["pTR02_Eta07_4060_WithoutUESub"], jetsAKTR02[0]); // Figure A.3
					}
					if(jetsAKTR02[0].pT() >= 60*GeV && jetsAKTR02[0].pT() < 80*GeV) {
						FillRadialDistribution(_p["pTR02_Eta07_6080_WithoutUESub"], jetsAKTR02[0]); // Figure A.3
					}
				}

				for(auto jet: jetsAKTR02) {
					_h["antiKTR02"]->fill(GetJetPtCorr(jet, rho02)); // Figure 3
					_h["ALICEvsMC_R02_Eta07"]->fill(GetJetPtCorr(jet, rho02)); // Figure 5
					_h["Ratio02_Numerator"]->fill(GetJetPtCorr(jet, rho02)); // Figure 6
				}

				// Anti-KT alg. - Resolution = 0.3, Eta = 0.6 (From 0.9 - 0.3)
				FastJets jetsAKTR03FJ = apply<FastJets>(event, "jetsAKTR03FJ");
				jetsAKTR03FJ.calc(ALICEparticles);
				Jets jetsAKTR03 = jetsAKTR03FJ.jetsByPt(Cuts::pT >= 20*GeV && Cuts::abseta < 0.6);
				double rho03 = 0;

				if (jetsAKTR03.size() != 0) {
					rho03 = GetRho(ALICEparticles, 0.3, jetsAKTR03[0]);
				}

				for(auto jet: jetsAKTR03) {
					_h["antiKTR03"]->fill(GetJetPtCorr(jet, rho03)); // Figure 3
				}

				// Anti-KT alg. - Resolution = 0.4, Eta = 0.5 (From 0.9 - 0.4)
				FastJets jetsAKTR04FJ = apply<FastJets>(event, "jetsAKTR04FJ");
				jetsAKTR04FJ.calc(ALICEparticles);
				Jets jetsAKTR04 = jetsAKTR04FJ.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.5);
				double rho04 = 0;
				if(jetsAKTR04.size() != 0) {
					rho04 = GetRho(ALICEparticles, 0.4, jetsAKTR04[0]);

					_p["mean_ALICEvsMC_R04_Eta_05"]->fill(GetJetPtCorr(jetsAKTR04[0], rho04), jetsAKTR04[0].particles().size()); // Figure 7
					_p["mean_ALICEvsMC_R04_Eta_05_WithoutUESub"]->fill(jetsAKTR04[0].pT()/GeV, jetsAKTR04[0].particles().size()); // Figure A.2
					FillRadius80PercPt(_p["avgpT_R04_Eta05"], jetsAKTR04[0]); // Figure 11

                                        if (inRange(GetJetPtCorr(jetsAKTR04[0], rho04), 20., 30.)) {
                                                _c["sow2030UESub"]->fill();
						FillRadialDistribution(_p["pTR04_Eta05_2030"], jetsAKTR04[0]); // Figure 9
                                                FillUEDist(ALICEparticles, 0.4, jetsAKTR04[0], _h["pTSpectraR04_Eta05_2030_ALICEvsMC_UE"], _h["pTSpectraR04_Eta05_2030_0to1_UE"], _h["pTSpectraR04_Eta05_2030_0to6_UE"]);
						for (auto p : jetsAKTR04[0].particles()) {
							_h["pTSpectraR04_Eta05_2030_ALICEvsMC"]->fill(p.pT()/GeV); // Figure 12
                                                        _h["pTSpectraR04_Eta05_2030_0to1"]->fill(p.pT()/jetsAKTR04[0].pT()); // Figure 13
                                                        _h["pTSpectraR04_Eta05_2030_0to6"]->fill(log(jetsAKTR04[0].pT()/p.pT())); // Figure 14
						}
					}
                                        else if(inRange(GetJetPtCorr(jetsAKTR04[0], rho04), 30., 40.)) {
                                                _c["sow3040UESub"]->fill();
						FillRadialDistribution(_p["pTR04_Eta05_3040"], jetsAKTR04[0]); // Figure 9
                                                FillUEDist(ALICEparticles, 0.4, jetsAKTR04[0], _h["pTSpectraR04_Eta05_3040_ALICEvsMC_UE"], _h["pTSpectraR04_Eta05_3040_0to1_UE"], _h["pTSpectraR04_Eta05_3040_0to6_UE"]);
						for (auto p : jetsAKTR04[0].particles()) {
							_h["pTSpectraR04_Eta05_3040_ALICEvsMC"]->fill(p.pT()/GeV); // Figure 12
                                                        _h["pTSpectraR04_Eta05_3040_0to1"]->fill(p.pT()/jetsAKTR04[0].pT()); // Figure 13
                                                        _h["pTSpectraR04_Eta05_3040_0to6"]->fill(log(jetsAKTR04[0].pT()/p.pT())); // Figure 14
						}
					}
                                        else if(inRange(GetJetPtCorr(jetsAKTR04[0], rho04), 40., 60.)) {
                                                _c["sow4060UESub"]->fill();
						FillRadialDistribution(_p["pTR04_Eta05_4060"], jetsAKTR04[0]); // Figure 9
                                                FillUEDist(ALICEparticles, 0.4, jetsAKTR04[0], _h["pTSpectraR04_Eta05_4060_ALICEvsMC_UE"], _h["pTSpectraR04_Eta05_4060_0to1_UE"], _h["pTSpectraR04_Eta05_4060_0to6_UE"]);
						for (auto p : jetsAKTR04[0].particles()) {
							_h["pTSpectraR04_Eta05_4060_ALICEvsMC"]->fill(p.pT()/GeV); // Figure 12
                                                        _h["pTSpectraR04_Eta05_4060_0to1"]->fill(p.pT()/jetsAKTR04[0].pT()); // Figure 13
                                                        _h["pTSpectraR04_Eta05_4060_0to6"]->fill(log(jetsAKTR04[0].pT()/p.pT())); // Figure 14
						}
					}
                                        else if(inRange(GetJetPtCorr(jetsAKTR04[0], rho04), 60., 80.)) {
                                                _c["sow6080UESub"]->fill();
						FillRadialDistribution(_p["pTR04_Eta05_6080"], jetsAKTR04[0]); // Figure 9
                                                FillUEDist(ALICEparticles, 0.4, jetsAKTR04[0], _h["pTSpectraR04_Eta05_6080_ALICEvsMC_UE"], _h["pTSpectraR04_Eta05_6080_0to1_UE"], _h["pTSpectraR04_Eta05_6080_0to6_UE"]);
						for (auto p : jetsAKTR04[0].particles()) {
							_h["pTSpectraR04_Eta05_6080_ALICEvsMC"]->fill(p.pT()/GeV); // Figure 12
                                                        _h["pTSpectraR04_Eta05_6080_0to1"]->fill(p.pT()/jetsAKTR04[0].pT()); // Figure 13
                                                        _h["pTSpectraR04_Eta05_6080_0to6"]->fill(log(jetsAKTR04[0].pT()/p.pT())); // Figure 14
						}
					}

                                        if(inRange(jetsAKTR04[0].pT(), 20., 30.)) {
						_c["sow2030"]->fill();
						FillRadialDistribution(_p["pTR04_Eta05_2030_WithoutUESub"], jetsAKTR04[0]); // Figure A.4
						for (auto p : jetsAKTR04[0].particles()) {
							_h["pTSpectraR04_Eta05_2030_ALICEvsMC_WithoutUESub"]->fill(p.pT()/GeV); // Figure A.6
                                                        _h["pTSpectraR04_Eta05_2030_0to1_WithoutUESub"]->fill(p.pT()/jetsAKTR04[0].pT()); // Figure A.7
                                                        _h["pTSpectraR04_Eta05_2030_0to6_WithoutUESub"]->fill(log(jetsAKTR04[0].pT()/p.pT())); // Figure A.8
						}
					}
                                        else if(inRange(jetsAKTR04[0].pT(), 30., 40.)) {
						_c["sow3040"]->fill();
						FillRadialDistribution(_p["pTR04_Eta05_3040_WithoutUESub"], jetsAKTR04[0]); // Figure A.4
						for (auto p : jetsAKTR04[0].particles()) {
							_h["pTSpectraR04_Eta05_3040_ALICEvsMC_WithoutUESub"]->fill(p.pT()/GeV); // Figure A.6
                                                        _h["pTSpectraR04_Eta05_3040_0to1_WithoutUESub"]->fill(p.pT()/jetsAKTR04[0].pT()); // Figure A.7
                                                        _h["pTSpectraR04_Eta05_3040_0to6_WithoutUESub"]->fill(log(jetsAKTR04[0].pT()/p.pT())); // Figure A.8
						}
					}
                                        else if(inRange(jetsAKTR04[0].pT(), 40., 60.)) {
						_c["sow4060"]->fill();
						FillRadialDistribution(_p["pTR04_Eta05_4060_WithoutUESub"], jetsAKTR04[0]); // Figure A.4
						for (auto p : jetsAKTR04[0].particles()) {
							_h["pTSpectraR04_Eta05_4060_ALICEvsMC_WithoutUESub"]->fill(p.pT()/GeV); // Figure A.6
                                                        _h["pTSpectraR04_Eta05_4060_0to1_WithoutUESub"]->fill(p.pT()/jetsAKTR04[0].pT()); // Figure A.7
                                                        _h["pTSpectraR04_Eta05_4060_0to6_WithoutUESub"]->fill(log(jetsAKTR04[0].pT()/p.pT())); // Figure A.8
						}
					}
                                        else if(inRange(jetsAKTR04[0].pT(), 60., 80.)) {
						_c["sow6080"]->fill();
						FillRadialDistribution(_p["pTR04_Eta05_6080_WithoutUESub"], jetsAKTR04[0]); // Figure A.4
						for (auto p : jetsAKTR04[0].particles()) {
							_h["pTSpectraR04_Eta05_6080_ALICEvsMC_WithoutUESub"]->fill(p.pT()/GeV); // Figure A.6
                                                        _h["pTSpectraR04_Eta05_6080_0to1_WithoutUESub"]->fill(p.pT()/jetsAKTR04[0].pT()); // Figure A.7
                                                        _h["pTSpectraR04_Eta05_6080_0to6_WithoutUESub"]->fill(log(jetsAKTR04[0].pT()/p.pT())); // Figure A.8
						}
					}
				}

				for(auto jet : jetsAKTR04) {
					_h["CrossSectionAntikT_R04"]->fill(GetJetPtCorr(jet, rho04)); // Figure 2
					_h["antiKTR04"]->fill(GetJetPtCorr(jet, rho04)); // Figure 3
					_h["antiKT_R04_ALICE_ATLAS"]->fill(GetJetPtCorr(jet, rho04)); // Figure 4
					_h["ALICEvsMC_R04_Eta05"]->fill(GetJetPtCorr(jet, rho04)); // Figure 5
					_h["Ratio04_Denominator"]->fill(GetJetPtCorr(jet, rho04)); // Figure 6
					_h["ALICEvsMC_R04_Eta05_WithoutUESub"]->fill(GetJetPtCorr(jet, rho04)); // Figure A.1
					_h["pTSpectraR04_Eta05_2030_0to1_WithoutUESub"]->fill(GetJetPtCorr(jet, rho04)); // Figure A.7
					_h["pTSpectraR04_Eta05_3040_0to1_WithoutUESub"]->fill(GetJetPtCorr(jet, rho04)); // Figure A.7
					_h["pTSpectraR04_Eta05_4060_0to1_WithoutUESub"]->fill(GetJetPtCorr(jet, rho04)); // Figure A.7
					_h["pTSpectraR04_Eta05_6080_0to1_WithoutUESub"]->fill(GetJetPtCorr(jet, rho04)); // Figure A.7
					_h["pTSpectraR04_Eta05_2030_0to6_WithoutUESub"]->fill(GetJetPtCorr(jet, rho04)); // Figure A.8
					_h["pTSpectraR04_Eta05_3040_0to6_WithoutUESub"]->fill(GetJetPtCorr(jet, rho04)); // Figure A.8
					_h["pTSpectraR04_Eta05_4060_0to6_WithoutUESub"]->fill(GetJetPtCorr(jet, rho04)); // Figure A.8
					_h["pTSpectraR04_Eta05_6080_0to6_WithoutUESub"]->fill(GetJetPtCorr(jet, rho04)); // Figure A.8
				}

				// KT alg. - Resolution = 0.4, Eta = 0.5 (From 0.9 - 0.4)
				FastJets jetsKTR04FJ = apply<FastJets>(event, "jetsKTR04FJ");
				jetsKTR04FJ.calc(ALICEparticles);
				Jets jetsKTR04 = jetsKTR04FJ.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.5);
				double rho04KT = 0;

				if (jetsKTR04.size() != 0) {
					rho04KT = GetRho(ALICEparticles, 0.4, jetsKTR04[0]);
				}
				for(auto jet : jetsKTR04) {
					_h["CrossSectionkT_R04"]->fill(GetJetPtCorr(jet, rho04KT)); // Figure 2
				}

 				// CONE alg. - Resolution = 0.4, Eta = 0.5 (From 0.9 - 0.4)
				FastJets jetsCONER04FJ = apply<FastJets>(event, "jetsCONER04FJ");
				jetsCONER04FJ.calc(ALICEparticles);
				Jets jetsCONER04 = jetsCONER04FJ.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.5);
                                double rho04CONE = 0;

                                if (jetsCONER04.size() != 0) {
					rho04CONE = GetRho(ALICEparticles, 0.4, jetsCONER04[0]);
				}

				for(auto jet : jetsCONER04){
					_h["CrossSectionSisCone_R04"]->fill(jet.pT()/GeV - (rho04CONE*M_PI*0.4*0.4)); // Figure 2
				}

				// Anti-KT alg. - Resolution = 0.6, Eta = 0.3 (From 0.9 - 0.6)
				FastJets jetsAKTR06FJ = apply<FastJets>(event, "jetsAKTR06FJ");
				jetsAKTR06FJ.calc(ALICEparticles);
				Jets jetsAKTR06 = jetsAKTR06FJ.jetsByPt(Cuts::pT >= 20*GeV && Cuts::abseta < 0.3);
				double rho06 = 0;

				if(jetsAKTR06.size() != 0) {
					rho06 = GetRho(ALICEparticles, 0.6, jetsAKTR06[0]);

					_p["mean_ALICEvsMC_R06_Eta_03"]->fill(GetJetPtCorr(jetsAKTR06[0], rho06), jetsAKTR06[0].particles().size()); // Figure 7
					_p["mean_ALICEvsMC_R06_Eta_03_WithoutUESub"]->fill(jetsAKTR06[0].pT()/GeV, jetsAKTR06[0].particles().size()); // Figure A.2
					FillRadius80PercPt(_p["avgpT_R06_Eta03"], jetsAKTR06[0]); // Figure 11

					if ((GetJetPtCorr(jetsAKTR06[0], rho06) >= 20*GeV) && (GetJetPtCorr(jetsAKTR06[0], rho06) < 30*GeV)) {
						FillRadialDistribution(_p["pTR06_Eta03_2030"], jetsAKTR06[0]); // Figure 10
					}
					if ((GetJetPtCorr(jetsAKTR06[0], rho06) >= 30*GeV) && (GetJetPtCorr(jetsAKTR06[0], rho06) < 40*GeV)) {
						FillRadialDistribution(_p["pTR06_Eta03_3040"], jetsAKTR06[0]); // Figure 10
					}
					if ((GetJetPtCorr(jetsAKTR06[0], rho06) >= 40*GeV) && (GetJetPtCorr(jetsAKTR06[0], rho06) < 60*GeV)) {
						FillRadialDistribution(_p["pTR06_Eta03_4060"], jetsAKTR06[0]); // Figure 10
					}
					if ((GetJetPtCorr(jetsAKTR06[0], rho06) >= 60*GeV) && (GetJetPtCorr(jetsAKTR06[0], rho06) < 80*GeV)) {
						FillRadialDistribution(_p["pTR06_Eta03_6080"], jetsAKTR06[0]); // Figure 10
					}

					if(jetsAKTR06[0].pT() >= 20*GeV && jetsAKTR06[0].pT() < 30*GeV) {
						FillRadialDistribution(_p["pTR06_Eta03_2030_WithoutUESub"], jetsAKTR06[0]); // Figure A.5
					}
					if(jetsAKTR06[0].pT() >= 30*GeV && jetsAKTR06[0].pT() < 40*GeV) {
						FillRadialDistribution(_p["pTR06_Eta03_3040_WithoutUESub"], jetsAKTR06[0]); // Figure A.5
					}
					if(jetsAKTR06[0].pT() >= 40*GeV && jetsAKTR06[0].pT() < 60*GeV) {
						FillRadialDistribution(_p["pTR06_Eta03_4060_WithoutUESub"], jetsAKTR06[0]); // Figure A.5
					}
					if(jetsAKTR06[0].pT() >= 60*GeV && jetsAKTR06[0].pT() < 80*GeV) {
						FillRadialDistribution(_p["pTR06_Eta03_6080_WithoutUESub"], jetsAKTR06[0]); // Figure A.5
					}
				}

				for(auto jet: jetsAKTR06) {
					_h["antiKTR06"]->fill(GetJetPtCorr(jet, rho06)); // Figure 3
					_h["antiKT_R06_ALICE_ATLAS"]->fill(GetJetPtCorr(jet, rho06)); // Figure 4
					_h["ALICEvsMC_R06_Eta03"]->fill(GetJetPtCorr(jet, rho06)); // Figure 5
					_h["Ratio06_Denominator"]->fill(GetJetPtCorr(jet, rho06)); // Figure 6
					_h["ALICEvsMC_R06_Eta03_WithoutUESub"]->fill(jet.pT()/GeV); // Figure A.1
				}
			}


			/// Normalise histograms etc., after the run
			void finalize() {
				// Figure 2
				_h["CrossSectionAntikT_R04"]->scaleW(crossSection()/(nanobarn*_c["sow"]->sumW()));
				_h["CrossSectionkT_R04"]->scaleW(crossSection()/(nanobarn*_c["sow"]->sumW()));
				_h["CrossSectionSisCone_R04"]->scaleW(crossSection()/(nanobarn*_c["sow"]->sumW()));

				// Figure 3
				_h["antiKTR02"]->scaleW(crossSection()/(nanobarn*1.4*_c["sow"]->sumW())); // 0.7 * 2 = 1.4
				_h["antiKTR03"]->scaleW(crossSection()/(nanobarn*1.2*_c["sow"]->sumW()));  // 0.6 * 2 = 1.2
				_h["antiKTR04"]->scaleW(crossSection()/(nanobarn*1.0*_c["sow"]->sumW()));  // 0.5 * 2 = 1
				_h["antiKTR06"]->scaleW(crossSection()/(nanobarn*0.6*_c["sow"]->sumW()));  // 0.3 * 2 = 0.6

				// Figure 4
				_h["antiKT_R04_ALICE_ATLAS"]->scaleW(crossSection()/(nanobarn*_c["sow"]->sumW()));
				_h["antiKT_R06_ALICE_ATLAS"]->scaleW(crossSection()/(nanobarn*0.6*_c["sow"]->sumW()));

				// Figure 5
				_h["ALICEvsMC_R02_Eta07"]->scaleW(crossSection()/(nanobarn*1.4*_c["sow"]->sumW()));
				_h["ALICEvsMC_R04_Eta05"]->scaleW(crossSection()/(nanobarn*_c["sow"]->sumW()));
				_h["ALICEvsMC_R06_Eta03"]->scaleW(crossSection()/(nanobarn*0.6*_c["sow"]->sumW()));

				// Figure 6
				_h["Ratio02_Numerator"]->scaleW(1/(1.4*_c["sow"]->sumW()));
				_h["Ratio04_Denominator"]->scaleW(1/_c["sow"]->sumW());
				_h["Ratio06_Denominator"]->scaleW(1/(0.6*_c["sow"]->sumW()));
				divide(_h["Ratio02_Numerator"], _h["Ratio04_Denominator"], _s["Ratio0204"]);
				divide(_h["Ratio02_Numerator"], _h["Ratio06_Denominator"], _s["Ratio0206"]);

				// Figure 7, 8, 9, 10, & 11 do not need to be normalized

				// Figure 12
                                _h["pTSpectraR04_Eta05_2030_ALICEvsMC"]->scaleW(1/_c["sow2030UESub"]->sumW());
				_h["pTSpectraR04_Eta05_3040_ALICEvsMC"]->scaleW(1/_c["sow3040UESub"]->sumW());
				_h["pTSpectraR04_Eta05_4060_ALICEvsMC"]->scaleW(1/_c["sow4060UESub"]->sumW());
				_h["pTSpectraR04_Eta05_6080_ALICEvsMC"]->scaleW(1/_c["sow6080UESub"]->sumW());
                                // Divided by 2 because there are two cones for UE calculation
                                _h["pTSpectraR04_Eta05_2030_ALICEvsMC_UE"]->scaleW(1./(2.*_c["sow2030UESub"]->sumW()));
                                _h["pTSpectraR04_Eta05_3040_ALICEvsMC_UE"]->scaleW(1./(2.*_c["sow3040UESub"]->sumW()));
                                _h["pTSpectraR04_Eta05_4060_ALICEvsMC_UE"]->scaleW(1./(2.*_c["sow4060UESub"]->sumW()));
                                _h["pTSpectraR04_Eta05_6080_ALICEvsMC_UE"]->scaleW(1./(2.*_c["sow6080UESub"]->sumW()));
                                Subtract(_h["pTSpectraR04_Eta05_2030_ALICEvsMC"], _h["pTSpectraR04_Eta05_2030_ALICEvsMC_UE"]);
                                Subtract(_h["pTSpectraR04_Eta05_3040_ALICEvsMC"], _h["pTSpectraR04_Eta05_3040_ALICEvsMC_UE"]);
                                Subtract(_h["pTSpectraR04_Eta05_4060_ALICEvsMC"], _h["pTSpectraR04_Eta05_4060_ALICEvsMC_UE"]);
                                Subtract(_h["pTSpectraR04_Eta05_6080_ALICEvsMC"], _h["pTSpectraR04_Eta05_6080_ALICEvsMC_UE"]);


				// Figure 13
				_h["pTSpectraR04_Eta05_2030_0to1"]->scaleW(1./_c["sow2030UESub"]->sumW());
				_h["pTSpectraR04_Eta05_3040_0to1"]->scaleW(1./_c["sow3040UESub"]->sumW());
				_h["pTSpectraR04_Eta05_4060_0to1"]->scaleW(1./_c["sow4060UESub"]->sumW());
				_h["pTSpectraR04_Eta05_6080_0to1"]->scaleW(1./_c["sow6080UESub"]->sumW());
                                // Divided by 2 because there are two cones for UE calculation
                                _h["pTSpectraR04_Eta05_2030_0to1_UE"]->scaleW(1./(2.*_c["sow2030UESub"]->sumW()));
				_h["pTSpectraR04_Eta05_3040_0to1_UE"]->scaleW(1./(2.*_c["sow3040UESub"]->sumW()));
				_h["pTSpectraR04_Eta05_4060_0to1_UE"]->scaleW(1./(2.*_c["sow4060UESub"]->sumW()));
				_h["pTSpectraR04_Eta05_6080_0to1_UE"]->scaleW(1./(2.*_c["sow6080UESub"]->sumW()));
                                Subtract(_h["pTSpectraR04_Eta05_2030_0to1"], _h["pTSpectraR04_Eta05_2030_0to1_UE"]);
                                Subtract(_h["pTSpectraR04_Eta05_3040_0to1"], _h["pTSpectraR04_Eta05_3040_0to1_UE"]);
                                Subtract(_h["pTSpectraR04_Eta05_4060_0to1"], _h["pTSpectraR04_Eta05_4060_0to1_UE"]);
                                Subtract(_h["pTSpectraR04_Eta05_6080_0to1"], _h["pTSpectraR04_Eta05_6080_0to1_UE"]);

				// Figure 14
				_h["pTSpectraR04_Eta05_2030_0to6"]->scaleW(1./_c["sow2030UESub"]->sumW());
				_h["pTSpectraR04_Eta05_3040_0to6"]->scaleW(1./_c["sow3040UESub"]->sumW());
				_h["pTSpectraR04_Eta05_4060_0to6"]->scaleW(1./_c["sow4060UESub"]->sumW());
				_h["pTSpectraR04_Eta05_6080_0to6"]->scaleW(1./_c["sow6080UESub"]->sumW());
                                // Divided by 2 because there are two cones for UE calculation
                                _h["pTSpectraR04_Eta05_2030_0to6_UE"]->scaleW(1./(2.*_c["sow2030UESub"]->sumW()));
				_h["pTSpectraR04_Eta05_3040_0to6_UE"]->scaleW(1./(2.*_c["sow3040UESub"]->sumW()));
				_h["pTSpectraR04_Eta05_4060_0to6_UE"]->scaleW(1./(2.*_c["sow4060UESub"]->sumW()));
				_h["pTSpectraR04_Eta05_6080_0to6_UE"]->scaleW(1./(2.*_c["sow6080UESub"]->sumW()));
                                Subtract(_h["pTSpectraR04_Eta05_2030_0to6"], _h["pTSpectraR04_Eta05_2030_0to6_UE"]);
                                Subtract(_h["pTSpectraR04_Eta05_3040_0to6"], _h["pTSpectraR04_Eta05_3040_0to6_UE"]);
                                Subtract(_h["pTSpectraR04_Eta05_4060_0to6"], _h["pTSpectraR04_Eta05_4060_0to6_UE"]);
                                Subtract(_h["pTSpectraR04_Eta05_6080_0to6"], _h["pTSpectraR04_Eta05_6080_0to6_UE"]);

				// Figure A.1
				_h["ALICEvsMC_R04_Eta05_WithoutUESub"]->scaleW(crossSection()/(nanobarn*_c["sow"]->sumW()));
				_h["ALICEvsMC_R06_Eta03_WithoutUESub"]->scaleW(crossSection()/(nanobarn*0.6*_c["sow"]->sumW()));

				// Figures A.2, A.3, A.4, & A.5 do not need to be normalized

				// Figure A.6
				_h["pTSpectraR04_Eta05_2030_ALICEvsMC_WithoutUESub"]->scaleW(1./_c["sow2030"]->sumW());
				_h["pTSpectraR04_Eta05_3040_ALICEvsMC_WithoutUESub"]->scaleW(1./_c["sow3040"]->sumW());
				_h["pTSpectraR04_Eta05_4060_ALICEvsMC_WithoutUESub"]->scaleW(1./_c["sow4060"]->sumW());
				_h["pTSpectraR04_Eta05_6080_ALICEvsMC_WithoutUESub"]->scaleW(1./_c["sow6080"]->sumW());

				// Figure A.7
				_h["pTSpectraR04_Eta05_2030_0to1_WithoutUESub"]->scaleW(1./_c["sow2030"]->sumW());
				_h["pTSpectraR04_Eta05_3040_0to1_WithoutUESub"]->scaleW(1./_c["sow3040"]->sumW());
				_h["pTSpectraR04_Eta05_4060_0to1_WithoutUESub"]->scaleW(1./_c["sow4060"]->sumW());
				_h["pTSpectraR04_Eta05_6080_0to1_WithoutUESub"]->scaleW(1./_c["sow6080"]->sumW());

				// Figure A.8
				_h["pTSpectraR04_Eta05_2030_0to6_WithoutUESub"]->scaleW(1./_c["sow2030"]->sumW());
				_h["pTSpectraR04_Eta05_3040_0to6_WithoutUESub"]->scaleW(1./_c["sow3040"]->sumW());
				_h["pTSpectraR04_Eta05_4060_0to6_WithoutUESub"]->scaleW(1./_c["sow4060"]->sumW());
				_h["pTSpectraR04_Eta05_6080_0to6_WithoutUESub"]->scaleW(1./_c["sow6080"]->sumW());

			}

			///@}


			/// @name Histograms
			///@{
			map<string, Histo1DPtr> _h;
			map<string, Profile1DPtr> _p;
			map<string, CounterPtr> _c;
			map<string, Scatter2DPtr> _s;

			fastjet::AreaDefinition *fjAreaDef02;
			fastjet::AreaDefinition *fjAreaDef03;
			fastjet::AreaDefinition *fjAreaDef04;
			fastjet::AreaDefinition *fjAreaDef04KT;
			fastjet::AreaDefinition *fjAreaDef06;
			///@}


	};


	RIVET_DECLARE_PLUGIN(ALICE_2015_I1328629);

}
