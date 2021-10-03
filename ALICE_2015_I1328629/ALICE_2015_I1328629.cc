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

			/// Book histograms and initialise projections before the run
			void init() {

				// Initialise and register projections

				// The basic final-state projection: all final-state particles within the given eta acceptance
				const FinalState fs(Cuts::abseta < 1.5);

				/* The final-state particles declared above are clustered using FastJet with the anti-kT 
				   algorithm and a jet-radius parameter 0.4 muons & neutrinos are excluded from the clustering */
				FastJets jetsAKT02FJ(fs, FastJets::ANTIKT, 0.2, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKT02FJ, "jetsAKT02FJ");

				FastJets jetsAKT03FJ(fs, FastJets::ANTIKT, 0.3, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKT03FJ, "jetsAKT03FJ");

				FastJets jetsAKTR04E05FJ(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKTR04E05FJ, "jetsAKTR04E05FJ");
				FastJets jetsKT04FJ(fs, FastJets::KT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsKT04FJ, "jetsKT04FJ");
				FastJets jetsCONE04FJ(fs, FastJets::SISCONE, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsCONE04FJ, "jetsCONE04FJ");

				FastJets jetsAKT06FJ(fs, FastJets::ANTIKT, 0.6, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKT06FJ, "jetsAKT06FJ");

				
				const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::abscharge > 0);
				declare(aprim, "aprim");

				book(_c["sow"], "sow");

				// Creates histograms 
				book(_h["CrossSectionAntikT_R04"], 1, 1, 1);  // Figure 2 Anti-kT
				book(_h["CrossSectionkT_R04"], 2, 1, 1);  // Figure 2 kT
				book(_h["CrossSectionSisCone_R04"], 3, 1, 1);  // Figure 2 SisCone

				book(_h["antiKTR02"], 4, 1, 1);  // Figure 3
				book(_h["antiKTR03"], 5, 1, 1);  // Figure 3
				book(_h["antiKTR04"], 6, 1, 1);  // Figure 3
				book(_h["antiKTR06"], 7, 1, 1);  // Figure 3

				book(_h["antiKT_R04_ALICE_ATLAS"], 8, 1, 1); // Figure 4
				book(_h["antiKT_R06_ALICE_ATLAS"], 9, 1, 1); // Figure 4

				book(_h["ALICEvsMC_R02_Eta07"], 10, 1, 1);  // Figure 5
				book(_h["ALICEvsMC_R04_Eta05"], 11, 1, 1);  // Figure 5
				book(_h["ALICEvsMC_R06_Eta03"], 12, 1, 1);  // Figure 5

				book(_h["ratio_AntiKT_R0204"], 13, 1, 1);  // Figure 6
				book(_h["ratio_AntiKT_R0206"], 14, 1, 1);  // Figure 6

				book(_h["mean_ALICEvsMC_R02_Eta_07"], 15, 1, 1);  // Figure 7
				book(_h["mean_ALICEvsMC_R04_Eta_05"], 16, 1, 1);  // Figure 7
				book(_h["mean_ALICEvsMC_R06_Eta_03"], 17, 1, 1);  // Figure 7

				book(_h["pTR02_Eta07_2030"], 18, 1, 1);  // Figure 8
				book(_h["pTR02_Eta07_3040"], 19, 1, 1);  // Figure 8
				book(_h["pTR02_Eta07_4060"], 20, 1, 1);  // Figure 8
				book(_h["pTR02_Eta07_6080"], 21, 1, 1);  // Figure 8

				book(_h["pTR04_Eta05_2030"], 22, 1, 1);  // Figure 9
				book(_h["pTR04_Eta05_3040"], 23, 1, 1);  // Figure 9
				book(_h["pTR04_Eta05_4060"], 24, 1, 1);  // Figure 9
				book(_h["pTR04_Eta05_6080"], 25, 1, 1);  // Figure 9

				book(_h["pTR04_Eta03_2030"], 26, 1, 1);  // Figure 10
				book(_h["pTR04_Eta03_3040"], 27, 1, 1);  // Figure 10
				book(_h["pTR04_Eta03_4060"], 28, 1, 1);  // Figure 10
				book(_h["pTR04_Eta03_6080"], 29, 1, 1);  // Figure 10

				book(_h["avgpT_R02_Eta07"], 30, 1, 1);  // Figure 11
				book(_h["avgpT_R04_Eta05"], 31, 1, 1);  // Figure 11
				book(_h["avgpT_R06_Eta03"], 32, 1, 1);  // Figure 11

				book(_h["pTSpectraR04_Eta05_2030_ALICEvsMC"], 33, 1, 1);  // Figure 12
				book(_h["pTSpectraR04_Eta05_3040_ALICEvsMC"], 34, 1, 1);  // Figure 12
				book(_h["pTSpectraR04_Eta05_4060_ALICEvsMC"], 35, 1, 1);  // Figure 12
				book(_h["pTSpectraR04_Eta05_6080_ALICEvsMC"], 36, 1, 1);  // Figure 12

				book(_h["pTSpectraR04_Eta05_2030_0to1"], 37, 1, 1);  // Figure 13
				book(_h["pTSpectraR04_Eta05_3040_0to1"], 38, 1, 1);  // Figure 13
				book(_h["pTSpectraR04_Eta05_4060_0to1"], 39, 1, 1);  // Figure 13
				book(_h["pTSpectraR04_Eta05_6080_0to1"], 40, 1, 1);  // Figure 13

				book(_h["pTSpectraR04_Eta05_2030_0to6"], 41, 1, 1);  // Figure 14
				book(_h["pTSpectraR04_Eta05_3040_0to6"], 42, 1, 1);  // Figure 14
				book(_h["pTSpectraR04_Eta05_4060_0to6"], 43, 1, 1);  // Figure 14
				book(_h["pTSpectraR04_Eta05_6080_0to6"], 44, 1, 1);  // Figure 14

			}


			/// Perform the per-event analysis
			void analyze(const Event& event) {
				//Increment the counter to keep track of the cross section
				_c["sow"]->fill();

				// Retrieve clustered jets, sorted by pT, with a minimum pT cut
				const ALICE::PrimaryParticles aprim = apply<ALICE::PrimaryParticles>(event, "aprim");
				const Particles ALICEparticles = aprim.particles();

				// Anti-KT alg. - Resolution = 0.2, Eta = 0.7 (From 0.9 - 0.2)
				FastJets jetsAKTR02E07FJ = apply<FastJets>(event, "jetsAKTR02E07FJ");
				jetsAKTR02E07FJ.calc(ALICEparticles);
				Jets jetsAKTR02E07 = jetsAKTR02E07FJ.jetsByPt(Cuts::pT >= 20*GeV && Cuts::abseta < 0.7);
				for(auto jet: jetsAKTR02E07FJ) {

				}

				// Anti-KT alg. - Resolution = 0.4, Eta = 0.5 (From 0.9 - 0.4)
				FastJets jetsAKTR04E05FJ = apply<FastJets>(event, "jetsAKTR04E05FJ");
				jetsAKTR04E05FJ.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
				Jets jetsAKTR04E05 = jetsAKTR04E05FJ.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.5); // abseta = 0.9 - R
				for(auto jet : jetsAKTR04E05FJ) {
					_h["CrossSectionAntikT_R04"]->fill(jet.pT()/GeV);
				}

				// KT alg. - Resolution = 0.4, Eta = 0.5 (From 0.9 - 0.4)
				FastJets jetsKT04FJ = apply<FastJets>(event, "jetsKT04FJ");
				jetsKT04FJ.calc(ALICEparticles); 
				Jets jetsKT04 = jetsKT04FJ.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.5);
				for(auto jet : jetsKT04FJ) {
					_h["CrossSectionkT_R04"]->fill(jet.pT()/GeV);
				}

 				// CONE alg. - Resolution = 0.4, Eta = 0.5 (From 0.9 - 0.4)
				FastJets jetsCONE04FJ = apply<FastJets>(event, "jetsCONE04FJ");
				jetsCONE04FJ.calc(ALICEparticles); 
				Jets jetsCONE04 = jetsCONE04FJ.jetsByPt(Cuts::pT >= 20.*GeV); 
				for(auto jet : jetsCONE04FJ){_h["CrossSectionSisCone_R04"]->fill(jet.pT()/GeV);}

				// Anti-KT alg. - Resolution = 0.6, Eta = 0.3 (From 0.9 - 0.6)
				FastJets jetsAKTR06E03FJ = apply<FastJets>(event, "jetsAKTR06E03FJ");
				jetsAKTR06E03FJ.calc(ALICEparticles);
				Jets jetsAKTR06E03 = jetsAKTR06E03FJ.jetsByPt(Cuts::pT >= 20*GeV && Cuts::abseta < 0.3);
				for(auto jet: jetsAKTR06E03FJ) {
					
				}
			}


			/// Normalise histograms etc., after the run
			void finalize() {
				 _h["CrossSectionAntikT_R04"]->scaleW(crossSection()/_c["sow"]->sumW());
				 _h["CrossSectionkT_R04"]->scaleW(crossSection()/_c["sow"]->sumW());
				 _h["CrossSectionSisCone_R04"]->scaleW(crossSection()/_c["sow"]->sumW());
			}

			///@}


			/// @name Histograms
			///@{
			map<string, Histo1DPtr> _h;
			map<string, Profile1DPtr> _p;
			map<string, CounterPtr> _c;
			///@}


	};


	RIVET_DECLARE_PLUGIN(ALICE_2015_I1328629);

}
