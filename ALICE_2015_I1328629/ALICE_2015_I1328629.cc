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
				FastJets jetsAKT02(fs, FastJets::ANTIKT, 0.2, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKT02, "jetsAKT02");


				FastJets jetsAKT03(fs, FastJets::ANTIKT, 0.3, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKT03, "jetsAKT03");

				FastJets jetsAKT04(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKT04, "jetsAKT04");
				FastJets jetsKT04(fs, FastJets::KT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsKT04, "jetsKT04");
				FastJets jetsCONE04(fs, FastJets::SISCONE, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsCONE04, "jetsCONE04");

				FastJets jetsAKT06(fs, FastJets::ANTIKT, 0.6, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKT06, "jetsAKT06");

				
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
				book(_h["antiKT_R04_06_ALICE_ATLAS"], 8, 1, 1); // Figure 4
				book(_h["antiKT_ALICEvsMCAlgs_R02_04_06_Eta07_05_03"], 10, 1, 1);  // Figure 5
				book(_h["ratio_AntiKT_R0204_0206_Eta03"], 13, 1, 1);  // Figure 6
				book(_h["mean_ALICEvsMC_R06_04_02_Eta_03_05_07"], 15, 1, 1);  // Figure 7
				book(_h["pTR02_Eta07Bins_2030_3040_4060_6080"], 18, 1, 1); // Figure 8
				book(_h["pTR04_Eta05Bins_2030_3040_4060_6080"], 22, 1, 1);  // Figure 9
				book(_h["pTR04_Eta03Bins_2030_3040_4060_6080"], 26, 1, 1);  // Figure 10
				book(_h["avgR80pT_R02_04_06_Eta07_05_03"], 30, 1, 1);  // Figure 11
				book(_h["pTSpectraR04_Eta05Bins_ALICEvsMC"], 33, 1, 1);  // Figure 12
				book(_h["pTSpectraR04_Eta05Bins0to1"], 37, 1, 1);  // Figure 13
				book(_h["pTSpectraR04_Eta05Bins0to6"], 41, 1, 1);  // Figure 14
			}


			/// Perform the per-event analysis
			void analyze(const Event& event) {
				//Increment the counter to keep track of the cross section
				_c["sow"]->fill();

				// Retrieve clustered jets, sorted by pT, with a minimum pT cut
				// Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

				const ALICE::PrimaryParticles aprim = apply<ALICE::PrimaryParticles>(event, "aprim");
				const Particles ALICEparticles = aprim.particles();


				// Jets jetsAKT02 = apply<FastJets>(event, "jetsAKT02").jetsByPt(Cuts::pT > 0.15*GeV && Cuts::abseta < 0.9);
				// Jets jetsKT02 = apply<FastJets>(event, "jetsKT02").jetsByPt(Cuts::pT > 0.15*GeV && Cuts::abseta < 0.9);
				// Jets jetsCONE02 = apply<FastJets>(event, "jetsCONE02").jetsByPt(Cuts::pT > 0.15*GeV && Cuts::abseta < 0.9);

				 Jets jetsAKT04FJ = apply<FastJets>(event, "jetsAKT04").jetsByPt(Cuts::pT > 0.00*GeV && Cuts::abseta < 0.9);
				 //The next two lines are only required to trick Rivet into using primary particles
				//jetsAKT04FJ.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
				//Jets jetsAKT04 = jetsAKT04FJ.jetsByPt(Cuts::pT >= 20.*GeV); //get jets (ordered by pT)
				for(auto jet : jetsAKT04FJ){_h["CrossSectionAntikT_R04"]->fill(jet.pT()/GeV);}


				 Jets jetsKT04FJ = apply<FastJets>(event, "jetsKT04").jetsByPt(Cuts::pT > 0.00*GeV && Cuts::abseta < 0.9);
				// jetsKT04FJ.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
				// Jets jetsKT04 = jetsKT04FJ.jetsByPt(Cuts::pT >= 20.*GeV); //get jets (ordered by pT)
				for(auto jet : jetsKT04FJ){_h["CrossSectionkT_R04"]->fill(jet.pT()/GeV);}

				 Jets jetsCONE04FJ = apply<FastJets>(event, "jetsCONE04").jetsByPt(Cuts::pT > 0.00*GeV && Cuts::abseta < 0.9);
				// jetsCONE04FJ.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
				 //Jets jetsCONE04 = jetsCONE04FJ.jetsByPt(Cuts::pT >= 20.*GeV); //get jets (ordered by pT)
				for(auto jet : jetsCONE04FJ){_h["CrossSectionSisCone_R04"]->fill(jet.pT()/GeV);}

				// Jets jetsAKT06 = apply<FastJets>(event, "jetsAKT06").jetsByPt(Cuts::pT > 0.15*GeV && Cuts::abseta < 0.9);
				// Jets jetsKT06 = apply<FastJets>(event, "jetsKT06").jetsByPt(Cuts::pT > 0.15*GeV && Cuts::abseta < 0.9);
				// Jets jetsCONE06 = apply<FastJets>(event, "jetsCONE06").jetsByPt(Cuts::pT > 0.15*GeV && Cuts::abseta < 0.9);
				
				// FJjets.calc(ALICEparticles); //give ALICE primary particles to FastJet projection

				// "sow" - Sum of Weights - Adding up equivalent proton-proton collisions that occurred
				// _c["sow"]->fill();
				// for(auto jet : jetsAKT04) { _h["all3JetFinders_R04"]->fill(jet.pT()/GeV); }

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
