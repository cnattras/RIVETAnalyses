// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include <stdio.h>
#include <vector>
#include <string>

namespace Rivet {


	/// @brief Add a short analysis description here
	class ALICE_2015_I1328629 : public Analysis {
		public:

			/// Constructor
			RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2015_I1328629);

			/// @name Analysis methods
			///@{
			vector <string> histogramNames;

			/// Book histograms and initialise projections before the run
			void init() {

				// Initialise and register projections

				// The basic final-state projection: all final-state particles within the given eta acceptance
				const FinalState fs(Cuts::abseta < 1.5);

				/* The final-state particles declared above are clustered using FastJet with the anti-kT 
				   algorithm and a jet-radius parameter 0.4 muons & neutrinos are excluded from the clustering */
				FastJets jetsAKT02(fs, FastJets::ANTIKT, 0.2, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKT02, "jetsAKT02");
				FastJets jetsKT02(fs, FastJets::KT, 0.2, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsKT02, "jetsKT02");
				FastJets jetsCONE02(fs, FastJets::SISCONE, 0.2, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsCONE02, "jetsCONE02");

				FastJets jetsAKT04(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKT04, "jetsAKT04");
				FastJets jetsKT04(fs, FastJets::KT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsKT04, "jetsKT04");
				FastJets jetsCONE04(fs, FastJets::SISCONE, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsCONE04, "jetsCONE04");

				FastJets jetsAKT06(fs, FastJets::ANTIKT, 0.6, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKT06, "jetsAKT06");
				FastJets jetsKT06(fs, FastJets::KT, 0.6, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsKT06, "jetsKT06");
				FastJets jetssCONE06(fs, FastJets::SISCONE, 0.6, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetssCONE06, "jetssCONE06");

				
				const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::abscharge > 0);
				declare(aprim, "aprim");

				book(_c["sow"], "sow");

				// Creates histograms for all Tables in .yoda file (73 tables in total, all in order)
				char name[100];
				for (int i = 1; i < 74; i++) {
					snprintf(name, sizeof(name), "Table%d", i);
					histogramNames.push_back(name)  // Keeps track of all table names
					book(_h[name], i, 1, 1);
				}

			}


			/// Perform the per-event analysis
			void analyze(const Event& event) {

				// Retrieve clustered jets, sorted by pT, with a minimum pT cut
				// Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

				const ALICE::PrimaryParticles aprim = apply<ALICE::PrimaryParticles>(event, "aprim");
				const Particles ALICEparticles = aprim.particles();

				Jets jetsAKT04 = apply<FastJets>(event, "jetsAKT04").jetsByPt(Cuts::pT >= 0*GeV && Cuts::abseta <0.5);
				Jets jetsKT04 = apply<FastJets>(event, "jetsKT04").jetsByPt(Cuts::pT >= 0*GeV && Cuts::abseta <0.5);
				Jets jetsCONE04 = apply<FastJets>(event, "jetsCONE04").jetsByPt(Cuts::pT >= 0*GeV && Cuts::abseta <0.5);
				
				// FJjets.calc(ALICEparticles); //give ALICE primary particles to FastJet projection

				// "sow" - Sum of Weights - Adding up equivalent proton-proton collisions that occurred
				// _c["sow"]->fill();
				// for(auto jet : jetsAKT04)
				// {
				// 	_h["Figure1"]->fill(jet.pT()/GeV);
				// }
			}


			/// Normalise histograms etc., after the run
			void finalize() {
				// _h["Figure1"]->scaleW(1./_c["sow"]->sumW());
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
