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

				// The basic final-state projection: all final-state particles within the given eta acceptance
				const FinalState fs(Cuts::abseta < 0.75);

				/*
					The final-state particles declared above are clustered using FastJet with the anti-kT 
					algorithm and a jet-radius parameter 0.4 muons & neutrinos are excluded from the clustering 
				*/
				FastJets jetsAKTR02FJ(fs, FastJets::ANTIKT, 0.2, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKTR02FJ, "jetsAKTR02FJ");

				FastJets jetsAKTR03FJ(fs, FastJets::ANTIKT, 0.3, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKTR03FJ, "jetsAKTR03FJ");

				FastJets jetsAKTR04FJ(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKTR04FJ, "jetsAKTR04FJ");

				FastJets jetsKTR04FJ(fs, FastJets::KT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsKTR04FJ, "jetsKTR04FJ");

				FastJets jetsCONER04FJ(fs, FastJets::SISCONE, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsCONER04FJ, "jetsCONER04FJ");

				FastJets jetsAKTR06FJ(fs, FastJets::ANTIKT, 0.6, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
				declare(jetsAKTR06FJ, "jetsAKTR06FJ");

				// Initialize ALICE primary particles
				const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::abscharge > 0);
				declare(aprim, "aprim");

				// Create counter
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

				string refname0204 = mkAxisCode(13, 1, 1);   // Figure 6
				const Scatter2D& refdata0204 = refData(refname0204);
				book(_h["Ratio02_Numerator"], refname0204 + "Ratio02_Numerator", refdata0204);
				book(_h["Ratio04_Denominator"], refname0204 + "Ratio04_Denominator", refdata0204);
				book(_s["Ratio0204"], refname0204);
				string refname0206 = mkAxisCode(14, 1, 1);   // Figure 6
				const Scatter2D& refdata0206 = refData(refname0206);
				book(_h["Ratio06_Denominator"], refname0206 + "Ratio06_Denominator", refdata0206);
				book(_s["Ratio0206"], refname0206);				
				
				book(_p["mean_ALICEvsMC_R02_Eta_07"], 15, 1, 1);  // Figure 7
				book(_p["mean_ALICEvsMC_R04_Eta_05"], 16, 1, 1);  // Figure 7
				book(_p["mean_ALICEvsMC_R06_Eta_03"], 17, 1, 1);  // Figure 7

				book(_h["pTR02_Eta07_2030"], 18, 1, 1);  // Figure 8
				book(_h["pTR02_Eta07_3040"], 19, 1, 1);  // Figure 8
				book(_h["pTR02_Eta07_4060"], 20, 1, 1);  // Figure 8
				book(_h["pTR02_Eta07_6080"], 21, 1, 1);  // Figure 8

				book(_h["pTR04_Eta05_2030"], 22, 1, 1);  // Figure 9
				book(_h["pTR04_Eta05_3040"], 23, 1, 1);  // Figure 9
				book(_h["pTR04_Eta05_4060"], 24, 1, 1);  // Figure 9
				book(_h["pTR04_Eta05_6080"], 25, 1, 1);  // Figure 9

				book(_h["pTR06_Eta03_2030"], 26, 1, 1);  // Figure 10
				book(_h["pTR06_Eta03_3040"], 27, 1, 1);  // Figure 10
				book(_h["pTR06_Eta03_4060"], 28, 1, 1);  // Figure 10
				book(_h["pTR06_Eta03_6080"], 29, 1, 1);  // Figure 10

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
				FastJets jetsAKTR02FJ = apply<FastJets>(event, "jetsAKTR02FJ");
				jetsAKTR02FJ.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
				Jets jetsAKTR02 = jetsAKTR02FJ.jetsByPt(Cuts::pT >= 20*GeV && Cuts::abseta < 0.7);

				if(jetsAKTR02.size() != 0)
					_p["mean_ALICEvsMC_R02_Eta_07"]->fill(jetsAKTR02[0].pT()/GeV, jetsAKTR02[0].particles().size());

				for(auto jet: jetsAKTR02) {
					_h["antiKTR02"]->fill(jet.pT()/GeV); // Figure 3
					_h["ALICEvsMC_R02_Eta07"]->fill(jet.pT()/GeV); // Figure 5
					_h["Ratio02_Numerator"]->fill(jet.pT()/GeV); // Figure 6
					// _h["mean_ALICEvsMC_R02_Eta_07"]->fill(jet.pT()/GeV); // Figure 7
					_h["pTR02_Eta07_2030"]->fill(jet.pT()/GeV); // Figure 8
					_h["pTR02_Eta07_3040"]->fill(jet.pT()/GeV); // Figure 8
					_h["pTR02_Eta07_4060"]->fill(jet.pT()/GeV); // Figure 8
					_h["pTR02_Eta07_6080"]->fill(jet.pT()/GeV); // Figure 8
					_h["avgpT_R02_Eta07"]->fill(jet.pT()/GeV); // Figure 11
				}

				// Anti-KT alg. - Resolution = 0.3, Eta = 0.6 (From 0.9 - 0.3)
				FastJets jetsAKTR03FJ = apply<FastJets>(event, "jetsAKTR03FJ");
				jetsAKTR03FJ.calc(ALICEparticles); 
				Jets jetsAKTR03 = jetsAKTR03FJ.jetsByPt(Cuts::pT >= 20*GeV && Cuts::abseta < 0.6);
				for(auto jet: jetsAKTR03) {
					_h["antiKTR03"]->fill(jet.pT()/GeV); // Figure 3

				}

				// Anti-KT alg. - Resolution = 0.4, Eta = 0.5 (From 0.9 - 0.4)
				FastJets jetsAKTR04FJ = apply<FastJets>(event, "jetsAKTR04FJ");
				jetsAKTR04FJ.calc(ALICEparticles); 
				Jets jetsAKTR04 = jetsAKTR04FJ.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.5); 

				if(jetsAKTR04.size() != 0)
					_p["mean_ALICEvsMC_R04_Eta_05"]->fill(jetsAKTR04[0].pT()/GeV, jetsAKTR04[0].particles().size());

				for(auto jet : jetsAKTR04) {
					_h["CrossSectionAntikT_R04"]->fill(jet.pT()/GeV); // Figure 2
					_h["antiKTR04"]->fill(jet.pT()/GeV); // Figure 3
					_h["antiKT_R04_ALICE_ATLAS"]->fill(jet.pT()/GeV); // Figure 4
					_h["ALICEvsMC_R04_Eta05"]->fill(jet.pT()/GeV); // Figure 5
					_h["Ratio04_Denominator"]->fill(jet.pT()/GeV); // Figure 6
					// _h["mean_ALICEvsMC_R04_Eta_05"]->fill(jet.pT()/GeV); // Figure 7
					_h["pTR04_Eta05_2030"]->fill(jet.pT()/GeV); // Figure 9
					_h["pTR04_Eta05_3040"]->fill(jet.pT()/GeV); // Figure 9
					_h["pTR04_Eta05_4060"]->fill(jet.pT()/GeV); // Figure 9
					_h["pTR04_Eta05_6080"]->fill(jet.pT()/GeV); // Figure 9
					_h["avgpT_R04_Eta05"]->fill(jet.pT()/GeV); // Figure 11
					_h["pTSpectraR04_Eta05_2030_ALICEvsMC"]->fill(jet.pT()/GeV); // Figure 12
					_h["pTSpectraR04_Eta05_3040_ALICEvsMC"]->fill(jet.pT()/GeV); // Figure 12
					_h["pTSpectraR04_Eta05_4060_ALICEvsMC"]->fill(jet.pT()/GeV); // Figure 12
					_h["pTSpectraR04_Eta05_6080_ALICEvsMC"]->fill(jet.pT()/GeV); // Figure 12
					_h["pTSpectraR04_Eta05_2030_0to1"]->fill(jet.pT()/GeV); // Figure 13
					_h["pTSpectraR04_Eta05_3040_0to1"]->fill(jet.pT()/GeV); // Figure 13
					_h["pTSpectraR04_Eta05_4060_0to1"]->fill(jet.pT()/GeV); // Figure 13
					_h["pTSpectraR04_Eta05_6080_0to1"]->fill(jet.pT()/GeV); // Figure 13
					_h["pTSpectraR04_Eta05_2030_0to6"]->fill(jet.pT()/GeV); // Figure 14
					_h["pTSpectraR04_Eta05_3040_0to6"]->fill(jet.pT()/GeV); // Figure 14
					_h["pTSpectraR04_Eta05_4060_0to6"]->fill(jet.pT()/GeV); // Figure 14
					_h["pTSpectraR04_Eta05_6080_0to6"]->fill(jet.pT()/GeV); // Figure 14

				}

				// KT alg. - Resolution = 0.4, Eta = 0.5 (From 0.9 - 0.4)
				FastJets jetsKTR04FJ = apply<FastJets>(event, "jetsKTR04FJ");
				jetsKTR04FJ.calc(ALICEparticles); 
				Jets jetsKTR04 = jetsKTR04FJ.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.5);
				for(auto jet : jetsKTR04) {
					_h["CrossSectionkT_R04"]->fill(jet.pT()/GeV); // Figure 2
				}

 				// CONE alg. - Resolution = 0.4, Eta = 0.5 (From 0.9 - 0.4)
				FastJets jetsCONER04FJ = apply<FastJets>(event, "jetsCONER04FJ");
				jetsCONER04FJ.calc(ALICEparticles); 
				Jets jetsCONER04 = jetsCONER04FJ.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.5); 
				for(auto jet : jetsCONER04){
					_h["CrossSectionSisCone_R04"]->fill(jet.pT()/GeV); // Figure 2
				}

				// Anti-KT alg. - Resolution = 0.6, Eta = 0.3 (From 0.9 - 0.6)
				FastJets jetsAKTR06FJ = apply<FastJets>(event, "jetsAKTR06FJ");
				jetsAKTR06FJ.calc(ALICEparticles);
				Jets jetsAKTR06 = jetsAKTR06FJ.jetsByPt(Cuts::pT >= 20*GeV && Cuts::abseta < 0.3);
				
				if(jetsAKTR06.size() != 0)
					_p["mean_ALICEvsMC_R06_Eta_03"]->fill(jetsAKTR06[0].pT()/GeV, jetsAKTR06[0].particles().size());

				for(auto jet: jetsAKTR06) {
					_h["antiKTR06"]->fill(jet.pT()/GeV); // Figure 3
					_h["antiKT_R06_ALICE_ATLAS"]->fill(jet.pT()/GeV); // Figure 4
					_h["ALICEvsMC_R06_Eta03"]->fill(jet.pT()/GeV); // Figure 5
					_h["Ratio06_Denominator"]->fill(jet.pT()/GeV); // Figure 6
					// _h["mean_ALICEvsMC_R06_Eta_03"]->fill(jet.pT()/GeV); // Figure 7
					_h["pTR06_Eta03_2030"]->fill(jet.pT()/GeV); // Figure 10
					_h["pTR06_Eta03_3040"]->fill(jet.pT()/GeV); // Figure 10
					_h["pTR06_Eta03_4060"]->fill(jet.pT()/GeV); // Figure 10
					_h["pTR06_Eta03_6080"]->fill(jet.pT()/GeV); // Figure 10
					_h["avgpT_R06_Eta03"]->fill(jet.pT()/GeV); // Figure 11
				}
			}


			/// Normalise histograms etc., after the run
			void finalize() {
				// Figure 2
				_h["CrossSectionAntikT_R04"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["CrossSectionkT_R04"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["CrossSectionSisCone_R04"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));

				// Figure 3
				_h["antiKTR02"]->scaleW(crossSection()/(millibarn*1.4*_c["sow"]->sumW())); // 0.7 * 2 = 1.4
				_h["antiKTR03"]->scaleW(crossSection()/(millibarn*1.2*_c["sow"]->sumW()));  // 0.6 * 2 = 1.2
				_h["antiKTR04"]->scaleW(crossSection()/(millibarn*1.0*_c["sow"]->sumW()));  // 0.5 * 2 = 1
				_h["antiKTR06"]->scaleW(crossSection()/(millibarn*0.6*_c["sow"]->sumW()));  // 0.3 * 2 = 0.6

				// Figure 4
				_h["antiKT_R04_ALICE_ATLAS"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW())); 
				_h["antiKT_R06_ALICE_ATLAS"]->scaleW(crossSection()/(millibarn*0.6*_c["sow"]->sumW()));

				// Figure 5
				_h["ALICEvsMC_R02_Eta07"]->scaleW(crossSection()/(millibarn*1.4*_c["sow"]->sumW()));
				_h["ALICEvsMC_R04_Eta05"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["ALICEvsMC_R06_Eta03"]->scaleW(crossSection()/(millibarn*0.6*_c["sow"]->sumW()));

				// Figure 6
				_h["Ratio02_Numerator"]->scaleW(1/(1.4*_c["sow"]->sumW()));
				_h["Ratio04_Denominator"]->scaleW(1/_c["sow"]->sumW());
				_h["Ratio06_Denominator"]->scaleW(1/(0.6*_c["sow"]->sumW()));
				divide(_h["Ratio02_Numerator"], _h["Ratio04_Denominator"], _s["Ratio0204"]);
				divide(_h["Ratio02_Numerator"], _h["Ratio06_Denominator"], _s["Ratio0206"]);

				// Figure 7 does not need to be normalized

				// Figure 8
				_h["pTR02_Eta07_2030"]->scaleW(crossSection()/(millibarn*1.4*_c["sow"]->sumW()));
				_h["pTR02_Eta07_3040"]->scaleW(crossSection()/(millibarn*1.4*_c["sow"]->sumW()));
				_h["pTR02_Eta07_4060"]->scaleW(crossSection()/(millibarn*1.4*_c["sow"]->sumW()));
				_h["pTR02_Eta07_6080"]->scaleW(crossSection()/(millibarn*1.4*_c["sow"]->sumW()));

				// Figure 9
				_h["pTR04_Eta05_2030"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTR04_Eta05_3040"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTR04_Eta05_4060"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTR04_Eta05_6080"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));

				// Figure 10
				_h["pTR06_Eta03_2030"]->scaleW(crossSection()/(millibarn*0.6*_c["sow"]->sumW()));
				_h["pTR06_Eta03_3040"]->scaleW(crossSection()/(millibarn*0.6*_c["sow"]->sumW()));
				_h["pTR06_Eta03_4060"]->scaleW(crossSection()/(millibarn*0.6*_c["sow"]->sumW()));
				_h["pTR06_Eta03_6080"]->scaleW(crossSection()/(millibarn*0.6*_c["sow"]->sumW()));

				// Figure 11 
				_h["avgpT_R02_Eta07"]->scaleW(crossSection()/(millibarn*1.4*_c["sow"]->sumW()));
				_h["avgpT_R04_Eta05"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["avgpT_R06_Eta03"]->scaleW(crossSection()/(millibarn*0.6*_c["sow"]->sumW()));

				// Figure 12
				_h["pTSpectraR04_Eta05_2030_ALICEvsMC"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTSpectraR04_Eta05_3040_ALICEvsMC"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTSpectraR04_Eta05_4060_ALICEvsMC"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTSpectraR04_Eta05_6080_ALICEvsMC"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));

				// Figure 13
				_h["pTSpectraR04_Eta05_2030_0to1"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTSpectraR04_Eta05_3040_0to1"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTSpectraR04_Eta05_4060_0to1"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTSpectraR04_Eta05_6080_0to1"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));

				// Figure 14
				_h["pTSpectraR04_Eta05_2030_0to6"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTSpectraR04_Eta05_3040_0to6"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTSpectraR04_Eta05_4060_0to6"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
				_h["pTSpectraR04_Eta05_6080_0to6"]->scaleW(crossSection()/(millibarn*_c["sow"]->sumW()));
			}

			///@}


			/// @name Histograms
			///@{
			map<string, Histo1DPtr> _h;
			map<string, Profile1DPtr> _p;
			map<string, CounterPtr> _c;
			map<string, Scatter2DPtr> _s;
			///@}


	};


	RIVET_DECLARE_PLUGIN(ALICE_2015_I1328629);

}
