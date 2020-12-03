// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "RHICCentrality.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class STAR_2010_I837075 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2010_I837075);

    /// Book histograms and initialise projections before the run
    void init() {
        w2 = 0;
        w_all = 0;

        declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

        // ----------------------------------------------------------
        // Cuts used to see if there is a trigger hit in the BBC
        // Refer to example/precedent from: 
        //   https://rivet.hepforge.org/analyses/ALICE_2010_I880049
        // ----------------------------------------------------------
        declare(ChargedFinalState((Cuts::eta >= 3.3 && Cuts::eta < 5.0) &&
                    Cuts::pT > 0.1*GeV), "BBC_East");
        declare(ChargedFinalState((Cuts::eta > -5.0 && Cuts::eta < -3.3) &&
                    Cuts::pT > 0.1*GeV), "BBC_West");
        // note: this cut probably isn't very restrictive in Cu+Cu, but more
        //       often in pp

        // ----------------------------------------------------------
        // Cuts used to get final state pi+, pi-, proton, pbar
        // Refer to example/precedent from: 
        //   https://rivet.hepforge.org/analyses/STAR_2006_S6500200.html
        // ----------------------------------------------------------
        const double min_pi_pT = 3*GeV; //3.*GeV;
        const double max_pi_pT = 10.*GeV;
        const double min_p_pT  = 3.*GeV;
        const double max_p_pT  = 6.*GeV;

        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_pi_pT  && Cuts::pT < max_pi_pT  && Cuts::pid == PID::PIPLUS),  "piplusFS");
        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_pi_pT  && Cuts::pT < max_pi_pT  && Cuts::pid == PID::PIMINUS), "piminusFS");
        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_p_pT   && Cuts::pT < max_p_pT   && Cuts::pid == PID::PROTON),  "protonFS");
        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_p_pT   && Cuts::pT < max_p_pT   && Cuts::pid == PID::PBAR),    "pbarFS");

        // counters
        array<string,4> cent {"_0_10","_10_20","_20_40","_40_60"};
        array<int,   4> i_minus { 1, 3, 5, 7 };
        array<int,   4> i_plus  { 2, 4, 6, 8 };

		book ( _c["pp"],"Nev_pp" );
		book ( _c["CuCu"],"Nev_CuCu" );
		 _c_dummy["pp"] = 0.;
		 _c_dummy["CuCu"] = 0.;

        for (int i{0};i<4;++i) {
            event_ctr[i] = 0.;

            book( _c4[i], "c"+cent[i]);
            _c4_dummy[i] = 0.;

            // FIG 1a : pi+ and pi- spectra
            book( _h4["PIminus"][i], 1, 1, i_minus[i] );
            book( _h4["PIplus" ][i], 1, 1, i_plus [i] );
            // FIG 1b : p and pbar spectra
            book( _h4["PBAR"   ][i], 2, 1, i_minus[i] );
            book( _h4["P"      ][i], 2, 1, i_plus [i] );

            // FIG 2a : ratio pi- to pi+
            string refname = mkAxisCode(3,1,i+1);
            book( _s4["PIminusOverPIplus"][i], refname, true);
            
            // FIG 2b : ratio pbar to p
            refname = mkAxisCode(4,1,i+1);
            book( _s4["PBARoverP"][i], refname, true);

            // FIG 3a : Pion Raa(pT) (range 3-6 GeV/c in data, although paper has 3-8 GeV/c)
            string refnameRaa_pion = mkAxisCode(5,1,i+1);
            const Scatter2D& refdataRaa_pion =refData(refnameRaa_pion);
            book(_h4["pion_pT_CuCu"][i], refnameRaa_pion + "_CuCu" + cent[i], refdataRaa_pion);
            // Required pp values
            if (i==0) {
                book(_h["pion_pT_pp"],   refnameRaa_pion + "_pp", refdataRaa_pion);
            }
            book(_s4["Raa_pion_pT"][i],  refnameRaa_pion, true);

            // FIG 4a : Proton Raa(pT)
            string refnameRaa_proton = mkAxisCode(7,1,i+1);
            const Scatter2D& refdataRaa_proton =refData(refnameRaa_proton);
            /* book(_h4["proton_pT_CuCu"][i], refnameRaa_proton + "_CuCu" + cent[i], refdataRaa_proton); */
            if (i==0) {
                book(_h["proton_pT_pp"],   refnameRaa_proton + "_pp", refdataRaa_proton);
            }
            book(_s4["Raa_proton"][i],     refnameRaa_proton);

            // FIG 4b -- no data for this figure

            // FIG 5a : p+p over pi+pi in 3-6 GeV range -- will use spectra from and FIG 1b for p+pbar
            //          and collect, newly, pions in 3-6 GeV
            refname = mkAxisCode(8,1,i+1);
            book (_s4["ratio_PtoPI"][i], refname, true);

            vector<double> bins { 3.,3.25,3.5,3.75,4.,4.5,5.,5.5,6. };
            book (_h4["PI_3_to_6"][i], "PI_3_to_6"+cent[i], bins );
        }
        // FIG 3b : Raa per centrality bin
        string refnameRaa_PIcent = mkAxisCode(6,1,2);
        const Scatter2D& refdataRaa_PIcent =refData(refnameRaa_PIcent);
        book(_h["pion_Cent_CuCu"], refnameRaa_PIcent + "_CuCu", refdataRaa_PIcent);
        book(_h["pion_Cent_pp"],   refnameRaa_PIcent + "_pp"  , refdataRaa_PIcent);
        book(_s["Raa_pion_Cent"],  refnameRaa_PIcent);

        // FIG 5b : Raa per centrality bin
		string refname_piToP_34 = mkAxisCode(9,1,2);
		const Scatter2D& refdata_piToP_34 =refData(refname_piToP_34);
		book(_h["PI_Cent_34"], "PI_Cent_34", refdata_piToP_34);
		book(_h["P_Cent_34"],  "P_Cent_34",  refdata_piToP_34);
		book(_s["PItoP_34"],   refname_piToP_34);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        /* cout << " event: " << w_all++ << endl; */
        if (debug > 2) {
            cout << " weights: ";
            for (auto& w : event.weights()) cout << " " << w;
            cout << endl;
        }

		int pid_p  = 2212;
		int NN_Cu  = 63;
		int pid_Cu = 1000290630;

		bool isCu;
		const ParticlePair& beam = beams();
		double NN = 0;
		string beamName = "Empty";
		if (beam.first.pid() == pid_Cu && beam.second.pid() == pid_Cu)
		{
			isCu = true;
			beamName = "CuCu";
			NN = NN_Cu;
			if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) beamName += "200GeV";
			else vetoEvent;
		} else if (beam.first.pid() == pid_p && beam.second.pid() == pid_p)
		{
			isCu = false;
			beamName = "pp";
			NN = 1.;
			if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) beamName += "200GeV";
			else vetoEvent;
		} else {
			vetoEvent;
		}

        const ChargedFinalState& bbc_east = applyProjection<ChargedFinalState>(event,"BBC_East");
        const ChargedFinalState& bbc_west = applyProjection<ChargedFinalState>(event,"BBC_West");
        if (!(bbc_east.size() || bbc_west.size())){
            vetoEvent;
        }

        // get final state particles:
        int k {0}; 
        double cent{-1};
        if (isCu) {
            // FIXME -- remove test_factor. 
            // It is used for testing because centraliy definition in Au+Au
            double const test_factor {0.5};
            cent = apply<CentralityProjection>(event,"CMULT")() * test_factor;
            if      (cent < 10.)  k = 0;
            else if (cent < 20.)  k = 1;
            else if (cent < 40.)  k = 2;
            else if (cent < 60.)  k = 3;
            else vetoEvent;
        }

        // All event vetos are passed. Get the particles
        const Particles& piplus  = apply<ALICE::PrimaryParticles>(event,"piplusFS").particles();
        const Particles& piminus = apply<ALICE::PrimaryParticles>(event,"piminusFS").particles();
        const Particles& proton  = apply<ALICE::PrimaryParticles>(event,"protonFS").particles();
        const Particles& pbar    = apply<ALICE::PrimaryParticles>(event,"pbarFS").particles();

		if (isCu) {
			_c_dummy["CuCu"] += 1.;
            _c4_dummy[k] += 1.;
            /* if (debug > 0 && k==2) */ 
                /* cout << " counter: _c4["<<k<<"]->sumW() : " << _c4[k]->sumW() << "  event: " << w2++ << endl; */

			// fill in the spectra
			for (auto& p : piminus) _h4["PIminus"][k]->fill(p.pT()/GeV);
			for (auto& p : piplus ) _h4["PIplus" ][k]->fill(p.pT()/GeV);
			for (auto& p : proton ) _h4["P"      ][k]->fill(p.pT()/GeV);
			for (auto& p : pbar   ) _h4["PBAR"   ][k]->fill(p.pT()/GeV);

			// fill in values for the RAA
			for (auto& p : piminus) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				/* cout << "w: " << w << endl; */
				if (pT < 6.) {
					_h4["pion_pT_CuCu"][k]->fill(pT,w);
					_h4["PI_3_to_6"][k]->fill(pT,w);
				}
				if (pT > 5. && pT < 8.) _h["pion_Cent_CuCu"]->fill(cent,1./Npart[k]);
				if (pT < 4. ) _h["PI_Cent_34"]->fill(cent); 
			}
			for (auto& p : piplus ) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				if (pT < 6.) {
					_h4["pion_pT_CuCu"][k]->fill(pT,w);
					_h4["PI_3_to_6"][k]->fill(pT,w);
				}
				if (pT > 5. && pT < 8.) _h["pion_Cent_CuCu"]->fill(cent,1./Npart[k]);
				if (pT < 4. ) _h["PI_Cent_34"]->fill(cent);
			}
			for (auto& p : pbar ) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				_h4["PBAR"][k]->fill(pT,w);
				if (pT < 4. ) _h["P_Cent_34"]->fill(cent);
			}
			for (auto& p : proton) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				_h4["P"][k]->fill(pT,w);
				if (pT < 4. ) _h["P_Cent_34"]->fill(cent);
			}
		} else { // is pp collision
			_c_dummy["pp"] += 1.;
			// fill in values for the RAA
			for (auto& p : piminus) {
				const double pT {p.pT()/GeV};
				double w  { 1. / (2*pi*pT) };
				if (pT < 6.) _h["pion_pT_pp"]->fill(pT,w);
				if (pT > 5. && pT < 8.) _h["pion_Cent_pp"]->fill(cent);
			}
			for (auto& p : piplus ) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				if (pT < 6.) _h["pion_pT_pp"]->fill(pT,w);
				if (pT > 5. && pT < 8.) _h["pion_Cent_pp"]->fill(cent);
			}
			for (auto& p : proton) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				_h["proton_pT_pp"]->fill(pT,w);
                // data missing for Raa for protons at 5-6 GeV/c
			}
			for (auto& p : pbar ) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				_h["proton_pT_pp"]->fill(pT,w);
                // data missing for Raa for protons at 5-6 GeV/c
			}
		}
    }


    /// Normalise histograms etc., after the run
    void finalize() {
        if (debug > 1) cout << " ZEBRA 2" << endl;

        bool has_CuCu { _c_dummy["CuCu"] != 0 }; //_c["CuCu"]->sumW() != 0 };
        bool has_pp   { _c_dummy["pp"] != 0 }; //_c["pp"]  ->sumW() != 0 };

        if (debug > 0) {
            cout << " has_CuCu: " << has_CuCu << endl;
            cout << " has_pp:   " << has_pp   << endl;
            
            cout << " has: " << _c["pp"]->sumW() << " pp events. " << endl;
            array<string, 4> i_str{"0_10","10_20","20_40","40_60"};
            if (debug > 2) {
                for (int i{0};i<4;++i) {
                    cout << " has: " << _c4[i]->sumW() 
                         << " Cu+Cu "<<i_str[i]<<" events." << endl;
                }
            }
            /* if (i==2) { */
                /* cout << " w2 weighting: " << weight << endl; */
            /* } */
        }

        if (has_pp) {
            _h["pion_pT_pp"  ]->scaleW(1./_c["pp"]->sumW());
        }

        if (has_CuCu) {
            array<string, 4> i_str{"0_10","10_20","20_40","40_60"};
            for (int i{0}; i<4; ++i) {
                //FIXME
                /* if (_c4[i]->sumW() == 0) { */
                if (_c4_dummy[i] == 0.) {
                    cout << " No events for centrality class " << i_str[0] << endl;
                    continue;
                }
                // FIG 1a : pi+ and pi- spectra
                _h4["PIminus"][i]->scaleW(1./_c4_dummy[i]);		
                _h4["PIplus" ][i]->scaleW(1./_c4_dummy[i]);		

                // FIG 1b : p and pbar spectra
                _h4["P"      ][i]->scaleW(1./_c4_dummy[i]);		
                _h4["PBAR"   ][i]->scaleW(1./_c4_dummy[i]);		

                // FIG 2a : ratio pi- to pi+
                divide (_h4["PIminus"][i], _h4["PIplus"][i], _s4["PIminusOverPIplus"][i]);

                // FIG 2b : ratio pbar to p
                divide (_h4["PBAR"][i], _h4["P"][i], _s4["PBARoverP"][i]);

                // FIG 3a : Pion Raa(pT)
                if (has_pp) {
                    _h4["pion_pT_CuCu"][i]->scaleW(1./_c4_dummy[i]);
                    divide(_h4["pion_pT_CuCu"][i], _h["pion_pT_pp"], _s4["Raa_pion_pT"][i]);
                    _s4["Raa_pion_pT"][i]->scaleY(1./Npart[i]);
                }

                // FIG 3b : (outside of centrality loop)
                // ----

                // need spectra of p + pbar for FIG 4a and FIG 5a
                string name = "PandPBAR" + i_str[i];
                Histo1D PandPBAR = Histo1D(*_h4["P"][i], name);
                PandPBAR += *_h4["PBAR"][i];

                // FIG 4a : Proton Raa(pT)
                if (has_pp) {
                    *_s4["Raa_proton"][i] = YODA::divide(PandPBAR, *_h["proton_pT_pp"]);
                }

                // FIG 4b -- no data for this figure

                // FIG 5a : p+p over pi+pi in 3-6 GeV range
                _h4["PI_3_to_6"][i]->scaleW(1./_c4_dummy[i]);
                *_s4["ratio_PtoPI"][i] = YODA::divide( PandPBAR, *_h4["PI_3_to_6"][i]);

                // FIG 5b (outside of centrality loop)
            }

            if (has_pp) {
                // FIG 3b : Raa per centrality bin
                _h["pion_Cent_CuCu"]->scaleW(1./_c["CuCu"]->sumW());
                _h["pion_Cent_pp"]  ->scaleW(1./_c["pp"  ]->sumW());
                divide(_h["pion_Cent_CuCu"], _h["pion_Cent_pp"], _s["Raa_pion_Cent"]);

                // FIG 5b : Raa per centrality bin
                _h["P_Cent_34"] ->scaleW(1./_c["CuCu"]->sumW());
                _h["PI_Cent_34"] ->scaleW(1./_c["pp"  ]->sumW());
                divide(_h["P_Cent_34"], _h["PI_Cent_34"], _s["PItoP_34"]);
            }
        }
    }

    //@{
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c; // Counters: number of p+p event, Cu+Cu events
    map<string, double> _c_dummy; // Counters: number of p+p event, Cu+Cu events

	array<double, 4> Npart{ 99., 74.6, 45.9, 21.5 };

    // There are a set of four indentical measurements for four centralities
    map<string, array<Histo1DPtr,4>> _h4; // Histograms at 4 bins of Cu+Cu centralities
    array<CounterPtr,4>              _c4; // Counts of 4 bins of Cu+Cu centralities
    array<double,4>              _c4_dummy; // Counts of 4 bins of Cu+Cu centralities
    array<double,4> event_ctr;

    map<string, array<Scatter2DPtr,4>> _s4; // Scatterplots of divisions at 4 bins of Cu+Cu centralities
    map<string, Scatter2DPtr> _s;           // Ratio plots for non-ratio events

    int debug { 1 };

    int w2;
    int w_all;

  };

  DECLARE_RIVET_PLUGIN(STAR_2010_I837075);
}

