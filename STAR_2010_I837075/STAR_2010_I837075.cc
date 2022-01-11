// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "../Centralities/RHICCentrality.hh"
#include <iostream>

namespace Rivet {
  /// @brief Add a short analysis description here
  class STAR_2010_I837075 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2010_I837075);

    /// Book histograms and initialise projections before the run
    void init() {

        beamOpt = getOption<string>("beam", "NONE");

    declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

        // ----------------------------------------------------------
        // Cuts used to see if there is a trigger hit in the BBC
        // Refer to example/precedent from:
        //   https://rivet.hepforge.org/analyses/ALICE_2010_I880049
        // ----------------------------------------------------------
        /* declare(ChargedFinalState((Cuts::eta >= 3.3 && Cuts::eta < 5.0) && */
        /*             Cuts::pT > 0.1*GeV), "BBC_East"); */
        /* declare(ChargedFinalState((Cuts::eta > -5.0 && Cuts::eta < -3.3) && */
        /*             Cuts::pT > 0.1*GeV), "BBC_West"); */
        // note: this cut probably isn't very restrictive in Cu+Cu, but more
        //       often in pp

        // ----------------------------------------------------------
        // Cuts used to get final state pi+, pi-, proton, pbar
        // Refer to example/precedent from:
        //   https://rivet.hepforge.org/analyses/STAR_2006_S6500200.html
        // ----------------------------------------------------------
        const double min_pi_pT =  3. * GeV;
        const double max_pi_pT = 10. * GeV;
        const double min_p_pT  =  3. * GeV;
        const double max_p_pT  =  6. * GeV;

        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_pi_pT  && Cuts::pT < max_pi_pT  && Cuts::pid == PID::PIPLUS),  "piplusFS");
        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_pi_pT  && Cuts::pT < max_pi_pT  && Cuts::pid == PID::PIMINUS), "piminusFS");
        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_p_pT   && Cuts::pT < max_p_pT   && Cuts::pid == PID::PROTON),  "protonFS");
        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_p_pT   && Cuts::pT < max_p_pT   && Cuts::pid == PID::PBAR),    "pbarFS");

        // counters
        array<int,   4> i_minus { 1, 3, 5, 7 };
        array<int,   4> i_plus  { 2, 4, 6, 8 };

		book ( _c["pp"],  "_Nev_pp" );
		book ( _c["CuCu"],"_Nev_CuCu" );

        for (int i{0};i<4;++i) {
            event_ctr[i] = 0.;

            book( _c4[i], "_c"+i_str[i]);
            /* _c4_dummy[i] = 0.; */

            // --------------------------------------------
            // each event is fill in analyze()  :  i.e. [A]
            //                    or finalize() :  i.e. (f)
            // --------------------------------------------

            // FIG 1a : pi+ and pi- spectra
            book( _h4["PIminus"][i], 1, 1, i_minus[i]) ; // 3-10 GeV/c [A]
            book( _h4["PIplus" ][i], 1, 1, i_plus [i]) ; // [A]

            // FIG 1b : p and pbar spectra
            book( _h4["PBAR"   ][i], 2, 1, i_minus[i]) ; // 3-6 GeV/c [A]
            book( _h4["P"      ][i], 2, 1, i_plus [i]) ; // [A]
            /* book( _h4["PBAR"   ][i], 2, 1, i_minus[i] ); // 3-6 GeV/c [A] */
            /* book( _h4["P"      ][i], 2, 1, i_plus [i] ); // [A] */

            // FIG 2a : ratio pi- to pi+
            book( _s4["PIminusOverPIplus"][i], mkAxisCode(3,1,i+1), true); // (f)

            // FIG 2b : ratio pbar to p
            book( _s4["PBARoverP"][i], mkAxisCode(4,1,i+1), true); // (f)

            // FIG 3a : Pion Raa(pT) (range 3-6 GeV/c in data, although paper has 3-8 GeV/c)
            book(_s4["Raa_pion_pT"][i], mkAxisCode(5,1,i+1), true); // (f)
            book(_h4["_pion_pT_CuCu"][i], "_pion_Pt_CuCu"+i_str[i], refData(5,1,i+1)); // [A]
            if (i==0) {
                book(_h["_pion_pT_pp"], "_pion_pT_pp", refData(5,1,i+1)); // [A]
            }

            // FIG 3b : Pion RAA (cent)
            if (i==0) {
                book(_h["pion_cent_CuCu_58"], 6,1,2); // [A]
                book(_c["_pion_cent_pp_58"],   "_c_pionCent");    // [A]
            }

            // FIG 4a : Proton Raa(pT)
            book(_s4["Raa_proton_pT"][i], mkAxisCode(7,1,i+1), true); // (f)
            book( _h4["_PPBAR"]  [i], "_ppbar"+i_str[i], refData(7,1,i+1)); //  (f)
            if (i==0) {
                book(_h["_proton_pT_pp"],   "_proton_pt_pp", refData(7,1,1));// [A]
            }

            // FIG 4b -- no data for this figure

            // FIG 5a : p+p over pi+pi in 3-6 GeV range -- will use spectra from and FIG 1b for p+pbar
            //          and collect, newly, pions in 3-6 GeV
            book (_s4["ratio_PtoPI"][i], mkAxisCode(8,1,i+1), true); // (f)
            // use _h4["PPBAR"] for protons
            book (_h4["_PI_3_to_6"][i], "_PI_3_to_6"+i_str[i], refData(8,1,i+1) ); // (f)

            // FIG 5b : ppbar and pions centrality
            if (i==0) {
                book(_s["PItoP_34"],  mkAxisCode(9,1,2), true); // (f)
                book(_h["_PI_Cent_34"], "_PI_Cent_34", refData(9,1,2)); // [A]
                book(_h["_P_Cent_34"],  "_P_Cent_34",  refData(9,1,2)); // [A]
            }
        }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ParticlePair& beam = beams();

      if (beamOpt == "NONE") {

        if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000290630)
        {
          isCu = true;
        } else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) {
          isCu = false;
        }
      }

      else if (beamOpt == "CUCU200") isCu = true;
      else if (beamOpt == "PP200") isCu = false;

        // move to init -- init uses data from the first event and fill a global flag
        // assumes all events are the same in analyze

        /* const ChargedFinalState& bbc_east = applyProjection<ChargedFinalState>(event,"BBC_East"); */
        /* const ChargedFinalState& bbc_west = applyProjection<ChargedFinalState>(event,"BBC_West"); */
        /* if (!bbc_east.size() && !bbc_west.size()){ */
        /*     cout << " n_events vetoed " << temp_nVetoEvent++ << endl; */
        /*     cout << " sizes: " << bbc_east.size() << " " << bbc_west.size() << endl; */
        /* } */

        // get final state particles:
        int k {0};
        double cent{-2.};
        if (isCu) {
            cent = apply<CentralityProjection>(event,"CMULT")();
            /* cout << cent << endl; */
            if      (cent < 10.)  k = 0;
            else if (cent < 20.)  k = 1;
            else if (cent < 40.)  k = 2;
            else if (cent < 60.)  k = 3;
            else vetoEvent;
        }

        // All event vetos are passed. Get the particles
        const Particles& piplus  = apply<ALICE::PrimaryParticles>(event,"piplusFS") .particles();
        const Particles& piminus = apply<ALICE::PrimaryParticles>(event,"piminusFS").particles();
        const Particles& proton  = apply<ALICE::PrimaryParticles>(event,"protonFS") .particles();
        const Particles& pbar    = apply<ALICE::PrimaryParticles>(event,"pbarFS")   .particles();

		if (isCu) {
            // ------------------
            //  Cu+Cu collisions
            // ------------------
            _c["CuCu"]->fill();
            _c4[k]->fill();

			// fill in values for the RAA
			for (auto& p : piminus) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				_h4["PIminus"][k]->fill(pT,w);
				if (pT < 6.) {
					_h4["_pion_pT_CuCu"][k]->fill(pT,w);
					_h4["_PI_3_to_6"][k]->fill(pT,w);
				}
				if (pT > 5. && pT < 8.) _h["pion_cent_CuCu_58"]->fill(cent,1./Npart[k]);
				if (pT < 4. ) _h["_PI_Cent_34"]->fill(cent);
			}
			for (auto& p : piplus ) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
			    _h4["PIplus" ][k]->fill(pT,w);
				if (pT < 6.) {
					_h4["_pion_pT_CuCu"][k]->fill(pT,w);
					_h4["_PI_3_to_6"][k]->fill(pT,w);
				}
				if (pT > 5. && pT < 8.) _h["pion_cent_CuCu_58"]->fill(cent,1./Npart[k]);
				if (pT < 4. ) _h["_PI_Cent_34"]->fill(cent);
			}
			for (auto& p : pbar ) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				_h4["PBAR"][k]->fill(pT,w);
				if (pT < 4. ) _h["_P_Cent_34"]->fill(cent);
			}
			for (auto& p : proton) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
			    _h4["P"][k]->fill(pT,w);
				if (pT < 4. ) _h["_P_Cent_34"]->fill(cent);
			}
		} else {
            // ---------------
            //  pp collisions
            // ---------------
            _c["pp"]->fill();
			// fill in values for the RAA
			for (auto& p : piminus) {
				const double pT {p.pT()/GeV};
				double w  { 1. / (2*pi*pT) };
				if (pT < 6.) _h["_pion_pT_pp"]->fill(pT,w);
				if (pT > 5. && pT < 8.) _c["_pion_cent_pp_58"]->fill();
			}
			for (auto& p : piplus ) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				if (pT < 6.) _h["_pion_pT_pp"]->fill(pT,w);
				if (pT > 5. && pT < 8.) _c["_pion_cent_pp_58"]->fill();
			}
			for (auto& p : proton) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				_h["_proton_pT_pp"]->fill(pT,w);
                // data missing for Raa for protons at 5-6 GeV/c
			}
			for (auto& p : pbar ) {
				const double pT {p.pT()/GeV};
				const double w  { 1. / (2*pi*pT) };
				_h["_proton_pT_pp"]->fill(pT,w);
                // data missing for Raa for protons at 5-6 GeV/c
			}
		}
    }


    /* /// Normalise histograms etc., after the run */
    void finalize() {


        /* cout << " n_events vetoed " << temp_nVetoEvent << endl; */

        // don't process any histograms unless there are both Cu+Cu and pp events
        for (int i{0}; i<4; ++i) {
            // FIG 1a : pi+ and pi- spectra
            _h4["PIminus"][i]->scaleW(1./_c4[i]->sumW());
            _h4["PIplus" ][i]->scaleW(1./_c4[i]->sumW());

            // FIG 1b : p and pbar spectra
            _h4["P"      ][i]->scaleW(1./_c4[i]->sumW());
            _h4["PBAR"   ][i]->scaleW(1./_c4[i]->sumW());

            // FIG 2a : ratio pi- to pi+
            divide (_h4["PIminus"][i], _h4["PIplus"][i], _s4["PIminusOverPIplus"][i]);

            // FIG 2b : ratio pbar to p
            divide (_h4["PBAR"][i], _h4["P"][i], _s4["PBARoverP"][i]);

            // FIG 3a : Pion Raa(pT)
            _h["_pion_pT_pp"]->scaleW(1./_c["pp"]->sumW());
            _h4["_pion_pT_CuCu"][i]->scaleW(1./_c4[i]->sumW()); //_c4_dummy[i]);
            divide(_h4["_pion_pT_CuCu"][i], _h["_pion_pT_pp"], _s4["Raa_pion_pT"][i]);
            _s4["Raa_pion_pT"][i]->scaleY(1./Npart[i]);

            // FIG 3b : Pion RAA in centrality bins
            _h["pion_cent_CuCu_58"]->scaleW(1./_c["CuCu"]->sumW());
            _c["_pion_cent_pp_58" ]->scaleW(1./_c["pp"  ]->sumW());
            _h["pion_cent_CuCu_58"]->scaleW(1./_c["_pion_cent_pp_58"]->sumW());


            // need spectra of p + pbar for FIG 4a and FIG 5a
            *_h4["_PPBAR"][i]  = *_h4["P"][i];
            *_h4["_PPBAR"][i] += *_h4["PBAR"][i];

            // FIG 4a : Proton Raa(pT)
            _h["_proton_pT_pp"]->scaleW(1./_c["pp"]->sumW());
            divide(_h4["_PPBAR"][i],_h["_proton_pT_pp"],_s4["Raa_proton_pT"][i]);

            // FIG 4b -- no data for this figure

            // FIG 5a : p+p over pi+pi in 3-6 GeV range
            _h4["_PI_3_to_6"][i]->scaleW(1./_c4[i]->sumW()); //_c4_dummy[i]);

            divide(_h4["_PPBAR"][i],_h4["_PI_3_to_6"][i],_s4["ratio_PtoPI"][i]);

            // FIG 5b (outside of centrality loop)
            _h["_P_Cent_34"]  ->scaleW(1./_c["CuCu"]->sumW());
            _h["_PI_Cent_34"] ->scaleW(1./_c["CuCu"]->sumW());
            divide(_h["_P_Cent_34"], _h["_PI_Cent_34"], _s["PItoP_34"]);
        }
    }

    //@{
    bool isCu = false;
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c; // Counters: number of p+p event, Cu+Cu events

	array<double,4> Npart{ 99., 74.6, 45.9, 21.5 };
    array<string,4>  i_str {"_0_10","_10_20","_20_40","_40_60"};

    // There are a set of four indentical measurements for four centralities
    map<string, array<Histo1DPtr,4>> _h4;   // Histograms at 4 bins of Cu+Cu centralities
    array<CounterPtr,4>              _c4;   // Counts of 4 bins of Cu+Cu centralities
    array<double,4>              _c4_dummy; // Counts of 4 bins of Cu+Cu centralities
    array<double,4> event_ctr;

    map<string, array<Scatter2DPtr,4>> _s4; // Scatterplots of divisions at 4 bins of Cu+Cu centralities
    map<string, Scatter2DPtr> _s;           // Ratio plots for ratio events

    string beamOpt = "NONE";

    int debug { 1 };

    /* double temp_nVetoEvent {0.}; */
  };

  DECLARE_RIVET_PLUGIN(STAR_2010_I837075);
}
