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


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
        declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");
        /* const FinalState fs(Cuts::abseta <= 0.5); // pg. 4 analsis */

        // ----------------------------------------------------------
        // Cuts used to see if there is a trigger hit in the BBC
        // Refer to example/precedent from: 
        //   https://rivet.hepforge.org/analyses/ALICE_2010_I880049
        // ----------------------------------------------------------
        /* declare(ALICE::PrimaryParticles( */
        /*     Cuts::abseta < 0.5 */ 
        /*  && Cuts::pT >= 200*MeV */
        /*  && Cuts::abscharge > 0), "TPC_PRIM"); */

        declare(ChargedFinalState((Cuts::eta >= 3.3 && Cuts::eta < 5.0) &&
                    Cuts::pT > 0.1*GeV), "BBC_East");
        declare(ChargedFinalState((Cuts::eta > -5.0 && Cuts::eta < -3.3) &&
                    Cuts::pT > 0.1*GeV), "BBC_West");

        // ----------------------------------------------------------
        // Cuts used to get final state pi+, pi-, proton, pbar
        // Refer to example/precedent from: 
        //   https://rivet.hepforge.org/analyses/STAR_2006_S6500200.html
        // ----------------------------------------------------------
        
        // Get the only particles that I care about: p/m_pi and p/pbar
        // FIXME: real pT requirements are 3-10 GeV
        /* const double min_pt { 0.2*GeV }; */
        /* const double max_pt { 10.*GeV }; */

        /* const double minpT_pi_fig1a_2a { 3. *GeV }; */
        /* const double maxpT_pi_fig1a_2a { 10.*GeV }; */

        /* const double minpT_p_fig1b_2b { 3. *GeV }; */
        /* const double maxpT_p_fig1b_2b { 6. *GeV }; */

        /* const double minpT_ppi_fig3a_4a_5a  { 3. *GeV }; */
        /* const double maxpT_ppi_fig3a_4a_5a  { 6. *GeV }; */

        /* const double minpT_pi_fig3b  { 5. *GeV }; */
        /* const double maxpT_pi_fig3b  { 8. *GeV }; */

        /* const double minpT_p_fig4b  { 5. *GeV }; */
        /* const double maxpT_p_fig4b  { 6. *GeV }; */

        /* const double minpT_ppi_fig5b  { 3. *GeV }; */
        /* const double maxpT_ppi_fig5b  { 4. *GeV }; */
         
        const double min_pi_pT = 3.;
        const double max_pi_pT = 10.;
        const double min_p_pT  = 3.;
        const double max_p_pT  = 10.;

        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_pi_pT  && Cuts::pT < max_pi_pT  && Cuts::pid == PID::PIPLUS),  "piplusFS");
        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_pi_pT  && Cuts::pT < max_pi_pT  && Cuts::pid == PID::PIMINUS), "piminusFS");
        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_p_pT   && Cuts::pT < max_p_pT   && Cuts::pid == PID::PROTON),  "protonFS");
        declare( ALICE::PrimaryParticles( Cuts::abseta < 0.5 && Cuts::pT > min_p_pT   && Cuts::pT < max_p_pT   && Cuts::pid == PID::PBAR),    "pbarFS");

        // counters
        array<string,4> cent {"_0_10","_10_20","_20_40","_40_60"};
        array<int,   4> i_minus { 1, 3, 5, 7 };
        array<int,   4> i_plus  { 2, 4, 6, 8 };
        for (int i{0};i<4;++i) {
            book( _c4[i], "c"+cent[i]);

            book( _h4["PIminus"][i], 1, 1, i_minus[i] );
            book( _h4["PIplus" ][i], 1, 1, i_plus [i] );
            book( _h4["PBAR"   ][i], 2, 1, i_minus[i] );
            book( _h4["P"      ][i], 2, 1, i_plus [i] );

            string refname = mkAxisCode(3,1,i+1);
            book( _s4["PIminusOverPIplus"][i], refname, true);

            refname = mkAxisCode(4,1,i+1);
            book( _s4["PBARoverP"][i], refname, true);
            book (_h4["RAA_PI"]     [i], 5, 1, i+1);
            /* book(_h["RAA_PIONS_pT_5to8"],    6,1,2); */
            book (_h4["RAA_P"]      [i], 7, 1, i+1);
            book (_h4["ratio_PItoP"][i], 8,1,i+1);
            /* book(_h["RATIO_PIONStoPROTONS_pT_3to4"],   9,1,1); */
        }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

        cout << " ZEBRA 0" << endl;
        const ChargedFinalState& bbc_east = applyProjection<ChargedFinalState>(event,"BBC_East");
        const ChargedFinalState& bbc_west = applyProjection<ChargedFinalState>(event,"BBC_West");
        if (!(bbc_east.size() || bbc_west.size())){
            vetoEvent;
        }

        // get final state particles:
        // FIXME -- remove test_factor
        double const test_factor {0.5};
        const double cent = apply<CentralityProjection>(event,"CMULT")() * test_factor;
        cout << "Centrality: " << cent << endl;
        int k = 0; // k is centrality index
        if      (cent < 10.)  k = 0;
        else if (cent < 20.)  k = 1;
        else if (cent < 40.)  k = 2;
        else if (cent < 60.)  k = 3;
        else vetoEvent;
        /* _c4[c]->fill(); */

        const Particles& piplus  = apply<ALICE::PrimaryParticles>(event,"piplusFS").particles();
        const Particles& piminus = apply<ALICE::PrimaryParticles>(event,"piminusFS").particles();
        const Particles& proton  = apply<ALICE::PrimaryParticles>(event,"protonFS").particles();
        const Particles& pbar    = apply<ALICE::PrimaryParticles>(event,"pbarFS").particles();

        for (auto& p : piminus) _h4["PIminus"][k]->fill(p.pT()/GeV);
        for (auto& p : piplus ) _h4["PIminus"][k]->fill(p.pT()/GeV);
        for (auto& p : proton ) _h4["P"      ][k]->fill(p.pT()/GeV);
        for (auto& p : pbar   ) _h4["PBAR"   ][k]->fill(p.pT()/GeV);

        /* // ok, fill all the histograms with appropriate spectra */
        /* if (cent < 10.) { */
        /*     for (auto& p : p_pi)  _h["PIplus_0_10cent"]->fill(p.pT()); // no need pT cut -- this is the max range; */
        /*     for (auto& p : n_pi)  _h["PIminus_0_10cent"]->fill(p.pT()); // no need pT cut -- this is the max range; */
        /*     for (auto& p : proton) { */
        /*         const double pT {p.pT()}; */ 
        /*         if (pT < 6.)  _h["P_0_10cent"]->fill(p.pT()); */
        /*     } */
        /*     for (auto& p : pbar) { */
        /*         const double pT {p.pT()}; */ 
        /*         if (pT < 6.)  _h["PBAR_0_10cent"]->fill(p.pT()); */
        /*     } */
        /* } else if (cent < 20.) { */ 
        /*     for (auto& p : p_pi)  _h["PIplus_10_20cent"]->fill(p.pT()); // no need pT cut -- this is the max range; */
        /*     for (auto& p : n_pi)  _h["PIminus_10_20cent"]->fill(p.pT()); // no need pT cut -- this is the max range; */
        /*     for (auto& p : proton) { */
        /*         const double pT {p.pT()}; */ 
        /*         if (pT < 6.)  _h["P_10_20cent"]->fill(p.pT()); */
        /*     } */
        /*     for (auto& p : pbar) { */
        /*         const double pT {p.pT()}; */ 
        /*         if (pT < 6.)  _h["PBAR_10_20cent"]->fill(p.pT()); */
        /*     } */

        /* } else if (cent < 40.) { */ 
        /* } else if (cent < 60.) { */ 
        /* } else vetoEvent; */

        





      /* // Retrieve dressed leptons, sorted by pT */
      /* vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons(); */

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      /* Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV); */

      /* // Remove all jets within dR < 0.2 of a dressed lepton */
      /* idiscardIfAnyDeltaRLess(jets, leptons, 0.2); */

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
      /* Jets bjets = filter_select(jets, [](const Jet& jet) { */
        /* return  jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5); */
      /* }); */

      // Veto event if there are no b-jets
      /* if (bjets.empty())  vetoEvent; */

      // Apply a missing-momentum cut
      /* if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent; */

      // Fill histogram with leading b-jet pT
      /* _h["XXXX"]->fill(bjets[0].pT()/GeV); */

    }


    /// Normalise histograms etc., after the run
    void finalize() {
        cout << " ZEBRA 1" << endl;

      /* normalize(_h["XXXX"]); // normalize to unity */
      /* normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts) */
      /* scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts) */

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr>   _c;

    // There are a set of four indentical measurements for four centralities
    map<string, array<Histo1DPtr,4>> _h4; // like _h but for four centralities
    array<CounterPtr,4>              _c4;
    /* map<string, Profile1DPtr> _p; */
    map<string, array<Scatter2DPtr,4>> _s4;
    map<string, Scatter2DPtr> _s;
    //@}


  };


  DECLARE_RIVET_PLUGIN(STAR_2010_I837075);

}
        // book all the histograms
        /* book(_h["PIminus_0_10cent"],1,1,1); */
        /* book(_h["PIplus_0_10cent"], 1,1,2); */
        /* book(_h["PIminus_10_20cent"],1,1,3); */
        /* book(_h["PIplus_10_20cent"], 1,1,4); */
        /* book(_h["PIminus_20_40cent"],1,1,5); */
        /* book(_h["PIplus_20_40cent"], 1,1,6); */
        /* book(_h["PIminus_40_60cent"],1,1,7); */
        /* book(_h["PIplus_40_60cent"], 1,1,8); */

        /* book(_h["PBAR_0_10cent"], 2,1,1); */
        /* book(_h["P_0_10cent"],    2,1,2); */
        /* book(_h["PBAR_10_20cent"],2,1,3); */
        /* book(_h["P_10_20cent"],   2,1,4); */
        /* book(_h["PBAR_20_40cent"],2,1,5); */
        /* book(_h["P_20_40cent"],   2,1,6); */
        /* book(_h["PBAR_40_60cent"],2,1,7); */
        /* book(_h["P_40_60cent"],   2,1,8); */

        /* string refname0 = mkAxisCode(3, 1, 1); */
        /* const Scatter2D& refdata0 = refData(refname0); */
        /* book(_h["PIminus_0_10"], refname0 + "_PIminus_0_10", refdata0); */
        /* book(_h["PIplus_0_10"],  refname0 + "_PIplus_0_10",  refdata0); */
        /* book(_s["PIminusOverPIplus_0_10"], refname0, true); */

        /* string refname1 = mkAxisCode(3, 1, 2); */
        /* const Scatter2D& refdata1 = refData(refname1); */
        /* book(_h["PIminus_10_20"], refname1 + "_PIminus_10_20", refdata1); */
        /* book(_h["PIplus_10_20"],  refname1 + "_PIplus_10_20",  refdata1); */
        /* book(_s["PIminusOverPIplus_10_20"], refname1, true); */

        /* string refname2 = mkAxisCode(3, 1, 3); */
        /* const Scatter2D& refdata2 = refData(refname2); */
        /* book(_h["PIminus_20_40"], refname2 + "_PIminus_20_40", refdata2); */
        /* book(_h["PIplus_20_40"],  refname2 + "_PIplus_20_40",  refdata2); */
        /* book(_s["PIminusOverPIplus_20_40"], refname2, true); */

        /* string refname4 = mkAxisCode(3, 1, 4); */
        /* const Scatter2D& refdata4 = refData(refname4); */
        /* book(_h["PIminus_40_60"], refname4 + "_PIminus_40_60", refdata4); */
        /* book(_h["PIplus_40_60"],  refname4 + "_PIplus_40_60",  refdata4); */
        /* book(_s["PIminusOverPIplus_40_60"], refname4, true); */

        /* string refname5 = mkAxisCode(4, 1, 1); */
        /* const Scatter2D& refdata5 = refData(refname5); */
        /* book(_h["PBAR_0_10"], refname5 + "_PBAR_0_10", refdata5); */
        /* book(_h["P_0_10"],  refname5 + "_P_0_10",  refdata5); */
        /* book(_s["PBAROverP_0_10"], refname5, true); */

        /* string refname6 = mkAxisCode(4, 1, 2); */
        /* const Scatter2D& refdata6 = refData(refname6); */
        /* book(_h["PBAR_10_20"], refname6 + "_PBAR_10_20", refdata6); */
        /* book(_h["P_10_20"],  refname6 + "_P_10_20",  refdata6); */
        /* book(_s["PBAROverP_10_20"], refname6, true); */

        /* string refname7 = mkAxisCode(4, 1, 3); */
        /* const Scatter2D& refdata7 = refData(refname7); */
        /* book(_h["PBAR_20_40"], refname7 + "_PBAR_20_40", refdata7); */
        /* book(_h["P_20_40"],  refname7 + "_P_20_40",  refdata7); */
        /* book(_s["PBAROverP_20_40"], refname7, true); */

        /* string refname8 = mkAxisCode(4, 1, 4); */
        /* const Scatter2D& refdata8 = refData(refname8); */
        /* book(_h["PBAR_40_60"], refname8 + "_PBAR_40_60", refdata8); */
        /* book(_h["P_40_60"],  refname8 + "_P_40_60",  refdata8); */
        /* book(_s["PBAROverP_40_60"], refname8, true); */

        /* book(_h["RAA_PIONS_0_10cent"],    5,1,1); */
        /* book(_h["RAA_PIONS_10_20cent"],   5,1,2); */
        /* book(_h["RAA_PIONS_20_40cent"],   5,1,3); */
        /* book(_h["RAA_PIONS_40_60cent"],   5,1,4); */

        /* book(_h["RAA_PROTONS_0_10cent"],    7,1,1); */
        /* book(_h["RAA_PROTONS_10_20cent"],   7,1,2); */
        /* book(_h["RAA_PROTONS_20_40cent"],   7,1,3); */
        /* book(_h["RAA_PROTONS_40_60cent"],   7,1,4); */

        /* book(_h["RATIO_PIONStoPROTONS_0_10cent"],    8,1,1); */
        /* book(_h["RATIO_PIONStoPROTONS_10_20cent"],   8,1,2); */
        /* book(_h["RATIO_PIONStoPROTONS_20_40cent"],   8,1,3); */
        /* book(_h["RATIO_PIONStoPROTONS_40_60cent"],   8,1,4); */

