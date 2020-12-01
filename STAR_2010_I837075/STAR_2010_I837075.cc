// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
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

        const FinalState fs(Cuts::abseta <= 1.0);

        book(_h["PIminus_0_10cent"],1,1,1);
        book(_h["PIplus_0_10cent"], 1,1,2);
        book(_h["PIminus_10_20cent"],1,1,3);
        book(_h["PIplus_10_20cent"], 1,1,4);
        book(_h["PIminus_20_40cent"],1,1,5);
        book(_h["PIplus_20_40cent"], 1,1,6);
        book(_h["PIminus_40_60cent"],1,1,7);
        book(_h["PIplus_40_60cent"], 1,1,8);

        book(_h["PBAR_0_10cent"], 2,1,1);
        book(_h["P_0_10cent"],    2,1,2);
        book(_h["PBAR_10_20cent"],2,1,3);
        book(_h["P_10_20cent"],   2,1,4);
        book(_h["PBAR_20_40cent"],2,1,5);
        book(_h["P_20_40cent"],   2,1,6);
        book(_h["PBAR_40_60cent"],2,1,7);
        book(_h["P_40_60cent"],   2,1,8);

        string refname0 = mkAxisCode(3, 1, 1);
        const Scatter2D& refdata0 = refData(refname0);
        book(_h["PIminus_0_10"], refname0 + "_PIminus_0_10", refdata0);
        book(_h["PIplus_0_10"],  refname0 + "_PIplus_0_10",  refdata0);
        book(_s["PIminusOverPIplus_0_10"], refname0, true);

        string refname1 = mkAxisCode(3, 1, 2);
        const Scatter2D& refdata1 = refData(refname1);
        book(_h["PIminus_10_20"], refname1 + "_PIminus_10_20", refdata1);
        book(_h["PIplus_10_20"],  refname1 + "_PIplus_10_20",  refdata1);
        book(_s["PIminusOverPIplus_10_20"], refname1, true);

        string refname2 = mkAxisCode(3, 1, 3);
        const Scatter2D& refdata2 = refData(refname2);
        book(_h["PIminus_20_40"], refname2 + "_PIminus_20_40", refdata2);
        book(_h["PIplus_20_40"],  refname2 + "_PIplus_20_40",  refdata2);
        book(_s["PIminusOverPIplus_20_40"], refname2, true);

        string refname4 = mkAxisCode(3, 1, 4);
        const Scatter2D& refdata4 = refData(refname4);
        book(_h["PIminus_40_60"], refname4 + "_PIminus_40_60", refdata4);
        book(_h["PIplus_40_60"],  refname4 + "_PIplus_40_60",  refdata4);
        book(_s["PIminusOverPIplus_40_60"], refname4, true);

        string refname5 = mkAxisCode(4, 1, 1);
        const Scatter2D& refdata5 = refData(refname5);
        book(_h["PBAR_0_10"], refname5 + "_PBAR_0_10", refdata5);
        book(_h["P_0_10"],  refname5 + "_P_0_10",  refdata5);
        book(_s["PBAROverP_0_10"], refname5, true);

        string refname6 = mkAxisCode(4, 1, 2);
        const Scatter2D& refdata6 = refData(refname6);
        book(_h["PBAR_10_20"], refname6 + "_PBAR_10_20", refdata6);
        book(_h["P_10_20"],  refname6 + "_P_10_20",  refdata6);
        book(_s["PBAROverP_10_20"], refname6, true);

        string refname7 = mkAxisCode(4, 1, 3);
        const Scatter2D& refdata7 = refData(refname7);
        book(_h["PBAR_20_40"], refname7 + "_PBAR_20_40", refdata7);
        book(_h["P_20_40"],  refname7 + "_P_20_40",  refdata7);
        book(_s["PBAROverP_20_40"], refname7, true);

        string refname8 = mkAxisCode(4, 1, 4);
        const Scatter2D& refdata8 = refData(refname8);
        book(_h["PBAR_40_60"], refname8 + "_PBAR_40_60", refdata8);
        book(_h["P_40_60"],  refname8 + "_P_40_60",  refdata8);
        book(_s["PBAROverP_40_60"], refname8, true);

        book(_h["RAA_PIONS_0_10cent"],    5,1,1);
        book(_h["RAA_PIONS_10_20cent"],   5,1,2);
        book(_h["RAA_PIONS_20_40cent"],   5,1,3);
        book(_h["RAA_PIONS_40_60cent"],   5,1,4);

        /* book(_h["RAA_PIONS_pT_5to8"],    6,1,2); */

        book(_h["RAA_PROTONS_0_10cent"],    7,1,1);
        book(_h["RAA_PROTONS_10_20cent"],   7,1,2);
        book(_h["RAA_PROTONS_20_40cent"],   7,1,3);
        book(_h["RAA_PROTONS_40_60cent"],   7,1,4);

        book(_h["RATIO_PIONStoPROTONS_0_10cent"],    8,1,1);
        book(_h["RATIO_PIONStoPROTONS_10_20cent"],   8,1,2);
        book(_h["RATIO_PIONStoPROTONS_20_40cent"],   8,1,3);
        book(_h["RATIO_PIONStoPROTONS_40_60cent"],   8,1,4);

        /* book(_h["RATIO_PIONStoPROTONS_pT_3to4"],   9,1,1); */




    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

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

      /* normalize(_h["XXXX"]); // normalize to unity */
      /* normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts) */
      /* scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts) */

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr>   _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr>   _c;
    map<string, Scatter2DPtr> _s;
    //@}


  };


  DECLARE_RIVET_PLUGIN(STAR_2010_I837075);

}
