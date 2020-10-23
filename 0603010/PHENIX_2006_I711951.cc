// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "RHICCentrality.hh"
#include <math.h>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES


namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2006_I711951 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2006_I711951);



    /// Book histograms and initialise projections before the run
    void init() {
      //Particles: pi^+, pi^-, k^+, k^-, p, p_bar
      //std::initializer_list<int> pdgIds = {221};

      const FinalState fs(Cuts::absrap<0.35&&Cuts::abscharge>0);
      declare(fs,"fs");

      beamOpt = getOption<string>("beam","NONE");
      if (beamOpt =="pp") collSys = pp;
      else if (beamOpt == "dAU200") collSys = dAu200;

      if (collSys != pp) {
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
      }


      //****Counters****
      //book(sow["sow_pp"],"sow_pp");
      //book(sow["sow_dAU20"],"sow_dAU20");
      //book(sow["sow_dAU40"],"sow_dAU40");
      //book(sow["sow_dAU60"],"sow_dAU60");
      //book(sow["sow_dAU88"],"sow_dAU88");
      book(sow,"sow");

      //****Book Histos****
      //book(hPionPlusPt,4,1,1);
      //book(sow,"sow");


      //for (i=0, i<4, i++){

      //Invariant Yield in pp (piplus)
      book(hPP_PiPlusPt["PiplusPP"], 1,1,1);

      //Invariant Yield in pp (piminus)
      book(hPP_PiMinusPt["PiminusPP"], 1,1,2);

      //Invariant Yield in pp (kplus)
      book(hPP_KPlusPt["KplusPP"], 5,1,1);

      //Invariant Yield in pp (kminus)
      book(hPP_KMinusPt["KminusPP"], 5,1,2);

      //Invariant Yield in pp (p)
      book(hPP_PPt["PPP"], 9,1,1);

      //Invariant Yield in pp (p_bar)
      book(hPP_P_BarPt["P_barPP"], 9,1,2);

      //Invariant Yield in dAU 0-20% (piplus)
      book(hdAu_PiPlusPt["PiplusC20"], 3,1,1);

      //Invariant Yield in dAU 0-20% (piminus)
      book(hdAu_PiMinusPt["PiminusC20"], 3,1,2);

      //Invariant Yield in dAU 0-20% (kplus)
      book(hdAu_KPlusPt["KplusC20"], 7,1,1);

      //Invariant Yield in dAU 0-20% (kminus)
      book(hdAu_KMinusPt["KminusC20"], 7,1,2);

      //Invariant Yield in dAU 0-20% (p)
      book(hdAu_PPt["PC20"], 10,1,1);

      //Invariant Yield in dAU 0-20% (p_bar)
      book(hdAu_P_BarPt["P_barC20"], 10,1,2);

      //Invariant Yield in dAU 20-40% (piplus)
      book(hdAu_PiPlusPt["PiplusC40"], 3,1,3);

      //Invariant Yield in dAU 20-40% (piminus)
      book(hdAu_PiMinusPt["PiminusC40"], 3,1,4);

      //Invariant Yield in dAU 20-40% (kplus)
      book(hdAu_KPlusPt["KplusC40"], 8,1,1);

      //Invariant Yield in dAU 20-40% (kminus)
      book(hdAu_KMinusPt["KminusC40"], 8,1,2);

      //Invariant Yield in dAU 20-40% (p)
      book(hdAu_PPt["PC40"], 10,1,3);

      //Invariant Yield in dAU 20-40% (p_bar)
      book(hdAu_P_BarPt["P_barC40"], 10,1,4);

      //Invariant Yield in dAU 40-60% (piplus)
      book(hdAu_PiPlusPt["PiplusC60"], 4,1,1);

      //Invariant Yield in dAU 40-60% (piminus)
      book(hdAu_PiMinusPt["PiminusC60"], 3,1,5);

      //Invariant Yield in dAU 40-60% (kplus)
      book(hdAu_KPlusPt["KplusC60"], 7,1,3);

      //Invariant Yield in dAU 40-60% (kminus)
      book(hdAu_KMinusPt["KminusC60"], 7,1,4);

      //Invariant Yield in dAU 40-60% (p)
      book(hdAu_PPt["PC60"], 10,1,5);

      //Invariant Yield in dAU 40-60% (p_bar)
      book(hdAu_P_BarPt["P_barC60"], 10,1,6);

      //Invariant Yield in dAU 60-88.5% (piplus)
      book(hdAu_PiPlusPt["PiplusC88"], 3,1,6);

      //Invariant Yield in dAU 60-88.5% (piminus)
      book(hdAu_PiMinusPt["PiminusC88"], 3,1,7);

      //Invariant Yield in dAU 60-88.5% (kplus)
      book(hdAu_KPlusPt["KplusC88"], 7,1,5);

      //Invariant Yield in dAU 60-88.5% (kminus)
      book(hdAu_KMinusPt["KminusC88"], 7,1,6);

      //Invariant Yield in dAU 60-88.5% (p)
      book(hdAu_PPt["PC88"], 10,1,7);

      //Invariant Yield in dAU 60-88.5% (p_bar)
      book(hdAu_P_BarPt["P_barC88"], 10,1,8);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      sow->fill();
      Particles chargedP = applyProjection<FinalState>(event,"fs").particles();

      if (collSys == pp) {
        //sow["sow_pp"]->fill();
        for (const Particle& p :chargedP) {
          double partPt = p.pT() / GeV;
          double pt_weight = 1. / (partPt*2.*M_PI);

          if (p.pid() == 211) {
            hPP_PiPlusPt["PiplusPP"]->fill(partPt, pt_weight);
          }
          if (p.pid() == -211) {
            hPP_PiMinusPt["PiminusPP"]->fill(partPt, pt_weight);
          }
          if (p.pid() == 321) {
            hPP_KPlusPt["KplusPP"]->fill(partPt, pt_weight);
          }
          if (p.pid() == -321) {
            hPP_KMinusPt["KminusPP"]->fill(partPt, pt_weight);
          }
          if (p.pid() == 2212) {
            if (!(p.hasAncestor(3122) ||                          //Lambda
                  p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                  p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                  p.hasAncestor(3334) )) {                        //Omega-
                    hPP_PPt["PPP"]->fill(partPt, pt_weight);
                  }
          }
          if (p.pid() == -2212) {
            if (!(p.hasAncestor(3122) ||                          //Lambda
                  p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                  p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                  p.hasAncestor(3334) )) {                        //Omega-
                    hPP_P_BarPt["P_barPP"]->fill(partPt, pt_weight);
                  }
          }
        }


      }

      if (collSys == dAu200) {

        const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
        const double cent = centProj();

        if (cent < 0. || cent > 88.5) vetoEvent;

        if (cent > 0. && cent < 20.) {
          //sow["sow_dAU20"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);

            if (p.pid() == 211) {
              hdAu_PiPlusPt["PiplusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == -211) {
              hdAu_PiMinusPt["PiminusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == 321) {
              hdAu_KPlusPt["KplusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == -321) {
              hdAu_KMinusPt["KminusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) ||                          //Lambda
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      hdAu_PPt["PC20"]->fill(partPt, pt_weight);
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) ||                          //Lambda
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      hdAu_P_BarPt["P_barC20"]->fill(partPt, pt_weight);
                    }
            }
          }

        }

        if (cent > 20. && cent < 40.) {
          //sow["sow_dAU40"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);

            if (p.pid() == 211) {
              hdAu_PiPlusPt["PiplusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == -211) {
              hdAu_PiMinusPt["PiminusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == 321) {
              hdAu_KPlusPt["KplusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == -321) {
              hdAu_KMinusPt["KminusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == 2212) {
              hdAu_PPt["PC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == -2212) {
              hdAu_P_BarPt["P_barC20"]->fill(partPt, pt_weight);
            }
          }


        }

        if (cent > 40. && cent < 60.) {
          //sow["sow_dAU60"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);

            if (p.pid() == 211) {
              hdAu_PiPlusPt["PiplusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == -211) {
              hdAu_PiMinusPt["PiminusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == 321) {
              hdAu_KPlusPt["KplusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == -321) {
              hdAu_KMinusPt["KminusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == 2212) {
              hdAu_PPt["PC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == -2212) {
              hdAu_P_BarPt["P_barC20"]->fill(partPt, pt_weight);
            }
          }

        }

        if (cent > 60. && cent < 88.5) {
          //sow["sow_dAU88"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);

            if (p.pid() == 211) {
              hdAu_PiPlusPt["PiplusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == -211) {
              hdAu_PiMinusPt["PiminusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == 321) {
              hdAu_KPlusPt["KplusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == -321) {
              hdAu_KMinusPt["KminusC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == 2212) {
              hdAu_PPt["PC20"]->fill(partPt, pt_weight);
            }
            if (p.pid() == -2212) {
              hdAu_P_BarPt["P_barC20"]->fill(partPt, pt_weight);
            }
          }

        }




        /**
        for(const Particle& p : chargedP)
        {
        if(p.pid() == 211)
        {
            double partPt = p.pT()/GeV;
                        double pt_weight = 1./(partPt*2.*M_PI);
            hPionPlusPt->fill(partPt, pt_weight);
        }

        }
        **/
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //****Scale Histos****
/**
      hdAu_PiPlusPt["PiplusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_PiMinusPt["PiminusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_KPlusPt["KplusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_KMinusPt["KminusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_PPt["PC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_P_BarPt["P_barC20"]->scaleW(1./sow["sow_dAU20"]->sumW());

      hdAu_PiPlusPt["PiplusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_PiMinusPt["PiminusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_KPlusPt["KplusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_KMinusPt["KminusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_PPt["PC40"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_P_BarPt["P_barC40"]->scaleW(1./sow["sow_dAU40"]->sumW());

      hdAu_PiPlusPt["PiplusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_PiMinusPt["PiminusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_KPlusPt["KplusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_KMinusPt["KminusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_PPt["PC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_P_BarPt["P_barC60"]->scaleW(1./sow["sow_dAU60"]->sumW());

      hdAu_PiPlusPt["PiplusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_PiMinusPt["PiminusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_KPlusPt["KplusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_KMinusPt["KminusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_PPt["PC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_P_BarPt["P_barC88"]->scaleW(1./sow["sow_dAU88"]->sumW());

      hPP_PiPlusPt["PiplusPP"]->scaleW(1./sow["sow_PP"]->sumW());
      hPP_PiMinusPt["PiminusPP"]->scaleW(1./sow["sow_PP"]->sumW());
      hPP_KPlusPt["KplusPP"]->scaleW(1./sow["sow_PP"]->sumW());
      hPP_KMinusPt["KminusPP"]->scaleW(1./sow["sow_PP"]->sumW());
      hPP_PPt["PPP"]->scaleW(1./sow["sow_PP"]->sumW());
      hPP_P_BarPt["P_barPP"]->scaleW(1./sow["sow_PP"]->sumW());

**/

      const double s = 1./sow->sumW();
      hdAu_PiPlusPt["PiplusC20"]->scaleW(s);
      hdAu_PiMinusPt["PiminusC20"]->scaleW(s);
      hdAu_KPlusPt["KplusC20"]->scaleW(s);
      hdAu_KMinusPt["KminusC20"]->scaleW(s);
      hdAu_PPt["PC20"]->scaleW(s);
      hdAu_P_BarPt["P_barC20"]->scaleW(s);

      hdAu_PiPlusPt["PiplusC40"]->scaleW(s);
      hdAu_PiMinusPt["PiminusC40"]->scaleW(s);
      hdAu_KPlusPt["KplusC40"]->scaleW(s);
      hdAu_KMinusPt["KminusC40"]->scaleW(s);
      hdAu_PPt["PC40"]->scaleW(s);
      hdAu_P_BarPt["P_barC40"]->scaleW(s);

      hdAu_PiPlusPt["PiplusC60"]->scaleW(s);
      hdAu_PiMinusPt["PiminusC60"]->scaleW(s);
      hdAu_KPlusPt["KplusC60"]->scaleW(s);
      hdAu_KMinusPt["KminusC60"]->scaleW(s);
      hdAu_PPt["PC60"]->scaleW(s);
      hdAu_P_BarPt["P_barC60"]->scaleW(s);

      hdAu_PiPlusPt["PiplusC88"]->scaleW(s);
      hdAu_PiMinusPt["PiminusC88"]->scaleW(s);
      hdAu_KPlusPt["KplusC88"]->scaleW(s);
      hdAu_KMinusPt["KminusC88"]->scaleW(s);
      hdAu_PPt["PC88"]->scaleW(s);
      hdAu_P_BarPt["P_barC88"]->scaleW(s);

      hPP_PiPlusPt["PiplusPP"]->scaleW(s);
      hPP_PiMinusPt["PiminusPP"]->scaleW(s);
      hPP_KPlusPt["KplusPP"]->scaleW(s);
      hPP_KMinusPt["KminusPP"]->scaleW(s);
      hPP_PPt["PPP"]->scaleW(s);
      hPP_P_BarPt["P_barPP"]->scaleW(s);

      //hPionPlusPt->scaleW(1./sow->sumW());
      //hPionMinusPt->scaleW(1./sow->sumW());
    }

    ///@}

    map<string, Histo1DPtr> hdAu_PiPlusPt;
    map<string, Histo1DPtr> hdAu_PiMinusPt;
    map<string, Histo1DPtr> hdAu_KPlusPt;
    map<string, Histo1DPtr> hdAu_KMinusPt;
    map<string, Histo1DPtr> hdAu_PPt;
    map<string, Histo1DPtr> hdAu_P_BarPt;

    map<string, Histo1DPtr> hPP_PiPlusPt;
    map<string, Histo1DPtr> hPP_PiMinusPt;
    map<string, Histo1DPtr> hPP_KPlusPt;
    map<string, Histo1DPtr> hPP_KMinusPt;
    map<string, Histo1DPtr> hPP_PPt;
    map<string, Histo1DPtr> hPP_P_BarPt;

    //map<string, CounterPtr> sow;
    string beamOpt;
    enum CollisionSystem { pp, dAu200 };
    CollisionSystem collSys;

    /// @name Histograms
    ///@{
    //Histo1DPtr hPionPlusPt;
    //Histo1DPtr hPionMinusPt;
    CounterPtr sow;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2006_I711951);

}
