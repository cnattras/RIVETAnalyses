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
#include "../Centralities/RHICCentrality.hh"
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

    bool getDeltaPt(YODA::Histo1D hist, double pT, double &deltaPt)
    {
      if(pT > hist.xMin() && pT < hist.xMax())
      {
        deltaPt = hist.bin(hist.binIndexAt(pT)).xMid();
        return true;
      }
      else return false;

    }

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
      book(sow["sow_pp"],"sow_pp");
      book(sow["sow_dAU20"],"sow_dAU20");
      book(sow["sow_dAU40"],"sow_dAU40");
      book(sow["sow_dAU60"],"sow_dAU60");
      book(sow["sow_dAU88"],"sow_dAU88");
      //book(sow,"sow");

      //****Book Histos****

      //for (i=0, i<4, i++){

      //Invariant Yield in pp (piplus)
      book(hPP_Yields["PiplusPP"], 1,1,1);

      //Invariant Yield in pp (piminus)
      book(hPP_Yields["PiminusPP"], 1,1,2);

      //Invariant Yield in pp (kplus)
      book(hPP_Yields["KplusPP"], 5,1,1);

      //Invariant Yield in pp (kminus)
      book(hPP_Yields["KminusPP"], 5,1,2);

      //Invariant Yield in pp (p)
      book(hPP_Yields["PPP"], 9,1,1);

      //Invariant Yield in pp (p_bar)
      book(hPP_Yields["P_barPP"], 9,1,2);

      //Invariant Yield in dAU 0-20% (piplus)
      book(hdAu_Yields["PiplusC20"], 3,1,1);

      //Invariant Yield in dAU 0-20% (piminus)
      book(hdAu_Yields["PiminusC20"], 3,1,2);

      //Invariant Yield in dAU 0-20% (kplus)
      book(hdAu_Yields["KplusC20"], 7,1,1);

      //Invariant Yield in dAU 0-20% (kminus)
      book(hdAu_Yields["KminusC20"], 7,1,2);

      //Invariant Yield in dAU 0-20% (p)
      book(hdAu_Yields["PC20"], 10,1,1);

      //Invariant Yield in dAU 0-20% (p_bar)
      book(hdAu_Yields["P_barC20"], 10,1,2);

      //Invariant Yield in dAU 20-40% (piplus)
      book(hdAu_Yields["PiplusC40"], 3,1,3);

      //Invariant Yield in dAU 20-40% (piminus)
      book(hdAu_Yields["PiminusC40"], 3,1,4);

      //Invariant Yield in dAU 20-40% (kplus)
      book(hdAu_Yields["KplusC40"], 8,1,1);

      //Invariant Yield in dAU 20-40% (kminus)
      book(hdAu_Yields["KminusC40"], 8,1,2);

      //Invariant Yield in dAU 20-40% (p)
      book(hdAu_Yields["PC40"], 10,1,3);

      //Invariant Yield in dAU 20-40% (p_bar)
      book(hdAu_Yields["P_barC40"], 10,1,4);

      //Invariant Yield in dAU 40-60% (piplus)
      book(hdAu_Yields["PiplusC60"], 4,1,1);

      //Invariant Yield in dAU 40-60% (piminus)
      book(hdAu_Yields["PiminusC60"], 3,1,5);

      //Invariant Yield in dAU 40-60% (kplus)
      book(hdAu_Yields["KplusC60"], 7,1,3);

      //Invariant Yield in dAU 40-60% (kminus)
      book(hdAu_Yields["KminusC60"], 7,1,4);

      //Invariant Yield in dAU 40-60% (p)
      book(hdAu_Yields["PC60"], 10,1,5);

      //Invariant Yield in dAU 40-60% (p_bar)
      book(hdAu_Yields["P_barC60"], 10,1,6);

      //Invariant Yield in dAU 60-88.5% (piplus)
      book(hdAu_Yields["PiplusC88"], 3,1,6);

      //Invariant Yield in dAU 60-88.5% (piminus)
      book(hdAu_Yields["PiminusC88"], 3,1,7);

      //Invariant Yield in dAU 60-88.5% (kplus)
      book(hdAu_Yields["KplusC88"], 7,1,5);

      //Invariant Yield in dAU 60-88.5% (kminus)
      book(hdAu_Yields["KminusC88"], 7,1,6);

      //Invariant Yield in dAU 60-88.5% (p)
      book(hdAu_Yields["PC88"], 10,1,7);

      //Invariant Yield in dAU 60-88.5% (p_bar)
      book(hdAu_Yields["P_barC88"], 10,1,8);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      //sow->fill();
      Particles chargedP = applyProjection<FinalState>(event,"fs").particles();

      if (collSys == pp) {
        sow["sow_pp"]->fill();
        for (const Particle& p :chargedP) {
          double partPt = p.pT() / GeV;
          double pt_weight = 1. / (partPt*2.*M_PI);
          double deltaPt = 0.;

          if (p.pid() == 211) {
            if (getDeltaPt(*hPP_Yields["PiplusPP"], partPt, deltaPt))
            {
              pt_weight /= deltaPt;
              hPP_Yields["PiplusPP"]->fill(partPt, pt_weight);
            }
          }
          if (p.pid() == -211) {
            if (getDeltaPt(*hPP_Yields["PiminusPP"], partPt, deltaPt))
            {
              pt_weight /= deltaPt;
              hPP_Yields["PiminusPP"]->fill(partPt, pt_weight);
            }
          }
          if (p.pid() == 321) {
            if (getDeltaPt(*hPP_Yields["KplusPP"], partPt, deltaPt))
            {
              pt_weight /= deltaPt;
              hPP_Yields["KplusPP"]->fill(partPt, pt_weight);
            }
          }
          if (p.pid() == -321) {
            if (getDeltaPt(*hPP_Yields["KminusPP"], partPt, deltaPt))
            {
              pt_weight /= deltaPt;
              hPP_Yields["KminusPP"]->fill(partPt, pt_weight);
            }
          }
          if (p.pid() == 2212) {
            if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                  p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                  p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                  p.hasAncestor(3334) )) {                        //Omega-
                    if (getDeltaPt(*hPP_Yields["PPP"], partPt, deltaPt))
                    {
                      pt_weight /= deltaPt;
                      hPP_Yields["PPP"]->fill(partPt, pt_weight);
                    }
                  }
          }
          if (p.pid() == -2212) {
            if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                  p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                  p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                  p.hasAncestor(3334) )) {                        //Omega-
                    if (getDeltaPt(*hPP_Yields["P_barPP"], partPt, deltaPt))
                    {
                      pt_weight /= deltaPt;
                      hPP_Yields["P_barPP"]->fill(partPt, pt_weight);
                    }
                  }
          }
        }
      }

      if (collSys == dAu200) {

        const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
        const double cent = centProj();

        if (cent < 0. || cent > 88.5) vetoEvent;

        if (cent > 0. && cent < 20.) {
          sow["sow_dAU20"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hdAu_Yields["PiplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC20"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hdAu_Yields["PiminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC20"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hdAu_Yields["KplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC20"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -321) {
              if (getDeltaPt(*hdAu_Yields["KminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC20"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hdAu_Yields["PC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC20"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hdAu_Yields["P_barC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC20"]->fill(partPt, pt_weight);
                      }
                    }
            }
          }

        }

        if (cent > 20. && cent < 40.) {
          sow["sow_dAU40"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hdAu_Yields["PiplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC40"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hdAu_Yields["PiminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC40"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hdAu_Yields["KplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC40"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -321) {
              if (getDeltaPt(*hdAu_Yields["KminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC40"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hdAu_Yields["PC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC40"]->fill(partPt, pt_weight);
                      }
              }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hdAu_Yields["P_barC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC40"]->fill(partPt, pt_weight);
                      }
              }
            }
          }


        }

        if (cent > 40. && cent < 60.) {
          sow["sow_dAU60"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hdAu_Yields["PiplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC60"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hdAu_Yields["PiminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC60"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hdAu_Yields["KplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC60"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -321) {
              if (getDeltaPt(*hdAu_Yields["KminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC60"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hdAu_Yields["PC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC60"]->fill(partPt, pt_weight);
                      }
              }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hdAu_Yields["P_barC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC60"]->fill(partPt, pt_weight);
                      }
              }
            }
          }

        }

        if (cent > 60. && cent < 88.5) {
          sow["sow_dAU88"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hdAu_Yields["PiplusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC88"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hdAu_Yields["PiminusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC88"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hdAu_Yields["KplusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC88"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -321) {
              if (getDeltaPt(*hdAu_Yields["KminusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC88"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hdAu_Yields["PC88"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC88"]->fill(partPt, pt_weight);
                      }
              }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hdAu_Yields["P_barC88"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC88"]->fill(partPt, pt_weight);
                      }
              }
            }
          }

        }

      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
//
      bool dAu200_available = false;
      bool pp_available = false;

      for (auto element : hdAu_Yields)
      {
        string name = element.second->name();
        if (name.find("pp") != std::string::npos)
        {
          if (element.second->numEntries() > 0) pp_available = true;
          else
          {
            pp_available = false;
            break;
          }
        }
        else if (name.find("C") != std::string::npos)
        {
          if (element.second->numEntries() > 0) dAu200_available = true;
          else
          {
            dAu200_available = false;
            break;
          }
        }
      }

      for (auto element : hPP_Yields)
      {
        string name = element.second->name();
        if (name.find("pp") != std::string::npos)
        {
          if (element.second->numEntries() > 0) pp_available = true;
          else
          {
            pp_available = false;
            break;
          }
        }
        else if (name.find("C") != std::string::npos)
        {
          if (element.second->numEntries() > 0) dAu200_available = true;
          else
          {
            dAu200_available = false;
            break;
          }
        }
      }

      //if (!(dAu200_available && pp_available)) return;
//

      //****Scale Histos****

      hdAu_Yields["PiplusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["PiminusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["KplusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["KminusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["PC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["P_barC20"]->scaleW(1./sow["sow_dAU20"]->sumW());

      hdAu_Yields["PiplusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_Yields["PiminusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_Yields["KplusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_Yields["KminusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_Yields["PC40"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["P_barC40"]->scaleW(1./sow["sow_dAU40"]->sumW());

      hdAu_Yields["PiplusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_Yields["PiminusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_Yields["KplusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_Yields["KminusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_Yields["PC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_Yields["P_barC60"]->scaleW(1./sow["sow_dAU60"]->sumW());

      hdAu_Yields["PiplusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_Yields["PiminusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_Yields["KplusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_Yields["KminusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_Yields["PC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_Yields["P_barC88"]->scaleW(1./sow["sow_dAU88"]->sumW());

      hPP_Yields["PiplusPP"]->scaleW(1./sow["sow_pp"]->sumW());
      hPP_Yields["PiminusPP"]->scaleW(1./sow["sow_pp"]->sumW());
      hPP_Yields["KplusPP"]->scaleW(1./sow["sow_pp"]->sumW());
      hPP_Yields["KminusPP"]->scaleW(1./sow["sow_pp"]->sumW());
      hPP_Yields["PPP"]->scaleW(1./sow["sow_pp"]->sumW());
      hPP_Yields["P_barPP"]->scaleW(1./sow["sow_pp"]->sumW());


/*
      const double s = 1./sow->sumW();
      hdAu_Yields["PiplusC20"]->scaleW(s);
      hdAu_Yields["PiminusC20"]->scaleW(s);
      hdAu_Yields["KplusC20"]->scaleW(s);
      hdAu_Yields["KminusC20"]->scaleW(s);
      hdAu_Yields["PC20"]->scaleW(s);
      hdAu_Yields["P_barC20"]->scaleW(s);

      hdAu_Yields["PiplusC40"]->scaleW(s);
      hdAu_Yields["PiminusC40"]->scaleW(s);
      hdAu_Yields["KplusC40"]->scaleW(s);
      hdAu_Yields["KminusC40"]->scaleW(s);
      hdAu_Yields["PC40"]->scaleW(s);
      hdAu_Yields["P_barC40"]->scaleW(s);

      hdAu_Yields["PiplusC60"]->scaleW(s);
      hdAu_Yields["PiminusC60"]->scaleW(s);
      hdAu_Yields["KplusC60"]->scaleW(s);
      hdAu_Yields["KminusC60"]->scaleW(s);
      hdAu_Yields["PC60"]->scaleW(s);
      hdAu_Yields["P_barC60"]->scaleW(s);

      hdAu_Yields["PiplusC88"]->scaleW(s);
      hdAu_Yields["PiminusC88"]->scaleW(s);
      hdAu_Yields["KplusC88"]->scaleW(s);
      hdAu_Yields["KminusC88"]->scaleW(s);
      hdAu_Yields["PC88"]->scaleW(s);
      hdAu_Yields["P_barC88"]->scaleW(s);

      hPP_Yields["PiplusPP"]->scaleW(s);
      hPP_Yields["PiminusPP"]->scaleW(s);
      hPP_Yields["KplusPP"]->scaleW(s);
      hPP_Yields["KminusPP"]->scaleW(s);
      hPP_Yields["PPP"]->scaleW(s);
      hPP_Yields["P_barPP"]->scaleW(s);

*/
    }

    ///@}

    map<string, Histo1DPtr> hdAu_Yields;
    map<string, Histo1DPtr> hPP_Yields;

/*  map<string, Histo1DPtr> hdAu_Yields;
    map<string, Histo1DPtr> hdAu_Yields;
    map<string, Histo1DPtr> hdAu_Yields;
    map<string, Histo1DPtr> hdAu_Yields;
    map<string, Histo1DPtr> hdAu_Yields;
    map<string, Histo1DPtr> hdAu_Yields;

    map<string, Histo1DPtr> hPP_Yields;
    map<string, Histo1DPtr> hPP_Yields;
    map<string, Histo1DPtr> hPP_Yields;
    map<string, Histo1DPtr> hPP_Yields;
    map<string, Histo1DPtr> hPP_Yields;
    map<string, Histo1DPtr> hPP_Yields;
*/

    map<string, CounterPtr> sow;
    string beamOpt;
    enum CollisionSystem { pp, dAu200 };
    CollisionSystem collSys;

    /// @name Histograms
    ///@{
    //Histo1DPtr hPionPlusPt;
    //Histo1DPtr hPionMinusPt;
    //CounterPtr sow;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2006_I711951);

}
