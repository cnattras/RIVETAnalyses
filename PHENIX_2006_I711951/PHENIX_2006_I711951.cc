// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
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
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2006_I711951);

    bool getDeltaPt(YODA::Histo1D hist, double pT, double &deltaPt)
    {
      if(pT > hist.xMin() && pT < hist.xMax())
      {
        deltaPt = hist.binAt(pT).xMid();
        return true;
      }
      else return false;

    }
//create binShift function
    void binShift(YODA::Histo1D& histogram) {
        const auto& binlist = histogram.bins();
        const auto& lastBin = histogram.bin(histogram.numBins());
        int n = 0;
        for (auto& bins : binlist) {
            double p_high = bins.xMax();
            double p_low = bins.xMin();
            //Now calculate f_corr
            if (bins.xMin() == binlist[0].xMin()) { //Check if we are working with first bin
                float b = 1 / (p_high - p_low) * log(binlist[0].sumW()/binlist[1].sumW());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
                n += 1;
            } else if (bins.xMin() == lastBin.xMin()){ //Check if we are working with last bin
                float b = 1 / (p_high - p_low) * log(binlist[binlist.size()-2].sumW() / lastBin.sumW());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
            } else { //Check if we are working with any middle bin
                float b = 1 / (p_high - p_low) * log(binlist[n-1].sumW() / binlist[n+1].sumW());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
                n += 1;
            }
        }
    }
    /// Book histograms and initialise projections before the run
    void init() {
      //Particles: pi^+, pi^-, k^+, k^-, p, p_bar
      //std::initializer_list<int> pdgIds = {221};

    // Alice projection? 
      const ALICE::PrimaryParticles cp(Cuts::absrap < 0.5);
      declare(cp,"cp");
   
    
      const ParticlePair& beam = beams();

      
      beamOpt = getOption<string>("beam","NONE");
      /*if (beamOpt =="pp") collSys = pp;
      else if (beamOpt == "DAu200") collSys = DAu200;
      else if (beamOpt == "AuAu200") collSys = AuAU200;*/


      if (beamOpt == "NONE") {
      if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
      {
        float NN = 1.;
        if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5e-3)) collSys = pp;
      }
      else if ((beam.first.pid() == 1000010020 && beam.second.pid() == 1000791970) || (beam.first.pid() == 1000791970 && beam.second.pid() == 1000010020))
      {
        if (fuzzyEquals(sqrtS()/GeV, 200*sqrt(197*2), 5e-3)) collSys = DAu200;
      }
      else if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000010020) collSys = DAu200;
      {
        if (fuzzyEquals(sqrtS()/GeV, 200*sqrt(197*2), 5e-3)) collSys = DAu200;
      }
      }


      if (beamOpt =="PP200") collSys = pp;
      else if (beamOpt == "DAU200") collSys = DAu200;
      
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
      


      //****Counters****
      book(sow["sow_pp"],"_sow_pp");
      book(sow["sow_DAu20"],"_sow_DAu20");
      book(sow["sow_DAu40"],"_sow_DAu40");
      book(sow["sow_DAu60"],"_sow_DAu60");
      book(sow["sow_DAu88"],"_sow_DAu88");
      //book(sow,"sow");

      //****Book Histos****

      //for (i=0, i<4, i++){

      //Invariant Yield in pp (piplus)
      book(hPP_Yields["PiplusPP"], 2,1,1);

      //Invariant Yield in pp (piminus)
      book(hPP_Yields["PiminusPP"], 2,1,2);

      //Invariant Yield in pp (kplus)
      book(hPP_Yields["KplusPP"], 4,1,1);

      //Invariant Yield in pp (kminus)
      book(hPP_Yields["KminusPP"], 4,1,2);

      //Invariant Yield in pp (p)
      book(hPP_Yields["PPP"], 6,1,1);

      //Invariant Yield in pp (p_bar)
      book(hPP_Yields["P_barPP"], 6,1,2);

      //Invariant Yield in DAu 0-20% (piplus)
      book(hDAu_Yields["PiplusC20"], 2,1,5);

      //Invariant Yield in DAu 0-20% (piminus)
      book(hDAu_Yields["PiminusC20"], 2,1,6);

      //Invariant Yield in DAu 0-20% (kplus)
      book(hDAu_Yields["KplusC20"], 4,1,5);

      //Invariant Yield in DAu 0-20% (kminus)
      book(hDAu_Yields["KminusC20"], 4,1,6);

      //Invariant Yield in DAu 0-20% (p)
      book(hDAu_Yields["PC20"], 6,1,5);

      //Invariant Yield in DAu 0-20% (p_bar)
      book(hDAu_Yields["P_barC20"], 6,1,6);

      //Invariant Yield in DAu 20-40% (piplus)
      book(hDAu_Yields["PiplusC40"], 2,1,7);

      //Invariant Yield in DAu 20-40% (piminus)
      book(hDAu_Yields["PiminusC40"], 2,1,8);

      //Invariant Yield in DAu 20-40% (kplus)
      book(hDAu_Yields["KplusC40"], 5,1,1);

      //Invariant Yield in DAu 20-40% (kminus)
      book(hDAu_Yields["KminusC40"], 5,1,2);

      //Invariant Yield in DAu 20-40% (p)
      book(hDAu_Yields["PC40"], 6,1,7);

      //Invariant Yield in DAu 20-40% (p_bar)
      book(hDAu_Yields["P_barC40"], 6,1,8);

      //Invariant Yield in DAu 40-60% (piplus)
      book(hDAu_Yields["PiplusC60"], 3,1,1);

      //Invariant Yield in DAu 40-60% (piminus)
      book(hDAu_Yields["PiminusC60"], 2,1,9);

      //Invariant Yield in DAu 40-60% (kplus)
      book(hDAu_Yields["KplusC60"], 4,1,7);

      //Invariant Yield in DAu 40-60% (kminus)
      book(hDAu_Yields["KminusC60"], 4,1,8);

      //Invariant Yield in DAu 40-60% (p)
      book(hDAu_Yields["PC60"], 6,1,9);

      //Invariant Yield in DAu 40-60% (p_bar)
      book(hDAu_Yields["P_barC60"], 6,1,10);

      //Invariant Yield in DAu 60-88.5% (piplus)
      book(hDAu_Yields["PiplusC88"], 2,1,10);

      //Invariant Yield in DAu 60-88.5% (piminus)
      book(hDAu_Yields["PiminusC88"], 2,1,11);

      //Invariant Yield in DAu 60-88.5% (kplus)
      book(hDAu_Yields["KplusC88"], 4,1,9);

      //Invariant Yield in DAu 60-88.5% (kminus)
      book(hDAu_Yields["KminusC88"], 4,1,10);

      //Invariant Yield in DAu 60-88.5% (p)
      book(hDAu_Yields["PC88"], 6,1,11);

      //Invariant Yield in DAu 60-88.5% (p_bar)
      book(hDAu_Yields["P_barC88"], 6,1,12);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      //sow->fill();
      Particles chargedP = apply<PrimaryParticles>(event,"cp").particles();

      
        
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
            if (!(p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) || // Lambda, Anti-Lambda
                 p.hasAncestorWith(Cuts::pid == 3212) || p.hasAncestorWith(Cuts::pid == -3212) || // Sigma0
                 p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) || // Sigmas
                 p.hasAncestorWith(Cuts::pid == 3112) || p.hasAncestorWith(Cuts::pid == -3112) || // Sigmas
                 p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) || // Cascades
                 p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) || // Cascades
                 p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334))) { // Omegas                        //Omega-
                    if (getDeltaPt(*hPP_Yields["PPP"], partPt, deltaPt))
                    {
                      pt_weight /= deltaPt;
                      hPP_Yields["PPP"]->fill(partPt, pt_weight);
                    }
                  }
          }
          if (p.pid() == -2212) {
            if (!(p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) || // Lambda, Anti-Lambda
                 p.hasAncestorWith(Cuts::pid == 3212) || p.hasAncestorWith(Cuts::pid == -3212) || // Sigma0
                 p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) || // Sigmas
                 p.hasAncestorWith(Cuts::pid == 3112) || p.hasAncestorWith(Cuts::pid == -3112) || // Sigmas
                 p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) || // Cascades
                 p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) || // Cascades
                 p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334))) { // Omegas
                    if (getDeltaPt(*hPP_Yields["P_barPP"], partPt, deltaPt))
                    {
                      pt_weight /= deltaPt;
                      hPP_Yields["P_barPP"]->fill(partPt, pt_weight);
                    }
                  }
          }
        }
      }

      if (collSys == DAu200) {

        const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
        const double cent = centProj();

        if (cent < 0. || cent > 88.5) vetoEvent;

        if (cent > 0. && cent < 20.) {
          sow["sow_DAu20"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hDAu_Yields["PiplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["PiplusC20"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hDAu_Yields["PiminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["PiminusC20"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hDAu_Yields["KplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["KplusC20"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -321) {
              if (getDeltaPt(*hDAu_Yields["KminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["KminusC20"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) || // Lambda, Anti-Lambda
                   p.hasAncestorWith(Cuts::pid == 3212) || p.hasAncestorWith(Cuts::pid == -3212) || // Sigma0
                   p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3112) || p.hasAncestorWith(Cuts::pid == -3112) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334))) { // Omegas
                      if (getDeltaPt(*hDAu_Yields["PC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hDAu_Yields["PC20"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) || // Lambda, Anti-Lambda
                   p.hasAncestorWith(Cuts::pid == 3212) || p.hasAncestorWith(Cuts::pid == -3212) || // Sigma0
                   p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3112) || p.hasAncestorWith(Cuts::pid == -3112) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334))) { // Omegas
                      if (getDeltaPt(*hDAu_Yields["P_barC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hDAu_Yields["P_barC20"]->fill(partPt, pt_weight);
                      }
                    }
            }
          }

        }

        if (cent > 20. && cent < 40.) {
          sow["sow_DAu40"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hDAu_Yields["PiplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["PiplusC40"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hDAu_Yields["PiminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["PiminusC40"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hDAu_Yields["KplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["KplusC40"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -321) {
              if (getDeltaPt(*hDAu_Yields["KminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["KminusC40"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) || // Lambda, Anti-Lambda
                   p.hasAncestorWith(Cuts::pid == 3212) || p.hasAncestorWith(Cuts::pid == -3212) || // Sigma0
                   p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3112) || p.hasAncestorWith(Cuts::pid == -3112) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334))) { // Omegas
                      if (getDeltaPt(*hDAu_Yields["PC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hDAu_Yields["PC40"]->fill(partPt, pt_weight);
                      }
              }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) || // Lambda, Anti-Lambda
                   p.hasAncestorWith(Cuts::pid == 3212) || p.hasAncestorWith(Cuts::pid == -3212) || // Sigma0
                   p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3112) || p.hasAncestorWith(Cuts::pid == -3112) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334))) { // Omegas
                      if (getDeltaPt(*hDAu_Yields["P_barC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hDAu_Yields["P_barC40"]->fill(partPt, pt_weight);
                      }
              }
            }
          }


        }

        if (cent > 40. && cent < 60.) {
          sow["sow_DAu60"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hDAu_Yields["PiplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["PiplusC60"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hDAu_Yields["PiminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["PiminusC60"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hDAu_Yields["KplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["KplusC60"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -321) {
              if (getDeltaPt(*hDAu_Yields["KminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["KminusC60"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) || // Lambda, Anti-Lambda
                   p.hasAncestorWith(Cuts::pid == 3212) || p.hasAncestorWith(Cuts::pid == -3212) || // Sigma0
                   p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3112) || p.hasAncestorWith(Cuts::pid == -3112) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334))) { // Omegas
                      if (getDeltaPt(*hDAu_Yields["PC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hDAu_Yields["PC60"]->fill(partPt, pt_weight);
                      }
              }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) || // Lambda, Anti-Lambda
                   p.hasAncestorWith(Cuts::pid == 3212) || p.hasAncestorWith(Cuts::pid == -3212) || // Sigma0
                   p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3112) || p.hasAncestorWith(Cuts::pid == -3112) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334))) { // Omegas
                      if (getDeltaPt(*hDAu_Yields["P_barC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hDAu_Yields["P_barC60"]->fill(partPt, pt_weight);
                      }
              }
            }
          }

        }

        if (cent > 60. && cent < 88.5) {
          sow["sow_DAu88"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (partPt*2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hDAu_Yields["PiplusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["PiplusC88"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hDAu_Yields["PiminusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["PiminusC88"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hDAu_Yields["KplusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["KplusC88"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -321) {
              if (getDeltaPt(*hDAu_Yields["KminusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hDAu_Yields["KminusC88"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) || // Lambda, Anti-Lambda
                   p.hasAncestorWith(Cuts::pid == 3212) || p.hasAncestorWith(Cuts::pid == -3212) || // Sigma0
                   p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3112) || p.hasAncestorWith(Cuts::pid == -3112) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334))) { // Omegas
                      if (getDeltaPt(*hDAu_Yields["PC88"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hDAu_Yields["PC88"]->fill(partPt, pt_weight);
                      }
              }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) || // Lambda, Anti-Lambda
                   p.hasAncestorWith(Cuts::pid == 3212) || p.hasAncestorWith(Cuts::pid == -3212) || // Sigma0
                   p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3112) || p.hasAncestorWith(Cuts::pid == -3112) || // Sigmas
                   p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) || // Cascades
                   p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334))) { // Omegas 
                      if (getDeltaPt(*hDAu_Yields["P_barC88"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hDAu_Yields["P_barC88"]->fill(partPt, pt_weight);
                      }
              }
            }
          }

        }

      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
//d01-x01-y01
        binShift(*hPP_Yields["PiplusPP"]);
	
	//d01-x01-y01
	
	      binShift(*hPP_Yields["PiminusPP"]);

	//d05-x01-y01
	
	      binShift(*hPP_Yields["KplusPP"]);
	
	//d05-x01-y02
	
	      binShift(*hPP_Yields["KminusPP"]);

	//d09-x01-y01
	
	      binShift(*hPP_Yields["PPP"]);
	
	//d09-x01-y02
	
	      binShift(*hPP_Yields["P_barPP"]);
	
	//d03-x01-y01
	
	      binShift(*hDAu_Yields["PiplusC20"]);
	
	//d03-x01-y02
	
	      binShift(*hDAu_Yields["PiminusC20"]);
	
	//d07-x01-y01
	
	      binShift(*hDAu_Yields["KplusC20"]);
	
	
	//d07-x01-y02
	
	      binShift(*hDAu_Yields["KminusC20"]);
	
	//d10-x01-y01
	
	      binShift(*hDAu_Yields["PC20"]);

	//d10-x01-y02
	
	      binShift(*hDAu_Yields["P_barC20"]);

	//d03-x01-y03

	      binShift(*hDAu_Yields["PiplusC40"]);
	
	
	//d03-x01-y04

	      binShift(*hDAu_Yields["PiminusC40"]);

	//d08-x01-y01

	      binShift(*hDAu_Yields["KplusC40"]);


	//d08-x01-y02

	      binShift(*hDAu_Yields["KminusC40"]);

	//d10-x01-y04

	      binShift(*hDAu_Yields["P_barC40"]);

	//d04-x01-y01
	
	      binShift(*hDAu_Yields["PiplusC60"]);

	//d03-x01-y05
	
	      binShift(*hDAu_Yields["PiminusC60"]);

	//d07-x01-y03
	
	      binShift(*hDAu_Yields["KplusC60"]);

	//d07-x01-y04
	
	      binShift(*hDAu_Yields["KminusC60"]);

	//d10-x01-y05
	
	      binShift(*hDAu_Yields["PC60"]);

	//d10-x01-y06
	
	      binShift(*hDAu_Yields["P_barC60"]);

	//d03-x01-y06
	
	      binShift(*hDAu_Yields["PiplusC88"]);

	//d03-x01-y07
  
	      binShift(*hDAu_Yields["PiminusC88"]);

	//d07-x01-y05
	
	      binShift(*hDAu_Yields["KplusC88"]);

	//d07-x01-y06
	
	      binShift(*hDAu_Yields["KminusC88"]);

	//d10-x01-y07
	
	      binShift(*hDAu_Yields["PC88"]);

	//d10-x01-y08
	
	      binShift(*hDAu_Yields["P_barC88"]);


      double xs = crossSection()/barn;
      double sf = 1.0 / (2 * M_PI );
      //****Scale Histos****
      hDAu_Yields["PiplusC20"]->scaleW(1./sow["sow_DAu20"]->sumW());
      hDAu_Yields["PiminusC20"]->scaleW(1./sow["sow_DAu20"]->sumW());
      hDAu_Yields["KplusC20"]->scaleW(1./sow["sow_DAu20"]->sumW());
      hDAu_Yields["KminusC20"]->scaleW(1./sow["sow_DAu20"]->sumW());
      hDAu_Yields["PC20"]->scaleW(1./sow["sow_DAu20"]->sumW());
      hDAu_Yields["P_barC20"]->scaleW(1./sow["sow_DAu20"]->sumW());

      hDAu_Yields["PiplusC40"]->scaleW(1./sow["sow_DAu40"]->sumW());
      hDAu_Yields["PiminusC40"]->scaleW(1./sow["sow_DAu40"]->sumW());
      hDAu_Yields["KplusC40"]->scaleW(1./sow["sow_DAu40"]->sumW());
      hDAu_Yields["KminusC40"]->scaleW(1./sow["sow_DAu40"]->sumW());
      hDAu_Yields["PC40"]->scaleW(1./sow["sow_DAu20"]->sumW());
      hDAu_Yields["P_barC40"]->scaleW(1./sow["sow_DAu40"]->sumW());

      hDAu_Yields["PiplusC60"]->scaleW(1./sow["sow_DAu60"]->sumW());
      hDAu_Yields["PiminusC60"]->scaleW(1./sow["sow_DAu60"]->sumW());
      hDAu_Yields["KplusC60"]->scaleW(1./sow["sow_DAu60"]->sumW());
      hDAu_Yields["KminusC60"]->scaleW(1./sow["sow_DAu60"]->sumW());
      hDAu_Yields["PC60"]->scaleW(1./sow["sow_DAu60"]->sumW());
      hDAu_Yields["P_barC60"]->scaleW(1./sow["sow_DAu60"]->sumW());

      hDAu_Yields["PiplusC88"]->scaleW(1./sow["sow_DAu88"]->sumW());
      hDAu_Yields["PiminusC88"]->scaleW(1./sow["sow_DAu88"]->sumW());
      hDAu_Yields["KplusC88"]->scaleW(1./sow["sow_DAu88"]->sumW());
      hDAu_Yields["KminusC88"]->scaleW(1./sow["sow_DAu88"]->sumW());
      hDAu_Yields["PC88"]->scaleW(1./sow["sow_DAu88"]->sumW());
      hDAu_Yields["P_barC88"]->scaleW(1./sow["sow_DAu88"]->sumW());

      hPP_Yields["PiplusPP"]->scaleW(xs*sf*1./sow["sow_pp"]->sumW());
      hPP_Yields["PiminusPP"]->scaleW(xs*sf*1./sow["sow_pp"]->sumW());
      hPP_Yields["KplusPP"]->scaleW(xs*sf*1./sow["sow_pp"]->sumW());
      hPP_Yields["KminusPP"]->scaleW(xs*sf*1./sow["sow_pp"]->sumW());
      hPP_Yields["PPP"]->scaleW(xs*sf*1./sow["sow_pp"]->sumW());
      hPP_Yields["P_barPP"]->scaleW(sf*1./sow["sow_pp"]->sumW());



    }

    ///@}

    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> hDAu_Yields;
    map<string, Histo1DPtr> hPP_Yields;

    map<string, CounterPtr> sow;
    string beamOpt;
    enum CollisionSystem { NONE, pp, AuAu200, DAu200 };
    CollisionSystem collSys;

    ///@}

  };


  RIVET_DECLARE_PLUGIN(PHENIX_2006_I711951);

}
