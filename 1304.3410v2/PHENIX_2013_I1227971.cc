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
  class PHENIX_2013_I1227971 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2013_I1227971);

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

      // Initialise and register projections

      //Particles: pi^+, pi^-, k^+, k^-, p, p_bar
      //std::initializer_list<int> pdgIds = {221};

      const FinalState fs(Cuts::absrap<0.35&&Cuts::abscharge>0);
      declare(fs,"fs");

      beamOpt = getOption<string>("beam","NONE");
      if (beamOpt == "dAU200") collSys = dAu200;
      else if (beamOpt == "AUAU200") collSys = AuAu200;

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      //****Counters****
      book(sow["sow_dAU20"],"sow_dAU20");
      book(sow["sow_dAU40"],"sow_dAU40");
      book(sow["sow_dAU60"],"sow_dAU60");
      book(sow["sow_dAU88"],"sow_dAU88");
      book(sow["sow_dAU100"], "sow_dAU100");
      book(sow["sow_AUAU10"],"sow_AUAU10");
      book(sow["sow_AUAU20"],"sow_AUAU20");
      book(sow["sow_AUAU40"],"sow_AUAU40");
      book(sow["sow_AUAU60"],"sow_AUAU60");
      book(sow["sow_AUAU92"],"sow_AUAU92");

      //invariant yield of kaons (AUAU)
      book(hAuAu_Yields["KminusC10"], 1,1,1);
      book(hAuAu_Yields["KminusC20"], 1,1,2);
      book(hAuAu_Yields["KminusC40"], 1,1,3);
      book(hAuAu_Yields["KminusC60"], 1,1,4);
      book(hAuAu_Yields["KminusC92"], 1,1,5);
      book(hAuAu_Yields["KplusC10"], 1,1,6);
      book(hAuAu_Yields["KplusC20"], 1,1,7);
      book(hAuAu_Yields["KplusC40"], 1,1,8);
      book(hAuAu_Yields["KplusC60"], 1,1,9);
      book(hAuAu_Yields["KplusC92"], 1,1,10);
      //invariant yield of kaons (dAU)
      book(hdAu_Yields["KminusC20"], 2,1,1);
      book(hdAu_Yields["KminusC100"], 2,1,2);
      book(hdAu_Yields["KminusC40"], 2,1,3);
      book(hdAu_Yields["KminusC60"], 2,1,4);
      book(hdAu_Yields["KminusC88"], 2,1,5);
      book(hdAu_Yields["KplusC20"], 2,1,6);
      book(hdAu_Yields["KplusC100"], 2,1,7);
      book(hdAu_Yields["KplusC40"], 2,1,8);
      book(hdAu_Yields["KplusC60"], 2,1,9);
      book(hdAu_Yields["KplusC88"], 2,1,10);
      //invariant yield of pions (AUAU)
      book(hAuAu_Yields["PiminusC10"], 3,1,1);
      book(hAuAu_Yields["PiminusC20"], 3,1,2);
      book(hAuAu_Yields["PiminusC40"], 3,1,3);
      book(hAuAu_Yields["PiminusC60"], 3,1,4);
      book(hAuAu_Yields["PiminusC92"], 3,1,5);
      book(hAuAu_Yields["PiplusC10"], 3,1,6);
      book(hAuAu_Yields["PiplusC20"], 3,1,7);
      book(hAuAu_Yields["PiplusC40"], 3,1,8);
      book(hAuAu_Yields["PiplusC60"], 3,1,9);
      book(hAuAu_Yields["PiplusC92"], 3,1,10);
      //invariant yield of pions (dAU)
      book(hdAu_Yields["PiminusC20"], 4,1,1);
      book(hdAu_Yields["PiminusC100"], 4,1,2);
      book(hdAu_Yields["PiminusC40"], 4,1,3);
      book(hdAu_Yields["PiminusC60"], 4,1,4);
      book(hdAu_Yields["PiminusC88"], 4,1,5);
      book(hdAu_Yields["PiplusC20"], 4,1,6);
      book(hdAu_Yields["PiplusC100"], 4,1,7);
      book(hdAu_Yields["PiplusC40"], 4,1,8);
      book(hdAu_Yields["PiplusC60"], 4,1,9);
      book(hdAu_Yields["PiplusC88"], 4,1,10);
      //invariant yield of protons (AUAU)
      book(hAuAu_Yields["P_barC10"], 5,1,1);
      book(hAuAu_Yields["P_barC20"], 5,1,2);
      book(hAuAu_Yields["P_barC40"], 5,1,3);
      book(hAuAu_Yields["P_barC60"], 5,1,4);
      book(hAuAu_Yields["P_barC92"], 5,1,5);
      book(hAuAu_Yields["PC10"], 5,1,6);
      book(hAuAu_Yields["PC20"], 5,1,7);
      book(hAuAu_Yields["PC40"], 5,1,8);
      book(hAuAu_Yields["PC60"], 5,1,9);
      book(hAuAu_Yields["PC92"], 5,1,10);
      //invariant yield of protons (dAU)
      book(hdAu_Yields["P_barC20"], 6,1,1);
      book(hdAu_Yields["P_barC100"], 6,1,2);
      book(hdAu_Yields["P_barC40"], 6,1,3);
      book(hdAu_Yields["P_barC60"], 6,1,4);
      book(hdAu_Yields["P_barC88"], 6,1,5);
      book(hdAu_Yields["PC20"], 6,1,6);
      book(hdAu_Yields["PC100"], 6,1,7);
      book(hdAu_Yields["PC40"], 6,1,8);
      book(hdAu_Yields["PC60"], 6,1,9);
      book(hdAu_Yields["PC88"], 6,1,10);

      //ratio of kaons (AUAU)
      book(hTemp_ratio_AuAu["KminusC10"], "KminusC10_AuAu", refData(7,1,1));
      book(hTemp_ratio_AuAu["KminusC20"], "KminusC20_AuAu", refData(7,1,2));
      book(hTemp_ratio_AuAu["KminusC40"], "KminusC40_AuAu", refData(7,1,3));
      book(hTemp_ratio_AuAu["KminusC60"], "KminusC60_AuAu", refData(7,1,4));
      book(hTemp_ratio_AuAu["KminusC92"], "KminusC92_AuAu", refData(7,1,5));
      book(hTemp_ratio_AuAu["KplusC10"], "KplusC10_AuAu", refData(7,1,1));
      book(hTemp_ratio_AuAu["KplusC20"], "KplusC20_AuAu", refData(7,1,2));
      book(hTemp_ratio_AuAu["KplusC40"], "KplusC40_AuAu", refData(7,1,3));
      book(hTemp_ratio_AuAu["KplusC60"], "KplusC60_AuAu", refData(7,1,4));
      book(hTemp_ratio_AuAu["KplusC92"], "KplusC92_AuAu", refData(7,1,5));
      book(RatioAuAu["KminustoKplusC10"], 7,1,1);
      book(RatioAuAu["KminustoKplusC20"], 7,1,2);
      book(RatioAuAu["KminustoKplusC40"], 7,1,3);
      book(RatioAuAu["KminustoKplusC60"], 7,1,4);
      book(RatioAuAu["KminustoKplusC92"], 7,1,5);
      //ratio of kaons (dAU)
      book(hTemp_ratio_dAu["KminusC20"], "KminusC20_dAu", refData(8,1,1));
      book(hTemp_ratio_dAu["KminusC100"], "KminusC100_dAu", refData(8,1,2));
      book(hTemp_ratio_dAu["KminusC40"], "KminusC40_dAu", refData(8,1,3));
      book(hTemp_ratio_dAu["KminusC60"], "KminusC60_dAu", refData(8,1,4));
      book(hTemp_ratio_dAu["KminusC88"], "KminusC88_dAu", refData(8,1,5));
      book(hTemp_ratio_dAu["KplusC20"], "KplusC20_dAu", refData(8,1,1));
      book(hTemp_ratio_dAu["KplusC100"], "KplusC100_dAu", refData(8,1,2));
      book(hTemp_ratio_dAu["KplusC40"], "KplusC40_dAu", refData(8,1,3));
      book(hTemp_ratio_dAu["KplusC60"], "KplusC60_dAu", refData(8,1,4));
      book(hTemp_ratio_dAu["KplusC88"], "KplusC88_dAu", refData(8,1,5));
      book(RatiodAu["KminustoKplusC20"], 8,1,1);
      book(RatiodAu["KminustoKplusC100"], 8,1,2);
      book(RatiodAu["KminustoKplusC40"], 8,1,3);
      book(RatiodAu["KminustoKplusC60"], 8,1,4);
      book(RatiodAu["KminustoKplusC88"], 8,1,5);
      //ratio of pions (AUAU)
      book(hTemp_ratio_AuAu["PiminusC10"], "PiminusC10_AuAu", refData(9,1,1));
      book(hTemp_ratio_AuAu["PiminusC20"], "PiminusC20_AuAu", refData(9,1,2));
      book(hTemp_ratio_AuAu["PiminusC40"], "PiminusC40_AuAu", refData(9,1,3));
      book(hTemp_ratio_AuAu["PiminusC60"], "PiminusC60_AuAu", refData(9,1,4));
      book(hTemp_ratio_AuAu["PiminusC92"], "PiminusC92_AuAu", refData(9,1,5));
      book(hTemp_ratio_AuAu["PiplusC10"], "PiplusC10_AuAu", refData(9,1,1));
      book(hTemp_ratio_AuAu["PiplusC20"], "PiplusC20_AuAu", refData(9,1,2));
      book(hTemp_ratio_AuAu["PiplusC40"], "PiplusC40_AuAu", refData(9,1,3));
      book(hTemp_ratio_AuAu["PiplusC60"], "PiplusC60_AuAu", refData(9,1,4));
      book(hTemp_ratio_AuAu["PiplusC92"], "PiplusC92_AuAu", refData(9,1,5));
      book(RatioAuAu["PiminustoPiplusC10"], 9,1,1);
      book(RatioAuAu["PiminustoPiplusC20"], 9,1,2);
      book(RatioAuAu["PiminustoPiplusC40"], 9,1,3);
      book(RatioAuAu["PiminustoPiplusC60"], 9,1,4);
      book(RatioAuAu["PiminustoPiplusC92"], 9,1,5);
      //ratio of pions (dAU)
      book(hTemp_ratio_dAu["PiminusC20"], "PiminusC20_dAu", refData(10,1,1));
      book(hTemp_ratio_dAu["PiminusC100"], "PiminusC100_dAu", refData(10,1,2));
      book(hTemp_ratio_dAu["PiminusC40"], "PiminusC40_dAu", refData(10,1,3));
      book(hTemp_ratio_dAu["PiminusC60"], "PiminusC60_dAu", refData(10,1,4));
      book(hTemp_ratio_dAu["PiminusC88"], "PiminusC88_dAu", refData(10,1,5));
      book(hTemp_ratio_dAu["PiplusC20"], "PiplusC20_dAu", refData(10,1,1));
      book(hTemp_ratio_dAu["PiplusC100"], "PiplusC100_dAu", refData(10,1,2));
      book(hTemp_ratio_dAu["PiplusC40"], "PiplusC40_dAu", refData(10,1,3));
      book(hTemp_ratio_dAu["PiplusC60"], "PiplusC60_dAu", refData(10,1,4));
      book(hTemp_ratio_dAu["PiplusC88"], "PiplusC88_dAu", refData(10,1,5));
      book(RatiodAu["PiminustoKplusC20"], 10,1,1);
      book(RatiodAu["PiminustoKplusC100"], 10,1,2);
      book(RatiodAu["PiminustoKplusC40"], 10,1,3);
      book(RatiodAu["PiminustoKplusC60"], 10,1,4);
      book(RatiodAu["PiminustoKplusC88"], 10,1,5);
      //ratio of protons (AUAU)
      book(hTemp_ratio_AuAu["P_barC10"], "P_barC10_AuAu", refData(11,1,1));
      book(hTemp_ratio_AuAu["P_barC20"], "P_barC20_AuAu", refData(11,1,2));
      book(hTemp_ratio_AuAu["P_barC40"], "P_barC40_AuAu", refData(11,1,3));
      book(hTemp_ratio_AuAu["P_barC60"], "P_barC60_AuAu", refData(11,1,4));
      book(hTemp_ratio_AuAu["P_barC92"], "P_barC92_AuAu", refData(11,1,5));
      book(hTemp_ratio_AuAu["PC10"], "PC10_AuAu", refData(11,1,1));
      book(hTemp_ratio_AuAu["PC20"], "PC20_AuAu", refData(11,1,2));
      book(hTemp_ratio_AuAu["PC40"], "PC40_AuAu", refData(11,1,3));
      book(hTemp_ratio_AuAu["PC60"], "PC60_AuAu", refData(11,1,4));
      book(hTemp_ratio_AuAu["PC92"], "PC92_AuAu", refData(11,1,5));
      book(RatioAuAu["P_bartoPC10"], 11,1,1);
      book(RatioAuAu["P_bartoPC20"], 11,1,2);
      book(RatioAuAu["P_bartoPC40"], 11,1,3);
      book(RatioAuAu["P_bartoPC60"], 11,1,4);
      book(RatioAuAu["P_bartoPC92"], 11,1,5);
      //ratio of protons (dAU)
      book(hTemp_ratio_dAu["P_barC20"], "P_barC20_dAu", refData(12,1,1));
      book(hTemp_ratio_dAu["P_barC100"], "P_barC100_dAu", refData(12,1,2));
      book(hTemp_ratio_dAu["P_barC40"], "P_barC40_dAu", refData(12,1,3));
      book(hTemp_ratio_dAu["P_barC60"], "P_barC60_dAu", refData(12,1,4));
      book(hTemp_ratio_dAu["P_barC88"], "P_barC88_dAu", refData(12,1,5));
      book(hTemp_ratio_dAu["PC20"], "PC20_dAu", refData(12,1,1));
      book(hTemp_ratio_dAu["PC100"], "PC100_dAu", refData(12,1,2));
      book(hTemp_ratio_dAu["PC40"], "PC40_dAu", refData(12,1,3));
      book(hTemp_ratio_dAu["PC60"], "PC60_dAu", refData(12,1,4));
      book(hTemp_ratio_dAu["PC88"], "PC88_dAu", refData(12,1,5));
      book(RatiodAu["P_bartoPC20"], 12,1,1);
      book(RatiodAu["P_bartoPC100"], 12,1,2);
      book(RatiodAu["P_bartoPC40"], 12,1,3);
      book(RatiodAu["P_bartoPC60"], 12,1,4);
      book(RatiodAu["P_bartoPC88"], 12,1,5);

      //ratio kaon/pion (AUAU)
      book(hTemp_ratio_AuAu["KC10"], "KC10_AuAu", refData(13,1,1));
      book(hTemp_ratio_AuAu["KC20"], "KC20_AuAu", refData(13,1,2));
      book(hTemp_ratio_AuAu["KC40"], "KC40_AuAu", refData(13,1,3));
      book(hTemp_ratio_AuAu["KC60"], "KC60_AuAu", refData(13,1,4));
      book(hTemp_ratio_AuAu["KC92"], "KC92_AuAu", refData(13,1,5));
      book(hTemp_ratio_AuAu["PiC10"], "PiC10_AuAu", refData(13,1,1));
      book(hTemp_ratio_AuAu["PiC20"], "PiC20_AuAu", refData(13,1,2));
      book(hTemp_ratio_AuAu["PiC40"], "PiC40_AuAu", refData(13,1,3));
      book(hTemp_ratio_AuAu["PiC60"], "PiC60_AuAu", refData(13,1,4));
      book(hTemp_ratio_AuAu["PiC92"], "PiC92_AuAu", refData(13,1,5));
      book(RatioAuAu["KtoPiC10"], 13,1,1);
      book(RatioAuAu["KtoPiC20"], 13,1,2);
      book(RatioAuAu["KtoPiC40"], 13,1,3);
      book(RatioAuAu["KtoPiC60"], 13,1,4);
      book(RatioAuAu["KtoPiC92"], 13,1,5);
      //ratio kaon/pion (dAU)
      book(hTemp_ratio_dAu["KC20"], "KC20_dAu", refData(14,1,1));
      book(hTemp_ratio_dAu["KC100"], "KC100_dAu", refData(14,1,2));
      book(hTemp_ratio_dAu["KC40"], "KC40_dAu", refData(14,1,3));
      book(hTemp_ratio_dAu["KC60"], "KC60_dAu", refData(14,1,4));
      book(hTemp_ratio_dAu["KC88"], "KC88_dAu", refData(14,1,5));
      book(hTemp_ratio_dAu["PiC20"], "PiC20_dAu", refData(14,1,1));
      book(hTemp_ratio_dAu["PiC100"], "PiC100_dAu", refData(14,1,2));
      book(hTemp_ratio_dAu["PiC40"], "PiC40_dAu", refData(14,1,3));
      book(hTemp_ratio_dAu["PiC60"], "PiC60_dAu", refData(14,1,4));
      book(hTemp_ratio_dAu["PiC88"], "PiC88_dAu", refData(14,1,5));
      book(RatiodAu["KtoPiC20"], 14,1,1);
      book(RatiodAu["KtoPiC100"], 14,1,2);
      book(RatiodAu["KtoPiC40"], 14,1,3);
      book(RatiodAu["KtoPiC60"], 14,1,4);
      book(RatiodAu["KtoPiC88"], 14,1,5);
      //ratio proton/pion (AUAU)
      book(hTemp_ratio_AuAu["PC10_2"], "PC10_2_AuAu", refData(15,1,1));
      book(hTemp_ratio_AuAu["PC20_2"], "PC20_2_AuAu", refData(15,1,2));
      book(hTemp_ratio_AuAu["PC40_2"], "PC40_2_AuAu", refData(15,1,3));
      book(hTemp_ratio_AuAu["PC60_2"], "PC60_2_AuAu", refData(15,1,4));
      book(hTemp_ratio_AuAu["PC92_2"], "PC92_2_AuAu", refData(15,1,5));
      book(hTemp_ratio_AuAu["PiC10_2"], "PiC10_2_AuAu", refData(15,1,1));
      book(hTemp_ratio_AuAu["PiC20_2"], "PiC20_2_AuAu", refData(15,1,2));
      book(hTemp_ratio_AuAu["PiC40_2"], "PiC40_2_AuAu", refData(15,1,3));
      book(hTemp_ratio_AuAu["PiC60_2"], "PiC60_2_AuAu", refData(15,1,4));
      book(hTemp_ratio_AuAu["PiC92_2"], "PiC92_2_AuAu", refData(15,1,5));
      book(RatioAuAu["PtoPiC10"], 15,1,1);
      book(RatioAuAu["PtoPiC20"], 15,1,2);
      book(RatioAuAu["PtoPiC40"], 15,1,3);
      book(RatioAuAu["PtoPiC60"], 15,1,4);
      book(RatioAuAu["PtoPiC92"], 15,1,5);
      //ratio proton/pion (dAU)
      book(hTemp_ratio_dAu["PC20_2"], "PC20_2_dAu", refData(16,1,1));
      book(hTemp_ratio_dAu["PC100_2"], "PC100_2_dAu", refData(16,1,2));
      book(hTemp_ratio_dAu["PC40_2"], "PC40_2_dAu", refData(16,1,3));
      book(hTemp_ratio_dAu["PC60_2"], "PC60_2_dAu", refData(16,1,4));
      book(hTemp_ratio_dAu["PC88_2"], "PC88_2_dAu", refData(16,1,5));
      book(hTemp_ratio_dAu["PiC20_2"], "PiC20_2_dAu", refData(16,1,1));
      book(hTemp_ratio_dAu["PiC100_2"], "PiC100_2_dAu", refData(16,1,2));
      book(hTemp_ratio_dAu["PiC40_2"], "PiC40_2_dAu", refData(16,1,3));
      book(hTemp_ratio_dAu["PiC60_2"], "PiC60_2_dAu", refData(16,1,4));
      book(hTemp_ratio_dAu["PiC88_2"], "PiC88_2_dAu", refData(16,1,5));
      book(RatiodAu["PtoPiC20"], 16,1,1);
      book(RatiodAu["PtoPiC100"], 16,1,2);
      book(RatiodAu["PtoPiC40"], 16,1,3);
      book(RatiodAu["PtoPiC60"], 16,1,4);
      book(RatiodAu["PtoPiC88"], 16,1,5);
      //rcp kaon

      //rcp pion

      //rcp proton

      //raa kaon

      //raa pion

      //raa proton

      //rda kaon

      //rda pion

      //rda proton

      //rpc_AuAu/dAu kaon

      //rpc_AuAu/dAu pion

      //rpc_AuAu/dAu proton





    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles chargedP = applyProjection<FinalState>(event,"fs").particles();

      if (collSys == AuAu200) {

        const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
        const double cent = centProj();

        if (cent < 0. || cent > 92) vetoEvent;

        if (cent > 0. && cent < 10.) {
          sow["sow_AUAU10"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC10"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC10_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC10_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiplusC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiplusC10"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["PiplusC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["PiplusC10"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC10"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC10_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC10_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiminusC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiminusC10"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["PiminusC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["PiminusC10"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hTemp_ratio_AuAu["KC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KC10"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["KplusC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KplusC10"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["KplusC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["KplusC10"]->fill(partPt, pt_weight);
              }
            }

            if (p.pid() == -321) {
              if (getDeltaPt(*hTemp_ratio_AuAu["KC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KC10"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["KminusC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KminusC10"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["KminusC10"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["KminusC10"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC10_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC10_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC10"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC10"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hAuAu_Yields["PC10"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hAuAu_Yields["PC10"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC10_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC10_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_AuAu["P_barC10"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["P_barC10"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hAuAu_Yields["P_barC10"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hAuAu_Yields["P_barC10"]->fill(partPt, pt_weight);
                      }
                    }
            }
          }

        }

        if (cent > 10. && cent < 20.) {
          sow["sow_AUAU20"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC20_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC20_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiplusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["PiplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["PiplusC20"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC20_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC20_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiminusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["PiminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["PiminusC20"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hTemp_ratio_AuAu["KC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["KplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KplusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["KplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["KplusC20"]->fill(partPt, pt_weight);
              }
            }

            if (p.pid() == -321) {
              if (getDeltaPt(*hTemp_ratio_AuAu["KC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["KminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KminusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["KminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["KminusC20"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC20_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC20_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC20"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hAuAu_Yields["PC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hAuAu_Yields["PC20"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC20_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC20_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_AuAu["P_barC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["P_barC20"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hAuAu_Yields["P_barC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hAuAu_Yields["P_barC20"]->fill(partPt, pt_weight);
                      }
                    }
            }
          }


        }

        if (cent > 20. && cent < 40.) {
          sow["sow_AUAU40"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC40_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC40_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiplusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["PiplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["PiplusC40"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC40_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC40_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiminusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["PiminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["PiminusC40"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hTemp_ratio_AuAu["KC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["KplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KplusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["KplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["KplusC40"]->fill(partPt, pt_weight);
              }
            }

            if (p.pid() == -321) {
              if (getDeltaPt(*hTemp_ratio_AuAu["KC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["KminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KminusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["KminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["KminusC40"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC40_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC40_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC40"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hAuAu_Yields["PC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hAuAu_Yields["PC40"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC40_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC40_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_AuAu["P_barC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["P_barC40"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hAuAu_Yields["P_barC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hAuAu_Yields["P_barC40"]->fill(partPt, pt_weight);
                      }
                    }
            }
          }

        }

        if (cent > 40. && cent < 60.) {
          sow["sow_AUAU60"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC60_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC60_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiplusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["PiplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["PiplusC60"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC60_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC60_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiminusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["PiminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["PiminusC60"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hTemp_ratio_AuAu["KC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["KplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KplusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["KplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["KplusC60"]->fill(partPt, pt_weight);
              }
            }

            if (p.pid() == -321) {
              if (getDeltaPt(*hTemp_ratio_AuAu["KC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["KminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KminusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["KminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["KminusC60"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC60_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC60_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC60"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hAuAu_Yields["PC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hAuAu_Yields["PC60"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC60_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC60_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_AuAu["P_barC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["P_barC60"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hAuAu_Yields["P_barC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hAuAu_Yields["P_barC60"]->fill(partPt, pt_weight);
                      }
                    }
            }
          }

        }

        if (cent > 60. && cent < 92.) {
          sow["sow_AUAU92"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC92"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC92_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC92_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiplusC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiplusC92"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["PiplusC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["PiplusC92"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC92"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiC92_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiC92_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["PiminusC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["PiminusC92"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["PiminusC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["PiminusC92"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hTemp_ratio_AuAu["KC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KC92"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["KplusC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KplusC92"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["KplusC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["KplusC92"]->fill(partPt, pt_weight);
              }
            }

            if (p.pid() == -321) {
              if (getDeltaPt(*hTemp_ratio_AuAu["KC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KC92"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_AuAu["KminusC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_AuAu["KminusC92"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hAuAu_Yields["KminusC92"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hAuAu_Yields["KminusC92"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC92_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC92_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC92"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC92"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hAuAu_Yields["PC92"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hAuAu_Yields["PC92"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_AuAu["PC92_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["PC92_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_AuAu["P_barC92"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_AuAu["P_barC92"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hAuAu_Yields["P_barC92"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hAuAu_Yields["P_barC92"]->fill(partPt, pt_weight);
                      }
                    }
            }
          }

        }

      }

      if (collSys == dAu200) {

        const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
        const double cent = centProj();

        if (cent < 0. || cent > 100.) vetoEvent;

        if (cent > 0. && cent < 20.) {
          sow["sow_dAU20"]->fill();
          sow["sow_dAU100"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hTemp_ratio_dAu["PiC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC20_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC20_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiplusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiplusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hTemp_ratio_dAu["PiC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC20_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC20_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiminusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiminusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hTemp_ratio_dAu["KC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KplusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KplusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KplusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC100"]->fill(partPt, pt_weight);
              }
            }

            if (p.pid() == -321) {
              if (getDeltaPt(*hTemp_ratio_dAu["KC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KminusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KminusC20"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC20"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KminusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_dAu["PC20_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC20_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC20"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["PC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC20"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["PC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC100"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_dAu["PC20_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC20_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["P_barC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["P_barC20"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["P_barC20"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC20"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["P_barC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["P_barC100"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["P_barC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC40"]->fill(partPt, pt_weight);
                      }
                    }
            }
          }

        }

        if (cent > 20. && cent < 40.) {
          sow["sow_dAU40"]->fill();
          sow["sow_dAU100"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hTemp_ratio_dAu["PiC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC40_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC40_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiplusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiplusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hTemp_ratio_dAu["PiC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC40_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC40_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiminusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiminusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hTemp_ratio_dAu["KC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KplusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KplusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KplusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC100"]->fill(partPt, pt_weight);
              }
            }

            if (p.pid() == -321) {
              if (getDeltaPt(*hTemp_ratio_dAu["KC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KminusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KminusC40"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC40"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KminusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_dAu["PC40_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC40_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC40"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["PC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC40"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["PC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC100"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_dAu["PC40_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC40_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["P_barC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["P_barC40"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["P_barC40"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC40"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["P_barC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["P_barC100"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["P_barC100"], partPt, deltaPt))
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
          sow["sow_dAU100"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hTemp_ratio_dAu["PiC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC60_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC60_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiplusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiplusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hTemp_ratio_dAu["PiC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC60_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC60_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiminusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiminusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hTemp_ratio_dAu["KC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KplusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KplusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KplusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC100"]->fill(partPt, pt_weight);
              }
            }

            if (p.pid() == -321) {
              if (getDeltaPt(*hTemp_ratio_dAu["KC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KminusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KminusC60"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC60"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KminusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_dAu["PC60_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC60_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC60"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["PC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC60"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["PC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC100"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_dAu["PC60_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC60_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["P_barC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["P_barC60"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["P_barC60"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC60"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["P_barC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["P_barC100"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["P_barC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC40"]->fill(partPt, pt_weight);
                      }
                    }
            }
          }

        }

        if (cent > 60. && cent < 88.5) {
          sow["sow_dAU88"]->fill();
          sow["sow_dAU100"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hTemp_ratio_dAu["PiC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC88_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC88_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiplusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiplusC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiplusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiplusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hTemp_ratio_dAu["PiC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC88_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC88_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiminusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiminusC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiminusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiminusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hTemp_ratio_dAu["KC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KplusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KplusC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KplusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KplusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC100"]->fill(partPt, pt_weight);
              }
            }

            if (p.pid() == -321) {
              if (getDeltaPt(*hTemp_ratio_dAu["KC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KminusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KminusC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KminusC88"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC88"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KminusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_dAu["PC88_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC88_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC88"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC88"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["PC88"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC88"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["PC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC100"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_dAu["PC88_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC88_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["P_barC88"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["P_barC88"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["P_barC88"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC88"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["P_barC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["P_barC100"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["P_barC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC40"]->fill(partPt, pt_weight);
                      }
                    }
            }
          }

        }

        if (cent > 88.5) {
          sow["sow_dAU100"]->fill();
          for (const Particle& p :chargedP) {
            double partPt = p.pT() / GeV;
            double pt_weight = 1. / (2.*M_PI);
            double deltaPt = 0.;

            if (p.pid() == 211) {
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiplusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiplusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == -211) {
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiC100_2"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiC100_2"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["PiminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["PiminusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["PiminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["PiminusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 321) {
              if (getDeltaPt(*hTemp_ratio_dAu["KC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KplusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KplusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KplusC100"]->fill(partPt, pt_weight);
              }
            }

            if (p.pid() == -321) {
              if (getDeltaPt(*hTemp_ratio_dAu["KC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hTemp_ratio_dAu["KminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hTemp_ratio_dAu["KminusC100"]->fill(partPt, pt_weight);
              }
              if (getDeltaPt(*hdAu_Yields["KminusC100"], partPt, deltaPt))
              {
                pt_weight /= deltaPt;
                hdAu_Yields["KminusC100"]->fill(partPt, pt_weight);
              }
            }
            if (p.pid() == 2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["PC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["PC100"]->fill(partPt, pt_weight);
                      }
                    }
            }
            if (p.pid() == -2212) {
              if (!(p.hasAncestor(3122) || p.hasAncestor(3122) ||   //Lambda+/-
                    p.hasAncestor(3212) || p.hasAncestor(3222) ||   //Sigma0, Sigma+
                    p.hasAncestor(3322) || p.hasAncestor(3312) ||   //Xi0, Xi-
                    p.hasAncestor(3334) )) {                        //Omega-
                      if (getDeltaPt(*hTemp_ratio_dAu["PC100_2"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["PC100_2"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hTemp_ratio_dAu["P_barC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hTemp_ratio_dAu["P_barC100"]->fill(partPt, pt_weight);
                      }
                      if (getDeltaPt(*hdAu_Yields["P_barC100"], partPt, deltaPt))
                      {
                        pt_weight /= deltaPt;
                        hdAu_Yields["P_barC40"]->fill(partPt, pt_weight);
                      }
                    }
            }
          }

        }
      }
    }



    /// Normalise histograms etc., after the run
    void finalize() {

      //****Scale Histos****

      bool dAu200_available = false;
      bool AuAu200_available = false;

      for (auto element : hdAu_Yields)
      {
        string name = element.second->name();
        if (name.find("AuAu") != std::string::npos)
        {
          if (element.second->numEntries() > 0) AuAu200_available = true;
          else
          {
            AuAu200_available = false;
            break;
          }
        }
        else if (name.find("dAu") != std::string::npos)
        {
          if (element.second->numEntries() > 0) dAu200_available = true;
          else
          {
            dAu200_available = false;
            break;
          }
        }
      }

      for (auto element : hAuAu_Yields)
      {
        string name = element.second->name();
        if (name.find("AuAu") != std::string::npos)
        {
          if (element.second->numEntries() > 0) AuAu200_available = true;
          else
          {
            AuAu200_available = false;
            break;
          }
        }
        else if (name.find("dAu") != std::string::npos)
        {
          if (element.second->numEntries() > 0) dAu200_available = true;
          else
          {
            dAu200_available = false;
            break;
          }
        }
      }

    //if (!(dAu200_available && AuAu200_available)) return;

      //Yields

      hAuAu_Yields["KminusC10"]->scaleW(1./sow["sow_AUAU10"]->sumW());
      hAuAu_Yields["KminusC20"]->scaleW(1./sow["sow_AUAU20"]->sumW());
      hAuAu_Yields["KminusC40"]->scaleW(1./sow["sow_AUAU40"]->sumW());
      hAuAu_Yields["KminusC60"]->scaleW(1./sow["sow_AUAU60"]->sumW());
      hAuAu_Yields["KminusC92"]->scaleW(1./sow["sow_AUAU92"]->sumW());
      hAuAu_Yields["KplusC10"]->scaleW(1./sow["sow_AUAU10"]->sumW());
      hAuAu_Yields["KplusC20"]->scaleW(1./sow["sow_AUAU20"]->sumW());
      hAuAu_Yields["KplusC40"]->scaleW(1./sow["sow_AUAU40"]->sumW());
      hAuAu_Yields["KplusC60"]->scaleW(1./sow["sow_AUAU60"]->sumW());
      hAuAu_Yields["KplusC92"]->scaleW(1./sow["sow_AUAU92"]->sumW());

      hAuAu_Yields["PiminusC10"]->scaleW(1./sow["sow_AUAU10"]->sumW());
      hAuAu_Yields["PiminusC20"]->scaleW(1./sow["sow_AUAU20"]->sumW());
      hAuAu_Yields["PiminusC40"]->scaleW(1./sow["sow_AUAU40"]->sumW());
      hAuAu_Yields["PiminusC60"]->scaleW(1./sow["sow_AUAU60"]->sumW());
      hAuAu_Yields["PiminusC92"]->scaleW(1./sow["sow_AUAU92"]->sumW());
      hAuAu_Yields["PiplusC10"]->scaleW(1./sow["sow_AUAU10"]->sumW());
      hAuAu_Yields["PiplusC20"]->scaleW(1./sow["sow_AUAU20"]->sumW());
      hAuAu_Yields["PiplusC40"]->scaleW(1./sow["sow_AUAU40"]->sumW());
      hAuAu_Yields["PiplusC60"]->scaleW(1./sow["sow_AUAU60"]->sumW());
      hAuAu_Yields["PiplusC92"]->scaleW(1./sow["sow_AUAU92"]->sumW());

      hAuAu_Yields["P_barC10"]->scaleW(1./sow["sow_AUAU10"]->sumW());
      hAuAu_Yields["P_barC20"]->scaleW(1./sow["sow_AUAU20"]->sumW());
      hAuAu_Yields["P_barC40"]->scaleW(1./sow["sow_AUAU40"]->sumW());
      hAuAu_Yields["P_barC60"]->scaleW(1./sow["sow_AUAU60"]->sumW());
      hAuAu_Yields["P_barC92"]->scaleW(1./sow["sow_AUAU92"]->sumW());
      hAuAu_Yields["PC10"]->scaleW(1./sow["sow_AUAU10"]->sumW());
      hAuAu_Yields["PC20"]->scaleW(1./sow["sow_AUAU20"]->sumW());
      hAuAu_Yields["PC40"]->scaleW(1./sow["sow_AUAU40"]->sumW());
      hAuAu_Yields["PC60"]->scaleW(1./sow["sow_AUAU60"]->sumW());
      hAuAu_Yields["PC92"]->scaleW(1./sow["sow_AUAU92"]->sumW());

      hdAu_Yields["KminusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["KminusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_Yields["KminusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_Yields["KminusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_Yields["KminusC100"]->scaleW(1./sow["sow_dAU100"]->sumW());
      hdAu_Yields["KplusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["KplusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_Yields["KplusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_Yields["KplusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_Yields["KplusC100"]->scaleW(1./sow["sow_dAU100"]->sumW());

      hdAu_Yields["PiminusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["PiminusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_Yields["PiminusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_Yields["PiminusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_Yields["PiminusC100"]->scaleW(1./sow["sow_dAU100"]->sumW());
      hdAu_Yields["PiplusC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["PiplusC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_Yields["PiplusC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_Yields["PiplusC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_Yields["PiplusC100"]->scaleW(1./sow["sow_dAU100"]->sumW());

      hdAu_Yields["P_barC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["P_barC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_Yields["P_barC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_Yields["P_barC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_Yields["P_barC100"]->scaleW(1./sow["sow_dAU100"]->sumW());
      hdAu_Yields["PC20"]->scaleW(1./sow["sow_dAU20"]->sumW());
      hdAu_Yields["PC40"]->scaleW(1./sow["sow_dAU40"]->sumW());
      hdAu_Yields["PC60"]->scaleW(1./sow["sow_dAU60"]->sumW());
      hdAu_Yields["PC88"]->scaleW(1./sow["sow_dAU88"]->sumW());
      hdAu_Yields["PC100"]->scaleW(1./sow["sow_dAU100"]->sumW());


      //Ratios

  /*    hTemp_ratio_AuAu["KminusC10"]
      hTemp_ratio_AuAu["KminusC20"]
      hTemp_ratio_AuAu["KminusC40"]
      hTemp_ratio_AuAu["KminusC60"]
      hTemp_ratio_AuAu["KminusC92"]
      hTemp_ratio_AuAu["KplusC10"]
      hTemp_ratio_AuAu["KplusC20"]
      hTemp_ratio_AuAu["KplusC40"]
      hTemp_ratio_AuAu["KplusC60"]
      hTemp_ratio_AuAu["KplusC92"] */
      divide(hTemp_ratio_AuAu["KminusC10"], hTemp_ratio_AuAu["KplusC10"], RatioAuAu["KminustoKplusC10"]);
      divide(hTemp_ratio_AuAu["KminusC20"], hTemp_ratio_AuAu["KplusC20"], RatioAuAu["KminustoKplusC20"]);
      divide(hTemp_ratio_AuAu["KminusC40"], hTemp_ratio_AuAu["KplusC40"], RatioAuAu["KminustoKplusC40"]);
      divide(hTemp_ratio_AuAu["KminusC60"], hTemp_ratio_AuAu["KplusC60"], RatioAuAu["KminustoKplusC60"]);
      divide(hTemp_ratio_AuAu["KminusC92"], hTemp_ratio_AuAu["KplusC92"], RatioAuAu["KminustoKplusC92"]);

  /*    hTemp_ratio_dAu["KminusC20"]
      hTemp_ratio_dAu["KminusC100"]
      hTemp_ratio_dAu["KminusC40"]
      hTemp_ratio_dAu["KminusC60"]
      hTemp_ratio_dAu["KminusC88"]
      hTemp_ratio_dAu["KplusC20"]
      hTemp_ratio_dAu["KplusC100"]
      hTemp_ratio_dAu["KplusC40"]
      hTemp_ratio_dAu["KplusC60"]
      hTemp_ratio_dAu["KplusC88"] */
      divide(hTemp_ratio_dAu["KminusC20"], hTemp_ratio_dAu["KplusC20"], RatiodAu["KminustoKplusC20"]);
      divide(hTemp_ratio_dAu["KminusC100"], hTemp_ratio_dAu["KplusC100"], RatiodAu["KminustoKplusC100"]);
      divide(hTemp_ratio_dAu["KminusC40"], hTemp_ratio_dAu["KplusC40"], RatiodAu["KminustoKplusC40"]);
      divide(hTemp_ratio_dAu["KminusC60"], hTemp_ratio_dAu["KplusC60"], RatiodAu["KminustoKplusC60"]);
      divide(hTemp_ratio_dAu["KminusC88"], hTemp_ratio_dAu["KplusC88"], RatiodAu["KminustoKplusC88"]);

  /*  hTemp_ratio_AuAu["PiminusC10"]
      hTemp_ratio_AuAu["PiminusC20"]
      hTemp_ratio_AuAu["PiminusC40"]
      hTemp_ratio_AuAu["PiminusC60"]
      hTemp_ratio_AuAu["PiminusC92"]
      hTemp_ratio_AuAu["PiplusC10"]
      hTemp_ratio_AuAu["PiplusC20"]
      hTemp_ratio_AuAu["PiplusC40"]
      hTemp_ratio_AuAu["PiplusC60"]
      hTemp_ratio_AuAu["PiplusC92"] */
      divide(hTemp_ratio_AuAu["PiminusC10"], hTemp_ratio_AuAu["PiplusC10"], RatioAuAu["PiminustoPiplusC10"]);
      divide(hTemp_ratio_AuAu["PiminusC20"], hTemp_ratio_AuAu["PiplusC20"], RatioAuAu["PiminustoPiplusC20"]);
      divide(hTemp_ratio_AuAu["PiminusC40"], hTemp_ratio_AuAu["PiplusC40"], RatioAuAu["PiminustoPiplusC40"]);
      divide(hTemp_ratio_AuAu["PiminusC60"], hTemp_ratio_AuAu["PiplusC60"], RatioAuAu["PiminustoPiplusC60"]);
      divide(hTemp_ratio_AuAu["PiminusC92"], hTemp_ratio_AuAu["PiplusC92"], RatioAuAu["PiminustoPiplusC92"]);

  /*  hTemp_ratio_dAu["PiminusC20"]
      hTemp_ratio_dAu["PiminusC100"]
      hTemp_ratio_dAu["PiminusC40"]
      hTemp_ratio_dAu["PiminusC60"]
      hTemp_ratio_dAu["PiminusC88"]
      hTemp_ratio_dAu["PiplusC20"]
      hTemp_ratio_dAu["PiplusC100"]
      hTemp_ratio_dAu["PiplusC40"]
      hTemp_ratio_dAu["PiplusC60"]
      hTemp_ratio_dAu["PiplusC88"] */
      divide(hTemp_ratio_dAu["PiminusC20"], hTemp_ratio_dAu["PiplusC20"], RatiodAu["PiminustoKplusC20"]);
      divide(hTemp_ratio_dAu["PiminusC100"], hTemp_ratio_dAu["PiplusC100"], RatiodAu["PiminustoKplusC100"]);
      divide(hTemp_ratio_dAu["PiminusC40"], hTemp_ratio_dAu["PiplusC40"], RatiodAu["PiminustoKplusC40"]);
      divide(hTemp_ratio_dAu["PiminusC60"], hTemp_ratio_dAu["PiplusC60"], RatiodAu["PiminustoKplusC60"]);
      divide(hTemp_ratio_dAu["PiminusC88"], hTemp_ratio_dAu["PiplusC88"], RatiodAu["PiminustoKplusC88"]);

   /* hTemp_ratio_AuAu["P_barC10"]
      hTemp_ratio_AuAu["P_barC20"]
      hTemp_ratio_AuAu["P_barC40"]
      hTemp_ratio_AuAu["P_barC60"]
      hTemp_ratio_AuAu["P_barC92"]
      hTemp_ratio_AuAu["PC10"]
      hTemp_ratio_AuAu["PC20"]
      hTemp_ratio_AuAu["PC40"]
      hTemp_ratio_AuAu["PC60"]
      hTemp_ratio_AuAu["PC92"] */
      divide(hTemp_ratio_AuAu["P_barC10"], hTemp_ratio_AuAu["PC10"], RatioAuAu["P_bartoPC10"]);
      divide(hTemp_ratio_AuAu["P_barC20"], hTemp_ratio_AuAu["PC20"], RatioAuAu["P_bartoPC20"]);
      divide(hTemp_ratio_AuAu["P_barC40"], hTemp_ratio_AuAu["PC40"], RatioAuAu["P_bartoPC40"]);
      divide(hTemp_ratio_AuAu["P_barC60"], hTemp_ratio_AuAu["PC60"], RatioAuAu["P_bartoPC60"]);
      divide(hTemp_ratio_AuAu["P_barC92"], hTemp_ratio_AuAu["PC92"], RatioAuAu["P_bartoPC92"]);

   /* hTemp_ratio_dAu["P_barC20"]
      hTemp_ratio_dAu["P_barC100"]
      hTemp_ratio_dAu["P_barC40"]
      hTemp_ratio_dAu["P_barC60"]
      hTemp_ratio_dAu["P_barC88"]
      hTemp_ratio_dAu["PC20"]
      hTemp_ratio_dAu["PC100"]
      hTemp_ratio_dAu["PC40"]
      hTemp_ratio_dAu["PC60"]
      hTemp_ratio_dAu["PC88"] */
      divide(hTemp_ratio_dAu["P_barC20"], hTemp_ratio_dAu["PC20"], RatiodAu["P_bartoPC20"]);
      divide(hTemp_ratio_dAu["P_barC100"], hTemp_ratio_dAu["PC100"], RatiodAu["P_bartoPC100"]);
      divide(hTemp_ratio_dAu["P_barC40"], hTemp_ratio_dAu["PC40"], RatiodAu["P_bartoPC40"]);
      divide(hTemp_ratio_dAu["P_barC60"], hTemp_ratio_dAu["PC60"], RatiodAu["P_bartoPC60"]);
      divide(hTemp_ratio_dAu["P_barC88"], hTemp_ratio_dAu["PC88"], RatiodAu["P_bartoPC88"]);

   /*  hTemp_ratio_AuAu["KC10"]
      hTemp_ratio_AuAu["KC20"]
      hTemp_ratio_AuAu["KC40"]
      hTemp_ratio_AuAu["KC60"]
      hTemp_ratio_AuAu["KC92"]
      hTemp_ratio_AuAu["PiC10"]
      hTemp_ratio_AuAu["PiC20"]
      hTemp_ratio_AuAu["PiC40"]
      hTemp_ratio_AuAu["PiC60"]
      hTemp_ratio_AuAu["PiC92"] */
      divide(hTemp_ratio_AuAu["KC10"], hTemp_ratio_AuAu["PiC10"], RatioAuAu["KtoPiC10"]);
      divide(hTemp_ratio_AuAu["KC20"], hTemp_ratio_AuAu["PiC20"], RatioAuAu["KtoPiC20"]);
      divide(hTemp_ratio_AuAu["KC40"], hTemp_ratio_AuAu["PiC40"], RatioAuAu["KtoPiC40"]);
      divide(hTemp_ratio_AuAu["KC60"], hTemp_ratio_AuAu["PiC60"], RatioAuAu["KtoPiC60"]);
      divide(hTemp_ratio_AuAu["KC92"], hTemp_ratio_AuAu["PiC92"], RatioAuAu["KtoPiC92"]);

   /* hTemp_ratio_dAu["KC20"]
      hTemp_ratio_dAu["KC100"]
      hTemp_ratio_dAu["KC40"]
      hTemp_ratio_dAu["KC60"]
      hTemp_ratio_dAu["KC88"]
      hTemp_ratio_dAu["PiC20"]
      hTemp_ratio_dAu["PiC100"]
      hTemp_ratio_dAu["PiC40"]
      hTemp_ratio_dAu["PiC60"]
      hTemp_ratio_dAu["PiC88"] */
      divide(hTemp_ratio_dAu["KC20"], hTemp_ratio_dAu["PiC20"], RatiodAu["KtoPiC20"]);
      divide(hTemp_ratio_dAu["KC100"], hTemp_ratio_dAu["PiC100"], RatiodAu["KtoPiC100"]);
      divide(hTemp_ratio_dAu["KC40"], hTemp_ratio_dAu["PiC40"], RatiodAu["KtoPiC40"]);
      divide(hTemp_ratio_dAu["KC60"], hTemp_ratio_dAu["PiC60"], RatiodAu["KtoPiC60"]);
      divide(hTemp_ratio_dAu["KC88"], hTemp_ratio_dAu["PiC88"], RatiodAu["KtoPiC88"]);

   /* hTemp_ratio_AuAu["PC10_2"]
      hTemp_ratio_AuAu["PC20_2"]
      hTemp_ratio_AuAu["PC40_2"]
      hTemp_ratio_AuAu["PC60_2"]
      hTemp_ratio_AuAu["PC92_2"]
      hTemp_ratio_AuAu["PiC10_2"]
      hTemp_ratio_AuAu["PiC20_2"]
      hTemp_ratio_AuAu["PiC40_2"]
      hTemp_ratio_AuAu["PiC60_2"]
      hTemp_ratio_AuAu["PiC92_2"] */
      divide(hTemp_ratio_AuAu["PC10_2"], hTemp_ratio_AuAu["PiC10_2"], RatioAuAu["PtoPiC10"]);
      divide(hTemp_ratio_AuAu["PC20_2"], hTemp_ratio_AuAu["PiC20_2"], RatioAuAu["PtoPiC20"]);
      divide(hTemp_ratio_AuAu["PC40_2"], hTemp_ratio_AuAu["PiC40_2"], RatioAuAu["PtoPiC40"]);
      divide(hTemp_ratio_AuAu["PC60_2"], hTemp_ratio_AuAu["PiC60_2"], RatioAuAu["PtoPiC60"]);
      divide(hTemp_ratio_AuAu["PC92_2"], hTemp_ratio_AuAu["PiC92_2"], RatioAuAu["PtoPiC92"]);

   /* hTemp_ratio_dAu["PC20_2"]
      hTemp_ratio_dAu["PC100_2"]
      hTemp_ratio_dAu["PC40_2"]
      hTemp_ratio_dAu["PC60_2"]
      hTemp_ratio_dAu["PC88_2"]
      hTemp_ratio_dAu["PiC20_2"]
      hTemp_ratio_dAu["PiC100_2"]
      hTemp_ratio_dAu["PiC40_2"]
      hTemp_ratio_dAu["PiC60_2"]
      hTemp_ratio_dAu["PiC88_2"] */
      divide(hTemp_ratio_dAu["PC20_2"], hTemp_ratio_dAu["PiC20_2"], RatiodAu["PtoPiC20"]);
      divide(hTemp_ratio_dAu["PC100_2"], hTemp_ratio_dAu["PiC100_2"], RatiodAu["PtoPiC100"]);
      divide(hTemp_ratio_dAu["PC40_2"], hTemp_ratio_dAu["PiC40_2"], RatiodAu["PtoPiC40"]);
      divide(hTemp_ratio_dAu["PC60_2"], hTemp_ratio_dAu["PiC60_2"], RatiodAu["PtoPiC60"]);
      divide(hTemp_ratio_dAu["PC88_2"], hTemp_ratio_dAu["PiC88_2"], RatiodAu["PtoPiC88"]);

    }


    map<string, Histo1DPtr> hdAu_Yields;
    map<string, Histo1DPtr> hAuAu_Yields;

    map<string, Histo1DPtr> hRcp;
    map<string, Histo1DPtr> hRAA;
    map<string, Histo1DPtr> hRdA;
    map<string, Histo1DPtr> hRpc_AuAu_dAu;


    map<string, Histo1DPtr> hTemp_ratio_dAu;
    map<string, Histo1DPtr> hTemp_ratio_AuAu;
    map<string, Scatter2DPtr> RatiodAu;
    map<string, Scatter2DPtr> RatioAuAu;

    map<string, CounterPtr> sow;

    string beamOpt;
    enum CollisionSystem { AuAu200, dAu200 };
    CollisionSystem collSys;


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2013_I1227971);

}
