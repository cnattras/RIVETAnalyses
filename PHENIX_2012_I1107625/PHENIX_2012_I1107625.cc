// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "../Centralities/RHICCentrality.hh" //external header for Centrality calculation
#include "Rivet/Projections/UnstableParticles.hh"
#include <math.h>
#define _USE_MATH_DEFINES

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2012_I1107625 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2012_I1107625);

//create binShift function
    void binShift(YODA::Histo1D& histogram) {
        std::vector<YODA::HistoBin1D> binlist = histogram.bins();
        int n = 0;
        for (YODA::HistoBin1D bins : binlist) {
            double p_high = bins.xMax();
            double p_low = bins.xMin();
            //Now calculate f_corr
            if (bins.xMin() == binlist[0].xMin()) { //Check if we are working with first bin
                float b = 1 / (p_high - p_low) * log(binlist[0].height()/binlist[1].height());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
                n += 1;
            } else if (bins.xMin() == binlist.back().xMin()){ //Check if we are working with last bin
                float b = 1 / (p_high - p_low) * log(binlist[binlist.size()-2].height() / binlist.back().height());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
            } else { //Check if we are working with any middle bin
                float b = 1 / (p_high - p_low) * log(binlist[n-1].height() / binlist[n+1].height());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
                n += 1;
            }
        }
    }
    
    void init() {
      //std::initializer_list<int> pdgIds = {111};  // Pion 0

      //const PrimaryParticles fs(pdgIds, Cuts::abseta < 0.35 && Cuts::abscharge == 0);
      //declare(fs, "fs");
      const UnstableParticles ufs(Cuts::abseta < 0.35 && Cuts::pT > 0.8*GeV && Cuts::pid == 111);
      declare(ufs, "ufs");

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      //Yields_________________
      book(hPion0Pt["ptyields39c10"], 1, 1, 1);
      book(hPion0Pt["ptyields39c20"], 1, 1, 2);
      book(hPion0Pt["ptyields39c40"], 2, 1, 1);
      book(hPion0Pt["ptyields39c60"], 1, 1, 3);
      book(hPion0Pt["ptyields39c86"], 1, 1, 4);
      book(hPion0Pt["ptyields39call"], 2, 1, 2);

      book(hPion0Pt["ptyields62c10"], 3, 1, 1);
      book(hPion0Pt["ptyields62c20"], 3, 1, 2);
      book(hPion0Pt["ptyields62c40"], 3, 1, 3);
      book(hPion0Pt["ptyields62c60"], 3, 1, 4);
      book(hPion0Pt["ptyields62c86"], 4, 1, 1);
      book(hPion0Pt["ptyields62call"], 3, 1, 5);

      book(sow["sow_AuAu39c10"],"_sow_AuAu39c10");
      book(sow["sow_AuAu39c20"],"_sow_AuAu39c20");
      book(sow["sow_AuAu39c40"],"_sow_AuAu39c40");
      book(sow["sow_AuAu39c60"],"_sow_AuAu39c60");
      book(sow["sow_AuAu39c86"],"_sow_AuAu39c86");
      book(sow["sow_AuAu39call"],"_sow_AuAu39call");

      book(sow["sow_AuAu62c10"],"_sow_AuAu62c10");
      book(sow["sow_AuAu62c20"],"_sow_AuAu62c20");
      book(sow["sow_AuAu62c40"],"_sow_AuAu62c40");
      book(sow["sow_AuAu62c60"],"_sow_AuAu62c60");
      book(sow["sow_AuAu62c86"],"_sow_AuAu62c86");
      book(sow["sow_AuAu62call"],"_sow_AuAu62call");

      //RAA _______________________________

      string refnameRaa1 = mkAxisCode(5,1,1);
      const Scatter2D& refdataRaa1 =refData(refnameRaa1);
      book(hPion0Pt["c10Pt_AuAu39"], refnameRaa1 + "_AuAu39", refdataRaa1);
      book(hPion0Pt["c10Pt39_pp"], refnameRaa1 + "_pp39", refdataRaa1);
      book(hRaa["Raa_c010_AuAu39"], refnameRaa1);

      string refnameRaa2 = mkAxisCode(5,1,2);
      const Scatter2D& refdataRaa2 =refData(refnameRaa2);
      book(hPion0Pt["c20Pt_AuAu39"], refnameRaa2 + "_AuAu39", refdataRaa2);
      book(hPion0Pt["c20Pt39_pp"], refnameRaa2 + "_pp39", refdataRaa2);
      book(hRaa["Raa_c1020_AuAu39"], refnameRaa2);

      string refnameRaa3 = mkAxisCode(6,1,1);
      const Scatter2D& refdataRaa3 =refData(refnameRaa3);
      book(hPion0Pt["c40Pt_AuAu39"], refnameRaa3 + "_AuAu39", refdataRaa3);
      book(hPion0Pt["c40Pt39_pp"], refnameRaa3 + "_pp39", refdataRaa3);
      book(hRaa["Raa_c2040_AuAu39"], refnameRaa3);

      string refnameRaa4 = mkAxisCode(5,1,3);
      const Scatter2D& refdataRaa4 =refData(refnameRaa4);
      book(hPion0Pt["c60Pt_AuAu39"], refnameRaa4 + "_AuAu39", refdataRaa4);
      book(hPion0Pt["c60Pt39_pp"], refnameRaa4 + "_pp39", refdataRaa4);
      book(hRaa["Raa_c4060_AuAu39"], refnameRaa4);

      string refnameRaa5 = mkAxisCode(5,1,4);
      const Scatter2D& refdataRaa5 =refData(refnameRaa5);
      book(hPion0Pt["c86Pt_AuAu39"], refnameRaa5 + "_AuAu39", refdataRaa5);
      book(hPion0Pt["c86Pt39_pp"], refnameRaa5 + "_pp39", refdataRaa5);
      book(hRaa["Raa_c6086_AuAu39"], refnameRaa5);

      string refnameRaa6 = mkAxisCode(6,1,2);
      const Scatter2D& refdataRaa6 =refData(refnameRaa6);
      book(hPion0Pt["callPt_AuAu39"], refnameRaa6 + "_AuAu39", refdataRaa6);
      book(hPion0Pt["callPt39_pp"], refnameRaa6 + "_pp39", refdataRaa6);
      book(hRaa["Raa_minbias_AuAu39"], refnameRaa6);

      string refnameRaa7 = mkAxisCode(7,1,1);
      const Scatter2D& refdataRaa7 =refData(refnameRaa7);
      book(hPion0Pt["c10Pt_AuAu62"], refnameRaa7 + "_AuAu62", refdataRaa7);
      book(hPion0Pt["c10Pt62_pp"], refnameRaa7 + "_pp62", refdataRaa7);
      book(hRaa["Raa_c010_AuAu62"], refnameRaa7);

      string refnameRaa8 = mkAxisCode(7,1,2);
      const Scatter2D& refdataRaa8 =refData(refnameRaa8);
      book(hPion0Pt["c20Pt_AuAu62"], refnameRaa8 + "_AuAu62", refdataRaa8);
      book(hPion0Pt["c20Pt62_pp"], refnameRaa8 + "_pp62", refdataRaa8);
      book(hRaa["Raa_c1020_AuAu62"], refnameRaa8);

      string refnameRaa9 = mkAxisCode(7,1,3);
      const Scatter2D& refdataRaa9 =refData(refnameRaa9);
      book(hPion0Pt["c40Pt_AuAu62"], refnameRaa9 + "_AuAu62", refdataRaa9);
      book(hPion0Pt["c40Pt62_pp"], refnameRaa9 + "_pp62", refdataRaa9);
      book(hRaa["Raa_c2040_AuAu62"], refnameRaa9);

      string refnameRaa10 = mkAxisCode(7,1,4);
      const Scatter2D& refdataRaa10 =refData(refnameRaa10);
      book(hPion0Pt["c60Pt_AuAu62"], refnameRaa10 + "_AuAu62", refdataRaa10);
      book(hPion0Pt["c60Pt62_pp"], refnameRaa10 + "_pp62", refdataRaa10);
      book(hRaa["Raa_c4060_AuAu62"], refnameRaa10);

      string refnameRaa11 = mkAxisCode(8,1,1);
      const Scatter2D& refdataRaa11 =refData(refnameRaa11);
      book(hPion0Pt["c86Pt_AuAu62"], refnameRaa11 + "_AuAu62", refdataRaa11);
      book(hPion0Pt["c86Pt62_pp"], refnameRaa11 + "_pp62", refdataRaa11);
      book(hRaa["Raa_c6086_AuAu62"], refnameRaa11);

      string refnameRaa12 = mkAxisCode(7,1,5);
      const Scatter2D& refdataRaa12 =refData(refnameRaa12);
      book(hPion0Pt["callPt_AuAu62"], refnameRaa12 + "_AuAu62", refdataRaa12);
      book(hPion0Pt["callPt62_pp"], refnameRaa12 + "_pp62", refdataRaa12);
      book(hRaa["Raa_minbias_AuAu62"], refnameRaa12);

      book(sow["sow_pp39"],"_sow_pp39");
      book(sow["sow_pp62"],"_sow_pp62");


//Centrality vs RAA

      string refnameRaa13 = mkAxisCode(9,1,1);
      const Scatter2D& refdataRaa13 =refData(refnameRaa13);
      book(hRaaNpart["Raa_pt46_AuAu39"], refnameRaa13, refdataRaa13);

      string refnameRaa14 = mkAxisCode(9,1,2);
      const Scatter2D& refdataRaa14 =refData(refnameRaa14);
      book(hRaaNpart["Raa_pt610_AuAu39"], refnameRaa14, refdataRaa14);

      string refnameRaa15 = mkAxisCode(9,1,3);
      const Scatter2D& refdataRaa15 =refData(refnameRaa15);
      book(hRaaNpart["Raa_pt46_AuAu62"], refnameRaa15, refdataRaa15);

      string refnameRaa16 = mkAxisCode(9,1,4);
      const Scatter2D& refdataRaa16 =refData(refnameRaa16);
      book(hRaaNpart["Raa_pt610_AuAu62"], refnameRaa16, refdataRaa16);

      centBins.insert(pair<string, int>("Raa_c010_AuAu39",0));
      centBins.insert(pair<string, int>("Raa_c1020_AuAu39",1));
      centBins.insert(pair<string, int>("Raa_c2040_AuAu39",2));
      centBins.insert(pair<string, int>("Raa_c4060_AuAu39",3));
      centBins.insert(pair<string, int>("Raa_c6086_AuAu39",4));
      centBins.insert(pair<string, int>("Raa_c010_AuAu62",0));
      centBins.insert(pair<string, int>("Raa_c1020_AuAu62",1));
      centBins.insert(pair<string, int>("Raa_c2040_AuAu62",2));
      centBins.insert(pair<string, int>("Raa_c4060_AuAu62",3));
      centBins.insert(pair<string, int>("Raa_c6086_AuAu62",4));

    }


    void analyze(const Event& event) {

      const ParticlePair& beam = beams();
      int NN = 0;

      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
          NN = 197.;
          if (fuzzyEquals(sqrtS()/GeV, 39*NN, 1E-3)) collSys = AuAu39;
          if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 1E-3)) collSys = AuAu62;
      }
      if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
      {
          if (fuzzyEquals(sqrtS()/GeV, 39, 1E-3)) collSys = pp39;
          if (fuzzyEquals(sqrtS()/GeV, 62.4, 1E-3)) collSys = pp62;
      }


            Particles neutralParticles = applyProjection<UnstableParticles>(event,"ufs").particles();

            if(collSys==pp39)
            {
                sow["sow_pp39"]->fill();
                for(Particle p : neutralParticles)
                {
                    hPion0Pt["c10Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c20Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c40Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c60Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c86Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["callPt39_pp"]->fill(p.pT()/GeV);

                }
                return;
            }

            if(collSys==pp62)
            {
                sow["sow_pp62"]->fill();
                for(Particle p : neutralParticles)
                {
                    hPion0Pt["c10Pt62_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c20Pt62_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c40Pt62_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c60Pt62_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c86Pt62_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["callPt62_pp"]->fill(p.pT()/GeV);

                }
                return;
            }

            const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
            const double c = cent();

            if ((c < 0.) || (c > 86.)) vetoEvent;

            if (collSys==AuAu39)
            {
                sow["sow_AuAu39call"]->fill();

                if((c >= 0.) && (c < 10.))
                {
                    sow["sow_AuAu39c10"]->fill();
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;
                        double pt_weight = 1./(partPt*2.*M_PI);
                        hPion0Pt["ptyields39c10"]->fill(partPt, pt_weight);
                        hPion0Pt["c10Pt_AuAu39"]->fill(p.pT()/GeV);
                        hPion0Pt["callPt_AuAu39"]->fill(p.pT()/GeV);
                        hPion0Pt["ptyields39call"]->fill(p.pT()/GeV, pt_weight);
                    }
                }
                else if((c >= 10.) && (c < 20.))
                {
                    sow["sow_AuAu39c20"]->fill();
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;
                        double pt_weight = 1./(partPt*2.*M_PI);
                        hPion0Pt["ptyields39c20"]->fill(partPt, pt_weight);
                        hPion0Pt["c20Pt_AuAu39"]->fill(p.pT()/GeV);
                        hPion0Pt["callPt_AuAu39"]->fill(p.pT()/GeV);
                        hPion0Pt["ptyields39call"]->fill(p.pT()/GeV, pt_weight);
                    }
                }
                else if((c >= 20.) && (c < 40.))
                {
                    sow["sow_AuAu39c40"]->fill();
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;
                        double pt_weight = 1./(partPt*2.*M_PI);
                        hPion0Pt["ptyields39c40"]->fill(partPt, pt_weight);
                        hPion0Pt["c40Pt_AuAu39"]->fill(p.pT()/GeV);
                        hPion0Pt["callPt_AuAu39"]->fill(p.pT()/GeV);
                        hPion0Pt["ptyields39call"]->fill(p.pT()/GeV, pt_weight);
                    }
                }
                else if((c >= 40.) && (c < 60.))
                {
                    sow["sow_AuAu39c60"]->fill();
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;
                        double pt_weight = 1./(partPt*2.*M_PI);
                        hPion0Pt["ptyields39c60"]->fill(partPt, pt_weight);
                        hPion0Pt["c60Pt_AuAu39"]->fill(p.pT()/GeV);
                        hPion0Pt["callPt_AuAu39"]->fill(p.pT()/GeV);
                        hPion0Pt["ptyields39call"]->fill(p.pT()/GeV, pt_weight);
                    }
                }
                else if((c >= 60.) && (c < 86.))
                {
                    sow["sow_AuAu39c86"]->fill();
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;
                        double pt_weight = 1./(partPt*2.*M_PI);
                        hPion0Pt["ptyields39c86"]->fill(partPt, pt_weight);
                        hPion0Pt["c86Pt_AuAu39"]->fill(p.pT()/GeV);
                        hPion0Pt["callPt_AuAu39"]->fill(p.pT()/GeV);
                        hPion0Pt["ptyields39call"]->fill(p.pT()/GeV, pt_weight);
                    }
                }
                return;
            }


            if (collSys==AuAu62)
            {
                sow["sow_AuAu62call"]->fill();

                if((c >= 0.) && (c < 10.))
                {
                    sow["sow_AuAu62c10"]->fill();
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;
                        double pt_weight = 1./(partPt*2.*M_PI);
                        hPion0Pt["ptyields62c10"]->fill(partPt, pt_weight);
                        hPion0Pt["c10Pt_AuAu62"]->fill(p.pT()/GeV);
                        hPion0Pt["callPt_AuAu62"]->fill(p.pT()/GeV);
                        hPion0Pt["ptyields62call"]->fill(p.pT()/GeV, pt_weight);
                    }
                }
                else if((c >= 10.) && (c < 20.))
                {
                    sow["sow_AuAu62c20"]->fill();
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;
                        double pt_weight = 1./(partPt*2.*M_PI);
                        hPion0Pt["ptyields62c20"]->fill(partPt, pt_weight);
                        hPion0Pt["c20Pt_AuAu62"]->fill(p.pT()/GeV);
                        hPion0Pt["callPt_AuAu62"]->fill(p.pT()/GeV);
                        hPion0Pt["ptyields62call"]->fill(p.pT()/GeV, pt_weight);
                    }
                }
                else if((c >= 20.) && (c < 40.))
                {
                    sow["sow_AuAu62c40"]->fill();
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;
                        double pt_weight = 1./(partPt*2.*M_PI);
                        hPion0Pt["ptyields62c40"]->fill(partPt, pt_weight);
                        hPion0Pt["c40Pt_AuAu62"]->fill(p.pT()/GeV);
                        hPion0Pt["callPt_AuAu62"]->fill(p.pT()/GeV);
                        hPion0Pt["ptyields62call"]->fill(p.pT()/GeV, pt_weight);
                    }
                }
                else if((c >= 40.) && (c < 60.))
                {
                    sow["sow_AuAu62c60"]->fill();
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;
                        double pt_weight = 1./(partPt*2.*M_PI);
                        hPion0Pt["ptyields62c60"]->fill(partPt, pt_weight);
                        hPion0Pt["c60Pt_AuAu62"]->fill(p.pT()/GeV);
                        hPion0Pt["callPt_AuAu62"]->fill(p.pT()/GeV);
                        hPion0Pt["ptyields62call"]->fill(p.pT()/GeV, pt_weight);
                    }
                }
                else if((c >= 60.) && (c < 86.))
                {
                    sow["sow_AuAu62c86"]->fill();
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;
                        double pt_weight = 1./(partPt*2.*M_PI);
                        hPion0Pt["ptyields62c86"]->fill(partPt, pt_weight);
                        hPion0Pt["c86Pt_AuAu62"]->fill(p.pT()/GeV);
                        hPion0Pt["callPt_AuAu62"]->fill(p.pT()/GeV);
                        hPion0Pt["ptyields62call"]->fill(p.pT()/GeV, pt_weight);
                    }
                }
                return;
            }

    }

    void finalize() {

        binShift(*hPion0Pt["ptyields39c10"]);
        binShift(*hPion0Pt["ptyields39c20"]);
        binShift(*hPion0Pt["ptyields39c40"]);
        binShift(*hPion0Pt["ptyields39c60"]);
        binShift(*hPion0Pt["ptyields39c86"]);
        binShift(*hPion0Pt["ptyields39call"]);

        binShift(*hPion0Pt["ptyields62c10"]);
        binShift(*hPion0Pt["ptyields62c20"]);
        binShift(*hPion0Pt["ptyields62c40"]);
        binShift(*hPion0Pt["ptyields62c60"]);
        binShift(*hPion0Pt["ptyields62c86"]);
        binShift(*hPion0Pt["ptyields62call"]);

        binShift(*hPion0Pt["c10Pt_AuAu39"]);
        binShift(*hPion0Pt["c20Pt_AuAu39"]);
        binShift(*hPion0Pt["c40Pt_AuAu39"]);
        binShift(*hPion0Pt["c60Pt_AuAu39"]);
        binShift(*hPion0Pt["c86Pt_AuAu39"]);
        binShift(*hPion0Pt["callPt_AuAu39"]);

        binShift(*hPion0Pt["c10Pt_AuAu62"]);
        binShift(*hPion0Pt["c20Pt_AuAu62"]);
        binShift(*hPion0Pt["c40Pt_AuAu62"]);
        binShift(*hPion0Pt["c60Pt_AuAu62"]);
        binShift(*hPion0Pt["c86Pt_AuAu62"]);
        binShift(*hPion0Pt["callPt_AuAu62"]);

        binShift(*hPion0Pt["c10Pt39_pp"]);
        binShift(*hPion0Pt["c20Pt39_pp"]);
        binShift(*hPion0Pt["c40Pt39_pp"]);
        binShift(*hPion0Pt["c60Pt39_pp"]);
        binShift(*hPion0Pt["c86Pt39_pp"]);
        binShift(*hPion0Pt["callPt39_pp"]);

        binShift(*hPion0Pt["c10Pt62_pp"]);
        binShift(*hPion0Pt["c20Pt62_pp"]);
        binShift(*hPion0Pt["c40Pt62_pp"]);
        binShift(*hPion0Pt["c60Pt62_pp"]);
        binShift(*hPion0Pt["c86Pt62_pp"]);
        binShift(*hPion0Pt["callPt62_pp"]);

 
/*
      bool AuAu39_available = false;
      bool AuAu62_available = false;
      bool pp39_available = false;
      bool pp62_available = false;

      for(auto element : sow)
      {
          string name = element.second->name();
          if(name.find("AuAu39") != std::string::npos)
          {
              if(element.second->sumW() > 0) AuAu39_available = true;
          }
          else if(name.find("AuAu62") != std::string::npos)
          {
              if(element.second->sumW() > 0) AuAu62_available = true;
          }
          else if(name.find("pp39") != std::string::npos)
          {
              if(element.second->sumW() > 0) pp39_available = true;
          }
          else if(name.find("pp62") != std::string::npos)
          {
              if(element.second->sumW() > 0) pp62_available = true;
          }
      }

      if(!(AuAu62_available && AuAu39_available && pp39_available && pp62_available)) return;
*/
//Yields_________________
      hPion0Pt["ptyields39c10"]->scaleW(1./sow["sow_AuAu39c10"]->sumW());
      hPion0Pt["ptyields39c20"]->scaleW(1./sow["sow_AuAu39c20"]->sumW());
      hPion0Pt["ptyields39c40"]->scaleW(1./sow["sow_AuAu39c40"]->sumW());
      hPion0Pt["ptyields39c60"]->scaleW(1./sow["sow_AuAu39c60"]->sumW());
      hPion0Pt["ptyields39c86"]->scaleW(1./sow["sow_AuAu39c86"]->sumW());
      hPion0Pt["ptyields39call"]->scaleW(1./sow["sow_AuAu39call"]->sumW());

      hPion0Pt["ptyields62c10"]->scaleW(1./sow["sow_AuAu62c10"]->sumW());
      hPion0Pt["ptyields62c20"]->scaleW(1./sow["sow_AuAu62c20"]->sumW());
      hPion0Pt["ptyields62c40"]->scaleW(1./sow["sow_AuAu62c40"]->sumW());
      hPion0Pt["ptyields62c60"]->scaleW(1./sow["sow_AuAu62c60"]->sumW());
      hPion0Pt["ptyields62c86"]->scaleW(1./sow["sow_AuAu62c86"]->sumW());
      hPion0Pt["ptyields62call"]->scaleW(1./sow["sow_AuAu62call"]->sumW());


//RAA _______________________________

      hPion0Pt["c10Pt_AuAu39"]->scaleW(1./sow["sow_AuAu39c10"]->sumW());
      hPion0Pt["c20Pt_AuAu39"]->scaleW(1./sow["sow_AuAu39c20"]->sumW());
      hPion0Pt["c40Pt_AuAu39"]->scaleW(1./sow["sow_AuAu39c40"]->sumW());
      hPion0Pt["c60Pt_AuAu39"]->scaleW(1./sow["sow_AuAu39c60"]->sumW());
      hPion0Pt["c86Pt_AuAu39"]->scaleW(1./sow["sow_AuAu39c86"]->sumW());
      hPion0Pt["callPt_AuAu39"]->scaleW(1./sow["sow_AuAu39call"]->sumW());

      hPion0Pt["c10Pt_AuAu62"]->scaleW(1./sow["sow_AuAu62c10"]->sumW());
      hPion0Pt["c20Pt_AuAu62"]->scaleW(1./sow["sow_AuAu62c20"]->sumW());
      hPion0Pt["c40Pt_AuAu62"]->scaleW(1./sow["sow_AuAu62c40"]->sumW());
      hPion0Pt["c60Pt_AuAu62"]->scaleW(1./sow["sow_AuAu62c60"]->sumW());
      hPion0Pt["c86Pt_AuAu62"]->scaleW(1./sow["sow_AuAu62c86"]->sumW());
      hPion0Pt["callPt_AuAu62"]->scaleW(1./sow["sow_AuAu62call"]->sumW());

      hPion0Pt["c10Pt39_pp"]->scaleW(1./sow["sow_pp39"]->sumW());
      hPion0Pt["c20Pt39_pp"]->scaleW(1./sow["sow_pp39"]->sumW());
      hPion0Pt["c40Pt39_pp"]->scaleW(1./sow["sow_pp39"]->sumW());
      hPion0Pt["c60Pt39_pp"]->scaleW(1./sow["sow_pp39"]->sumW());
      hPion0Pt["c86Pt39_pp"]->scaleW(1./sow["sow_pp39"]->sumW());
      hPion0Pt["callPt39_pp"]->scaleW(1./sow["sow_pp39"]->sumW());

      hPion0Pt["c10Pt62_pp"]->scaleW(1./sow["sow_pp62"]->sumW());
      hPion0Pt["c20Pt62_pp"]->scaleW(1./sow["sow_pp62"]->sumW());
      hPion0Pt["c40Pt62_pp"]->scaleW(1./sow["sow_pp62"]->sumW());
      hPion0Pt["c60Pt62_pp"]->scaleW(1./sow["sow_pp62"]->sumW());
      hPion0Pt["c86Pt62_pp"]->scaleW(1./sow["sow_pp62"]->sumW());
      hPion0Pt["callPt62_pp"]->scaleW(1./sow["sow_pp62"]->sumW());

      divide(hPion0Pt["c10Pt_AuAu39"],hPion0Pt["c10Pt39_pp"],hRaa["Raa_c010_AuAu39"]);
      divide(hPion0Pt["c20Pt_AuAu39"],hPion0Pt["c20Pt39_pp"],hRaa["Raa_c1020_AuAu39"]);
      divide(hPion0Pt["c40Pt_AuAu39"],hPion0Pt["c40Pt39_pp"],hRaa["Raa_c2040_AuAu39"]);
      divide(hPion0Pt["c60Pt_AuAu39"],hPion0Pt["c60Pt39_pp"],hRaa["Raa_c4060_AuAu39"]);
      divide(hPion0Pt["c86Pt_AuAu39"],hPion0Pt["c86Pt39_pp"],hRaa["Raa_c6086_AuAu39"]);
      divide(hPion0Pt["callPt_AuAu39"],hPion0Pt["callPt39_pp"],hRaa["Raa_minbias_AuAu39"]);

      divide(hPion0Pt["c10Pt_AuAu62"],hPion0Pt["c10Pt62_pp"],hRaa["Raa_c010_AuAu62"]);
      divide(hPion0Pt["c20Pt_AuAu62"],hPion0Pt["c20Pt62_pp"],hRaa["Raa_c1020_AuAu62"]);
      divide(hPion0Pt["c40Pt_AuAu62"],hPion0Pt["c40Pt62_pp"],hRaa["Raa_c2040_AuAu62"]);
      divide(hPion0Pt["c60Pt_AuAu62"],hPion0Pt["c60Pt62_pp"],hRaa["Raa_c4060_AuAu62"]);
      divide(hPion0Pt["c86Pt_AuAu62"],hPion0Pt["c86Pt62_pp"],hRaa["Raa_c6086_AuAu62"]);
      divide(hPion0Pt["callPt_AuAu62"],hPion0Pt["callPt62_pp"],hRaa["Raa_minbias_AuAu62"]);

      hRaa["Raa_c010_AuAu39"]->scaleY(1./777.2);
      hRaa["Raa_c1020_AuAu39"]->scaleY(1./496.7);
      hRaa["Raa_c2040_AuAu39"]->scaleY(1./253.6);
      hRaa["Raa_c4060_AuAu39"]->scaleY(1./81.81);
      hRaa["Raa_c6086_AuAu39"]->scaleY(1./13.88);
      hRaa["Raa_minbias_AuAu39"]->scaleY(1./401.6);

      hRaa["Raa_c010_AuAu62"]->scaleY(1./843.0);
      hRaa["Raa_c1020_AuAu62"]->scaleY(1./535.8);
      hRaa["Raa_c2040_AuAu62"]->scaleY(1./270.5);
      hRaa["Raa_c4060_AuAu62"]->scaleY(1./85.71);
      hRaa["Raa_c6086_AuAu62"]->scaleY(1./14.29);
      hRaa["Raa_minbias_AuAu62"]->scaleY(1./426.9);

//Centrality vs RAA_______________________________

      for(auto element : hRaa)
      {
          string name = element.first;
          if(name.find("Raa_minbias_AuAu") != std::string::npos) continue;

          YODA::Scatter2D yodaRaa = *element.second;

          double averageRaaPt46 = 0.;
          double averageRaaPt610 = 0.;
          double averageRaaPt46Error = 0.;
          double averageRaaPt610Error = 0.;
          int nbinsPt46 = 0;
          int nbinsPt610 = 0;

          for(auto &point : yodaRaa.points())
          {

              if(point.x() > 4. && point.x() < 6.)
              {
                  if(!isnan(point.y()))
                  {
                      averageRaaPt46 += point.y();
                      averageRaaPt46Error += point.yErrPlus()*point.yErrPlus();
                      nbinsPt46++;
                  }

              }
              if(point.x() > 6.)
              {
                  if(!isnan(point.y()))
                  {
                      averageRaaPt610 += point.y();
                      averageRaaPt610Error += point.yErrPlus()*point.yErrPlus();
                      nbinsPt610++;
                  }
              }
          }


          if(name.find("39") != std::string::npos)
          {
              hRaaNpart["Raa_pt46_AuAu39"]->point(centBins[name]).setY(averageRaaPt46/nbinsPt46);
              hRaaNpart["Raa_pt46_AuAu39"]->point(centBins[name]).setYErr(sqrt(averageRaaPt46Error)/nbinsPt46);
              hRaaNpart["Raa_pt610_AuAu39"]->point(centBins[name]).setY(averageRaaPt610/nbinsPt610);
              hRaaNpart["Raa_pt610_AuAu39"]->point(centBins[name]).setYErr(sqrt(averageRaaPt610Error)/nbinsPt610);
          }
          else if(name.find("62") != std::string::npos)
          {
              hRaaNpart["Raa_pt46_AuAu62"]->point(centBins[name]).setY(averageRaaPt46/nbinsPt46);
              hRaaNpart["Raa_pt46_AuAu62"]->point(centBins[name]).setYErr(sqrt(averageRaaPt46Error)/nbinsPt46);
              hRaaNpart["Raa_pt610_AuAu62"]->point(centBins[name]).setY(averageRaaPt610/nbinsPt610);
              hRaaNpart["Raa_pt610_AuAu62"]->point(centBins[name]).setYErr(sqrt(averageRaaPt610Error)/nbinsPt610);
          }

      }


    }


    map<string, Histo1DPtr> hPion0Pt;
    map<string, Scatter2DPtr> hRaa;
    map<string, CounterPtr> sow;
    map<string, Scatter2DPtr> hRaaNpart;
    map<string, int> centBins;
    string beamOpt;
    enum CollisionSystem {pp39, pp62, AuAu39, AuAu62};
    CollisionSystem collSys;

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2012_I1107625);

}
