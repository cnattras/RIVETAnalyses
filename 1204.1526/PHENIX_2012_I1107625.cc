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
#include "Centrality/RHICCentrality.hh" //external header for Centrality calculation
#include <math.h>
#define _USE_MATH_DEFINES

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2012_I1107625 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2012_I1107625);

    void init() {
      std::initializer_list<int> pdgIds = {111};  // Pion 0

      const PrimaryParticles fs(pdgIds, Cuts::abseta < 0.35 && Cuts::abscharge == 0);
      declare(fs, "fs");

      beamOpt = getOption<string>("beam","NONE");


      if(beamOpt=="PP") collSys = pp;
        else if(beamOpt=="AUAU39") collSys = AuAu39;
        else if(beamOpt=="AUAU62") collSys = AuAu62;
        else if(beamOpt=="AUAU200") collSys = AuAu200;


      if(!(collSys == pp)) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      //Yields_________________
      book(ptyields39c10, 1, 1, 1);
      book(ptyields39c20, 1, 1, 2);
      book(ptyields39c40, 2, 1, 1);
      book(ptyields39c60, 1, 1, 3);
      book(ptyields39c86, 1, 1, 4);
      book(ptyields39call, 2, 1, 2);

      book(ptyields62c10, 3, 1, 1);
      book(ptyields62c20, 3, 1, 2);
      book(ptyields62c40, 3, 1, 3);
      book(ptyields62c60, 3, 1, 4);
      book(ptyields62c86, 4, 1, 1);
      book(ptyields62call, 3, 1, 5);

      book(sow["sow_39c10"],"sow_39c10");
      book(sow["sow_39c20"],"sow_39c20");
      book(sow["sow_39c40"],"sow_39c40");
      book(sow["sow_39c60"],"sow_39c60");
      book(sow["sow_39c86"],"sow_39c86");
      book(sow["sow_39call"],"sow_39call");

      book(sow["sow_62c10"],"sow_62c10");
      book(sow["sow_62c20"],"sow_62c20");
      book(sow["sow_62c40"],"sow_62c40");
      book(sow["sow_62c60"],"sow_62c60");
      book(sow["sow_62c86"],"sow_62c86");
      book(sow["sow_39call"],"sow_39call");

      //RAA _______________________________
      book(ptRAA39c10, 5, 1, 1);
      book(ptRAA39c20, 5, 1, 2);
      book(ptRAA39c40, 6, 1, 1);
      book(ptRAA39c60, 5, 1, 3);
      book(ptRAA39c86, 5, 1, 4);
      book(ptRAA39call, 6, 1, 2);

      book(ptRAA62c10, 7, 1, 1);
      book(ptRAA62c20, 7, 1, 2);
      book(ptRAA62c40, 7, 1, 3);
      book(ptRAA62c60, 7, 1, 4);
      book(ptRAA62c86, 8, 1, 1);
      book(ptRAA62call, 7, 1, 5);

      string refnameRaa1 = mkAxisCode(5,1,1);
            const Scatter2D& refdataRaa1 =refData(refnameRaa1);
      book(hPion0Pt["c10Pt_AuAu39"], refnameRaa1 + "_AuAu39", refdataRaa1);
      book(hPion0Pt["c10Pt39_pp"], refnameRaa1 + "_pp", refdataRaa1);
      book(hRaa1, refnameRaa1);

      string refnameRaa2 = mkAxisCode(5,1,2);
            const Scatter2D& refdataRaa2 =refData(refnameRaa2);
      book(hPion0Pt["c20Pt_AuAu39"], refnameRaa2 + "_AuAu39", refdataRaa2);
      book(hPion0Pt["c20Pt39_pp"], refnameRaa2 + "_pp", refdataRaa2);
      book(hRaa2, refnameRaa2);

      string refnameRaa3 = mkAxisCode(6,1,1);
            const Scatter2D& refdataRaa3 =refData(refnameRaa3);
      book(hPion0Pt["c40Pt_AuAu39"], refnameRaa3 + "_AuAu39", refdataRaa3);
      book(hPion0Pt["c40Pt39_pp"], refnameRaa3 + "_pp", refdataRaa3);
      book(hRaa3, refnameRaa3);

      string refnameRaa4 = mkAxisCode(5,1,3);
            const Scatter2D& refdataRaa4 =refData(refnameRaa4);
      book(hPion0Pt["c60Pt_AuAu39"], refnameRaa4 + "_AuAu39", refdataRaa4);
      book(hPion0Pt["c60Pt39_pp"], refnameRaa4 + "_pp", refdataRaa4);
      book(hRaa4, refnameRaa4);

      string refnameRaa5 = mkAxisCode(5,1,4);
            const Scatter2D& refdataRaa5 =refData(refnameRaa5);
      book(hPion0Pt["c86Pt_AuAu39"], refnameRaa5 + "_AuAu39", refdataRaa5);
      book(hPion0Pt["c86Pt39_pp"], refnameRaa5 + "_pp", refdataRaa5);
      book(hRaa5, refnameRaa5);

      string refnameRaa6 = mkAxisCode(6,1,2);
            const Scatter2D& refdataRaa6 =refData(refnameRaa6);
      book(hPion0Pt["callPt_AuAu39"], refnameRaa6 + "_AuAu39", refdataRaa6);
      book(hPion0Pt["callPt39_pp"], refnameRaa6 + "_pp", refdataRaa6);
      book(hRaa6, refnameRaa6);

      string refnameRaa7 = mkAxisCode(7,1,1);
            const Scatter2D& refdataRaa7 =refData(refnameRaa2);
      book(hPion0Pt["c10Pt_AuAu62"], refnameRaa7 + "_AuAu62", refdataRaa7);
      book(hPion0Pt["c10Pt62_pp"], refnameRaa7 + "_pp", refdataRaa7);
      book(hRaa7, refnameRaa7);

      string refnameRaa8 = mkAxisCode(7,1,2);
            const Scatter2D& refdataRaa8 =refData(refnameRaa8);
      book(hPion0Pt["c20Pt_AuAu62"], refnameRaa8 + "_AuAu62", refdataRaa8);
      book(hPion0Pt["c20Pt62_pp"], refnameRaa8 + "_pp", refdataRaa8);
      book(hRaa8, refnameRaa8);

      string refnameRaa9 = mkAxisCode(7,1,3);
            const Scatter2D& refdataRaa9 =refData(refnameRaa29);
      book(hPion0Pt["c40Pt_AuAu62"], refnameRaa9 + "_AuAu62", refdataRaa9);
      book(hPion0Pt["c40Pt62_pp"], refnameRaa9 + "_pp", refdataRaa9);
      book(hRaa9, refnameRaa9);

      string refnameRaa10 = mkAxisCode(7,1,4);
            const Scatter2D& refdataRaa10 =refData(refnameRaa10);
      book(hPion0Pt["c60Pt_AuAu62"], refnameRaa10 + "_AuAu62", refdataRaa10);
      book(hPion0Pt["c60Pt62_pp"], refnameRaa10 + "_pp", refdataRaa10);
      book(hRaa10, refnameRaa10);

      string refnameRaa11 = mkAxisCode(8,1,1);
            const Scatter2D& refdataRaa11 =refData(refnameRaa11);
      book(hPion0Pt["c86Pt_AuAu62"], refnameRaa11 + "_AuAu62", refdataRaa11);
      book(hPion0Pt["c86Pt62_pp"], refnameRaa11 + "_pp", refdataRaa11);
      book(hRaa11, refnameRaa11);

      string refnameRaa12 = mkAxisCode(7,1,5);
            const Scatter2D& refdataRaa12 =refData(refnameRaa12);
      book(hPion0Pt["callPt_AuAu62"], refnameRaa12 + "_AuAu62", refdataRaa12);
      book(hPion0Pt["callPt62_pp"], refnameRaa12 + "_pp", refdataRaa12);
      book(hRaa12, refnameRaa12);

      book(sow["sow_pp"],"sow_pp");


//Centrality vs RAA
      book(centRAA39, 9, 1, 1);
      book(centRAA62, 9, 1, 3);
      }


    void analyze(const Event& event) {


            Particles neutralParticles = applyProjection<PrimaryParticles>(event,"fs").particles();

            if(collSys==pp)
            {
                sow["sow_pp"]->fill();
                for(Particle p : neutralParticles)
                {
                    hPion0Pt["c10Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c20Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c40Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c60Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c86Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["callPt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c10Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c20Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c40Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c60Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["c86Pt39_pp"]->fill(p.pT()/GeV);
                    hPion0Pt["callPt39_pp"]->fill(p.pT()/GeV);

                }
                return;
            }
            const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
            const double c = cent();
            if (collSys==AuAu39)
            {

              if ((c < 0.) || (c > 86.)) vetoEvent;

              sow["sow_AuAu39"]->fill();

              for(const Particle& p : neutralParticles)
              {
                  double partPt = p.pT()/GeV;
                  double pt_weight = 1./(partPt*2.*M_PI);

                  sow["sow_39call"]->fill();
                  hPion0Pt["callPt_AuAu39"]->fill(p.pT()/GeV);

                  if ((c >= 0.) && (c < 10.))
                  {
                    sow["sow_39c10"]->fill();
                    ptyields39c10->fill(partPt, pt_weight);
                    hPion0Pt["c10Pt_AuAu39"]->fill(p.pT()/GeV);
                  }
                  if ((c >= 10.) && (c < 20.))
                  {
                    sow["sow_39c20"]->fill();
                    ptyields39c20->fill(partPt, pt_weight);
                    hPion0Pt["c20Pt_AuAu39"]->fill(p.pT()/GeV);
                  }
                  if ((c >= 20.) && (c < 40.))
                  {
                    sow["sow_39c40"]->fill();
                    ptyields39c40->fill(partPt, pt_weight);
                    hPion0Pt["c40Pt_AuAu39"]->fill(p.pT()/GeV);
                  }
                  if ((c >= 40.) && (c < 60.))
                  {
                    sow["sow_39c60"]->fill();
                    ptyields39c60->fill(partPt, pt_weight);
                    hPion0Pt["c60Pt_AuAu39"]->fill(p.pT()/GeV);
                  }
                  if ((c >= 60.) && (c < 86.))
                  {
                    sow["sow_39c86"]->fill();
                    ptyields39c86->fill(partPt, pt_weight);
                    hPion0Pt["c86Pt_AuAu39"]->fill(p.pT()/GeV);
                  }

                }
              }

              if (collSys==AuAu62)
              {
                if ((c < 0.) || (c > 86.)) vetoEvent;

                sow["sow_AuAu62"]->fill();

                for(const Particle& p : neutralParticles)
                {
                    double partPt = p.pT()/GeV;
                    double pt_weight = 1./(partPt*2.*M_PI);

                    sow["sow_62call"]->fill();
                    hPion0Pt["callPt_AuAu62"]->fill(p.pT()/GeV);

                    if ((c >= 0.) && (c < 10.))
                    {
                      sow["sow_62c10"]->fill();
                      ptyields62c10->fill(partPt, pt_weight);
                      hPion0Pt["c10Pt_AuAu39"]->fill(p.pT()/GeV);
                    }
                    if ((c >= 10.) && (c < 20.))
                    {
                      sow["sow_62c20"]->fill();
                      ptyields62c20->fill(partPt, pt_weight);
                      hPion0Pt["c20Pt_AuAu39"]->fill(p.pT()/GeV);
                    }
                    if ((c >= 20.) && (c < 40.))
                    {
                      sow["sow_62c40"]->fill();
                      ptyields62c40->fill(partPt, pt_weight);
                      hPion0Pt["c40Pt_AuAu39"]->fill(p.pT()/GeV);
                    }
                    if ((c >= 40.) && (c < 60.))
                    {
                      sow["sow_62c60"]->fill();
                      ptyields62c60->fill(partPt, pt_weight);
                      hPion0Pt["c60Pt_AuAu39"]->fill(p.pT()/GeV);
                    }
                    if ((c >= 60.) && (c < 86.))
                    {
                      sow["sow_62c86"]->fill();
                      ptyields62c86->fill(partPt, pt_weight);
                      hPion0Pt["c86Pt_AuAu39"]->fill(p.pT()/GeV);
                    }
                  }
                }
    }

    void finalize() {

//Yields_________________
ptyields39c10->scaleW(1./sow["sow_39c10"]->sumW());
ptyields39c20->scaleW(1./sow["sow_39c20"]->sumW());
ptyields39c40->scaleW(1./sow["sow_39c40"]->sumW());
ptyields39c60->scaleW(1./sow["sow_39c60"]->sumW());
ptyields39c86->scaleW(1./sow["sow_39c86"]->sumW());
ptyields39call->scaleW(1./sow["sow_39call"]->sumW());

ptyields62c10->scaleW(1./sow["sow_62c10"]->sumW());
ptyields62c20->scaleW(1./sow["sow_62c20"]->sumW());
ptyields62c40->scaleW(1./sow["sow_62c40"]->sumW());
ptyields62c60->scaleW(1./sow["sow_62c60"]->sumW());
ptyields62c86->scaleW(1./sow["sow_62c86"]->sumW());
ptyields39call->scaleW(1./sow["sow_39call"]->sumW())


//RAA _______________________________
      bool AuAu39_available = false;
      bool AuAu62_available = false;
      bool AuAu200_available = false;
      bool pp200_available = false;

      for(auto element : hPion0Pt)
      {
          string name = element.second->name();
          if(name.find("AuAu") != std::string::npos)
          {
              if(element.second->numEntries() > 0) AuAu200_available = true;
          }
          else if(name.find("pp") != std::string::npos)
          {
              if(element.second->numEntries() > 0) pp200_available = true;
          }
      }

      if(!(AuAu200_available && pp200_available)) return;

      hPion0Pt["c10Pt_AuAu39"]->scaleW(1./sow["sow_AuAu39"]->sumW());
      hPion0Pt["c20Pt_AuAu39"]->scaleW(1./sow["sow_AuAu39"]->sumW());
      hPion0Pt["c40Pt_AuAu39"]->scaleW(1./sow["sow_AuAu39"]->sumW());
      hPion0Pt["c60Pt_AuAu39"]->scaleW(1./sow["sow_AuAu39"]->sumW());
      hPion0Pt["c86Pt_AuAu39"]->scaleW(1./sow["sow_AuAu39"]->sumW());
      hPion0Pt["callPt_AuAu39"]->scaleW(1./sow["sow_AuAu39"]->sumW());

      hPion0Pt["c10Pt_AuAu62"]->scaleW(1./sow["sow_AuAu62"]->sumW());
      hPion0Pt["c20Pt_AuAu62"]->scaleW(1./sow["sow_AuAu62"]->sumW());
      hPion0Pt["c40Pt_AuAu62"]->scaleW(1./sow["sow_AuAu62"]->sumW());
      hPion0Pt["c60Pt_AuAu62"]->scaleW(1./sow["sow_AuAu62"]->sumW());
      hPion0Pt["c86Pt_AuAu62"]->scaleW(1./sow["sow_AuAu62"]->sumW());
      hPion0Pt["callPt_AuAu62"]->scaleW(1./sow["sow_AuAu62"]->sumW());

      hPion0Pt["c10Pt39_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      hPion0Pt["c20Pt39_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      hPion0Pt["c40Pt39_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      hPion0Pt["c60Pt39_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      hPion0Pt["c86Pt39_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      hPion0Pt["callPt39_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      hPion0Pt["c10Pt62_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      hPion0Pt["c20Pt62_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      hPion0Pt["c40Pt62_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      hPion0Pt["c60Pt62_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      hPion0Pt["c86Pt62_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      hPion0Pt["callPt62_pp"]->scaleW(1./sow["sow_pp"]->sumW());

      divide(hPion0Pt["c10Pt_AuAu39"],hPion0Pt["c10Pt39_pp"],hRaa1);
      divide(hPion0Pt["c20Pt_AuAu39"],hPion0Pt["c20Pt39_pp"],hRaa2);
      divide(hPion0Pt["c40Pt_AuAu39"],hPion0Pt["c40Pt39_pp"],hRaa3);
      divide(hPion0Pt["c60Pt_AuAu39"],hPion0Pt["c60Pt39_pp"],hRaa4);
      divide(hPion0Pt["c86Pt_AuAu39"],hPion0Pt["c86Pt39_pp"],hRaa5);
      divide(hPion0Pt["callPt_AuAu39"],hPion0Pt["callPt39_pp"],hRaa6);
      divide(hPion0Pt["c10Pt_AuAu62"],hPion0Pt["c10Pt62_pp"],hRaa7);
      divide(hPion0Pt["c20Pt_AuAu62"],hPion0Pt["c20Pt62_pp"],hRaa8);
      divide(hPion0Pt["c40Pt_AuAu62"],hPion0Pt["c40Pt62_pp"],hRaa9);
      divide(hPion0Pt["c60Pt_AuAu62"],hPion0Pt["c60Pt62_pp"],hRaa10);
      divide(hPion0Pt["c86Pt_AuAu62"],hPion0Pt["c86Pt62_pp"],hRaa11);
      divide(hPion0Pt["callPt_AuAu62"],hPion0Pt["callPt62_pp"],hRaa12);

      hRaa1->scaleY(1./1051.3);
      hRaa2->scaleY(1./1051.3);
      hRaa3->scaleY(1./1051.3);
      hRaa4->scaleY(1./1051.3);
      hRaa5->scaleY(1./1051.3);
      hRaa6->scaleY(1./1051.3);
      hRaa7->scaleY(1./1051.3);
      hRaa8->scaleY(1./1051.3);
      hRaa9->scaleY(1./1051.3);
      hRaa10->scaleY(1./1051.3);
      hRaa11->scaleY(1./1051.3);
      hRaa12->scaleY(1./1051.3);
//_______________________________







    }


    map<string, Histo1DPtr> hPion0Pt;
    Scatter2DPtr hRaa;
    map<string, CounterPtr> sow;
    string beamOpt;
    enum CollisionSystem {pp, AuAu200};
    CollisionSystem collSys;


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2012_I1107625);

}
