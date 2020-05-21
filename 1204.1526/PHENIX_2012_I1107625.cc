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
      CentralityBins = {10., 20., 40., 60., 86.};
      std::initializer_list<int> pdgIds = {111};  // Pion 0

      const PrimaryParticles fs(pdgIds, Cuts::abseta < 0.35 && Cuts::abscharge == 0);
      declare(fs, "fs");

      beamOpt = getOption<string>("beam","NONE");


      if(beamOpt=="PP") collSys = pp;
        else if(beamOpt=="AUAU39") collSys = AuAu39;
        else if(beamOpt=="AUAU62") collSys = AuAu62;
        else if(beamOpt=="AUAU200") collSys = AuAu200;


      if(!(collSys == pp)) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      //Yields_________________have to book individually bc of inconsistent data tables
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

      for(int i = 0, N = CentralityBins.size();i < N; ++i) {

            book(sow[CentralityBins[i]],"sow"+to_string(i+1));
      }


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

      //Dont know how to change this for this analyze
      string refnameRaa = mkAxisCode(1,1,1);
            const Scatter2D& refdataRaa =refData(refnameRaa);
      book(hPion0Pt["Pion0Pt_AuAu"], refnameRaa + "_AuAu", refdataRaa);
      book(hPion0Pt["Pion0Pt_pp"], refnameRaa + "_pp", refdataRaa);
      book(hRaa, refnameRaa);

      book(sow["sow_pp"],"sow_pp");
      book(sow["sow_AuAu39"],"sow_AuAu39");
      book(sow["sow_AuAu62"],"sow_AuAu62");


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
                    hPion0Pt["Pion0Pt_pp"]->fill(p.pT()/GeV);
                }
                return;
            }

            if (collSys==AuAu39)
            {

            const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
            const double c = cent();

              if ((c < 0.) || (c > 86.)) vetoEvent;

              sow["sow_AuAu39"]->fill();

              for(const Particle& p : neutralParticles)
              {
                  double partPt = p.pT()/GeV;
                  double pt_weight = 1./(partPt*2.*M_PI);
                  hPion0Pt["Pion0Pt_AuAu"]->fill(p.pT()/GeV);

                  if ((c >= 0.) && (c < 10.))
                  {

                  }
                  if ((c >= 10.) && (c < 20.))
                  {

                  }
                  if ((c >= 20.) && (c < 40.))
                  {

                  }
                  if ((c >= 40.) && (c < 60.))
                  {

                  }
                  if ((c >= 60.) && (c < 86.))
                  {

                  }

                }
              }

              if (collSys==AuAu62)
              {

              }


    }

    void finalize() {

        for (int i = 0, N = CentralityBins.size();i < N; ++i)
        {
            ptyields39c[CentralityBins[i]]->scaleW(1./sow[CentralityBins[i]]->sumW());
            ptyields62c[CentralityBins[i]]->scaleW(1./sow[CentralityBins[i]]->sumW());
            }
            ptyields39call->scaleW(1./sow[all->sumW());
            ptyields62call->scaleW(1./sowall->sumW());

//RAA _______________________________
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

      hPion0Pt["Pion0Pt_AuAu"]->scaleW(1./sow["sow_AuAu"]->sumW());
      hPion0Pt["Pion0Pt_pp"]->scaleW(1./sow["sow_pp"]->sumW());

      divide(hPion0Pt["Pion0Pt_AuAu"],hPion0Pt["Pion0Pt_pp"],hRaa);
      hRaa->scaleY(1./1051.3);
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
