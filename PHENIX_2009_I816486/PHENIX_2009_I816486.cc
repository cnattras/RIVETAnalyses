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
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2009_I816486 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2009_I816486);

    void init() {
      std::initializer_list<int> pdgIds = {111};  // Pion 0

      const PrimaryParticles fs(pdgIds, Cuts::abseta < 0.35 && Cuts::abscharge == 0);
      declare(fs, "fs");


      const ParticlePair& beam = beams();
      beamOpt = getOption<string>("beam","NONE");

    
      if (beamOpt == "NONE"){
      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu200;
      else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = pp;
      }

      if(beamOpt=="PP") collSys = pp;
      else if(beamOpt=="AUAU200") collSys = AuAu200;

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");



      //v_2_________________NEEDS TO BE CORRECTED

      book(hPion0Pt["ptv2c0005"], 1, 1, 1);
      book(hPion0Pt["ptv2c0510"], 1, 1, 2);

      for(int i = 0, N = CentralityBins.size();i < N-2; ++i)
      {
        book(hPion0Pt["ptv2c" + std::to_string(CentralityBins[i+2])], 1, 1, i+3);
      }

      book(hPion0Pt["ptv2c0020"], 2, 1, 1);
      book(hPion0Pt["ptv2c2040"], 2, 1, 2);
      book(hPion0Pt["ptv2c4060"], 2, 1, 3);
      book(hPion0Pt["ptv2c0092"], 2, 1, 4);


      //RAA _______________________________
      for(int i = 0, N = CentralityBins.size();i < N-2; ++i)
      {
        string refnameRaa=mkAxisCode(3,1,i+1);
        const Scatter2D& refdataRaa =refData(refnameRaa);
        book(hPion0Pt["c" + std::to_string(CentralityBins[i+2]) + "pt_AuAu200"],"_" + refnameRaa + "_AuAu200", refdataRaa);
        book(hPion0Pt["c" + std::to_string(CentralityBins[i+2]) + "pt_pp"], "_" + refnameRaa + "_pp", refdataRaa);
        book(hRaa["Raa" + std::to_string(CentralityBins[i+2])], refnameRaa);

        book(sow["sow_c" + std::to_string(CentralityBins[i+2])], "_sow_c" + std::to_string(CentralityBins[i+2]));
      }

        book(sow["sow_pp"],"_sow_pp");



      //dphi vs RAA_______________ NEEDS TO BE CORRECTED
      int h;

      for(int j = 0, N = PtBins.size();j < N; ++j)
      {
        for(int i = 0, M = CentralityBins.size();i < M-3; ++i)
        {
          h=i+(j*5);
          string refnameRaa = mkAxisCode(4,1,h+1);
          const Scatter2D& refdataRaa =refData(refnameRaa);
          book(hRaadphi["Raa_c"+ std::to_string(CentralityBins[i+3]) + "_pt" + std::to_string(PtBins[j]) + std::to_string(PtBins[j+1]) + "_AuAu200"], refnameRaa, refdataRaa);
        }
      }
    }

    void analyze(const Event& event) {
      Particles neutralParticles = applyProjection<PrimaryParticles>(event,"fs").particles();

            if(collSys==pp)
            {
              sow["sow_pp"]->fill();
              for(Particle p : neutralParticles)
              {
                for(int i = 0, N = CentralityBins.size();i < N-2; ++i)
                {
                  hPion0Pt["c" + std::to_string(CentralityBins[i+2]) + "pt_pp"]->fill(p.pT()/GeV);
                }

              }
              return;
            }

            const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
            const double c = cent();

            if ((c < 0.) || (c > 92.)) vetoEvent;

            if (collSys==AuAu200)
            {
                if((c >= 0.) && (c < 5.))
                {
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;

                        hPion0Pt["c10pt_AuAu200"]->fill(partPt);
                    }
                }
                else if((c >= 5.) && (c < 10.))
                {
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;

                        hPion0Pt["c10pt_AuAu200"]->fill(partPt);
                    }
                }
                else if((c >= 10.) && (c < 20.))
                {
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;

                        hPion0Pt["c20pt_AuAu200"]->fill(partPt);
                    }
                }
                else if((c >= 20.) && (c < 30.))
                {
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;

                        hPion0Pt["c30pt_AuAu200"]->fill(partPt);
                    }
                }
                else if((c >= 30.) && (c < 40.))
                {
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;

                        hPion0Pt["c40pt_AuAu200"]->fill(partPt);
                    }
                }
                else if((c >= 40.) && (c < 50.))
                {
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;

                        hPion0Pt["c50pt_AuAu200"]->fill(partPt);
                    }
                }
                else if((c >= 50.) && (c < 60.))
                {
                    for(const Particle& p : neutralParticles)
                    {
                        double partPt = p.pT()/GeV;

                        hPion0Pt["c60pt_AuAu200"]->fill(partPt);
                    }
                }

                return;
            }


    }

    void finalize() {


      //v_2_________________NEEDS TO BE IMPLEMENTED




      //RAA _______________________________

      for(int i = 0, N = CentralityBins.size();i < N-2; ++i)
      {
        hPion0Pt["c" + std::to_string(CentralityBins[i+2]) + "pt_AuAu200"]->scaleW(1./sow["sow_c" + std::to_string(CentralityBins[i+2])]->sumW());
        divide(hPion0Pt["c" + std::to_string(CentralityBins[i+2]) + "pt_AuAu200"],hPion0Pt["c" + std::to_string(CentralityBins[i+2]) + "pt_pp"],hRaa["Raa" + std::to_string(CentralityBins[i+2])]);
      }

      //need to be fixed
      hRaa["Raa10"]->scaleY(1./777.2);
      hRaa["Raa20"]->scaleY(1./777.2);
      hRaa["Raa30"]->scaleY(1./777.2);
      hRaa["Raa40"]->scaleY(1./777.2);
      hRaa["Raa50"]->scaleY(1./777.2);
      hRaa["Raa60"]->scaleY(1./777.2);


      //dphi vs RAA_______________________________NEEDS TO BE IMPLEMENTED



    }


    map<string, Histo1DPtr> hPion0Pt;
    map<string, Scatter2DPtr> hRaa;
    map<string, CounterPtr> sow;
    map<string, Scatter2DPtr> hRaadphi;
    string beamOpt;
    enum CollisionSystem {pp, AuAu200};
    CollisionSystem collSys;
    vector<int> CentralityBins {0, 5, 10, 20, 30, 40, 50, 60};
    vector<double> PtBins {1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

  };



  DECLARE_RIVET_PLUGIN(PHENIX_2009_I816486);

}
