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
//#include "Centrality/RHICCentrality.hh" //external header for Centrality calculation
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
      //V2CentralityBins = {0005., 0510., 0010., 1020., 2030., 3040., 4050., 5060.};
      RaaCentralityBins = {10., 20., 30., 40., 50., 60.};
      //PtBins = {1.5., 2.0., 2.5., 3.0., 3.5., 4.0., 5.0., 6.0., 7.0., 8.0., 9.0, 10.0.};

      const PrimaryParticles fs(pdgIds, Cuts::abseta < 0.35 && Cuts::abscharge == 0);
      declare(fs, "fs");

      beamOpt = getOption<string>("beam","NONE");


      if(beamOpt=="PP") collSys = pp;
      else if(beamOpt=="AUAU200") collSys = AuAu200;


    //  if(!(collSys == pp)) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");



      //v_2_________________

      //RAA _______________________________
      for(int i = 0, N = RaaCentralityBins.size();i < N; ++i)
      {
        string refnameRaa=mkAxisCode(3,1,i+1);
        const Scatter2D& refdataRaa =refData(refnameRaa);
        book(hPion0Pt[RaaCentralityBins[i] + "pt_AuAu200"], refnameRaa + "_AuAu200", refdataRaa);
        book(hPion0Pt[RaaCentralityBins[i] + "pt_pp"], refnameRaa + "_pp", refdataRaa);
        book(hRaa["Raa" + RaaCentralityBins[i]], refnameRaa);
      }

        book(sow["sow_pp"],"sow_pp");


      //dphi vs RAA


    }


    void analyze(const Event& event) {



    }

    void finalize() {


      //v_2_________________




      //RAA _______________________________



      //dphi vs RAA_______________________________



    }


    map<string, Histo1DPtr> hPion0Pt;
    map<string, Scatter2DPtr> hRaa;
    map<string, CounterPtr> sow;
    map<string, Scatter2DPtr> hRaaNpart;
    map<string, int> centBins;
    string beamOpt;
    enum CollisionSystem {pp, AuAu200};
    CollisionSystem collSys;
    vector<double> RaaCentralityBins;

  };



  DECLARE_RIVET_PLUGIN(PHENIX_2009_I816486);

}
