i/ -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "../Centralities/RHICCentrality.hh" //external header for Centrality calculation
#define _USE_MATH_DEFINES
namespace Rivet {
  /// @brief Add a short analysis description here
  class STAR_2003_I619063 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2003_I619063);
    /// @name Analysis methods
    //@{
    bool getBinCenter(YODA::Histo1D hist, double pT, double &binCenter)
    {
        if(pT > hist.xMin() && pT < hist.xMax())
        {
            binCenter = hist.bin(hist.binIndexAt(pT)).xMid();
            return true;
        }
        else return false;
    }

    /// Book histograms and initialise projections before the run
    void init() {

      beamOpt = getOption<string>("beam", "NONE");

      // Initialise and register projections

      std::initializer_list<int> pdgIds = {PID::PROTON, PID::PIPLUS, PID::KPLUS};
		  const PrimaryParticles cpp(pdgIds, Cuts::absrap < 0.5);
      declare(cpp,"cpp");
      //Ratios using AuAu at 130 GeV
      string refnameEnergyRatio05 = mkAxisCode(2,1,1);
      const Scatter2D& refdataEnergyRatio05 =refData(refnameEnergyRatio05);

      string refnameEnergyRatio2030 = mkAxisCode(2,1,2);
      const Scatter2D& refdataEnergyRatio2030 =refData(refnameEnergyRatio2030);

      string refnameEnergyRatio3040 = mkAxisCode(2,1,3);
      const Scatter2D& refdataEnergyRatio3040 =refData(refnameEnergyRatio3040);

      string refnameEnergyRatio4060 = mkAxisCode(2,1,4);
      const Scatter2D& refdataEnergyRatio4060 =refData(refnameEnergyRatio4060);

      //Ratios using pp at 200 GeV
      string refnameppRatio05 = mkAxisCode(3,1,1);
      const Scatter2D& refdatappRatio05 =refData(refnameppRatio05);

      string refnameppRatio1020 = mkAxisCode(3,1,2);
      const Scatter2D& refdatappRatio1020 =refData(refnameppRatio1020);

      string refnameppRatio2030 = mkAxisCode(3,1,3);
      const Scatter2D& refdatappRatio2030 =refData(refnameppRatio2030);

      string refnameppRatio3040 = mkAxisCode(3,1,4);
      const Scatter2D& refdatappRatio3040 =refData(refnameppRatio3040);

      string refnameppRatio4060 = mkAxisCode(3,1,5);
      const Scatter2D& refdatappRatio4060 =refData(refnameppRatio4060);

      string refnameppRatio6080 = mkAxisCode(3,1,6);
      const Scatter2D& refdatappRatio6080 =refData(refnameppRatio6080);

      //Centrality
      declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

      //pp at 200 GeV
      book(chSpectrum["chSpectrum0_5_pp"], refnameppRatio05 + "_pp", refdatappRatio05);
      book(R_AB["R_200_pp_05"], refnameppRatio05);

      book(chSpectrum["chSpectrum10_20_pp"], refnameppRatio1020 + "_pp", refdatappRatio1020);
      book(R_AB["R_200_pp_1020"], refnameppRatio1020);

      book(chSpectrum["chSpectrum20_30_pp"], refnameppRatio2030 + "_pp", refdatappRatio2030);
      book(R_AB["R_200_pp_2030"], refnameppRatio2030);

      book(chSpectrum["chSpectrum30_40_pp"], refnameppRatio3040 + "_pp", refdatappRatio3040);
      book(R_AB["R_200_pp_3040"], refnameppRatio3040);

      book(chSpectrum["chSpectrum40_60_pp"], refnameppRatio4060 + "_pp", refdatappRatio4060);
      book(R_AB["R_200_pp_4060"], refnameppRatio4060);

      book(chSpectrum["chSpectrum60_80_pp"], refnameppRatio6080 + "_pp", refdatappRatio6080);
      book(R_AB["R_200_pp_6080"], refnameppRatio6080);

      book(sow["sow_pp"], "sow_pp");

      //AuAu at 200 GeV
      book(chSpectrum["chSpectrum0_5"], 1, 1, 1);
      book(chSpectrum["chSpectrum5_10"], 1, 1, 2);
      book(chSpectrum["chSpectrum10_20"], 1, 1, 3);
      book(chSpectrum["chSpectrum20_30"], 1, 1, 4);
      book(chSpectrum["chSpectrum30_40"], 1, 1, 5);
      book(chSpectrum["chSpectrum40_60"], 1, 1, 6);
      book(chSpectrum["chSpectrum60_80"], 1, 1, 7);
      book(sow["sow0_5"], "sow0_5");
      book(sow["sow5_10"], "sow5_10");
      book(sow["sow10_20"], "sow10_20");
      book(sow["sow20_30"], "sow20_30");
      book(sow["sow30_40"], "sow30_40");
      book(sow["sow40_60"], "sow40_60");
      book(sow["sow60_80"], "sow60_80");

      string refnameCentRatio05_4060 = mkAxisCode(4,1,1);
      book(R_AB["Rcp0_5_over_40_60"], refnameCentRatio05_4060);

      string refnameCentRatio05_6080 = mkAxisCode(4,1,2);
      book(R_AB["Rcp0_5_over_60_80"], refnameCentRatio05_6080);

      book(chSpectrum["chSpectrum_C05_S200130"], refnameEnergyRatio05 + "_200", refdataEnergyRatio05);
      book(chSpectrum["chSpectrum_C2030_S200130"], refnameEnergyRatio2030 + "_200", refdataEnergyRatio2030);
      book(chSpectrum["chSpectrum_C3040_S200130"], refnameEnergyRatio3040 + "_200", refdataEnergyRatio3040);
      book(chSpectrum["chSpectrum_C4060_S200130"], refnameEnergyRatio4060 + "_200", refdataEnergyRatio4060);

      book(chSpectrum["chSpectrum_C05_S200pp"], refnameppRatio05 + "_200", refdatappRatio05);
      book(chSpectrum["chSpectrum_C1020_S200pp"], refnameppRatio1020 + "_200", refdatappRatio1020);
      book(chSpectrum["chSpectrum_C2030_S200pp"], refnameppRatio2030 + "_200", refdatappRatio2030);
      book(chSpectrum["chSpectrum_C3040_S200pp"], refnameppRatio3040 + "_200", refdatappRatio3040);
      book(chSpectrum["chSpectrum_C4060_S200pp"], refnameppRatio4060 + "_200", refdatappRatio4060);
      book(chSpectrum["chSpectrum_C6080_S200pp"], refnameppRatio6080 + "_200", refdatappRatio6080);

      //AuAu at 130 GeV
      book(sow["sow0_5_130"], "sow0_5_130");
      book(sow["sow20_30_130"], "sow20_30_130");
      book(sow["sow30_40_130"], "sow30_40_130");
      book(sow["sow40_60_130"], "sow40_60_130");


      book(chSpectrum["chSpectrum0_5_130"], refnameEnergyRatio05 + "_130", refdataEnergyRatio05);
      book(R_AB["R_200_130_05"], refnameEnergyRatio05);

      book(chSpectrum["chSpectrum20_30_130"], refnameEnergyRatio2030 + "_130", refdataEnergyRatio2030);
      book(R_AB["R_200_130_2030"], refnameEnergyRatio2030);

      book(chSpectrum["chSpectrum30_40_130"], refnameEnergyRatio3040 + "_130", refdataEnergyRatio3040);
      book(R_AB["R_200_130_3040"], refnameEnergyRatio3040);

      book(chSpectrum["chSpectrum40_60_130"], refnameEnergyRatio4060 + "_130", refdataEnergyRatio4060);
      book(R_AB["R_200_130_4060"], refnameEnergyRatio4060);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here

      const ParticlePair& beam = beams();
      int NN = 0;

      if(beamOpt == "NONE")
      {
              if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
              {
                  NN = 197.;
                  if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSys = AuAu200;
                  if (fuzzyEquals(sqrtS()/GeV, 130*NN, 1E-3)) collSys = AuAu130;
              }
              if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
              {
                  collSys = pp;
              }
      }
      else if(beamOpt == "PP200") collSys = pp;
      else if(beamOpt == "AUAU200") collSys = AuAu200;
      else if(beamOpt == "AUAU130") collSys = AuAu130;



      PrimaryParticles cpp = applyProjection<PrimaryParticles>(event,"cpp");
      Particles particles = cpp.particles();

      if(collSys == pp)
      {
          sow["sow_pp"]->fill();
          for(Particle p : particles)
          {
              chSpectrum["chSpectrum0_5_pp"]->fill(p.pT()/GeV);
              chSpectrum["chSpectrum10_20_pp"]->fill(p.pT()/GeV);
              chSpectrum["chSpectrum20_30_pp"]->fill(p.pT()/GeV);
              chSpectrum["chSpectrum30_40_pp"]->fill(p.pT()/GeV);
              chSpectrum["chSpectrum40_60_pp"]->fill(p.pT()/GeV);
              chSpectrum["chSpectrum60_80_pp"]->fill(p.pT()/GeV);
          }
          return;
      }

      //Get cenrality
      const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
      const double c = cent();

      if(c >= 80) vetoEvent;

      double binCenter = 0.;

      if(c < 5.)
      {

          if(collSys == AuAu200)
          {
              for(Particle p : particles)
              {
                  if(getBinCenter(*chSpectrum["chSpectrum0_5"], p.pT()/GeV, binCenter))
                  {
                      chSpectrum["chSpectrum0_5"]->fill(p.pT()/GeV, 1./binCenter);
                  }
                  chSpectrum["chSpectrum_C05_S200130"]->fill(p.pT()/GeV);
                  chSpectrum["chSpectrum_C05_S200pp"]->fill(p.pT()/GeV);
              }
              sow["sow0_5"]->fill();
          }
          else if(collSys == AuAu130)
          {
              for(Particle p : particles)
              {
                  chSpectrum["chSpectrum0_5_130"]->fill(p.pT()/GeV);
              }
              sow["sow0_5_130"]->fill();
          }
      }
      else if(c >= 5 && c < 10)
      {
          if(collSys == AuAu200)
          {
              for(Particle p : particles)
              {
                  if(getBinCenter(*chSpectrum["chSpectrum5_10"], p.pT()/GeV, binCenter))
                  {
                      chSpectrum["chSpectrum5_10"]->fill(p.pT()/GeV, 1./binCenter);
                  }
              }
              sow["sow5_10"]->fill();
          }
      }
      else if(c >= 10 && c < 20)
      {
          if(collSys == AuAu200)
          {
              for(Particle p : particles)
              {
                  if(getBinCenter(*chSpectrum["chSpectrum10_20"], p.pT()/GeV, binCenter))
                  {
                      chSpectrum["chSpectrum10_20"]->fill(p.pT()/GeV, 1./binCenter);
                  }
                  chSpectrum["chSpectrum_C1020_S200pp"]->fill(p.pT()/GeV);
              }
              sow["sow10_20"]->fill();
          }
      }
      else if(c >= 20 && c < 30)
      {
          if(collSys == AuAu200)
          {
              for(Particle p : particles)
              {
                  if(getBinCenter(*chSpectrum["chSpectrum20_30"], p.pT()/GeV, binCenter))
                  {
                      chSpectrum["chSpectrum20_30"]->fill(p.pT()/GeV, 1./binCenter);
                  }
                  chSpectrum["chSpectrum_C2030_S200130"]->fill(p.pT()/GeV);
                  chSpectrum["chSpectrum_C2030_S200pp"]->fill(p.pT()/GeV);
              }
              sow["sow20_30"]->fill();
          }
          else if(collSys == AuAu130)
          {
              for(Particle p : particles)
              {
                  chSpectrum["chSpectrum20_30_130"]->fill(p.pT()/GeV);
              }
              sow["sow20_30_130"]->fill();
          }
      }
      else if(c >= 30 && c < 40)
      {
          if(collSys == AuAu200)
          {
              for(Particle p : particles)
              {
                  if(getBinCenter(*chSpectrum["chSpectrum30_40"], p.pT()/GeV, binCenter))
                  {
                      chSpectrum["chSpectrum30_40"]->fill(p.pT()/GeV, 1./binCenter);
                  }
                  chSpectrum["chSpectrum_C3040_S200130"]->fill(p.pT()/GeV);
                  chSpectrum["chSpectrum_C3040_S200pp"]->fill(p.pT()/GeV);
              }
              sow["sow30_40"]->fill();
          }
          else if(collSys == AuAu130)
          {
              for(Particle p : particles)
              {
                  chSpectrum["chSpectrum30_40_130"]->fill(p.pT()/GeV);
              }
              sow["sow30_40_130"]->fill();
          }

      }
      else if(c >= 40 && c < 60)
      {
          if(collSys == AuAu200)
          {
              for(Particle p : particles)
              {
                  if(getBinCenter(*chSpectrum["chSpectrum40_60"], p.pT()/GeV, binCenter))
                  {
                      chSpectrum["chSpectrum40_60"]->fill(p.pT()/GeV, 1./binCenter);
                  }
                  chSpectrum["chSpectrum_C4060_S200130"]->fill(p.pT()/GeV);
                  chSpectrum["chSpectrum_C4060_S200pp"]->fill(p.pT()/GeV);
              }
              sow["sow40_60"]->fill();
          }
          else if(collSys == AuAu130)
          {
              for(Particle p : particles)
              {
                  chSpectrum["chSpectrum40_60_130"]->fill(p.pT()/GeV);
              }
              sow["sow40_60_130"]->fill();
          }
      }
      else if(c >= 60 && c < 80)
      {
          if(collSys == AuAu200)
          {
              for(Particle p : particles)
              {
                  if(getBinCenter(*chSpectrum["chSpectrum60_80"], p.pT()/GeV, binCenter))
                  {
                      chSpectrum["chSpectrum60_80"]->fill(p.pT()/GeV, 1./binCenter);
                  }
                  chSpectrum["chSpectrum_C6080_S200pp"]->fill(p.pT()/GeV);
              }
              sow["sow60_80"]->fill();
          }
      }



    }


    /// Normalise histograms etc., after the run
    void finalize() {
      //These lines normalize per event and pT bin width.
      //You need to also scale by 1.0/(2.0*3.14159*[eta acceptance])

      bool AuAu130_available = false;
      bool AuAu200_available = false;
      bool pp200_available = false;

      for(auto element : chSpectrum)
      {
          string name = element.second->name();
          if(name.find("200") != std::string::npos)
          {
              if(element.second->numEntries() > 0) AuAu200_available = true;
          }
          else if(name.find("130") != std::string::npos)
          {
              if(element.second->numEntries() > 0) AuAu130_available = true;
          }
          else if(name.find("pp") != std::string::npos)
          {
              if(element.second->numEntries() > 0) pp200_available = true;
          }
      }

      if(!(AuAu200_available && AuAu130_available && pp200_available)) return;

      chSpectrum["chSpectrum0_5"]->scaleW(1./(sow["sow0_5"]->sumW()*4.0*M_PI));
      chSpectrum["chSpectrum5_10"]->scaleW(1./(sow["sow5_10"]->sumW()*4.0*M_PI));
      chSpectrum["chSpectrum10_20"]->scaleW(1./(sow["sow10_20"]->sumW()*4.0*M_PI));
      chSpectrum["chSpectrum20_30"]->scaleW(1./(sow["sow20_30"]->sumW()*4.0*M_PI));
      chSpectrum["chSpectrum30_40"]->scaleW(1./(sow["sow30_40"]->sumW()*4.0*M_PI));
      chSpectrum["chSpectrum40_60"]->scaleW(1./(sow["sow40_60"]->sumW()*4.0*M_PI));
      chSpectrum["chSpectrum60_80"]->scaleW(1./(sow["sow60_80"]->sumW()*4.0*M_PI));

      chSpectrum["chSpectrum_C05_S200130"]->scaleW(1./sow["sow0_5"]->sumW());
      chSpectrum["chSpectrum_C2030_S200130"]->scaleW(1./sow["sow20_30"]->sumW());
      chSpectrum["chSpectrum_C3040_S200130"]->scaleW(1./sow["sow30_40"]->sumW());
      chSpectrum["chSpectrum_C4060_S200130"]->scaleW(1./sow["sow40_60"]->sumW());

      chSpectrum["chSpectrum_C05_S200pp"]->scaleW(1./sow["sow0_5"]->sumW());
      chSpectrum["chSpectrum_C1020_S200pp"]->scaleW(1./sow["sow10_20"]->sumW());
      chSpectrum["chSpectrum_C2030_S200pp"]->scaleW(1./sow["sow20_30"]->sumW());
      chSpectrum["chSpectrum_C3040_S200pp"]->scaleW(1./sow["sow30_40"]->sumW());
      chSpectrum["chSpectrum_C4060_S200pp"]->scaleW(1./sow["sow40_60"]->sumW());
      chSpectrum["chSpectrum_C6080_S200pp"]->scaleW(1./sow["sow60_80"]->sumW());

      divide(chSpectrum["chSpectrum0_5"], chSpectrum["chSpectrum40_60"], R_AB["Rcp0_5_over_40_60"]);
      R_AB["Rcp0_5_over_40_60"]->scaleY(93.6/1051.3);

      divide(chSpectrum["chSpectrum0_5"], chSpectrum["chSpectrum60_80"], R_AB["Rcp0_5_over_60_80"]);
      R_AB["Rcp0_5_over_60_80"]->scaleY(21.2/1051.3);

      chSpectrum["chSpectrum0_5_130"]->scaleW(1./sow["sow0_5_130"]->sumW());
      chSpectrum["chSpectrum20_30_130"]->scaleW(1./sow["sow20_30_130"]->sumW());
      chSpectrum["chSpectrum30_40_130"]->scaleW(1./sow["sow30_40_130"]->sumW());
      chSpectrum["chSpectrum40_60_130"]->scaleW(1./sow["sow40_60_130"]->sumW());

      divide(chSpectrum["chSpectrum_C05_S200130"],chSpectrum["chSpectrum0_5_130"],R_AB["R_200_130_05"]);
      divide(chSpectrum["chSpectrum_C2030_S200130"],chSpectrum["chSpectrum20_30_130"],R_AB["R_200_130_2030"]);
      divide(chSpectrum["chSpectrum_C3040_S200130"],chSpectrum["chSpectrum30_40_130"],R_AB["R_200_130_3040"]);
      divide(chSpectrum["chSpectrum_C4060_S200130"],chSpectrum["chSpectrum40_60_130"],R_AB["R_200_130_4060"]);

      chSpectrum["chSpectrum0_5_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      chSpectrum["chSpectrum10_20_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      chSpectrum["chSpectrum20_30_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      chSpectrum["chSpectrum30_40_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      chSpectrum["chSpectrum40_60_pp"]->scaleW(1./sow["sow_pp"]->sumW());
      chSpectrum["chSpectrum60_80_pp"]->scaleW(1./sow["sow_pp"]->sumW());

      divide(chSpectrum["chSpectrum_C05_S200pp"],chSpectrum["chSpectrum0_5_pp"],R_AB["R_200_pp_05"]);
      R_AB["R_200_pp_05"]->scaleY(1./1051.3);

      divide(chSpectrum["chSpectrum_C1020_S200pp"],chSpectrum["chSpectrum10_20_pp"],R_AB["R_200_pp_1020"]);
      R_AB["R_200_pp_1020"]->scaleY(1./591.3);

      divide(chSpectrum["chSpectrum_C2030_S200pp"],chSpectrum["chSpectrum20_30_pp"],R_AB["R_200_pp_2030"]);
      R_AB["R_200_pp_2030"]->scaleY(1./368.6);

      divide(chSpectrum["chSpectrum_C3040_S200pp"],chSpectrum["chSpectrum30_40_pp"],R_AB["R_200_pp_3040"]);
      R_AB["R_200_pp_3040"]->scaleY(1./220.2);

      divide(chSpectrum["chSpectrum_C4060_S200pp"],chSpectrum["chSpectrum40_60_pp"],R_AB["R_200_pp_4060"]);
      R_AB["R_200_pp_4060"]->scaleY(1./93.6);

      divide(chSpectrum["chSpectrum_C6080_S200pp"],chSpectrum["chSpectrum60_80_pp"],R_AB["R_200_pp_6080"]);
      R_AB["R_200_pp_6080"]->scaleY(1./21.2);

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> chSpectrum;
    map<string, CounterPtr> sow;
    map<string, Scatter2DPtr> R_AB;
    enum CollisionSystem {pp, AuAu130, AuAu200};
    CollisionSystem collSys;
    string beamOpt = "NONE";
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(STAR_2003_I619063);


}
