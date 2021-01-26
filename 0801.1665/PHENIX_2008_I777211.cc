// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/Cuts.hh"
#include "../Centralities/RHICCentrality.hh" //external header for Centrality calculation
#include <math.h>
#define _USE_MATH_DEFINES


namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2008_I777211 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2008_I777211);

    void init() {

      const UnstableParticles pi0(Cuts::absrap < 0.35 && Cuts::pT > 1*GeV && Cuts::abspid == 111 );
      declare(pi0, "pi0");

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      string refnameRaa = mkAxisCode(1,1,1);
      const Scatter2D& refdataRaa =refData(refnameRaa);

      book(hPion0Pt["Pion0Pt_AuAu"], refnameRaa + "_AuAu", refdataRaa);
      book(hPion0Pt["Pion0Pt_pp"], refnameRaa + "_pp", refdataRaa);
      book(hRaa, refnameRaa);

      book(sow["sow_AuAu"],"sow_AuAu");
      book(sow["sow_pp"],"sow_pp");

    }


    void analyze(const Event& event) {

      const ParticlePair& beam = beams();

      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu200;
      else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = pp;

      Particles neutralParticles = applyProjection<UnstableParticles>(event,"pi0").particles();

      if(collSys==pp)
      {
          sow["sow_pp"]->fill();
          for(Particle p : neutralParticles)
          {
              hPion0Pt["Pion0Pt_pp"]->fill(p.pT()/GeV);
          }
          return;
      }

      const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
      const double c = cent();

      if (c > 5.) vetoEvent;
      sow["sow_AuAu"]->fill();

      for(const Particle& p : neutralParticles)
      {
          hPion0Pt["Pion0Pt_AuAu"]->fill(p.pT()/GeV);
      }


    }

    void finalize() {

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

    }

    map<string, Histo1DPtr> hPion0Pt;
    Scatter2DPtr hRaa;
    map<string, CounterPtr> sow;
    enum CollisionSystem {pp, AuAu200};
    CollisionSystem collSys;


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2008_I777211);

}
