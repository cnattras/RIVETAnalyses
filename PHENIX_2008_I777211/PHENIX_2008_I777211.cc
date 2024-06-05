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

      beamOpt = getOption<string>("beam", "NONE");

      const UnstableParticles pi0(Cuts::absrap < 0.35 && Cuts::pT > 1*GeV && Cuts::abspid == 111 );
      declare(pi0, "pi0");

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      string refnameRaa = mkAxisCode(1,1,1);
      const Scatter2D& refdataRaa =refData(refnameRaa);

      book(hPion0Pt["Pion0Pt_AuAu"], refnameRaa + "_AuAu", refdataRaa);
      book(hPion0Pt["Pion0Pt_pp"], refnameRaa + "_pp", refdataRaa);
      book(hRaa, refnameRaa);

      book(sow["sow_AuAu"],"_sow_AuAu");
      book(sow["sow_pp"],"_sow_pp");

    }


    void analyze(const Event& event) {

      const ParticlePair& beam = beams();

      if (beamOpt == "NONE") {

        if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu200;
        else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = pp;
      }

      else if (beamOpt == "PP200") collSys = pp;
      else if (beamOpt == "AUAU200") collSys = AuAu200;

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
      binShift(*hPion0Pt["Pion0Pt_AuAu"]);
      binShift(*hPion0Pt["Pion0Pt_pp"]);
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
    string beamOpt;

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2008_I777211);

}
