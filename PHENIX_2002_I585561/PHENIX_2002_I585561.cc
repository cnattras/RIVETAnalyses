// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "../Centralities/RHICCentrality.hh"

#include <cmath>

namespace Rivet {


  class PHENIX_2002_I585561 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2002_I585561);


    void init() override {
      declare(FinalState(), "FS");

      const RHICCentrality rhicCent("PHENIX");
      declare(rhicCent, "RHIC_CENT");
      declareCentrality(rhicCent,
                        "RHIC_2019_CentralityCallibration:exp=PHENIX",
                        "CMULT", "CENT");

      // Figure 2 spectra from the reference data.
      book(_hLambdaMB, 1, 1, 1);
      book(_hLambdaBarMB, 1, 1, 2);
      book(_hLambdaCentral, 1, 1, 3);
      book(_hLambdaBarCentral, 1, 1, 4);

      // Keep the temporary ratio histograms on the same binning as the reference plots.
      book(_hRatioPtNum, "TMP/LambdaRatioPtNum", refData(2, 1, 1));
      book(_hRatioPtDen, "TMP/LambdaRatioPtDen", refData(2, 1, 1));
      book(_sRatioPt, 2, 1, 1);

      book(_hRatioNpartNum, "TMP/LambdaRatioNpartNum", refData(3, 1, 1));
      book(_hRatioNpartDen, "TMP/LambdaRatioNpartDen", refData(3, 1, 1));
      book(_sRatioNpart, 3, 1, 1);

      book(_cAllEvents, "TMP/AllEvents");
      book(_cCentralEvents, "TMP/CentralEvents");
    }


    void analyze(const Event& event) override {
      const CentralityProjection& centProj = apply<CentralityProjection>(event, "CENT");
      //if (!centProj.isValueSet()) vetoEvent;

      const double cent = centProj();
      const bool isCentral = (cent >= 0.0 && cent < 5.0);
      cout<<"centrality "<<cent<<endl;

       // const RHICCentrality& rhicCent = apply<RHICCentrality>(event, "RHIC_CENT");
      //const double npart = rhicCent.npart();

      _cAllEvents->fill();
      if (isCentral) _cCentralEvents->fill();

      for (const Particle& p : apply<FinalState>(event, "FS").particles()) {
        if (p.abspid() != PID::LAMBDA) continue;
        if (std::abs(p.rapidity()) > 0.5) continue;
        cout<<"Found Lambda"<<endl;

        const double pt = p.pT() / GeV;
        if (!std::isfinite(pt)) continue;

        if (p.pid() == PID::LAMBDA) {
          _hLambdaMB->fill(pt);
          if (isCentral) _hLambdaCentral->fill(pt);

          _hRatioPtDen->fill(pt);
          //if (std::isfinite(npart)) _hRatioNpartDen->fill(npart);
        } else {
          _hLambdaBarMB->fill(pt);
          if (isCentral) _hLambdaBarCentral->fill(pt);

          _hRatioPtNum->fill(pt);
          //if (std::isfinite(npart)) _hRatioNpartNum->fill(npart);
        }
      }
    }


    void finalize() override {
      const double nAll = _cAllEvents->sumW();
      const double nCentral = _cCentralEvents->sumW();

      if (nAll > 0.0) {
        scaleInvariantYield(_hLambdaMB, nAll);
        scaleInvariantYield(_hLambdaBarMB, nAll);
      }

      if (nCentral > 0.0) {
        scaleInvariantYield(_hLambdaCentral, nCentral);
        scaleInvariantYield(_hLambdaBarCentral, nCentral);
      }

      fillRatio(_sRatioPt, _hRatioPtNum, _hRatioPtDen);
      fillRatio(_sRatioNpart, _hRatioNpartNum, _hRatioNpartDen);
    }


  private:

    static constexpr double TWOPI = 6.28318530717958647693;
    static constexpr double DY = 1.0;

    void scaleInvariantYield(Histo1DPtr hist, double nEvents) const {
      if (!hist || nEvents <= 0.0) return;

      for (size_t i = 0; i < hist->numBins(); ++i) {
        auto& bin = hist->bin(i);
        const double pt = bin.xMid();
        const double dpt = bin.xWidth();

        if (pt <= 0.0 || dpt <= 0.0) continue;

        const double factor = 1.0 / (nEvents * TWOPI * pt * dpt * DY);
        if (std::isfinite(factor) && factor > 0.0) {
          bin.scaleW(factor);
        }
      }
    }


    void fillRatio(Scatter2DPtr scatter, Histo1DPtr num, Histo1DPtr den) const {
      if (!scatter || !num || !den) return;

      scatter->reset();

      for (size_t i = 0; i < num->numBins(); ++i) {
        const auto& numBin = num->bin(i);
        const auto& denBin = den->bin(i);

        const double numerator = numBin.sumW();
        const double denominator = denBin.sumW();
        if (numerator <= 0.0 || denominator <= 0.0) continue;

        const double x = numBin.xMid();
        const double dx = 0.5 * numBin.xWidth();
        const double ratio = numerator / denominator;
        const double err = ratioError(numerator, numBin.sumW2(),
                                      denominator, denBin.sumW2());

        if (!std::isfinite(x) || !std::isfinite(dx)) continue;
        if (!std::isfinite(ratio) || !std::isfinite(err)) continue;

        scatter->addPoint(x, ratio, dx, err);
      }
    }


    double ratioError(double numerator, double numeratorErr2,
                      double denominator, double denominatorErr2) const {
      if (numerator <= 0.0 || denominator <= 0.0) return 0.0;
      if (numeratorErr2 < 0.0 || denominatorErr2 < 0.0) return 0.0;

      const double ratio = numerator / denominator;
      const double relErr2 = numeratorErr2 / sqr(numerator)
                           + denominatorErr2 / sqr(denominator);
      const double err = ratio * std::sqrt(relErr2);
      return std::isfinite(err) ? err : 0.0;
    }


    Histo1DPtr _hLambdaMB;
    Histo1DPtr _hLambdaBarMB;
    Histo1DPtr _hLambdaCentral;
    Histo1DPtr _hLambdaBarCentral;

    Histo1DPtr _hRatioPtNum;
    Histo1DPtr _hRatioPtDen;
    Histo1DPtr _hRatioNpartNum;
    Histo1DPtr _hRatioNpartDen;

    Scatter2DPtr _sRatioPt;
    Scatter2DPtr _sRatioNpart;

    CounterPtr _cAllEvents;
    CounterPtr _cCentralEvents;
  };


  RIVET_DECLARE_PLUGIN(PHENIX_2002_I585561);

}
