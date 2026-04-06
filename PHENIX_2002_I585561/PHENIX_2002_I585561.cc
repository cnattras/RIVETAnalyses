// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

#include <cmath>

namespace Rivet {


  class PHENIX_2002_I585561 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2002_I585561);


    void init() override {
      declare(UnstableParticles(), "UFS");

      // Minimum-bias spectra
      book(_hLambdaMB, 1, 1, 1);
      book(_hLambdaBarMB, 1, 1, 2);

      // Ratio vs pT
      book(_hRatioPtNum, "TMP/LambdaRatioPtNum", refData(1, 1, 3));
      book(_hRatioPtDen, "TMP/LambdaRatioPtDen", refData(1, 1, 3));
      book(_sRatioPt, 1, 1, 3);

      book(_cAllEvents, "TMP/AllEvents");
    }


    void analyze(const Event& event) override {
      _cAllEvents->fill();

      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        if (p.abspid() != PID::LAMBDA) continue;
        if (std::abs(p.rapidity()) > 0.5) continue;

        const double pt = p.pT() / GeV;
        if (!std::isfinite(pt)) continue;

        if (p.pid() == PID::LAMBDA) {
          _hLambdaMB->fill(pt);
          _hRatioPtDen->fill(pt);
        } else if (p.pid() == -PID::LAMBDA) {
          _hLambdaBarMB->fill(pt);
          _hRatioPtNum->fill(pt);
        }
      }
    }


    void finalize() override {
      const double nAll = _cAllEvents->sumW();
      if (nAll <= 0.0) return;

      scaleInvariantYield(_hLambdaMB, nAll);
      scaleInvariantYield(_hLambdaBarMB, nAll);
      scaleInvariantYield(_hRatioPtNum, nAll);
      scaleInvariantYield(_hRatioPtDen, nAll);

      binShift(_hLambdaMB);
      binShift(_hLambdaBarMB);
      binShift(_hRatioPtNum);
      binShift(_hRatioPtDen);

      fillRatio(_sRatioPt, _hRatioPtNum, _hRatioPtDen);
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
        if (std::isfinite(factor) && factor > 0.0) bin.scaleW(factor);
      }
    }


    double binHeight(Histo1DPtr hist, size_t i) const {
      const auto& bin = hist->bin(i);
      const double width = bin.xWidth();
      return (width > 0.0) ? bin.sumW() / width : 0.0;
    }


    void binShift(Histo1DPtr hist) const {
      if (!hist || hist->numBins() < 2) return;

      for (size_t i = 0; i < hist->numBins(); ++i) {
        const auto& bin = hist->bin(i);
        const double pHigh = bin.xMax();
        const double pLow = bin.xMin();
        const double width = pHigh - pLow;
        if (width <= 0.0) continue;

        const double yThis = binHeight(hist, i);
        if (yThis <= 0.0) continue;

        double yLo = yThis;
        double yHi = yThis;

        if (i == 0) {
          yHi = binHeight(hist, i + 1);
        } else if (i + 1 == hist->numBins()) {
          yLo = binHeight(hist, i - 1);
        } else {
          yLo = binHeight(hist, i - 1);
          yHi = binHeight(hist, i + 1);
        }

        if (yLo <= 0.0 || yHi <= 0.0) continue;

        const double b = std::log(yLo / yHi) / width;
        if (!std::isfinite(b)) continue;

        const double expLow = std::exp(-b * pLow);
        const double expHigh = std::exp(-b * pHigh);
        const double denom = expHigh - expLow;
        if (fuzzyEquals(denom, 0.0)) continue;

        const double fCorr =
          -b * width * std::exp(-b * 0.5 * (pHigh + pLow)) / denom;

        if (std::isfinite(fCorr) && fCorr > 0.0 && fCorr < 10.0) {
          hist->bin(i).scaleW(fCorr);
        }
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


    void fillRatio(Scatter2DPtr scatter, Histo1DPtr num, Histo1DPtr den) const {
      if (!scatter || !num || !den) return;

      scatter->reset();

      for (size_t i = 0; i < num->numBins(); ++i) {
        const auto& numBin = num->bin(i);
        const auto& denBin = den->bin(i);

        const double numerator = numBin.sumW();
        const double denominator = denBin.sumW();
        const double x = numBin.xMid();
        const double dx = 0.5 * numBin.xWidth();

        double ratio = 0.0;
        double err = 0.0;

        if (numerator > 0.0 && denominator > 0.0) {
          ratio = numerator / denominator;
          err = ratioError(numerator, numBin.sumW2(), denominator, denBin.sumW2());
        }

        if (!std::isfinite(x) || !std::isfinite(dx)) continue;
        if (!std::isfinite(ratio)) ratio = 0.0;
        if (!std::isfinite(err)) err = 0.0;

        scatter->addPoint(x, ratio, dx, err);
      }
    }


    Histo1DPtr _hLambdaMB;
    Histo1DPtr _hLambdaBarMB;

    Histo1DPtr _hRatioPtNum;
    Histo1DPtr _hRatioPtDen;
    Scatter2DPtr _sRatioPt;

    CounterPtr _cAllEvents;
  };


  RIVET_DECLARE_PLUGIN(PHENIX_2002_I585561);

}
