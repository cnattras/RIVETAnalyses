// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/ImpactParameterProjection.hh"
#include "RHICCentrality.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Tools/BeamConstraint.hh"
#include "Rivet/Projections/GeneratedPercentileProjection.hh"
#include "Rivet/Projections/UserCentEstimate.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#define _USE_MATH_DEFINES

namespace Rivet {

  /// @brief Centrality projection for PHENIX AuAu.
  class RHIC_2019_CentralityCalibration : public Analysis {

  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(RHIC_2019_CentralityCalibration);

    /// Book histograms and initialise projections before the run
    void init() {
      // One projection for the actual observable, and one for the
      // generated impact parameter.

      string experiment = getOption<string>("exp","STAR");
      set<string> done;
      MSG_INFO("RHIC Experiment: " << experiment);
      declare(RHICCentrality(experiment), "Centrality");
      declare(ImpactParameterProjection(), "IMP");

      // The calibration histogram:
      book(_calib, "CMULT", 100, 0.0, 2000.0);

      // If histogram was pre-loaded, the calibration is done.
      //_done = ( _calib->numEntries() > 0 );

      // The alternative histogram based on impact parameter. Note that
      // it MUST be named the same as the histogram for the experimental
      // observable with an added _IMP suffix for the Pecentile<>
      // binning to work properly.
      book(_impcalib, "CMULT_IMP", 400, 0.0, 20.0);


    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // The alternative centrality based on generated impact
      // parameter, assumes that the generator does not describe the
      // full final state, and should therefore be filled even if the
      // event is not triggered.
        _impcalib->fill(apply<SingleValueProjection>(event, "IMP")());

        _calib->fill(apply<SingleValueProjection>(event, "Centrality")());

    }

    /// Finalize
    void finalize() {

      _calib->normalize();
      _impcalib->normalize();

    }

    /// The calibration histograms.
    Histo1DPtr _calib;
    Histo1DPtr _impcalib;
    };

    // The hook for the plugin system
    RIVET_DECLARE_PLUGIN(RHIC_2019_CentralityCalibration);
}

