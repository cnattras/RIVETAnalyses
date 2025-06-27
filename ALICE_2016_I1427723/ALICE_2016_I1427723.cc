// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <fstream>
#include <iostream>
#include <math.h>
#define _USE_MATH_DEFINES

using namespace std;

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2016_I1427723 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2016_I1427723);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // the basic final-state projection: all final-state particles within the given eta acceptance
        
      // Declare centrality projection
      const FinalState fs(Cuts::abseta < 0.6);
      declare(fs, "FS");
      const ChargedFinalState cfs(Cuts::abseta < 0.6);
      declare(cfs, "CFS");

      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");
      
      book(_h["ET"], 1, 1, 2);
      book(_h["EToverNpart"], 1, 1, 3);
      
      string refnameETChpart = mkAxisCode(3,1,1);
      const Estimate1D& refdataETChpart = refData(refnameETChpart);
      book(_h["mean_ET_ratio"],refnameETChpart + "_ET", refdataETChpart.xEdges());
      book(_h["mean_Chpart_ratio"],refnameETChpart + "_Chpart", refdataETChpart.xEdges());
      book(MeanETChpart,refnameETChpart);

      // The centrality bins upper bin edges.
      centralityBins = { 2.5,5.0,7.5,10., 12.5,15.0,17.5, 20.,22.5,25.0,27.5,30., 32.5,35.0,37.5,40., 45.0,50., 55.0,60., 65.0,70., 75.0,80.};
      npart = {396.1,369.3,342.4,316.5, 292.5,269.9,248.8,229.3, 210.8,193.5,177,161.9, 147.6,134.5,121.9,110.2, 94.11,75.3, 59.24,45.58, 34.33,25.21, 17.96,12.58};
      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const FinalState& fs = apply<FinalState>(event, "FS");
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

      // Prepare centrality projection and value
      const CentralityProjection& centrProj = apply<CentralityProjection>(event, "V0M");
      double centr = centrProj();

      // Veto event for too large centralities since those are not used in the analysis at all
      if ((centr < 0.) || (centr >= 80.)){
        vetoEvent;
      }

      //calculate ET in the event
      
      int centbin = -1;
      if(centr < 40.) //Just a machinery to take the right npart accordingly to centrality
      {
          centbin = ((int)centr/2.5);
      }
      else
      {
          centbin = 16 + ((int)(centr-40.)/5.);
      }
      
      double mynpart = npart[centbin];
            
      double ET = 0.;
      for(const Particle& p : fs.particles())
      {
          ET += p.Et()/GeV;
      }
      double mult = cfs.size();
      
      _h["ET"]->fill(centr,ET);
      _h["EToverNpart"]->fill(centr,ET/(mynpart/2.));
      _h["mean_ET_ratio"]->fill(centr,ET);
      _h["mean_Chpart_ratio"]->fill(centr,mult);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

        _h["ET"]->scale(1, 1./1.2);
        _h["EToverNpart"]->scale(1, 1./1.2);
        _h["mean_ET_ratio"]->scale(1, 1./1.2);
        _h["mean_Chpart_ratio"]->scale(1, 1./1.2);
        divide(_h["mean_ET_ratio"],_h["mean_Chpart_ratio"],MeanETChpart);
    }

    //@}

    vector<double> centralityBins;
    vector<double> npart;
    /// @name Histograms
    //@{
    map<string, Profile1DPtr> _h;
    Estimate1DPtr MeanETChpart;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALICE_2016_I1427723);


}
