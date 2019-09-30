
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

  class ALICE_2013_I1222333 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2013_I1222333);


    /// name Analysis methods

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // the basic final-state projection: all final-state particles within the given eta acceptance

      const ChargedFinalState cfs(Cuts::abseta < 5 && Cuts::pT > 100*MeV);
      declare(cfs, "CFS");

      // Declare centrality projection
      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Book histograms
			    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

      // Prepare centrality projection and value
      const CentralityProjection& centrProj = apply<CentralityProjection>(event, "V0M");
      double centr = centrProj();

      // Veto event for too large centralities since those are not used in the analysis at all
      if ((centr < 0.) || (centr > 80.)){
        vetoEvent;
      }


    // loop over charged final state particles
      for(const Particle& p : cfs.particles()) {

	       // protections against mc generators decaying long-lived particles
	        if (!( p.hasAncestor(310) || p.hasAncestor(-310)  ||     // K0s
	           p.hasAncestor(130)  || p.hasAncestor(-130)  ||     // K0l
	           p.hasAncestor(3322) || p.hasAncestor(-3322) ||     // Xi0
	           p.hasAncestor(3122) || p.hasAncestor(-3122) ||     // Lambda
	           p.hasAncestor(3222) || p.hasAncestor(-3222) ||     // Sigma+/-
	           p.hasAncestor(3312) || p.hasAncestor(-3312) ||     // Xi-/+
	           p.hasAncestor(3334) || p.hasAncestor(-3334) ))    // Omega-/+
	        {
          if (centr >= 0. && centr <= 5.) {
            switch (p.pid()) {
              	case 211: // pi+
              	  {

              	    break;
              	  }
              	case -211: // pi-
              	  {

              	    break;
              	  }
              	case 2212: // proton
              	  {

              	    break;
              	  }
              	case -2212: // anti-proton
                  {

              	    break;
              	  }
              	case 321: // K+
              	  {

              	    break;
              	  }
              	case -321: // K-
              	  {

              	    break;
              	  }
              }
          }

          if (centr > 5. && centr <= 10.) {
            	switch (p.pid()) {
                      	case 211: // pi+
                      	  {

                      	    break;
                      	  }
                      	case -211: // pi-
                      	  {

                      	    break;
                      	  }
                      	case 2212: // proton
                      	  {

                      	    break;
                      	  }
                      	case -2212: // anti-proton
                          {

                      	    break;
                      	  }
                      	case 321: // K+
                      	  {

                      	    break;
                      	  }
                      	case -321: // K-
                      	  {

                      	    break;
                      	  }
                }
            }

          if (centr > 10. && centr <= 20.) {
            switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }

          if (centr > 20. && centr <= 30.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }
          if (centr > 30. && centr <= 40.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }
          if (centr > 40. && centr <= 50.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }
          if (centr > 50. && centr <= 60.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }
          if (centr > 60. && centr <= 70.) {
            switch (p.pid()) {
              case 211: // pi+
                {

                  break;
                }
              case -211: // pi-
                {

                  break;
                }
              case 2212: // proton
                {

                  break;
                }
              case -2212: // anti-proton
                {

                  break;
                }
              case 321: // K+
                {

                  break;
                }
              case -321: // K-
                {

                  break;
                }
            }
          }
          if (centr > 70. && centr <= 80.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }
          if (centr > 80. && centr <= 90.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }


	          // particle switch
         } // primaty pi, K, p only

      } // particle loop

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double norm = sumOfWeights() *2.*M_PI;
      //scale(_h["0111"], 1./norm);

    }

    //@}


    /// name Histograms
    //{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2013_I1222333);


}
