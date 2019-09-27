
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

  class PHENIX_2013_I1207323 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2013_I1207323);


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
				 book(_h["0111"], 1, 1, 1);
	 book(_h["0112"], 1, 1, 2);
	 book(_h["0211"], 2, 1, 1);
	 book(_h["0212"], 2, 1, 2);
	 book(_h["0311"], 3, 1, 1);
	 book(_h["0312"], 3, 1, 2);
	 book(_h["0411"], 4, 1, 1);
	 book(_h["0412"], 4, 1, 2);
	 book(_h["0511"], 5, 1, 1);
	 book(_h["0512"], 5, 1, 2);
	 book(_h["0611"], 6, 1, 1);
	 book(_h["0612"], 6, 1, 2);
	 book(_h["0711"], 7, 1, 1);
	 book(_h["0712"], 7, 1, 2);
	 book(_h["0811"], 8, 1, 1);
	 book(_h["0812"], 8, 1, 2);
	 book(_h["0911"], 9, 1, 1);
	 book(_h["0912"], 9, 1, 2);
	 book(_h["1011"], 10, 1, 1);
	 book(_h["1012"], 10, 1, 2);
	 book(_h["1111"], 11, 1, 1);
	 book(_h["1112"], 11, 1, 2);
	 book(_h["1211"], 12, 1, 1);
	 book(_h["1212"], 12, 1, 2);
	 book(_h["1311"], 13, 1, 1);
	 book(_h["1312"], 13, 1, 2);
	 book(_h["1411"], 14, 1, 1);
	 book(_h["1412"], 14, 1, 2);
	 book(_h["1511"], 15, 1, 1);
	 book(_h["1512"], 15, 1, 2);
	 book(_h["1611"], 16, 1, 1);
	 book(_h["1612"], 16, 1, 2);
	 book(_h["1711"], 17, 1, 1);
	 book(_h["1712"], 17, 1, 2);
	 book(_h["1811"], 18, 1, 1);
	 book(_h["1812"], 18, 1, 2);
	 book(_h["1911"], 19, 1, 1);
	 book(_h["1912"], 19, 1, 2);
	 book(_h["2011"], 20, 1, 1);
	 book(_h["2012"], 20, 1, 2);
	 book(_h["2111"], 21, 1, 1);
	 book(_h["2112"], 21, 1, 2);
	 book(_h["2211"], 22, 1, 1);
	 book(_h["2212"], 22, 1, 2);
	 book(_h["2311"], 23, 1, 1);
	 book(_h["2312"], 23, 1, 2);
	 book(_h["2411"], 24, 1, 1);
	 book(_h["2412"], 24, 1, 2);
	 book(_h["2511"], 25, 1, 1);
	 book(_h["2512"], 25, 1, 2);
	 book(_h["2611"], 26, 1, 1);
	 book(_h["2612"], 26, 1, 2);
	 book(_h["2711"], 27, 1, 1);
	 book(_h["2712"], 27, 1, 2);
	 book(_h["2811"], 28, 1, 1);
	 book(_h["2812"], 28, 1, 2);
	 book(_h["2911"], 29, 1, 1);
	 book(_h["2912"], 29, 1, 2);
	 book(_h["3011"], 30, 1, 1);
	 book(_h["3012"], 30, 1, 2);
	 book(_h["3111"], 31, 1, 1);
	 book(_h["3112"], 31, 1, 2);
	 book(_h["3113"], 31, 1, 3);
	 book(_h["3114"], 31, 1, 4);
	 book(_h["3115"], 31, 1, 5);
	 book(_h["3116"], 31, 1, 6);
	 book(_h["3117"], 31, 1, 7);
	 book(_h["3211"], 32, 1, 1);
	 book(_h["3212"], 32, 1, 2);
	 book(_h["3213"], 32, 1, 3);
	 book(_h["3214"], 32, 1, 4);
	 book(_h["3215"], 32, 1, 5);
	 book(_h["3216"], 32, 1, 6);
	 book(_h["3311"], 33, 1, 1);
	 book(_h["3312"], 33, 1, 2);
	 book(_h["3411"], 34, 1, 1);
	 book(_h["3412"], 34, 1, 2);
	 book(_h["3511"], 35, 1, 1);
	 book(_h["3512"], 35, 1, 2);
	 book(_h["3611"], 36, 1, 1);
	 book(_h["3612"], 36, 1, 2);
	 book(_h["3711"], 37, 1, 1);
	 book(_h["3712"], 37, 1, 2);
	 book(_h["3811"], 38, 1, 1);
	 book(_h["3911"], 39, 1, 1);
	 book(_h["4011"], 40, 1, 1);
	 book(_h["4111"], 41, 1, 1);
	 book(_h["4211"], 42, 1, 1);
	 book(_h["4311"], 43, 1, 1);
	 book(_h["4411"], 44, 1, 1);
	 book(_h["4511"], 45, 1, 1);
	 book(_h["4611"], 46, 1, 1);
	 book(_h["4711"], 47, 1, 1);
	 book(_h["4811"], 48, 1, 1);
	 book(_h["4911"], 49, 1, 1);
	 book(_h["5011"], 50, 1, 1);
	 book(_h["5111"], 51, 1, 1);
	 book(_h["5211"], 52, 1, 1);
	 book(_h["5311"], 53, 1, 1);
	 book(_h["5411"], 54, 1, 1);
	 book(_h["5511"], 55, 1, 1);
	 book(_h["5611"], 56, 1, 1);
	 book(_h["5711"], 57, 1, 1);
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
  DECLARE_RIVET_PLUGIN(PHENIX_2013_I1207323);


}