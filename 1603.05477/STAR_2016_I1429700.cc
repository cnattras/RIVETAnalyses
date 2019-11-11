
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
static const int numTrigPtBins = 8;
static const float pTTrigBins[] = {2.0,2.5,3.0,3.5,4.0,4.5,5.0};
static const int numAssocPtBins = 5;
static const float pTAssocBins[] = {1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0};
static const int numCentBins = 5;
static const int centBins[] = {0,12,20,30,40,50,60,70,80,90};

using namespace std;
namespace Rivet {

  class STAR_2016_I1429700 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2016_I1429700);


    /// name Analysis methods

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // the basic final-state projection: all final-state particles within the given eta acceptance

      const ChargedFinalState cfs(Cuts::abseta < 5 && Cuts::pT > 100*MeV);
      declare(cfs, "CFS");

      // Declare centrality projection
      //declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Book histograms 
	// book(_h["0111"], 1, 1, 1);
	//(Fig.2)
	//correlations functions for charged baryons and charged mesons for 0-20% Cu-Cu collisions
	//and 40-80% Au+Au collisions  
	 book(_h["0211"], 2, 1, 1);
	 book(_h["0212"], 2, 1, 2);
	 book(_h["0213"], 2, 1, 3);
	 book(_h["0214"], 2, 1, 4);
	 book(_h["0215"], 2, 1, 5);
	 book(_h["0216"], 2, 1, 6);
	//(Fig.3)
	//ratio of baryons/mesons in 0-60% Cu+Cu collisions
	 book(_h["0311"], 3, 1, 1);
	//(Fig.4)
	//jet-like yeild as a function of trigger particles for charged mesons and charged baryons
	//in minimum bias d+Au collisions 
	 book(_h["0411"], 4, 1, 1);
	 book(_h["0412"], 4, 1, 2);
	//(Fig.4)
	//jet-like yeild as a function of trigger particles for charged mesons and charged baryons
	//in 0-60% Cu+Cu collision
	 book(_h["0511"], 5, 1, 1);
	 book(_h["0512"], 5, 1, 2);
	//(Fig.4)
	//jet-like yeild as a function of trigger particles for charged mesons and charged baryons
	//in 40-80% Au+Au collision
	 book(_h["0611"], 6, 1, 1);
	 book(_h["0612"], 6, 1, 2);
	//(Fig.4)
	//jet-like yeild as a function of trigger particles for charged mesons and charged baryons
	//in 0-12% Au+Au collision
	 book(_h["0711"], 7, 1, 1);
	//(Fig.5)
	//centrality dependence of the jet-like yield of charged mesons and charged baryons in Cu+Cu
	//collisions 
	 book(_h["0811"], 8, 1, 1);
	//(Fig.5)
	//centrality dependence of the jet-like yield of charged mesons and charged baryons in Au+Au
	//collisions
	 book(_h["0911"], 9, 1, 1);
	//(Fig.5)
	//centrality dependence of the jet-like yield of charged mesons and charged baryons in d+Au
	//collisions
	 book(_h["1011"], 10, 1, 1);
	//(Fig.6)
	//jet-like yeild as a function of associated trigger particles for charged mesons and charged
	//baryons in d+Au and 0-60% Cu+Cu collisions  
	 book(_h["1111"], 11, 1, 1);
	 book(_h["1112"], 11, 1, 2);
	 book(_h["1113"], 11, 1, 3);
	 book(_h["1114"], 11, 1, 4);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const ChargedFinalState& cfsTrig = apply<ChargedFinalState>(event, "CFSTrig");

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
  DECLARE_RIVET_PLUGIN(STAR_2016_I1429700);


}
