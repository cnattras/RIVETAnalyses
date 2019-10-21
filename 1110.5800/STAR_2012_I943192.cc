
// -*- C++ -*-
#include "Rivet/Analysis.hh"
//#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector> 
#define _USE_MATH_DEFINES
static const int numTrigPtBins = 7;
static const float pTTrigBins[] = {2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0};
static const int numAssocPtBins = 9;
static const float pTAssocBins[] = {1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0};

using namespace std;
namespace Rivet {

  class STAR_2012_I943192 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2012_I943192);


    /// name Analysis methods

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // the basic final-state projection: all final-state particles within the given eta acceptance

      const ChargedFinalState cfs(Cuts::abseta < 1.0 && Cuts::pT > 1*GeV);
      declare(cfs, "CFS");
      const ChargedFinalState cfsTrig(Cuts::abseta < 1.0 && Cuts::pT > 2*GeV);
      declare(cfsTrig, "CFSTrig");

      // Declare centrality projection
      //LATER FIX TO USE STAR
      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Book histograms
      //efficiency - this was just a plot to demonstrate the method
      //book(_h["0111"], 1, 1, 1);
      //correlation functions in Au-Au before and after track merging, demonstrates the method
      //book(_h["0211"], 2, 1, 1);
      //book(_h["0212"], 2, 1, 2);
      //sample correlation functions in delta eta before background subtraction
      book(_h["0311"], 3, 1, 1);
      book(_h["0312"], 3, 1, 2);
	 book(_h["0313"], 3, 1, 3);
	 book(_h["0314"], 3, 1, 4);
	 book(_h["0315"], 3, 1, 5);
	 book(_h["0316"], 3, 1, 6);
	 //sample correlation functions in delta phi before background subtraction
	 book(_h["0411"], 4, 1, 1);
	 book(_h["0412"], 4, 1, 2);
	 book(_h["0413"], 4, 1, 3);
	 book(_h["0414"], 4, 1, 4);
	 book(_h["0415"], 4, 1, 5);
	 book(_h["0416"], 4, 1, 6);
	 //sample correlation functions in delta eta after background subtraction	 
	 book(_h["0511"], 5, 1, 1);
	 book(_h["0512"], 5, 1, 2);
	 book(_h["0513"], 5, 1, 3);
	 book(_h["0514"], 5, 1, 4);
	 book(_h["0515"], 5, 1, 5);
	 book(_h["0516"], 5, 1, 6);
	 //sample correlation functions in delta phi before background subtraction
	 book(_h["0611"], 6, 1, 1);
	 book(_h["0612"], 6, 1, 2);
	 book(_h["0613"], 6, 1, 3);
	 book(_h["0614"], 6, 1, 4);
	 book(_h["0615"], 6, 1, 5);
	 book(_h["0616"], 6, 1, 6);
	 //
	 book(_h["0711"], 7, 1, 1);
	 book(_h["0811"], 8, 1, 1);
	 book(_h["0911"], 9, 1, 1);
	 book(_h["1011"], 10, 1, 1);
	 book(_h["1111"], 11, 1, 1);
	 book(_h["1211"], 12, 1, 1);
	 book(_h["1212"], 12, 1, 2);
	 book(_h["1213"], 12, 1, 3);
	 book(_h["1214"], 12, 1, 4);
	 book(_h["1215"], 12, 1, 5);
	 book(_h["1216"], 12, 1, 6);
	 book(_h["1217"], 12, 1, 7);
	 book(_h["1218"], 12, 1, 8);
	 book(_h["1311"], 13, 1, 1);
	 book(_h["1312"], 13, 1, 2);
	 book(_h["1313"], 13, 1, 3);
	 book(_h["1314"], 13, 1, 4);
	 book(_h["1315"], 13, 1, 5);
	 book(_h["1316"], 13, 1, 6);
	 book(_h["1317"], 13, 1, 7);
	 book(_h["1318"], 13, 1, 8);
	 book(_h["1411"], 14, 1, 1);
	 book(_h["1412"], 14, 1, 2);
	 book(_h["1413"], 14, 1, 3);
	 book(_h["1414"], 14, 1, 4);
	 book(_h["1415"], 14, 1, 5);
	 book(_h["1416"], 14, 1, 6);
	 book(_h["1417"], 14, 1, 7);
	 book(_h["1418"], 14, 1, 8);
	 book(_h["1511"], 15, 1, 1);
	 book(_h["1512"], 15, 1, 2);
	 book(_h["1513"], 15, 1, 3);
	 book(_h["1514"], 15, 1, 4);
	 book(_h["1515"], 15, 1, 5);
	 book(_h["1516"], 15, 1, 6);
	 book(_h["1517"], 15, 1, 7);
	 book(_h["1518"], 15, 1, 8);
	 book(_h["1611"], 16, 1, 1);
	 book(_h["1711"], 17, 1, 1);
	 book(_h["1811"], 18, 1, 1);
	 book(_h["1911"], 19, 1, 1);
	 book(_h["2011"], 20, 1, 1);
	 book(_h["2111"], 21, 1, 1);
	 book(_h["2112"], 21, 1, 2);
	 book(_h["2113"], 21, 1, 3);
	 book(_h["2114"], 21, 1, 4);
	 book(_h["2115"], 21, 1, 5);
	 book(_h["2116"], 21, 1, 6);
	 book(_h["2117"], 21, 1, 7);
	 book(_h["2118"], 21, 1, 8);
	 book(_h["2211"], 22, 1, 1);
	 book(_h["2212"], 22, 1, 2);
	 book(_h["2213"], 22, 1, 3);
	 book(_h["2214"], 22, 1, 4);
	 book(_h["2215"], 22, 1, 5);
	 book(_h["2216"], 22, 1, 6);
	 book(_h["2217"], 22, 1, 7);
	 book(_h["2218"], 22, 1, 8);
	 book(_h["2311"], 23, 1, 1);
	 book(_h["2411"], 24, 1, 1);
	 book(_h["2511"], 25, 1, 1);
	 book(_h["2611"], 26, 1, 1);
	 book(_h["2711"], 27, 1, 1);
	 book(_h["2811"], 28, 1, 1);
	 book(_h["2911"], 29, 1, 1);
	 book(_h["3011"], 30, 1, 1);
	 book(_h["3111"], 31, 1, 1);
	 book(_h["3211"], 32, 1, 1);
	 book(_h["3311"], 33, 1, 1);
	 book(_h["3411"], 34, 1, 1);
	 book(_h["3511"], 35, 1, 1);
	 book(_h["3611"], 36, 1, 1);
	 book(_h["3711"], 37, 1, 1);
	 book(_h["3811"], 38, 1, 1);
	 book(_h["3911"], 39, 1, 1);
	 book(_h["4011"], 40, 1, 1);
	 //Add declaration of histograms for correlation functions with all correlation functions
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
      cout<<"Hi"<<endl;

      vector<int> triggerPhis[numTrigPtBins];
    // loop over charged final state particles
      for(const Particle& pTrig : cfsTrig.particles()) {

	       // protections against mc generators decaying long-lived particles
	        if (!( pTrig.hasAncestor(310) || pTrig.hasAncestor(-310)  ||     // K0s
	           pTrig.hasAncestor(130)  || pTrig.hasAncestor(-130)  ||     // K0l
	           pTrig.hasAncestor(3322) || pTrig.hasAncestor(-3322) ||     // Xi0
	           pTrig.hasAncestor(3122) || pTrig.hasAncestor(-3122) ||     // Lambda
	           pTrig.hasAncestor(3222) || pTrig.hasAncestor(-3222) ||     // Sigma+/-
	           pTrig.hasAncestor(3312) || pTrig.hasAncestor(-3312) ||     // Xi-/+
	           pTrig.hasAncestor(3334) || pTrig.hasAncestor(-3334) ))    // Omega-/+
	        {
		  if( abs(pTrig.pid())==211 || abs(pTrig.pid())==2212 || abs(pTrig.pid())==321){



		    for(const Particle& pAssoc : cfs.particles()) {
		      
		      // protections against mc generators decaying long-lived particles
		      if (!( pAssoc.hasAncestor(310) || pAssoc.hasAncestor(-310)  ||     // K0s
			     pAssoc.hasAncestor(130)  || pAssoc.hasAncestor(-130)  ||     // K0l
			     pAssoc.hasAncestor(3322) || pAssoc.hasAncestor(-3322) ||     // Xi0
			     pAssoc.hasAncestor(3122) || pAssoc.hasAncestor(-3122) ||     // Lambda
			     pAssoc.hasAncestor(3222) || pAssoc.hasAncestor(-3222) ||     // Sigma+/-
			     pAssoc.hasAncestor(3312) || pAssoc.hasAncestor(-3312) ||     // Xi-/+
			     pAssoc.hasAncestor(3334) || pAssoc.hasAncestor(-3334) ))    // Omega-/+
			{
			  if( abs(pAssoc.pid())==211 || abs(pAssoc.pid())==2212 || abs(pAssoc.pid())==321){

			    //https://rivet.hepforge.org/code/dev/structRivet_1_1DeltaPhiInRange.html

			    double dPhi = deltaPhi(pTrig, pAssoc, true);//this does NOT rotate the delta phi to be in a given range
			    double dEta = deltaEta(pTrig, pAssoc);
			    cout<<" pTTrig "<<pTrig.momentum().pt()<<" pTassoc "<<pAssoc.momentum().pt()<<"dphi "<<dPhi<<" dEta "<<dEta<<endl;
			    //for this analysis we divide by the acceptance in delta eta as well when we fill
			    double acceptance  = -dEta/2.0 + 1.0;
			    //fill histogram with 1.0/acceptance
			  }
			}
		    }


		    //centr is centrality in percentile
		    //triggerPhis.push_back(p.momentum().pt()); 


		  }	          // particle switch
         } // primary pi, K, p only

      } // particle loop

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      //normalize correlation histograms by scaling by 1.0/(Ntrig*binwidthphi*binwidtheta) in each bin BUT also be careful when rebinning.  Probably best to FIRST add histograms for correlation functions THEN normalize
      //do background subtraction ala zyam
      //calculate yields

      //double norm = sumOfWeights() *2.*M_PI;
      //scale(_h["0111"], 1./norm);

    }

    //@}


    /// name Histograms
    //{
    map<string, Histo1DPtr> _h;
    //map<string, Profile1DPtr> _p;
    //map<string, CounterPtr> _c;
    //}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(STAR_2012_I943192);


}
