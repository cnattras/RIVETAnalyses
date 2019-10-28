
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector> 
#define _USE_MATH_DEFINES
// Delta_phi(rad) bins
static const int numDelta_phiBins = 29;
static const float Delta_phiBins[] = {-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6};
// zT bins
static const int numzTBins = 8;
static const float zTBins[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
// zT_corr bins 
static const int numzT_corrBins = 12;
static const float zT_corrBins[] = {0.09,0.14,0.19,0.24,0.29,0.34,0.39,0.44,0.49,0.54,0.59,0.64,0.69};
// pT_trig(GeV/c) bins 
static const int numTrigPtBins = 3;
static const float pTTrigBins[] = {9.0,12.0,15.0,18.0}; 
// pT_assoc(GeV/c) bins 
static const int numAssocPtBins = 3;    
static const float pTAssocBins[] = {0.0,3.0,6.0,9.0};  
// centrality bins
static const int numCentBins = 9;
static const int centBins[] = {0,5,10,20,30,40,50,60,70,80};

using namespace std;
namespace Rivet {

  class STAR_2016_I1442357 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2016_I1442357);


    /// name Analysis methods
int GetTrigBin(float trigpT){
      if(trigpT<pTTrigBins[0]){
  cerr<<"Warning: trigpT "<<trigpT<<" is less than the minimum trigger momentum!"<<endl;
  return -1;
      }
      for(int i=0;i<numTrigPtBins;i++){
  if(trigpT<pTTrigBins[i+1]) return i;
      }
      cerr<<"Warning: trigpT "<<trigpT<<" is greater than the maximum trigger momentum!"<<endl;
      return -1;
    }
    int GetAssocBin(float assocpT){
      if(assocpT<pTAssocBins[0]){
  cerr<<"Warning: assocpT "<<assocpT<<" is less than the minimum associated momentum! "<<endl;
  return -1;
      }
      for(int i=0;i<numAssocPtBins;i++){
        if(assocpT<pTAssocBins[i+1]) return i;
      }
      cerr<<"Warning: assocpT "<<assocpT<<" is greater than the maximum associated momentum!"<<endl;
      return -1;
    }
    int GetCentBin(float cent){
      if(cent<centBins[0]){
  cerr<<"Warning: cent "<<cent<<" is less than the minimum centrality bin! "<<endl;
  return -1;
      }
      for(int i=0;i<numCentBins;i++){
  if(cent<centBins[i+1]) return i;
      }
      cerr<<"Warning: cent "<<cent<<" is greater than the maximum centrality bin!"<<endl;
      return -1;
    }
    /// Book histograms and initialise projections before the run
    void init() {              /// ******** ADD PROJECTION FOR PHOTONS ***********************

      // Initialise and register projections

      // the basic final-state projection: all final-state particles within the given eta acceptance

      const ChargedFinalState cfs(Cuts::abseta < 5 && Cuts::pT > 100*MeV); /////// LOOK IN PAPER TO FIGURE OUT CUTS**********
      declare(cfs, "CFS");
      ////const ChargedFinalState cfs(Cuts::abseta < 1.0 && Cuts::pT > 1*GeV);
      ////declare(cfs, "CFS");
      ////const ChargedFinalState cfsTrig(Cuts::abseta < 1.0 && Cuts::pT > 2*GeV);
      ////declare(cfsTrig, "CFSTrig");
      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      // pi0 also


      // Declare centrality projection
      //LATER FIX TO USE STAR
      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M"); /// *****CHANGE ALICE... *******

      // Book histograms

	 //*****FIGURE 1, Delta_phi(rad) // The Azimuthal correlation functions of charged hadrons per trigger,12 < pT^trig < 20 GeV/c
   //Data from Fig. 1a
   book(_h["0111"], 1, 1, 1); //AuAu, 1.2<pT^assoc<3 GeV/c, gamma
   book(_h["0211"], 2, 1, 1); //AuAu, 1.2<pT^assoc<3 GeV/c, pi0
   //Data from Fig. 1b
	 book(_h["0311"], 3, 1, 1); //AuAu, 3<pT^assoc< 5 GeV/c, gamma
   book(_h["0411"], 4, 1, 1); //AuAu, 3<pT^assoc<5 GeV/c, pi0
   //Data from Fig. 1c
   book(_h["0511"], 5, 1, 1); //pp, 1.2<pT^assoc<3 GeV/c, gamma
   book(_h["0611"], 6, 1, 1); //pp, 1.2<pT^assoc<3 GeV/c, pi0
   //Data from Fig. 1d
   book(_h["0711"], 7, 1, 1); //pp, 3<pT^assoc<5 GeV/c, gamma
   book(_h["0811"], 8, 1, 1); //pp, 3<pT^assoc<5 GeV/c, pi0
   
   //*****FIGURE 2, zT //pi0-hadron correlations**************************************************
   //Data from Fig. 2a
   book(_h["0911"], 9, 1, 1);  //Away-side, Au+Au
   book(_h["1011"], 10, 1, 1); //Away-side, p+p
   //Data from Fig. 2b
   book(_h["1111"], 11, 1, 1); //Near-side, Au+Au
   book(_h["1211"], 12, 1, 1); //Near-side, p+p
   
   //*****FIGURE 3, zT //Direct Photon-hadron correlations****************************************
   book(_h["1311"], 13, 1, 1); //Away-side, AuAu 
   book(_h["1411"], 14, 1, 1); //Away-side, pp
   
   //*****FIGURE 4, zT_corr **********************************************************************
   book(_h["1511"], 15, 1, 1); //Pi0-hadron correlations, Away-side, pp
   
   //*****FIGURE 5, zT ***************************************************************************
   book(_h["1611"], 16, 1, 1); //gamma, Away-side, I_AA
   book(_h["1711"], 17, 1, 1); //pi0, Away-side, I_AA
   
   //*****FIGURE 6, zT ***************************************************************************
   book(_h["1811"], 18, 1, 1); //gamma
   book(_h["1911"], 19, 1, 1); //pi0
   
   //*****FIGURE 7 *******************************************************************************
   book(_h["2011"], 20, 1, 1); //gamma, Away-side, I_AA (pT^trig)
   book(_h["2111"], 21, 1, 1); //gamma, Away-side, I_AA (pT^assoc)
   
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      ///////const ChargedFinalState& cfsTrig = apply<ChargedFinalState>(event, "CFSTrig");


      /// DO NOT MESS WITH THIS
      // Prepare centrality projection and value
      const CentralityProjection& centrProj = apply<CentralityProjection>(event, "V0M");
      double centr = centrProj();

      // Veto event for too large centralities since those are not used in the analysis at all
      if ((centr < 0.) || (centr > 80.)){
        vetoEvent;
      }
      cout<<"Hi"<<endl;


      //////////vector<int> triggerPhis[numTrigPtBins]
      

      ///////////////////////////////// STOPPING POINT ///////////////////////////////////////
    

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
  DECLARE_RIVET_PLUGIN(STAR_2016_I1442357);


}
