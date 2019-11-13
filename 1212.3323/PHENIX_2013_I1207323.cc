
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include "Rivet/Projections/RHICCentrality.hh"
#define _USE_MATH_DEFINES

static const int numDelPhiBins = 10;
static const float DelPhiBins[] = {0.16,0.47,0.79,1.10,1.41,1.73,2.04,2.36,2.67,2.98};

static const int numXiBins = 6;
static const float XiBins[] = {0.10,0.50,1.10,1.50,1.90,2.30};

static const int numTrigPtBins = 1;
static const float pTTrigBins[] = {5.0,9.0};

using namespace std;
namespace Rivet {

  class Correlator {
  
  public:
    
    /// Constructor
    Correlator(int index) {
      _index = index;
    }

    void SetCollSystemAndEnergy(string s){ _collSystemAndEnergy = s; }
    void SetCentrality(double cmin, double cmax){ _centrality = make_pair(cmin, cmax); }
    void SetTriggerRange(double tmin, double tmax){ _triggerRange = make_pair(tmin, tmax); }
    void SetAssociatedRange(double amin, double amax){ _associatedRange = make_pair(amin, amax); }
    void SetXiRange(double xmin, double xmax){ _xiRange = make_pair(xmin, xmax); }
    
    string GetCollSystemAndEnergy(){ return _collSystemAndEnergy; }
    pair<double,double> GetCentrality(){ return _centrality; }
    double GetCentralityMin(){ return _centrality.first; }
    double GetCentralityMax(){ return _centrality.second; }
    pair<double,double> GetTriggerRange(){ return _triggerRange; }
    double GetTriggerRangeMin(){ return _triggerRange.first; }
    double GetTriggerRangeMax(){ return _triggerRange.second; }
    pair<double,double> GetAssociatedRange(){ return _associatedRange; }
    double GetAssociatedRangeMin(){ return _associatedRange.first; }
    double GetAssociatedRangeMax(){ return _associatedRange.second; }
    pair<double,double> GetXiRange(){ return _xiRange; }
    double GetXiRangeMin(){ return _xiRange.first; }
    double GerXiRangeMax(){ return _xiRange.second; }
    
    int GetIndex(){ return _index; }
    
    bool CheckCollSystemAndEnergy(string s){ return _collSystemAndEnergy.compare(s) == 0 ? true : false; }
    bool CheckCentrality(double cent){ return (cent>_centrality.first && cent<_centrality.second) ? true : false; }
    bool CheckTriggerRange(double tpt){ return (tpt>_triggerRange.first && tpt<_triggerRange.second) ? true : false; }
    bool CheckAssociatedRange(double apt){ return (apt>_associatedRange.first && apt<_associatedRange.second) ? true : false; }
    bool CheckAssociatedRangeMaxTrigger(double apt, double tpt){ return (apt>_associatedRange.first && apt<tpt) ? true : false; }
    bool CheckXiRange(double xpt){ return (xpt>_xiRange.first && xpt<_xiRange.second) ? true : false; }
    bool CheckXiConditions(string s, double cent, double tpt, double apt, double xpt)
    {
        if(!CheckConditions(s, cent, tpt)) return false;
        if(!CheckAssociatedRange(apt)) return false;
        if(!CheckXiRange(xpt)) return false;
        
        return true;
        
    }
    bool CheckConditions(string s, double cent, double tpt)
    {
        if(!CheckConditions(s, cent)) return false;
        if(!CheckTriggerRange(tpt)) return false;
        
        return true;
        
    }
    bool CheckConditions(string s, double cent)
    {
        if(!CheckCollSystemAndEnergy(s)) return false;
        if(!CheckCentrality(cent)) return false;
        
        return true;
        
    }
    bool CheckConditionsMaxTrigger(string s, double cent, double tpt, double apt)
    {
        if(!CheckCollSystemAndEnergy(s)) return false;
        if(!CheckCentrality(cent)) return false;
        if(!CheckTriggerRange(tpt)) return false;
        if(!CheckAssociatedRangeMaxTrigger(apt,tpt)) return false;
        
        return true;
        
    }
    
    int _index;
    string _collSystemAndEnergy;
    pair<double,double> _centrality;
    pair<double,double> _triggerRange;
    pair<double,double> _associatedRange;
    pair<double,double> _xiRange; 

  
  };


  class PHENIX_2013_I1207323 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2013_I1207323);


    /// name Analysis methods
    bool isSameParticle(const Particle& p1, const Particle& p2)
    {
        //if pT, eta and phi are equal, they are the same particle
        if(p1.pt() != p2.pt()) return false;
        if(p1.eta() != p2.eta()) return false;
        if(p1.phi() != p2.phi()) return false;
    
        return true;
    }
    
    bool isSecondary(Particle p)
    {
        //return true if is secondary
        if (( p.hasAncestor(310) || p.hasAncestor(-310)  ||     // K0s
            p.hasAncestor(130)  || p.hasAncestor(-130)  ||     // K0l
            p.hasAncestor(3322) || p.hasAncestor(-3322) ||     // Xi0
            p.hasAncestor(3122) || p.hasAncestor(-3122) ||     // Lambda
            p.hasAncestor(3222) || p.hasAncestor(-3222) ||     // Sigma+/-
            p.hasAncestor(3312) || p.hasAncestor(-3312) ||     // Xi-/+
            p.hasAncestor(3334) || p.hasAncestor(-3334) ))    // Omega-/+
        return true;
        else return false;
        
    }

    /// Book histograms and initialise projections before the run
    void init() {
      // the basic final-state projection: all final-state particles within the given eta acceptance

      const ChargedFinalState cfs(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV);
      declare(cfs, "CFS");
      const ChargedFinalState cfsTrig(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV);
      declare(cfsTrig, "CFSTrig");

      // FinalState of prompt photons and bare muons and electrons in the event
      const ChargedFinalState pfsTrig(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV);
      //const PromptFinalState pfsTrig(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV);
      //const PromptFinalState pfsTrig(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV && Cuts::abspid == PID::PHOTON);
      declare(pfsTrig, "PFSTrig"); 
      
      //PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      //declare(pfsTrig, "PFSTrig")


      // Declare centrality projection
      // declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      //==================================================
      // Create one correlator for each set of Collisions System / Beam Energy / Centrality Interval / Trigger pT interval / Associated pT interval
      // The number used in the constructor is the y-axis index of the histograms (d00-x00-y00) of the paper
      // Ex.: Correlator c1(1); -> is the correlator for histograms _h["0311"], _h["0411"], etc
      // Ex.: Correlator c2(2); -> is the correlator for histograms _h["0312"], _h["0412"], etc
       //==================================================
      
      Correlator c1(1);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 40.);
      c1.SetXiRange(0., 0.4);
      c1.SetTriggerRange(5., 9.);
      c1.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c1);
      
      Correlator c2(2);
      c2.SetCollSystemAndEnergy("pp200GeV");
      c2.SetCentrality(0., 40.);
      c2.SetXiRange(0., 0.4);
      c2.SetTriggerRange(5., 9.);
      c2.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c2);
      
      Correlator c3(3);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0., 40.);
      c3.SetXiRange(0.4, 0.8);
      c3.SetTriggerRange(5., 9.);
      c3.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c3);
      
      Correlator c4(4);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetCentrality(0., 40.);
      c4.SetXiRange(0.4, 0.8);
      c4.SetTriggerRange(5., 9.);
      c4.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c4);
      
      Correlator c5(5);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0., 40.);
      c3.SetXiRange(0.8, 1.2);
      c3.SetTriggerRange(5., 9.);
      c3.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c5);
      
      Correlator c6(6);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetCentrality(0., 40.);
      c4.SetXiRange(0.8, 1.2);
      c4.SetTriggerRange(5., 9.);
      c4.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c6);

      Correlator c7(7);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 40.);
      c1.SetXiRange(1.2, 1.6);
      c1.SetTriggerRange(5., 9.);
      c1.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c7);
      
      Correlator c8(8);
      c2.SetCollSystemAndEnergy("pp200GeV");
      c2.SetCentrality(0., 40.);
      c2.SetXiRange(1.2, 1.6);
      c2.SetTriggerRange(5., 9.);
      c2.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c8);
      
      Correlator c9(9);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0., 40.);
      c3.SetXiRange(1.6, 2.0);
      c3.SetTriggerRange(5., 9.);
      c3.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c9);
      
      Correlator c10(10);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetCentrality(0., 40.);
      c4.SetXiRange(1.6, 2.0);
      c4.SetTriggerRange(5., 9.);
      c4.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c10);
      
      Correlator c11(11);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0., 40.);
      c3.SetXiRange(2.0, 2.4);
      c3.SetTriggerRange(5., 9.);
      c3.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c11);
      
      Correlator c12(12);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetCentrality(0., 40.);
      c4.SetXiRange(2.0, 2.4);
      c4.SetTriggerRange(5., 9.);
      c4.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c12);

      Correlator c13(13);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0., 40.);
      c3.SetXiRange(0.2, 2.2);
      c3.SetTriggerRange(5., 9.);
      c3.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c13);
      
      Correlator c14(14);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetCentrality(0., 40.);
      c4.SetXiRange(0.2, 2.2);
      c4.SetTriggerRange(5., 9.);
      c4.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c14);


      // Book histograms
      // Order of plots may be off as they appear in HEP Data 12 11 10 1-9
      // Correlation functions
      // AuAu 0<xi<0.4
 	    book(_h["0111"], 1, 1, 1);
       //pp 0<xi<0.4
    	 book(_h["0211"], 2, 1, 1);
       //AuAu 0.4<xi<0.8
    	 book(_h["0311"], 3, 1, 1);
       //pp 0.4<xi<0.8
    	 book(_h["0411"], 4, 1, 1);
       //AuAu 0.8<xi<1.2
    	 book(_h["0511"], 5, 1, 1);
       //pp 0.8<xi<1.2
    	 book(_h["0611"], 6, 1, 1);
       //AuAu 1.2<xi<1.6
    	 book(_h["0711"], 7, 1, 1);
       //pp 1.2<xi<1.6
    	 book(_h["0811"], 8, 1, 1);
       //AuAu 1.6<xi<2.0
    	 book(_h["0911"], 9, 1, 1);
       //pp 1.6<xi<2.0
    	 book(_h["1011"], 10, 1, 1);
       //AuAu 2.0<xi<2.4
    	 book(_h["1111"], 11, 1, 1);
       //pp 2.0<xi<2.4
    	 book(_h["1211"], 12, 1, 1);
       //trigger yields for 0-40% Au+Au
    	 book(_h["1311"], 13, 1, 1);
       //trigger yields for 0-40% pp
    	 book(_h["1411"], 14, 1, 1);
       //I_AA for 0-40% Au+Au/p+p
    	 book(_h["1511"], 15, 1, 1);
       //I_AA of 0-40% Au+Au/p+p for |Dphi-pi| < pi/2
    	 book(_h["1611"], 16, 1, 1);
       //I_AA of 0-40% Au+Au/p+p for |Dphi-pi| < pi/3
    	 book(_h["1711"], 17, 1, 1);
       //I_AA of 0-40% Au+Au/p+p for |Dphi-pi| < pi/6
    	 book(_h["1811"], 18, 1, 1);
       //Ratio of |Dphi-pi| < pi/2 I_AA to |Dphi-pi| < pi/6 I_AA for 0-40% Au+Au/p+p
    	 book(_h["1911"], 19, 1, 1);
      //Don't know 
      book(_h["2011"], 20, 1, 1);

      for(unsigned int i = 1; i<= Correlators.size(); i++)
      {
          if (i<10){
           book(sow[i],"sow" + to_string(i));
           book(_h["0" + to_string(i) + "11"], i, 1, 1);
          }
          else {
           book(sow[i],"sow" + to_string(i));
           book(_h[to_string(i) + "11"], i, 1, 1);
          }
      }
      
      nEvents.assign(Correlators.size()+1, 0); 
      nTriggers.assign(Correlators.size()+1, 0); 

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // Prepare centrality projection and value
      const CentralityProjection& centrProj = apply<CentralityProjection>(event, "CMULT");
      double centr = centrProj();

      // Veto event for too large centralities since those are not used in the analysis at all
      if ((centr < 0.) || (centr > 40.)){
        vetoEvent;
      }


      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

      //const PromptFinalState& pfsTrig = apply<PromptFinalState>(event, "PFSTrig");
      const ChargedFinalState& pfsTrig = apply<ChargedFinalState>(event, "PFSTrig");
      // loop over charged final state particles
      for(const Particle& p : pfsTrig.particles()) {
	       // protections against mc generators decaying long-lived particles
	        if (!( p.hasAncestor(310) || p.hasAncestor(-310)  ||     // K0s
	           p.hasAncestor(130)  || p.hasAncestor(-130)  ||     // K0l
	           p.hasAncestor(3322) || p.hasAncestor(-3322) ||     // Xi0
	           p.hasAncestor(3122) || p.hasAncestor(-3122) ||     // Lambda
	           p.hasAncestor(3222) || p.hasAncestor(-3222) ||     // Sigma+/-
	           p.hasAncestor(3312) || p.hasAncestor(-3312) ||     // Xi-/+
	           p.hasAncestor(3334) || p.hasAncestor(-3334) ))    // Omega-/+
	        {
         
         }
      } // particle loop

    }


    /// Normalise histograms etc., after the run
    void finalize() {

     // double norm = sumOfWeights() *2.*M_PI;
      //scale(_h["0111"], 1./norm);

    }

    //@}


    
    
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    
    map<int, CounterPtr> sow;
    bool fillTrigger = true;
    vector<int> nTriggers;
    vector<int> nEvents;
    vector<Correlator> Correlators;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PHENIX_2013_I1207323);


}
