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
#include "RHICCentrality.hh"
#define _USE_MATH_DEFINES

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
    
    int GetIndex(){ return _index; }
    
    bool CheckCollSystemAndEnergy(string s){ return _collSystemAndEnergy.compare(s) == 0 ? true : false; }
    bool CheckCentrality(double cent){ return (cent>_centrality.first && cent<_centrality.second) ? true : false; }
    bool CheckTriggerRange(double tpt){ return (tpt>_triggerRange.first && tpt<_triggerRange.second) ? true : false; }
    bool CheckAssociatedRange(double apt){ return (apt>_associatedRange.first && apt<_associatedRange.second) ? true : false; }
    bool CheckAssociatedRangeMaxTrigger(double apt, double tpt){ return (apt>_associatedRange.first && apt<tpt) ? true : false; }
    bool CheckConditions(string s, double cent, double tpt, double apt)
    {
        if(!CheckConditions(s, cent, tpt)) return false;
        if(!CheckAssociatedRange(apt)) return false;
        
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

  
  };

  class STAR_2012_I943192 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2012_I943192);


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

      const ChargedFinalState cfs(Cuts::pT > 1*GeV); //Not cutting in eta, so no need to correct for pair acceptance
      declare(cfs, "CFS");
      const ChargedFinalState cfsTrig(Cuts::abseta < 1.0 && Cuts::pT > 2*GeV);
      declare(cfsTrig, "CFSTrig");

      // Declare centrality projection
      //LATER FIX TO USE STAR
      //declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");
      declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

      //==================================================
      // Create one correlator for each set of Collisions System / Beam Energy / Centrality Interval / Trigger pT interval / Associated pT interval
      // The number used in the constructor is the y-axis index of the histograms (d00-x00-y00) of the paper
      // Ex.: Correlator c1(1); -> is the correlator for histograms _h["0311"], _h["0411"], etc
      // Ex.: Correlator c2(2); -> is the correlator for histograms _h["0312"], _h["0412"], etc
       //==================================================
      
      Correlator c1(1);
      c1.SetCollSystemAndEnergy("CuCu62GeV");
      c1.SetCentrality(0., 60.);
      c1.SetTriggerRange(3., 6.);
      c1.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c1);
      
      Correlator c2(2);
      c2.SetCollSystemAndEnergy("AuAu62GeV");
      c2.SetCentrality(0., 80.);
      c2.SetTriggerRange(3., 6.);
      c2.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c2);
      
      Correlator c3(3);
      c3.SetCollSystemAndEnergy("dAu200GeV");
      c3.SetCentrality(0., 95.);
      c3.SetTriggerRange(3., 6.);
      c3.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c3);
      
      Correlator c4(4);
      c4.SetCollSystemAndEnergy("CuCu200GeV");
      c4.SetCentrality(0., 60.);
      c4.SetTriggerRange(3., 6.);
      c4.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c4);
      
      Correlator c5(5);
      c5.SetCollSystemAndEnergy("AuAu200GeV");
      c5.SetCentrality(40., 80.);
      c5.SetTriggerRange(3., 6.);
      c5.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c5);
      
      Correlator c6(6);
      c6.SetCollSystemAndEnergy("AuAu200GeV");
      c6.SetCentrality(0., 12.);
      c6.SetTriggerRange(3., 6.);
      c6.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c6);
      
      /*
      book(_h["0311"], 3, 1, 1);
      book(_h["0312"], 3, 1, 2);
      book(_h["0313"], 3, 1, 3);
      book(_h["0314"], 3, 1, 4);
      book(_h["0315"], 3, 1, 5);
      book(_h["0316"], 3, 1, 6);
      //sample correlation functions in delta phi before background subtraction (Fig. 6)
      book(_h["0411"], 4, 1, 1);
      book(_h["0412"], 4, 1, 2);
      book(_h["0413"], 4, 1, 3);
      book(_h["0414"], 4, 1, 4);
      book(_h["0415"], 4, 1, 5);
      book(_h["0416"], 4, 1, 6);
      //sample correlation functions in delta eta after background subtraction	  (Fig. 7)
      book(_h["0511"], 5, 1, 1);
      book(_h["0512"], 5, 1, 2);
      book(_h["0513"], 5, 1, 3);
      book(_h["0514"], 5, 1, 4);
      book(_h["0515"], 5, 1, 5);
      book(_h["0516"], 5, 1, 6);
      //sample correlation functions in delta phi before background subtraction (Fig. 7)
      book(_h["0611"], 6, 1, 1);
      book(_h["0612"], 6, 1, 2);
      book(_h["0613"], 6, 1, 3);
      book(_h["0614"], 6, 1, 4);
      book(_h["0615"], 6, 1, 5);
      book(_h["0616"], 6, 1, 6);
      //yield vs Npart (Fig. 8)
      book(_h["0711"], 7, 1, 1);
      book(_h["0811"], 8, 1, 1);
      book(_h["0911"], 9, 1, 1);
      book(_h["1011"], 10, 1, 1);
      book(_h["1111"], 11, 1, 1);
      //yield vs pTtrig (Fig. 9)
      book(_h["1211"], 12, 1, 1);
      book(_h["1212"], 12, 1, 2);
      book(_h["1213"], 12, 1, 3);
      book(_h["1214"], 12, 1, 4);
      book(_h["1215"], 12, 1, 5);
      book(_h["1216"], 12, 1, 6);
      book(_h["1217"], 12, 1, 7);
      book(_h["1218"], 12, 1, 8);
      //Yield vs pTassoc (Fig. 10).  Commented out histograms are PYTHIA
      //book(_h["1311"], 13, 1, 1);
      book(_h["1312"], 13, 1, 2);
      book(_h["1313"], 13, 1, 3);
      //book(_h["1314"], 13, 1, 4);
      book(_h["1315"], 13, 1, 5);
      book(_h["1316"], 13, 1, 6);
      book(_h["1317"], 13, 1, 7);
      book(_h["1318"], 13, 1, 8);
      //Delta Phis Widths vs pTtrigger (Fig. 11a).  PYTHIA commended out.
      //book(_h["1411"], 14, 1, 1);
      book(_h["1412"], 14, 1, 2);
      book(_h["1413"], 14, 1, 3);
      //book(_h["1414"], 14, 1, 4);
      book(_h["1415"], 14, 1, 5);
      book(_h["1416"], 14, 1, 6);
      book(_h["1417"], 14, 1, 7);
      book(_h["1418"], 14, 1, 8);
      //Delta Phis Widths vs pTassoc (Fig. 11b).  PYTHIA commended out.
      //book(_h["1511"], 15, 1, 1);
      book(_h["1512"], 15, 1, 2);
      book(_h["1513"], 15, 1, 3);
      //book(_h["1514"], 15, 1, 4);
      book(_h["1515"], 15, 1, 5);
      book(_h["1516"], 15, 1, 6);
      book(_h["1517"], 15, 1, 7);
      book(_h["1518"], 15, 1, 8);
      //Delta Phis Widths vs Npart (Fig. 11c).
      book(_h["1611"], 16, 1, 1);
      book(_h["1711"], 17, 1, 1);
      book(_h["1811"], 18, 1, 1);
      book(_h["1911"], 19, 1, 1);
      book(_h["2011"], 20, 1, 1);
      //Delta Eta Widths vs pTtrig (Fig. 11d).  PYTHIA commended out.
      //book(_h["2111"], 21, 1, 1);
      book(_h["2112"], 21, 1, 2);
      book(_h["2113"], 21, 1, 3);
      //book(_h["2114"], 21, 1, 4);
      book(_h["2115"], 21, 1, 5);
      book(_h["2116"], 21, 1, 6);
      book(_h["2117"], 21, 1, 7);
      book(_h["2118"], 21, 1, 8);
      //Delta Eta Widths vs pTassoc (Fig. 11e).  PYTHIA commended out.
      //book(_h["2211"], 22, 1, 1);
      book(_h["2212"], 22, 1, 2);
      book(_h["2213"], 22, 1, 3);
      //book(_h["2214"], 22, 1, 4);
      book(_h["2215"], 22, 1, 5);
      book(_h["2216"], 22, 1, 6);
      book(_h["2217"], 22, 1, 7);
      book(_h["2218"], 22, 1, 8);
      //Delta Eta Widths vs Npart (Fig. 11f)
      book(_h["2311"], 23, 1, 1);
      book(_h["2411"], 24, 1, 1);
      book(_h["2511"], 25, 1, 1);
      book(_h["2611"], 26, 1, 1);
      book(_h["2711"], 27, 1, 1);
      //Figure 12 - ridge yields, obsolete, not implementing!
      //book(_h["2811"], 28, 1, 1);
      //book(_h["2911"], 29, 1, 1);
      //book(_h["3011"], 30, 1, 1);
      //book(_h["3111"], 31, 1, 1);
      //Figure 13 - ridge/jet yields, obsolete, not implementing!
      //book(_h["3211"], 32, 1, 1);
      //book(_h["3311"], 33, 1, 1);
      //book(_h["3411"], 34, 1, 1);
      //book(_h["3511"], 35, 1, 1);
      //Fig. 14 v3^2/v2^2, mostly obsolete but not easy to implement anyways - not implementing!
      //book(_h["3611"], 36, 1, 1);
      //book(_h["3711"], 37, 1, 1);
      //book(_h["3811"], 38, 1, 1);
      //book(_h["3911"], 39, 1, 1);
      //book(_h["4011"], 40, 1, 1);
      //Add declaration of histograms for correlation functions with all correlation functions
      */
      
      for(unsigned int i = 1; i<= Correlators.size(); i++)
      {
          book(sow[i],"sow" + to_string(i));
          book(_h["031" + to_string(i)], 3, 1, i);
          book(_h["041" + to_string(i)], 4, 1, i);
      }
      
      nEvents.assign(Correlators.size()+1, 0); 
      nTriggers.assign(Correlators.size()+1, 0); 
      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const ChargedFinalState& cfsTrig = apply<ChargedFinalState>(event, "CFSTrig");

      //==================================================
      // Select the histograms accordingly to the collision system, beam energy and centrality
      // WARNING: Still not implemented for d-Au
      //==================================================
      double nNucleons = 0.;
      string CollSystem = "Empty";
      const ParticlePair& beam = beams();
      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
          CollSystem = "AuAu";
          nNucleons = 197.;
      }
      else if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000290630) 
      {
          CollSystem = "CuCu";
          nNucleons = 63.;
      }
      else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) 
      {
          CollSystem = "pp";
          nNucleons = 1.;
      }
      //if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000010020) CollSystem = "dAu";
      //if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000290630) CollSystem = "dAu";
      if(CollSystem.compare("Empty") == 0) return;
      
      string cmsEnergy = "Empty";
      if (fuzzyEquals(sqrtS()/GeV, 200*nNucleons, 1E-3)) cmsEnergy = "200GeV";
      if (fuzzyEquals(sqrtS()/GeV, 62.3*nNucleons, 1E-3)) cmsEnergy = "62GeV";
      if(cmsEnergy.compare("Empty") == 0) return;
      
      string SysAndEnergy = CollSystem + cmsEnergy;

      // Prepare centrality projection and value
      const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
      double centr = centProj();
    
    double triggerptMin = 999.;
    double triggerptMax = -999.;
    double associatedptMin = 999.;
    double associatedptMax = -999.;
    
    bool isVeto = true;
    
    for(Correlator& corr : Correlators)
    {
        if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
        if(!corr.CheckCentrality(centr)) continue;
        
        //If event is accepted for the correlator, fill event weights
        sow[corr.GetIndex()]->fill();
        nEvents[corr.GetIndex()]++;
        
        isVeto = false;
        
        //Check min and max of the trigger and associated particles in order to speed up the particle loops
        if(corr.GetTriggerRangeMin() < triggerptMin) triggerptMin = corr.GetTriggerRangeMin();
        if(corr.GetTriggerRangeMax() > triggerptMax) triggerptMax = corr.GetTriggerRangeMax();
        
        if(corr.GetAssociatedRangeMin() < associatedptMin) associatedptMin = corr.GetAssociatedRangeMin();
        if(corr.GetAssociatedRangeMax() > associatedptMax) associatedptMax = corr.GetAssociatedRangeMax();
    }
    
    if(isVeto) vetoEvent;
    
    // loop over charged final state particles
      for(const Particle& pTrig : cfsTrig.particles()) {

	       if(pTrig.pt()/GeV < triggerptMin || pTrig.pt()/GeV > triggerptMax) continue;

          //Check if is secondary
          if(isSecondary(pTrig)) continue;
          
		  if( abs(pTrig.pid())==211 || abs(pTrig.pid())==2212 || abs(pTrig.pid())==321){

            for(Correlator& corr : Correlators)
            {
                if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
                nTriggers[corr.GetIndex()]++;
            }

		    for(const Particle& pAssoc : cfs.particles()) {
                
                if(pAssoc.pt()/GeV < associatedptMin || pAssoc.pt()/GeV > pTrig.pt()/GeV) continue;
                
                //Check if Trigger and Associated are the same particle
                if(isSameParticle(pTrig,pAssoc)) continue;
                
                //Check if is secondary
                if(isSecondary(pAssoc)) continue;
                
			  if( abs(pAssoc.pid())==211 || abs(pAssoc.pid())==2212 || abs(pAssoc.pid())==321){
			    //int mybin = GetTrigBin(pTrig.pt());
			    //int mybina = GetAssocBin(pAssoc.pt());

			    //https://rivet.hepforge.org/code/dev/structRivet_1_1DeltaPhiInRange.html
			    double dPhi = deltaPhi(pTrig, pAssoc, true);//this does NOT rotate the delta phi to be in a given range
			    double dEta = deltaEta(pTrig, pAssoc);
			    
                for(Correlator& corr : Correlators)
                {
                    if(!corr.CheckConditionsMaxTrigger(SysAndEnergy, centr, pTrig.pt()/GeV, pAssoc.pt()/GeV)) continue;
                    
                    if(abs(dPhi) < 0.78)
                    {
                        _h["031" + to_string(corr.GetIndex())]->fill(-abs(dEta), 0.5);
                    }
                    
                    if(abs(dEta) < 0.78)
                    {
                        _h["041" + to_string(corr.GetIndex())]->fill(-abs(dPhi), 0.5);
                    }
                    
                } //end of correlators loop 
                
			  } // associated hadrons
		    } // end of loop over associated particles
		  } // trigger hadrons
      } // particle loop

    }


    /// Normalise histograms etc., after the run
    void finalize() {
        
        for(unsigned int i = 1; i <= Correlators.size(); i++)
        {
            if(nTriggers[i] > 0)
            {
                _h["031" + to_string(i)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                _h["041" + to_string(i)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
            }
            
        }
        
      //normalize correlation histograms by scaling by 1.0/(Ntrig*binwidthphi*binwidtheta) in each bin BUT also be careful when rebinning.  Probably best to FIRST add histograms for correlation functions THEN normalize
      //do background subtraction ala zyam
      //calculate yields

      //double norm = sumOfWeights() *2.*M_PI;
      //scale(_h["0111"], 1./norm);

    }

    //Histograms and variables
    map<string, Histo1DPtr> _h;
    map<int, CounterPtr> sow;
    bool fillTrigger = true;
    vector<int> nTriggers;
    vector<int> nEvents;
    vector<Correlator> Correlators;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(STAR_2012_I943192);


}
