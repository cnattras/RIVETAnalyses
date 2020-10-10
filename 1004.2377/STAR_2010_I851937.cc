// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include "../Centralities/RHICCentrality.hh"

#define _USE_MATH_DEFINES
static const int numTrigPtBins = 4;
static const float pTTrigBins[] = {2.5,3.0,4.0,6.0,10.0};
static const int numAssocPtBins = 5;
static const float pTAssocBins[] = {0.25,0.5,1.0,1.5,2.5,4.0};
static const int numCentBins = 5;
static const float CentBins[] = {0.0,12.0,20.0,40.0,60.0,80.0};
using namespace std;
namespace Rivet {

class Correlator {
      
    private:
      int _index;
      int _subindex;
      int _subsubindex;
      string _collSystemAndEnergy;
      pair<double,double> _centrality;
      pair<double,double> _triggerRange;
      pair<double,double> _associatedRange;
      vector<int> _pid;
      bool _NoPTassociated = false;

    public:
    
      /// Constructor
      Correlator(int index, int subindex, int subsubindex) {
        _index = index;
        _subindex = subindex;
        _subsubindex = subsubindex;
      }
      void SetCollSystemAndEnergy(string s){ _collSystemAndEnergy = s; }
      void SetCentrality(double cmin, double cmax){ _centrality = make_pair(cmin, cmax); }
      void SetNoPTassociated(){_NoPTassociated = true; }
      void SetTriggerRange(double tmin, double tmax){ _triggerRange = make_pair(tmin, tmax); }
      void SetAssociatedRange(double amin, double amax){ _associatedRange = make_pair(amin, amax); }
      void SetPID(std::initializer_list<int> pid){ _pid = pid; }
    
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
      vector<int> GetPID(){ return _pid; }
    
      int GetIndex(){ return _index; }
      int GetSubIndex(){ return _subindex; }
      int GetSubSubIndex(){ return _subsubindex; }
      string GetFullIndex()
      {
          string fullIndex = to_string(GetIndex()) + to_string(GetSubIndex()) + to_string(GetSubSubIndex());
          return fullIndex;
      }
    
      bool CheckCollSystemAndEnergy(string s){ return _collSystemAndEnergy.compare(s) == 0 ? true : false; }
      bool CheckCentrality(double cent){ return (cent>_centrality.first && cent<_centrality.second) ? true : false; }
      bool CheckTriggerRange(double tpt){ return (tpt>_triggerRange.first && tpt<_triggerRange.second) ? true : false; }
      bool CheckAssociatedRange(double apt){ return (apt>_associatedRange.first && apt<_associatedRange.second || _NoPTassociated == true) ? true : false; }
      bool CheckAssociatedRangeMaxTrigger(double apt, double tpt){ return (apt>_associatedRange.first && apt<tpt) ? true : false; }
      bool CheckPID(std::initializer_list<int> pid)
      {
          
          bool inList = false;
          
          for(int id : pid)
          {
              auto it = std::find(_pid.begin(), _pid.end(), id);
              
              if(it != _pid.end())
              {
                  inList = true;
                  break;
              }
          }
          
          return inList;
          
      }
      
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
  };
  
  /// @brief Add a short analysis description here
  class STAR_2010_I851937 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2010_I851937);


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

    double CalculateVn(YODA::Histo1D& hist, int nth)
    {
        int nBins = hist.numBins();

        double integral = 0.;
        
        double Vn = 0.;

        for (int i = 0; i < nBins; i++)
        {
            integral += hist.bin(i).sumW();
            Vn += hist.bin(i).sumW()*cos(nth*hist.bin(i).xMid());
        }

        Vn /= integral;
        return Vn;
    }
    
    int FindBinAtMinimum(YODA::Histo1D& hist, double bmin, double bmax)
    {
        int minBin = -1;
        double minVal = 999.;
        
        for(unsigned int i = 0; i < hist.numBins(); i++)
        {
            if(hist.bin(i).xMid() < bmin || hist.bin(i).xMid() > bmax) continue;
            if( (hist.bin(i).sumW()/hist.bin(i).xWidth()) < minVal )
            {
                minVal = hist.bin(i).sumW()/hist.bin(i).xWidth();
                minBin = i;
            }
        }
        
        return minBin;
        
    }
    
    void SubtractBackground(YODA::Histo1D& fullHist, YODA::Histo1D& hist, vector<int> n, double bmin, double bmax)
    {
        vector<double> Vn(n.size(), 0);
        for(unsigned int i = 0; i < n.size(); i++)
        {
            Vn[i] = CalculateVn(fullHist, n[i]);
        }
        
        double bmod = 1.;
        int minBin = FindBinAtMinimum(fullHist, bmin, bmax);
        
        for(unsigned int i = 0; i < Vn.size(); i++)
        {
            bmod += 2*Vn[i]*cos(n[i]*fullHist.bin(minBin).xMid());
        }
        
        double b = (fullHist.bin(minBin).sumW()/fullHist.bin(minBin).xWidth())/bmod; //Divided by bin width in order to generalize it and enable it to be used for histograms with different binning
                
        for(unsigned int ibin = 0; ibin < hist.numBins(); ibin++)
        {
            double modulation = 1;
            for(unsigned int i = 0; i < Vn.size(); i++)
            {
                modulation += 2*Vn[i]*cos(n[i]*hist.bin(ibin).xMid());
            }
            modulation *= b;
            hist.bin(ibin).scaleW(1 - (modulation/(hist.bin(ibin).sumW()/hist.bin(ibin).xWidth()))); //Divided by bin width to compensate the calculation of "b"
        }
        
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

        const ChargedFinalState cfs(Cuts::abseta < 0.35);
        declare(cfs, "CFS");
        
        const PrimaryParticles pp(pdgPi0, Cuts::abseta < 0.35);
        declare(pp, "PP");
        
        const PromptFinalState pfs(Cuts::abseta < 0.35 && Cuts::pid == 22);
        declare(pfs, "PFS"); 
     
/*/SAVE ME TO BE ON THE SAFE SIDE
      Correlator c1(1,1);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0.0, 12.0);
      c1.SetTriggerRange(2.5, 3);
      c1.SetAssociatedRange(1.0, 2.5);
      //c1.SetPID(pdgPi0);
      Correlators.push_back(c1);

	  Correlator c2(1,1);
      c1.SetCollSystemAndEnergy("dEta200GeV");
      c1.SetCentrality(0.0, 12.0);
      c1.SetTriggerRange(2.5, 3);
      c1.SetAssociatedRange(1.0, 2.5);
      //c1.SetPID(pdgPi0);
      Correlators.push_back(c2);
*/  
      //NOTES: 

      int i = 0; //raw, bksb, dEta
      int ptt = 0; //2.5-3,3-4,...
      int cb = 0; //Centralities
      for (i = 0; i < 3; i++){
      	for (ptt = 0; ptt < numTrigPtBins - 1; ptt++){
      		for (cb = 0; cb < numCentBins - 1; cb++){
      			Correlator c1 (i,ptt,cb);
      			c1.SetCollSystemAndEnergy("AuAu200GeV");
      			c1.SetCentrality(CentBins[cb], CentBins[cb+1]);
      			c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt + 1]);
      			c1.SetNoPTassociated();
      			//c1.SetPID(pdgPi0)
      			Correlators.push_back(c1);
      		}
      	}
      }
  	}	


    /* FIX ME: LOOP OVER THE YODA FILE
         //Booking the histograms
    
    	for (Correlator& corr : Correlators){
    		string name = "010101";
    		book(_h[name],01,01,01);
    		book(sow[name], "sow" + name);
    		nTriggers[name] = 0;
    	  	  }
     	}

	*/

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const PrimaryParticles& ppTrigPi0 = apply<PrimaryParticles>(event, "PP");
      const PromptFinalState& pfsTrigPhotons = apply<PromptFinalState>(event, "PFS");
      //==================================================
      // Select the histograms accordingly to the collision system, beam energy and centrality
      // WARNING: Still not implemented for d-Au
      //==================================================
      double nNucleons = 0.;
      string CollSystem = "Empty";
      const ParticlePair& beam = beams();
          CollSystem = "pp";
          nNucleons = 1.;
      //if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000010020) CollSystem = "dAu";
      //if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000290630) CollSystem = "dAu";
     
      string cmsEnergy = "Empty";
      if (fuzzyEquals(sqrtS()/GeV, 200*nNucleons, 1E-3)) cmsEnergy = "200GeV";
     
      string SysAndEnergy = CollSystem + cmsEnergy;
    
      double triggerptMin = 999.;
      double triggerptMax = -999.;
      double associatedptMin = 999.;
      double associatedptMax = -999.;
    
      bool isVeto = true;
    
      for(Correlator& corr : Correlators)
      {
        if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
        //if(!corr.CheckCentrality(centr)) continue;
        
        //If event is accepted for the correlator, fill event weights
        sow[corr.GetFullIndex()]->fill();
        
        isVeto = false;
        
        //Check min and max of the trigger and associated particles in order to speed up the particle loops
        if(corr.GetTriggerRangeMin() < triggerptMin) triggerptMin = corr.GetTriggerRangeMin();
        if(corr.GetTriggerRangeMax() > triggerptMax) triggerptMax = corr.GetTriggerRangeMax();
        
        if(corr.GetAssociatedRangeMin() < associatedptMin) associatedptMin = corr.GetAssociatedRangeMin();
        if(corr.GetAssociatedRangeMax() > associatedptMax) associatedptMax = corr.GetAssociatedRangeMax();
      }
    
      if(isVeto) vetoEvent;
    }

    /// Normalise histograms etc., after the run
    void finalize() {
	}
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> sow;
    map<string, Histo1DPtr> _DeltaPhixE;
    map<int, Histo1DPtr> _DeltaPhiSub;
    map<string, int> nTriggers;

    vector<Correlator> Correlators;
    
    std::initializer_list<int> pdgPi0 = {111, -111};  // Pion 0
    std::initializer_list<int> pdgPhoton = {22};  // Pion 0
    //@}


  };


  DECLARE_RIVET_PLUGIN(STAR_2010_I851937);

}