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
static const int numTrigPtBins = 3;
static const float pTTrigBins[] = {4.0,6.0,8.0,15.0};
static const int numAssocPtBins = 3;
static const float pTAssocBins[] = {3.0,4.0,6.0,10};
static const int numzTBins = 7;
static const float zTBins[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
static const int numCentBins = 6;
static const float centBins[] = {0.0,5.0,10,20,30,40,80};
static const double pi = 2*acos(0.0);

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
      pair<double,double> _zTRange;
      vector<int> _pid;

    public:
    
      /// Constructor
      Correlator(int index, int subindex, int subsubindex) {
        _index = index;
        _subindex = subindex;
        _subsubindex = subsubindex;
      }

      void SetCollSystemAndEnergy(string s){ _collSystemAndEnergy = s; }
      void SetCentrality(double cmin, double cmax){ _centrality = make_pair(cmin, cmax); }
      void SetTriggerRange(double tmin, double tmax){ _triggerRange = make_pair(tmin, tmax); }
      void SetAssociatedRange(double amin, double amax){ _associatedRange = make_pair(amin, amax); }
      void SetzTRange(double amin, double amax){ _associatedRange = make_pair(amin, amax); }
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
      pair<double,double> GetzTRange(){ return _zTRange; }
      double GetzTRangeMin(){ return _zTRange.first; }
      double GetzTRangeMax(){ return _zTRange.second; }
      vector<int> GetPID(){ return _pid; }
    
      int GetIndex(){ return _index; }
      int GetSubIndex(){ return _subindex; }
      int GetSubSubIndex(){ return _subsubindex; }
      string GetFullIndex()
      {
          string fullIndex = to_string(GetIndex()) + to_string(GetSubIndex())+ to_string(GetSubSubIndex());
          return fullIndex;
      }
    
      bool CheckCollSystemAndEnergy(string s){ return _collSystemAndEnergy.compare(s) == 0 ? true : false; }
      bool CheckCentrality(double cent){ return (cent>_centrality.first && cent<_centrality.second) ? true : false; }
      bool CheckTriggerRange(double tpt){ return (tpt>_triggerRange.first && tpt<_triggerRange.second) ? true : false; }
      bool CheckAssociatedRange(double apt){ return (apt>_associatedRange.first && apt<_associatedRange.second) ? true : false; }
      bool CheckAssociatedRangeMaxTrigger(double apt, double tpt){ return (apt>_associatedRange.first && apt<tpt) ? true : false; }
      bool CheckzTRange(double apt){ return (apt>_zTRange.first && apt<_zTRange.second) ? true : false; }
      bool CheckzTRangeMaxTrigger(double apt, double tpt){ return (apt>_zTRange.first && apt<tpt) ? true : false; }
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
  class STAR_2006_I715470 : public Analysis {
  public:
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
   
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2006_I715470);
    void init() {
      const ChargedFinalState cfs(Cuts::abseta < 0.35);
      declare(cfs, "CFS");

      const PrimaryParticles pp(pdgPi0, Cuts::abseta < 0.35);
      declare(pp, "PP");

      const PromptFinalState pfs(Cuts::abseta < 0.35 && Cuts::pid == 22);
      declare(pfs, "PFS");

int ndPhiBins = 72;
double lowedge = -pi/2.0;
double highedge = 3.0*pi/2.0;


        for(int ncb=0;ncb<numCentBins;ncb++){
      for(int ntrig=0;ntrig<numTrigPtBins;ntrig++){
          for(int nassoc=0;nassoc<numAssocPtBins;nassoc++){
      //Correlators
            Correlator c1(ncb,ntrig,nassoc);
            c1.SetCollSystemAndEnergy("AuAu200GeV");
            c1.SetCentrality(centBins[ncb], centBins[ncb+1]);
            c1.SetTriggerRange(pTTrigBins[ntrig], pTTrigBins[ntrig+1]);
            c1.SetAssociatedRange(pTAssocBins[nassoc],pTAssocBins[nassoc+1]);
            c1.SetzTRange(0,1);
            Correlators.push_back(c1);
            string name = c1.GetCollSystemAndEnergy()+c1.GetFullIndex();
            book(_h[name], name, ndPhiBins,lowedge,highedge);
            book(sow[c1.GetFullIndex()],"sow" + c1.GetFullIndex());
            if(ncb==0){

              Correlator c2(100,ntrig,nassoc);
              c2.SetCollSystemAndEnergy("dAu200GeV");
              c2.SetCentrality(0,80);
              c2.SetTriggerRange(pTTrigBins[ntrig], pTTrigBins[ntrig+1]);
              c2.SetAssociatedRange(pTAssocBins[nassoc],pTAssocBins[nassoc+1]);
              c2.SetzTRange(0,1); 
              Correlators.push_back(c2);
            string name2 = c2.GetCollSystemAndEnergy()+c2.GetFullIndex();
            book(_h[name2], name2, ndPhiBins,lowedge,highedge);
            book(sow[c2.GetFullIndex()],"sow" + c2.GetFullIndex());
            }
          }


        }//end loop over trigger pT bins

          for(int nzT=0;nzT<numzTBins;nzT++){

      //Correlators
            Correlator c1(ncb,1000,1000);
            c1.SetCollSystemAndEnergy("AuAu200GeV");
            c1.SetCentrality(centBins[ncb], centBins[ncb+1]);
            c1.SetTriggerRange(8,15);
            c1.SetzTRange(zTBins[nzT],zTBins[nzT+1]);
            c1.SetAssociatedRange(3,15);
            Correlators.push_back(c1);
            string name = c1.GetCollSystemAndEnergy()+c1.GetFullIndex();
            book(_h[name], name, ndPhiBins,lowedge,highedge);
            book(sow[c1.GetFullIndex()],"sow" + c1.GetFullIndex());
            if(ncb==0){

              Correlator c2(1000,1000,1000);
              c2.SetCollSystemAndEnergy("dAu200GeV");
              c2.SetCentrality(0,80);
              c2.SetTriggerRange(8,15);
              c2.SetzTRange(zTBins[nzT],zTBins[nzT+1]);
              c2.SetAssociatedRange(3,15);
              Correlators.push_back(c2);
            string name2 = c2.GetCollSystemAndEnergy()+c2.GetFullIndex();
            book(_h[name2], name2, ndPhiBins,lowedge,highedge);
            book(sow[c2.GetFullIndex()],"sow" + c2.GetFullIndex());
            }
          }//End of the zT loop

      }//end loop over centrality bins

//AJ - book histograms which were in the .yoda file
      book(_h["010101"],1,1,1);
      book(_h["010102"],1,1,2);
      book(_h["010103"], 1, 1, 3);
      book(_h["010104"], 1, 1, 4);
      book(_h["010105"], 1, 1, 5);
      book(_h["010106"], 1, 1, 6);
      book(_h["020101"], 2, 1, 1);
      book(_h["020102"], 2, 1, 2);
      book(_h["030101"], 3, 1, 1);
      book(_h["030102"], 3, 1, 2);
      book(_h["030103"], 3, 1, 3);
      book(_h["030104"], 3, 1, 4);
      book(_h["030105"], 3, 1, 5);
      book(_h["030101"], 3, 1, 1);
      book(_h["040101"], 3, 1, 1);
      //Add counters to keep track of trigger particles - something like this
      //book(sow[corr.GetFullIndex()],"sow" + corr.GetFullIndex());

    }
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      //const PrimaryParticles& ppTrigPi0 = apply<PrimaryParticles>(event, "PP");
      //const PromptFinalState& pfsTrigPhotons = apply<PromptFinalState>(event, "PFS");
      //==================================================
      // Select the histograms accordingly to the collision system, beam energy and centrality
      // WARNING: Still not implemented for d-Au
      //==================================================
      double nNucleons = 0.;
      string CollSystem = "Empty";
      const ParticlePair& beam = beams();
          CollSystem = "AuAu";
          nNucleons = 197;
      if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000010020){ 
        CollSystem = "dAu";
          nNucleons = 1;
    }
      if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000290630){

       CollSystem = "dAu";
          nNucleons = 1;//needs checking
     }
     
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
    // loop over charged final state particles - PI
     for(const Particle& pTrig : cfs.particles()){

        //Check if is secondary
        //if(isSecondary(pTrig)) continue;
          
        
        for(Correlator& corr : Correlators)
        {//This loops over all correlators and asks if the particle is in the trigger range.  If so, add it to the list of triggers.
            //if(!corr.CheckPID(pdgPi0)) continue;
            if(!corr.CheckTriggerRange(pTrig.pT()/GeV)) continue;  
            nTriggers[corr.GetFullIndex()]++;
            cout<<"I have a trigger particle!"<<endl;
        }
        // Hadron loop
        for(const Particle& pAssoc : cfs.particles())
        {//This loops over all particles and considers them associated particles
                
            //Check if Trigger and Associated are the same particle
            if(isSameParticle(pTrig,pAssoc)) continue;
                
            //Check if is secondary
            //if(isSecondary(pAssoc)) continue;

            //https://rivet.hepforge.org/code/dev/structRivet_1_1DeltaPhiInRange.html
            double dPhi = deltaPhi(pTrig, pAssoc, true);//this does NOT rotate the delta phi to be in a given range
                        
            //double xE = GetXE(pTrig,pAssoc);
            
            for(Correlator& corr : Correlators)
            {
                //if(!corr.CheckPID(pdgPi0)) continue;
                
                if(!corr.CheckTriggerRange(pTrig.pT()/GeV)) continue;
                //AJ add a spot where you fill histograms
                //See Nora's code in 08014545 around line 1168

              string name = corr.GetCollSystemAndEnergy()+corr.GetFullIndex();
                  //if(corr.GetSubSubIndex()==-1){
                   //string name = "58010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  _h[name]->fill(dPhi);
                  //} 
                
            } //end of correlators loop 
                
        } // end of loop over associated particles

    } // particle loop


    }
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


  };


  DECLARE_RIVET_PLUGIN(STAR_2006_I715470);

}
