
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
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

using namespace std;
namespace Rivet {
    
  class Correlator {
      
    private:

      int _index;
      int _subindex;
      string _collSystemAndEnergy;
      pair<double,double> _centrality;
      pair<double,double> _triggerRange;
      pair<double,double> _associatedRange;
      vector<int> _pid;

    public:
    
      /// Constructor
      Correlator(int index, int subindex) {
        _index = index;
        _subindex = subindex;
      }

      void SetCollSystemAndEnergy(string s){ _collSystemAndEnergy = s; }
      void SetCentrality(double cmin, double cmax){ _centrality = make_pair(cmin, cmax); }
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
      string GetFullIndex()
      {
          string fullIndex = to_string(GetIndex()) + to_string(GetSubIndex());
          return fullIndex;
      }
    
      bool CheckCollSystemAndEnergy(string s){ return _collSystemAndEnergy.compare(s) == 0 ? true : false; }
      bool CheckCentrality(double cent){ return (cent>_centrality.first && cent<_centrality.second) ? true : false; }
      bool CheckTriggerRange(double tpt){ return (tpt>_triggerRange.first && tpt<_triggerRange.second) ? true : false; }
      bool CheckAssociatedRange(double apt){ return (apt>_associatedRange.first && apt<_associatedRange.second) ? true : false; }
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

  class PHENIX_2010_I857187 : public Analysis {
  
    public:

     /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2010_I857187);


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
        if (( p.hasAncestorWith(Cuts::pid == 310) || p.hasAncestorWith(Cuts::pid == -310)  ||     // K0s
          p.hasAncestorWith(Cuts::pid == 130)  || p.hasAncestorWith(Cuts::pid == -130)  ||     // K0l
          p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) ||     // Xi0
          p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) ||     // Lambda
          p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) ||     // Sigma+/-
          p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) ||     // Xi-/+
          p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334) ))    // Omega-/+
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
    
    double GetXE(Particle pTrig, Particle pAssoc)
    {
        double xE = -(pAssoc.pT()/pTrig.pT())*cos(pTrig.phi()-pAssoc.phi());
                
        return xE;
    }
    
    bool FindHistoXE(YODA::Histo1D hist, double xE, string &histoID)
    {
                
        for(auto &bin : hist.bins())
        {
            if(xE < bin.xMin() || xE > bin.xMax()) continue;
            else
            {
                histoID += to_string(bin.xMin()) + "_" + to_string(bin.xMax());
                return true;
            }
        }
        
        return false;
        
    }
    
    //vmin is included in the integral, but vmax is not. Range = [vmin, vmax[
    //Give the min and max values to calculate the integral
    double GetYieldInUserRange(YODA::Histo1D& hist, double vmin, double vmax, double &n)
    {        
        double integral = 0.;
        double entries = 0.;
        
        if(vmin < hist.bin(0).xMin() || vmax > hist.bin((int)hist.numBins()-1).xMax())
        {
            MSG_ERROR("Out of range!");
            return 0.;
        }
                
        int bmin = hist.indexAt(vmin);
        int bmax = hist.indexAt(vmax);
        if(bmax < 0) bmax = (int)hist.numBins()-1;
        
        for(int i = bmin; i <= bmax; i++)
        {
            integral += hist.bin(i).sumW();
            entries += hist.bin(i).numEntries();
        }
        
        n = entries;
        
        return integral;
        
    }
    
    double GetYieldInUserRangeZYAM(YODA::Histo1D& hist, double vmin, double vmax, double &n)
    {        
        double integral = 0.;
        double entries = 0.;
        
        if(vmin < hist.bin(0).xMin() || vmax > hist.bin((int)hist.numBins()-1).xMax())
        {
            MSG_ERROR("Out of range!");
            return 0.;
        }
                
        int bmin = hist.indexAt(vmin);
        int bmax = hist.indexAt(vmax);
        if(bmax < 0) bmax = (int)hist.numBins()-1;
        
        double nbins = 0.;
        
        for(int i = bmin; i <= bmax; i++)
        {
            integral += hist.bin(i).sumW();
            entries += hist.bin(i).numEntries();
            nbins += 1.;
        }
        
        n = entries;
        
        double minValue = sqrt(-2);
        
        for(auto &bin : hist.bins())
        {
            if(std::isnan(minValue)) minValue = bin.sumW();
            if(bin.sumW() < minValue) minValue = bin.sumW();
        }
        
        if(std::isnan(minValue))
        {
            MSG_ERROR("Not possible to apply ZYAM! Returning integral without underlying event subtraction!");
            return integral;
        }
                
        integral -= nbins*minValue;
        
        return integral;
        
    }
    
    double GetPout(Particle pTrig, Particle pAssoc)
    {
        double pout = (pAssoc.pT()/GeV)*sin(pTrig.phi()-pAssoc.phi());
                
        return abs(pout);
    }

      /// Book histograms and initialise projections before the run
      void init() {

        // Initialise and register projections

        // the basic final-state projection: all final-state particles within the given eta acceptance
        
        const ChargedFinalState cfs(Cuts::abseta < 0.35);
        declare(cfs, "CFS");
        
        const PrimaryParticles pp({111, -111}, Cuts::abseta < 0.35);
        declare(pp, "PP");
        
        const PromptFinalState pfs(Cuts::abseta < 0.35 && Cuts::pid == 22);
        declare(pfs, "PFS");
        
        // the basic final-state projection: all final-state photon within the given eta acceptance


        // Declare centrality projection
        //declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      //==================================================
      // Create one correlator for each set of Collisions System / Beam Energy / Centrality Interval / Trigger pT interval / Associated pT interval
      // The number used in the constructor is the y-axis index of the histograms (d00-x00-y00) of the paper
      // Ex.: Correlator c1(1); -> is the correlator for histograms _h["0311"], _h["0411"], etc
      // Ex.: Correlator c2(2); -> is the correlator for histograms _h["0312"], _h["0412"], etc
      //==================================================
      
      //Correlators
      Correlator c1(1,1);
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(5., 7.);
      //c1.SetAssociatedRange(0., 999.);
      c1.SetPID({111, -111});
      Correlators.push_back(c1);
      
      Correlator c2(1,2);
      c2.SetCollSystemAndEnergy("pp200GeV");
      c2.SetCentrality(0., 100.);
      c2.SetTriggerRange(5., 7.);
      //c2.SetAssociatedRange(0., 7.);
      c2.SetPID({22});
      Correlators.push_back(c2);
      
      Correlator c3(2,1);
      c3.SetCollSystemAndEnergy("pp200GeV");
      c3.SetCentrality(0., 100.);
      c3.SetTriggerRange(7., 9.);
      //c3.SetAssociatedRange(0., 999.);
      c3.SetPID({111, -111});
      Correlators.push_back(c3);
      
      Correlator c4(2,2);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetCentrality(0., 100.);
      c4.SetTriggerRange(7., 9.);
      //c4.SetAssociatedRange(0., 9.);
      c4.SetPID({22});
      Correlators.push_back(c4);
      
      Correlator c5(3,1);
      c5.SetCollSystemAndEnergy("pp200GeV");
      c5.SetCentrality(0., 100.);
      c5.SetTriggerRange(9., 12.);
      //c5.SetAssociatedRange(0., 999.);
      c5.SetPID({111, -111});
      Correlators.push_back(c5);
      
      Correlator c6(3,2);
      c6.SetCollSystemAndEnergy("pp200GeV");
      c6.SetCentrality(0., 100.);
      c6.SetTriggerRange(9., 12.);
      //c6.SetAssociatedRange(0., 12.);
      c6.SetPID({22});
      Correlators.push_back(c6);
      
      Correlator c7(4,1);
      c7.SetCollSystemAndEnergy("pp200GeV");
      c7.SetCentrality(0., 100.);
      c7.SetTriggerRange(12., 15.);
      //c7.SetAssociatedRange(0., 999.);
      c7.SetPID({111, -111});
      Correlators.push_back(c7);
      
      Correlator c8(4,2);
      c8.SetCollSystemAndEnergy("pp200GeV");
      c8.SetCentrality(0., 100.);
      c8.SetTriggerRange(12., 15.);
      //c8.SetAssociatedRange(0., 15.);
      c8.SetPID({22});
      Correlators.push_back(c8);

      for(Correlator& corr : Correlators)
      {
          book(sow[corr.GetFullIndex()],"sow" + corr.GetFullIndex());
          book(_h["0" + to_string(corr.GetIndex()) + "1" + to_string(corr.GetSubIndex())], corr.GetIndex(), 1, corr.GetSubIndex());
          book(_h["0" + to_string(corr.GetIndex()+4) + "1" + to_string(corr.GetSubIndex())], corr.GetIndex()+4, 1, corr.GetSubIndex());
          
          string refname = mkAxisCode(corr.GetIndex(), 1, corr.GetSubIndex());
          const Estimate1D& refdata = refData(refname);
          for(auto &bin : refdata.bins())
          {
              book(_DeltaPhixE["0" + to_string(corr.GetIndex()) + "1" + to_string(corr.GetSubIndex()) + "xE_" + to_string(bin.xMin()) + "_" + to_string(bin.xMax())], "DeltaPhi_0" + to_string(corr.GetIndex()) + "1" + to_string(corr.GetSubIndex()) + "xE_" + to_string(bin.xMin()) + "_" + to_string(bin.xMax()), 24, 0, M_PI);
          }
          nTriggers[corr.GetFullIndex()] = 0;
      }
      
    } // End of init

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
          CollSystem = "pp";
          nNucleons = 1.;
      //if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000010020) CollSystem = "dAu";
      //if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000290630) CollSystem = "dAu";
      
      string cmsEnergy = "Empty";
      if (fuzzyEquals(sqrtS()/GeV, 200*nNucleons, 5)) cmsEnergy = "200GeV";
      
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
    
    // loop over charged final state particles - PI0
    for(const Particle& pTrig : ppTrigPi0.particles())
    {
        //Check if is secondary
        if(isSecondary(pTrig)) continue;
          
        
        for(Correlator& corr : Correlators)
        {
            if(!corr.CheckPID({111, -111})) continue;
            if(!corr.CheckTriggerRange(pTrig.pT()/GeV)) continue;  
            nTriggers[corr.GetFullIndex()]++;
        }

        // Hadron loop
        for(const Particle& pAssoc : cfs.particles())
        {
                
            //Check if Trigger and Associated are the same particle
            if(isSameParticle(pTrig,pAssoc)) continue;
                
            //Check if is secondary
            if(isSecondary(pAssoc)) continue;

            //https://rivet.hepforge.org/code/dev/structRivet_1_1DeltaPhiInRange.html
            double dPhi = deltaPhi(pTrig, pAssoc, true);//this does NOT rotate the delta phi to be in a given range
                        
            double xE = GetXE(pTrig,pAssoc);
            
            for(Correlator& corr : Correlators)
            {
                if(!corr.CheckPID({111, -111})) continue;
                
                if(!corr.CheckTriggerRange(pTrig.pT()/GeV)) continue;
                
                if(xE > 0.)
                {
                    string histoID = "xE_";
                    if(!FindHistoXE(*_h["0" + to_string(corr.GetIndex()) + "1" + to_string(corr.GetSubIndex())], xE, histoID)) continue;
                    _DeltaPhixE["0" + to_string(corr.GetIndex()) + "1" + to_string(corr.GetSubIndex()) + histoID]->fill(abs(dPhi));
                }
                
                if(pAssoc.pT()/GeV > 2. && pAssoc.pT()/GeV < 10.)
                {
                    double pout = GetPout(pTrig,pAssoc);
                    _h["0" + to_string(corr.GetIndex()+4) + "1" + to_string(corr.GetSubIndex())]->fill(pout);
                }
                
            } //end of correlators loop 
                
        } // end of loop over associated particles

    } // particle loop
    
    
    // loop over charged final state particles - PHOTONS
    for(const Particle& pTrig : pfsTrigPhotons.particles())
    {
        //Check if is secondary
        if(isSecondary(pTrig)) continue;
          
        
        for(Correlator& corr : Correlators)
        {
            if(!corr.CheckPID({22})) continue;
            if(!corr.CheckTriggerRange(pTrig.pT()/GeV)) continue;  
            nTriggers[corr.GetFullIndex()]++;
        }
        
        // Hadron loop
        for(const Particle& pAssoc : cfs.particles())
        {
                
            //Check if Trigger and Associated are the same particle
            if(isSameParticle(pTrig,pAssoc)) continue;
                
            //Check if is secondary
            if(isSecondary(pAssoc)) continue;

            //https://rivet.hepforge.org/code/dev/structRivet_1_1DeltaPhiInRange.html
            double dPhi = deltaPhi(pTrig, pAssoc, true);//this does NOT rotate the delta phi to be in a given range
                        
            double xE = GetXE(pTrig,pAssoc);
            
            
            
            for(Correlator& corr : Correlators)
            {
                if(!corr.CheckPID({22})) continue;
                
                if(!corr.CheckTriggerRange(pTrig.pT()/GeV)) continue;
                
                if(xE > 0.)
                {
                    //cout << "xE: " << xE << endl;
                    //cout << "DeltaPhi: " << dPhi << endl;
                    string histoID = "xE_";
                    if(!FindHistoXE(*_h["0" + to_string(corr.GetIndex()) + "1" + to_string(corr.GetSubIndex())], xE, histoID)) continue;
                    _DeltaPhixE["0" + to_string(corr.GetIndex()) + "1" + to_string(corr.GetSubIndex()) + histoID]->fill(abs(dPhi));
                }
                
                
                if(pAssoc.pT()/GeV > 2. && pAssoc.pT()/GeV < 10.)
                {
                    double pout = GetPout(pTrig,pAssoc);
                    _h["0" + to_string(corr.GetIndex()+4) + "1" + to_string(corr.GetSubIndex())]->fill(pout);
                        
                }
                
                
            } //end of correlators loop 
                
        } // end of loop over associated particles

    } // particle loop
    
    
    
    
    
  }
    /// Normalise histograms etc., after the run
    void finalize() {
        
        for(Correlator& corr : Correlators)
        {
            string refname = mkAxisCode(corr.GetIndex(), 1, corr.GetSubIndex());
            const Estimate1D& refdata = refData(refname);
            for(auto &bin : refdata.bins())
            {
                _DeltaPhixE["0" + to_string(corr.GetIndex()) + "1" + to_string(corr.GetSubIndex()) + "xE_" + to_string(bin.xMin()) + "_" + to_string(bin.xMax())]->scaleW(sow[corr.GetFullIndex()]->numEntries()/(nTriggers[corr.GetFullIndex()]*sow[corr.GetFullIndex()]->sumW()));
                double entries = 0.;
                double yield = GetYieldInUserRangeZYAM(*_DeltaPhixE["0" + to_string(corr.GetIndex()) + "1" + to_string(corr.GetSubIndex()) + "xE_" + to_string(bin.xMin()) + "_" + to_string(bin.xMax())], M_PI/2., M_PI, entries);
                string hname = "0" + to_string(corr.GetIndex()) + "1" + to_string(corr.GetSubIndex());
                _h[hname]->fill(bin.xMid(),yield/entries);
                    
            }
            
            _h["0" + to_string(corr.GetIndex()+4) + "1" + to_string(corr.GetSubIndex())]->scaleW(sow[corr.GetFullIndex()]->numEntries()/(nTriggers[corr.GetFullIndex()]*sow[corr.GetFullIndex()]->sumW()));
            
        }
        

    }


    //Histograms and variables

    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> sow;
    map<string, Histo1DPtr> _DeltaPhixE;
    map<int, Histo1DPtr> _DeltaPhiSub;
    map<string, int> nTriggers;
    vector<Correlator> Correlators;


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(PHENIX_2010_I857187);
  }
