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
      string _collSystemAndEnergy;
      pair<double,double> _centrality;
      pair<double,double> _triggerRange;
      pair<double,double> _associatedRange;
      vector<int> _pid;
      bool _NoPTassociated = false;

    public:
    
      /// Constructor
      Correlator(int index, int subindex) {
        _index = index;
        _subindex = subindex;
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
      string GetFullIndex()
      {
          string fullIndex = to_string(GetIndex()) + to_string(GetSubIndex());
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

double GetDeltaPhi(Particle pAssoc, Particle pTrig)
    {
        //https://rivet.hepforge.org/code/dev/structRivet_1_1DeltaPhiInRange.html
        double dPhi = deltaPhi(pTrig, pAssoc, true);//this does NOT rotate the delta phi to be in a given range
                    
        if(dPhi < -M_PI/2.)
        {
            dPhi += 2.*M_PI;
        }
        else if(dPhi > 3.*M_PI/2.)
        {
            dPhi -= 2*M_PI;
        }
                    
        return dPhi;   
      }
      
      Histo1DPtr SubtractBackgroundZYAM(Histo1DPtr histo)
      {
        
        YODA::Histo1D hist = *histo;
        
        double minValue = sqrt(-2);
        double binWidth = 0.;
        int minValueEntries = 0.;
        
        for(auto &bin : hist.bins())
        {
            if(std::isnan(minValue))
            {
                minValue = bin.sumW();
                binWidth = bin.width();
                minValueEntries = bin.numEntries();
            }
            if(bin.sumW()/bin.width() < minValue/binWidth)
            {
                minValue = bin.sumW();
                binWidth = bin.width();
                minValueEntries = bin.numEntries();
            }
        }
                
        hist.reset();
        
        for(auto &bin : hist.bins())
        {
            bin.fillBin((minValue*bin.width())/(minValueEntries*binWidth), minValueEntries);
        }
        
        *histo = YODA::subtract(*histo, hist);
        
        return histo;
                
      }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

        const ChargedFinalState cfs(Cuts::abseta < 1.0);
        declare(cfs, "CFS");

        const ChargedFinalState cfsEta(Cuts::abseta < 0.7);
        declare(cfsEta, "CFSETA");
        
        beamOpt = getOption<string>("beam", "NONE");

        if (beamOpt == "DAU") collSys = dAu;
        else if (beamOpt == "AUAU") collSys = AuAu;
     
        declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");
       
        
      //===========================================================
      //===========================================================
      //                   Figure 2   : FIX ME: add in the nTriggers
        //   for the raw und eta. book (sow) for all but sub.
      //===========================================================
      //===========================================================
      for (int ptt = 0; ptt < numTrigPtBins - 1; ptt++)
      {
        for (int cb = 0; cb < numCentBins - 1; cb++)
        {
            Correlator c1 (ptt,cb);
            c1.SetCollSystemAndEnergy("AuAu200GeV");
            c1.SetCentrality(CentBins[cb], CentBins[cb+1]);
            c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt + 1]);
            c1.SetNoPTassociated();
            Correlators.push_back(c1);
        }
      }

      for(Correlator corr : Correlators)
      {
          if(corr.GetIndex() <= 1)
          {
              //raw |eta| < 1
              string name_raw = "raw_d" + to_string((corr.GetIndex()*2)+1) + "x1y" + to_string((corr.GetSubIndex()*2)+1);
              book(_h[name_raw], (corr.GetIndex()*2)+1, 1, (corr.GetSubIndex()*2)+1);
              nTriggers[name_raw]=0;
              book(sow[name_raw], "sow" + name_raw);

              //limited eta acceptance |eta| < 0.7
              string name_eta = "eta_d" + to_string((corr.GetIndex()*2)+1) + "x1y" + to_string((corr.GetSubIndex()*2)+2);
              book(_h[name_eta], (corr.GetIndex()*2)+1, 1, (corr.GetSubIndex()*2)+2);
              nTriggers[name_eta]=0;
              book(sow[name_eta], "sow" + name_eta);

              //Background subtracted |eta| < 1
              string name_sub = "sub_d" + to_string((corr.GetIndex()*2)+2) + "x1y" + to_string(corr.GetSubIndex()+1);
              book(_h[name_sub], (corr.GetIndex()*2)+2, 1, corr.GetSubIndex()+1);

          }
        
        if(corr.GetIndex() == 2)
        {
            if(corr.GetSubIndex() <= 2)
            {
                //raw |eta| < 1
                string name_raw = "raw_d" + to_string((corr.GetIndex()*2)+1) + "x1y" + to_string((corr.GetSubIndex()*2)+1);
                book(_h[name_raw], (corr.GetIndex()*2)+1, 1, (corr.GetSubIndex()*2)+1);
                nTriggers[name_raw]=0;
                book(sow[name_raw], "sow" + name_raw);
                //limited eta acceptance |eta| < 0.7
                string name_eta = "eta_d" + to_string((corr.GetIndex()*2)+1) + "x1y" + to_string((corr.GetSubIndex()*2)+2);
                book(_h[name_eta], (corr.GetIndex()*2)+1, 1, (corr.GetSubIndex()*2)+2);
                nTriggers[name_eta]=0;
                book(sow[name_eta], "sow" + name_eta);
            }
            else
            {
                //raw |eta| < 1
                string name_raw = "raw_d" + to_string((corr.GetIndex()*2)+2) + "x1y1";
                book(_h[name_raw], (corr.GetIndex()*2)+2, 1, 1);
                nTriggers[name_raw]=0;
                book(sow[name_raw], "sow" + name_raw);
                //limited eta acceptance |eta| < 0.7
                string name_eta = "eta_d" + to_string((corr.GetIndex()*2)+2) + "x1y1";
                book(_h[name_eta], (corr.GetIndex()*2)+2, 1, 1);
                nTriggers[name_eta]=0;
                book(sow[name_eta], "sow" + name_eta);

            }
            
            //Background subtracted |eta| < 1
            string name_sub = "sub_d" + to_string((corr.GetIndex()*2)+3) + "x1y" + to_string(corr.GetSubIndex()+1);
            book(_h[name_sub], (corr.GetIndex()*2)+3, 1, corr.GetSubIndex()+1);
            
        }
        
        
        if(corr.GetIndex() == 3)
        {
            //raw |eta| < 1
            string name_raw = "raw_d" + to_string((corr.GetIndex()*2)+2) + "x1y" + to_string((corr.GetSubIndex()*2)+1);
            book(_h[name_raw], (corr.GetIndex()*2)+2, 1, (corr.GetSubIndex()*2)+1);
            nTriggers[name_raw]=0;
            book(sow[name_raw], "sow" + name_raw);
            //limited eta acceptance |eta| < 0.7
            string name_eta = "eta_d" + to_string((corr.GetIndex()*2)+2) + "x1y" + to_string((corr.GetSubIndex()*2)+2);
            book(_h[name_eta], (corr.GetIndex()*2)+2, 1, (corr.GetSubIndex()*2)+2);
            nTriggers[name_eta]=0;
            book(sow[name_eta], "sow" + name_eta);
            //Background subtracted |eta| < 1
            string name_sub = "sub_d" + to_string((corr.GetIndex()*2)+3) + "x1y" + to_string(corr.GetSubIndex()+1);
            book(_h[name_sub], (corr.GetIndex()*2)+3, 1, corr.GetSubIndex()+1);
        }
      }

    
    //YODA FILE EXPLANATION
    //d01,x01,y01 is raw 0-12%
    //d01,x01,y02 is deta 0-12%
    //d01,x01,y03 is raw 20-40% etc.
    //Only for figure 2 and 3 when the title is not background subtracted.

  

      //===================================================================
      //===================================================================
      //Figure 3 FIX ME:
      // 
      //===================================================================
      //===================================================================
      for (int ptt = 0; ptt < numTrigPtBins; ptt++)
      {
        for (int pta = 1; pta < numAssocPtBins; pta++)
        {
            Correlator c2 (ptt,pta);
            c2.SetCollSystemAndEnergy("AuAu200GeV");
            c2.SetCentrality(pTAssocBins[pta], pTAssocBins[pta + 1]);
            c2.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt + 1]);
            c2.SetAssociatedRange(1,6);
            Correlators3.push_back(c2);
        }
      }

    for(Correlator corr : Correlators3)
      {
     	/// Yoda d10-x01-y01 to d10-x01-y16: all Bksub for Fig 3. 1st 4 are pTassoc .5-1
         // 2nd 4 are pTassoc 1-1.5, 3rd 4 are pTassoc 1.5-2.5, last 4 are 2.5-4 pTassoc. However
         // 1,2,3,4 of each set are different pTTrig.
     	string name_sub = "sub_d10x01y" + to_string(corr.GetIndex() + ((corr.GetSubIndex()-1)*4)+1);
        book(_h[name_sub], 10, 1, corr.GetIndex() + ((corr.GetSubIndex() - 1) * 4) + 1);
        //cout << "name_sub: " << name_sub << endl; //Debugging

        //Yoda d11-x01-y01 to d11-x01-y16: all AuAu Raw for Fig 3. Same logic as above
        string name_AuAuRaw = "AuAuRaw_d11x1y" + to_string(corr.GetIndex() + ((corr.GetSubIndex()-1)*4)+1);
      	book(_h[name_AuAuRaw], 11, 1, corr.GetIndex() + ((corr.GetSubIndex()-1)*4)+1);
      	nTriggers[name_AuAuRaw]=0;
        book(sow[name_AuAuRaw], "sow" + name_AuAuRaw);

    	// Yoda d12-x01-y01 to d12-x01-y16: all dAu for Fig 3. Same logic as above
         string name_dAu = "dAu_d12x1y" + to_string(corr.GetIndex() + ((corr.GetSubIndex()-1)*4)+1);
         book(_h[name_dAu], 12, 1, corr.GetIndex() + ((corr.GetSubIndex()-1)*4)+1);
         nTriggers[name_dAu]=0;
         book(sow[name_dAu], "sow" + name_dAu);

        // Yoda d13-x01-y01 to d13-x01-y16: same logic
         string name_dEta = "dEta_d13x1y" + to_string(corr.GetIndex() + ((corr.GetSubIndex()-1)*4)+1);
         book(_h[name_dEta], 13, 1, corr.GetIndex() + ((corr.GetSubIndex()-1)*4)+1);
         nTriggers[name_dEta]=0;
         book(sow[name_dEta], "sow" + name_dEta);
 
      }

     ////////////////////////////////////////////////////////////////////////////////////////////////////
     ////////////////////////////////////////////////////////////////////////////////////////////////////
     //                                 Figure 6                                                       //
     ////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////
      for (int ptt = 0; ptt < numTrigPtBins - 1; ptt++)
      {
      	for (int cb = 0; cb < 1; cb++)
		{
			Correlator c6 (ptt,cb);
			c6.SetCollSystemAndEnergy("dAu200GeV");
			c6.SetCentrality(CentBins[cb], CentBins[cb+1]);
			c6.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
			c6.SetNoPTassociated();
			Correlators6.push_back(c6);
		}
      }



     ////////////////////////////////////////////////////////////////////////////////////////////////////
     ////////////////////////////////////////////////////////////////////////////////////////////////////
     //                                  Figure 7    : add in book(sow)                                                  //
     ////////////////////////////////////////////////////////////////////////////////////////////////////
     ////////////////////////////////////////////////////////////////////////////////////////////////////

  	for (int ptt = 1; ptt < numTrigPtBins; ptt++)
      {
        for (int cb = 0; cb < 1; cb++)
        {
            Correlator c7 (ptt,cb);
            c7.SetCollSystemAndEnergy("AuAu200GeV");
            c7.SetCentrality(CentBins[cb], CentBins[cb+1]);
            c7.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt + 1]);
            c7.SetNoPTassociated();
            Correlators7.push_back(c7);
        }
      }

  for(Correlator corr : Correlators7)
      {
        if (corr.GetIndex() == 1)
        {
        	string name_dPhi = "dPhi_d22x1y1";
        	book(_h[name_dPhi], 22 , 1 , 1);
        	nTriggers[name_dPhi]=0;
        	book(sow[name_dPhi], "sow" + name_dPhi);
        }
        else if (corr.GetIndex() == 2)
        {
       		string name_dPhi2 = "dPhi2_d22x1y2";
        	book(_h[name_dPhi2], 22 , 1 , 2);
        	nTriggers[name_dPhi2]=0;
        	book(sow[name_dPhi2], "sow" + name_dPhi2);
        }

        else if (corr.GetIndex() == 3)
        {
        	string name_dPhi3 = "dPhi3_d22x1y1";
        	book(_h[name_dPhi3], 23 , 1 , 1);
        	nTriggers[name_dPhi3]=0;
       		book(sow[name_dPhi3], "sow" + name_dPhi3);
        }
      }
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     //                                  Figure 8 : FIX ME ANTONIO PLEASE                                                  //
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////////////////////////////////////

      for (int ptt = 1; ptt < numTrigPtBins - 1; ptt++)
      {
        for (int cb = 0; cb < numCentBins - 1; cb++)
        {
            Correlator c8 (ptt,cb);
            c8.SetCollSystemAndEnergy("AuAu200GeV");
            c8.SetCentrality(CentBins[cb], CentBins[cb+1]);
            c8.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
            c8.SetNoPTassociated();
            Correlators8.push_back(c8);
        }
      }

      for (Correlator corr : Correlators8)
      {
        string name_AuAu = "AuAu_d24x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        book(_h[name_AuAu], 24, 1, 1);
        book(_h[name_AuAu], 24, 1, 2);

        string name_AuAu2 = "AuAu2_d25x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        book(_h[name_AuAu2], 25, 1, 1);

        string name_AuAu3 = "AuAu3_d26x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        book(_h[name_AuAu3], 26, 1, 1);
        book(_h[name_AuAu3], 26, 1, 2);

        string name_AuAu4 = "AuAu4_d27x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        book(_h[name_AuAu4], 27,1, 1);

        string name_AuAu5 = "AuAu5_d28x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        book(_h[name_AuAu5], 28,1, 1);
        book(_h[name_AuAu5], 28,1, 2);

        string name_AuAu6 = "AuAu5_d29x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        book(_h[name_AuAu6], 29,1, 1);

        string name_dAuAu = "dAuAu_d30x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        book(_h[name_dAuAu], 30,1, 1);

      }

    } //ends the init



    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
      const double c = cent();

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      
      //==================================================
      // Select the histograms accordingly to the collision system, beam energy and centrality
      // WARNING: Still not implemented for d-Au
      //==================================================
     
      string SysAndEnergy = "";

      if (collSys == dAu) SysAndEnergy = "dAu200GeV";
      else if (collSys == AuAu) SysAndEnergy = "AuAu200GeV";


      //SysAndEnergy = CollSystem + cmsEnergy;

      double triggerptMin = 999.;
      double triggerptMax = -999.;
      double associatedptMin = 999.;
      double associatedptMax = -999.;
    
      bool isVeto = true;

      //INCOMPLETE: Fill Histograms for figure 2
      for(Correlator corr : Correlators)
      {
          if(!corr.CheckCentrality(c)) continue;
          if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
          
          if(corr.GetIndex() <=1)
          {
          string name_raw = "raw_d" + to_string((corr.GetIndex()*2)+1) + "x1y" + to_string((corr.GetSubIndex()*2)+1);
          sow[name_raw]->fill();

          string name_eta = "eta_d" + to_string((corr.GetIndex()*2)+1) + "x1y" + to_string((corr.GetSubIndex()*2)+2);
          sow[name_eta]->fill();

          string name_sub = "sub_d" + to_string((corr.GetIndex()*2)+2) + "x1y" + to_string(corr.GetSubIndex()+1);
          sow[name_sub]->fill();
      	  }

      	  if (corr.GetIndex() == 2)
      	  {
      	  	if (corr.GetSubIndex() <= 2)
      	  	{
      	  		string name_raw = "raw_d" + to_string((corr.GetIndex()*2)+1) + "x1y" + to_string((corr.GetSubIndex()*2)+1);
      	  		sow[name_raw]->fill();

      	  		string name_eta = "eta_d" + to_string((corr.GetIndex()*2)+1) + "x1y" + to_string((corr.GetSubIndex()*2)+2);
      	  		sow[name_eta]->fill();
      	  	}
      	  	else
      	  	{
      	  		string name_raw = "raw_d" + to_string((corr.GetIndex()*2)+2) + "x1y1";
      	  		sow[name_raw]->fill();

      	  		string name_eta = "eta_d" + to_string((corr.GetIndex()*2)+2) + "x1y1";
      	  		sow[name_raw]->fill();
      	  	}

      	  	string name_sub = "sub_d" + to_string((corr.GetIndex()*2)+3) + "x1y" + to_string(corr.GetSubIndex()+1);
      	  	sow[name_sub]->fill();
      	  }

      	  if (corr.GetIndex() == 3)
      	  {
      	  	string name_raw = "raw_d" + to_string((corr.GetIndex()*2)+2) + "x1y" + to_string((corr.GetSubIndex()*2)+1);
      	  	sow[name_raw]->fill();

      	  	string name_eta = "eta_d" + to_string((corr.GetIndex()*2)+2) + "x1y" + to_string((corr.GetSubIndex()*2)+2);
      	  	sow[name_eta]->fill();

      	  	string name_sub = "sub_d" + to_string((corr.GetIndex()*2)+3) + "x1y" + to_string(corr.GetSubIndex()+1);
      	  	sow[name_sub]->fill();
      	  }
      }

      //Fill Histograms for figure 3
      for (Correlator corr : Correlators3)
      {
        
        if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;

        string name_sub = "sub_d10x01y" + to_string(corr.GetIndex() + ((corr.GetSubIndex()-1)*4)+1); 
        sow[name_sub]->fill();

        string name_AuAuRaw = "AuAuRaw_d11x1y" + to_string(corr.GetIndex() + ((corr.GetSubIndex()-1)*4)+1);
        sow[name_AuAuRaw]->fill();

        string name_dAu = "dAu_d12x1y" + to_string(corr.GetIndex() + ((corr.GetSubIndex()-1)*4)+1);
        sow[name_dAu]->fill();

        string name_dEta = "dEta_d13x1y" + to_string(corr.GetIndex() + ((corr.GetSubIndex()-1)*4)+1);
        sow[name_dEta]->fill();

      }

      //Fill Histograms for figure 6

      //Fill Histograms for figure 7
      for (Correlator corr : Correlators7)
      {
        if (!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;

        string name_dPhi = "dPhi_d22x1y1";
        sow[name_dPhi]->fill();

        string name_dPhi2 = "dPhi2_d22x1y2";
        sow[name_dPhi2]->fill();

        string name_dPhi3 = "dPhi3_d22x1y1";
        sow[name_dPhi3]->fill();
      }

      //Fill Histograms for figure 8; FIX ME TO THE FIG 8 LOGIC
      for (Correlator corr : Correlators8)
      {
        if (!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;

        string name_AuAu = "AuAu_d24x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        sow[name_AuAu]->fill();

        string name_AuAu2 = "AuAu2_d25x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        sow[name_AuAu2]->fill();

        string name_AuAu3 = "AuAu3_d26x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        sow[name_AuAu3]->fill();

        string name_AuAu4 = "AuAu4_d27x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        sow[name_AuAu4]->fill();

        string name_AuAu5 = "AuAu5_d28x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        sow[name_AuAu5]->fill();

        string name_AuAu6 = "AuAu5_d29x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        sow[name_AuAu6]->fill();

        string name_dAuAu = "dAuAu_d30x1y" + to_string(corr.GetIndex() + corr.GetSubIndex());
        sow[name_dAuAu]->fill();
      }

Particles chargedParticles = cfs.particles();
        
        for(Particle pTrig : chargedParticles)
        {
            for(Correlator corr : Correlators)
            {
                if(!corr.CheckCentrality(c)) continue;
                if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
                
                string name_raw = "raw_d" + to_string((corr.GetIndex()*2)+1) + "x1y" + to_string((corr.GetSubIndex()*2)+1);
                nTriggers[name_raw]++;
                
            }
            
            for(Particle pAssoc : chargedParticles)
            {
                if(isSameParticle(pTrig, pAssoc)) continue;
                
                double DeltaPhi = GetDeltaPhi(pAssoc, pTrig);
            
                for(Correlator corr : Correlators)
                {
                    if(!corr.CheckCentrality(c)) continue;
                    if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                    if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
                    
                    if(corr.GetIndex() <= 1)
                    {
                        string name_raw = "raw_d" + to_string((corr.GetIndex()*2)+1) + "x1y" + to_string((corr.GetSubIndex()*2)+1);
                        _h[name_raw]->fill(DeltaPhi);
                    }
                    else if(corr.GetIndex() == 2)
                    {
                        
                        
                    }
                    
        
                }
            
            }
        
        
        }
    
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
      bool AuAu200_available = false;
      bool dAu_available = false;

      for (auto element : _h)
      {
        string name = element.second->name();
        if (name.find("AuAu") != std::string::npos)
        {
          if (element.second->numEntries()>0) AuAu200_available=true;
          else
          {
            AuAu200_available=false;
            break;
          }

        }
         else if (name.find("dAu") != std::string::npos)
        {
          if (element.second->numEntries()>0) dAu_available=true;
          else
          {
            dAu_available=false;
            break;
          }
          
        }
      }
      
      //if((!AuAu200_available) || (!dAu_available)) return;
      
      for(Correlator corr : Correlators) //Finish with rest of logic
      {
          string name_raw = "raw_d" + to_string((corr.GetIndex()*2)+1) + "x1y" + to_string((corr.GetSubIndex()*2)+1);
          _h[name_raw]->scaleW(sow[name_raw]->numEntries()/(nTriggers[name_raw]*sow[name_raw]->sumW()));
          
          string name_sub = "sub_d" + to_string((corr.GetIndex()*2)+3) + "x1y" + to_string(corr.GetSubIndex()+1);
          _h[name_sub] = SubtractBackgroundZYAM(_h[name_raw]);
      } 
  }


    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> sow;
    map<string, Histo1DPtr> _DeltaPhixE;
    map<int, Histo1DPtr> _DeltaPhiSub;
    map<string, int> nTriggers;

    vector<Correlator> Correlators;
    vector<Correlator> Correlators3;
    vector<Correlator> Correlators6;
    vector<Correlator> Correlators7;
    vector<Correlator> Correlators8;
    
    string beamOpt;
    enum CollisionSystem {dAu, AuAu};
    CollisionSystem collSys;
    //@}


  };


  DECLARE_RIVET_PLUGIN(STAR_2010_I851937);

}
