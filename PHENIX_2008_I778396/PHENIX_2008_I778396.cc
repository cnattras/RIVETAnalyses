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
static const float pTTrigBins[] = {2.0,3.0,4.0,5.0,10.0};
static const int numAssocPtBins = 6;
static const float pTAssocBins[] = {0.4,1.0,2.0,3.0,4.0,5.0,10.0};
static const int numCentBins = 5;
static const float CentBins[] = {0.0,12.0,20.0,40.0,60.0,92.0};
static const int numDeltaPhiBins = 26;
static const float DeltaPhiBins[] = {-1.5048,-1.275,-0.98,-0.685,-0.415,-0.22,-0.1,0,0.1,0.22,0.415,0.685,0.98,
1.275,1.595,1.965,2.305,2.55,2.75,2.945,3.14,3.44,3.535,3.73,3.935,4.32,4.7012};
static const int numpTAssocBins2 = 4;
static const float pTAssocBins2[] = {1.0,1.5,2.0,2.5,3.0};
//For figure 12
static const int numCentBins12 = 4;
static const float CentBins12[] = {0.0,20.0,40.0,60.0,92.0};
static const int numpTAssocBins12 = 8;
static const float pTAssocBins12[] = {0.345,0.915,1.415,1.915,2.415,3.055,4.095,5.285,5.88};
//For figure 25 
static const int numCentBins25 = 10;
static const float CentBins25[] = {0.0,5.0,10.0,20.0,30.0,40.04,50.0,60.0,70.0,92.0,100};

//For figure 26 
static const int numCentBins26 = 5;
static const float CentBins26[] = {0.0,20.0,40.0,60.0,70.0,92.0};

static const int numCentBins31 = 4;
static const float CentBins31[] = {0,20,40,60,92};
static const float Ncoll[] = {779,297,90.6,14.5};

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
      bool _noCentrality = false;
      bool _noAssoc = false;
    public:
   
      /// Constructor
      Correlator(int index, int subindex, int subsubindex) {
        _index = index;
        _subindex = subindex;
        _subsubindex = subsubindex;
      }
      void SetCollSystemAndEnergy(string s){ _collSystemAndEnergy = s; }
      void SetCentrality(double cmin, double cmax){ _centrality = make_pair(cmin, cmax); }
      void SetNoCentrality(){ _noCentrality = true; }
      void SetNoAssoc(){ _noAssoc = true; }
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
          string fullIndex = to_string(GetIndex()) + to_string(GetSubIndex())+ to_string(GetSubIndex());
          return fullIndex;
      }
   
      bool CheckCollSystemAndEnergy(string s){ return _collSystemAndEnergy.compare(s) == 0 ? true : false; }
     bool CheckCentrality(double cent){ return ((cent>_centrality.first && cent<_centrality.second) || _noCentrality == true) ? true : false; }
      bool CheckTriggerRange(double tpt){ return (tpt>_triggerRange.first && tpt<_triggerRange.second) ? true : false; }
      bool CheckAssociatedRange(double apt){ return ((apt>_associatedRange.first && apt<_associatedRange.second) || _noAssoc == true) ? true : false; }
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
  class PHENIX_2008_I778396 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2008_I778396);
    
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

    
    double getYieldRangeUser(Histo1DPtr histo, double xmin, double xmax, double &fraction)
    {
        //This will include bins partially covered by the user range
        
        YODA::Histo1D hist = *histo;
        
        double integral = 0.;
        
        if(xmax < xmin) throw RangeError("Error: xmin > xmax");
        if(xmin < hist.bin(0).xMin()) throw RangeError("xmin is out of range");
        if(xmax > hist.bin(hist.numBins()-1).xMax()) throw RangeError("xmax is out of range");
        
        for(auto &bin : hist.bins())
        {
            if((bin.xMin() > xmin) && (bin.xMax() < xmax))
            {
                integral += bin.sumW();
                fraction += bin.numEntries();
            }
            else if((bin.xMin() < xmin) && (bin.xMax() > xmin))
            {
                double perc = bin.xMax() - xmin;
                integral += perc*bin.sumW();
                fraction += perc*bin.numEntries();
                
            }
            else if((bin.xMin() < xmax) && (bin.xMax() > xmax))
            {
                double perc = xmax - bin.xMin();
                integral += perc*bin.sumW();
                fraction += perc*bin.numEntries();
            }
        }
        
        return integral;
        
    }

    /*
    double getTruncatedMean(Histo1DPtr histo, double xmin, double xmax, double &fraction)
    {
        //This will include bins partially covered by the user range
        
        YODA::Histo1D hist = *histo;
        
        double integral = 0.;
        
        if(xmax < xmin) throw RangeError("Error: xmin > xmax");
        if(xmin < hist.bin(0).xMin()) throw RangeError("xmin is out of range");
        if(xmax > hist.bin(hist.numBins()-1).xMax()) throw RangeError("xmax is out of range");
        
        for(auto &bin : hist.bins())
        {
            if((bin.xMin() > xmin) && (bin.xMax() < xmax))
            {
                integral += bin.sumW()/bin.width();
                fraction += bin.numEntries()/bin.width();
            }
            else if((bin.xMin() < xmin) && (bin.xMax() > xmin))
            {
                double perc = (bin.xMax() - xmin)/bin.width();
                integral += perc*(bin.sumW()/bin.width());
                fraction += perc*(bin.numEntries()/bin.width());
                
            }
            else if((bin.xMin() < xmax) && (bin.xMax() > xmax))
            {
                double perc = (xmax - bin.xMin())/bin.width();
                integral += perc*(bin.sumW()/bin.width());
                fraction += perc*(bin.numEntries()/bin.width());
            }
        }
        
        return integral;
        
    }
    */



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
/*
    double GetDeltaEta(Particle pAssoc, Particle pTrig)
    {
           double dEta = deltaEta(pTrig, pAssoc);

          if(dEta < -M_PI/2.)
          {
            dEta += 2.*M_PI;
          }
          else if(dEta > 3.*M_PI/2.)
          {
            dEta -= 2*M_PI;
          }
                    
        return dEta;   
    }
   */
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

    void init() {
        const ChargedFinalState cfs(Cuts::abseta < 0.35 && Cuts::abscharge > 0);
        declare(cfs, "CFS");
        
        const PrimaryParticles pp(pdgPi0, Cuts::abseta < 0.35);
        declare(pp, "PP");
        
        const PromptFinalState pfs(Cuts::abseta < 0.35 && Cuts::pid == 22);
        declare(pfs, "PFS");

        beamOpt = getOption<string>("beam", "NONE");

        if (beamOpt == "PP") collSys = pp0;
        else if (beamOpt == "AUAU") collSys = AuAu;
        
        declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");


      int pta;
        int ptt;
        int pta2;
        int cb;
        int i;

//*****************************************************************************
// The following will book the histograms for Figure 36-38
for(pta = 0; pta<numAssocPtBins; pta++){
	for(ptt = 0; ptt<numTrigPtBins; ptt++){

      Correlator c1(pta,ptt,-1);
      c1.SetCollSystemAndEnergy("pp200GeV");
     	c1.SetNoCentrality();
     	c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
     	c1.SetAssociatedRange(pTAssocBins[pta],pTAssocBins[pta+1]);
     	//c1.SetPID(pdgPi0);
     	Correlators38.push_back(c1);

      Correlator c2(pta,ptt,1);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
     	c2.SetCentrality(0,20);
     	c2.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
     	c2.SetAssociatedRange(pTAssocBins[pta],pTAssocBins[pta+1]);
     	//c1.SetPID(pdgPi0);
     	Correlators38.push_back(c2);

      Correlator c3(pta,ptt,2);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(20,40);
      c3.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c3.SetAssociatedRange(pTAssocBins[pta],pTAssocBins[pta+1]);
      //c1.SetPID(pdgPi0);
      Correlators38.push_back(c3);

      Correlator c4(pta,ptt,3);
      c4.SetCollSystemAndEnergy("AuAu200GeV");
      c4.SetCentrality(60,92);
      c4.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c4.SetAssociatedRange(pTAssocBins[pta],pTAssocBins[pta+1]);
      //c1.SetPID(pdgPi0);
      Correlators38.push_back(c4);
	
	}
}
for(Correlator& corr : Correlators38)
  {
    if(corr.GetSubSubIndex()==-1){
      string name = "58010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
        book(_h[name], 58,01,(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    } 
    else if(corr.GetSubSubIndex()==1){
      string name = "53010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
        book(_h[name], 53,01,(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    } 
    else if(corr.GetSubSubIndex()==2){
      string name = "55010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
        book(_h[name], 55,01,(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0; 
    }
    else if(corr.GetSubSubIndex()==3){
      string name = "57010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
        book(_h[name], 57,01,(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0; 
    }
        
  }

 
//*****************************************************************************
// The following will book the histograms for Figure 31
for(ptt = 0; ptt<numTrigPtBins; ptt++){
  	for(pta = 0; pta < numpTAssocBins12; pta++){
  	for(cb = 0; cb<numCentBins31; cb++){
        Correlator c1(pta,ptt,cb);
        c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(CentBins31[cb],CentBins31[cb+1]);
      c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c1.SetAssociatedRange(pTAssocBins12[pta],pTAssocBins12[pta+1]);
      //c1.SetPID(pdgPi0);
      Correlators31.push_back(c1);
  	}
    Correlator c2(pta,ptt,4); // pp
        c2.SetCollSystemAndEnergy("pp200GeV");
      c2.SetNoCentrality();
      c2.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c2.SetAssociatedRange(pTAssocBins12[pta],pTAssocBins12[pta+1]);
      //c1.SetPID(pdgPi0);
      Correlators31.push_back(c2);
  }
}

for(Correlator& corr : Correlators31)
  {
        string name = "Fig31CorrFunc_" + to_string(49 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
        book(_h[name], name, 49 + corr.GetSubSubIndex(), -M_PI/2, 3*M_PI/2); // Changed from 34, FIXME
        //string name = "Fig12CorrFunc_pp_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
        //book(_h[name], (35 - corr.GetSubSubIndex()),01,(1 + corr.GetSubIndex()+1));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
  }

 for(ptt = 0; ptt<numTrigPtBins; ptt++){

	    string name1 = "49010" + to_string(1+ptt);
        book(_h[name1], 49,01,(1+ptt));
        book(sow[name1],"sow" + name1);
        nTriggers[name1] = 0;

      string name2 = "50010" + to_string(1+ptt);
        book(_h[name2], 50,01,(1+ptt));
        book(sow[name2],"sow" + name2);
        nTriggers[name2] = 0;

      string name3 = "51010" + to_string(1+ptt);
        book(_h[name3], 51,01,(1+ptt));
        book(sow[name3],"sow" + name3);
        nTriggers[name3] = 0; 

      string name4 = "52010" + to_string(1+ptt);
        book(_h[name4], 52,01,(1+ptt));
        book(sow[name4],"sow" + name4);
        nTriggers[name4] = 0; 

    
  }


/*for(ptt = 0; ptt<numTrigPtBins; ptt++){
  for(cb = 0; cb<numCentBins31; cb++){
    Correlator c2(pta,ptt,4); // pp
        c2.SetCollSystemAndEnergy("pp200GeV");
      c2.SetNoCentrality();
      c2.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c2.SetAssociatedRange(pTAssocBins12[pta],pTAssocBins12[pta+1]);
      //c1.SetPID(pdgPi0);
      Correlators31pp.push_back(c2);
  }
}

  for(Correlator& corr : Correlators31pp)
  {
        string name = "Fig31CorrFunc_pp_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
        book(_h[name], name, 53, -M_PI/2, 3*M_PI/2); // Changed from 34, FIXME
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
  }
  	string name = "Fig31CorrFunc_pp_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex()); // FIXME This needs to be 53 
*/


/*for(Correlator& corr : Correlators31)
  {
        if(corr.GetSubSubIndex()==0){
      string name = "49010" + to_string(1+corr.GetSubIndex());
        book(_h[name], 49,01,(1+corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    } 
    else if(corr.GetSubSubIndex()==1){
      string name = "50010" + to_string(1+corr.GetSubIndex());
        book(_h[name], 50,01,(1+corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    } 
    else if(corr.GetSubSubIndex()==2){
      string name = "51010" + to_string(1+corr.GetSubIndex());
        book(_h[name], 51,01,(1+corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0; 
    }
    else if(corr.GetSubSubIndex()==3){
      string name = "52010" + to_string(1+corr.GetSubIndex());
        book(_h[name], 52,01,(1+corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0; 
    }
  } */

 
//*****************************************************************************
// The following will book the histograms for Figures 28-30

for(ptt = 0; ptt<numTrigPtBins; ptt++){
  	for(pta = 0; pta < numpTAssocBins12; pta++){
  	for(cb = 0; cb<numCentBins31; cb++){
        Correlator c1(pta,ptt,cb);
        c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(CentBins31[cb],CentBins31[cb+1]);
      c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c1.SetAssociatedRange(pTAssocBins12[pta],pTAssocBins12[pta+1]);
      //c1.SetPID(pdgPi0);
      Correlators30.push_back(c1);
  	}
    Correlator c2(pta,ptt,4); // pp
        c2.SetCollSystemAndEnergy("pp200GeV");
      c2.SetNoCentrality();
      c2.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c2.SetAssociatedRange(pTAssocBins12[pta],pTAssocBins12[pta+1]);
      //c1.SetPID(pdgPi0);
      Correlators30.push_back(c2);
  }
}

for(Correlator& corr : Correlators30)
  {
        string name = "Fig30CorrFunc_" + to_string(45 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
        book(_h[name], name, 45 + corr.GetSubSubIndex(), -M_PI/2, 3*M_PI/2); // Changed from 34, FIXME
        //string name = "Fig12CorrFunc_pp_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
        //book(_h[name], (35 - corr.GetSubSubIndex()),01,(1 + corr.GetSubIndex()+1));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
  }

 for(ptt = 0; ptt<numTrigPtBins; ptt++){

	    string name1 = "45010" + to_string(1+ptt);
        book(_h[name1], 45,01,(1+ptt));
        book(sow[name1],"sow" + name1);
        nTriggers[name1] = 0;

      string name2 = "46010" + to_string(1+ptt);
        book(_h[name2], 46,01,(1+ptt));
        book(sow[name2],"sow" + name2);
        nTriggers[name2] = 0;

      string name3 = "47010" + to_string(1+ptt);
        book(_h[name3], 47,01,(1+ptt));
        book(sow[name3],"sow" + name3);
        nTriggers[name3] = 0; 

      string name4 = "48010" + to_string(1+ptt);
        book(_h[name4], 48,01,(1+ptt));
        book(sow[name4],"sow" + name4);
        nTriggers[name4] = 0; 

    
  }




/*for(ptt = 0; ptt<numTrigPtBins; ptt++){
  for(cb = 0; cb<numCentBins; cb++){
        Correlator c1(pta,ptt,cb);
        c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(CentBins[cb],CentBins[cb+1]);
      c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c1.SetAssociatedRange(pTAssocBins[pta],pTAssocBins[pta+1]);
      //c1.SetPID(pdgPi0);
      Correlators30.push_back(c1);
  }
}

for(Correlator& corr : Correlators30)
  {
        if(corr.GetSubSubIndex()==0){
      string name = "45010" + to_string(1+corr.GetSubIndex());
        book(_h[name], 45,01,(1+corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    } 
    else if(corr.GetSubSubIndex()==1){
      string name = "46010" + to_string(1+corr.GetSubIndex());
        book(_h[name], 46,01,(1+corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    } 
    else if(corr.GetSubSubIndex()==2){
      string name = "47010" + to_string(1+corr.GetSubIndex());
        book(_h[name], 47,01,(1+corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0; 
    }
    else if(corr.GetSubSubIndex()==3){
      string name = "48010" + to_string(1+corr.GetSubIndex());
        book(_h[name], 48,01,(1+corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0; 
    }
  }*/ 
 
//*****************************************************************************
// The following will book the histograms for Figure 26, where pta2 
for(pta2 = 0; pta2<numpTAssocBins2; pta2++){
	for(cb = 0; cb < numCentBins26; cb++){ 
      Correlator c1(cb,0,pta2);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetNoCentrality();
      c1.SetTriggerRange(3, 4);
      c1.SetAssociatedRange(pTAssocBins2[pta2],7);
      //c1.SetPID(pdgPi0);
      Correlators26.push_back(c1);

      Correlator c2(cb,1,pta2);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetNoCentrality();
      c2.SetTriggerRange(4, 5);
      c2.SetAssociatedRange(pTAssocBins2[pta2],7);
      //c1.SetPID(pdgPi0);
      Correlators26.push_back(c2);
  }
}
for(Correlator& corr : Correlators26)
  {
	  string name = "Fig26CorrFunc_" + to_string(43 + corr.GetSubIndex()) + "_1_" + to_string(1+corr.GetSubSubIndex())  + "_Centrality_" + to_string(1+corr.GetIndex());        
        book(_p[name], name, 1, CentBins26[corr.GetIndex()], CentBins26[corr.GetIndex()+1]);
        //book(_h[name], (43 + corr.GetSubIndex()),01,(1+corr.GetSubSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
  }
  for(pta = 0; pta<numpTAssocBins2; pta++){

      string name2 = "43010" + to_string(pta+1);
      book(_h[name2],43, 1, pta+1);
      // Near side 
      string name3 = "44010" + to_string(pta+1);
      book(_h[name3],44, 1, pta+1);

  }

 
//*****************************************************************************
// The following will book the histograms for Figure 25 
for(ptt = 0; ptt<numTrigPtBins; ptt++){
	for(cb = 0; cb < numCentBins25; cb++){
      Correlator c1(cb,ptt,0);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(CentBins25[cb],CentBins25[cb + 1]);
      c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c1.SetAssociatedRange(1,5);
      //c1.SetPID(pdgPi0);
      Correlators25.push_back(c1);
      /*
      Correlator c2(cb,ptt,1);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetNoCentrality();
      c2.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c2.SetAssociatedRange(1,5);
      //c1.SetPID(pdgPi0);
      Correlators25.push_back(c2);

      Correlator c3(cb,ptt,2);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetNoCentrality();
      c3.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c3.SetAssociatedRange(1,5);
      //c1.SetPID(pdgPi0);
      Correlators25.push_back(c3);
      */
  }
}

/*for(Correlator& corr : Correlators25) 
  {
        string name = (42 - corr.GetSubSubIndex()) + "01" + to_string(1+corr.GetSubIndex());
        book(_h[name], (42 - corr.GetSubSubIndex()),01,(1+corr.GetSubIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
      
  }*/ 

  for(Correlator& corr : Correlators25)
  {
        //string name = "Fig25CorrFunc_" + to_string(42 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
        //string name = "Fig25CorrFunc_42_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
        //book(_p[name], name, 10, 0.,100.);
        //book(_p[name], name, 1, CentBins25[corr.GetIndex()], CentBins25[corr.GetIndex()+1]);
        //book(_h[name], name, 42 - corr.GetSubSubIndex(), -M_PI/2, 3*M_PI/2);
        //string name = "Fig12CorrFunc_pp_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
        //book(_h[name], (35 - corr.GetSubSubIndex()),01,(1 + corr.GetSubIndex()+1));
        //book(_p[name], name, 10, 0.,100.);
        string name1 = "Fig25CorrFunc_40_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
        book(_p[name1], name1, 1, CentBins25[corr.GetIndex()], CentBins25[corr.GetIndex()+1]);
        book(sow[name1],"sow" + name1);
        nTriggers[name1] = 0;
        
        string name2 = "Fig25CorrFunc_41_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
        book(_p[name2], name2, 1, CentBins25[corr.GetIndex()], CentBins25[corr.GetIndex()+1]);
        book(sow[name2],"sow" + name2);
        nTriggers[name2] = 0;

        string name3 = "Fig25CorrFunc_42_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
        book(_p[name3], name3, 1, CentBins25[corr.GetIndex()], CentBins25[corr.GetIndex()+1]);
        book(sow[name3],"sow" + name3);
        nTriggers[name3] = 0;

  }

 for(ptt = 0; ptt<numTrigPtBins; ptt++){
      // Head 
      string name1 = "Figure25_AuAu_40_1_" + to_string(ptt+1);
      book(_h[name1],40, 1, ptt+1);
      // Shoulder 
      string name2 = "Figure25_AuAu_41_1_" + to_string(ptt+1);
      book(_h[name2],41, 1, ptt+1);
      // Near side 
      string name3 = "Figure25_AuAu_42_1_" + to_string(ptt+1);
      book(_h[name3],42, 1, ptt+1);

  }

 
//*****************************************************************************
// The following will book the histograms for Figure 24 
for(i = 0; i<2; i++){
      Correlator c1(1,1,i);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0,20);
      c1.SetTriggerRange(2, 3);
      c1.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators24.push_back(c1);

      Correlator c2(2,2,i);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetCentrality(0,20);
      c2.SetTriggerRange(3, 4);
      c2.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators24.push_back(c2);

      Correlator c3(3,3,i);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0,20);
      c3.SetTriggerRange(3, 4);
      c3.SetAssociatedRange(3,5);
      //c1.SetPID(pdgPi0);
      Correlators24.push_back(c3);

      Correlator c4(4,4,i);
      c4.SetCollSystemAndEnergy("AuAu200GeV");
      c4.SetCentrality(0,20);
      c4.SetTriggerRange(5,10);
      c4.SetAssociatedRange(5,10);
      //c1.SetPID(pdgPi0);
      Correlators24.push_back(c4);

}
for(Correlator& corr : Correlators24)
  {
      	string name = to_string(39 - corr.GetSubSubIndex()) + "010" + to_string(corr.GetIndex());
        book(_h[name], (39 - corr.GetSubSubIndex()),01,(corr.GetIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
      
  }

//*****************************************************************************
// The following will book the histograms for Figure 23 
for(ptt = 0; ptt<1; ptt++){
      Correlator c1(1,ptt,0);
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetNoCentrality();
      c1.SetTriggerRange(2, 3);
      c1.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators23.push_back(c1);

      Correlator c2(2,ptt,0);
      c2.SetCollSystemAndEnergy("pp200GeV");
      c2.SetNoCentrality();
      c2.SetTriggerRange(3, 4);
      c2.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators23.push_back(c2);

      Correlator c3(3,ptt,0);
      c3.SetCollSystemAndEnergy("pp200GeV");
      c3.SetNoCentrality();
      c3.SetTriggerRange(3, 4);
      c3.SetAssociatedRange(3,5);
      //c1.SetPID(pdgPi0);
      Correlators23.push_back(c3);

      Correlator c4(4,ptt,0);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetNoCentrality();
      c4.SetTriggerRange(5,10);
      c4.SetAssociatedRange(5,10);
      //c1.SetPID(pdgPi0);
      Correlators23.push_back(c4);

      Correlator c5(1,5,1);
      c5.SetCollSystemAndEnergy("AuAu200GeV");
      c5.SetCentrality(0,20);
      c5.SetTriggerRange(2, 3);
      c5.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators23.push_back(c5);

      Correlator c6(2,6,1);
      c6.SetCollSystemAndEnergy("AuAu200GeV");
      c6.SetCentrality(0,20);
      c6.SetTriggerRange(3, 4);
      c6.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators23.push_back(c6);

      Correlator c7(3,7,1);
      c7.SetCollSystemAndEnergy("AuAu200GeV");
      c7.SetCentrality(0,20);
      c7.SetTriggerRange(3, 4);
      c7.SetAssociatedRange(3,5);
      //c1.SetPID(pdgPi0);
      Correlators23.push_back(c7);

      Correlator c8(4,8,1);
      c8.SetCollSystemAndEnergy("AuAu200GeV");
      c8.SetCentrality(0,20);
      c8.SetTriggerRange(5,10);
      c8.SetAssociatedRange(5,10);
      //c1.SetPID(pdgPi0);
      Correlators23.push_back(c8);

}
for(Correlator& corr : Correlators23)
  {
        //string name = to_string(37-corr.GetSubSubIndex()) + "01" + to_string(corr.GetIndex()+1);
        //book(_h[name], (37 - corr.GetSubSubIndex()),1,(corr.GetIndex()));
        //book(sow[name],"sow" + name);
        //nTriggers[name] = 0;
        //string name = to_string(37 - corr.GetSubSubIndex()) + "010" + to_string(corr.GetIndex());
        //book(_h[name], (37 - corr.GetSubSubIndex()),01,(corr.GetIndex()));
        //book(sow[name],"sow" + name);
        //nTriggers[name] = 0;
        if(corr.GetSubIndex() > 4)
        {
          string name = "36010" + to_string(corr.GetIndex());
          book(_h[name], 36,01,(corr.GetIndex()));
          book(sow[name],"sow" + name);
          nTriggers[name] = 0;
        }
        //else
          /* FIXME pp, bin definitions acting up
        {
          string name = "35010" + to_string(corr.GetIndex());
          book(_h[name], 37,01,(corr.GetIndex()));
          book(sow[name],"sow" + name);
          nTriggers[name] = 0;
        }
        */

  }

 
//*****************************************************************************
// The following will book the histograms for Figure 18 

for(ptt = 0; ptt<numTrigPtBins; ptt++){
   for(pta = 0; pta < numpTAssocBins12; pta++){
      Correlator c1(pta,ptt,1);
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetNoCentrality();
      c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c1.SetAssociatedRange(pTAssocBins12[pta],pTAssocBins12[pta+1]);
      //c1.SetPID(pdgPi0);
      Correlators18.push_back(c1);

      Correlator c2(pta,ptt,0);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetCentrality(0,20);
      c2.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c2.SetAssociatedRange(pTAssocBins12[pta],pTAssocBins12[pta+1]);  
      //c1.SetPID(pdgPi0);
      Correlators18.push_back(c2);
  }
}

/*for(Correlator& corr : Correlators18)
  {
  	if(corr.GetSubSubIndex() == 0){
        //string name = "Fig18CorrFunc_AuAu_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
        string name = "Fig18CorrFunc_AuAu_35_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
        book(_h[name], name, 35, -M_PI/2, 3*M_PI/2);
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    }
    else if(corr.GetSubSubIndex() == 1){
    	//string name = "Fig18CorrFunc_pp_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
        string name = "Fig18CorrFunc_pp_34_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
        book(_h[name], name, 34, -M_PI/2, 3*M_PI/2);
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    }
  }*/

for(Correlator& corr : Correlators18)
  {
        string name = "Fig18CorrFunc_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
        book(_h[name], name, 35 - corr.GetSubSubIndex(), -M_PI/2, 3*M_PI/2); // Changed from 34, FIXME
        //string name = "Fig12CorrFunc_pp_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
        //book(_h[name], (35 - corr.GetSubSubIndex()),01,(1 + corr.GetSubIndex()+1));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
  }

 for(ptt = 0; ptt<numTrigPtBins; ptt++){
      //AuAu
      string name = "Figure18_AuAu_35_1_" + to_string(ptt+1);
      book(_h[name],35, 1, ptt+1);
      //pp
      string nameNS = "Figure18_pp_34_1_" + to_string(ptt+1);
      book(_h[nameNS],34, 1, ptt+1);
  }


 
//*****************************************************************************
// The following will book the histograms for Figure 17 
for(ptt = 0; ptt<1; ptt++){
      Correlator c1(1,ptt,1);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetNoCentrality();
      c1.SetTriggerRange(2, 3);
      c1.SetAssociatedRange(0.4,1);
      //c1.SetPID(pdgPi0);
      Correlators17.push_back(c1);

      Correlator c2(2,ptt,1);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetNoCentrality();
      c2.SetTriggerRange(2,3);
      c2.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators17.push_back(c2);

      Correlator c3(3,ptt,1);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetNoCentrality();
      c3.SetTriggerRange(3, 4);
      c3.SetAssociatedRange(3,4);
      //c1.SetPID(pdgPi0);
      Correlators17.push_back(c3);

      Correlator c4(4,ptt,1);
      c4.SetCollSystemAndEnergy("AuAu200GeV");
      c4.SetNoCentrality();
      c4.SetTriggerRange(4,5);
      c4.SetAssociatedRange(4,5);
      //c1.SetPID(pdgPi0);
      Correlators17.push_back(c4);

      Correlator c5(5,ptt,1);
      c5.SetCollSystemAndEnergy("AuAu200GeV");
      c5.SetNoCentrality();
      c5.SetTriggerRange(5,10);
      c5.SetAssociatedRange(5,10);
      //c1.SetPID(pdgPi0);
      Correlators17.push_back(c5);

}
for(Correlator& corr : Correlators17)
  {
        string name = "33010" + to_string(corr.GetIndex());
        book(_h[name], 33,01,corr.GetIndex());
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
  }

//*****************************************************************************
// The following will book the histograms for Figure 16 

for(ptt = 0; ptt<numTrigPtBins-1; ptt++){
      Correlator c1(-1,ptt,1);
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetNoCentrality();
      c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c1.SetNoAssoc();
      //c1.SetPID(pdgPi0);
      Correlators16.push_back(c1);

      Correlator c2(1,ptt,0);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetCentrality(0,20);
      c2.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c2.SetNoAssoc(); 
      //c1.SetPID(pdgPi0);
      Correlators16.push_back(c2);
}
for(Correlator& corr : Correlators16)
  {
        string name = to_string(32 - corr.GetSubSubIndex()) + "010" + to_string(corr.GetSubIndex()+1);
        book(_h[name], (35 - corr.GetSubSubIndex()),01,(1 + corr.GetSubIndex()+1));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
  }
 
//*****************************************************************************
// The following will book the histograms for Figure 12 
for(ptt = 0; ptt<numTrigPtBins; ptt++){
  for(pta = 0; pta < numpTAssocBins12; pta++){
      Correlator c1(pta,ptt,4);
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetNoCentrality();
      c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c1.SetAssociatedRange(pTAssocBins12[pta],pTAssocBins12[pta+1]); 
      //c1.SetPID(pdgPi0);
      Correlators12.push_back(c1);

    for(cb = 0; cb<numCentBins12; cb++){
        Correlator c1(pta,ptt,cb);
        c1.SetCollSystemAndEnergy("AuAu200GeV");
        c1.SetCentrality(CentBins12[cb],CentBins12[cb+1]);
        c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
        c1.SetAssociatedRange(pTAssocBins12[pta],pTAssocBins12[pta+1]);
        //c1.SetPID(pdgPi0);
        Correlators12.push_back(c1);
   }
  }
}
for(Correlator& corr : Correlators12)
  {
  	if(corr.GetSubSubIndex() < 4){
        string name = "Fig12CorrFunc_AuAu_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
        book(_h[name], name, 36, -M_PI/2, 3*M_PI/2);
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    }
    else if(corr.GetSubSubIndex()==4){
        string name = "Fig12CorrFunc_pp_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
        book(_h[name], name, 36, -M_PI/2, 3*M_PI/2);
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    }
  }
  
  for(ptt = 0; ptt<numTrigPtBins; ptt++){
      //AuAu
      for(cb = 0; cb<numCentBins12; cb++){
          string name = "Figure12_AuAu_" + to_string(21+cb) + "_1_" + to_string(ptt+1);
          book(_h[name],21+cb, 1, ptt+1);
          string nameSR = "Figure12_AuAu_" + to_string(26+cb) + "_1_" + to_string(ptt+1);
          book(_h[nameSR],26+cb, 1, ptt+1);
      }
      //pp
      string name = "Figure12_pp_25_1_" + to_string(ptt+1);
      book(_h[name],25, 1, ptt+1);
      string nameSR = "Figure12_pp_" + to_string(26+cb) + "_1_" + to_string(ptt+1);
      book(_h[nameSR],26+cb, 1, ptt+1);
  }

//*****************************************************************************
// The following will book the histograms for Figure 11
for(i=0;i<2;i++){
  for(ptt = 0; ptt<numTrigPtBins-2; ptt++){
      Correlator c1(1,ptt,i);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0,20);
      c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c1.SetNoAssoc();
      //c1.SetPID(pdgPi0);
      Correlators11.push_back(c1);
  }
}
for(Correlator& corr : Correlators11)
  {

	if(corr.GetSubSubIndex()==0){
        string name = (18 + corr.GetSubIndex()) + "0101" + to_string(corr.GetIndex());
        book(_h[name], (18 + corr.GetSubIndex()),01,(corr.GetIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    }
    else if(corr.GetSubSubIndex()==1){
        string name = (15 + corr.GetSubIndex()) + "01" + to_string(corr.GetIndex());
        book(_h[name], (15 + corr.GetSubIndex()),01,(corr.GetIndex()));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
    }
  }
 //
//*****************************************************************************
//The following will book the histograms for Figure 10, where i represents FIT number 
for(i=0;i<3;i++){
      Correlator c1(i,1,1);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetNoCentrality();
      c1.SetTriggerRange(2,3);
      c1.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators10.push_back(c1);
}
for(Correlator& corr : Correlators10)
{
        string name = to_string(corr.GetSubSubIndex()) + "010" + to_string(corr.GetIndex()+1);
        book(_h[name], corr.GetSubSubIndex(),01,corr.GetIndex()+1);
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
  }
 
//*****************************************************************************
// The following will book the histograms for Figure 9 
for(i=0;i<1;i++){
      Correlator c1(1,1,-1);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0,5);
      c1.SetTriggerRange(2,3);
      c1.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators9.push_back(c1);
}
for(Correlator& corr : Correlators9)
  {
        string name = "130101";
        book(_h[name], 13,01,01);
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
  }
 
//*****************************************************************************
// The following will book the histograms for Figure 8 
for(cb = 0; cb < numCentBins26; cb++){
      Correlator c1(cb,1,-1);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(CentBins26[cb],CentBins26[cb + 1]);
      c1.SetTriggerRange(2,3);
      c1.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators8.push_back(c1);
}

for(Correlator& corr : Correlators8)
  {
        string namehead = "Fig8CorrFunc_H_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
        book(_h[namehead], namehead, 12, -M_PI/2, 3*M_PI/2);
        book(sow[namehead],"sow" + namehead);
        nTriggers[namehead] = 0;

        string nameshoulder = "Fig8CorrFunc_S_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
        book(_h[nameshoulder], nameshoulder, 12, -M_PI/2, 3*M_PI/2);
        book(sow[nameshoulder],"sow" + nameshoulder);
        nTriggers[nameshoulder] = 0;
  }

  string name = "120101";
  book(_h[name], 12, 1 , 1);
  book(sow[name],"sow" + name);
  nTriggers[name] = 0;
 
//*****************************************************************************
// The following will book the histograms for Figure 7 

for(ptt = 0; ptt<numTrigPtBins; ptt++){
      Correlator c1(1,ptt,1);
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetNoCentrality();
      c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c1.SetNoAssoc();
      //c1.SetPID(pdgPi0);
      Correlators7.push_back(c1);

      Correlator c2(1,ptt,0);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetCentrality(0,20);
      c2.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c2.SetNoAssoc(); 
      //c1.SetPID(pdgPi0);
      Correlators7.push_back(c2);
}
for(Correlator& corr : Correlators7)
  {
        string name = to_string(9-corr.GetSubSubIndex()) + "010" + to_string(corr.GetSubIndex() + 1);
        book(_h[name], (9-corr.GetSubSubIndex()),01,(corr.GetSubIndex() + 1));
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
  }
 
//*****************************************************************************
// The following will book the histograms for Figure 6 
for(i=0;i<1;i++){
   for(pta = 0; pta<4; pta++){

      Correlator c1(pta,1,0);
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetNoCentrality();
      c1.SetTriggerRange(3,4);
      c1.SetAssociatedRange(pTAssocBins[pta],pTAssocBins[pta+1]);
      //c1.SetPID(pdgPi0);
      Correlators6.push_back(c1);

    }

      Correlator c5(5,2,0);
      c5.SetCollSystemAndEnergy("AuAu200GeV");
      c5.SetNoCentrality();
      c5.SetTriggerRange(5,10);
      c5.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators6.push_back(c5);

      Correlator c6(6,2,0);
      c6.SetCollSystemAndEnergy("AuAu200GeV");
      c6.SetNoCentrality();
      c6.SetTriggerRange(4,5);
      c6.SetAssociatedRange(4,5);
      //c1.SetPID(pdgPi0);
      Correlators6.push_back(c6);

      Correlator c7(7,2,0);
      c7.SetCollSystemAndEnergy("AuAu200GeV");
      c7.SetNoCentrality();
      c7.SetTriggerRange(5,10);
      c7.SetAssociatedRange(3,5);
      //c1.SetPID(pdgPi0);
      Correlators6.push_back(c7);

      Correlator c8(8,2,0);
      c8.SetCollSystemAndEnergy("AuAu200GeV");
      c8.SetNoCentrality();
      c8.SetTriggerRange(5,10);
      c8.SetAssociatedRange(5,10);
      //c1.SetPID(pdgPi0);
      Correlators6.push_back(c8);

  for(pta = 0; pta<4; pta++){
   for(ptt = 0; ptt<4; ptt++){

      Correlator c9(pta,1,1);
      c9.SetCollSystemAndEnergy("pp200GeV");
      c9.SetCentrality(0,20);
      c9.SetTriggerRange(3,4);
      c9.SetAssociatedRange(pTAssocBins[pta],pTAssocBins[pta+1]);
      //c1.SetPID(pdgPi0);
      Correlators6.push_back(c9);

    }
   } 

      Correlator c11(5,2,1);
      c11.SetCollSystemAndEnergy("AuAu200GeV");
      c11.SetCentrality(0,20);
      c11.SetTriggerRange(5,10);
      c11.SetAssociatedRange(2,3);
      //c1.SetPID(pdgPi0);
      Correlators6.push_back(c11);

      Correlator c12(6,2,1);
      c12.SetCollSystemAndEnergy("AuAu200GeV");
      c12.SetCentrality(0,20);
      c12.SetTriggerRange(4,5);
      c12.SetAssociatedRange(4,5);
      //c1.SetPID(pdgPi0);
      Correlators6.push_back(c12);

      Correlator c13(7,2,1);
      c13.SetCollSystemAndEnergy("AuAu200GeV");
      c13.SetCentrality(0,20);
      c13.SetTriggerRange(5,10);
      c13.SetAssociatedRange(3,5);
      //c1.SetPID(pdgPi0);
      Correlators6.push_back(c13);

      Correlator c14(8,2,1);
      c14.SetCollSystemAndEnergy("AuAu200GeV");
      c14.SetCentrality(0,20);
      c14.SetTriggerRange(5,10);
      c14.SetAssociatedRange(5,10);
      //c1.SetPID(pdgPi0);
      Correlators6.push_back(c14);
}
for(Correlator& corr : Correlators6)
  {
  	if(corr.GetSubSubIndex()==0){
  		if(corr.GetSubIndex()==1){
  			string name = "07010" + to_string(corr.GetIndex()+1);
  			book(_h[name], 7,01,corr.GetIndex()+1);
        	book(sow[name],"sow" + name);
        	nTriggers[name] = 0;
  		}
  		else if(corr.GetSubIndex()==2){
  			string name = "06010" + to_string(corr.GetIndex()-4);
  			book(_h[name], 6,01,corr.GetIndex()-4);
        	book(sow[name],"sow" + name);
        	nTriggers[name] = 0;
  		}
  	}
  	else if(corr.GetSubSubIndex()==1){
  		if(corr.GetSubIndex()==1){
  			string name = "08010" + to_string(corr.GetIndex()+1);
  			book(_h[name], 8,01,corr.GetIndex()+1);
        	book(sow[name],"sow" + name);
        	nTriggers[name] = 0;
  		}
  		else if(corr.GetSubIndex()==2){
  			string name = "09010" + to_string(corr.GetIndex()-4);
  			book(_h[name], 9,01,corr.GetIndex()-4);
        	book(sow[name],"sow" + name);
        	nTriggers[name] = 0;
  		}
  	}

  }

//*****************************************************************************
// The following will book the histograms for Figure 4, where i is the NS,AS,AH 
for(i=0;i<3;i++){
  for(ptt = 0; ptt<numTrigPtBins-1; ptt++){
      Correlator c1(i,ptt,1);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0,20);
      c1.SetTriggerRange(pTTrigBins[ptt], pTTrigBins[ptt+1]);
      c1.SetNoAssoc();
      //c1.SetPID(pdgPi0);
      Correlators4.push_back(c1);
  }
}
for(Correlator& corr : Correlators4)
  {
        string name = to_string(corr.GetSubIndex()+1) + "010" + to_string(corr.GetIndex()+1);
        book(_h[name], corr.GetSubIndex()+1,01,corr.GetIndex()+1);
        book(sow[name],"sow" + name);
        nTriggers[name] = 0;
  }


//for(Correlator& corr : CorrelatorsB)
  //{
    //    book(_h["0" + to_string(corr.GetIndex()) + "11"], corr.GetIndex(), 1, 1);
    //    book(sow[corr.GetIndex()],"sow" + to_string(corr.GetIndex()));
    //    nTriggers[corr.GetIndex()] = 0;
   // }


    }
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const PrimaryParticles& ppTrigPi0 = apply<PrimaryParticles>(event, "PP");
      const PromptFinalState& pfsTrigPhotons = apply<PromptFinalState>(event, "PFS");
      //==================================================
      // Select the histograms accordingly to the collision system, beam energy and centrality
      // WARNING: Still not implemented for d-Au
      //==================================================
      string SysAndEnergy = "";
      //string cmsEnergy = "200GeV";

      if(collSys == pp0) SysAndEnergy = "pp200GeV";
      else if(collSys == AuAu) SysAndEnergy = "AuAu200GeV";
    
      //SysAndEnergy = collSys + cmsEnergy;
   
      double triggerptMin = 999.;
      double triggerptMax = -999.;
      double associatedptMin = 999.;
      double associatedptMax = -999.;
   
      bool isVeto = true;

      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();
      
       /*
      for(Correlator& corr : Correlators) // Think about adding for Correlators38 etc. FIXME 
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
      */
      //*****************************************************************************
      // The following will fill the sow for Figure 36-8 
      for(Correlator& corr : Correlators38)
        {
                  
                  
                  if(!corr.CheckCentrality(c)) continue;
                  
                  if(corr.GetSubSubIndex()==-1){
                  string name = "58010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  sow[name]->fill();
                  } 
                  else if(corr.GetSubSubIndex()==1){
                  string name = "53010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  sow[name]->fill();
                  } 
                  else if(corr.GetSubSubIndex()==2){
                  string name = "55010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  sow[name]->fill();
                  }
                  //FIXME This one causes an error
                  else if(corr.GetSubSubIndex()==3){
                  string name = "57010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  sow[name]->fill();
                  }
                  

                  //string name = "58010" + to_string((corr.GetIndex()+1) + corr.GetSubIndex());
                  
                  //sow[name]->fill();
                  
        }
      //*****************************************************************************
      // The following will fill the sow for Figure 6
        //FIXME ... the figure is for pT 
        
        for(Correlator& corr : Correlators6)
        {

                  if(!corr.CheckCentrality(c)) continue;
                  
                  
                if(corr.GetSubSubIndex()==0){
                  if(corr.GetSubIndex()==1){
                   string name = "07010" + to_string(corr.GetIndex()+1);
                   sow[name]->fill();
                  }
                 else if(corr.GetSubIndex()==2){
                  string name = "06010" + to_string(corr.GetIndex()-4);
                  sow[name]->fill();
                 }
                }
                else if(corr.GetSubSubIndex()==1){
                 if(corr.GetSubIndex()==1){
                   string name = "08010" + to_string(corr.GetIndex()+1);
                   sow[name]->fill();
                  }
                 else if(corr.GetSubIndex()==2){
                  string name = "09010" + to_string(corr.GetIndex()-4);
                  sow[name]->fill();
                }
              }
        }
          

      //*****************************************************************************
      // The following will fill the sow for Figure 12 
      for(Correlator& corr : Correlators12)
      {
          if(!corr.CheckCentrality(c)) continue;

          if(corr.GetSubSubIndex() < 4){
              string name = "Fig12CorrFunc_AuAu_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
              sow[name]->fill();
          }
          else if(corr.GetSubSubIndex() == 4){
              string name = "Fig12CorrFunc_pp_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
              sow[name]->fill();
          }

      }

      //*****************************************************************************
      // The following will fill the sow for Figure 18 
      
      for(Correlator& corr : Correlators18)
      {

          if(!corr.CheckCentrality(c)) continue;
          string name = "Fig18CorrFunc_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
			//string name = "Fig18CorrFunc_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
			sow[name]->fill();
          /*if(corr.GetSubSubIndex() == 0){
              //string name = "Fig18CorrFunc_AuAu_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
              string name = "Fig18CorrFunc_AuAu_35_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
              sow[name]->fill();
          }
          else if(corr.GetSubSubIndex() == 1){
              //string name = "Fig18CorrFunc_AuAu_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
              string name = "Fig18CorrFunc_pp_34_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
              sow[name]->fill();
          } */

      } //FIXME


      //*****************************************************************************
      // The following will fill the sow for Figure 25 
      
      for(Correlator& corr : Correlators25)
      {
          if(!corr.CheckCentrality(c)) continue;


              string name1 = "Fig25CorrFunc_40_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
              sow[name1]->fill();
				
              string name2 = "Fig25CorrFunc_41_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
              sow[name2]->fill();

              string name3 = "Fig25CorrFunc_42_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
              sow[name3]->fill();
      }


      //*****************************************************************************
      // The following will fill the sow for Figure 8 

      for(Correlator& corr : Correlators8)
      {
          if(!corr.CheckCentrality(c)) continue;


        string namehead = "Fig8CorrFunc_H_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
        sow[namehead]->fill();

        string nameshoulder = "Fig8CorrFunc_S_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
        sow[nameshoulder]->fill();

      }



      //*****************************************************************************
      // The following will fill the sow for Figure 26 
      
      for(Correlator& corr : Correlators26)
      {
          if(!corr.CheckCentrality(c)) continue;


			string name = "Fig26CorrFunc_" + to_string(43 + corr.GetSubIndex()) + "_1_" + to_string(1+corr.GetSubSubIndex())  + "_Centrality_" + to_string(1+corr.GetIndex());
            sow[name]->fill();
      }
        
      

      for(Correlator& corr : Correlators31)
      {
          if(!corr.CheckCentrality(c)) continue;

         string name = "Fig31CorrFunc_" + to_string(49 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
			sow[name]->fill();
    

    }


    for(Correlator& corr : Correlators30)
      {
          if(!corr.CheckCentrality(c)) continue;

         string name = "Fig30CorrFunc_" + to_string(45 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
			sow[name]->fill();
    

    }



     for(const Particle& pTrig : cfs.particles())
     {
         //Check if is secondary    FIXME, omitted for run time 
         //if(isSecondary(pAssoc)) continue;
          
         //Trigger counting Figure 38
         for(Correlator& corr : Correlators38)
         {
             if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
             if(!corr.CheckCentrality(c)) continue;
              
             if(corr.GetSubSubIndex()==-1){
                  string name = "58010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  nTriggers[name]++;
             } 
             else if(corr.GetSubSubIndex()==1){
                  string name = "53010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  nTriggers[name]++;
             } 
             else if(corr.GetSubSubIndex()==2){
                  string name = "55010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  nTriggers[name]++;
             }
             else if(corr.GetSubSubIndex()==3){
                  string name = "57010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  nTriggers[name]++;
             }
         }

         //Trigger counting Figure 12
         for(Correlator& corr : Correlators12)
         {
             if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
             if(!corr.CheckCentrality(c)) continue;

             if(corr.GetSubSubIndex() < 4){
               string name = "Fig12CorrFunc_AuAu_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
               nTriggers[name]++;
                  }
             else if(corr.GetSubSubIndex() == 4){
               string name = "Fig12CorrFunc_pp_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
                nTriggers[name]++;
                  }

          }

          //Trigger counting Figure 18
          
         for(Correlator& corr : Correlators18)
         {
             if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
             if(!corr.CheckCentrality(c)) continue;

             string name = "Fig18CorrFunc_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
             nTriggers[name]++;
             /*if(corr.GetSubSubIndex() == 0){
               string name = "Fig18CorrFunc_AuAu_35_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
               nTriggers[name]++;
                  }
             else if(corr.GetSubSubIndex() == 1){
               string name = "Fig18CorrFunc_pp_34_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
                nTriggers[name]++;
                  }*/

          } //FIXME 

          //Trigger counting Figure 25
          
         for(Correlator& corr : Correlators25)
         {
             if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
             if(!corr.CheckCentrality(c)) continue;

             string name1 = "Fig25CorrFunc_40_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
             nTriggers[name1]++;

             string name2 = "Fig25CorrFunc_41_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
             nTriggers[name2]++;

             string name3 = "Fig25CorrFunc_42_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
             nTriggers[name3]++;

          } //FIXME 



          //Trigger counting Figure 8 
          
         for(Correlator& corr : Correlators8)
         {
             if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
             if(!corr.CheckCentrality(c)) continue;

        string namehead = "Fig8CorrFunc_H_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
        nTriggers[namehead]++;

        string nameshoulder = "Fig8CorrFunc_S_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
        nTriggers[nameshoulder]++;

          } 




          //Trigger counting Figure 26
          
         for(Correlator& corr : Correlators26)
         {
             if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
             if(!corr.CheckCentrality(c)) continue;

             string name = "Fig26CorrFunc_" + to_string(43 + corr.GetSubIndex()) + "_1_" + to_string(1+corr.GetSubSubIndex())  + "_Centrality_" + to_string(1+corr.GetIndex());
             nTriggers[name]++;

          } //FIXME 



          for(Correlator& corr : Correlators31)
         {
             if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
             if(!corr.CheckCentrality(c)) continue;

             string name = "Fig31CorrFunc_" + to_string(49 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
				nTriggers[name]++;
    			


          } 


          for(Correlator& corr : Correlators30)
         {
             if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
             if(!corr.CheckCentrality(c)) continue;

             string name = "Fig30CorrFunc_" + to_string(45 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
				nTriggers[name]++;
    			


          } 

          		
          
          
        //Trigger counting Figure 6
         for(Correlator& corr : Correlators6)
         {
             if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
             if(!corr.CheckCentrality(c)) continue;

            if(corr.GetSubSubIndex()==0){
                  if(corr.GetSubIndex()==1){
                   string name = "07010" + to_string(corr.GetIndex()+1);
                   nTriggers[name]++;
                  }
                 else if(corr.GetSubIndex()==2){
                  string name = "06010" + to_string(corr.GetIndex()-4);
                  nTriggers[name]++;
                 }
                }
                else if(corr.GetSubSubIndex()==1){
                 if(corr.GetSubIndex()==1){
                   string name = "08010" + to_string(corr.GetIndex()+1);
                  nTriggers[name]++;
                  }
                 else if(corr.GetSubIndex()==2){
                  string name = "09010" + to_string(corr.GetIndex()-4);
                  nTriggers[name]++;
                } 
               }

        }





          for(const Particle& pTAssoc : cfs.particles())
          {
              //Check if Trigger and Associated are the same particle
              if(isSameParticle(pTrig,pTAssoc)) continue; //I changed this FIXME
                
              //Check if is secondary
              //if(isSecondary(pAssoc)) continue;
              
            //*****************************************************************************
            // The following will fill the histograms for Figure 38 
              for(Correlator& corr : Correlators38)
              {
                  if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                  
                  if(!corr.CheckAssociatedRange(pTAssoc.pt()/GeV)) continue;
                  
                  if(!corr.CheckCentrality(c)) continue;
                  
                  double DeltaPhi = GetDeltaPhi(pTrig, pTAssoc);

                  if(corr.GetSubSubIndex()==-1){
                   string name = "58010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  _h[name]->fill(DeltaPhi);
                  } 
                  else if(corr.GetSubSubIndex()==1){
                  string name = "53010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  _h[name]->fill(DeltaPhi);
                  } 
                  else if(corr.GetSubSubIndex()==2){
                  string name = "55010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  _h[name]->fill(DeltaPhi);
                  }
                  else if(corr.GetSubSubIndex()==3){
                  string name = "57010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
                  _h[name]->fill(DeltaPhi);
                  }
                  
                  
                  //string name = "58010" + to_string((corr.GetIndex()+1) + corr.GetSubIndex());
                  //_h[name]->fill(DeltaPhi);
                
                
                    //nTriggers[name]++;
  
                  
              }
              //*****************************************************************************
              // The following will fill the histograms for Figure 6
              //FIXME ... the figure is for pT 
              for(Correlator& corr : Correlators6)
              {
                  if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                  
                  if(!corr.CheckAssociatedRange(pTAssoc.pt()/GeV)) continue;
                  
                  if(!corr.CheckCentrality(c)) continue;
                  
                  double DeltaPhi = GetDeltaPhi(pTrig, pTAssoc);

                if(corr.GetSubSubIndex()==0){
                  if(corr.GetSubIndex()==1){
                   string name = "07010" + to_string(corr.GetIndex()+1);
                   _h[name]->fill(DeltaPhi);
                  }
                 else if(corr.GetSubIndex()==2){
                  string name = "06010" + to_string(corr.GetIndex()-4);
                  _h[name]->fill(DeltaPhi);
                 }
                }
                else if(corr.GetSubSubIndex()==1){
                 if(corr.GetSubIndex()==1){
                   string name = "08010" + to_string(corr.GetIndex()+1);
                   _h[name]->fill(DeltaPhi);
                  }
                 else if(corr.GetSubIndex()==2){
                  string name = "09010" + to_string(corr.GetIndex()-4);
                  _h[name]->fill(DeltaPhi);
                } 
               }
              }
          //*****************************************************************************
              // The following will fill the histograms for Figure 12
          for(Correlator& corr : Correlators12)
          {
                  if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                  
                  if(!corr.CheckAssociatedRange(pTAssoc.pt()/GeV)) continue;
                  
                  if(!corr.CheckCentrality(c)) continue;


                  
                  double DeltaPhi = GetDeltaPhi(pTrig, pTAssoc);

                  if(corr.GetSubSubIndex() < 4){
                     string name = "Fig12CorrFunc_AuAu_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
                    _h[name]->fill(DeltaPhi);
                  }
                  else if(corr.GetSubSubIndex() == 4){
                   string name = "Fig12CorrFunc_pp_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
                    _h[name]->fill(DeltaPhi);
                  }

          }

          //*****************************************************************************
              // The following will fill the histograms for Figure 18
          
          for(Correlator& corr : Correlators18)
          {
                  if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                  
                  if(!corr.CheckAssociatedRange(pTAssoc.pt()/GeV)) continue;
                  
                  if(!corr.CheckCentrality(c)) continue;
                  
                  double DeltaPhi = GetDeltaPhi(pTrig, pTAssoc);
					
				string name = "Fig18CorrFunc_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
				_h[name]->fill(DeltaPhi);
                  /*if(corr.GetSubSubIndex() == 0){
                     string name = "Fig18CorrFunc_AuAu_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
                     //string name = "Fig18CorrFunc_AuAu_35_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
                    _h[name]->fill(DeltaPhi);
                  }
                  else if(corr.GetSubSubIndex() == 1){
                   string name = "Fig18CorrFunc_pp_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
                   //string name = "Fig18CorrFunc_pp_34_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
                    _h[name]->fill(DeltaPhi);
                  }*/


          }  //FIXME 



          //*****************************************************************************
              // The following will fill the histograms for Figure 25
          
          for(Correlator& corr : Correlators25)
          {
                  if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                  
                  if(!corr.CheckAssociatedRange(pTAssoc.pt()/GeV)) continue;
                  
                  if(!corr.CheckCentrality(c)) continue;
                  
                  double DeltaPhi = GetDeltaPhi(pTrig, pTAssoc);

                  if(DeltaPhi < (M_PI + M_PI/6) && DeltaPhi > (M_PI - M_PI/6)){
                  	string name1 = "Fig25CorrFunc_40_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
					_p[name1]->fill(corr.GetIndex(), pTAssoc.pT()/GeV);
					//_h[name]->fill(DeltaPhi);
                  }

                  if((DeltaPhi < M_PI-M_PI/6. && DeltaPhi > M_PI/2.) || (DeltaPhi < 3.*M_PI/2. && DeltaPhi > M_PI+M_PI/6.)){
                  	string name2 = "Fig25CorrFunc_41_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
					_p[name2]->fill(corr.GetIndex(), pTAssoc.pT()/GeV);
					//_h[name]->fill(DeltaPhi);
                  }

                  if(DeltaPhi < M_PI/3 && DeltaPhi > -M_PI/3){
                  	string name3 = "Fig25CorrFunc_42_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
					_p[name3]->fill(corr.GetIndex(), pTAssoc.pT()/GeV);
					//_h[name]->fill(DeltaPhi);
                  }


          }  //FIXME 



          //*****************************************************************************
              // The following will fill the histograms for Figure 26 
          
          for(Correlator& corr : Correlators26)
          {
                  if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                  
                  if(!corr.CheckAssociatedRange(pTAssoc.pt()/GeV)) continue;
                  
                  if(!corr.CheckCentrality(c)) continue;
                  
                  double DeltaPhi = GetDeltaPhi(pTrig, pTAssoc);

                  if(DeltaPhi < (M_PI + M_PI/6) && DeltaPhi > (M_PI - M_PI/6)){
                  	string name = "Fig26CorrFunc_" + to_string(43 + corr.GetSubIndex()) + "_1_" + to_string(1+corr.GetSubSubIndex())  + "_Centrality_" + to_string(1+corr.GetIndex());
					_p[name]->fill(corr.GetIndex(), pTAssoc.pT()/GeV);
					//_h[name]->fill(DeltaPhi);
                  }


          }  



          //*****************************************************************************
              // The following will fill the histograms for Figure 8 
          
          for(Correlator& corr : Correlators8)
          {
                  if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                  
                  if(!corr.CheckAssociatedRange(pTAssoc.pt()/GeV)) continue;
                  
                  if(!corr.CheckCentrality(c)) continue;
                  
                  double DeltaPhi = GetDeltaPhi(pTrig, pTAssoc);

                    string namehead = "Fig8CorrFunc_H_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
                    _h[namehead]->fill(corr.GetIndex(), pTAssoc.pT()/GeV);
                  
                    string nameshoulder = "Fig8CorrFunc_S_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
                    _h[nameshoulder]->fill(corr.GetIndex(), pTAssoc.pT()/GeV);
                  /*if(DeltaPhi < (M_PI + M_PI/6) && DeltaPhi > (M_PI - M_PI/6)){
                    string namehead = "Fig8CorrFunc_H_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
                    _h[namehead]->fill(corr.GetIndex(), pTAssoc.pT()/GeV);
                  }

                  if((DeltaPhi < M_PI-M_PI/6. && DeltaPhi > M_PI/2.) || (DeltaPhi < 3.*M_PI/2. && DeltaPhi > M_PI+M_PI/6.)){
                    string nameshoulder = "Fig8CorrFunc_S_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
                    _h[nameshoulder]->fill(corr.GetIndex(), pTAssoc.pT()/GeV);
                  }*/



          } 


          
			

              //*****************************************************************************
              // The following will fill the histograms for Figure 23

              for(Correlator& corr : Correlators23)
              {
                  if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                  
                  if(!corr.CheckAssociatedRange(pTAssoc.pt()/GeV)) continue;
                  
                  if(!corr.CheckCentrality(c)) continue;
                  
                  double DeltaEta = pTrig.eta() - pTAssoc.eta();
                      // Name is only for AuAu, see above FIXME 
                   string name = "36010" + to_string(corr.GetIndex());
                   _h[name]->fill(DeltaEta);

              }
              
              for(Correlator& corr : Correlators24)
              {
                  if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                  
                  if(!corr.CheckAssociatedRange(pTAssoc.pt()/GeV)) continue;
                  
                  if(!corr.CheckCentrality(c)) continue;
                  
                  double DeltaEta = pTrig.eta() - pTAssoc.eta();
                  double DeltaPhi = GetDeltaPhi(pTrig, pTAssoc);
                      // Name is only for AuAu, see above FIXME 

                  	string name = to_string(39 - corr.GetSubSubIndex()) + "010" + to_string(corr.GetIndex());
                   	if(corr.GetSubSubIndex()==1) {
                   		_h[name]->fill(DeltaEta);
                  	}
                  else if(corr.GetSubSubIndex()==0) {
                   	_h[name]->fill(DeltaPhi);
                  }
                  
                  /*if(corr.GetSubSubIndex()==1){
                  	string name = "38010" + to_string(corr.GetIndex());
                   	_h[name]->fill(DeltaEta);
                  }
                  
                  else if(corr.GetSubSubIndex()==0){
                  	string name = "39010" + to_string(corr.GetIndex());
                   	_h[name]->fill(DeltaEta);
                  }*/


              }

               for(Correlator& corr : Correlators25)
     		 {
			//string name = "Fig25CorrFunc_42_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
			//string name1 = "Figure25_AuAu_40_1_" + to_string(ptt+1); // + 41, 42 

            /*string name = "Fig25CorrFunc_42_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
            _p[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            _p[name] = SubtractBackgroundZYAM(_p[name]);*/
            //cerr << pTAssoc << endl;
            /*
            string name1 = "Figure25_AuAu_40_1_" + to_string(corr.GetSubIndex()+1);
            double fraction = 0.;
            double yield1 = getYieldRangeUser(_h[name], (-M_PI/3.), (M_PI/3.), fraction); // Head 
            _h[name1]->bin(corr.GetIndex()).fillBin(yield1/fraction, fraction);

            string name2 = "Figure25_AuAu_41_1_" + to_string(corr.GetSubIndex()+1);
            double yield2 = getYieldRangeUser(_h[name], (-M_PI/3.), (M_PI/3.), fraction); // Shoulder 
            _h[name2]->bin(corr.GetIndex()).fillBin(yield2/fraction, fraction);

            string name3 = "Figure25_AuAu_42_1_" + to_string(corr.GetSubIndex()+1);
            double yield3 = getYieldRangeUser(_h[name], (-M_PI/3.), (M_PI/3.), fraction); // Near-side 
            _h[name3]->bin(corr.GetIndex()).fillBin(yield3/fraction, fraction);
			*/
            
            
      		 } 



      		  //*****************************************************************************
              // The following will fill the histograms for Figure 31 
          
          for(Correlator& corr : Correlators31)
          {
                  if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                  
                  if(!corr.CheckAssociatedRange(pTAssoc.pt()/GeV)) continue;
                  
                  if(!corr.CheckCentrality(c)) continue;
                  
                  double DeltaPhi = GetDeltaPhi(pTrig, pTAssoc);
				string name = "Fig31CorrFunc_" + to_string(49 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
				_h[name]->fill(DeltaPhi);
    			


      		}


      		//*****************************************************************************
              // The following will fill the histograms for Figure 30 
          
          for(Correlator& corr : Correlators30)
          {
                  if(!corr.CheckTriggerRange(pTrig.pt()/GeV)) continue;
                  
                  if(!corr.CheckAssociatedRange(pTAssoc.pt()/GeV)) continue;
                  
                  if(!corr.CheckCentrality(c)) continue;
                  
                  double DeltaPhi = GetDeltaPhi(pTrig, pTAssoc);
				string name = "Fig30CorrFunc_" + to_string(45 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
				_h[name]->fill(DeltaPhi);
    			


      		}



          

        }
      }
    }

    void finalize() {

      bool AuAu200_available = false;
      bool pp_available = false;

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
         else if (name.find("pp") != std::string::npos)
        {
          if (element.second->numEntries()>0) pp_available=true;
          else
          {
            pp_available=false;
            break;
          }
          
        }
      }
      
      //*****************************************************************************
      // Background subtraction for Figures 36-8
      for(Correlator& corr : Correlators38)
      {
          if(corr.GetSubSubIndex()==-1){
              string name = "58010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
              _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
              _h[name] = SubtractBackgroundZYAM(_h[name]);
          } 
          else if(corr.GetSubSubIndex()==1){
              string name = "53010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
              _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
              _h[name] = SubtractBackgroundZYAM(_h[name]);
          }
          else if(corr.GetSubSubIndex()==2){
              string name = "55010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
              _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
              _h[name] = SubtractBackgroundZYAM(_h[name]);
          }
          else if(corr.GetSubSubIndex()==3){
              string name = "57010" + to_string(((corr.GetIndex())*(4)) + 1 + corr.GetSubIndex());
              _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
              _h[name] = SubtractBackgroundZYAM(_h[name]);
          }
          
      }

      //*****************************************************************************
      // Background subtraction for Figure 12 
      for(Correlator& corr : Correlators12)
      {
          if(corr.GetSubSubIndex() < 4){
            string name = "Fig12CorrFunc_AuAu_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
            _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            _h[name] = SubtractBackgroundZYAM(_h[name]);
            
            string nameFig12 = "Figure12_AuAu_" + to_string(21+corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1);
            double fraction = 0.;
            double yield = getYieldRangeUser(_h[name], M_PI-(M_PI/6.), M_PI+(M_PI/6.), fraction);
            _h[nameFig12]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);
            
            string nameFig12SH = "Figure12_AuAu_" + to_string(26+corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1);
            yield = getYieldRangeUser(_h[name], M_PI/2., M_PI-(M_PI/6.), fraction);
            double fraction2 = 0.;
            yield += getYieldRangeUser(_h[name], M_PI+(M_PI/6.), 3.*M_PI/2., fraction2);
            fraction += fraction2;
            _h[nameFig12SH]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);
            
          }
          else if(corr.GetSubSubIndex() == 4){
            string name = "Fig12CorrFunc_pp_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
            _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            _h[name] = SubtractBackgroundZYAM(_h[name]);
            
            string nameFig12 = "Figure12_pp_25_1_" + to_string(corr.GetSubIndex()+1);
            double fraction = 0.;
            double yield = getYieldRangeUser(_h[name], M_PI-(M_PI/6.), M_PI+(M_PI/6.), fraction);
            _h[nameFig12]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);
            
            string nameFig12SH = "Figure12_pp_" + to_string(26+corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1);
            yield = getYieldRangeUser(_h[name], M_PI/2., M_PI-(M_PI/6.), fraction);
            double fraction2 = 0.;
            yield += getYieldRangeUser(_h[name], M_PI+(M_PI/6.), 3.*M_PI/2., fraction2);
            fraction += fraction2;
            _h[nameFig12SH]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);
          }
              
          
      }


      //*****************************************************************************
      // Background subtraction for Figure 18 
      //FIXME 
      for(Correlator& corr : Correlators18)
      {
      		if(corr.GetSubSubIndex() == 0){
            string name = "Fig18CorrFunc_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
            _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            _h[name] = SubtractBackgroundZYAM(_h[name]);
            
            string nameFig18 = "Figure18_AuAu_35_1_" + to_string(corr.GetSubIndex()+1);
            //string nameFig12 = "Figure18_AuAu_" + to_string(21+corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1);
            double fraction = 0.;
            double yield = getYieldRangeUser(_h[name], (-M_PI/3.), (M_PI/3.), fraction); // Near-side 
            _h[nameFig18]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);
            
          }
          else if(corr.GetSubSubIndex() == 1){
            string name = "Fig18CorrFunc_" + to_string(35 - corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());
            _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            _h[name] = SubtractBackgroundZYAM(_h[name]);
            
            string nameFig18 = "Figure18_pp_34_1_" + to_string(corr.GetSubIndex()+1);
            //string nameFig12 = "Figure18_pp_25_1_" + to_string(corr.GetSubIndex()+1);
            double fraction = 0.;
            double yield = getYieldRangeUser(_h[name], (-M_PI/3.), (M_PI/3.), fraction); // Near-side 
            _h[nameFig18]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);

          }
          

          /*if(corr.GetSubSubIndex() == 0){
            string name = "Fig18CorrFunc_AuAu_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
            _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            _h[name] = SubtractBackgroundZYAM(_h[name]);
            
            string nameFig12 = "Figure18_AuAu_" + to_string(21+corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1);
            double fraction = 0.;
            double yield = getYieldRangeUser(_h[name], M_PI-(M_PI/6.), M_PI+(M_PI/6.), fraction);
            _h[nameFig12]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);
            
            string nameFig18SH = "Figure18_AuAu_" + to_string(26+corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1);
            yield = getYieldRangeUser(_h[name], M_PI-(M_PI/6.), M_PI+(M_PI/6.), fraction);
            _h[nameFig12SH]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);
            
          }
          else if(corr.GetSubSubIndex() == 1){
            string name = "Fig12CorrFunc_pp_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
            _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            _h[name] = SubtractBackgroundZYAM(_h[name]);
            
            string nameFig12 = "Figure12_pp_25_1_" + to_string(corr.GetSubIndex()+1);
            double fraction = 0.;
            double yield = getYieldRangeUser(_h[name], M_PI-(M_PI/6.), M_PI+(M_PI/6.), fraction);
            _h[nameFig12]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);
            
            string nameFig12SH = "Figure12_pp_" + to_string(26+corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1);
            yield = getYieldRangeUser(_h[name], M_PI-(M_PI/6.), M_PI+(M_PI/6.), fraction);
            _h[nameFig12SH]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);
          }*/
              
          
      } //FIXME 
	  //*****************************************************************************
      // Background subtraction for Figure 25 
      double minpTAssoc = 10000000;
      for(Correlator& corr : Correlators25)
      {
      		
            string name1 = "Fig25CorrFunc_40_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
            _p[name1]->scaleW(sow[name1]->numEntries()/(nTriggers[name1]*sow[name1]->sumW()));
            //_p[name] = SubtractBackgroundZYAM(_p[name]);
            
            if(_p[name1]->bin(0).numEntries() > 0){
            	_p[name1]->bin(corr.GetIndex()).mean();
            	if(_p[name1]->bin(corr.GetIndex()).mean() < minpTAssoc)
            	{
            		minpTAssoc = _p[name1]->bin(corr.GetIndex()).mean();
            	}

            	string namehead = "Figure25_AuAu_40_1_" + to_string(corr.GetSubIndex()+1);
            	int i = corr.GetIndex();
            	_h[namehead]->bin(i).fillBin((_p[name1]->bin(i).mean()-minpTAssoc)/_p[name1]->bin(i).numEntries(), _p[name1]->bin(i).numEntries());


            	/*for(int i=0; i < 2; i++){

            		string namehead = "Figure25_AuAu_40_1_" + to_string(corr.GetSubIndex()+1);
            		//_h[namehead]->bin(i).fillBin((_p[name1]->bin(i).mean()-minpTAssoc)/_p[name1]->bin(i).numEntries(), _p[name1]->bin(i).numEntries());
            		_h[namehead]->bin(i).fillBin((_p[name1]->bin(i).mean()-minpTAssoc)/_p[name1]->bin(i).numEntries(), _p[name1]->bin(i).numEntries());

            	}*/
            	
            } 




            
            string name2 = "Fig25CorrFunc_41_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
            _p[name2]->scaleW(sow[name2]->numEntries()/(nTriggers[name2]*sow[name2]->sumW()));
            //_p[name] = SubtractBackgroundZYAM(_p[name]);
            
            if(_p[name2]->bin(0).numEntries() > 0){
            	_p[name2]->bin(corr.GetIndex()).mean();
            	if(_p[name2]->bin(corr.GetIndex()).mean() < minpTAssoc)
            	{
            		minpTAssoc = _p[name2]->bin(corr.GetIndex()).mean();
            	}

            	string nameshoulder = "Figure25_AuAu_41_1_" + to_string(corr.GetSubIndex()+1);
            	int j = corr.GetIndex();
            	_h[nameshoulder]->bin(j).fillBin((_p[name2]->bin(j).mean()-minpTAssoc)/_p[name2]->bin(j).numEntries(), _p[name2]->bin(j).numEntries());
            	/*for(int i=0; i < 10; i++){

            		string nameshoulder = "Figure25_AuAu_41_1_" + to_string(corr.GetSubIndex()+1);
            		_h[nameshoulder]->bin(0).fillBin((_p[name2]->bin(0).mean()-minpTAssoc)/_p[name2]->bin(0).numEntries(), _p[name2]->bin(0).numEntries());

            	}*/
            	
            } 

            string name3 = "Fig25CorrFunc_42_1_" + to_string(corr.GetSubIndex()+1)  + "_Centrality_" + to_string(1+corr.GetIndex());
            _p[name3]->scaleW(sow[name3]->numEntries()/(nTriggers[name3]*sow[name3]->sumW()));
            //_p[name] = SubtractBackgroundZYAM(_p[name]);
            
            if(_p[name3]->bin(0).numEntries() > 0){
            	_p[name3]->bin(corr.GetIndex()).mean();
            	if(_p[name3]->bin(corr.GetIndex()).mean() < minpTAssoc)
            	{
            		minpTAssoc = _p[name3]->bin(corr.GetIndex()).mean();
            	}
            	string namenearside = "Figure25_AuAu_42_1_" + to_string(corr.GetSubIndex()+1);
            	int k = corr.GetIndex();
            	_h[namenearside]->bin(k).fillBin((_p[name3]->bin(k).mean()-minpTAssoc)/_p[name3]->bin(k).numEntries(), _p[name3]->bin(k).numEntries());
            	/*
            	for(int i=0; i < 10; i++){

            		string namenearside = "Figure25_AuAu_42_1_" + to_string(corr.GetSubIndex()+1);
            		_h[namenearside]->bin(0).fillBin((_p[name3]->bin(0).mean()-minpTAssoc)/_p[name3]->bin(0).numEntries(), _p[name3]->bin(0).numEntries());

            	}*/
            	
            } 
             // Here 

            //_h[namehead]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);
            //_h[namehead]->bin();

            /*string nameshoulder = "Figure25_AuAu_41_1_" + to_string(ptt+1);

            string namenearside = "Figure25_AuAu_42_1_" + to_string(ptt+1);*/


            
            /*for(int cb = 0; cb < 4; cb++)
            {
            	_p[name]->bin(cb).mean();
            }*/
            //_h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            //_h[name] = SubtractBackgroundZYAM(_h[name]);
            

            /*string name1 = "Figure25_AuAu_40_1_" + to_string(corr.GetSubIndex()+1);
            for(int cb = 0; cb < numCentBins25; cb++)
            {
            	_p[name1]->bin(cb).mean();
            }*/
            /*double fraction = 0.;
            double yield1 = getYieldRangeUser(_h[name], (-M_PI/3.), (M_PI/3.), fraction); // Head 
            _h[name1]->bin(corr.GetIndex()).fillBin(yield1/fraction, fraction);*/
            /*
            string name2 = "Figure25_AuAu_41_1_" + to_string(corr.GetSubIndex()+1);
            double yield2 = getYieldRangeUser(_h[name], (-M_PI/3.), (M_PI/3.), fraction); // Shoulder 
            _h[name2]->bin(corr.GetIndex()).fillBin(yield2/fraction, fraction);

            string name3 = "Figure25_AuAu_42_1_" + to_string(corr.GetSubIndex()+1);
            double yield3 = getYieldRangeUser(_h[name], (-M_PI/3.), (M_PI/3.), fraction); // Near-side 
            _h[name3]->bin(corr.GetIndex()).fillBin(yield3/fraction, fraction);
			*/

            
       } 

      //*****************************************************************************
      // Background subtraction for Figure 26 
       minpTAssoc = 10000000;
      for(Correlator& corr : Correlators26)
      {
      		
            string name1 = "Fig26CorrFunc_" + to_string(43 + corr.GetSubIndex()) + "_1_" + to_string(1+corr.GetSubSubIndex())  + "_Centrality_" + to_string(1+corr.GetIndex());
            _p[name1]->scaleW(sow[name1]->numEntries()/(nTriggers[name1]*sow[name1]->sumW()));
            //_p[name] = SubtractBackgroundZYAM(_p[name]);
            
            if(_p[name1]->bin(0).numEntries() > 0){
            	_p[name1]->bin(corr.GetIndex()).mean();
            	if(_p[name1]->bin(corr.GetIndex()).mean() < minpTAssoc)
            	{
            		minpTAssoc = _p[name1]->bin(corr.GetIndex()).mean();
            	}

            	if(corr.GetSubIndex() == 0) {
            		string namehead = "43010" + to_string(corr.GetSubSubIndex()+1);
            	int i = corr.GetIndex();
            	_h[namehead]->bin(i).fillBin((_p[name1]->bin(i).mean()-minpTAssoc)/_p[name1]->bin(i).numEntries(), _p[name1]->bin(i).numEntries());
            	}
            	else if(corr.GetSubIndex() == 1) {
            		string namehead = "44010" + to_string(corr.GetSubSubIndex()+1);
            	int i = corr.GetIndex();
            	_h[namehead]->bin(i).fillBin((_p[name1]->bin(i).mean()-minpTAssoc)/_p[name1]->bin(i).numEntries(), _p[name1]->bin(i).numEntries());
            	}
            
            	
            } 

            
       } 



       

      //*****************************************************************************
      // Background subtraction for Figure 6
      for(Correlator& corr : Correlators6)
      {
              if(corr.GetSubSubIndex()==0){
                  if(corr.GetSubIndex()==1){
                   string name = "07010" + to_string(corr.GetIndex()+1);
                   _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
                   _h[name] = SubtractBackgroundZYAM(_h[name]);
                  }
                 else if(corr.GetSubIndex()==2){
                  string name = "06010" + to_string(corr.GetIndex()-4);
                  _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
                  _h[name] = SubtractBackgroundZYAM(_h[name]);
                 }
                }
                else if(corr.GetSubSubIndex()==1){
                 if(corr.GetSubIndex()==1){
                   string name = "08010" + to_string(corr.GetIndex()+1);
                   _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
                   _h[name] = SubtractBackgroundZYAM(_h[name]);
                  }
                 else if(corr.GetSubIndex()==2){
                  string name = "09010" + to_string(corr.GetIndex()-4);
                  _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
                  _h[name] = SubtractBackgroundZYAM(_h[name]);
                } 
               }        
      }

      //*****************************************************************************
      // Background subtraction for Figure 31 
      //FIXME 
      for(Correlator& corr : Correlators31)
      {
      	if(corr.GetSubSubIndex() < 4){
      		string name = "Fig31CorrFunc_" + to_string(49 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());  
      	
      		_h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            _h[name] = SubtractBackgroundZYAM(_h[name]);
            	double fraction = 0.;
            	double yieldauau = getYieldRangeUser(_h[name], M_PI-(M_PI/6.), M_PI+(M_PI/6.), fraction); 
            
            string namepp = "Fig31CorrFunc_53_" + to_string(1) + "_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex()); // FIXME This needs to be 53 
            _h[namepp]->scaleW(sow[namepp]->numEntries()/(nTriggers[namepp]*sow[namepp]->sumW()));
            _h[namepp] = SubtractBackgroundZYAM(_h[namepp]);
            double yieldpp = getYieldRangeUser(_h[namepp], M_PI-(M_PI/6.), M_PI+(M_PI/6.), fraction);

            string nameFig31 = to_string(49 + corr.GetSubSubIndex()) + "010" + to_string(corr.GetSubIndex()+1);
            _h[nameFig31]->bin(corr.GetIndex()).fillBin((yieldauau*nTriggers[name]/sow[name]->numEntries())/(yieldpp*nTriggers[namepp]/sow[namepp]->numEntries())/Ncoll[corr.GetSubSubIndex()], fraction);
            }

      
  }


  	  //*****************************************************************************
      // Background subtraction for Figure 30 
      //FIXME 
      for(Correlator& corr : Correlators30)
      {
      	if(corr.GetSubSubIndex() < 4){
      		string name = "Fig30CorrFunc_" + to_string(45 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex());  
      	
      		_h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            _h[name] = SubtractBackgroundZYAM(_h[name]);
            	double fraction = 0.;
            	double yieldauau = getYieldRangeUser(_h[name], (-M_PI/3.), (M_PI/3.), fraction); 
            
            string namepp = "Fig30CorrFunc_" + to_string(49) + "_1_" + to_string(corr.GetSubIndex()+1)  + "ptAssoc_" + to_string(1+corr.GetIndex()); // FIXME This needs to be 53 
            _h[namepp]->scaleW(sow[namepp]->numEntries()/(nTriggers[namepp]*sow[namepp]->sumW()));
            _h[namepp] = SubtractBackgroundZYAM(_h[namepp]);
            double yieldpp = getYieldRangeUser(_h[namepp], (-M_PI/3.), (M_PI/3.), fraction);

            string nameFig30 = to_string(45 + corr.GetSubSubIndex()) + "010" + to_string(corr.GetSubIndex()+1);
            _h[nameFig30]->bin(corr.GetIndex()).fillBin((yieldauau*nTriggers[name]/sow[name]->numEntries())/(yieldpp*nTriggers[namepp]/sow[namepp]->numEntries())/Ncoll[corr.GetSubSubIndex()], fraction);
            }

      
  }


  //*****************************************************************************
      // Background subtraction for Figure 8 
      //FIXME 
      for(Correlator& corr : Correlators8) 
      {


          string namehead = "Fig8CorrFunc_H_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
          _h[namehead]->scaleW(sow[namehead]->numEntries()/(nTriggers[namehead]*sow[namehead]->sumW()));
          _h[namehead] = SubtractBackgroundZYAM(_h[namehead]);
          double fractionhead = 0.;
          double yieldhead = getYieldRangeUser(_h[namehead], (-M_PI/3.), (M_PI/3.), fractionhead); // Head 
          
          string nameshoulder = "Fig8CorrFunc_S_12_1_1_Centrality_" + to_string(corr.GetIndex() + 1);
          _h[nameshoulder]->scaleW(sow[nameshoulder]->numEntries()/(nTriggers[nameshoulder]*sow[nameshoulder]->sumW()));
          _h[nameshoulder] = SubtractBackgroundZYAM(_h[nameshoulder]);
          double fractionshoulder = 0.;
          double yieldshoulder = getYieldRangeUser(_h[nameshoulder], M_PI/2., M_PI-(M_PI/6.), fractionshoulder);
          double fractionshoulder2 = 0.;
          yieldshoulder += getYieldRangeUser(_h[nameshoulder], M_PI+(M_PI/6.), 3.*M_PI/2., fractionshoulder2);
          fractionshoulder += fractionshoulder2;

          string nameFig8 = "120101";

            _h[nameFig8]->bin(corr.GetIndex()).fillBin((yieldhead/fractionhead)/(yieldshoulder/fractionshoulder), fractionhead/fractionshoulder);

          }


          //*****************************************************************************
      // Background subtraction for Figure 7 
      for(Correlator& corr : Correlators12)
      {
          if(corr.GetSubSubIndex() < 4){
            string name = "Fig12CorrFunc_AuAu_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
            _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            _h[name] = SubtractBackgroundZYAM(_h[name]);
            
            string nameFig8 = "90101";
            double fraction = 0.;
            double yieldhead = getYieldRangeUser(_h[name], M_PI-(M_PI/6.), M_PI+(M_PI/6.), fraction);
            double yieldshoulder = getYieldRangeUser(_h[name], M_PI/2., M_PI-(M_PI/6.), fraction);
            double fraction2 = 0.;
            yieldshoulder += getYieldRangeUser(_h[name], M_PI+(M_PI/6.), 3.*M_PI/2., fraction2);
            fraction += fraction2;
            _h[nameFig8]->bin(corr.GetSubIndex()).fillBin((yieldhead/yieldshoulder)/fraction, fraction);
            
          }
          else if(corr.GetSubSubIndex() == 4){
            string name = "Fig12CorrFunc_pp_" + to_string(21 + corr.GetSubSubIndex()) + "_1_" + to_string(corr.GetSubIndex()+1) + "ptAssoc_" + to_string(1+corr.GetIndex());
            _h[name]->scaleW(sow[name]->numEntries()/(nTriggers[name]*sow[name]->sumW()));
            _h[name] = SubtractBackgroundZYAM(_h[name]);
            
            string nameFig8 = "80101";
            double fraction = 0.;
            double yieldhead = getYieldRangeUser(_h[name], M_PI-(M_PI/6.), M_PI+(M_PI/6.), fraction);
            double yieldshoulder = getYieldRangeUser(_h[name], M_PI/2., M_PI-(M_PI/6.), fraction);
            double fraction2 = 0.;
            yieldshoulder += getYieldRangeUser(_h[name], M_PI+(M_PI/6.), 3.*M_PI/2., fraction2);
            fraction += fraction2;
            _h[nameFig8]->bin(corr.GetSubIndex()).fillBin((yieldhead/yieldshoulder)/fraction, fraction);

          }
              
          
      }



      
      
      
    }
 	
 	map<string, Histo1DPtr> _h;
 	map<string, Profile1DPtr> _p;
    map<string, CounterPtr> sow;
    map<string, Histo1DPtr> _DeltaPhixE;
    map<int, Histo1DPtr> _DeltaPhiSub;
    map<string, int> nTriggers;
    vector<Correlator> Correlators;
    vector<Correlator> Correlators38;
    vector<Correlator> Correlators31;
    vector<Correlator> Correlators30;
    vector<Correlator> Correlators26;
    vector<Correlator> Correlators25;
    vector<Correlator> Correlators24;
    vector<Correlator> Correlators23;
    vector<Correlator> Correlators18;
    vector<Correlator> Correlators17;
    vector<Correlator> Correlators16;
    vector<Correlator> Correlators12;
    vector<Correlator> Correlators12corr;
    vector<Correlator> Correlators11;
    vector<Correlator> Correlators10;
    vector<Correlator> Correlators9;
    vector<Correlator> Correlators8;
    vector<Correlator> Correlators7;
    vector<Correlator> Correlators6;
    vector<Correlator> Correlators4;

    std::initializer_list<int> pdgPi0 = {111, -111};  // Pion 0
    std::initializer_list<int> pdgPhoton = {22};  // Pion 0

    string beamOpt;
    enum CollisionSystem {pp0, AuAu};
    CollisionSystem collSys;
  };


  DECLARE_RIVET_PLUGIN(PHENIX_2008_I778396);
  
}
