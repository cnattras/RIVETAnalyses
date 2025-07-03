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
static const int numAssocPtBins = 3;
static const float pTAssocBins[] = {3.0,4.0,6.0,10};
static const int numzTBins = 7;
static const float zTBins[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
static const int numCentBins = 6;
static const float centBins[] = {0.0,5.0,10,20,30,40,80};

using namespace std;
namespace Rivet {

    
  class Correlator {
      
    private:
		std::vector<int> _indices; //added
      int _index;
      int _subindex;
      string _collSystemAndEnergy;
      pair<double,double> _centrality;
      pair<double,double> _triggerRange;
      pair<double,double> _associatedRange;
      pair<double,double> _zTRange;
      vector<int> _pid;
	  Histo1DPtr _deltaPhi; //added
	  CounterPtr _counter; //added
	  CounterPtr _cTriggers; //added

    public:
    
      /// Constructor
      // Correlator(int index, int subindex) {
      //   _index = index;
      //   _subindex = subindex;
      // }

	  Correlator(int index0, int index1, int index2) {
		  _indices = { index0, index1, index2 };
	  } //added

	  Correlator(int index0, int index1) {
		  _indices = { index0, index1 };
	  }

	  Correlator(int index0) {
		  _indices = { index0 };
	  } //added

	  Correlator(std::vector<int> vindex) {
		  _indices = vindex;
	  } //added





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
            string GetFullIndex()
      {
              string fullIndex = "";
              for(int index : _indices)
              {
                      fullIndex += to_string(index);
              }

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
                binWidth = bin.xWidth();
                minValueEntries = bin.numEntries();
            }
            if(bin.sumW()/bin.xWidth() < minValue/binWidth)
            {
                minValue = bin.sumW();
                binWidth = bin.xWidth();
                minValueEntries = bin.numEntries();
            }
        }
                
        hist.reset();
        
        for(auto &bin : hist.bins())
        {
            hist.fill(bin.xWidth(), (minValue*bin.xWidth())/(minValueEntries*binWidth));
        }
        
        *histo -= hist;
        
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
                integral += bin.sumW()/bin.xWidth();
                fraction += bin.numEntries()/bin.xWidth();
            }
            else if((bin.xMin() < xmin) && (bin.xMax() > xmin))
            {
                double perc = (bin.xMax() - xmin)/bin.xWidth();
                integral += perc*(bin.sumW()/bin.xWidth());
                fraction += perc*(bin.numEntries()/bin.xWidth());
                
            }
            else if((bin.xMin() < xmax) && (bin.xMax() > xmax))
            {
                double perc = (xmax - bin.xMin())/bin.xWidth();
                integral += perc*(bin.sumW()/bin.xWidth());
                fraction += perc*(bin.numEntries()/bin.xWidth());
            }
        }
        
        return integral;
        
    }
   
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(STAR_2006_I715470);
    void init() {
      const ChargedFinalState cfs(Cuts::abseta < 0.35);
      declare(cfs, "CFS");

      beamOpt = getOption<string>("beam", "NONE");

      if (beamOpt == "AUAU200") collSys = AuAu200;
      else if (beamOpt == "dAU200") collSys = dAu200;
      
      declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

int ndPhiBins = 72*4;
double lowedge = -pi/2.0;
double highedge = 3.0*pi/2.0;


        for(int ncb=0;ncb<numCentBins;ncb++){
          for(int nassoc=0;nassoc<numAssocPtBins;nassoc++){
      //Correlators
            Correlator c1(ncb,nassoc);
            c1.SetCollSystemAndEnergy("AuAu200");
            c1.SetCentrality(centBins[ncb], centBins[ncb+1]);
            c1.SetTriggerRange(8., 15.);
            c1.SetAssociatedRange(pTAssocBins[nassoc],pTAssocBins[nassoc+1]);
            c1.SetzTRange(0,1);
            Correlators.push_back(c1);
            string name = c1.GetCollSystemAndEnergy()+c1.GetFullIndex();
            //string name3 = c1.GetCollSystemAndEnergy()+c1.GetFullIndex()+"BkgdSubtracted";
            book(_h[name], name, ndPhiBins,lowedge,highedge);
            cout<<"395 making "<<name<<endl;
            //book(_h[name3], name3, ndPhiBins,lowedge,highedge);
            book(sow[c1.GetFullIndex()],"sow" + c1.GetFullIndex());
            if(ncb==0){

              Correlator c2(100,nassoc);
              c2.SetCollSystemAndEnergy("dAu200");
              c2.SetCentrality(0,80);
              c2.SetTriggerRange(8., 15.);
              c2.SetAssociatedRange(pTAssocBins[nassoc],pTAssocBins[nassoc+1]);
              c2.SetzTRange(0,1); 
              Correlators.push_back(c2);
            string name2 = c2.GetCollSystemAndEnergy()+c2.GetFullIndex();
            //string name4 = c2.GetCollSystemAndEnergy()+c2.GetFullIndex()+"BkgdSubtracted";
            book(_h[name2], name2, ndPhiBins,lowedge,highedge);
            //book(_h[name4], name4, ndPhiBins,lowedge,highedge);
            book(sow[c2.GetFullIndex()],"sow" + c2.GetFullIndex());
            cout<<"412 making* "<<name2<<endl;
            }
          }



          for(int nzT=0;nzT<numzTBins;nzT++){
cout<<"This loop repeats some names for histograms which are made above.  I did not find the logic error.  6 Aug 2021 CN"<<endl;
      //Correlators
            Correlator c1(ncb,1000);
            c1.SetCollSystemAndEnergy("AuAu200");
            c1.SetCentrality(centBins[ncb], centBins[ncb+1]);
            c1.SetTriggerRange(8,15);
            c1.SetzTRange(zTBins[nzT],zTBins[nzT+1]);
            c1.SetAssociatedRange(3,15);
            Correlators.push_back(c1);
            string name = c1.GetCollSystemAndEnergy()+c1.GetFullIndex();
            //string name3 = c1.GetCollSystemAndEnergy()+c1.GetFullIndex()+"BkgdSubtracted";
            book(_h[name], name, ndPhiBins,lowedge,highedge);
            //book(_h[name3], name3, ndPhiBins,lowedge,highedge);
            book(sow[c1.GetFullIndex()],"sow" + c1.GetFullIndex());
            cout<<"433 making* "<<name<<endl;
            if(ncb==0){

              Correlator c2(1000,1000);
              c2.SetCollSystemAndEnergy("dAu200");
              c2.SetCentrality(0,80);
              c2.SetTriggerRange(8,15);
              c2.SetzTRange(zTBins[nzT],zTBins[nzT+1]);
              c2.SetAssociatedRange(3,15);
              Correlators.push_back(c2);
            string name2 = c2.GetCollSystemAndEnergy()+c2.GetFullIndex();
            //string name4 = c2.GetCollSystemAndEnergy()+c2.GetFullIndex()+"BkgdSubtracted";
            book(_h[name2], name2, ndPhiBins,lowedge,highedge);
            //book(_h[name4], name4, ndPhiBins,lowedge,highedge);
            book(sow[c2.GetFullIndex()],"sow" + c2.GetFullIndex());
            cout<<"448 making* "<<name2<<endl;
            }
          }//End of the zT loop

      }//end loop over centrality bins

        //Fig 3 - yield per trigger vs centrality
        book(_h["010101"], 1, 1, 1);//pTassoc = 3-4 NS
        book(_h["010102"], 1, 1, 2);//pTassoc = 4-6 NS
        book(_h["010103"], 1, 1, 3);//pTassoc > 6 NS
        book(_h["020101"], 2, 1, 1);//pTassoc = 3-4 AS
        book(_h["020102"], 2, 1, 2);//pTassoc = 4-6 AS 
        book(_h["020103"], 2, 1, 3);//pTassoc > 6 AS

        //Fig 4 - D(z) vs z_T
        book(_h["030101"], 3, 1, 1);//d+Au NS   
		book(_h["040101"], 4, 1, 1);//Au+Au 0-5% NS
		book(_h["040102"], 4, 1, 2);//Au+Au 20-40% NS *can't fill yet

        book(_h["050101"], 5, 1, 1);//d+Au AS
        book(_h["050102"], 5, 1, 2);//0-5% AS
		book(_h["050103"], 5, 1, 3);//Au+Au 20-40% AS *can't fill yet
    }
    void analyze(const Event& event) {
        
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      
      //==================================================
      // Select the histograms accordingly to the collision system, beam energy and centrality
      // WARNING: Still not implemented for d-Au
      //==================================================
      string SysAndEnergy = "";

      if(collSys == AuAu200) SysAndEnergy = "AuAu200";
      else if(collSys == dAu200) SysAndEnergy = "dAu200";

      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();
    
      double triggerptMin = 999.;
      double triggerptMax = -999.;
      double associatedptMin = 999.;
      double associatedptMax = -999.;
    
      bool isVeto = true;
    
      for(Correlator& corr : Correlators)
      {
        if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
        if(!corr.CheckCentrality(c)) continue;
        
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
            if(!corr.CheckCentrality(c)) continue;//May need to check that the centrality for d+Au matches the centrality in the correlator
            if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
            if(!corr.CheckTriggerRange(pTrig.pT()/GeV)) continue;  
            nTriggers[corr.GetFullIndex()]++;//Antonio is this OK? Nora uses strings.
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
            //double dPhi = deltaPhi(pTrig, pAssoc, true);//this does NOT rotate the delta phi to be in a given range
            double dPhi = GetDeltaPhi(pTrig, pAssoc);
                        
            //double xE = GetXE(pTrig,pAssoc);
            
            for(Correlator& corr : Correlators)
            {
                if(!corr.CheckCentrality(c)) continue;
                if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
                if(!corr.CheckAssociatedRange(pTrig.pT()/GeV)) continue;
                
                
                
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
        for(Correlator& corr : Correlators)
        {
            
            string name = corr.GetCollSystemAndEnergy()+corr.GetFullIndex();
            //AJ add some cross checks to make sure you do not divide by zero here.

            if(nTriggers[corr.GetFullIndex()] > 0) _h[name]->scaleW(sow[corr.GetFullIndex()]->numEntries()/(nTriggers[corr.GetFullIndex()]*sow[corr.GetFullIndex()]->sumW()));
            //string name2 = corr.GetCollSystemAndEnergy()+corr.GetFullIndex()+"BkgdSubtracted";
                  _h[name] = SubtractBackgroundZYAM(_h[name]);
                  //You will need something like this

  //            double fraction = 0.;
                  //AJ work on this
                  //1.  Just calculate the yield but adjust the range so that it matched the paper.
                  //2.  Try to fill the histogram as below, except you'll need to figure out what the name is.
      //      _h[nameFig12]->bin(corr.GetIndex()).fillBin(yield/fraction, fraction);
            //  _h[nameFig8]->bin(corr.GetSubIndex()).fillBin((yieldhead/yieldshoulder)/fraction, fraction);
            //Histogram 1: all correlators with pTassoc = 3-4 GeV, trigger 8-15


if(corr.CheckConditions("AuAu200",2.5)){//Central Au+Au
//cout<<"I am AuAu200GeV 0-5%"<<endl;
}
if(corr.CheckConditions("AuAu200",35)){//Central Au+Au
//cout<<"I am AuAu200GeV 20-40%"<<endl;
}
if(corr.CheckConditions("dAu200",2.5)){//Central Au+Au
//cout<<"I am dAu200GeV"<<endl;
}
        }
        
        
    }
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> sow;
    map<string, int> nTriggers;
    vector<Correlator> Correlators;
    
    string beamOpt;
    enum CollisionSystem {AuAu200, dAu200};
    CollisionSystem collSys;


  };


  RIVET_DECLARE_PLUGIN(STAR_2006_I715470);

}
