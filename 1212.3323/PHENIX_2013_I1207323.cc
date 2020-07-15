
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

  private:

    int _index;
    int _subindex;
    string _collSystemAndEnergy;
    pair<double,double> _centrality;
    pair<double,double> _triggerRange;
    pair<double,double> _associatedRange;
    pair<double,double> _xiRange;
  
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
    //Instead of directly taking the calculated Delta_eta, this method takes the center of the bin associated to the Delta_eta.
    //This is done to avoid values of zero efficiency when |Delta_eta| is close to maxDeltaEta
    double EtaEffCorrection(double deltaEta, YODA::Histo1D& hist)
    {
        double maxDeltaEta = 2.;
        
        int binEta = hist.binIndexAt(deltaEta);
        
        if(binEta < 0)
        {
            MSG_INFO("Eta is out of bounds!");
            return 0.;
        }
        
        double binCenterEta = hist.bin(binEta).xMid();
        
        return 1. - (1./maxDeltaEta)*abs(binCenterEta);
    }
    
    //bmin and bmax are included in the integral. Range = [bmin, bmax]
    //Give the min and max bins to calculate the integral
    double GetYieldInBinRange(YODA::Histo1D& hist, int bmin, int bmax)
    {
        double integral = 0.;
        
        if(bmin < 0 || bmax >= (int)hist.numBins())
        {
            MSG_ERROR("Out of range!");
            return 0.;
        }
        
        for(int i = bmin; i <= bmax; i++)
        {
            integral += hist.bin(i).sumW();
        }
        
        return integral;
        
    }
    
    //vmin is included in the integral, but vmax is not. Range = [vmin, vmax[
    //Give the min and max values to calculate the integral
    double GetYieldInUserRange(YODA::Histo1D& hist, double vmin, double vmax)
    {
        double integral = 0.;
        
        if(hist.binIndexAt(vmin) < 0 || hist.binIndexAt(vmax) < 0)
        {
            MSG_ERROR("Out of range!");
            return 0.;
        }
        
        for(int i = hist.binIndexAt(vmin); i < hist.binIndexAt(vmax); i++)
        {
            integral += hist.bin(i).sumW();
        }
        
        return integral;
        
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
      
      Correlator c1(1,1);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 40.);
      c1.SetXiRange(0., 0.4);
      c1.SetTriggerRange(5., 9.);
      c1.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c1);
      
      Correlator c2(1,2);
      c2.SetCollSystemAndEnergy("pp200GeV");
      c2.SetCentrality(0., 40.);
      c2.SetXiRange(0., 0.4);
      c2.SetTriggerRange(5., 9.);
      c2.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c2);
      
      Correlator c3(1,3);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0., 40.);
      c3.SetXiRange(0.4, 0.8);
      c3.SetTriggerRange(5., 9.);
      c3.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c3);
      
      Correlator c4(1,4);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetCentrality(0., 40.);
      c4.SetXiRange(0.4, 0.8);
      c4.SetTriggerRange(5., 9.);
      c4.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c4);
      
      Correlator c5(1,5);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0., 40.);
      c3.SetXiRange(0.8, 1.2);
      c3.SetTriggerRange(5., 9.);
      c3.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c5);
      
      Correlator c6(1,6);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetCentrality(0., 40.);
      c4.SetXiRange(0.8, 1.2);
      c4.SetTriggerRange(5., 9.);
      c4.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c6);

      Correlator c7(1,7);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 40.);
      c1.SetXiRange(1.2, 1.6);
      c1.SetTriggerRange(5., 9.);
      c1.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c7);
      
      Correlator c8(1,8);
      c2.SetCollSystemAndEnergy("pp200GeV");
      c2.SetCentrality(0., 40.);
      c2.SetXiRange(1.2, 1.6);
      c2.SetTriggerRange(5., 9.);
      c2.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c8);
      
      Correlator c9(1,9);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0., 40.);
      c3.SetXiRange(1.6, 2.0);
      c3.SetTriggerRange(5., 9.);
      c3.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c9);
      
      Correlator c10(1,10);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetCentrality(0., 40.);
      c4.SetXiRange(1.6, 2.0);
      c4.SetTriggerRange(5., 9.);
      c4.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c10);
      
      Correlator c11(1,11);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0., 40.);
      c3.SetXiRange(2.0, 2.4);
      c3.SetTriggerRange(5., 9.);
      c3.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c11);
      
      Correlator c12(1,12);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetCentrality(0., 40.);
      c4.SetXiRange(2.0, 2.4);
      c4.SetTriggerRange(5., 9.);
      c4.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c12);

      Correlator c13(1,13);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0., 40.);
      c3.SetXiRange(0.2, 2.2);
      c3.SetTriggerRange(5., 9.);
      c3.SetAssociatedRange(0.5, 7.0);
      Correlators.push_back(c13);
      
      Correlator c14(1,14);
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
 	    // book(_h["0111"], 1, 1, 1);
      //  //pp 0<xi<0.4
    	 // book(_h["0211"], 2, 1, 1);
      //  //AuAu 0.4<xi<0.8
    	 // book(_h["0311"], 3, 1, 1);
      //  //pp 0.4<xi<0.8
    	 // book(_h["0411"], 4, 1, 1);
      //  //AuAu 0.8<xi<1.2
    	 // book(_h["0511"], 5, 1, 1);
      //  //pp 0.8<xi<1.2
    	 // book(_h["0611"], 6, 1, 1);
      //  //AuAu 1.2<xi<1.6
    	 // book(_h["0711"], 7, 1, 1);
      //  //pp 1.2<xi<1.6
    	 // book(_h["0811"], 8, 1, 1);
      //  //AuAu 1.6<xi<2.0
    	 // book(_h["0911"], 9, 1, 1);
      //  //pp 1.6<xi<2.0
    	 // book(_h["1011"], 10, 1, 1);
      //  //AuAu 2.0<xi<2.4
    	 // book(_h["1111"], 11, 1, 1);
      //  //pp 2.0<xi<2.4
    	 // book(_h["1211"], 12, 1, 1);
      //  //trigger yields for 0-40% Au+Au
    	 // book(_h["1311"], 13, 1, 1);
      //  //trigger yields for 0-40% pp
    	 // book(_h["1411"], 14, 1, 1);
      //  //I_AA for 0-40% Au+Au/p+p
    	 // book(_h["1511"], 15, 1, 1);
      //  //I_AA of 0-40% Au+Au/p+p for |Dphi-pi| < pi/2
    	 // book(_h["1611"], 16, 1, 1);
      //  //I_AA of 0-40% Au+Au/p+p for |Dphi-pi| < pi/3
    	 // book(_h["1711"], 17, 1, 1);
      //  //I_AA of 0-40% Au+Au/p+p for |Dphi-pi| < pi/6
    	 // book(_h["1811"], 18, 1, 1);
      //  //Ratio of |Dphi-pi| < pi/2 I_AA to |Dphi-pi| < pi/6 I_AA for 0-40% Au+Au/p+p
    	 // book(_h["1911"], 19, 1, 1);
      // //Don't know 
      // book(_h["2011"], 20, 1, 1);

      for(Correlator& corr : Correlators)
      {
        book(sow[corr.GetFullIndex()],"sow" + corr.GetFullIndex());
        book(_h["011" + to_string(corr.GetSubIndex())], 1, 1, corr.GetSubIndex());

        string refname = mkAxisCode(1, 1, corr.GetSubIndex());
          const Histo1D& refdata = refData(refname);
          for(auto &bin : refdata.bins())
          {
              book(_DeltaPhixi["011" + to_string(corr.GetSubIndex()) + "xi_" + to_string(bin.xMin()) + "_" + to_string(bin.xMax())], "DeltaPhi_011" + to_string(corr.GetSubIndex()) + "xi_" + to_string(bin.xMin()) + "_" + to_string(bin.xMax()), 24, 0, M_PI);          
          }
          nTriggers[corr.GetFullIndex()] = 0;
      }
      

    } // End of init


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
      
      if (beam.first.pid() == 2212 && beam.second.pid() == 2212) 
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
        //sow[corr.GetIndex()]->fill();
        //nEvents[corr.GetIndex()]++;
        
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
                //nTriggers[corr.GetIndex()]++;
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
       
                for(Correlator& corr : Correlators)
                {
                    if(!corr.CheckConditionsMaxTrigger(SysAndEnergy, centr, pTrig.pt()/GeV, pAssoc.pt()/GeV)) continue;
                    
                } //end of correlators loop 
                
     } // associated hadrons
      } // end of loop over associated particles
    } // trigger hadrons
      } // particle loop

    }


    /// Normalise histograms etc., after the run
    /*void finalize() {
        
        
        for(unsigned int i = 1; i <= Correlators.size() + 5; i++)
         {
          
                 _h["011" + to_string(i)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
            
         }
        
      //normalize correlation histograms by scaling by 1.0/(Ntrig*binwidthphi*binwidtheta) in each bin BUT also be careful when rebinning.  Probably best to FIRST add histograms for correlation functions THEN normalize
      //do background subtraction ala zyam
      //calculate yields

      double norm = sumOfWeights() *2.*M_PI;
      scale(_h["0111"], 1./norm);

    }*/
    
    // Histograms and variables
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Histo1DPtr> _DeltaPhi;
    map<string, Histo1DPtr> _DeltaPhixi;
    map<int, Histo1DPtr> _DeltaPhiSub;
    map<string, CounterPtr> sow;
    bool fillTrigger = true;
    map<string, int> nTriggers;
    //vector<int> nEvents;
    vector<Correlator> Correlators;

   };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PHENIX_2013_I1207323);

}
