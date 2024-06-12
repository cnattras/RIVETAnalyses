
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
#include "../Centralities/RHICCentrality.hh"

#define _USE_MATH_DEFINES

static const int numTrigPtBins = 3;
static const float pTTrigBins[] = {2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0};
static const int numAssocPtBins = 10;
static const float pTAssocBins[] = {1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0};
static const int numCentBins = 8;
static const int centBins[] = {0,5,10,20, 30,40, 50,60};

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










  class PHENIX_2009_I815824 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2009_I815824);


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

      // Initialise and register projections

      // the basic final-state projection: all final-state particles within the given eta acceptance

      const ChargedFinalState cfs(Cuts::abseta < 5 && Cuts::pT > 100*MeV);
      declare(cfs, "CFS");
	  
	  const ChargedFinalState cfsTrig(Cuts::abseta < 1.0 && Cuts::pT > 2*GeV);
      declare(cfsTrig, "CFSTrig");

      // Declare centrality projection
	  
	  //Later Fix to Use STAR(/PHENIX)?
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
	  
	  //==================================================
      // Create one correlator for each set of Collisions System / Beam Energy / Centrality Interval / Trigger pT interval / Associated pT interval
      // The number used in the constructor is the y-axis index of the histograms (d00-x00-y00) of the paper
     
	 
	 // Need 24 total from fig 4
	 // +6 for fig 6
	 // redo these....................
	 
	 //AuAu200GeV
	 
		//5-7
	 Correlator c1(1);
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0,20.);
      c1.SetTriggerRange(5., 7.);
      c1.SetAssociatedRange(1., 2.);
      Correlators.push_back(c1);
	  
	  Correlator c2(2);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetCentrality(0,20.);
      c2.SetTriggerRange(5., 7.);
      c2.SetAssociatedRange(2., 3.);
      Correlators.push_back(c2);
	  
	  Correlator c3(3);
      c3.SetCollSystemAndEnergy("AuAu200GeV");
      c3.SetCentrality(0,20.);
      c3.SetTriggerRange(5., 7.);
      c3.SetAssociatedRange(3., 5.);
      Correlators.push_back(c3);
	  
		//7-9
	  Correlator c4(4);
      c4.SetCollSystemAndEnergy("AuAu200GeV");
      c4.SetCentrality(0,20.);
      c4.SetTriggerRange(7., 9.);
      c4.SetAssociatedRange(1., 2.);
      Correlators.push_back(c4);
	  
	  Correlator c5(5);
      c5.SetCollSystemAndEnergy("AuAu200GeV");
      c5.SetCentrality(0,20.);
      c5.SetTriggerRange(7., 9.);
      c5.SetAssociatedRange(2., 3.);
      Correlators.push_back(c5);
	  
	  Correlator c6(6);
      c6.SetCollSystemAndEnergy("AuAu200GeV");
      c6.SetCentrality(0,20.);
      c6.SetTriggerRange(7., 9.);
      c6.SetAssociatedRange(3., 5.);
      Correlators.push_back(c6);
	  
	  
	  //9-12
	  Correlator c7(7);
      c7.SetCollSystemAndEnergy("AuAu200GeV");
      c7.SetCentrality(0,20.);
      c7.SetTriggerRange(9., 12.);
      c7.SetAssociatedRange(1., 2.);
      Correlators.push_back(c7);
	  
	  Correlator c8(8);
      c8.SetCollSystemAndEnergy("AuAu200GeV");
      c8.SetCentrality(0,20.);
      c8.SetTriggerRange(9., 12.);
      c8.SetAssociatedRange(2., 3.);
      Correlators.push_back(c8);
	  
	  Correlator c9(9);
      c9.SetCollSystemAndEnergy("AuAu200GeV");
      c9.SetCentrality(0,20.);
      c9.SetTriggerRange(9., 12.);
      c9.SetAssociatedRange(3., 5.);
      Correlators.push_back(c9);
	 
		//12-15
	 Correlator c10(10);
      c10.SetCollSystemAndEnergy("AuAu200GeV");
      c10.SetCentrality(0,20.);
      c10.SetTriggerRange(12., 15.);
      c10.SetAssociatedRange(1., 2.);
      Correlators.push_back(c10);
	  
	  Correlator c11(11);
      c11.SetCollSystemAndEnergy("AuAu200GeV");
      c11.SetCentrality(0,20.);
      c11.SetTriggerRange(12., 15.);
      c11.SetAssociatedRange(2., 3.);
      Correlators.push_back(c11);
	  
	  Correlator c12(12);
      c12.SetCollSystemAndEnergy("AuAu200GeV");
      c12.SetCentrality(0,20.);
      c12.SetTriggerRange(12., 15.);
      c12.SetAssociatedRange(3., 5.);
      Correlators.push_back(c12);
	 
	 //pp200GeV
	 		//5-7
	 Correlator c13(13);
      c13.SetCollSystemAndEnergy("pp200GeV");
      c13.SetCentrality(0,20.);
      c13.SetTriggerRange(5., 7.);
      c13.SetAssociatedRange(1., 2.);
      Correlators.push_back(c13);
	  
	  Correlator c14(14);
      c14.SetCollSystemAndEnergy("pp200GeV");
      c14.SetCentrality(0,20.);
      c14.SetTriggerRange(5., 7.);
      c14.SetAssociatedRange(2., 3.);
      Correlators.push_back(c14);
	  
	  Correlator c15(15);
      c15.SetCollSystemAndEnergy("pp200GeV");
      c15.SetCentrality(0,20.);
      c15.SetTriggerRange(5., 7.);
      c15.SetAssociatedRange(3., 5.);
      Correlators.push_back(c15);
	  
		//7-9
	  Correlator c16(16);
      c16.SetCollSystemAndEnergy("pp200GeV");
      c16.SetCentrality(0,20.);
      c16.SetTriggerRange(7., 9.);
      c16.SetAssociatedRange(1., 2.);
      Correlators.push_back(c16);
	  
	  Correlator c17(17);
      c17.SetCollSystemAndEnergy("pp200GeV");
      c17.SetCentrality(0,20.);
      c17.SetTriggerRange(7., 9.);
      c17.SetAssociatedRange(2., 3.);
      Correlators.push_back(c17);
	  
	  Correlator c18(18);
      c18.SetCollSystemAndEnergy("pp200GeV");
      c18.SetCentrality(0,20.);
      c18.SetTriggerRange(7., 9.);
      c18.SetAssociatedRange(3., 5.);
      Correlators.push_back(c18);
	  
	  
	  //9-12
	  Correlator c19(19);
      c19.SetCollSystemAndEnergy("pp200GeV");
      c19.SetCentrality(0,20.);
      c19.SetTriggerRange(9., 12.);
      c19.SetAssociatedRange(1., 2.);
      Correlators.push_back(c19);
	  
	  Correlator c20(20);
      c20.SetCollSystemAndEnergy("pp200GeV");
      c20.SetCentrality(0,20.);
      c20.SetTriggerRange(9., 12.);
      c20.SetAssociatedRange(2., 3.);
      Correlators.push_back(c20);
	  
	  Correlator c21(21);
      c21.SetCollSystemAndEnergy("pp200GeV");
      c21.SetCentrality(0,20.);
      c21.SetTriggerRange(9., 12.);
      c21.SetAssociatedRange(3., 5.);
      Correlators.push_back(c21);
	 
		//12-15
	 Correlator c22(22);
      c22.SetCollSystemAndEnergy("pp200GeV");
      c22.SetCentrality(0,20.);
      c22.SetTriggerRange(12., 15.);
      c22.SetAssociatedRange(1., 2.);
      Correlators.push_back(c22);
	  
	  Correlator c23(23);
      c23.SetCollSystemAndEnergy("pp200GeV");
      c23.SetCentrality(0,20.);
      c23.SetTriggerRange(12., 15.);
      c23.SetAssociatedRange(2., 3.);
      Correlators.push_back(c23);
	  
	  Correlator c24(24);
      c24.SetCollSystemAndEnergy("pp200GeV");
      c24.SetCentrality(0,20.);
      c24.SetTriggerRange(12., 15.);
      c24.SetAssociatedRange(3., 5.);
      Correlators.push_back(c24);
	  
	  	 
	  
	  ///////////////fig 6
	  
	  // 0-20%
	 Correlator c25(25);
      c25.SetCollSystemAndEnergy("y-h200GeV");
      c25.SetCentrality(0,20.);
      c25.SetTriggerRange(5., 15.);
      c25.SetAssociatedRange(3., 5.);
      Correlators.push_back(c25);
	  
	  Correlator c26(26);
      c26.SetCollSystemAndEnergy("hh200GeV");
      c26.SetCentrality(0,20.);
      c26.SetTriggerRange(5., 10.);
      c26.SetAssociatedRange(3., 4.);
      Correlators.push_back(c26);
	  
	  
	  //20-40%
	  Correlator c27(27);
      c27.SetCollSystemAndEnergy("y-h200GeV");
      c27.SetCentrality(20.,40.);
      c27.SetTriggerRange(5., 15.);
      c27.SetAssociatedRange(3., 5.);
      Correlators.push_back(c27); 
	  
	  
	 
	  Correlator c28(28);
      c28.SetCollSystemAndEnergy("hh200GeV");
      c28.SetCentrality(20.,40.);
      c28.SetTriggerRange(5., 10.);
      c28.SetAssociatedRange(3., 4.);
      Correlators.push_back(c28);
	  
	  
	  //40-60%
	  Correlator c29(29);
      c29.SetCollSystemAndEnergy("y-h200GeV");
      c29.SetCentrality(40.,60.);
      c29.SetTriggerRange(5., 15.);
      c29.SetAssociatedRange(3., 5.);
      Correlators.push_back(c29);
	  
	  Correlator c30(30);
      c30.SetCollSystemAndEnergy("hh200GeV");
      c30.SetCentrality(40.,60.);
      c30.SetTriggerRange(5., 10.);
      c30.SetAssociatedRange(3., 4.);
      Correlators.push_back(c30);
	  
      
	  
		
	  
	  

      // Book histograms
		//Head Region Yields 0-20% Centrality 5-7
	 book(_h["0111"], 1, 1, 1);
	 book(_h["0112"], 1, 1, 2);
		//Head Region Yields 0-20% Centrality 7-9
	 book(_h["0211"], 2, 1, 1);
	 book(_h["0212"], 2, 1, 2);
		//Head Region Yields 0-20% Centrality 9-12
	 book(_h["0311"], 3, 1, 1);
	 book(_h["0312"], 3, 1, 2);
		//Head Region Yields 0-20% Centrality 12-15
	 book(_h["0411"], 4, 1, 1);
	 book(_h["0412"], 4, 1, 2);
		//Head Region Yields 20-40% Centrality 5-7
	 book(_h["0511"], 5, 1, 1);
	 book(_h["0512"], 5, 1, 2);
		//Head Region Yields 20-40% Centrality 7-9
	 book(_h["0611"], 6, 1, 1);
	 book(_h["0612"], 6, 1, 2);
		//Head Region Yields 20-40% Centrality 9-12
	 book(_h["0711"], 7, 1, 1);
	 book(_h["0712"], 7, 1, 2);
		//Head Region Yields 20-40% Centrality 12-15
	 book(_h["0811"], 8, 1, 1);
	 book(_h["0812"], 8, 1, 2);
		//Head Region Yields 40-60% Centrality 5-7
	 book(_h["0911"], 9, 1, 1);
	 book(_h["0912"], 9, 1, 2);
		//Head Region Yields 40-60% Centrality 7-9
	 book(_h["1011"], 10, 1, 1);
	 book(_h["1012"], 10, 1, 2);
		//Head Region Yields 40-60% Centrality 9-12
	 book(_h["1111"], 11, 1, 1);
	 book(_h["1112"], 11, 1, 2);
		//Head Region Yields 40-60% Centrality 12-15
	 book(_h["1211"], 12, 1, 1);
	 book(_h["1212"], 12, 1, 2);
		//I_{AA} Ratio fig5a
	 book(_h["1311"], 13, 1, 1);
		//I_{AA} Ratio fig5b
	 book(_h["1411"], 14, 1, 1);
		//I_{AA} Ratio fig5c
	 book(_h["1511"], 15, 1, 1);
		//I_{AA} Ratio fig5d
	 book(_h["1611"], 16, 1, 1);
		//I_{AA}(Z_T) Ratio fig8a
	 book(_h["1711"], 17, 1, 1);
		//I_{AA}(Z_T) Ratio fig8b
	 book(_h["1811"], 18, 1, 1);
		//I_{AA}(Z_T) Ratio fig8c
	 book(_h["1911"], 19, 1, 1);
		//I_{AA}(Z_T) Ratio fig8d
	 book(_h["2011"], 20, 1, 1);
	 
	 
	 
	 
	  for(unsigned int i = 1; i<= Correlators.size(); i++)
      {
		  book(sow[i],"sow" + to_string(i));
		  if (i < 10 ){
			//book(sow[i],"sow" + to_string(i));
			//book(_h["0" + to_string(i) + "11"], i, 1, 1);
			//book(_h["0" + to_string(i) + "11"], i, 1, 2);
		  }
		  else {
			//book(sow[i],"sow" + to_string(i));
			
			//book(_h[to_string(i) + "11"], i, 1, 1);
		  }
      }
      
      nEvents.assign(Correlators.size()+1, 0); 
      nTriggers.assign(Correlators.size()+1, 0); 
      
    }
	 
	 
	 
    


    /// Perform the per-event analysis
    void analyze(const Event& event) {
		

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
	  
	  const ChargedFinalState& cfsTrig = apply<ChargedFinalState>(event, "CFSTrig");

      // Prepare centrality projection and value
      //const CentralityProjection& centrProj = apply<CentralityProjection>(event, "V0M");
      //double centr = centrProj();
	  
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
                    
                    
					double etaCorrection = EtaEffCorrection(abs(dEta), *_DeltaEta[corr.GetIndex()]);
					
					if(abs(dPhi) < 0.78)
                    {
                        _h["031" + to_string(corr.GetIndex())]->fill(-abs(dEta), 0.5);
                    }
                    
                    if(abs(dEta) < 0.78)
                    {
                        _h["041" + to_string(corr.GetIndex())]->fill(-abs(dPhi), 0.5);
						_h["061" + to_string(corr.GetIndex())]->fill(-abs(dPhi), 0.5/etaCorrection);
                        _DeltaPhi[corr.GetIndex()]->fill(abs(dPhi), 0.5/etaCorrection);
                        _DeltaPhiSub[corr.GetIndex()]->fill(abs(dPhi), 0.5/etaCorrection);
                    }
                    
                } //end of correlators loop 
                
			  } // associated hadrons
		    } // end of loop over associated particles
		  } // trigger hadrons
      } // particle loop

    }


      
    /// Normalise histograms etc., after the run
    void finalize() {

      //double norm = sumOfWeights() *2.*M_PI;
      //scale(_h["0111"], 1./norm);


		for(unsigned int i = 1; i <= Correlators.size(); i++)
        {
            if(nTriggers[i] > 0)
            {
                _h["031" + to_string(i)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                _h["041" + to_string(i)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
				_h["061" + to_string(i)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                
                _DeltaPhi[i]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                _DeltaPhiSub[i]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                
                vector<int> n{2,3};
                SubtractBackground(*_DeltaPhi[i], *_h["061" + to_string(i)], n, 0.63, 2.51);
                SubtractBackground(*_DeltaPhi[i], *_DeltaPhiSub[i], n, 0.63, 2.51);
            }
				
				
				
            
            
        }


    }

    //@}


    /// name Histograms and Variables
    //{
    map<string, Histo1DPtr> _h;
   // map<string, Profile1DPtr> _p;
   // map<string, CounterPtr> _c;
    //}
	map<int, CounterPtr> sow;
	map<int, Histo1DPtr> _DeltaPhi;
    map<int, Histo1DPtr> _DeltaEta;
    map<int, Histo1DPtr> _DeltaPhiSub;
    bool fillTrigger = true;
    vector<int> nTriggers;
    vector<int> nEvents;
    vector<Correlator> Correlators;
    string beamOpt;
    enum CollisionSystem { pp, AuAu, CuCu};
    CollisionSystem collSys;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PHENIX_2009_I815824);


}