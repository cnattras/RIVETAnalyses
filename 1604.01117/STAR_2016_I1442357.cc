
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

using namespace std;
namespace Rivet {
   //------------------------------------------------------------------------------------------------------------------------ 
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



  class STAR_2016_I1442357 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2016_I1442357);

 //Histograms and variables

    map<string, Histo1DPtr> _h;
    map<int, CounterPtr> sow;
    map<int, Histo1DPtr> _DeltaPhi;
    map<int, Histo1DPtr> _DeltaEta;
    map<int, Histo1DPtr> _DeltaPhiSub;
    bool fillTrigger = true;
    vector<int> nTriggers;
    vector<int> nEvents;
    vector<Correlator> Correlators;





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

      const ChargedFinalState cfs(Cuts::abseta < 1.0 && Cuts::pT > 1*GeV);
      declare(cfs, "CFS");
      const ChargedFinalState cfsTrig(Cuts::abseta < 1.0 && Cuts::pT > 2*GeV);
      declare(cfsTrig, "CFSTrig");
      // FinalState of prompt photons and bare muons and electrons in the event
      const PromptFinalState pfsTrig(Cuts::abspid == PID::PHOTON);
      declare(pfsTrig, "PFSTrig");
      // pi0 also


      // Declare centrality projection
      declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");


      //==================================================
      // Create one correlator for each set of Collisions System / Beam Energy / Centrality Interval / Trigger pT interval / Associated pT interval
      // The number used in the constructor is the y-axis index of the histograms (d00-x00-y00) of the paper
      // Ex.: Correlator c1(1); -> is the correlator for histograms _h["0311"], _h["0411"], etc
      // Ex.: Correlator c2(2); -> is the correlator for histograms _h["0312"], _h["0412"], etc
       //==================================================

      Correlator c1(1); 
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 12.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.1, 0.2); // zT
      Correlators.push_back(c1);

      Correlator c2(2); 
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 12.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.2, 0.3); // zT
      Correlators.push_back(c2);

      Correlator c3(3); 
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 12.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.3, 0.4); // zT
      Correlators.push_back(c3);

      Correlator c4(4); 
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 12.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.4, 0.5); // zT
      Correlators.push_back(c4);

      Correlator c5(5); 
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 12.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.5, 0.6); // zT
      Correlators.push_back(c5);

      Correlator c6(6); 
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 12.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.6, 0.7); // zT
      Correlators.push_back(c6);
      
      Correlator c7(7); 
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 12.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.7, 0.8); // zT
      Correlators.push_back(c7);

      Correlator c8(8); 
      c1.SetCollSystemAndEnergy("AuAu200GeV");
      c1.SetCentrality(0., 12.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.8, 0.9); // zT
      Correlators.push_back(c8);
// pp
      Correlator c9(9); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.1, 0.2); // zT
      Correlators.push_back(c9);

      Correlator c10(10); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.2, 0.3); // zT
      Correlators.push_back(c10);

      Correlator c11(11); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.3, 0.4); // zT
      Correlators.push_back(c11);

      Correlator c12(12); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.4, 0.5); // zT
      Correlators.push_back(c12);

      Correlator c13(13); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.5, 0.6); // zT
      Correlators.push_back(c13);

      Correlator c14(14); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.6, 0.7); // zT
      Correlators.push_back(c14);
      
      Correlator c15(15); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.7, 0.8); // zT
      Correlators.push_back(c15);

      Correlator c16(16); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.8, 0.9); // zT
      Correlators.push_back(c16);

      Correlator c17(17); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.10, 0.18); // zT
      Correlators.push_back(c17);

      Correlator c18(18); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.18, 0.20); // zT
      Correlators.push_back(c18);

      Correlator c19(19); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.20, 0.30); // zT
      Correlators.push_back(c19);

      Correlator c20(20); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.30, 0.40); // zT
      Correlators.push_back(c20);
      
      Correlator c21(21); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.40, 0.50); // zT
      Correlators.push_back(c21);

      Correlator c22(22); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.50, 0.59); // zT
      Correlators.push_back(c22);

      Correlator c23(23); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.59, 0.65); // zT
      Correlators.push_back(c23);

      Correlator c24(24); 
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(12., 20.);  //pT_trig
      c1.SetAssociatedRange(0.65, 0.70); // zT
      Correlators.push_back(c24);



      for(unsigned int i = 1; i<= Correlators.size(); i++)
      {
        book(sow[i],"sow" + to_string(i));
        book(_DeltaPhi[i], "DeltaPhi" + to_string(i), 24, 0, M_PI);
        book(_DeltaPhiSub[i], "DeltaPhiSub" + to_string(i), 24, 0, M_PI);
      }

      for(unsigned int i = 1; i<= 21; i++)
      {
          book(_h[to_string(i) + "11"], i, 1, 1); 
      }
      
      nEvents.assign(Correlators.size()+1, 0); 
      nTriggers.assign(Correlators.size()+1, 0); 
      
    }
      // Book histograms
/*
	 // *****FIGURE 1, Delta_phi(rad) // The Azimuthal correlation functions of charged hadrons per trigger,12 < pT^trig < 20 GeV/c
   ////// DO NOT NEED THIS FIGURE ///////////
   //per trigger Yields (rad^-1)
   //Data from Fig. 1a
   book(_h["0111"], 1, 1, 1); //AuAu, 1.2<pT^assoc<3 GeV/c, gamma
   book(_h["0211"], 2, 1, 1); //AuAu, 1.2<pT^assoc<3 GeV/c, pi0
   //Data from Fig. 1b
	 book(_h["0311"], 3, 1, 1); //AuAu, 3<pT^assoc<5 GeV/c, gamma
   book(_h["0411"], 4, 1, 1); //AuAu, 3<pT^assoc<5 GeV/c, pi0
   //Data from Fig. 1c
   book(_h["0511"], 5, 1, 1); //pp, 1.2<pT^assoc<3 GeV/c, gamma
   book(_h["0611"], 6, 1, 1); //pp, 1.2<pT^assoc<3 GeV/c, pi0
   //Data from Fig. 1d
   book(_h["0711"], 7, 1, 1); //pp, 3<pT^assoc<5 GeV/c, gamma
   book(_h["0811"], 8, 1, 1); //pp, 3<pT^assoc<5 GeV/c, pi0
   
   // *****FIGURE 2, zT //pi0-hadron correlations**************************************************
   //D(zT)
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   //Data from Fig. 2a
   book(_h["0911"], 9, 1, 1);  //Away-side, Au+Au
   book(_h["1011"], 10, 1, 1); //Away-side, p+p
   //Data from Fig. 2b
   book(_h["1111"], 11, 1, 1); //Near-side, Au+Au
   book(_h["1211"], 12, 1, 1); //Near-side, p+p
   
   // *****FIGURE 3, zT //Direct Photon-hadron correlations****************************************
   //D(zT)
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   book(_h["1311"], 13, 1, 1); //Away-side, AuAu 
   book(_h["1411"], 14, 1, 1); //Away-side, pp
   
   // *****FIGURE 4, zT_corr **********************************************************************
   //D(zT_Corr)
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   book(_h["1511"], 15, 1, 1); //Pi0-hadron correlations, Away-side, pp
   
   // *****FIGURE 5, zT ***************************************************************************
   //AuAu
   //I_AA
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   book(_h["1611"], 16, 1, 1); //gamma, Away-side, I_AA
   book(_h["1711"], 17, 1, 1); //pi0, Away-side, I_AA
   
   // *****FIGURE 6, zT ***************************************************************************
   //AuAu
   //Ratio
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   book(_h["1811"], 18, 1, 1); //gamma
   book(_h["1911"], 19, 1, 1); //pi0
   
   // *****FIGURE 7 *******************************************************************************
   //I_AA
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   book(_h["2011"], 20, 1, 1); //gamma, Away-side, I_AA (pT^trig)
   book(_h["2111"], 21, 1, 1); //gamma, Away-side, I_AA (pT^assoc)
   
   */
  //}

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const PromptFinalState& pfsTrig = apply<PromptFinalState>(event, "PFSTrig");
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
      if(cmsEnergy.compare("Empty") == 0) return;
      
      string SysAndEnergy = CollSystem + cmsEnergy;


      /// DO NOT MESS WITH THIS
      // Prepare centrality projection and value
      const CentralityProjection& centrProj = apply<CentralityProjection>(event, "CMULT");
      double centr = centrProj();

      // Veto event for too large centralities since those are not used in the analysis at all
      if ((centr < 0.) || (centr > 12.)){
        vetoEvent;
      }
    

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
          
                for(Correlator& corr : Correlators)
                {
                    if(!corr.CheckConditionsMaxTrigger(SysAndEnergy, centr, pTrig.pt()/GeV, pAssoc.pt()/GeV)) continue;
                    
                    //double etaCorrection = EtaEffCorrection(abs(dEta), *_DeltaEta[corr.GetIndex()]);
                
                        _h[to_string(corr.GetIndex()) + "11"]->fill(-abs(dPhi),0.5); 
                        //_DeltaPhi[corr.GetIndex()]->fill(abs(dPhi), 0.5/etaCorrection);
                        //_DeltaPhiSub[corr.GetIndex()]->fill(abs(dPhi), 0.5/etaCorrection);
                    
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
                _h[to_string(i) + "11"]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                
                _DeltaPhi[i]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                _DeltaPhiSub[i]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                
                vector<int> n{2,3};
               // SubtractBackground(*_DeltaPhi[i], *_h["061" + to_string(i)], n, 0.63, 2.51);
               // SubtractBackground(*_DeltaPhi[i], *_DeltaPhiSub[i], n, 0.63, 2.51);
            }    
        }
        
      //normalize correlation histograms by scaling by 1.0/(Ntrig*binwidthphi*binwidtheta) in each bin BUT also be careful when rebinning.  Probably best to FIRST add histograms for correlation functions THEN normalize
      //do background subtraction ala zyam
      //calculate yields

    } //end void finalize

};
  DECLARE_RIVET_PLUGIN(STAR_2016_I1442357);
}

