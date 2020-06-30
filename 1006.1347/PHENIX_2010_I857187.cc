
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include "RHICCentrality.hh"
#define _USE_MATH_DEFINES
static const int numTrigPtBins = 4;
static const float pTTrigBins[] = {5.0,7.0,9.0,12.0,15.0};
static const int numXeBins = 14;
static const float XeBins[] = {0.1,0.16,0.2,0.26,0.3,0.36,0.4,0.46,0.5,0.56,0.71,0.85,0.95,1.2,1.5};
static const int numAssocPtBins = 4;
static const float pTAssocBins[] = {5.0,7.0,9.0,12.0,15.0};
static const int numPoutBins = 8;
static const float PoutBins[] = {0.2,0.3,0.8,1.5,2.0,3.0,4.5,6.5,9.0};

using namespace std;
namespace Rivet {
    
  class Correlator {
      
    private:

      int _index;
      string _collSystemAndEnergy;
      pair<double,double> _centrality;
      pair<double,double> _triggerRange;
      pair<double,double> _associatedRange;

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
  };

  class PHENIX_2010_I857187 : public Analysis {
  
    public:

     /// Constructor
      DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2010_I857187);


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

      /// Book histograms and initialise projections before the run
      void init() {

        // Initialise and register projections

        // the basic final-state projection: all final-state particles within the given eta acceptance

        const ChargedFinalState cfs(Cuts::abseta < 0.35);
        declare(cfs, "CFS");
        const ChargedFinalState cfsTrig(Cuts::abseta < 0.35);
        declare(cfsTrig, "CFSTrig");

        // the basic final-state projection: all final-state photon within the given eta acceptance


        // Declare centrality projection
        declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      //==================================================
      // Create one correlator for each set of Collisions System / Beam Energy / Centrality Interval / Trigger pT interval / Associated pT interval
      // The number used in the constructor is the y-axis index of the histograms (d00-x00-y00) of the paper
      // Ex.: Correlator c1(1); -> is the correlator for histograms _h["0311"], _h["0411"], etc
      // Ex.: Correlator c2(2); -> is the correlator for histograms _h["0312"], _h["0412"], etc
      //==================================================
      
      //Real correlators
      /*Correlator c1(1);
      c1.SetCollSystemAndEnergy("pp200Gev");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(5., 7.);
      c1.SetAssociatedRange(5., 7.);
      Correlators.push_back(c1);
      
      Correlator c2(2);
      c2.SetCollSystemAndEnergy("pp200Gev");
      c2.SetCentrality(0., 100.);
      c2.SetTriggerRange(7., 9.);
      c2.SetAssociatedRange(7., 9.);
      Correlators.push_back(c2);
      
      Correlator c3(3);
      c3.SetCollSystemAndEnergy("pp200Gev");
      c3.SetCentrality(0., 100.);
      c3.SetTriggerRange(9., 12.);
      c3.SetAssociatedRange(9., 12.);
      Correlators.push_back(c3);
      
      Correlator c4(4);
      c4.SetCollSystemAndEnergy("pp200Gev");
      c4.SetCentrality(0., 100.);
      c4.SetTriggerRange(12., 15.);
      c4.SetAssociatedRange(12., 15.);
      Correlators.push_back(c4);
      
      Correlator c5(5);
      c5.SetCollSystemAndEnergy("pp200Gev");
      c5.SetCentrality(0., 100.);
      c5.SetTriggerRange(5., 7.);
      c5.SetAssociatedRange(5., 7.);
      Correlators.push_back(c5);
      
      Correlator c6(6);
      c6.SetCollSystemAndEnergy("pp200Gev");
      c6.SetCentrality(0., 100.);
      c6.SetTriggerRange(7., 9.);
      c6.SetAssociatedRange(7., 9.);
      Correlators.push_back(c6);

      Correlator c7(7);
      c7.SetCollSystemAndEnergy("pp200Gev");
      c7.SetCentrality(0., 100.);
      c7.SetTriggerRange(9., 12.);
      c7.SetAssociatedRange(9., 12.);
      Correlators.push_back(c7);

      Correlator c8(8);
      c8.SetCollSystemAndEnergy("pp200Gev");
      c8.SetCentrality(0., 100.);
      c8.SetTriggerRange(12., 15.);
      c8.SetAssociatedRange(12., 15.);
      Correlators.push_back(c8);*/
      
      //Debug correlators
      Correlator c1(1);
      c1.SetCollSystemAndEnergy("pp200GeV");
      c1.SetCentrality(0., 100.);
      c1.SetTriggerRange(0., 7.);
      c1.SetAssociatedRange(0., 7.);
      Correlators.push_back(c1);
      
      Correlator c2(2);
      c2.SetCollSystemAndEnergy("pp200GeV");
      c2.SetCentrality(0., 100.);
      c2.SetTriggerRange(0., 9.);
      c2.SetAssociatedRange(0., 9.);
      Correlators.push_back(c2);
      
      Correlator c3(3);
      c3.SetCollSystemAndEnergy("pp200GeV");
      c3.SetCentrality(0., 100.);
      c3.SetTriggerRange(0., 12.);
      c3.SetAssociatedRange(0., 12.);
      Correlators.push_back(c3);
      
      Correlator c4(4);
      c4.SetCollSystemAndEnergy("pp200GeV");
      c4.SetCentrality(0., 100.);
      c4.SetTriggerRange(0., 15.);
      c4.SetAssociatedRange(0., 15.);
      Correlators.push_back(c4);
      
      Correlator c5(5);
      c5.SetCollSystemAndEnergy("pp200GeV");
      c5.SetCentrality(0., 100.);
      c5.SetTriggerRange(0., 7.);
      c5.SetAssociatedRange(0., 7.);
      Correlators.push_back(c5);
      
      Correlator c6(6);
      c6.SetCollSystemAndEnergy("pp200GeV");
      c6.SetCentrality(0., 100.);
      c6.SetTriggerRange(0., 9.);
      c6.SetAssociatedRange(0., 9.);
      Correlators.push_back(c6);

      Correlator c7(7);
      c7.SetCollSystemAndEnergy("pp200GeV");
      c7.SetCentrality(0., 100.);
      c7.SetTriggerRange(0., 12.);
      c7.SetAssociatedRange(0., 12.);
      Correlators.push_back(c7);

      Correlator c8(8);
      c8.SetCollSystemAndEnergy("pp200GeV");
      c8.SetCentrality(0., 100.);
      c8.SetTriggerRange(0., 15.);
      c8.SetAssociatedRange(0., 15.);
      Correlators.push_back(c8);

 /* Old booking code! 
      // Book histograms
      //Pi0 - Hadron 5 < pT,trigger < 7 GeV/c
  book(_h["0111"], 1, 1, 1);
   // Pi0 - Hadron 7 < pT,trigger < 9 GeV/c
	 book(_h["0211"], 2, 1, 1);
   //Pi0 - Hadron 9 < pT,trigger < 12 GeV/c
	 book(_h["0311"], 3, 1, 1);
   //Pi0 - Hadron 12 < pT,trigger < 15 GeV/c
	 book(_h["0411"], 4, 1, 1);
   //Isolated Direct Photon - Hadron 5 < pT,trigger < 7 GeV/c
	 book(_h["0511"], 5, 1, 1);
   //Isolated Direct Photon - Hadron 7 < pT,trigger < 9 GeV/c
	 book(_h["0611"], 6, 1, 1);
   //Isolated Direct Photon - Hadron 9 < pT,trigger < 12 GeV/c
	 book(_h["0711"], 7, 1, 1);
   //Isolated Direct Photon - Hadron 12 < pT,trigger < 15 GeV/c
	 book(_h["0811"], 8, 1, 1);
   //Pi0 - Hadron 5 < pT,t < 7
	 book(_h["0911"], 9, 1, 1);
   //Pi0 - Hadron 7 < pT,t < 9
	 book(_h["1011"], 10, 1, 1);
   //Pi0 - Hadron 9 < pT,t < 12
	 book(_h["1111"], 11, 1, 1);
   //Pi0 - Hadron 12 < pT,t < 15
	 book(_h["1211"], 12, 1, 1);
   //Isolated Direct Photon - Hadron 5 < pT,t < 7
	 book(_h["1311"], 13, 1, 1);
   //Isolated Direct Photon - Hadron 7 < pT,t < 9
	 book(_h["1411"], 14, 1, 1);
   //Isolated Direct Photon - Hadron 9 < pT,t < 12
	 book(_h["1511"], 15, 1, 1);
   //Isolated Direct Photon - Hadron 12 < pT,t < 15
	 book(_h["1611"], 16, 1, 1);
    }
  */

      for(unsigned int i = 1; i<= 8; i++)
      {
        book(sow[i],"sow" + to_string(i));
        switch (i) {
          case 1:
            for(int ii=1; ii<= 2; ii++){
              book(_h["0" + to_string(i) + "1" + to_string(ii)], i, 1, ii);
            }
            break;
          case 2:
            for(int ii=1; ii<= 2; ii++){
              book(_h["0" + to_string(i) + "1" + to_string(ii)], i, 1, ii);
            }
            break;
          case 3:
            for(int ii=1; ii<= 2; ii++){
              book(_h["0" + to_string(i) + "1" + to_string(ii)], i, 1, ii);
            }
            break;
          case 4:
            for(int ii=1; ii<= 2; ii++){
              book(_h["0" + to_string(i) + "1" + to_string(ii)], i, 1, ii);
            }
            break;
          case 5:
            for(int ii=1; ii<= 2; ii++){
              book(_h["0" + to_string(i) + "1" + to_string(ii)], i, 1, ii);
            }
            break;
          case 6:
            for(int ii=1; ii<= 2; ii++){
              book(_h["0" + to_string(i) + "1" + to_string(ii)], i, 1, ii);
            }
            break;
          case 7:
            for(int ii=1; ii<= 2; ii++){
              book(_h["0" + to_string(i) + "1" + to_string(ii)], i, 1, ii);
            }
            break;
          case 8:
            for(int ii=1; ii<= 2; ii++){
              book(_h["0" + to_string(i) + "1" + to_string(ii)], i, 1, ii);
            }
            break;
          //book(_h[to_string(i) + "11"], i, 1, 1);
        }
          book(_DeltaPhi[i], "DeltaPhi" + to_string(i), 24, 0, M_PI);
          book(_DeltaPhiSub[i], "DeltaPhiSub" + to_string(i), 24, 0, M_PI);
      }
      
      nEvents.assign(Correlators.size()+1, 0); 
      nTriggers.assign(Correlators.size()+1, 0); 
      
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
          CollSystem = "pp";
          nNucleons = 1.;
      //if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000010020) CollSystem = "dAu";
      //if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000290630) CollSystem = "dAu";
      
      string cmsEnergy = "Empty";
      if (fuzzyEquals(sqrtS()/GeV, 200*nNucleons, 1E-3)) cmsEnergy = "200GeV";
      
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
        //if(!corr.CheckCentrality(centr)) continue;
        
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
        //cout << "trigger loop" << endl;
        if(pTrig.pt()/GeV < triggerptMin || pTrig.pt()/GeV > triggerptMax) continue;
        //cout << "trigger pt: " << pTrig.pT()/GeV << endl;
          //Check if is secondary
        if(isSecondary(pTrig)) continue;
          
        if( abs(pTrig.pid())==211 || abs(pTrig.pid())==2212 || abs(pTrig.pid())==321){

          for(Correlator& corr : Correlators)
          {
            //if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
            nTriggers[corr.GetIndex()]++;
          }

          // Hadron loop
          for(const Particle& pAssoc : cfs.particles()) {
                //cout << "Assoc loop" << endl;
            if(pAssoc.pt()/GeV < associatedptMin || pAssoc.pt()/GeV > pTrig.pt()/GeV) continue;
                
            //Check if Trigger and Associated are the same particle
            if(isSameParticle(pTrig,pAssoc)) continue;
                
            //Check if is secondary
            if(isSecondary(pAssoc)) continue;
                
            if(abs(pAssoc.pid())==211 || abs(pAssoc.pid())==2212 || abs(pAssoc.pid())==321) {
            //int mybin = GetTrigBin(pTrig.pt());
            //int mybina = GetAssocBin(pAssoc.pt());

            //https://rivet.hepforge.org/code/dev/structRivet_1_1DeltaPhiInRange.html
            double dPhi = deltaPhi(pTrig, pAssoc, true);//this does NOT rotate the delta phi to be in a given range
            double dEta = deltaEta(pTrig, pAssoc);
          
            //cout << "dPhi: " << dPhi << " dEta: " << dEta << endl;
            for(Correlator& corr : Correlators)
            {
              //if(!corr.CheckConditionsMaxTrigger(SysAndEnergy, centr, pTrig.pt()/GeV, pAssoc.pt()/GeV)) continue;
                    
              if(abs(dPhi) < 0.78)
              {
                switch(corr.GetIndex()) {
                  
                  case 1:
                    for(int ii=1; ii<= 2; ii++){
                        //cout << "Filling histogram" << endl;
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dEta), 0.5);
                    }
                    break;
                  case 2:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dEta), 0.5);
                    }
                    break;
                  case 3:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dEta), 0.5);
                    }
                    break;
                  case 4:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dEta), 0.5);
                    }
                    break;
                  case 5:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dEta), 0.5);
                    }
                    break;
                  case 6:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dEta), 0.5);
                    }
                    break;
                  case 7:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dEta), 0.5);
                    }
                    break;
                  case 8:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dEta), 0.5);
                    }
                    break;
                }
              }
                switch(corr.GetIndex()) {
                  
                  case 1:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dPhi), 0.5);
                    }
                    break;
                  case 2:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dPhi), 0.5);
                    }
                    break;
                  case 3:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dPhi), 0.5);
                    }
                    break;
                  case 4:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dPhi), 0.5);
                    }
                    break;
                  case 5:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dPhi), 0.5);
                    }
                    break;
                  case 6:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dPhi), 0.5);
                    }
                    break;
                  case 7:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dPhi), 0.5);
                    }
                    break;
                  case 8:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(corr.GetIndex()) + "1" + to_string(ii)]->fill(abs(dPhi), 0.5);
                    }
                    break;
                }
              
            
                    
            } //end of correlators loop 
                
          } // associated hadrons
        } // end of loop over associated particles
      } // trigger hadrons

    } // particle loop
  }
    /// Normalise histograms etc., after the run
    void finalize() {
        
        for(unsigned int i = 1; i <= 8; i++) /*I set the the range to i <= 10 at first, but then I ran the
                                              code and there were issues. So I made it leass than 10, i.e. 9,
                                              but it gives me nan for each entry in the Rivet.yoda file.*/
        {
           switch(i) {
                  
                  case 1:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(i) + "1" + to_string(ii)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                    }
                    break;
                  case 2:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(i) + "1" + to_string(ii)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                    }
                    break;
                  case 3:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(i) + "1" + to_string(ii)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                    }
                    break;
                  case 4:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(i) + "1" + to_string(ii)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                    }
                    break;
                  case 5:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(i) + "1" + to_string(ii)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                    }
                    break;
                  case 6:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(i) + "1" + to_string(ii)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                    }
                    break;
                  case 7:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(i) + "1" + to_string(ii)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                    }
                    break;
                  case 8:
                    for(int ii=1; ii<= 2; ii++){
                      _h["0" + to_string(i) + "1" + to_string(ii)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                    }
                    break;
              }
          //This is the old code. I switched it with the code above, but running an analysis with the code above game me some issues.
            /*if(nTriggers[i] > 0)
            {
              if (i <= 9) {
                //_h["0" + to_string(i) + "11"]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
              }
              else {
                //_h[to_string(i) + "11"]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
              }*/
              _DeltaPhi[i]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
              _DeltaPhiSub[i]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));

              //vector<int> n{2,3};
              //SubtractBackground(*_DeltaPhi[i], *_DeltaPhiSub[i], n, 0.63, 2.51);
            
            
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
    map<int, Histo1DPtr> _DeltaPhi;
    map<int, Histo1DPtr> _DeltaPhiSub;
    bool fillTrigger = true;
    vector<int> nTriggers;
    vector<int> nEvents;
    vector<Correlator> Correlators;

    /* 
    map<string, Histo1DPtr> _h;
    map<int, CounterPtr> sow;
    bool fillTrigger = true;
    vector<int> nTriggers;
    vector<int> nEvents;
    vector<Correlator> Correlators;
    */

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PHENIX_2010_I857187);
  }
