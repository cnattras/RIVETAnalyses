// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

#include <fstream>
#include <iostream>
#include <math.h>
#include <vector> 
#include "Rivet/Projections/RHICCentrality.hh"
#define _USE_MATH_DEFINES

using namespace std;

namespace Rivet {
    
  class Correlator {
  
    private:
      
      int _index;
      string _collSystemAndEnergy;
      pair<double, double> _centrality;
      pair<double, double> _triggerRange;
      pair<double, double> _associatedRange;
      pair<double, double> _azimuthalRange;
  
    public:
    
      /// Constructor
      Correlator(int index) {
        _index = index;
      }

      void SetCollSystemAndEnergy(string s){ _collSystemAndEnergy = s; }
      void SetCentrality(double cmin, double cmax){ _centrality = make_pair(cmin, cmax); }
      void SetTriggerRange(double tmin, double tmax){ _triggerRange = make_pair(tmin, tmax); }
      void SetAssociatedRange(double amin, double amax){ _associatedRange = make_pair(amin, amax); }
      
      // Sets azimuthal angle with respect to the reaction plane. 
      void SetAzimuthalRange(double pmin, double pmax) { _azimuthalRange = make_pair(pmin, pmax); }
      
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
      
      double GetAzimuthalRangeMin() { return _azimuthalRange.first; }
      double GetAzimuthalRangeMax() { return _azimuthalRange.second; }
    
      int GetIndex(){ return _index; }
    
      bool CheckCollSystemAndEnergy(string s){ return _collSystemAndEnergy.compare(s) == 0 ? true : false; }
      bool CheckCentrality(double cent){ return (cent>_centrality.first && cent<_centrality.second) ? true : false; }
      bool CheckTriggerRange(double tpt){ return (tpt>_triggerRange.first && tpt<_triggerRange.second) ? true : false; }
      bool CheckAssociatedRange(double apt){ return (apt>_associatedRange.first && apt<_associatedRange.second) ? true : false; }
      bool CheckAssociatedRangeMaxTrigger(double apt, double tpt){ return (apt>_associatedRange.first && apt<tpt) ? true : false; }
      
      bool CheckAzimuthalRange(double azi) { return (azi > _azimuthalRange.first && azi < _azimuthalRange.second) ? true : false; }
      
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
  
  }; // End of Correlator class

  /// @brief Add a short analysis description here
  class PHENIX_2011_I872172 : public Analysis {

    private:
      
      //Histograms and variables
      map<string, Histo1DPtr> _h;
      map<int, CounterPtr> sow;
      map<int, Histo1DPtr> _DeltaPhi;
      map<int, Histo1DPtr> _DeltaEta;
      map<int, Histo1DPtr> _DeltaPhiSub;
      map<int, int> nTriggers;
      map<int, int> nEvents;
      bool fillTrigger = true;
      vector<Correlator> Correlators;


    public:

      /// Constructor
      DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2011_I872172);
      
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

      void init() {

         // Initialise and register projections

        // the basic final-state projection: 
        // all final-state particles within 
        // the given eta acceptance
        const FinalState fs(Cuts::abseta < 4.9);

        // the final-state particles declared above are clustered using FastJet with
        // the anti-kT algorithm and a jet-radius parameter 0.4
        // muons and neutrinos are excluded from the clustering
        FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
        declare(jetfs, "jets");

        // FinalState of prompt photons and bare muons and electrons in the event
        //PromptFinalState photons(Cuts::abspid == PID::PHOTON);
        ChargedFinalState photons(Cuts::abspid == PID::PHOTON);
        declare(photons, "photons");
        PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

        // dress the prompt bare leptons with prompt photons within dR < 0.1
        // apply some fiducial cuts on the dressed leptons
        Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
        DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
        declare(dressed_leps, "leptons");

        // missing momentum
        declare(MissingMomentum(fs), "MET");


        // the basic final-state projection: all final-state particles within the given eta acceptance
        const ChargedFinalState cfs(Cuts::pT > 1*GeV); //Not cutting in eta, so no need to correct for pair acceptance
        declare(cfs, "CFS");
        const ChargedFinalState cfsTrig(Cuts::abseta < 1.0 && Cuts::pT > 2*GeV);
        declare(cfsTrig, "CFSTrig");
        
        declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

        //==================================================
        // Create one correlator for each set of Collisions System / Beam Energy / Centrality Interval / Trigger pT interval / Associated pT interval
        // The number used in the constructor is the y-axis index of the histograms (d00-x00-y00) of the paper
        // Ex.: Correlator c1(1); -> is the correlator for histograms _h["0311"], _h["0411"], etc
        // Ex.: Correlator c2(2); -> is the correlator for histograms _h["0312"], _h["0412"], etc
        //==================================================
        
        /*
        // Fig 2. "cfs" phi 0-15 and 75-90
        Correlator d1x1y1(1);
        d1x1y1.SetCollSystemAndEnergy("AuAu200GeV");
        d1x1y1.SetCentrality(0.0, 20.0);
        d1x1y1.SetTriggerRange(4.0, 7.0);
        d1x1y1.SetAssociatedRange(3.0, 4.0);
        d1x1y1.SetAzimuthalRange(0.0, 15.0);
        Correlators.push_back(d1x1y1);
        
        Correlator d1x1y2(2);
        d1x1y2.SetCollSystemAndEnergy("AuAu200GeV");
        d1x1y2.SetCentrality(20.0, 60.0);
        d1x1y2.SetTriggerRange(4.0, 7.0);
        d1x1y2.SetAssociatedRange(3.0, 4.0);
        d1x1y2.SetAzimuthalRange(0.0, 15.0);
        Correlators.push_back(d1x1y2);
        
        Correlator d1x1y3(3);
        d1x1y3.SetCollSystemAndEnergy("AuAu200GeV");
        d1x1y3.SetCentrality(0.0, 20.0);
        d1x1y3.SetTriggerRange(4.0, 7.0);
        d1x1y3.SetAssociatedRange(3.0, 4.0);
        Correlators.push_back(d1x1y3);
        
        Correlator d1x1y4(4);
        d1x1y4.SetCollSystemAndEnergy("AuAu200GeV");
        d1x1y4.SetCentrality(0.0, 20.0);
        d1x1y4.SetTriggerRange(4.0, 7.0);
        d1x1y4.SetAssociatedRange(3.0, 4.0);
        Correlators.push_back(d1x1y4);  
        
        // Fig 3. "cfs_ave" phi 0-90
        Correlator d2x1y1(1);
        d2x1y1.SetCollSystemAndEnergy("AuAu200GeV");
        d2x1y1.SetCentrality(0.0, 20.0);
        d2x1y1.SetTriggerRange(4.0, 7.0);
        d2x1y1.SetAssociatedRange(3.0, 4.0);
        d2x1y2.SetAzimuthalRange(0.0, 90.0);
        Correlators.push_back(d2x1y1);
        
        Correlator d2x1y2(2);
        d2x1y2.SetCollSystemAndEnergy("AuAu200GeV");
        d2x1y2.SetCentrality(20.0, 60.0);
        d2x1y2.SetTriggerRange(4.0, 7.0);
        d2x1y2.SetAssociatedRange(3.0, 4.0);
        d2x1y2.SetAzimuthalRange(0.0, 90.0);
        Correlators.push_back(d2x1y2);
        */
        // Fig 13. "cfsall_00-20 3-4"
        Correlator d3x1y1(311);
        d3x1y1.SetCollSystemAndEnergy("AuAu200GeV");
        d3x1y1.SetCentrality(0.0, 100.0);
        d3x1y1.SetTriggerRange(0.0, 7.0);
        d3x1y1.SetAssociatedRange(0.0, 7.0);
        d3x1y1.SetAzimuthalRange(0.0, 15.0);
        Correlators.push_back(d3x1y1);
        
        Correlator d3x1y2(312);
        d3x1y2.SetCollSystemAndEnergy("AuAu200GeV");
        d3x1y2.SetCentrality(0.0, 20.0);
        d3x1y2.SetTriggerRange(4.0, 7.0);
        d3x1y2.SetAssociatedRange(3.0, 4.0);
        d3x1y2.SetAzimuthalRange(15.0, 30.0);
        Correlators.push_back(d3x1y2);
        
        Correlator d3x1y3(313);
        d3x1y3.SetCollSystemAndEnergy("AuAu200GeV");
        d3x1y3.SetCentrality(0.0, 20.0);
        d3x1y3.SetTriggerRange(4.0, 7.0);
        d3x1y3.SetAssociatedRange(3.0, 4.0);
        d3x1y3.SetAzimuthalRange(30.0, 45.0);
        Correlators.push_back(d3x1y3);
        
        Correlator d3x1y4(314);
        d3x1y4.SetCollSystemAndEnergy("AuAu200GeV");
        d3x1y4.SetCentrality(0.0, 20.0);
        d3x1y4.SetTriggerRange(4.0, 7.0);
        d3x1y4.SetAssociatedRange(3.0, 4.0);
        d3x1y4.SetAzimuthalRange(45.0, 60.0);
        Correlators.push_back(d3x1y4);
        
        Correlator d3x1y5(315);
        d3x1y5.SetCollSystemAndEnergy("AuAu200GeV");
        d3x1y5.SetCentrality(0.0, 20.0);
        d3x1y5.SetTriggerRange(4.0, 7.0);
        d3x1y5.SetAssociatedRange(3.0, 4.0);
        d3x1y5.SetAzimuthalRange(60.0, 75.0);
        Correlators.push_back(d3x1y5);
        
        Correlator d3x1y6(316);
        d3x1y6.SetCollSystemAndEnergy("AuAu200GeV");
        d3x1y6.SetCentrality(0.0, 20.0);
        d3x1y6.SetTriggerRange(4.0, 7.0);
        d3x1y6.SetAssociatedRange(3.0, 4.0);
        d3x1y6.SetAzimuthalRange(75.0, 90.0);
        Correlators.push_back(d3x1y6);
        
        // Fig 14. cfsall_20-60 3-4
        Correlator d4x1y1(411);
        d4x1y1.SetCollSystemAndEnergy("AuAu200GeV");
        d4x1y1.SetCentrality(0.0, 20.0);
        d4x1y1.SetTriggerRange(4.0, 7.0);
        d4x1y1.SetAssociatedRange(4.0, 5.0);
        d4x1y1.SetAzimuthalRange(0.0, 15.0);
        Correlators.push_back(d4x1y1);
        
        Correlator d4x1y2(412);
        d4x1y2.SetCollSystemAndEnergy("AuAu200GeV");
        d4x1y2.SetCentrality(0.0, 20.0);
        d4x1y2.SetTriggerRange(4.0, 7.0);
        d4x1y2.SetAssociatedRange(4.0, 5.0);
        d4x1y2.SetAzimuthalRange(15.0, 30.0);
        Correlators.push_back(d4x1y2);
        
        Correlator d4x1y3(413);
        d4x1y3.SetCollSystemAndEnergy("AuAu200GeV");
        d4x1y3.SetCentrality(0.0, 20.0);
        d4x1y3.SetTriggerRange(4.0, 7.0);
        d4x1y3.SetAssociatedRange(4.0, 5.0);
        d4x1y3.SetAzimuthalRange(30.0, 45.0);
        Correlators.push_back(d4x1y3);
        
        Correlator d4x1y4(414);
        d4x1y4.SetCollSystemAndEnergy("AuAu200GeV");
        d4x1y4.SetCentrality(0.0, 20.0);
        d4x1y4.SetTriggerRange(4.0, 7.0);
        d4x1y4.SetAssociatedRange(4.0, 5.0);
        d4x1y4.SetAzimuthalRange(45.0, 60.0);
        Correlators.push_back(d4x1y4);
        
        Correlator d4x1y5(415);
        d4x1y5.SetCollSystemAndEnergy("AuAu200GeV");
        d4x1y5.SetCentrality(0.0, 20.0);
        d4x1y5.SetTriggerRange(4.0, 7.0);
        d4x1y5.SetAssociatedRange(4.0, 5.0);
        d4x1y5.SetAzimuthalRange(60.0, 75.0);
        Correlators.push_back(d4x1y5);

        Correlator d4x1y6(416);
        d4x1y6.SetCollSystemAndEnergy("AuAu200GeV");
        d4x1y6.SetCentrality(0.0, 20.0);
        d4x1y6.SetTriggerRange(4.0, 7.0);
        d4x1y6.SetAssociatedRange(4.0, 5.0);
        d4x1y6.SetAzimuthalRange(75.0, 90.0);
        Correlators.push_back(d4x1y6);
        
        // Fig 13. cfsall_00-20 5-7
        Correlator d5x1y1(511);
        d5x1y1.SetCollSystemAndEnergy("AuAu200GeV");
        d5x1y1.SetCentrality(0.0, 20.0);
        d5x1y1.SetTriggerRange(4.0, 7.0);
        d5x1y1.SetAssociatedRange(5.0, 7.0);
        d5x1y1.SetAzimuthalRange(0.0, 15.0);
        Correlators.push_back(d5x1y1);

        Correlator d5x1y2(512);
        d5x1y2.SetCollSystemAndEnergy("AuAu200GeV");
        d5x1y2.SetCentrality(0.0, 20.0);
        d5x1y2.SetTriggerRange(4.0, 7.0);
        d5x1y2.SetAssociatedRange(5.0, 7.0);
        d5x1y2.SetAzimuthalRange(15.0, 30.0);
        Correlators.push_back(d5x1y2);

        Correlator d5x1y3(513);
        d5x1y3.SetCollSystemAndEnergy("AuAu200GeV");
        d5x1y3.SetCentrality(0.0, 20.0);
        d5x1y3.SetTriggerRange(4.0, 7.0);
        d5x1y3.SetAssociatedRange(5.0, 7.0);
        d5x1y3.SetAzimuthalRange(30.0, 45.0);
        Correlators.push_back(d5x1y3);

        Correlator d5x1y4(514);
        d5x1y4.SetCollSystemAndEnergy("AuAu200GeV");
        d5x1y4.SetCentrality(0.0, 20.0);
        d5x1y4.SetTriggerRange(4.0, 7.0);
        d5x1y4.SetAssociatedRange(5.0, 7.0);
        d5x1y4.SetAzimuthalRange(45.0, 60.0);
        Correlators.push_back(d5x1y4);

        Correlator d5x1y5(515);
        d5x1y5.SetCollSystemAndEnergy("AuAu200GeV");
        d5x1y5.SetCentrality(0.0, 20.0);
        d5x1y5.SetTriggerRange(4.0, 7.0);
        d5x1y5.SetAssociatedRange(5.0, 7.0);
        d5x1y5.SetAzimuthalRange(60.0, 75.0);
        Correlators.push_back(d5x1y5);

        Correlator d5x1y6(516);
        d5x1y6.SetCollSystemAndEnergy("AuAu200GeV");
        d5x1y6.SetCentrality(0.0, 20.0);
        d5x1y6.SetTriggerRange(4.0, 7.0);
        d5x1y6.SetAssociatedRange(5.0, 7.0);
        d5x1y6.SetAzimuthalRange(75.0, 90.0);
        Correlators.push_back(d5x1y6);
        
        // Fig 13. cfsall_00-20 4-5
        Correlator d6x1y1(611);
        d6x1y1.SetCollSystemAndEnergy("AuAu200GeV");
        d6x1y1.SetCentrality(20.0, 60.0);
        d6x1y1.SetTriggerRange(4.0, 7.0);
        d6x1y1.SetAssociatedRange(3.0, 4.0);
        d6x1y1.SetAzimuthalRange(0.0, 15.0);
        Correlators.push_back(d6x1y1);

        Correlator d6x1y2(612);
        d6x1y2.SetCollSystemAndEnergy("AuAu200GeV");
        d6x1y2.SetCentrality(20.0, 60.0);
        d6x1y2.SetTriggerRange(4.0, 7.0);
        d6x1y2.SetAssociatedRange(3.0, 4.0);
        d6x1y2.SetAzimuthalRange(15.0, 30.0);
        Correlators.push_back(d6x1y2);

        Correlator d6x1y3(613);
        d6x1y3.SetCollSystemAndEnergy("AuAu200GeV");
        d6x1y3.SetCentrality(20.0, 60.0);
        d6x1y3.SetTriggerRange(4.0, 7.0);
        d6x1y3.SetAssociatedRange(3.0, 4.0);
        d6x1y3.SetAzimuthalRange(30.0, 45.0);
        Correlators.push_back(d6x1y3);

        Correlator d6x1y4(614);
        d6x1y4.SetCollSystemAndEnergy("AuAu200GeV");
        d6x1y4.SetCentrality(20.0, 60.0);
        d6x1y4.SetTriggerRange(4.0, 7.0);
        d6x1y4.SetAssociatedRange(3.0, 4.0);
        d6x1y4.SetAzimuthalRange(45.0, 60.0);
        Correlators.push_back(d6x1y4);

        Correlator d6x1y5(615);
        d6x1y5.SetCollSystemAndEnergy("AuAu200GeV");
        d6x1y5.SetCentrality(20.0, 60.0);
        d6x1y5.SetTriggerRange(4.0, 7.0);
        d6x1y5.SetAssociatedRange(3.0, 4.0);
        d6x1y5.SetAzimuthalRange(60.0, 75.0);
        Correlators.push_back(d6x1y5);

        Correlator d6x1y6(616);
        d6x1y6.SetCollSystemAndEnergy("AuAu200GeV");
        d6x1y6.SetCentrality(20.0, 60.0);
        d6x1y6.SetTriggerRange(4.0, 7.0);
        d6x1y6.SetAssociatedRange(3.0, 4.0);
        d6x1y6.SetAzimuthalRange(75.0, 90.0);
        Correlators.push_back(d6x1y6);
        
        // Fig 14. cfs_20-60 4-5
        Correlator d7x1y1(711);
        d7x1y1.SetCollSystemAndEnergy("AuAu200GeV");
        d7x1y1.SetCentrality(20.0, 60.0);
        d7x1y1.SetTriggerRange(4.0, 7.0);
        d7x1y1.SetAssociatedRange(4.0, 5.0);
        d7x1y1.SetAzimuthalRange(0.0, 15.0);
        Correlators.push_back(d7x1y1);

        Correlator d7x1y2(712);
        d7x1y2.SetCollSystemAndEnergy("AuAu200GeV");
        d7x1y2.SetCentrality(20.0, 60.0);
        d7x1y2.SetTriggerRange(4.0, 7.0);
        d7x1y2.SetAssociatedRange(4.0, 5.0);
        d7x1y2.SetAzimuthalRange(15.0, 30.0);
        Correlators.push_back(d7x1y2);

        Correlator d7x1y3(713);
        d7x1y3.SetCollSystemAndEnergy("AuAu200GeV");
        d7x1y3.SetCentrality(20.0, 60.0);
        d7x1y3.SetTriggerRange(4.0, 7.0);
        d7x1y3.SetAssociatedRange(4.0, 5.0);
        d7x1y3.SetAzimuthalRange(30.0, 45.0);
        Correlators.push_back(d7x1y3);

        Correlator d7x1y4(714);
        d7x1y4.SetCollSystemAndEnergy("AuAu200GeV");
        d7x1y4.SetCentrality(20.0, 60.0);
        d7x1y4.SetTriggerRange(4.0, 7.0);
        d7x1y4.SetAssociatedRange(4.0, 5.0);
        d7x1y4.SetAzimuthalRange(45.0, 60.0);
        Correlators.push_back(d7x1y4);

        Correlator d7x1y5(715);
        d7x1y5.SetCollSystemAndEnergy("AuAu200GeV");
        d7x1y5.SetCentrality(20.0, 60.0);
        d7x1y5.SetTriggerRange(4.0, 7.0);
        d7x1y5.SetAssociatedRange(4.0, 5.0);
        d7x1y5.SetAzimuthalRange(60.0, 75.0);
        Correlators.push_back(d7x1y5);

        Correlator d7x1y6(716);
        d7x1y6.SetCollSystemAndEnergy("AuAu200GeV");
        d7x1y6.SetCentrality(20.0, 60.0);
        d7x1y6.SetTriggerRange(4.0, 7.0);
        d7x1y6.SetAssociatedRange(4.0, 5.0);
        d7x1y6.SetAzimuthalRange(75.0, 90.0);
        Correlators.push_back(d7x1y6);
        
        // Fig 14. cfs_20-60 5-7
        Correlator d8x1y1(811);
        d8x1y1.SetCollSystemAndEnergy("AuAu200GeV");
        d8x1y1.SetCentrality(20.0, 60.0);
        d8x1y1.SetTriggerRange(4.0, 7.0);
        d8x1y1.SetAssociatedRange(5.0, 7.0);
        d8x1y1.SetAzimuthalRange(0.0, 15.0);
        Correlators.push_back(d8x1y1);

        Correlator d8x1y2(812);
        d8x1y2.SetCollSystemAndEnergy("AuAu200GeV");
        d8x1y2.SetCentrality(20.0, 60.0);
        d8x1y2.SetTriggerRange(4.0, 7.0);
        d8x1y2.SetAssociatedRange(5.0, 7.0);
        d8x1y2.SetAzimuthalRange(15.0, 30.0);
        Correlators.push_back(d8x1y2);

        Correlator d8x1y3(813);
        d8x1y3.SetCollSystemAndEnergy("AuAu200GeV");
        d8x1y3.SetCentrality(20.0, 60.0);
        d8x1y3.SetTriggerRange(4.0, 7.0);
        d8x1y3.SetAssociatedRange(5.0, 7.0);
        d8x1y3.SetAzimuthalRange(30.0, 45.0);
        Correlators.push_back(d8x1y3);

        Correlator d8x1y4(814);
        d8x1y4.SetCollSystemAndEnergy("AuAu200GeV");
        d8x1y4.SetCentrality(20.0, 60.0);
        d8x1y4.SetTriggerRange(4.0, 7.0);
        d8x1y4.SetAssociatedRange(5.0, 7.0);
        d8x1y4.SetAzimuthalRange(45.0, 60.0);
        Correlators.push_back(d8x1y4);

        Correlator d8x1y5(815);
        d8x1y5.SetCollSystemAndEnergy("AuAu200GeV");
        d8x1y5.SetCentrality(20.0, 60.0);
        d8x1y5.SetTriggerRange(4.0, 7.0);
        d8x1y5.SetAssociatedRange(5.0, 7.0);
        d8x1y5.SetAzimuthalRange(60.0, 75.0);
        Correlators.push_back(d8x1y5);

        Correlator d8x1y6(816);
        d8x1y6.SetCollSystemAndEnergy("AuAu200GeV");
        d8x1y6.SetCentrality(20.0, 60.0);
        d8x1y6.SetTriggerRange(4.0, 7.0);
        d8x1y6.SetAssociatedRange(5.0, 7.0);
        d8x1y6.SetAzimuthalRange(75.0, 90.0);
        Correlators.push_back(d8x1y6);

        for(unsigned int i = 1; i <= Correlators.size(); i++)
        {
          //book(sow[i],"sow" + to_string(i));
          book(sow[Correlators[i].GetIndex()], "sow" + to_string(Correlators[i].GetIndex()));
          //book(_h[to_string(Correlators[i].GetIndex())]);
          //book(_h["031" + to_string(i)], 3, 1, i);
          //book(_h["041" + to_string(i)], 4, 1, i);
        }
        
        // Fig 13
        book(_h["0311"], 3, 1, 1);
        book(_h["0312"], 3, 1, 2);
        book(_h["0313"], 3, 1, 3);
        book(_h["0314"], 3, 1, 4);
        book(_h["0315"], 3, 1, 5);
        book(_h["0316"], 3, 1, 6);
        book(_h["0411"], 4, 1, 1);
        book(_h["0412"], 4, 1, 2);
        book(_h["0413"], 4, 1, 3);
        book(_h["0414"], 4, 1, 4);
        book(_h["0415"], 4, 1, 5);
        book(_h["0416"], 4, 1, 6);
        book(_h["0511"], 5, 1, 1);
        book(_h["0512"], 5, 1, 2);
        book(_h["0513"], 5, 1, 3);
        book(_h["0514"], 5, 1, 4);
        book(_h["0515"], 5, 1, 5);
        book(_h["0516"], 5, 1, 6);
        book(_h["0611"], 6, 1, 1);
        book(_h["0612"], 6, 1, 2);
        book(_h["0613"], 6, 1, 3);
        book(_h["0614"], 6, 1, 4);
        book(_h["0615"], 6, 1, 5);
        book(_h["0616"], 6, 1, 6);
        book(_h["0711"], 7, 1, 1);
        book(_h["0712"], 7, 1, 2);
        book(_h["0713"], 7, 1, 3);
        book(_h["0714"], 7, 1, 4);
        book(_h["0715"], 7, 1, 5);
        book(_h["0716"], 7, 1, 6);
        book(_h["0811"], 8, 1, 1);
        book(_h["0812"], 8, 1, 2);
        book(_h["0813"], 8, 1, 3);
        book(_h["0814"], 8, 1, 4);
        book(_h["0815"], 8, 1, 5);
        book(_h["0816"], 8, 1, 6);
        
        for (vector<Correlator>::size_type i = 0; i < Correlators.size(); i++) {
          nEvents[Correlators[i].GetIndex()] = 0;
          nTriggers[Correlators[i].GetIndex()] = 0;
        }

      }

      void analyze(const Event& event) {

        const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
        const ChargedFinalState& cfsTrig = apply<ChargedFinalState>(event, "CFSTrig");
        
        // Determine the beam system and energy being used. 
        double nNucleons = 0;
        string collSystem;
        const ParticlePair& beam = beams();
        
        if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) {
          collSystem = "AuAu";
          nNucleons = 197.0;
        } else if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000290630) {
          collSystem = "CuCu";
          nNucleons = 63.0;
        } else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) {
          collSystem = "pp";
          nNucleons = 1.0;
        }
        
        if(collSystem.empty()) return;
      
        string cmsEnergy;
        if (fuzzyEquals(sqrtS()/GeV, 200*nNucleons, 1E-3)) cmsEnergy = "200GeV";
        if (fuzzyEquals(sqrtS()/GeV, 62.3*nNucleons, 1E-3)) cmsEnergy = "62GeV";
        if(cmsEnergy.empty()) return;
      
        string SysAndEnergy = collSystem + cmsEnergy;
        if (SysAndEnergy != "AuAu200GeV") {
          cerr << "This is the wrong beam type!" << endl;
          vetoEvent;
        }
        
        double triggerptMin = 999.;
        double triggerptMax = -999.;
        double associatedptMin = 999.;
        double associatedptMax = -999.;
        
        // Prepare centrality projection and value
        const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
        double centr = centProj();
        
        // Check if the event is in the centrality range we want to search in.
        if (centr < 0.0 && centr > 60.0) {
          vetoEvent;
        }
        
        bool isVeto = true;
    
        // Validate the events for the correlators. 
        for (Correlator& corr : Correlators) {
          if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
          if(!corr.CheckCentrality(centr)) continue;
        
          //If event is accepted for the correlator, fill event weights
          sow[corr.GetIndex()]->fill();
          nEvents[corr.GetIndex()]++;
        
          isVeto = false;
        
          //Check min and max of the trigger and associated particles in order to speed up the particle loops
          if (corr.GetTriggerRangeMin() < triggerptMin) triggerptMin = corr.GetTriggerRangeMin();
          if (corr.GetTriggerRangeMax() > triggerptMax) triggerptMax = corr.GetTriggerRangeMax();
        
          if (corr.GetAssociatedRangeMin() < associatedptMin) associatedptMin = corr.GetAssociatedRangeMin();
          if (corr.GetAssociatedRangeMax() > associatedptMax) associatedptMax = corr.GetAssociatedRangeMax();
        }
        
        if (isVeto) vetoEvent;
        
        // loop over charged final state particles
        for (const Particle& pTrig : cfsTrig.particles()) {
        
          if (pTrig.pt()/GeV < triggerptMin || pTrig.pt()/GeV > triggerptMax) continue;
          if (isSecondary(pTrig)) continue;
          
          if (abs(pTrig.pid()) == 211 || abs(pTrig.pid()) == 2212 || abs(pTrig.pid()) == 321) {
          
            for (Correlator& corr : Correlators) {
              if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
              nTriggers[corr.GetIndex()]++;
            }
            
            for (const Particle& pAssoc : cfs.particles()) {
              if(pAssoc.pt()/GeV < associatedptMin || pAssoc.pt()/GeV > pTrig.pt()/GeV) continue;
              if(isSameParticle(pTrig,pAssoc)) continue;
              if(isSecondary(pAssoc)) continue;
              
              // https://home.fnal.gov/~mrenna/lutp0613man2/node44.html
              // 211 = pi+, 2212 = p+, 321 = K+
              if(abs(pAssoc.pid()) == 211 || abs(pAssoc.pid()) == 2212 || abs(pAssoc.pid()) == 321) {
                //int mybin = GetTrigBin(pTrig.pt());
                //int mybina = GetAssocBin(pAssoc.pt());

                //https://rivet.hepforge.org/code/dev/structRivet_1_1DeltaPhiInRange.html
                double dPhi = deltaPhi(pTrig, pAssoc, true); //this does NOT rotate the delta phi to be in a given range
                double dEta = deltaEta(pTrig, pAssoc);
              
                for (Correlator& corr : Correlators) {
                  if(!corr.CheckConditionsMaxTrigger(SysAndEnergy, centr, pTrig.pt()/GeV, pAssoc.pt()/GeV)) continue;
                    
                  if (abs(dEta) < 1.78) {
                    _h[to_string(corr.GetIndex())]->fill(-abs(dPhi), 0.5);
                    _DeltaPhi[corr.GetIndex()]->fill(abs(dPhi), 0.5);
                    _DeltaPhiSub[corr.GetIndex()]->fill(abs(dPhi), 0.5);
                  }
                }
              }
            } // End of pAssoc loop.
          }
        } // End of pTrig loop.
      } // End of analysis()

      void finalize() {
      
        for (unsigned int i = 1; i <= Correlators.size(); i++) {
          int index = Correlators[i].GetIndex();
          if (nTriggers[index] > 0) {
            _h["0" + to_string(index)]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));
            _DeltaPhi[index]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));
            _DeltaPhiSub[index]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));
          }
            
        }

      } // End of finalize()

    
  }; // End of PHENIX_2011_I872172

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PHENIX_2011_I872172);
}
