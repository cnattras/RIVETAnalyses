// Andrew Bryant -- abryan39
// 6 December 2019
//
// This code implements an analysis for the 1010.1521 paper. Much of the code
// is adapted from the 1010.5800 analysis and tweaked for use in this analysis. 
// https://github.com/cnattras/RIVETAnalyses/tree/master/1110.5800

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
      // Unlike the code found in the 1110.5800 analysis, we use a map instead of a vector.
      // This is because y-axis values are not unique and using the full histogram index is easier.
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
        // Create a correlator for each histogram. A different correlator must be created for each beam, centrality, and range.
        // The number used in the constructor is a representation of the histograms d00-x00-y00 format.
        // Ex. index 311 corresponds to histogram d03-x01-y01.      
        //==================================================
        
        // Fig 13. "cfsall_00-20 3-4"
        for (int index = 311; index < 317; index++) {
          for (double azimuthal = 0.0; azimuthal < 90.0; azimuthal += 15.0) {
            Correlator corr(index);
            corr.SetCollSystemAndEnergy("AuAu200GeV");
            corr.SetCentrality(0.0, 20.0);
            corr.SetTriggerRange(4.0, 7.0);
            corr.SetAssociatedRange(3.0, 4.0);
            corr.SetAzimuthalRange(azimuthal, azimuthal+15.0);
            Correlators.push_back(corr);
          }
        }      
        
        // Fig 13. cfsall_00-20 4-5
        for (int index = 411; index < 417; index++) {
          for (double azimuthal = 0.0; azimuthal < 90.0; azimuthal += 15.0) {
            Correlator corr(index);
            corr.SetCollSystemAndEnergy("AuAu200GeV");
            corr.SetCentrality(0.0, 20.0);
            corr.SetTriggerRange(4.0, 7.0);
            corr.SetAssociatedRange(4.0, 5.0);
            corr.SetAzimuthalRange(azimuthal, azimuthal+15.0);
            Correlators.push_back(corr);
          }
        }
        
        // Fig 13. cfsall_00-20 5-7
        for (int index = 511; index < 517; index++) {
          for (double azimuthal = 0.0; azimuthal < 90.0; azimuthal += 15.0) {
            Correlator corr(index);
            corr.SetCollSystemAndEnergy("AuAu200GeV");
            corr.SetCentrality(0.0, 20.0);
            corr.SetTriggerRange(4.0, 7.0);
            corr.SetAssociatedRange(5.0, 7.0);
            corr.SetAzimuthalRange(azimuthal, azimuthal+15.0);
            Correlators.push_back(corr);
          }
        }
        
        // Fig 14. cfsall_20-60 3-4
        for (int index = 611; index < 617; index++) {
          for (double azimuthal = 0.0; azimuthal < 90.0; azimuthal += 15.0) {
            Correlator corr(index);
            corr.SetCollSystemAndEnergy("AuAu200GeV");
            corr.SetCentrality(20.0, 60.0);
            corr.SetTriggerRange(4.0, 7.0);
            corr.SetAssociatedRange(3.0, 4.0);
            corr.SetAzimuthalRange(azimuthal, azimuthal+15.0);
            Correlators.push_back(corr);
          }
        }
        
        // Fig 14. cfs_20-60 4-5
        for (int index = 711; index < 717; index++) {
          for (double azimuthal = 0.0; azimuthal < 90.0; azimuthal += 15.0) {
            Correlator corr(index);
            corr.SetCollSystemAndEnergy("AuAu200GeV");
            corr.SetCentrality(20.0, 60.0);
            corr.SetTriggerRange(4.0, 7.0);
            corr.SetAssociatedRange(4.0, 5.0);
            corr.SetAzimuthalRange(azimuthal, azimuthal+15.0);
            Correlators.push_back(corr);
          }
        }
        
        // Fig 14. cfs_20-60 5-7        
        for (int index = 811; index < 817; index++) {
          for (double azimuthal = 0.0; azimuthal < 90.0; azimuthal += 15.0) {
            Correlator corr(index);
            corr.SetCollSystemAndEnergy("AuAu200GeV");
            corr.SetCentrality(20.0, 60.0);
            corr.SetTriggerRange(4.0, 7.0);
            corr.SetAssociatedRange(5.0, 7.0);
            corr.SetAzimuthalRange(azimuthal, azimuthal+15.0);
            Correlators.push_back(corr);
          }
        }
        
        // Create sum_of_weights and delta phi books that can be filled.
        for(vector<Correlator>::size_type i = 0; i < Correlators.size(); i++) {
          int index = Correlators[i].GetIndex();
          book(sow[index], "sow" + to_string(index));
          book(_DeltaPhi[index], "DeltaPhi" + to_string(index), 24, 0, M_PI);
        }
        
        //*******************************
        // Book histrograms to be filled.
        //*******************************
        
        // y-axis ranges of the histograms to be booked.
        // The position of the pair in the vector + 1 is the data set.
        // histogram_ranges[0] corresponds to d01-x01-y01 - d01-x01-y04.
        vector<pair<int, int>> histogram_ranges = { make_pair(1, 4), make_pair(1, 2), make_pair(1, 6), 
                                                    make_pair(1, 6), make_pair(1, 6), make_pair(1, 6),
                                                    make_pair(1, 6), make_pair(1, 6), make_pair(1, 2),
                                                    make_pair(1, 2), make_pair(1, 4), make_pair(1, 4),
                                                    make_pair(1, 4), make_pair(1, 2), make_pair(1, 4),
                                                    make_pair(1, 4), make_pair(1, 6), make_pair(1, 6),
                                                    make_pair(1, 6), make_pair(1, 6), make_pair(1, 6),
                                                    make_pair(1, 6)
                                                  };
       
        for (vector<pair<int,int>>::size_type d = 0; d < histogram_ranges.size(); d++) {
          for (int y = histogram_ranges[d].first; y <= histogram_ranges[d].second; y++) {
            book(_h[to_string(static_cast<int>(d)+1) + "1" + to_string(y)], static_cast<int>(d)+1, 1, y);
          }
        }
        
        // Initialize the nEvents and nTriggers maps.
        for (vector<Correlator>::size_type i = 0; i < Correlators.size(); i++) {
          nEvents[Correlators[i].GetIndex()] = 0;
          nTriggers[Correlators[i].GetIndex()] = 0;
        }

      } // End of init()

      void analyze(const Event& event) {
      
        const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
        const ChargedFinalState& cfsTrig = apply<ChargedFinalState>(event, "CFSTrig");
        
        double triggerptMin = 999.;
        double triggerptMax = -999.;
        double associatedptMin = 999.;
        double associatedptMax = -999.;
        
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
        
        // Prepare centrality projection and value.
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
        
          //If the event is accepted for the correlator, fill event weights.
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
          
          // https://home.fnal.gov/~mrenna/lutp0613man2/node44.html
          // 211 = pi+, 2212 = p+, 321 = K+
          if (abs(pTrig.pid()) == 211 || abs(pTrig.pid()) == 2212 || abs(pTrig.pid()) == 321) {
          
            for (Correlator& corr : Correlators) {
              if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
              nTriggers[corr.GetIndex()]++;
            }
            
            for (const Particle& pAssoc : cfs.particles()) {
              if(pAssoc.pt()/GeV < associatedptMin || pAssoc.pt()/GeV > pTrig.pt()/GeV) continue;
              if(isSameParticle(pTrig,pAssoc)) continue;
              if(isSecondary(pAssoc)) continue;
              
              if(abs(pAssoc.pid()) == 211 || abs(pAssoc.pid()) == 2212 || abs(pAssoc.pid()) == 321) {

                //https://rivet.hepforge.org/code/dev/structRivet_1_1DeltaPhiInRange.html
                double dPhi = deltaPhi(pTrig, pAssoc, true); //this does NOT rotate the delta phi to be in a given range
                double dEta = deltaEta(pTrig, pAssoc);
              
                for (Correlator& corr : Correlators) {
                  if(!corr.CheckConditionsMaxTrigger(SysAndEnergy, centr, pTrig.pt()/GeV, pAssoc.pt()/GeV)) continue;

                  if (abs(dEta) < 1.78) {
                    _h[to_string(corr.GetIndex())]->fill(-abs(dPhi), 0.5);
                    _DeltaPhi[corr.GetIndex()]->fill(abs(dPhi), 0.5);
                  }
                }
              }
            } // End of pAssoc loop.
          }
        } // End of pTrig loop.    
      } // End of analysis()
      
      void finalize() { 
        for (vector<Correlator>::size_type i = 0; i < Correlators.size(); i++) {
          int index = Correlators[i].GetIndex();

          if (nTriggers[index] > 0) {
            _h[to_string(index)]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));
            _DeltaPhi[index]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));
          }  
        }
      } // End of finalize()    
    
  }; // End of PHENIX_2011_I872172

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PHENIX_2011_I872172);
}
