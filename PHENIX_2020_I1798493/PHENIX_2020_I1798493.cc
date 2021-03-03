// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

namespace Rivet {
  class Correlator {

    private:
      std::vector<int> _indices;
      string _collSystemAndEnergy;
      pair<double,double> _centrality;
      pair<double,double> _triggerRange;
      pair<double,double> _associatedRange;
      vector<int> _pid;
      bool _noCentrality = false;
      bool _noAssoc = false;
			int _nTriggers = 0;
			int _nEvents = 0;
			Histo1DPtr _deltaPhi;
			CounterPtr _counter;
    public:

      /// Constructor
      Correlator(int index0, int index1, int index2) {
        _indices = {index0, index1, index2};
      }

			Correlator(int index0, int index1) {
        _indices = {index0, index1};
      }

			Correlator(int index0) {
        _indices = {index0};
      }

			Correlator(std::vector<int> vindex) {
        _indices = vindex;
      }

      void SetCollSystemAndEnergy(string s){ _collSystemAndEnergy = s; }
      void SetCentrality(double cmin, double cmax){ _centrality = make_pair(cmin, cmax); }
      void SetNoCentrality(){ _noCentrality = true; }
      void SetNoAssoc(){ _noAssoc = true; }
      void SetTriggerRange(double tmin, double tmax){ _triggerRange = make_pair(tmin, tmax); }
      void SetAssociatedRange(double amin, double amax){ _associatedRange = make_pair(amin, amax); }
      void SetPID(std::initializer_list<int> pid){ _pid = pid; }
			void SetCorrelationFunction(Histo1DPtr cf){ _deltaPhi = cf; }
			void SetCounter(CounterPtr c){ _counter = c; }

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
			double GetWeight(){ return _counter->sumW(); }
			Histo1DPtr GetCorrelationFunction(){ return _deltaPhi; }
			CounterPtr GetCounter(){ return _counter; }

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

			void AddCorrelation(Particle pTrig, Particle pAssoc)
			{
					double dPhi = GetDeltaPhi(pTrig, pAssoc);

					_deltaPhi->fill(dPhi);
			}

			void AddWeight()
			{
					_counter->fill();
					_nEvents++;
			}

      int GetIndex(int i){ return _indices[i]; }
      string GetFullIndex()
      {
          string fullIndex = "";
					for(int index : _indices)
					{
							fullIndex += to_string(index);
					}
          return fullIndex;
      }

			void AddTrigger()
			{
					_nTriggers++;
			}

			void Normalize(double weight = 1.)
			{
					if(_nTriggers*_counter->sumW() > 0) _deltaPhi->scaleW((weight*_nEvents)/(_nTriggers*_counter->sumW()));
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
  class PHENIX_2020_I1798493 : public Analysis {
    public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2020_I1798493);
	

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV);
      declare(cfs, "CFS");
	    declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
	
	    // const FinalState fs(Cuts::abseta < 4.9);


      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");
      // Book histograms
      // specify custom binning
      //book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
      //book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      //book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["PerTriggerVsdphiAUAU20to24"], 1, 1, 1);
    	book(_h["PerTriggerVsdphiAUAU16to20"], 1, 1, 2);
 	    book(_h["PerTriggerVsdphiAUAU12to16"], 1, 1, 3);
    	book(_h["PerTriggerVsdphiAUAU8to12"], 1, 1, 4);
    	book(_h["PerTriggerVsdphiAUAU4to8"], 1, 1, 5);
    	book(_h["PerTriggerVsdphiAUAU0to4"], 1, 1, 6);
    	book(_h["PerTriggerVsdphidAU20to24"], 2, 1, 1);
	    book(_h["PerTriggerVsdphidAU16to20"], 2, 1, 2);
    	book(_h["PerTriggerVsdphidAU12to16"], 2, 1, 3);
	    book(_h["PerTriggerVsdphidAU8to12"], 2, 1,  4);
    	book(_h["PerTriggerVsdphidAU4to8"], 2, 1, 5);
	    book(_h["GammaDirhPertriggerVsXiAUAU"], 3, 1, 1);
    	book(_h["GammaDirhPertriggerVsXidAU"], 3, 1, 2);
	    book(_h["IAA"], 4, 1 , 1 );
	    book(_h["IdA"], 4,  1,  2);
	    book(_h["IAAVsXiDirectPhoton5to7GevC"], 5, 1, 1);
	    book(_h["IAAVsXiDirectPhoton7to9GevC"], 6, 1, 1);
    	book(_h["IAAVsXiDirectPhoton9to12GevC"], 7, 1, 1);
    	book(_h["IAAVsXiDirectPhoton5to7PiOver2"], 8, 1, 1);
    	book(_h["IAAVsXiDirectPhoton5to7PiOver3"], 8, 1, 2);
    	book(_h["IAAVsXiDirectPhoton5to7PiOver6"], 8, 1, 3);
    	book(_h["IAAVsXiDirectPhoton7to9PiOver2"], 9, 1, 1);
    	book(_h["IAAVsXiDirectPhoton7to9PiOver3"], 9, 1, 2);
    	book(_h["IAAVsXiDirectPhoton7to9PiOver6"], 9, 1, 3);
    	book(_h["IAAVsXiDirectPhoton9to12PiOver2"], 10, 1, 1);
    	book(_h["IAAVsXiDirectPhoton9to12PiOver3"], 10, 1, 2);
    	book(_h["IAAVsXiDirectPhoton9to12PiOver6"], 10, 1, 3);
    	book(_h["RatiosOfIAAVsDirectPhotonPtPiOver2"], 11, 1, 1);
    	book(_h["RatiosOfIAAVsDirectPhotonPtPiOver3"], 11, 1, 2);
    	book(_h["RatiosOfIAAVsDirectPhotonPtPiOver6"], 11, 1, 3);
    	book(_h["IAAVsXiDirectPhoton5to7GevC"], 12, 1, 1);
    	book(_h["IAAVsXiDirectPhoton7to9GevC"], 13, 1, 1);
    	book(_h["IAAVsXiDirectPhoton9to12GevC"], 14, 1, 1);
    	book(_c["sow_AuAu200_PerTriggerVsdphiAUAU20to24"], "sow_AuAu200_PerTriggerVsdphiAUAU20to24");
      book(_c["sow_AuAu200_PerTriggerVsdphiAUAU16to20"], "sow_AuAu200_PerTriggerVsdphiAUAU16to20");
      book(_c["sow_AuAu200_PerTriggerVsdphiAUAU12to16"], "sow_AuAu200_PerTriggerVsdphiAUAU12to16");
      book(_c["sow_AuAu200_PerTriggerVsdphiAUAU8to12"], "sow_AuAu200_PerTriggerVsdphiAUAU8to12");
      book(_c["sow_AuAu200_PerTriggerVsdphiAUAU4to8"], "sow_AuAu200_PerTriggerVsdphiAUAU4to8");
      book(_c["sow_AuAu200_PerTriggerVsdphiAUAU0to4"], "sow_AuAu200_PerTriggerVsdphiAUAU0to4");
      vector<Correlator> Correlators;
      Correlator corr1(1);
    	corr1.SetCollSystemAndEnergy("AuAu200GeV");
    	corr1.SetCentrality(0.,40.);
    	corr1.SetTriggerRange(5., 9.);
	    corr1.SetAssociatedRange(0.67, 1.);
    	corr1.SetCorrelationFunction(_h["PerTriggerVsdphiAUAU20to24"]);
	    corr1.SetCounter(_c["sow_AuAu200_PerTriggerVsdphiAUAU20to24"]);
    	Correlators.push_back(corr1); 

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve dressed leptons, sorted by pT
      vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

      // Remove all jets within dR < 0.2 of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
      Jets bjets = filter_select(jets, [](const Jet& jet) {
        return  jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
      });

      // Veto event if there are no b-jets
      if (bjets.empty())  vetoEvent;

      // Apply a missing-momentum cut
      if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent;

      // Fill histogram with leading b-jet pT
      // _h["XXXX"]->fill(bjets[0].pT()/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // normalize(_h["XXXX"]); // normalize to unity
      // normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
      // scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}

  vector<Correlator> Correlators;

  enum CollisionSystem {pp, AuAu, dAu};
  CollisionSystem collSys;
  };



  DECLARE_RIVET_PLUGIN(PHENIX_2020_I1798493);

}

