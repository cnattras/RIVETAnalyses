// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"

namespace Rivet {


  class Correlator {

    private:
      std::vector<int> _indices;
      string _collSystemAndEnergy;
      pair<double,double> _centrality;
      pair<double,double> _triggerRange;
      pair<double,double> _associatedRange;
      vector<int> _pid;
	int _underlyingEvent = 0;
      bool _noCentrality = false;
      bool _noAssoc = false;
	int _nTriggers = 0;
	int _nEvents = 0;
	Histo1DPtr _deltaPhi;
	CounterPtr _counter;
	

    public:

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
	void SetUnderlyingEvent(int underEvent){ _underlyingEvent = underEvent; } 

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
  class PHENIX_2018_I1658594 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2018_I1658594);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const ChargedFinalState cfs(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV);
      declare(cfs, "CFS");
      
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);
     
	//fig 6 
      book(_h["XXXXX"], "myh1", 20, 0.0, 100.0);
      book(_h["pi0Sepctra0To10V2"], 1, 1, 1);
      book(_h["pi0Spectra10to20V2"],1,1,2);
      book(_h["pi0Spectra20to30V2"],1,1,3);
      book(_h["pi0Spectra30to40V2"],1,1,4);
      book(_h["pi0Spectra40to50V2"],1,1,5);
      book(_h["pi0Sepctra0To10V3"],1,1,6);
      book(_h["pi0Spectra10to20V3"],1,1,7);
      book(_h["pi0Spectra20to30V3"],1,1,8);
      book(_h["pi0Spectra30to40V3"],1,1,9);
      book(_h["pi0Spectra40to50V3"],1,1,10);
      book(_h["pi0Sepctra0to10V4"], 1, 1, 11);
      book(_h["pi0Spectra10to20V4"],1,1,12);
      book(_h["pi0Spectra20to30V4"],1,1,13);
      book(_h["pi0Spectra30to40V4"],1,1,14);
      book(_h["pi0Spectra40to50V4"],1,1,15);
      book(_h["pi0Sepctra0To10V4psi2"],1,1,16);
      book(_h["pi0Spectra10to20V4psi2"],1,1,17);
      book(_h["pi0Spectra20to30V4psi2"],1,1,18);
      book(_h["pi0Spectra30to40V4psi2"],1,1,19);
      book(_h["pi0Spectra40to50V4psi2"],1,1,20);
	//fig 12
      book(_h["pi0Sepctra0to10_4-10x0.5-1"], 2, 1, 1);
      book(_h["pi0Spectra0to10_4-10x1-2"],2,1,2);
      book(_h["pi0Spectra0to10_4-10x2-4"],2,1,3);
      book(_h["pi0Spectra0to10_4-10x4-10"],2,1,4);
      book(_h["pi0Sepctra10to20_4-10x0.5-1"], 2, 1, 5);
      book(_h["pi0Spectra10to20_4-10x1-2"],2,1,6);
      book(_h["pi0Spectra10to20_4-10x2-4"],2,1,7);
      book(_h["pi0Spectra10to20_4-10x4-10"],2,1,8);
      book(_h["pi0Sepctra20to30_4-10x0.5-1"], 2, 1, 9);
      book(_h["pi0Spectra20to30_4-10x1-2"],2,1,10);
      book(_h["pi0Spectra20to30_4-10x2-4"],2,1,11);
      book(_h["pi0Spectra20to30_4-10x4-10"],2,1,12);
      book(_h["pi0Sepctra30to40_4-10x0.5-1"], 2, 1, 13);
      book(_h["pi0Spectra30to40_4-10x1-2"],2,1,14);
      book(_h["pi0Spectra30to40_4-10x2-4"],2,1,15);
      book(_h["pi0Spectra30to40_4-10x4-10"],2,1,16);
      book(_h["pi0Sepctra40to50_4-10x0.5-1"], 2, 1, 17);
      book(_h["pi0Spectra40to50_4-10x1-2"],2,1,18);
      book(_h["pi0Spectra40to50_4-10x2-4"],2,1,19);
      book(_h["pi0Spectra40to50_4-10x4-10"],2,1,20);
	//fig 15
      book(_h["pi0Sepctra0to10_1-2x0.5-1"], 3, 1, 1);
      book(_h["pi0Spectra0to10_1-2x1-2"],3,1,2);
      book(_h["pi0Spectra0to10_2-4x0.5-1"],3,1,3);
      book(_h["pi0Spectra0to10_2-4x1-2"],3,1,4);
      book(_h["pi0Sepctra0to10_2-4x2-4"], 3, 1, 5);
      book(_h["pi0Spectra10to20_1-2x0.5-1"],3,1,6);
      book(_h["pi0Spectra10to20_1-2x1-2"],3,1,7);
      book(_h["pi0Spectra10to20_2-4x0.5-1"],3,1,8);
      book(_h["pi0Sepctra10to20_2-4x1-2"], 3, 1, 9);
      book(_h["pi0Spectra10to20_2-4x2-4"],3,1,10);
      book(_h["pi0Spectra20to30_1-2x0.5-1"],3,1,11);
      book(_h["pi0Spectra20to30_1-2x1-2"],3,1,12);
      book(_h["pi0Sepctra20to30_2-4x0.5-1"], 3, 1, 13);
      book(_h["pi0Spectra20to30_2-4to1-2"],3,1,14);
      book(_h["pi0Spectra20to30_2-4x2-4"],3,1,15);
      book(_h["pi0Spectra30to40_1-2x0.5-1"],3,1,16);
      book(_h["pi0Sepctra30to40_1-2x1-2"], 3, 1, 17);
      book(_h["pi0Spectra30to40_2-4x2-4"],3,1,18);
      book(_h["pi0Spectra30to40_2-4x1-2"],3,1,19);
      book(_h["pi0Spectra30to40_2-4x2-4"],3,1,20);
      book(_h["pi0Spectra40to50_1-2x0.5-1"],3,1,21);
      book(_h["pi0Sepctra40to50_1-2x1-2"], 3, 1, 22);
      book(_h["pi0Spectra40to50_2-4to2-4"],3,1,23);
      book(_h["pi0Spectra40to50_2-4x1-2"],3,1,24);
      book(_h["pi0Spectra40to50_2-4x2-4"],3,1,25);
	//fig 18
      book(_h["pi0Sepctra0to10NearSide&2-4x1-2"], 4, 1, 1);
      book(_h["pi0Spectra0to10NearSide&2-4to2-4"],4,1,2);
      book(_h["pi0Spectra0to10NearSide&4-10x2-4"],4,1,3);
      book(_h["pi0Sepctra10to20NearSide&2-4x1-2"],4, 1, 4);
      book(_h["pi0Spectra10to20NearSide&2-4x2-4"],4,1,5);
      book(_h["pi0Spectra10to20NearSide&4-10x2-4"],4,1,6);
      book(_h["pi0Sepctra20to30NearSide&2-4x1-2"], 4, 1,7);
      book(_h["pi0Spectra20to30NearSide&2-4x2-4"],4,1,8);
      book(_h["pi0Spectra20to30NearSide&4-10x2-4"],4,1,9);
      book(_h["pi0Sepctra30to40NearSide&2-4x1-2"], 4, 1, 10);
      book(_h["pi0Spectra30to40NearSide&2-4x2-4"],4,1,11);
      book(_h["pi0Spectra30to40NearSide&4-10x2-4"],4,1,12);
      book(_h["pi0Sepctra40to50NearSide&2-4x1-2"], 4, 1, 13);
      book(_h["pi0Spectra40to50NearSide&2-4x2-4"],4,1,14);
      book(_h["pi0Spectra40to50NearSide&4-10x2-4"],4,1,15);
      book(_h["pi0Sepctra0to10FarSide&2-4x1-2"], 4, 1, 16);
      book(_h["pi0Spectra0to10FarSide&2-4to2-4"],4,1,17);
      book(_h["pi0Spectra0to10FarSide&4-10x2-4"],4,1,18);
      book(_h["pi0Sepctra10to20FarSide&2-4x1-2"],4, 1, 19);
      book(_h["pi0Spectra10to20FarSide&2-4x2-4"],4,1,20);
      book(_h["pi0Spectra10to20FarSide&4-10x2-4"],4,1,21);
      book(_h["pi0Sepctra20to30FarSide&2-4x1-2"], 4, 1, 22);
      book(_h["pi0Spectra20to30FarSide&2-4x2-4"],4,1,23);
      book(_h["pi0Spectra20to30FarSide&4-10x2-4"],4,1,24);
      book(_h["pi0Sepctra30to40FarSide&2-4x1-2"], 4, 1, 25);
      book(_h["pi0Spectra30to4FarSide&2-4x2-4"],4,1,26);
      book(_h["pi0Spectra30to40FarSide&4-10x2-4"],4,1,27);
      book(_h["pi0Sepctra40to50FarSide&2-4x1-2"], 4, 1, 28);
      book(_h["pi0Spectra40to50FarSide&2-4x2-4"],4,1,29);
      book(_h["pi0Spectra40to50FarSide&4-10x2-4"],4,1,30);

	//fig 12
	int minCent=0, maxCent=40, min_pT=4, max_pT=4, 
		min_pA=0.5, max_pA=4, minV=1, maxV=1;
	Correlator corrFig12(20);
	for(int a=minCent; a<=maxCent; a+=10){
		for(int b=min_pT; b<=max_pT; b*=2){
			for(int c=min_pA; c<=max_pA; c*=2){
				for(int d=minV; d<=maxV; d++){
					corrFig12.SetCollSystemAndEnergy("AuAu200GeV");
					corrFig12.SetCentrality(a,a+10);
					if(b==0.5 || b==1 || b==2) corrFig12.SetTriggerRange(b,b*2);
					else corrFig12.SetTriggerRange(b,10);
					if(c==0.5 || c==1 || c==2) corrFig12.SetAssociatedRange(c,c*2);
					else corrFig12.SetAssociatedRange(b,10);
					corrFig12.SetUnderlyingEvent(d);
					//corrFig12.SetCorrelatorFunction(_h["???"])
					//correFig12.SetCounter(_c["???"]);
				}
			}
		}
	}
	Correlators.push_back(corrFig12);
			

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

    }


    /// Normalise histograms etc., after the run
    void finalize() {
/*
  
*/
    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}

    vector<Correlator> Correlators;

    enum CollisionSystem {pp, AuAu};
    CollisionSystem collSys;

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2018_I1658594);

}
