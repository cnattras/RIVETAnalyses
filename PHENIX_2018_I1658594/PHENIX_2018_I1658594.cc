// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include <stdio.h>

namespace Rivet {


  class Correlator {

    private:
      std::vector<int> _indices;
      string _collSystemAndEnergy;
      pair<double,double> _centrality;
      pair<double,double> _triggerRange;
      pair<double,double> _associatedRange;
	pair<int,int> _RxnPlaneAngleRange;
      vector<int> _pid;
	int _RxnPlaneAngle = 0;
	int _eventPlaneMethod = 0;
      bool _noCentrality = false;
      bool _noAssoc = false;
	Histo1DPtr _deltaPhi;
	CounterPtr _counter;
        CounterPtr _cTriggers;


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
	void SetRxnPlaneAngle(int rxnMin, int rxnMax){ _RxnPlaneAngleRange = make_pair(rxnMin,rxnMax); }
	void SetEventPlaneMethod (int eventPlane) { _eventPlaneMethod = eventPlane; }
        void SetTriggerCounter(CounterPtr c){ _cTriggers = c; }

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
	pair<int,int> GetRxnPlaneAngle() { return _RxnPlaneAngleRange; }
	int GetEventPlaneMethod() { return _eventPlaneMethod; }
      vector<int> GetPID(){ return _pid; }
	double GetWeight(){ return _counter->sumW(); }
	Histo1DPtr GetCorrelationFunction(){ return _deltaPhi; }
	CounterPtr GetCounter(){ return _counter; }

      double GetDeltaPhi(Particle pAssoc, Particle pTrig)
  	   {
	double dPhi = deltaPhi(pTrig, pAssoc, true);
		//this does NOT rotate the delta phi to be in a given range

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
		_cTriggers->fill();
	}

      void Normalize(double weight = 1.)
	{
		if(_cTriggers->effNumEntries()*_counter->sumW() > 0) _deltaPhi->scaleW((weight*_counter->effNumEntries())/(_cTriggers->effNumEntries()*_counter->sumW()));
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
      //const FinalState fs(Cuts::abseta < 4.9);

	//fig 6
      //book(_h["XXXXX"], "myh1", 20, 0.0, 100.0);
      /*book(_h["pi0Sepctra0To10V2"], 1, 1, 1);
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
      book(_h["pi0Spectra40to50V4psi2"],1,1,20);*/
	//fig 12
      /*book(_h["pi0Sepctra0to10_4-10x0.5-1"], 2, 1, 1);
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
      book(_h["pi0Spectra40to50_4-10x4-10"],2,1,20);*/
	//fig 15
      /*book(_h["pi0Sepctra0to10_1-2x0.5-1"], 3, 1, 1);
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
      book(_h["pi0Spectra40to50_2-4x2-4"],3,1,25);*/
	//fig 18
      /*book(_h["pi0Sepctra0to10NearSide&2-4x1-2"], 4, 1, 1);
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
      book(_h["pi0Spectra40to50FarSide&4-10x2-4"],4,1,30);*/
	//fig 20
	/*book(_h["NearSide0to10_2-4x1-2"],5,1,1);
	book(_h["NearSide0to10_2x4x2-4"],5,1,2);
	book(_h["NearSide0to10_4-10x2-4"],5,1,3);
	book(_h["NearSide10to20_2-4x1-2"],5,1,4);
        book(_h["NearSide10to20_2x4x2-4"],5,1,5);
        book(_h["NearSide10to20_4-10x2-4"],5,1,6);
	book(_h["NearSide20to30_2-4x1-2"],5,1,7);
        book(_h["NearSide20to30_2x4x2-4"],5,1,8);
        book(_h["NearSide20to30_4-10x2-4"],5,1,9);
	book(_h["NearSide30to40_2-4x1-2"],5,1,10);
        book(_h["NearSide30to40_2x4x2-4"],5,1,11);
        book(_h["NearSide30to40_4-10x2-4"],5,1,12);
	book(_h["NearSide40to50_2-4x1-2"],5,1,13);
        book(_h["NearSide40to50_2x4x2-4"],5,1,14);
        book(_h["NearSide40to50_4-10x2-4"],5,1,15);
	book(_h["FarSide0to10_2-4x1-2"],5,1,16);
        book(_h["FarSide0to10_2x4x2-4"],5,1,17);
        book(_h["FarSide0to10_4-10x2-4"],5,1,18);
        book(_h["FarSide10to20_2-4x1-2"],5,1,19);
        book(_h["FarSide10to20_2x4x2-4"],5,1,20);
        book(_h["FarSide10to20_4-10x2-4"],5,1,21);
        book(_h["FarSide20to30_2-4x1-2"],5,1,22);
        book(_h["FarSide20to30_2x4x2-4"],5,1,23);
        book(_h["FarSide20to30_4-10x2-4"],5,1,24);
        book(_h["FarSide30to40_2-4x1-2"],5,1,25);
        book(_h["FarSide30to40_2x4x2-4"],5,1,26);
        book(_h["FarSide30to40_4-10x2-4"],5,1,27);
        book(_h["FarSide40to50_2-4x1-2"],5,1,28);
        book(_h["FarSide40to50_2x4x2-4"],5,1,29);
        book(_h["FarSide40to50_4-10x2-4"],5,1,30);*/

	//initialize iterators and varibles
	int minCent=0, maxCent=0, minPlane=0, maxPlane=0, a=0, e=0, iterator=1;
	float min_pT=0, max_pT=0, min_pA=0, max_pA=0, minV=0, maxV=0, b=0, c=0, d=0;
	char corrName[200], bookName[200], corrNameTrigger[200];

	//fig 6
	minCent=0, maxCent=40, minPlane=2, maxPlane=5;
	for(e=minPlane; e<=maxPlane; e++){
		for(a=minCent; a<=maxCent; a+=10){
			if(e!=5){
				snprintf(corrName,200,"CounterFig6EventPlane%iCent%iTo%i",e,a,a+10);
                                snprintf(corrNameTrigger,200,"CounterFig6EventPlane%iCent%iTo%i%s",e,a,a+10,"_Triggers");
				snprintf(bookName,200,"Fig6EventPlane%iCent%iTo%i",e,a,a+10);
			}
			else{
				snprintf(corrName,200,"CounterFig6EventPlane%i{Psi2}Cent%iTo%i",e,a,a+10);
                                snprintf(corrNameTrigger,200,"CounterFig6EventPlane%i{Psi2}Cent%iTo%i%s",e,a,a+10,"_Triggers");
				snprintf(bookName,200,"Fig6EventPlane%i{Psi2}Cent%iTo%i",e,a,a+10);
			}
			book(_h[bookName],1,1,iterator);
			book(_c[corrName], corrName);
                        book(_c[corrNameTrigger], corrNameTrigger);

			Correlator corrFig6(a,e);
                        corrFig6.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig6.SetCentrality(a,a+10);
                        corrFig6.SetEventPlaneMethod(e);
                        corrFig6.SetCorrelationFunction(_h[bookName]);
                        corrFig6.SetCounter(_c[corrName]);
                        corrFig6.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig6);
                        iterator++;
		}
	}
	iterator=1;

	//fig 12
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=0.5, max_pA=4, minV=1, maxV=1;
	b=2;
	for(a=minCent; a<=maxCent; a+=10){
		for(c=min_pA; c<=max_pA; c*=2){
			snprintf(corrName,200,"CounterFig12Cent%iTo%iPtA%2.1fTo%2.1f",a,a+10,c,c*2);
                        snprintf(corrNameTrigger,200,"CounterFig12Cent%iTo%iPtA%2.1fTo%2.1f%s",a,a+10,c,c*2, "_Triggers");
			snprintf(bookName,200,"Fig12Cent%iTo%iPtA%2.1fTo%2.1f",a,a+10,c,c*2);
			book(_h[bookName],2,1,iterator);
			book(_c[corrName], corrName);
                        book(_c[corrNameTrigger], corrNameTrigger);

			Correlator corrFig12(a,(int) c*10);
	 		corrFig12.SetCollSystemAndEnergy("AuAu200GeV");
	 		corrFig12.SetCentrality(a,a+10);
	 		corrFig12.SetTriggerRange(b,(b==4)?10:b*2);
	 		corrFig12.SetAssociatedRange(c,(c==4)?10:c*2);
	 		corrFig12.SetCorrelationFunction(_h[bookName]);
	 		corrFig12.SetCounter(_c[corrName]);
                        corrFig12.SetTriggerCounter(_c[corrNameTrigger]);
			Correlators.push_back(corrFig12);
			iterator++;
		}
	}
	iterator=1;


	//fig 15
	minCent=0, maxCent=40, min_pT=1, max_pT=2, min_pA=0.5, max_pA=2, minV=1, maxV=1;
	for(a=minCent; a<=maxCent; a+=10){
		for(b=min_pT; b<=max_pT; b*=2){
			for(c=min_pA; c<=max_pA && c!=b ; c*=2){
				snprintf(corrName,200,"CounterFig15Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
					a,a+10,b,b*2,c,c*2);
				snprintf(bookName,200,"Fig15Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
					a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
				snprintf(corrNameTrigger,200,"CounterFig15Cent%iTo%iPtA%2.1fTo%2.1f%s",a,a+10,c,c*2, "_Triggers");
				book(_h[bookName],3,1,iterator);
				book(_c[corrName], corrName);
				book(_c[corrNameTrigger], corrNameTrigger);

				Correlator corrFig15(a,(int)b*10,(int)c*10);
				corrFig15.SetCollSystemAndEnergy("AuAu200GeV");
				corrFig15.SetCentrality(a,a+10);
				corrFig15.SetTriggerRange(b,(b==4)?10:b*2);
				corrFig15.SetAssociatedRange(c,(c==4)?10:c*2);
				corrFig15.SetCorrelationFunction(_h[bookName]);
				corrFig15.SetCounter(_c[corrName]);
				corrFig15.SetTriggerCounter(_c[corrNameTrigger]);
				Correlators.push_back(corrFig15);
				iterator++;
			}
		}
	}
	iterator=1;

	//Fig 18
	minCent=0, maxCent=40, min_pT=2, max_pT=4, min_pA=1, max_pA=2, minV=1, maxV=1;
	for(int sides=0; sides<2; sides++){
		for(a=minCent; a<=maxCent; a+=10){
                	for(b=min_pT; b<=max_pT; b*=2){
                        	for((b==4) ? c=min_pA*2 : c=min_pA; c<=max_pA; c*=2){
                                	if(sides==0){
						snprintf(bookName,200,
							"NearSideFig18Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
							a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrName,200,
                                                	"CounterNearSideFig18Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
	                                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrNameTrigger,200,
							"CounterNearSideFig18Cent%iTo%iPtA%2.1fTo%2.1f%s",a,a+10,c,c*2, "_Triggers");
					}
					else{
						snprintf(bookName,200,
							"FarSideFig18Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
							a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrName,200,
	                                                "CounterFarSideFig18Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
        	                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrNameTrigger,200,
							"CounterFarSideFig12Cent%iTo%iPtA%2.1fTo%2.1f%s",a,a+10,c,c*2, "_Triggers");
					}
                                	book(_h[bookName],4,1,iterator);
					book(_c[corrName], corrName);
					book(_c[corrNameTrigger], corrNameTrigger);

	                                Correlator corrFig18(a,(int)b*10,(int)c*10);
        	                        corrFig18.SetCollSystemAndEnergy("AuAu200GeV");
                	                corrFig18.SetCentrality(a,a+10);
	                                if(b==0.5 || b==1 || b==2) corrFig18.SetTriggerRange(b,b*2);
        	                        else corrFig18.SetTriggerRange(b,10);
                	                if(c==0.5 || c==1 || c==2) corrFig18.SetAssociatedRange(c,c*2);
                        	        else corrFig18.SetAssociatedRange(c,10);
					corrFig18.SetCorrelationFunction(_h[bookName]);
	                                corrFig18.SetCounter(_c[corrName]);
					corrFig18.SetTriggerCounter(_c[corrNameTrigger]);
        	                        Correlators.push_back(corrFig18);
                	                iterator++;
				}
			}
		}
	}
	iterator=1;

	//Fig 20
	minCent=0, maxCent=40, min_pT=2, max_pT=4, min_pA=1, max_pA=2, minV=1, maxV=1;
        for(int sides=0; sides<2; sides++){
                for(a=minCent; a<=maxCent; a+=10){
                        for(b=min_pT; b<=max_pT; b*=2){
                                for((b==4) ? c=min_pA*2 : c=min_pA; c<=max_pA; c*=2){
                                	snprintf(corrName,200,
						"CounterFig20Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
                                        	a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
                                	if(sides==0){
						snprintf(bookName,200,
							"NearSideFig20Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
							a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrName,200,
	                                                "CounterNearSideFig20Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
        	                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrNameTrigger,200,
							"CounterNearSideFig20Cent%iTo%iPtA%2.1fTo%2.1f%s",a,a+10,c,c*2, "_Triggers");
					}
                                	else{
						snprintf(bookName,200,
							"FarSideFig20Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
							a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrName,200,
	                                                "CounterFarSideFig20Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
        	                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrNameTrigger,200,
							"CounterFarSideFig20Cent%iTo%iPtA%2.1fTo%2.1f%s",a,a+10,c,c*2, "_Triggers");

					}
                                        book(_h[bookName],5,1,iterator);
					book(_c[corrName], corrName);
					book(_c[corrNameTrigger], corrNameTrigger);

                                        Correlator corrFig20(a,(int)b*10,(int)c*10);
                                        corrFig20.SetCollSystemAndEnergy("AuAu200GeV");
                                        corrFig20.SetCentrality(a,a+10);
                                        if(b==0.5 || b==1 || b==2) corrFig20.SetTriggerRange(b,b*2);
                                        else corrFig20.SetTriggerRange(b,10);
                                        if(c==0.5 || c==1 || c==2) corrFig20.SetAssociatedRange(c,c*2);
                                        else corrFig20.SetAssociatedRange(c,10);
					corrFig20.SetCorrelationFunction(_h[bookName]);
                                        corrFig20.SetCounter(_c[corrName]);
					corrFig20.SetTriggerCounter(_c[corrNameTrigger]);
                                        Correlators.push_back(corrFig20);
                                        iterator++;
				}
			}
		}
	}
	iterator=1;

	//Fig 21
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=0.5, max_pA=4, minV=1, maxV=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(b=min_pT; b<=max_pT; b*=2){
                        for(c=min_pA; c<=max_pA; c*=2){
                                snprintf(corrName,200,"CounterFig21Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
				snprintf(bookName,200,
					"Fig21Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
					a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2);
				snprintf(corrNameTrigger,200,"CounterFig21Cent%iTo%iPtA%2.1fTo%2.1f%s",a,a+10,c,c*2, "_Triggers");
                                book(_h[bookName],6,1,iterator);
				book(_c[corrName], corrName);
				book(_c[corrNameTrigger], corrNameTrigger);

                                Correlator corrFig21(a,(int)c*10);
                                corrFig21.SetCollSystemAndEnergy("AuAu200GeV");
                                corrFig21.SetCentrality(a,a+10);
				corrFig21.SetTriggerRange(b,(b==4)?10:b*2);
				corrFig21.SetAssociatedRange(c,(c==4)?10:c*2);
				corrFig21.SetCorrelationFunction(_h[bookName]);
                                corrFig21.SetCounter(_c[corrName]);
				corrFig21.SetTriggerCounter(_c[corrNameTrigger]);
                                Correlators.push_back(corrFig21);
				iterator++;
			}
		}
	}
	iterator=1;

	//Fig22
	minCent=0, maxCent=40, min_pT=1, max_pT=2, min_pA=0.5, max_pA=2, minV=1, maxV=1;
	for(a=minCent; a<=maxCent; a+=10){
        	for(b=min_pT; b<=max_pT; b*=2){
        		for(c=min_pA; c<=max_pA && !(c>b); c*=2){
                                snprintf(corrName,200,
					"CounterFig22Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
                                snprintf(bookName,200,
                                        "Fig22Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
                                        a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2);
				snprintf(corrNameTrigger,200,"CounterFig22Cent%iTo%iPtA%2.1fTo%2.1f%s",a,a+10,c,c*2, "_Triggers");
                                book(_h[bookName],7,1,iterator);
				book(_c[corrName], corrName);
				book(_c[corrNameTrigger], corrNameTrigger);

                                Correlator corrFig22(a,(int)c*10);
                                corrFig22.SetCollSystemAndEnergy("AuAu200GeV");
                                corrFig22.SetCentrality(a,a+10);
                                corrFig22.SetTriggerRange(b,(b==4)?10:b*2);
                                corrFig22.SetAssociatedRange(c,(c==4)?10:c*2);
				corrFig22.SetCorrelationFunction(_h[bookName]);
                                corrFig22.SetCounter(_c[corrName]);
				corrFig22.SetTriggerCounter(_c[corrNameTrigger]);
                                Correlators.push_back(corrFig22);
                                iterator++;
                        }
                }
        }
        iterator=1;

	//Fig23
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
	for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
			snprintf(corrName,200,
				"CounterFig23Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig23Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,"CounterFig23Cent%iTo%iPtA%2.1fTo%2.1f%s",a,a+10,c,c*2, "_Triggers");
                        book(_h[bookName],8,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig23(a,(int)d);
                        corrFig23.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig23.SetCentrality(a,a+10);
                        corrFig23.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig23.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig23.SetRxnPlaneAngle(d,d+1);
			corrFig23.SetCounter(_c[corrName]);
			corrFig23.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig23);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 24
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig24Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig24Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],9,1,iterator);
                        Correlator corrFig24(a,d);
                        corrFig24.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig24.SetCentrality(a,a+10);
                        corrFig24.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig24.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig24.SetRxnPlaneAngle(d,d+1);
                        corrFig24.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig24);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 25
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig25Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig25Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],10,1,iterator);
                        Correlator corrFig25(a,d);
                        corrFig25.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig25.SetCentrality(a,a+10);
                        corrFig25.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig25.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig25.SetRxnPlaneAngle(d,d+1);
                        corrFig25.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig25);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 26
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig26Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig26Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],11,1,iterator);
                        Correlator corrFig26(a,d);
                        corrFig26.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig26.SetCentrality(a,a+10);
                        corrFig26.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig26.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig26.SetRxnPlaneAngle(d,d+1);
                        corrFig26.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig26);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 27
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig27Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig27Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],12,1,iterator);
                        Correlator corrFig27(a,d);
                        corrFig27.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig27.SetCentrality(a,a+10);
                        corrFig27.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig27.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig27.SetRxnPlaneAngle(d,d+1);
                        corrFig27.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig27);
                        iterator++;
                }
        }
        iterator=1;


	//Fig 28
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig28Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig28Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],13,1,iterator);
                        Correlator corrFig28(a,d);
                        corrFig28.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig28.SetCentrality(a,a+10);
                        corrFig28.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig28.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig28.SetRxnPlaneAngle(d,d+1);
                        corrFig28.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig28);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 29
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig29Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig29Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],14,1,iterator);
                        Correlator corrFig29(a,d);
                        corrFig29.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig29.SetCentrality(a,a+10);
                        corrFig29.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig29.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig29.SetRxnPlaneAngle(d,d+1);
                        corrFig29.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig29);
                        iterator++;
                }
        }
        iterator=1;


	//Fig 30
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig30Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig30Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],15,1,iterator);
                        Correlator corrFig30(a,d);
                        corrFig30.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig30.SetCentrality(a,a+10);
                        corrFig30.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig30.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig30.SetRxnPlaneAngle(d,d+1);
                        corrFig30.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig30);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 31
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig31Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig31Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],16,1,iterator);
                        Correlator corrFig31(a,d);
                        corrFig31.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig31.SetCentrality(a,a+10);
                        corrFig31.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig31.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig31.SetRxnPlaneAngle(d,d+1);
                        corrFig31.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig31);
                        iterator++;
                }
        }
        iterator=1;


	//Fig 32
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig32Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig32Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],17,1,iterator);
                        Correlator corrFig32(a,d);
                        corrFig32.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig32.SetCentrality(a,a+10);
                        corrFig32.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig32.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig32.SetRxnPlaneAngle(d,d+1);
                        corrFig32.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig32);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 33
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig33Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig33Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],18,1,iterator);
                        Correlator corrFig33(a,d);
                        corrFig33.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig33.SetCentrality(a,a+10);
                        corrFig33.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig33.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig33.SetRxnPlaneAngle(d,d+1);
                        corrFig33.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig33);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 34
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig34Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig34Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],19,1,iterator);
                        Correlator corrFig34(a,d);
                        corrFig34.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig34.SetCentrality(a,a+10);
                        corrFig34.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig34.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig34.SetRxnPlaneAngle(d,d+1);
                        corrFig34.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig34);
                        iterator++;
                }
        }
        iterator=1;

	//Psi2_Cor2-4x1-2GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig35Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig35Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],20,1,iterator);
                        Correlator corrFig35(a,d);
                        corrFig35.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig35.SetCentrality(a,a+10);
                        corrFig35.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig35.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig35.SetRxnPlaneAngle(d,d+1);
                        corrFig35.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig35);
                        iterator++;
                }
        }
        iterator=1;


	//Psi2_Cor2-4x2-4GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig36Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig36Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],21,1,iterator);
                        Correlator corrFig36(a,d);
                        corrFig36.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig36.SetCentrality(a,a+10);
                        corrFig36.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig36.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig36.SetRxnPlaneAngle(d,d+1);
                        corrFig36.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig36);
                        iterator++;
                }
        }
        iterator=1;

	//Psi2_Cor4-10x2-4GeV
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig37Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig37Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],22,1,iterator);
                        Correlator corrFig37(a,d);
                        corrFig37.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig37.SetCentrality(a,a+10);
                        corrFig37.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig37.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig37.SetRxnPlaneAngle(d,d+1);
                        corrFig37.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig37);
                        iterator++;
                }
        }
        iterator=1;


	//Psi2_PTY2-4x1-2GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig38Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig38Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],23,1,iterator);
                        Correlator corrFig38(a,d);
                        corrFig38.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig38.SetCentrality(a,a+10);
                        corrFig38.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig38.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig38.SetRxnPlaneAngle(d,d+1);
                        corrFig38.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig38);
                        iterator++;
                }
        }
        iterator=1;

	//Psi2_PTY2-4x2-4GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig38Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig38Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],24,1,iterator);
                        Correlator corrFig38(a,d);
                        corrFig38.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig38.SetCentrality(a,a+10);
                        corrFig38.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig38.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig38.SetRxnPlaneAngle(d,d+1);
                        corrFig38.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig38);
                        iterator++;
                }
        }
        iterator=1;


	//Psi2_PTY4-10x2-4GeV
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig39Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig39Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],25,1,iterator);
                        Correlator corrFig39(a,d);
                        corrFig39.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig39.SetCentrality(a,a+10);
                        corrFig39.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig39.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig39.SetRxnPlaneAngle(d,d+1);
                        corrFig39.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig39);
                        iterator++;
                }
        }
        iterator=1;

	//Psi3_Cor2-4x1-2GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(40)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig(40)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],26,1,iterator);
                        Correlator corrFig40(a,d);
                        corrFig40.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig40.SetCentrality(a,a+10);
                        corrFig40.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig40.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig40.SetRxnPlaneAngle(d,d+1);
                        corrFig40.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig40);
                        iterator++;
                }
        }
        iterator=1;

	//Psi3_Cor2-4x2-4GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(41)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig(41)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],27,1,iterator);
                        Correlator corrFig41(a,d);
                        corrFig41.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig41.SetCentrality(a,a+10);
                        corrFig41.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig41.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig41.SetRxnPlaneAngle(d,d+1);
                        corrFig41.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig41);
                        iterator++;
                }
        }
        iterator=1;

	//Psi3_Cor4-10x2-4GeV
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(42)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig(42)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],28,1,iterator);
                        Correlator corrFig42(a,d);
                        corrFig42.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig42.SetCentrality(a,a+10);
                        corrFig42.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig42.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig42.SetRxnPlaneAngle(d,d+1);
                        corrFig42.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig42);
                        iterator++;
                }
        }
        iterator=1;


	//Psi3_PTY2-4x1-2GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(43)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig(43)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],29,1,iterator);
                        Correlator corrFig43(a,d);
                        corrFig43.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig43.SetCentrality(a,a+10);
                        corrFig43.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig43.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig43.SetRxnPlaneAngle(d,d+1);
                        corrFig43.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig43);
                        iterator++;
                }
        }
        iterator=1;

	//Psi3_PTY2-4x2-4GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(44)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig(44)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],30,1,iterator);
                        Correlator corrFig44(a,d);
                        corrFig44.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig44.SetCentrality(a,a+10);
                        corrFig44.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig44.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig44.SetRxnPlaneAngle(d,d+1);
                        corrFig44.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig44);
                        iterator++;
                }
        }
        iterator=1;


	//Psi3_PTY4-10x2-4GeV
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(45)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,d,d+1);
                        snprintf(bookName,200,
                                "Fig(45)Cent%i-%iPtT%2.1f-%2.1fPtA%2.1f-%2.1fAngle%i/pi<phi-Psi<%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
                        book(_h[bookName],31,1,iterator);
                        Correlator corrFig45(a,d);
                        corrFig45.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig45.SetCentrality(a,a+10);
                        corrFig45.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig45.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig45.SetRxnPlaneAngle(d,d+1);
                        corrFig45.SetCounter(_c[corrName]);
                        Correlators.push_back(corrFig45);
                        iterator++;
                }
        }
        iterator=1;

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

    const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      //add calculation of reaction plane angle
      const double c = cent();
      const ParticlePair& beam = beams();
      string CollSystem = "Empty";

      //add pp collision eventually
	if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
          CollSystem = "AUAU200GeV";
          //if(fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) CollSystem += "200GeV";
      }

      if(CollSystem == "AUAU200GeV" && c > 50)
      {
        vetoEvent;
      }

      //cout << c << endl;
      for(Correlator& corr : Correlators)
      {
        if(!corr.CheckCollSystemAndEnergy(CollSystem)) continue;
        if(!corr.CheckCentrality(c)) continue;
        corr.AddWeight();
      }

      //Correlator corr = Correlators[0];

      for(auto pTrig : cfs.particles())
      {
        for (Correlator &corr : Correlators)
        {
          if (!corr.CheckCollSystemAndEnergy(CollSystem)) continue;
          if (!corr.CheckCentrality(c)) continue;
          //cout << "hi" << '\n';
          if (!corr.CheckTriggerRange(pTrig.pT() / GeV)) continue;
            ///Add a function to check reaction plane angle and make sure that the
            //if (!corr.CheckXiRange(log(pTrig.pT()/ pAssoc.pT()))) continue;
          //somthing here is stoping the pp
          //cout << c << '\n';
          //cout << "hi" << '\n';
          corr.AddTrigger();
        }
        for (auto pAssoc : cfs.particles())
        {
          for (Correlator &corr : Correlators)
          {
            if (!corr.CheckCollSystemAndEnergy(CollSystem)) continue;
            if (!corr.CheckCentrality(c)) continue;
            if (!corr.CheckTriggerRange(pTrig.pT() / GeV)) continue;
            if (!corr.CheckAssociatedRange(pAssoc.pT() / GeV)) continue;
            ///Add a function to check reaction plane angle and make sure that the
            //if (!corr.CheckXiRange(log(pTrig.pT()/ pAssoc.pT()))) continue;
            corr.AddCorrelation(pTrig, pAssoc);
          }
        }
      }



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
