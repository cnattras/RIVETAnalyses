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
#include <string>

namespace Rivet {
  class Correlator {

    private:
      std::vector<int> _indices;
      string _collSystemAndEnergy;
      pair<double,double> _centrality;
      pair<double,double> _triggerRange;
      pair<double,double> _associatedRange;
      pair<double,double> _XiRange;
      vector<int> _pid;
      bool _noCentrality = false;
      bool _noAssoc = false;
      bool _noXi = false;
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
      void SetNoXi(){ _noXi = true; }
      void SetTriggerRange(double tmin, double tmax){ _triggerRange = make_pair(tmin, tmax); }
      void SetAssociatedRange(double amin, double amax){ _associatedRange = make_pair(amin, amax); }
      void SetXiRange(double ximin, double ximax){ _XiRange = make_pair(ximin, ximax); }
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
      pair<double,double> GetXiRange(){ return _XiRange; }
      double GetXiRangeMin(){ return _XiRange.first; }
      double GetXiRangeMax(){ return _XiRange.second; }      
      vector<int> GetPID(){ return _pid; }
			double GetWeight(){ return _counter->sumW(); }
			Histo1DPtr GetCorrelationFunction(){ return _deltaPhi; }
			CounterPtr GetCounter(){ return _counter; }

			double GetDeltaPhi(Particle pAssoc, Particle pTrig)
	    {
	        // need to work on this dphi range
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
      bool CheckXiRange(double apt){ return ((apt>_XiRange.first && apt<_XiRange.second) || _noXi == true) ? true : false; }
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
        if(!CheckXiRange(apt/tpt)) return false;
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
      // fig 4 a
	    book(_h["GammaDirhPertriggerVsXiAUAU"], 3, 1, 1);
    	book(_h["GammaDirhPertriggerVsXidAU"], 3, 1, 2);

      // fig 4 b
      book(_h["IAA"], 4, 1 , 1 );
	    book(_h["IdA"], 4,  1,  2);

      // fig 5
      book(_h["IAAVsXiDirectPhoton5to7GevC"], 5, 1, 1);
	    book(_h["IAAVsXiDirectPhoton7to9GevC"], 6, 1, 1);
    	book(_h["IAAVsXiDirectPhoton9to12GevC"], 7, 1, 1);
      
      // fig 7
    	book(_h["RatiosOfIAAVsDirectPhotonPtPiOver2"], 11, 1, 1);
    	book(_h["RatiosOfIAAVsDirectPhotonPtPiOver3"], 11, 1, 2);
    	book(_h["RatiosOfIAAVsDirectPhotonPtPiOver6"], 11, 1, 3);
    	
      // fig 8
      book(_h["IAAVsXiDirectPhoton5to7GevC"], 12, 1, 1);
    	book(_h["IAAVsXiDirectPhoton7to9GevC"], 13, 1, 1);
    	book(_h["IAAVsXiDirectPhoton9to12GevC"], 14, 1, 1);
      
      //corrolators
      //fig 4 a
    	for (int i=0;i<6;i++){
        float xilower= 0+i*.4;
        float xiupper= xilower+.4;
        string corra = "sow_AUAU200_GammaDirhPertriggerVsXiAUAU" + to_string(xilower) + "to" + to_string(xiupper);
        string corrd = "sow_dAU200_GammaDirhPertriggerVsXidAU" + to_string(xilower) + "to" + to_string(xiupper);
        book(_c[corra], corra);
    	  book(_c[corrd], corrd);
        Correlator corrfi4a(1);
        corrfi4a.SetCollSystemAndEnergy("AUAU200GeV");
    	  corrfi4a.SetCentrality(0.,40.);
    	  corrfi4a.SetTriggerRange(5., 9.);
	      corrfi4a.SetAssociatedRange(.5, 7.);
        corrfi4a.SetXiRange(xilower,xiupper);
    	  corrfi4a.SetCorrelationFunction(_h[corra+"h"]);
	      corrfi4a.SetCounter(_c[corra]);
        Correlators.push_back(corrfi4a);
	      Correlator corrfi4ad(1);
        corrfi4ad.SetCollSystemAndEnergy("dAU200GeV");
    	  corrfi4ad.SetCentrality(0.,40.);
    	  corrfi4ad.SetTriggerRange(5., 9.);
	      corrfi4ad.SetAssociatedRange(.5, 7.);
        corrfi4ad.SetXiRange(xilower,xiupper);
    	  corrfi4ad.SetCorrelationFunction(_h[corrd+"h"]);
	      corrfi4ad.SetCounter(_c[corrd]);
        Correlators.push_back(corrfi4ad);
      }
      //fig 4 b
      for (int i=0;i<6;i++){
        float xilower= 0+i*.4;
        float xiupper= xilower+.4;
        string corra = "sow_AUAU200_IAA" + to_string(xilower) + "to" + to_string(xiupper);
        string corrd = "sow_dAU200_IdA" + to_string(xilower) + "to" + to_string(xiupper);
        book(_c[corra], corra);
    	  book(_c[corrd], corrd);
        Correlator corrfi4b(1);
        corrfi4b.SetCollSystemAndEnergy("AUAU200GeV");
    	  corrfi4b.SetCentrality(0.,40.);
    	  corrfi4b.SetTriggerRange(5., 9.);
	      corrfi4b.SetAssociatedRange(.5, 7.);
        corrfi4b.SetXiRange(xilower,xiupper);
    	  corrfi4b.SetCorrelationFunction(_h[corra+"h"]);
	      corrfi4b.SetCounter(_c[corra]);
        Correlators.push_back(corrfi4b);
        Correlator corrfi4bd(1);
        corrfi4bd.SetCollSystemAndEnergy("dAU200GeV");
    	  corrfi4bd.SetCentrality(0.,40.);
    	  corrfi4bd.SetTriggerRange(5., 9.);
	      corrfi4bd.SetAssociatedRange(.5, 7.);
        corrfi4bd.SetXiRange(xilower,xiupper);
    	  corrfi4bd.SetCorrelationFunction(_h[corrd+"h"]);
	      corrfi4bd.SetCounter(_c[corrd]);
        Correlators.push_back(corrfi4bd);
      }
      
      //fig 5 same as pi/2's of fig 6

      //fig 7
      for (int i=0;i<3;i++){
        int upper = 0;
        int lower = 0;

            if (i==0){
              upper = 7; 
              lower = 5;
            }
            if (i==1){
              upper = 9; 
              lower = 7;
            }
            if (i==2){
              upper = 12; 
              lower = 9;
            }
        string corrsless = "sow_AUAU200_RatiosOfIAAVsDirectPhotonPtLessThan1.2" + to_string(lower) + "to" + to_string(upper);
        string corrsmore = "sow_AUAU200_RatiosOfIAAVsDirectPhotonPtmoreThan1.2" + to_string(lower) + "to" + to_string(upper);
    	  book(_c[corrsless], corrsless);
    	  book(_c[corrsmore], corrsmore);
        Correlator corrfig7less(1);
        corrfig7less.SetCollSystemAndEnergy("AUAU200GeV");
    	  corrfig7less.SetCentrality(0.,40.);
    	  corrfig7less.SetTriggerRange(lower, upper);
	      corrfig7less.SetAssociatedRange(.5, 7.);
        corrfig7less.SetXiRange(log(5./7.),1.2);
    	  corrfig7less.SetCorrelationFunction(_h[corrsless+"h"]);
	      corrfig7less.SetCounter(_c[corrsless]);
        Correlators.push_back(corrfig7less);
        Correlator corrfig7more(1);
        corrfig7more.SetCollSystemAndEnergy("AUAU200GeV");
    	  corrfig7more.SetCentrality(0.,40.);
    	  corrfig7more.SetTriggerRange(lower, upper);
	      corrfig7more.SetAssociatedRange(.5, 7.);
        corrfig7more.SetXiRange(1.2,log(12./.5));
    	  corrfig7more.SetCorrelationFunction(_h[corrsmore+"h"]);
	      corrfig7more.SetCounter(_c[corrsmore]);
        Correlators.push_back(corrfig7more);
      }

      //fig 8 same as fig 5
      
      //histogram and corralotor loops
      //fig 2
      for (int i=0;i<6;i++){
        float tlow = 5;
        float tup = 9;
        float xilow = ((2.)-i*.4);
        float xiup = ((2.4)-i*.4);
        float aup  = exp(xiup)/tlow;
        float alow = exp(xilow)/tup;
        //string corrnum = "corr" + to_string(i+1);
        string books= "PerTriggerVsdphiAUAU" + to_string(((2.)-i*.4)) + "to" + to_string(((2.4)-i*.4));
        string corrs= "sow_AUAU200_" + books;
        book(_h[books], 1, 1, i+1);
        book(_c[corrs], corrs);
        Correlator corrfig2(i);
        corrfig2.SetCollSystemAndEnergy("AUAU200GeV");
    	  corrfig2.SetCentrality(0.,40.);
    	  corrfig2.SetTriggerRange(tlow, tup);
	      corrfig2.SetAssociatedRange(alow, aup);
        corrfig2.SetXiRange(xilow,xiup);
    	  corrfig2.SetCorrelationFunction(_h[corrs+"h"]);
	      corrfig2.SetCounter(_c[corrs]);
        Correlators.push_back(corrfig2);
      };

      //fig 3
      for (int i=0;i<5;i++){
        float tlow = 5;
        float tup = 9;
        float xilow = ((2.)-i*.4);
        float xiup = ((2.4)-i*.4);
        float aup  = exp(xiup)/tlow;
        float alow = exp(xilow)/tup;
        string corrnum = "corr" + to_string(i+1);
        string books= "PerTriggerVsdphidAU" + to_string(((2.)-i*.4)) + "to" + to_string(((2.4)-i*.4));
        string corrs= "sow_dAU200_" + books;
        book(_h[books], 2, 1, i+1);
        book(_c[corrs], corrs);
        Correlator corrfig3(i);
        corrfig3.SetCollSystemAndEnergy("dAU200GeV");
    	  corrfig3.SetCentrality(0.,40.);
    	  corrfig3.SetTriggerRange(tlow, tup);
	      corrfig3.SetAssociatedRange(alow, aup);
        corrfig3.SetXiRange(xilow,xiup);
    	  corrfig3.SetCorrelationFunction(_h[corrs+"h"]);
	      corrfig3.SetCounter(_c[corrs]);
        Correlators.push_back(corrfig3);
      };

      //fig 6 
      for (int i=0;i<3;i++){
        int upper = 0;
        int lower = 0;
        int pi = 0;
        int binnum = 0;
        float xiupper = 0;
        float xilower = 0;
        float alow = 0.5;
        float aup = 7.0;
        if (i==0){
          upper = 7;
          lower = 5;
          binnum = 5;
        }
        if (i==1){
          upper = 9;
          lower = 7;
          binnum = 6;
        }
        if (i==2){
          upper = 12;
          lower = 9;
          binnum = 6;
        }
        string forcor = "IAAVsXiDirectPhoton" + to_string(lower) + "to" + to_string(upper);
        for (int j = 0; j < 3; j++){
          if (j==0){
            pi = 2;
          }
          if (j==1){
            pi = 3;
          }
          if (j==2){
            pi = 6;
          }
          string books = forcor + "PiOver" + to_string(pi);
          book(_h[books], 8+i, 1, 1+j);
        }
        for(int k=0;k<binnum;k++){
          if (i==0){
            xilower = (-0.1) +k*.4;
            xiupper = xilower+.4;
          }
          if (i==1){
            xilower = (-0.1)+k*.4;
            xiupper = xilower+.4;
          }
          if (i==2){
            xilower = 0.3+k*.4;
            xiupper = xilower+.4;
          }
          if (xilower < 0){
            xilower = 0;
            xiupper = 0.3;
          }
          string corrs = "sow_AuAu200_" + forcor + "AndXi" + to_string(xilower) + "to" + to_string(xiupper);
          //string corfunc= forcor + to_string(xilower) + "to" + to_string(xiupper);
          book(_c[corrs], corrs);
          Correlator corrfi6(i+k+1);
          corrfi6.SetCollSystemAndEnergy("AUAU200GeV");
          corrfi6.SetCentrality(0., 40.);
          corrfi6.SetTriggerRange(lower, upper);
          corrfi6.SetAssociatedRange(alow, aup);
          corrfi6.SetXiRange(xilower, xiupper);
          corrfi6.SetCorrelationFunction(_h[corrs+"h"]);
          corrfi6.SetCounter(_c[corrs]);
          Correlators.push_back(corrfi6);
        }
      }

    };


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const CentralityProjection &cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();

      if (c > 40){
        vetoEvent;
      }
      /*
      _c["CCCC"]->fill();

      Particles fsParticles = applyProjection<FinalState>(event, "fs").particles();

      for (const Particle &p : fsParticles)
      {
        if (p.pid() == 321){
          _h["AAAA"]->fill(p.pT() / GeV);
        }
      }*/
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

