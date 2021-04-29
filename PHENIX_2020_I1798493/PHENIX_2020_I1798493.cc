// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
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
      Histo1DPtr _deltaPhi;
      CounterPtr _counter;
      CounterPtr _cTriggers;

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

        void AddCorrelation(Particle pTrig, Particle pAssoc, bool is0toPI = false)
	{
                double dPhi = GetDeltaPhi(pTrig, pAssoc);
                if(is0toPI)
                {
                        dPhi = mapAngle0ToPi(dPhi);
                        _deltaPhi->fill(dPhi, 0.5);
                }
                else
                {
                        _deltaPhi->fill(dPhi);
                }

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
        if(!CheckXiRange(log(tpt/apt))) return false;
        // should the above be log(apt/tpt) or apt/tpt
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
  void DivideScatter2D(Scatter2DPtr s1, Scatter2DPtr s2, Scatter2DPtr s)
  {
    for(unsigned int i = 0; i < s2->numPoints(); i++)
    {
      if(s2->point(i).y() == 0)
      {
        s->addPoint(s2->point(i).x(), std::numeric_limits<double>::quiet_NaN());
        continue;
      }
      double yErr = (s1->point(i).y()/s2->point(i).y())*std::sqrt(std::pow(s1->point(i).yErrPlus()/s1->point(i).y(), 2) + std::pow(s2->point(i).yErrPlus()/s2->point(i).y(), 2));
      s->addPoint(s2->point(i).x(), s1->point(i).y()/s2->point(i).y(), s1->point(i).xErrPlus(), yErr);
    }
  }
	/// @brief Add a short analysis description here
  class PHENIX_2020_I1798493 : public Analysis {
    public:

    /// Constructor
      DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2020_I1798493);

      Histo1DPtr SubtractBackgroundZYAM(Histo1DPtr histo)
      {

        YODA::Histo1D hist = *histo;

        double minValue = sqrt(-2);
        double binWidth = 0.;
        int minValueEntries = 0.;

        for (auto &bin : hist.bins())
        {
          if (std::isnan(minValue))
          {
            minValue = bin.sumW();
            binWidth = bin.width();
            minValueEntries = bin.numEntries();
          }
          if (bin.sumW() / bin.width() < minValue / binWidth)
          {
            minValue = bin.sumW();
            binWidth = bin.width();
            minValueEntries = bin.numEntries();
          }
        }
        if (minValue == 0 || minValueEntries == 0) return histo;

        hist.reset();

        for (auto &bin : hist.bins())
        {
          bin.fillBin((minValue * bin.width()) / (minValueEntries * binWidth), minValueEntries);
        }

        *histo = YODA::subtract(*histo, hist);

        return histo;
      }

      double getYieldRangeUser(Histo1DPtr histo, double xmin, double xmax, double &fraction)
      {
        //This will include bins partially covered by the user range

        YODA::Histo1D hist = *histo;

        double integral = 0.;

        if (xmax < xmin) throw RangeError("Error: xmin > xmax");
        if (xmin < hist.bin(0).xMin()) throw RangeError("xmin is out of range");
        if (xmax > hist.bin(hist.numBins() - 1).xMax()) throw RangeError("xmax is out of range");

        for (auto &bin : hist.bins())
        {
          if ((bin.xMin() > xmin) && (bin.xMax() < xmax))
          {
            integral += bin.sumW();
            fraction += bin.numEntries();
          }
          else if ((bin.xMin() < xmin) && (bin.xMax() > xmin))
          {
            double perc = bin.xMax() - xmin;
            integral += perc * bin.sumW();
            fraction += perc * bin.numEntries();
          }
          else if ((bin.xMin() < xmax) && (bin.xMax() > xmax))
          {
            double perc = xmax - bin.xMin();
            integral += perc * bin.sumW();
            fraction += perc * bin.numEntries();
          }
        }

        return integral;
      }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV);
      declare(cfs, "CFS");
	    const PromptFinalState pfs(Cuts::abseta < 0.35 && Cuts::pid == 22);
      declare(pfs, "pfs");

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      // Book histograms
      // fig 4 a
	    book(_h["GammaDirhPertriggerVsXiAUAU"], 3, 1, 1);
    	book(_h["GammaDirhPertriggerVsXidAU"], 3, 1, 2);

      //fig 4 b
      string refname = mkAxisCode(4,1,1);
      const Scatter2D& refdata = refData(refname);
      string refname1 = mkAxisCode(4,1,2);
      const Scatter2D& refdata1 = refData(refname1);
      book(_h["IAA_AuAu_4.b"], refname + "_AuAu", refdata);
      book(_h["IAA_pp_4.b"], refname + "_pp", refdata);
      book(_h["IdA_dAu_4.b"], refname1 + "_dAu", refdata1);
      book(_s["IdA_4.b"], refname1);
      book(_s["IAA_4.b"], refname);

      // fig 7
    	for (int i=0;i<3;i++){
        string pt = "0";
        if (i==0){pt="2";};
        if (i==1){pt="3";};
        if (i==2){pt="6";};
        string refname = mkAxisCode(11,1,1+i);
        const Scatter2D& refdata = refData(refname);
        string hAA = "IAA_AuAu_7.less." + pt;
        string hpp = "IAA_pp_7.less." + pt;
        string siaa = "IAA_7.less."+ pt;
        book(_h[hAA], refname + "_AuAu_less", refdata);
        book(_h[hpp], refname + "_pp_less", refdata);
        book(_s[siaa], refname + "_IAA_less");
        hAA = "IAA_AuAu_7.over." + pt;
        hpp = "IAA_pp_7.over." + pt;
        siaa = "IAA_7.over."+ pt;
        book(_h[hAA], refname + "_AuAu_over", refdata);
        book(_h[hpp], refname + "_pp_over", refdata);
        book(_s[siaa], refname + "_IAA_over");
        siaa = "IAA_7."+ pt;
        book(_s[siaa], refname);
      }

      //fig 6 and 5 and 8
      for (int i=0;i<3;i++){
        for (int k=0;k<3;k++){
          string pt = "0";
          string pii = "0";
          if (i==0){pt="a";};
          if (i==1){pt="b";};
          if (i==2){pt="c";};
          if (k==0){pii = ".2";};
          if (k==1){pii = ".3";};
          if (k==2){pii = ".6";};
          string hAA = "IAA_AuAu_6." + pt + pii;
          string hpp = "IAA_pp_6." + pt + pii;
          string siaa = "IAA_6."+ pt + pii;
          string refname = mkAxisCode(8+i,1,1+k);
          const Scatter2D& refdata = refData(refname);
          book(_h[hAA], refname + "_AuAu", refdata);
          book(_h[hpp], refname + "_pp", refdata);
          book(_s[siaa], refname);
          //fig 5 and 8
          if (k==0){
            string siaa5 = "IAA_5." + pt;
            string siaa8 = "IAA_8." + pt;
            book(_s[siaa5], 5+i, 1, 1);
            book(_s[siaa8], 12+i, 1, 1);
          };
        }
      }

      //corrolators
      int dphibinNum = 60;
      int dphibinNumCor = 20;
      //fig 4 a
      for (int i=0;i<6;i++){
        float xilower= 0+i*.4;
        float xiupper= xilower+.4;
        char buffXiLow [5];
        char buffXiUpp [5];
        snprintf(buffXiLow,5,"%2.1f", xilower);
        snprintf(buffXiUpp,5,"%2.1f", xiupper);
        string Xilow = buffXiLow;
        string Xiupp = buffXiUpp;
        string corra = "sow_AUAU200_GammaDirhPertriggerVsXiAUAU" + Xilow + "to" + Xiupp;
        string corrd = "sow_dAU200_GammaDirhPertriggerVsXidAU" + Xilow + "to" + Xiupp;
        string corra2 = "dphi_AUAU200_GammaDirhPertriggerVsXiAUAU" + Xilow + "to" + Xiupp;
        string corrd2 = "dphi_dAU200_GammaDirhPertriggerVsXidAU" + Xilow + "to" + Xiupp;
        book(_c[corra], corra);
        book(_c[corra+"_Triggers"], corra+"_Triggers");
        book(_c[corrd], corrd);
        book(_c[corrd+"_Triggers"], corrd+"_Triggers");
        book(_h[corra2], corra2, dphibinNum, -M_PI/2., 1.5*M_PI);
        book(_h[corrd2], corrd2, dphibinNum, -M_PI/2., 1.5*M_PI);
        Correlator corrfi4a(4010);
        corrfi4a.SetCollSystemAndEnergy("AUAU200GeV");
    	  corrfi4a.SetCentrality(0.,40.);
    	  corrfi4a.SetTriggerRange(5., 9.);
	      corrfi4a.SetAssociatedRange(.5, 7.);
        corrfi4a.SetXiRange(xilower,xiupper);
    	  corrfi4a.SetCorrelationFunction(_h[corra2]);
	      corrfi4a.SetCounter(_c[corra]);
        corrfi4a.SetTriggerCounter(_c[corra+"_Triggers"]);
        Correlators.push_back(corrfi4a);
	      Correlator corrfi4ad(4020);
        corrfi4ad.SetCollSystemAndEnergy("dAU200GeV");
    	  corrfi4ad.SetNoCentrality();
    	  corrfi4ad.SetTriggerRange(5., 9.);
	      corrfi4ad.SetAssociatedRange(.5, 7.);
        corrfi4ad.SetXiRange(xilower,xiupper);
    	  corrfi4ad.SetCorrelationFunction(_h[corrd2]);
	      corrfi4ad.SetCounter(_c[corrd]);
        corrfi4ad.SetTriggerCounter(_c[corrd+"_Triggers"]);
        Correlators.push_back(corrfi4ad);
      }

      //fig 4 b
      for (int i=0;i<6;i++){
        float xilower= 0+i*.4;
        float xiupper= xilower+.4;
        char buffXiLow [5];
        char buffXiUpp [5];
        snprintf(buffXiLow,5,"%2.1f", xilower);
        snprintf(buffXiUpp,5,"%2.1f", xiupper);
        string Xilow = buffXiLow;
        string Xiupp = buffXiUpp;
        string corra = "sow_AUAU200_IAA" + Xilow + "to" + Xiupp;
        string corrd = "sow_dAU200_IdA" + Xilow + "to" + Xiupp;
        string corra2 = "dphi_AUAU200_IAA" + Xilow + "to" + Xiupp;
        string corrd2 = "dphi_dAU200_IdA" + Xilow + "to" + Xiupp;
        string corrp = "sow_pp200_IdA" + Xilow + "to" + Xiupp;
        string corrp2 = "dphi_pp200_IdA" + Xilow + "to" + Xiupp;
        book(_c[corra], corra);
        book(_c[corra+"_Triggers"], corra+"_Triggers");
        book(_c[corrd], corrd);
        book(_c[corrd+"_Triggers"], corrd+"_Triggers");
        book(_c[corrp], corrp);
        book(_c[corrp+"_Triggers"], corrp+"_Triggers");
        book(_h[corra2], corra2, dphibinNum, -M_PI/2., 1.5*M_PI);
        book(_h[corrd2], corrd2, dphibinNum, -M_PI/2., 1.5*M_PI);
        book(_h[corrp2], corrp2, dphibinNum, -M_PI/2., 1.5*M_PI);
        Correlator corrfi4b(4110);
        corrfi4b.SetCollSystemAndEnergy("AUAU200GeV");
    	  corrfi4b.SetCentrality(0.,40.);
    	  corrfi4b.SetTriggerRange(5., 9.);
	      corrfi4b.SetAssociatedRange(.5, 7.);
        corrfi4b.SetXiRange(xilower,xiupper);
    	  corrfi4b.SetCorrelationFunction(_h[corra2]);
	      corrfi4b.SetCounter(_c[corra]);
        corrfi4b.SetTriggerCounter(_c[corra+"_Triggers"]);
        Correlators.push_back(corrfi4b);
        Correlator corrfi4bd(4120);
        corrfi4bd.SetCollSystemAndEnergy("dAU200GeV");
    	  corrfi4bd.SetNoCentrality();
    	  corrfi4bd.SetTriggerRange(5., 9.);
	      corrfi4bd.SetAssociatedRange(.5, 7.);
        corrfi4bd.SetXiRange(xilower,xiupper);
    	  corrfi4bd.SetCorrelationFunction(_h[corrd2]);
	      corrfi4bd.SetCounter(_c[corrd]);
        corrfi4bd.SetTriggerCounter(_c[corrd+"_Triggers"]);
        Correlators.push_back(corrfi4bd);
        Correlator corrfi4bp(4130);
        corrfi4bp.SetCollSystemAndEnergy("pp200GeV");
    	  corrfi4bp.SetNoCentrality();
    	  corrfi4bp.SetTriggerRange(5., 9.);
	      corrfi4bp.SetAssociatedRange(.5, 7.);
        corrfi4bp.SetXiRange(xilower,xiupper);
    	  corrfi4bp.SetCorrelationFunction(_h[corrp2]);
	      corrfi4bp.SetCounter(_c[corrp]);
        corrfi4bp.SetTriggerCounter(_c[corrp+"_Triggers"]);
        Correlators.push_back(corrfi4bp);
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
    	  string corrsless2 = "dphi_AUAU200_RatiosOfIAAVsDirectPhotonPtLessThan1.2" + to_string(lower) + "to" + to_string(upper);
        string corrsmore2 = "dphi_AUAU200_RatiosOfIAAVsDirectPhotonPtmoreThan1.2" + to_string(lower) + "to" + to_string(upper);
        book(_c[corrsless], corrsless);
        book(_c[corrsless+"_Triggers"], corrsless+"_Triggers");
        book(_c[corrsmore], corrsmore);
        book(_c[corrsmore+"_Triggers"], corrsmore+"_Triggers");
        book(_h[corrsless2], corrsless2, dphibinNum, -M_PI/2., 1.5*M_PI);
        book(_h[corrsmore2], corrsmore2, dphibinNum, -M_PI/2., 1.5*M_PI);
        string corrslessp = "sow_pp200_RatiosOfIAAVsDirectPhotonPtLessThan1.2" + to_string(lower) + "to" + to_string(upper);
        string corrsmorep = "sow_pp200_RatiosOfIAAVsDirectPhotonPtmoreThan1.2" + to_string(lower) + "to" + to_string(upper);
    	  string corrslessp2 = "dphi_pp200_RatiosOfIAAVsDirectPhotonPtLessThan1.2" + to_string(lower) + "to" + to_string(upper);
        string corrsmorep2 = "dphi_pp200_RatiosOfIAAVsDirectPhotonPtmoreThan1.2" + to_string(lower) + "to" + to_string(upper);
        book(_c[corrslessp], corrslessp);
        book(_c[corrslessp+"_Triggers"], corrslessp+"_Triggers");
        book(_c[corrsmorep], corrsmorep);
        book(_c[corrsmorep+"_Triggers"], corrsmorep+"_Triggers");
        book(_h[corrslessp2], corrslessp2, dphibinNum, -M_PI/2., 1.5*M_PI);
        book(_h[corrsmorep2], corrsmorep2, dphibinNum, -M_PI/2., 1.5*M_PI);
        Correlator corrfig7less(7010);
        corrfig7less.SetCollSystemAndEnergy("AUAU200GeV");
    	  corrfig7less.SetCentrality(0.,40.);
    	  corrfig7less.SetTriggerRange(lower, upper);
	      corrfig7less.SetAssociatedRange(.5, 7.);
        corrfig7less.SetXiRange(log(lower/7.),1.2);
    	  corrfig7less.SetCorrelationFunction(_h[corrsless2]);
	      corrfig7less.SetCounter(_c[corrsless]);
        corrfig7less.SetTriggerCounter(_c[corrsless+"_Triggers"]);
        Correlators.push_back(corrfig7less);
        Correlator corrfig7more(7110);
        corrfig7more.SetCollSystemAndEnergy("AUAU200GeV");
    	  corrfig7more.SetCentrality(0.,40.);
    	  corrfig7more.SetTriggerRange(lower, upper);
	      corrfig7more.SetAssociatedRange(.5, 7.);
        corrfig7more.SetXiRange(1.2,log(upper/.5));
    	  corrfig7more.SetCorrelationFunction(_h[corrsmore2]);
	      corrfig7more.SetCounter(_c[corrsmore]);
        corrfig7more.SetTriggerCounter(_c[corrsmore+"_Triggers"]);
        Correlators.push_back(corrfig7more);
        Correlator corrfig7lessp(7030);
        corrfig7lessp.SetCollSystemAndEnergy("pp200GeV");
    	  corrfig7lessp.SetNoCentrality();
    	  corrfig7lessp.SetTriggerRange(lower, upper);
	      corrfig7lessp.SetAssociatedRange(.5, 7.);
        corrfig7lessp.SetXiRange(log(5./7.),1.2);
    	  corrfig7lessp.SetCorrelationFunction(_h[corrslessp2]);
	      corrfig7lessp.SetCounter(_c[corrslessp]);
        corrfig7lessp.SetTriggerCounter(_c[corrslessp+"_Triggers"]);
        Correlators.push_back(corrfig7lessp);
        Correlator corrfig7morep(7130);
        corrfig7morep.SetCollSystemAndEnergy("pp200GeV");
    	  corrfig7morep.SetNoCentrality();
    	  corrfig7morep.SetTriggerRange(lower, upper);
	      corrfig7morep.SetAssociatedRange(.5, 7.);
        corrfig7morep.SetXiRange(1.2,log(12./.5));
    	  corrfig7morep.SetCorrelationFunction(_h[corrsmorep2]);
	      corrfig7morep.SetCounter(_c[corrsmorep]);
        corrfig7morep.SetTriggerCounter(_c[corrsmorep+"_Triggers"]);
        Correlators.push_back(corrfig7morep);
      }

      //fig 8 same as fig 5

      //histogram and corralotor loops
      //fig 2
      for (int i=0;i<6;i++){
        float tlow = 5;
        float tup = 9;
        float xilow = ((2.)-i*.4);
        float xiup = ((2.4)-i*.4);
        float aup = tup/exp(xiup);
        float alow  = tlow/exp(xilow);
        char buffXiLow [5];
        char buffXiUpp [5];
        snprintf(buffXiLow,5,"%2.1f", xilow);
        snprintf(buffXiUpp,5,"%2.1f", xiup);
        string Xilow = buffXiLow;
        string Xiupp = buffXiUpp;
        string books = "PerTriggerVsdphiAUAU" + Xilow + "to" + Xiupp;
        string corrs = "sow_AUAU200_" + books;
        string corrs2 = "dphi_AUAU200_" + books;
        book(_h[books], 1, 1, i+1);
        book(_c[corrs], corrs);
        book(_c[corrs+"_Triggers"], corrs+"_Triggers");
        //book(_h[corrs2], corrs2, dphibinNumCor, 0, M_PI);
        Correlator corrfig2(2000);
        corrfig2.SetCollSystemAndEnergy("AUAU200GeV");
    	  corrfig2.SetCentrality(0.,40.);
    	  corrfig2.SetTriggerRange(tlow, tup);
	      corrfig2.SetAssociatedRange(alow, aup);
        corrfig2.SetXiRange(xilow,xiup);
    	  corrfig2.SetCorrelationFunction(_h[books]);
	      corrfig2.SetCounter(_c[corrs]);
        corrfig2.SetTriggerCounter(_c[corrs+"_Triggers"]);
        Correlators.push_back(corrfig2);
        /*
        string corrsp = "sow_pp200_" + books;
        string corrs2p = "dphi_pp200_" + books;
        book(_c[corrsp], corrsp);
        book(_c[corrsp+"_Triggers"], corrsp+"_Triggers");
        book(_h[corrs2p], corrs2p, dphibinNum, -M_PI/2., 1.5*M_PI);
        Correlator corrfig2(i);
        corrfig2p.SetCollSystemAndEnergy("pp200GeV");
    	  corrfig2p.SetNoCentrality();
    	  corrfig2p.SetTriggerRange(tlow, tup);
	      corrfig2p.SetAssociatedRange(alow, aup);
        corrfig2p.SetXiRange(xilow,xiup);
    	  corrfig2p.SetCorrelationFunction(_h[corrs2p]);
	      corrfig2p.SetCounter(_c[corrsp]);
        corrfig2p.SetTriggerCounter(_c[corrsp+"_Triggers"]);
        Correlators.push_back(corrfig2p);
        */
      };

      //fig 3
      for (int i=0;i<5;i++){
        float tlow = 5;
        float tup = 9;
        float xilow = ((2.)-i*.4);
        float xiup = ((2.4)-i*.4);
        float aup = tup/exp(xiup);
        float alow  = tlow/exp(xilow);
        char buffXiLow [5];
        char buffXiUpp [5];
        snprintf(buffXiLow,5,"%2.1f", xilow);
        snprintf(buffXiUpp,5,"%2.1f", xiup);
        string Xilower = buffXiLow;
        string XiUpper = buffXiUpp;
        string books = "PerTriggerVsdphidAU" + Xilower + "to" + XiUpper;
        string corrs = "sow_dAU200_" + books;
        string corrs2 = "dphi_dAU200_" + books;
        book(_h[books], 2, 1, i+1);
        book(_c[corrs], corrs);
        book(_c[corrs+"_Triggers"], corrs+"_Triggers");
        //book(_h[corrs2], corrs2, dphibinNumCor, -0, M_PI);
        Correlator corrfig3(3000);
        corrfig3.SetCollSystemAndEnergy("dAU200GeV");
    	  corrfig3.SetNoCentrality();
    	  corrfig3.SetTriggerRange(tlow, tup);
	      corrfig3.SetAssociatedRange(alow, aup);
        corrfig3.SetXiRange(xilow,xiup);
    	  corrfig3.SetCorrelationFunction(_h[books]);
	      corrfig3.SetCounter(_c[corrs]);
        corrfig3.SetTriggerCounter(_c[corrs+"_Triggers"]);
        Correlators.push_back(corrfig3);
        /*
        string corrs = "sow_pp200_" + books;
        string corrs2 = "dphi_pp200_" + books;
        book(_c[corrsp], corrs);
        book(_c[corrsp+"_Triggers"], corrsp+"_Triggers");
        book(_h[corrs2p], corrs2p, dphibinNum, -M_PI/2., 1.5*M_PI);
        Correlator corrfig3p(i+1);
        corrfig3p.SetCollSystemAndEnergy("pp200GeV");
    	  corrfig3p.SetNoCentrality();
    	  corrfig3p.SetTriggerRange(tlow, tup);
	      corrfig3p.SetAssociatedRange(alow, aup);
        corrfig3p.SetXiRange(xilow,xiup);
    	  corrfig3p.SetCorrelationFunction(_h[corrs2p]);
	      corrfig3p.SetCounter(_c[corrsp]);
        corrfig3p.SetTriggerCounter(_c[corrsp+"_Triggers"]);
        Correlators.push_back(corrfig3p);
        */
      };

      //fig 6
      for (int i=0;i<3;i++){
        int upper = 0;
        int lower = 0;
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
        for(int k=0;k<binnum;k++){
          if (i==0){
            xilower = (-0.1)+k*.4;
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
          char buffXiLow5 [5];
          char buffXiUpp5 [5];
          snprintf(buffXiLow5,5,"%2.1f", xilower);
          snprintf(buffXiUpp5,5,"%2.1f", xiupper);
          string corrs = "sow_AUAU200_" + forcor + "AndXi" + buffXiLow5 + "to" + buffXiUpp5;
          string corrs2 = "dphi_AUAU200_" + forcor + "AndXi" + buffXiLow5 + "to" + buffXiUpp5;
          string corrsp = "sow_pp200_" + forcor + "AndXi" + buffXiLow5 + "to" + buffXiUpp5;
          string corrsp2 = "dphi_pp200_" + forcor + "AndXi" + buffXiLow5 + "to" + buffXiUpp5;
          book(_c[corrsp], corrsp);
          book(_c[corrsp+"_Triggers"], corrsp+"_Triggers");
          book(_h[corrsp2], corrsp2, dphibinNum, -M_PI/2., 1.5*M_PI);
          book(_c[corrs], corrs);
          book(_c[corrs+"_Triggers"], corrs+"_Triggers");
          book(_h[corrs2], corrs2, dphibinNum, -M_PI/2., 1.5*M_PI);
          Correlator corrfi6(6010+(i*100));
          corrfi6.SetCollSystemAndEnergy("AUAU200GeV");
          corrfi6.SetCentrality(0., 40.);
          corrfi6.SetTriggerRange(lower, upper);
          corrfi6.SetAssociatedRange(alow, aup);
          corrfi6.SetXiRange(xilower, xiupper);
          corrfi6.SetCorrelationFunction(_h[corrs2]);
          corrfi6.SetCounter(_c[corrs]);
          corrfi6.SetTriggerCounter(_c[corrs+"_Triggers"]);
          Correlators.push_back(corrfi6);
          Correlator corrfi6p(6030+(i*100));
          corrfi6p.SetCollSystemAndEnergy("pp200GeV");
          corrfi6p.SetNoCentrality();
          corrfi6p.SetTriggerRange(lower, upper);
          corrfi6p.SetAssociatedRange(alow, aup);
          corrfi6p.SetXiRange(xilower, xiupper);
          corrfi6p.SetCorrelationFunction(_h[corrsp2]);
          corrfi6p.SetCounter(_c[corrsp]);
          corrfi6p.SetTriggerCounter(_c[corrsp+"_Triggers"]);
          Correlators.push_back(corrfi6p);
        }
      }
    };


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const PromptFinalState& pfs = apply<PromptFinalState>(event, "pfs");
      const double c = cent();
      const ParticlePair& beam = beams();
      string CollSystem = "Empty";

      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
        CollSystem = "AUAU200GeV";
        //if(fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) CollSystem += "200GeV";
      }
      if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
      {
        CollSystem = "pp200GeV";
        //if(fuzzyEquals(sqrtS()/GeV, 200., 1E-3)) CollSystem += "200GeV";
      }
      if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000791970)
      {
        CollSystem = "dAU200GeV";
        //if(fuzzyEquals(sqrtS()/GeV, 200., 1E-3)) CollSystem += "200GeV";
      }
      if(CollSystem == "AUAU200GeV" && c > 40)
      {
        vetoEvent;
      }

      for(Correlator& corr : Correlators)
      {
	      if(!corr.CheckCollSystemAndEnergy(CollSystem)) continue;
	      if(!corr.CheckCentrality(c)) continue;
	      corr.AddWeight();
      }

      for(auto pTrig : pfs.particles())
      //for(auto pTrig : cfs.particles())
      {
        for (Correlator &corr : Correlators)
        {
          if (!corr.CheckCollSystemAndEnergy(CollSystem)) continue;
          if (!corr.CheckCentrality(c)) continue;
          if (!corr.CheckTriggerRange(pTrig.pT() / GeV)) continue;
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
            if (!corr.CheckXiRange(log(pTrig.pT()/ pAssoc.pT()))) continue;
            if(corr.GetIndex(0) == 2000 || corr.GetIndex(0) == 3000) corr.AddCorrelation(pTrig, pAssoc, true);
            else corr.AddCorrelation(pTrig, pAssoc);
          }
        }
      }
    }

    // for possable later use: mapAngle0ToPi() changing 02pi to 0pi
    /// Normalise histograms etc., after the run
    void finalize() {


      for(Correlator& corr : Correlators)
      {
        corr.Normalize();
        Histo1DPtr h = corr.GetCorrelationFunction();
        h = SubtractBackgroundZYAM(h);
        double fraction = 0.;
        double yield = getYieldRangeUser(h, M_PI/2., 3.*M_PI/2., fraction);
        double yield3 = getYieldRangeUser(h, 2.*M_PI/3., 4.*M_PI/3., fraction);
        double yield6 = getYieldRangeUser(h, 5.*M_PI/6., 7.*M_PI/6., fraction);

        //figure 4.b
        if(corr.GetIndex(0) == 4110) _h["IAA_AuAu_4.b"]->bin(_h["IAA_AuAu_4.b"]->binIndexAt( (corr.GetXiRangeMin()+corr.GetXiRangeMax())/2. )).fillBin(yield/fraction, fraction);
        if(corr.GetIndex(0) == 4120) _h["IdA_dAu_4.b"]->bin(_h["IdA_dAu_4.b"]->binIndexAt( (corr.GetXiRangeMin()+corr.GetXiRangeMax())/2. )).fillBin(yield/fraction, fraction);
        if(corr.GetIndex(0) == 4130) _h["IAA_pp_4.b"]->bin(_h["IAA_pp_4.b"]->binIndexAt( (corr.GetXiRangeMin()+corr.GetXiRangeMax())/2. )).fillBin(yield/fraction, fraction);

        //figure 6
        for (int i=0;i<3;i++){
          for (int k=0;k<3;k++){
            double y = 0;
            string pt = "0";
            string pii = "0";
            int aaid = (6010+i*100);
            int ppid = (6030+i*100);
            if (i==0){pt="a";};
            if (i==1){pt="b";};
            if (i==2){pt="c";};
            if (k==0){y=yield; pii = ".2";};
            if (k==1){y=yield3; pii = ".3";};
            if (k==2){y=yield6; pii = ".6";};
            string hAA = "IAA_AuAu_6." + pt + pii;
            string hpp = "IAA_pp_6." + pt + pii;
            if(corr.GetIndex(0) == aaid) _h[hAA]->bin(_h[hAA]->binIndexAt( (corr.GetXiRangeMin()+corr.GetXiRangeMax())/2. )).fillBin(y/fraction, fraction);
            if(corr.GetIndex(0) == ppid) _h[hpp]->bin(_h[hpp]->binIndexAt( (corr.GetXiRangeMin()+corr.GetXiRangeMax())/2. )).fillBin(y/fraction, fraction);
          }
        }

        //figure 7
        for (int i=0;i<3;i++){
          string pt = "0";
          double y = 0;
          if (i==0){pt="2"; y=yield;};
          if (i==1){pt="3"; y=yield3;};
          if (i==2){pt="6"; y=yield6;};
          //<1.2
          int aaid = (7010);
          int ppid = (7030);
          string hAA = "IAA_AuAu_7.less." + pt;
          string hpp = "IAA_pp_7.less." + pt;
          if(corr.GetIndex(0) == aaid) _h[hAA]->bin(_h[hAA]->binIndexAt( (corr.GetTriggerRangeMin()+corr.GetTriggerRangeMax())/2. )).fillBin(y/fraction, fraction);
          if(corr.GetIndex(0) == ppid) _h[hpp]->bin(_h[hpp]->binIndexAt( (corr.GetTriggerRangeMin()+corr.GetTriggerRangeMax())/2. )).fillBin(y/fraction, fraction);
          //>1.2
          aaid = (7110);
          ppid = (7130);
          hAA = "IAA_AuAu_7.over." + pt;
          hpp = "IAA_pp_7.over." + pt;
          if(corr.GetIndex(0) == aaid) _h[hAA]->bin(_h[hAA]->binIndexAt( (corr.GetTriggerRangeMin()+corr.GetTriggerRangeMax())/2. )).fillBin(y/fraction, fraction);
          if(corr.GetIndex(0) == ppid) _h[hpp]->bin(_h[hpp]->binIndexAt( (corr.GetTriggerRangeMin()+corr.GetTriggerRangeMax())/2. )).fillBin(y/fraction, fraction);
        }

      }

      //fill mod charts
      //fig 4.b
      divide(_h["IAA_AuAu_4.b"], _h["IAA_pp_4.b"], _s["IAA_4.b"]);
      divide(_h["IdA_dAu_4.b"], _h["IAA_pp_4.b"], _s["IdA_4.b"]);

      //fig 6
      for (int i=0;i<3;i++){
        for (int k=0;k<3;k++){
          string pt = "0";
          string pii = "0";
          if (i==0){pt="a";};
          if (i==1){pt="b";};
          if (i==2){pt="c";};
          if (k==0){pii = ".2";};
          if (k==1){pii = ".3";};
          if (k==2){pii = ".6";};
          string hAA = "IAA_AuAu_6." + pt + pii;
          string hpp = "IAA_pp_6." + pt + pii;
          string siaa = "IAA_6."+ pt + pii;
          divide(_h[hAA], _h[hpp], _s[siaa]);
          //fig 5 and 8
          if (k==0){
            string siaa5 = "IAA_5." + pt;
            string siaa8 = "IAA_8." + pt;
            divide(_h[hAA], _h[hpp], _s[siaa5]);
            divide(_h[hAA], _h[hpp], _s[siaa8]);
          };
        }
      }

      //fig 7
      for (int i=0;i<3;i++){
        string pt = "0";
        if (i==0){pt="2";};
        if (i==1){pt="3";};
        if (i==2){pt="6";};
        string hAA = "IAA_AuAu_7.less." + pt;
        string hpp = "IAA_pp_7.less." + pt;
        string siaal = "IAA_7.less."+ pt;
        divide(_h[hAA], _h[hpp], _s[siaal]);
        hAA = "IAA_AuAu_7.over." + pt;
        hpp = "IAA_pp_7.over." + pt;
        string siaam = "IAA_7.over."+ pt;
        divide(_h[hAA], _h[hpp], _s[siaam]);
        string siaa = "IAA_7."+ pt;
        DivideScatter2D(_s[siaal], _s[siaam], _s[siaa]);
      }

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;
    //@}

    vector<Correlator> Correlators;

    enum CollisionSystem {pp, AuAu, dAu};
    CollisionSystem collSys;
  };



  DECLARE_RIVET_PLUGIN(PHENIX_2020_I1798493);

}
