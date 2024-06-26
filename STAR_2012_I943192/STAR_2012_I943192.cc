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

//Christine was here
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
    void SetTriggerBins(vector<double> tbins)
    {
        int imax = tbins.size() - 1;
        for(int i = 0; i < imax; i++)
        {
            _triggerBins.push_back(make_pair(tbins.at(i), tbins.at(i+1)));
        }

        SetTriggerRange(_triggerBins.front().first, _triggerBins.back().second);
    }
    void SetAssiciatedBins(vector<double> abins)
    {
        int imax = abins.size() - 1;
        for(int i = 0; i < imax; i++)
        {
            _associatedBins.push_back(make_pair(abins.at(i), abins.at(i+1)));
        }

        SetAssociatedRange(_associatedBins.front().first, _associatedBins.back().second);
    }

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
    vector<pair<double,double>> GetTriggerBins(){ return _triggerBins; }
    vector<pair<double,double>> GetAssociatedBins(){ return _associatedBins; }

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
    bool CheckTriggerBin(double tpt){ return (tpt>_triggerBins.front().first && tpt<_triggerBins.back().second) ? true : false; }
    int GetTriggerBinIndex(double tpt)
    {
        if(!CheckTriggerBin(tpt)) return -1;
        for(unsigned int i = 0; i < _triggerBins.size(); i++)
        {
            if(tpt < _triggerBins.at(i).second) return i;
        }

        return -1;
    }
    bool CheckAssociatedBin(double apt){ return (apt>_associatedBins.front().first && apt<_associatedBins.back().second) ? true : false; }
    int GetAssociatedBinIndex(double apt)
    {
        if(!CheckAssociatedBin(apt)) return -1;
        for(unsigned int i = 0; i < _associatedBins.size(); i++)
        {
            if(apt < _associatedBins.at(i).second) return i;
        }

        return -1;
    }





    int _index;
    string _collSystemAndEnergy;
    pair<double,double> _centrality;
    pair<double,double> _triggerRange;
    pair<double,double> _associatedRange;
    vector<pair<double,double>> _triggerBins;
    vector<pair<double,double>> _associatedBins;


  };

  class STAR_2012_I943192 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2012_I943192);


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
        double minVal = 1.e99;

        for(unsigned int i = 0; i < hist.numBins(); i++)
        {
            if(hist.bin(i).numEntries() == 0) continue;
            if(hist.bin(i).xMin() < bmin || hist.bin(i).xMax() > bmax) continue;
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

        if(minBin < 0) return;

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
        double etaCorrFactor = 1.;

        double maxDeltaEta = 2.;

        int binEta = hist.binIndexAt(deltaEta);

        if(binEta < 0 && abs(deltaEta) < 2)
        {
            etaCorrFactor = 1. - (1./maxDeltaEta)*abs(deltaEta);
        }
        else
        {
            double binCenterEta = hist.bin(binEta).xMid();
            etaCorrFactor = 1. - (1./maxDeltaEta)*abs(binCenterEta);
        }

        return etaCorrFactor;
    }

    //bmin and bmax are included in the integral. Range = [bmin, bmax]
    //Give the min and max bins to calculate the integral
    double GetYieldInBinRange(YODA::Histo1D& hist, int bmin, int bmax, double &n)
    {
        double integral = 0.;
        double entries = 0.;

        if(bmin < 0 || bmax > (int)hist.numBins())
        {
            MSG_ERROR("Out of range!");
            return 0.;
        }

        for(int i = bmin; i <= bmax; i++)
        {
            integral += hist.bin(i).sumW();
            entries += hist.bin(i).numEntries();
        }

        n = entries;

        return integral;

    }

    //vmin is included in the integral, but vmax is not. Range = [vmin, vmax[
    //Give the min and max values to calculate the integral
    double GetYieldInUserRange(YODA::Histo1D& hist, double vmin, double vmax, double &n)
    {
        double integral = 0.;
        double entries = 0.;

        if(vmin < hist.bin(0).xMin() || vmax > hist.bin((int)hist.numBins()-1).xMax())
        {
            MSG_ERROR("Out of range!");
            return 0.;
        }

        int bmin = hist.binIndexAt(vmin);
        int bmax = hist.binIndexAt(vmax);
        if(bmax < 0) bmax = (int)hist.numBins()-1;

        for(int i = bmin; i <= bmax; i++)
        {
            integral += hist.bin(i).sumW();
            entries += hist.bin(i).numEntries();
        }

        n = entries;

        return integral;

    }

    /// Book histograms and initialise projections before the run
    void init() {

      // the basic final-state projection: all final-state particles within the given eta acceptance

      const ChargedFinalState cfs(Cuts::abseta < 1.0 && Cuts::pT > 1*GeV && Cuts::abscharge > 0);
      declare(cfs, "CFS");
      const ChargedFinalState cfsTrig(Cuts::abseta < 1.0 && Cuts::pT > 2*GeV && Cuts::abscharge > 0);
      declare(cfsTrig, "CFSTrig");

      // Declare centrality projection
      //LATER FIX TO USE STAR
      //declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");
      declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

      //==================================================
      // Create one correlator for each set of Collisions System / Beam Energy / Centrality Interval / Trigger pT interval / Associated pT interval
      // The number used in the constructor is the y-axis index of the histograms (d00-x00-y00) of the paper
      // Ex.: Correlator c1(1); -> is the correlator for histograms _h["0311"], _h["0411"], etc
      // Ex.: Correlator c2(2); -> is the correlator for histograms _h["0312"], _h["0412"], etc
       //==================================================

      Correlator c1(1);
      c1.SetCollSystemAndEnergy("CuCu62GeV");
      c1.SetCentrality(0., 60.);
      //c1.SetTriggerRange(3., 6.);
      //c1.SetAssociatedRange(1.5, 999.);
      c1.SetAssiciatedBins(_assocBins);
      c1.SetTriggerBins(_trigBins);
      Correlators.push_back(c1);

      Correlator c2(2);
      c2.SetCollSystemAndEnergy("AuAu62GeV");
      c2.SetCentrality(0., 80.);
      //c2.SetTriggerRange(3., 6.);
      //c2.SetAssociatedRange(1.5, 999.);
      c2.SetAssiciatedBins(_assocBins);
      c2.SetTriggerBins(_trigBins);
      Correlators.push_back(c2);

      Correlator c3(3);
      c3.SetCollSystemAndEnergy("dAu200GeV");
      c3.SetCentrality(0., 95.);
      //c3.SetTriggerRange(3., 6.);
      //c3.SetAssociatedRange(1.5, 999.);
      c3.SetAssiciatedBins(_assocBins);
      c3.SetTriggerBins(_trigBins);
      Correlators.push_back(c3);

      Correlator c4(4);
      c4.SetCollSystemAndEnergy("CuCu200GeV");
      c4.SetCentrality(0., 60.);
      //c4.SetTriggerRange(3., 6.);
      //c4.SetAssociatedRange(1.5, 999.);
      c4.SetAssiciatedBins(_assocBins);
      c4.SetTriggerBins(_trigBins);
      Correlators.push_back(c4);

      Correlator c5(5);
      c5.SetCollSystemAndEnergy("AuAu200GeV");
      c5.SetCentrality(40., 80.);
      c5.SetTriggerRange(3., 6.);
      c5.SetAssociatedRange(1.5, 999.);
      c5.SetAssiciatedBins(_assocBins);
      c5.SetTriggerBins(_trigBins);
      Correlators.push_back(c5);

      Correlator c6(6);
      c6.SetCollSystemAndEnergy("AuAu200GeV");
      c6.SetCentrality(0., 12.);
      c6.SetTriggerRange(3., 6.);
      c6.SetAssociatedRange(1.5, 999.);
      c6.SetAssiciatedBins(_assocBins);
      c6.SetTriggerBins(_trigBins);
      Correlators.push_back(c6);



      /*
      book(_h["0311"], 3, 1, 1);
      book(_h["0312"], 3, 1, 2);
      book(_h["0313"], 3, 1, 3);
      book(_h["0314"], 3, 1, 4);
      book(_h["0315"], 3, 1, 5);
      book(_h["0316"], 3, 1, 6);
      //sample correlation functions in delta phi before background subtraction (Fig. 6)
      book(_h["0411"], 4, 1, 1);
      book(_h["0412"], 4, 1, 2);
      book(_h["0413"], 4, 1, 3);
      book(_h["0414"], 4, 1, 4);
      book(_h["0415"], 4, 1, 5);
      book(_h["0416"], 4, 1, 6);
      //sample correlation functions in delta eta after background subtraction	  (Fig. 7)
      book(_h["0511"], 5, 1, 1);
      book(_h["0512"], 5, 1, 2);
      book(_h["0513"], 5, 1, 3);
      book(_h["0514"], 5, 1, 4);
      book(_h["0515"], 5, 1, 5);
      book(_h["0516"], 5, 1, 6);
      //sample correlation functions in delta phi before background subtraction (Fig. 7)
      book(_h["0611"], 6, 1, 1);
      book(_h["0612"], 6, 1, 2);
      book(_h["0613"], 6, 1, 3);
      book(_h["0614"], 6, 1, 4);
      book(_h["0615"], 6, 1, 5);
      book(_h["0616"], 6, 1, 6);
      //yield vs Npart (Fig. 8)
      book(_h["0711"], 7, 1, 1);
      book(_h["0811"], 8, 1, 1);
      book(_h["0911"], 9, 1, 1);
      book(_h["1011"], 10, 1, 1);
      book(_h["1111"], 11, 1, 1);
      //yield vs pTtrig (Fig. 9)
      book(_h["1211"], 12, 1, 1);
      book(_h["1212"], 12, 1, 2);
      book(_h["1213"], 12, 1, 3);
      book(_h["1214"], 12, 1, 4);
      book(_h["1215"], 12, 1, 5);
      book(_h["1216"], 12, 1, 6);
      book(_h["1217"], 12, 1, 7);
      book(_h["1218"], 12, 1, 8);
      //Yield vs pTassoc (Fig. 10).  Commented out histograms are PYTHIA
      //book(_h["1311"], 13, 1, 1);
      book(_h["1312"], 13, 1, 2);
      book(_h["1313"], 13, 1, 3);
      //book(_h["1314"], 13, 1, 4);
      book(_h["1315"], 13, 1, 5);
      book(_h["1316"], 13, 1, 6);
      book(_h["1317"], 13, 1, 7);
      book(_h["1318"], 13, 1, 8);
      //Delta Phis Widths vs pTtrigger (Fig. 11a).  PYTHIA commended out.
      //book(_h["1411"], 14, 1, 1);
      book(_h["1412"], 14, 1, 2);
      book(_h["1413"], 14, 1, 3);
      //book(_h["1414"], 14, 1, 4);
      book(_h["1415"], 14, 1, 5);
      book(_h["1416"], 14, 1, 6);
      book(_h["1417"], 14, 1, 7);
      book(_h["1418"], 14, 1, 8);
      //Delta Phis Widths vs pTassoc (Fig. 11b).  PYTHIA commended out.
      //book(_h["1511"], 15, 1, 1);
      book(_h["1512"], 15, 1, 2);
      book(_h["1513"], 15, 1, 3);
      //book(_h["1514"], 15, 1, 4);
      book(_h["1515"], 15, 1, 5);
      book(_h["1516"], 15, 1, 6);
      book(_h["1517"], 15, 1, 7);
      book(_h["1518"], 15, 1, 8);
      //Delta Phis Widths vs Npart (Fig. 11c).
      book(_h["1611"], 16, 1, 1);
      book(_h["1711"], 17, 1, 1);
      book(_h["1811"], 18, 1, 1);
      book(_h["1911"], 19, 1, 1);
      book(_h["2011"], 20, 1, 1);
      //Delta Eta Widths vs pTtrig (Fig. 11d).  PYTHIA commended out.
      //book(_h["2111"], 21, 1, 1);
      book(_h["2112"], 21, 1, 2);
      book(_h["2113"], 21, 1, 3);
      //book(_h["2114"], 21, 1, 4);
      book(_h["2115"], 21, 1, 5);
      book(_h["2116"], 21, 1, 6);
      book(_h["2117"], 21, 1, 7);
      book(_h["2118"], 21, 1, 8);
      //Delta Eta Widths vs pTassoc (Fig. 11e).  PYTHIA commended out.
      //book(_h["2211"], 22, 1, 1);
      book(_h["2212"], 22, 1, 2);
      book(_h["2213"], 22, 1, 3);
      //book(_h["2214"], 22, 1, 4);
      book(_h["2215"], 22, 1, 5);
      book(_h["2216"], 22, 1, 6);
      book(_h["2217"], 22, 1, 7);
      book(_h["2218"], 22, 1, 8);
      //Delta Eta Widths vs Npart (Fig. 11f)
      book(_h["2311"], 23, 1, 1);
      book(_h["2411"], 24, 1, 1);
      book(_h["2511"], 25, 1, 1);
      book(_h["2611"], 26, 1, 1);
      book(_h["2711"], 27, 1, 1);
      //Figure 12 - ridge yields, obsolete, not implementing!
      //book(_h["2811"], 28, 1, 1);
      //book(_h["2911"], 29, 1, 1);
      //book(_h["3011"], 30, 1, 1);
      //book(_h["3111"], 31, 1, 1);
      //Figure 13 - ridge/jet yields, obsolete, not implementing!
      //book(_h["3211"], 32, 1, 1);
      //book(_h["3311"], 33, 1, 1);
      //book(_h["3411"], 34, 1, 1);
      //book(_h["3511"], 35, 1, 1);
      //Fig. 14 v3^2/v2^2, mostly obsolete but not easy to implement anyways - not implementing!
      //book(_h["3611"], 36, 1, 1);
      //book(_h["3711"], 37, 1, 1);
      //book(_h["3811"], 38, 1, 1);
      //book(_h["3911"], 39, 1, 1);
      //book(_h["4011"], 40, 1, 1);
      //Add declaration of histograms for correlation functions with all correlation functions
      */



      for(Correlator& corr : Correlators)
      {
          int index = corr.GetIndex();

          book(sow[index],"sow" + to_string(index));
          book(_h["031" + to_string(index)], 3, 1, index);
          book(_h["041" + to_string(index)], 4, 1, index);
          book(_h["061" + to_string(index)], 6, 1, index);
          book(_DeltaPhi[index], "DeltaPhi" + to_string(index), 24, 0, M_PI);
          book(_DeltaPhiSub[index], "DeltaPhiSub" + to_string(index), 24, 0, M_PI);
          book(_DeltaEta[index], "DeltaEta" + to_string(index), 20, 0, 2);



          for(unsigned int itr = 0; itr < corr.GetTriggerBins().size(); itr++)
          {
              book(_DeltaEtaForYieldsTriggerBins[index][itr], "DeltaEtaForYields" + to_string(index) + "_TriggerBin" + to_string(itr), 20, -2, 0);
              book(_DeltaEtaForYieldsTriggerBinsCorr[index][itr], "DeltaEtaForYieldsCorr" + to_string(index) + "_TriggerBin" + to_string(itr), 20, -2, 0);
          }

          for(unsigned int ias = 0; ias < corr.GetAssociatedBins().size(); ias++)
          {
              book(_DeltaEtaForYieldsAssociatedBins[index][ias], "DeltaEtaForYields" + to_string(index) + "_AssociatedBin" + to_string(ias), 20, -2, 0);
              book(_DeltaEtaForYieldsAssociatedBinsCorr[index][ias], "DeltaEtaForYieldsCorr" + to_string(index) + "_AssociatedBin" + to_string(ias), 20, -2, 0);
          }

      }

      //Yields. Declared here because indeces of correlators do not match hepdata
      book(_YieldsDeltaEtaTriggerBins[c1.GetIndex()], 12, 1, 2);
      book(_YieldsDeltaEtaTriggerBins[c2.GetIndex()], 12, 1, 3);
      book(_YieldsDeltaEtaTriggerBins[c3.GetIndex()], 12, 1, 5);
      book(_YieldsDeltaEtaTriggerBins[c4.GetIndex()], 12, 1, 6);
      book(_YieldsDeltaEtaTriggerBins[c5.GetIndex()], 12, 1, 7);
      book(_YieldsDeltaEtaTriggerBins[c6.GetIndex()], 12, 1, 8);

      nEvents.assign(Correlators.size()+1, 0);
      nTriggers.assign(Correlators.size()+1, 0);

      vector<int> v;
      v.assign(7, 0);
      nTriggersPerTriggerBin.assign(Correlators.size()+1, v);


    }


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
                nTriggersPerTriggerBin[corr.GetIndex()][corr.GetTriggerBinIndex(pTrig.pt()/GeV)]++;
            }

		    for(const Particle& pAssoc : cfs.particles()) {

                if(pAssoc.pt()/GeV < associatedptMin || pAssoc.pt()/GeV > associatedptMax) continue;

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

                    double etaCorrection = 1.;

                    //if(abs(dEta) < 1.78) etaCorrection = EtaEffCorrection(-abs(dEta), *_DeltaEta[corr.GetIndex()]);

                    //if(pTrig.pt()/GeV > 3. && pAssoc.pt()/GeV > 1.5)
                    //{
                    if(abs(dPhi) < 0.78)
                    {
                        if(pAssoc.pt()/GeV > 1.5 && pAssoc.pt()/GeV < pTrig.pt()/GeV)
                        {
                            _DeltaEtaForYieldsTriggerBins[corr.GetIndex()][corr.GetTriggerBinIndex(pTrig.pt()/GeV)]->fill(-abs(dEta), 0.5/etaCorrection);
                            _DeltaEtaForYieldsTriggerBinsCorr[corr.GetIndex()][corr.GetTriggerBinIndex(pTrig.pt()/GeV)]->fill(-abs(dEta), 0.5/etaCorrection);
                            if(pTrig.pt()/GeV > 3.)
                            {
                                _h["031" + to_string(corr.GetIndex())]->fill(-abs(dEta), 0.5/etaCorrection);
                            }
                        }

                        if(pTrig.pt()/GeV > 3.)
                        {
                            _DeltaEtaForYieldsAssociatedBins[corr.GetIndex()][corr.GetAssociatedBinIndex(pAssoc.pt()/GeV)]->fill(-abs(dEta), 0.5/etaCorrection);
                            _DeltaEtaForYieldsAssociatedBinsCorr[corr.GetIndex()][corr.GetAssociatedBinIndex(pAssoc.pt()/GeV)]->fill(-abs(dEta), 0.5/etaCorrection);
                        }

                    }

                    if(abs(dEta) < 1.78)
                    {
                        if(pTrig.pt()/GeV > 3. && pAssoc.pt()/GeV > 1.5 && pAssoc.pt()/GeV < pTrig.pt()/GeV)
                        {
                            _h["041" + to_string(corr.GetIndex())]->fill(-abs(dPhi), 0.5);
                            _h["061" + to_string(corr.GetIndex())]->fill(-abs(dPhi), 0.5/etaCorrection);
                            _DeltaPhi[corr.GetIndex()]->fill(abs(dPhi), 0.5/etaCorrection);
                            _DeltaPhiSub[corr.GetIndex()]->fill(abs(dPhi), 0.5/etaCorrection);
                        }

                    }

                    //}





                } //end of correlators loop

			  } // associated hadrons
		    } // end of loop over associated particles
		  } // trigger hadrons
      } // particle loop

    }


    /// Normalise histograms etc., after the run
    void finalize() {



        //for(unsigned int i = 1; i <= Correlators.size(); i++)
        for(Correlator& corr : Correlators)
        {
            int index = corr.GetIndex();

            if(nTriggers[index] > 0)
            {
                _h["031" + to_string(index)]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));
                _h["041" + to_string(index)]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));
                _h["061" + to_string(index)]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));

                _DeltaPhi[index]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));
                _DeltaPhiSub[index]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));

                vector<int> n{2,3};
                SubtractBackground(*_DeltaPhi[index], *_h["061" + to_string(index)], n, 0.63, 2.51);
                SubtractBackground(*_DeltaPhi[index], *_DeltaPhiSub[index], n, 0.63, 2.51);

                for(unsigned int itr = 0; itr < corr.GetTriggerBins().size(); itr++)
                {
                    //TO DO: Correct number of triggers has to be checked
                    //_DeltaEtaForYieldsTriggerBins[index][itr]->scaleW((double)nEvents[index]/(nTriggersPerTriggerBin[index][itr]*sow[index]->sumW()));
                    //_DeltaEtaForYieldsTriggerBinsCorr[index][itr]->scaleW((double)nEvents[index]/(nTriggersPerTriggerBin[index][itr]*sow[index]->sumW()));

                    _DeltaEtaForYieldsTriggerBins[index][itr]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));
                    _DeltaEtaForYieldsTriggerBinsCorr[index][itr]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));

                    SubtractBackground(*_DeltaEtaForYieldsTriggerBins[index][itr], *_DeltaEtaForYieldsTriggerBinsCorr[index][itr], n, -2, -0.78);


                    double entries = 0.;
                    double yields = GetYieldInUserRange(*_DeltaEtaForYieldsTriggerBinsCorr[index][itr],-0.78,0., entries);


                    (*_YieldsDeltaEtaTriggerBins[index]).bin(itr).fillBin(yields*2);



                }

                for(unsigned int ias = 0; ias < corr.GetAssociatedBins().size(); ias++)
                {
                    _DeltaEtaForYieldsAssociatedBins[index][ias]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));
                    _DeltaEtaForYieldsAssociatedBinsCorr[index][ias]->scaleW((double)nEvents[index]/(nTriggers[index]*sow[index]->sumW()));

                    SubtractBackground(*_DeltaEtaForYieldsAssociatedBins[index][ias], *_DeltaEtaForYieldsAssociatedBinsCorr[index][ias], n, -2, -0.78);
                }


            }

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
    map<int, Histo1DPtr> _DeltaEta;
    map<int, Histo1DPtr> _DeltaPhiSub;

    map<int, map<int, Histo1DPtr>> _DeltaEtaForYieldsTriggerBins;
    map<int, map<int, Histo1DPtr>> _DeltaEtaForYieldsAssociatedBins;

    map<int, map<int, Histo1DPtr>> _DeltaEtaForYieldsTriggerBinsCorr;
    map<int, map<int, Histo1DPtr>> _DeltaEtaForYieldsAssociatedBinsCorr;

    map<int, Histo1DPtr> _YieldsDeltaEtaTriggerBins;

    bool fillTrigger = true;
    vector<int> nTriggers;
    vector<int> nEvents;
    vector<vector<int>> nTriggersPerTriggerBin;
    vector<Correlator> Correlators;
    vector<double> _trigBins{2., 2.5, 3., 3.5, 4., 4.5, 5., 6.};
    vector<double> _assocBins{1., 1.5, 2., 3.};

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(STAR_2012_I943192);


}
