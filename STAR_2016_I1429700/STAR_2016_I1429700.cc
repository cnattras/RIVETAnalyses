
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

#include "../Centralities/RHICCentrality.hh"

#define _USE_MATH_DEFINES

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
    void SetDeltaEtaRange(double emin, double emax) { _etarange = make_pair(emin, emax); }
    
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
    pair<double,double> _etarange;

  
  }; 




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
	if (( p.hasAncestorWith(Cuts::pid == 310) || p.hasAncestorWith(Cuts::pid == -310)  ||     // K0s
		   p.hasAncestorWith(Cuts::pid == 130)  || p.hasAncestorWith(Cuts::pid == -130)  ||     // K0l
	           p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) ||     // Xi0
	           p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) ||     // Lambda
	           p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) ||     // Sigma+/-
	           p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) ||     // Xi-/+
	           p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334) ))    // Omega-/+
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
        
        int binEta = hist.indexAt(deltaEta);
        
        if(binEta < 0)
        {
            return 0.;
        }
        
        double binCenterEta = hist.bin(binEta).xMid();
        
        return 1. - (1./maxDeltaEta)*abs(binCenterEta);
    }

 class STAR_2016_I1429700 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(STAR_2016_I1429700);

/// Book histograms and initialise projections before the ruvoid init() {

      // the basic final-state projection: all final-state particles within the given eta acceptance
     
     void init() {
      const ChargedFinalState cfs(Cuts::pT > 1*GeV); //Not cutting in eta, so no need to correct for pair acceptance
      declare(cfs, "CFS");
      const ChargedFinalState cfsTrig(Cuts::abseta < 1.0 && Cuts::pT > 2*GeV);
      declare(cfsTrig, "CFSTrig");
      

      // Declare centrality projection
      //LATER FIX TO USE STAR
      //declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");
      declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");
      



/*    //Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // the basic final-state projection: all final-state particles within the given eta acceptance

      const ChargedFinalState cfs(Cuts::abseta < 5 && Cuts::pT > 100*MeV);
      declare(cfs, "CFS");
      const ChargedFinalState cfsTrig(Cuts::abseta < 5 && Cuts::pT > 2*GeV);
      declare(cfsTrig, "CFSTrig");

      // Declare centrality projection
      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");
*/


//==================================================
// Create one correlator for each set of Collisions System / Beam Energy / Centrality Interval / Trigger pT interval / Associated pT interval
// The number used in the constructor is the y-axis index of the histograms (d00-x00-y00) of the paper
// Ex.: Correlator c1(1); -> is the correlator for histograms _h["0311"], _h["0411"], etc
// Ex.: Correlator c2(2); -> is the correlator for histograms _h["0312"], _h["0412"], etc
//==================================================

      //Figure 3
      //Cu+Cu baryons/mesons ratio
      Correlator c1(1);
      c1.SetCollSystemAndEnergy("CuCu200GeV");
      c1.SetCentrality(0., 20.);
      c1.SetTriggerRange(3., 6.);
      c1.SetAssociatedRange(2., 3.);
      Correlators.push_back(c1);

//Figure 2
//d+Au Charged baryons
      Correlator c2(2);
      c2.SetCollSystemAndEnergy("dAu200GeV");
      c2.SetCentrality(40., 80.);
      c2.SetTriggerRange(3., 6.);
      c2.SetAssociatedRange(1.5, 999.);
      c2.SetDeltaEtaRange(-1.5, -1.);
      Correlators.push_back(c2);

      Correlator c3(1);
      c3.SetCollSystemAndEnergy("dAu200GeV");
      c3.SetCentrality(40., 80.);
      c3.SetTriggerRange(3., 6.);
      c3.SetAssociatedRange(1.5, 999.);
      c3.SetDeltaEtaRange(-1., -0.5);
      Correlators.push_back(c3);

      Correlator c4(1);
      c4.SetCollSystemAndEnergy("dAu200GeV");
      c4.SetCentrality(40., 80.);
      c4.SetTriggerRange(3., 6.);
      c4.SetAssociatedRange(1.5, 999.);
      c4.SetDeltaEtaRange(-0.5, 0.);
      Correlators.push_back(c4);

      Correlator c5(1);
      c5.SetCollSystemAndEnergy("dAu200GeV");
      c5.SetCentrality(40., 80.);
      c5.SetTriggerRange(3., 6.);
      c5.SetAssociatedRange(1.5, 999.);
      c5.SetDeltaEtaRange(0., 0.5);
      Correlators.push_back(c5);

      Correlator c6(1);
      c6.SetCollSystemAndEnergy("dAu200GeV");
      c6.SetCentrality(40., 80.);
      c6.SetTriggerRange(3., 6.);
      c6.SetAssociatedRange(1.5, 999.);
      c6.SetDeltaEtaRange(0.5, 1.);
      Correlators.push_back(c6);

      Correlator c7(1);
      c7.SetCollSystemAndEnergy("dAu200GeV");
      c7.SetCentrality(40., 80.);
      c7.SetTriggerRange(3., 6.);
      c7.SetAssociatedRange(1.5, 999.);
      c7.SetDeltaEtaRange(1., 1.5);
      Correlators.push_back(c7);

      //Figure 2
      //Cu+Cu charged baryons
      Correlator c8(3);
      c8.SetCollSystemAndEnergy("CuCu200GeV");
      c8.SetCentrality(0., 20.);
      c8.SetTriggerRange(3., 6.);
      c8.SetAssociatedRange(1.5, 999.);
      c8.SetDeltaEtaRange(-1.5, -1.);
      Correlators.push_back(c8);

      Correlator c9(2);
      c9.SetCollSystemAndEnergy("CuCu200GeV");
      c9.SetCentrality(0., 20.);
      c9.SetTriggerRange(3., 6.);
      c9.SetAssociatedRange(1.5, 999.);
      c9.SetDeltaEtaRange(-1., -0.5);
      Correlators.push_back(c9);

      Correlator c10(2);
      c10.SetCollSystemAndEnergy("CuCu200GeV");
      c10.SetCentrality(0., 20.);
      c10.SetTriggerRange(3., 6.);
      c10.SetAssociatedRange(1.5, 999.);
      c10.SetDeltaEtaRange(-0.5, 0.);
      Correlators.push_back(c10);

      Correlator c11(2);
      c11.SetCollSystemAndEnergy("CuCu200GeV");
      c11.SetCentrality(0., 20.);
      c11.SetTriggerRange(3., 6.);
      c11.SetAssociatedRange(1.5, 999.);
      c11.SetDeltaEtaRange(0., 0.5);
      Correlators.push_back(c11);

      Correlator c12(2);
      c12.SetCollSystemAndEnergy("CuCu200GeV");
      c12.SetCentrality(0., 20.);
      c12.SetTriggerRange(3., 6.);
      c12.SetAssociatedRange(1.5, 999.);
      c12.SetDeltaEtaRange(0.5, 1.);
      Correlators.push_back(c12);

      Correlator c13(2);
      c13.SetCollSystemAndEnergy("CuCu200GeV");
      c13.SetCentrality(0., 20.);
      c13.SetTriggerRange(3., 6.);
      c13.SetAssociatedRange(1.5, 999.);
      c13.SetDeltaEtaRange(1., 1.5);
      Correlators.push_back(c13);

      //Figure 2
      //Au+Au charged baryons
      Correlator c14(3);
      c14.SetCollSystemAndEnergy("AuAu200GeV");
      c14.SetCentrality(40., 80.);
      c14.SetTriggerRange(3., 6.);
      c14.SetAssociatedRange(1.5, 999.);
      c14.SetDeltaEtaRange(-1.5, -1.);
      Correlators.push_back(c14);

      Correlator c15(3);
      c15.SetCollSystemAndEnergy("AuAu200GeV");
      c15.SetCentrality(40., 80.);
      c15.SetTriggerRange(3., 6.);
      c15.SetAssociatedRange(1.5, 999.);
      c15.SetDeltaEtaRange(-1., -0.5);
      Correlators.push_back(c15);

      Correlator c16(3);
      c16.SetCollSystemAndEnergy("AuAu200GeV");
      c16.SetCentrality(40., 80.);
      c16.SetTriggerRange(3., 6.);
      c16.SetAssociatedRange(1.5, 999.);
      c16.SetDeltaEtaRange(-0.5, 0.);
      Correlators.push_back(c16);

      Correlator c17(3);
      c17.SetCollSystemAndEnergy("AuAu200GeV");
      c17.SetCentrality(40., 80.);
      c17.SetTriggerRange(3., 6.);
      c17.SetAssociatedRange(1.5, 999.);
      c17.SetDeltaEtaRange(0., 0.5);
      Correlators.push_back(c17);

      Correlator c18(3);
      c18.SetCollSystemAndEnergy("AuAu200GeV");
      c18.SetCentrality(40., 80.);
      c18.SetTriggerRange(3., 6.);
      c18.SetAssociatedRange(1.5, 999.);
      c18.SetDeltaEtaRange(0.5, 1.);
      Correlators.push_back(c18);

      Correlator c19(3);
      c19.SetCollSystemAndEnergy("AuAu200GeV");
      c19.SetCentrality(40., 80.);
      c19.SetTriggerRange(3., 6.);
      c19.SetAssociatedRange(1.5, 999.);
      c19.SetDeltaEtaRange(1., 1.5);
      Correlators.push_back(c19);

      //Figure 2
      // d+Au Charged mesons 
      Correlator c20(4);
      c20.SetCollSystemAndEnergy("dAu200GeV");
      c20.SetCentrality(40., 80.);
      c20.SetTriggerRange(3., 6.);
      c20.SetAssociatedRange(1.5, 999.);
      c20.SetDeltaEtaRange(-1.5, -1.);
      Correlators.push_back(c20);

      Correlator c21(4);
      c21.SetCollSystemAndEnergy("dAu200GeV");
      c21.SetCentrality(40., 80.);
      c21.SetTriggerRange(3., 6.);
      c21.SetAssociatedRange(1.5, 999.);
      c21.SetDeltaEtaRange(-1., -0.5);
      Correlators.push_back(c21);

      Correlator c22(4);
      c22.SetCollSystemAndEnergy("dAu200GeV");
      c22.SetCentrality(40., 80.);
      c22.SetTriggerRange(3., 6.);
      c22.SetAssociatedRange(1.5, 999.);
      c22.SetDeltaEtaRange(-0.5, 0.);
      Correlators.push_back(c22);

      Correlator c23(4);
      c23.SetCollSystemAndEnergy("dAu200GeV");
      c23.SetCentrality(40., 80.);
      c23.SetTriggerRange(3., 6.);
      c23.SetAssociatedRange(1.5, 999.);
      c23.SetDeltaEtaRange(0., 0.5);
      Correlators.push_back(c23);

      Correlator c24(4);
      c24.SetCollSystemAndEnergy("dAu200GeV");
      c24.SetCentrality(40., 80.);
      c24.SetTriggerRange(3., 6.);
      c24.SetAssociatedRange(1.5, 999.);
      c24.SetDeltaEtaRange(0.5, 1.);
      Correlators.push_back(c24);

      Correlator c25(4);
      c25.SetCollSystemAndEnergy("dAu200GeV");
      c25.SetCentrality(40., 80.);
      c25.SetTriggerRange(3., 6.);
      c25.SetAssociatedRange(1.5, 999.);
      c25.SetDeltaEtaRange(1., 1.5);
      Correlators.push_back(c25);

//Figure 2
//Cu+Cu charged mesons
      Correlator c26(5);
      c26.SetCollSystemAndEnergy("CuCu200GeV");
      c26.SetCentrality(0., 20.);
      c26.SetTriggerRange(3., 6.);
      c26.SetAssociatedRange(1.5, 999.);
      c26.SetDeltaEtaRange(-1.5, -1.);
      Correlators.push_back(c26);

      Correlator c27(5);
      c27.SetCollSystemAndEnergy("CuCu200GeV");
      c27.SetCentrality(0., 20.);
      c27.SetTriggerRange(3., 6.);
      c27.SetAssociatedRange(1.5, 999.);
      c27.SetDeltaEtaRange(-1., -0.5);
      Correlators.push_back(c27);

      Correlator c28(5);
      c28.SetCollSystemAndEnergy("CuCu200GeV");
      c28.SetCentrality(0., 20.);
      c28.SetTriggerRange(3., 6.);
      c28.SetAssociatedRange(1.5, 999.);
      c28.SetDeltaEtaRange(-0.5, 0.);
      Correlators.push_back(c28);

      Correlator c29(5);
      c29.SetCollSystemAndEnergy("CuCu200GeV");
      c29.SetCentrality(0., 20.);
      c29.SetTriggerRange(3., 6.);
      c29.SetAssociatedRange(1.5, 999.);
      c29.SetDeltaEtaRange(0., 0.5);
      Correlators.push_back(c29);

      Correlator c30(5);
      c30.SetCollSystemAndEnergy("CuCu200GeV");
      c30.SetCentrality(0., 20.);
      c30.SetTriggerRange(3., 6.);
      c30.SetAssociatedRange(1.5, 999.);
      c30.SetDeltaEtaRange(0.5, 1.);
      Correlators.push_back(c30);

      Correlator c31(5);
      c31.SetCollSystemAndEnergy("CuCu200GeV");
      c31.SetCentrality(0., 20.);
      c31.SetTriggerRange(3., 6.);
      c31.SetAssociatedRange(1.5, 999.);
      c31.SetDeltaEtaRange(1., 1.5);
      Correlators.push_back(c31);
//Figure 2
//Au+Au charged mesons
      Correlator c32(6);
      c32.SetCollSystemAndEnergy("AuAu200GeV");
      c32.SetCentrality(40., 80.);
      c32.SetTriggerRange(3., 6.);
      c32.SetAssociatedRange(1.5, 999.);
      c32.SetDeltaEtaRange(-1.5, -1.);
      Correlators.push_back(c32);

      Correlator c33(6);
      c33.SetCollSystemAndEnergy("AuAu200GeV");
      c33.SetCentrality(40., 80.);
      c33.SetTriggerRange(3., 6.);
      c33.SetAssociatedRange(1.5, 999.);
      c33.SetDeltaEtaRange(-1., -0.5);
      Correlators.push_back(c33);

      Correlator c34(6);
      c34.SetCollSystemAndEnergy("AuAu200GeV");
      c34.SetCentrality(40., 80.);
      c34.SetTriggerRange(3., 6.);
      c34.SetAssociatedRange(1.5, 999.);
      c34.SetDeltaEtaRange(-0.5, 0.);
      Correlators.push_back(c34);

      Correlator c35(6);
      c35.SetCollSystemAndEnergy("AuAu200GeV");
      c35.SetCentrality(40., 80.);
      c35.SetTriggerRange(3., 6.);
      c35.SetAssociatedRange(1.5, 999.);
      c35.SetDeltaEtaRange(0., 0.5);
      Correlators.push_back(c35);

      Correlator c36(6);
      c36.SetCollSystemAndEnergy("AuAu200GeV");
      c36.SetCentrality(40., 80.);
      c36.SetTriggerRange(3., 6.);
      c36.SetAssociatedRange(1.5, 999.);
      c36.SetDeltaEtaRange(0.5, 1.);
      Correlators.push_back(c36);

      Correlator c37(6);
      c37.SetCollSystemAndEnergy("AuAu200GeV");
      c37.SetCentrality(40., 80.);
      c37.SetTriggerRange(3., 6.);
      c37.SetAssociatedRange(1.5, 999.);
      c37.SetDeltaEtaRange(1., 1.5);
      Correlators.push_back(c37);

//Figure 4
////d+Au charged mesons
 Correlator c38(1);
      c38.SetCollSystemAndEnergy("dAu200GeV");
      c38.SetCentrality(40., 80.);
      c38.SetTriggerRange(3., 3.5);
      c38.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c38);
//Figure 4
//d+Au charged baryons
      Correlator c39(2);
      c39.SetCollSystemAndEnergy("dAu200GeV");
      c39.SetCentrality(40., 80.);
      c39.SetTriggerRange(3.5, 4.);
      c39.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c39);

// Figure 4
//Cu+Cu charged mesons(empty)
      Correlator c40(1);
      c40.SetCollSystemAndEnergy("CuCu200GeV");
      c40.SetCentrality(0., 60.);
      c40.SetTriggerRange(2., 2.5);
      c40.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c40);

      Correlator c41(1);
      c41.SetCollSystemAndEnergy("CuCu200GeV");
      c41.SetCentrality(0., 60.);
      c41.SetTriggerRange(2.5, 3.);
      c41.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c41);

      Correlator c42(1);
      c42.SetCollSystemAndEnergy("CuCu200GeV");
      c42.SetCentrality(0., 60.);
      c42.SetTriggerRange(3., 3.5);
      c42.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c42);

      Correlator c43(1);
      c43.SetCollSystemAndEnergy("CuCu200GeV");
      c43.SetCentrality(0., 60.);
      c43.SetTriggerRange(3.5, 4.);
      c43.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c43);

      // Figure 4
      // Cu+Cu charged baryons(empty)       
      Correlator c44(2);
      c44.SetCollSystemAndEnergy("CuCu200GeV");
      c44.SetCentrality(0., 60.);
      c44.SetTriggerRange(0., 2.5);
      c44.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c44);

      Correlator c45(2);
      c45.SetCollSystemAndEnergy("CuCu200GeV");
      c45.SetCentrality(0., 60.);
      c45.SetTriggerRange(2.5, 3.);
      c45.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c45);

      Correlator c46(2);
      c46.SetCollSystemAndEnergy("CuCu200GeV");
      c46.SetCentrality(0., 60.);
      c46.SetTriggerRange(3., 3.5);
      c46.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c46);

      Correlator c47(2);
      c47.SetCollSystemAndEnergy("CuCu200GeV");
      c47.SetCentrality(0., 60.);
      c47.SetTriggerRange(3.5, 4.);
      c47.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c47);  

//Figure 4
//Au+Au charged mesons
      Correlator c48(1);
      c48.SetCollSystemAndEnergy("AuAu200GeV");
      c48.SetCentrality(40., 80.);
      c48.SetTriggerRange(2., 2.5);
      c48.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c48);

      Correlator c49(1);
      c49.SetCollSystemAndEnergy("AuAu200GeV");
      c49.SetCentrality(40., 80.);
      c49.SetTriggerRange(2.5, 3.);
      c49.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c49);

      Correlator c50(1);
      c50.SetCollSystemAndEnergy("AuAu200GeV");
      c50.SetCentrality(40., 60.);
      c50.SetTriggerRange(3., 3.5);
      c50.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c50);

      Correlator c51(1);
      c51.SetCollSystemAndEnergy("AuAu200GeV");
      c51.SetCentrality(40., 80.);
      c51.SetTriggerRange(3.5, 4.);
      c51.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c51);

      Correlator c52(1);
      c52.SetCollSystemAndEnergy("AuAu200GeV");
      c52.SetCentrality(40., 80.);
      c52.SetTriggerRange(4.5, 5.);
      c52.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c52);
//Figure 4
//Au+Au charged baryons
      Correlator c53(2);
      c53.SetCollSystemAndEnergy("AuAu200GeV");
      c53.SetCentrality(40., 80.);
      c53.SetTriggerRange(2.5, 3.);
      c53.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c53);

      Correlator c54(2);
      c54.SetCollSystemAndEnergy("AuAu200GeV");
      c54.SetCentrality(40., 80.);
      c54.SetTriggerRange(3.5, 4.);
      c54.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c54);

      Correlator c55(2);
      c55.SetCollSystemAndEnergy("AuAu200GeV");
      c55.SetCentrality(40., 80.);
      c55.SetTriggerRange(4.5, 5.);
      c55.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c55);

//There is no HEPDATA for centralities 0-12% Au+Au collisions 
    /*  // Figure 4
      // Au+Au charged mesons
      Correlator c1(1);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetCentrality(0., 12.);
      c2.SetTriggerRange(3.5, 4.);
      c2.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c2); 

      Correlator c1(1);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetCentrality(0., 12.);
      c2.SetTriggerRange(4.5, 5.);
      c2.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c2);

      //Figure 4
      //Au+Au charged baryons
      Correlator c1(1);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetCentrality(0., 12.);
      c2.SetTriggerRange(3., 3.5);
      c2.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c2);

      Correlator c1(1);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetCentrality(0., 12.);
      c2.SetTriggerRange(3.5, 4.);
      c2.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c2);

      Correlator c1(1);
      c2.SetCollSystemAndEnergy("AuAu200GeV");
      c2.SetCentrality(0., 12.);
      c2.SetTriggerRange(4.5, 5.);
      c2.SetAssociatedRange(1.5, 999.);
      Correlators.push_back(c2);
    */

  
     //Figure 6
     //Cu+Cu charged mesons
      Correlator c56(1);
      c56.SetCollSystemAndEnergy("CuCu200GeV");
      c56.SetCentrality(0., 60.);
      c56.SetTriggerRange(3., 6.);
      c56.SetAssociatedRange(1.2, 1.4);
      Correlators.push_back(c56);

      Correlator c57(1);
      c57.SetCollSystemAndEnergy("CuCu200GeV");
      c57.SetCentrality(0., 60.);
      c57.SetTriggerRange(3., 6.);
      c57.SetAssociatedRange(1.6, 1.8);
      Correlators.push_back(c57);

      Correlator c58(1);
      c58.SetCollSystemAndEnergy("CuCu200GeV");
      c58.SetCentrality(0., 60.);
      c58.SetTriggerRange(3., 6.);
      c58.SetAssociatedRange(2.2, 2.4);
      Correlators.push_back(c58);


      //Figure 6
      //d+Au charged mesons
      Correlator c59(2);
      c59.SetCollSystemAndEnergy("dAu200GeV");
      //c59.SetCentrality(0., 60.);
      c59.SetTriggerRange(3., 6.);
      c59.SetAssociatedRange(1.2, 1.4);
      Correlators.push_back(c59);

      Correlator c60(2);
      c60.SetCollSystemAndEnergy("dAu200GeV");
      //c60.SetCentrality(0., 60.);
      c60.SetTriggerRange(3., 6.);
      c60.SetAssociatedRange(1.6, 1.8);
      Correlators.push_back(c60);

      Correlator c61(2);
      c61.SetCollSystemAndEnergy("dAu200GeV");
      //c61.SetCentrality(0., 60.);
      c61.SetTriggerRange(3., 6.);
      c61.SetAssociatedRange(2.2, 2.4);
      Correlators.push_back(c61);
 
      //Figure 6
      //Cu+Cu charged baryons
      Correlator c62(3);
      c62.SetCollSystemAndEnergy("CuCu200GeV");
      c62.SetCentrality(0., 60.);
      c62.SetTriggerRange(3., 6.);
      c62.SetAssociatedRange(1.2, 1.4);
      Correlators.push_back(c62);

      Correlator c63(3);
      c63.SetCollSystemAndEnergy("CuCu200GeV");
      c63.SetCentrality(0., 60.);
      c63.SetTriggerRange(3., 6.);
      c63.SetAssociatedRange(1.6, 1.8);
      Correlators.push_back(c63);

      //Figure 6
      //d+Au charged baryons
      Correlator c64(4);
      c64.SetCollSystemAndEnergy("dAu200GeV");
      c64.SetCentrality(40., 80.);
      c64.SetTriggerRange(3., 6.);
      c64.SetAssociatedRange(1.2, 1.4);
      Correlators.push_back(c64);

      Correlator c65(4);
      c65.SetCollSystemAndEnergy("dAu200GeV");
      c65.SetCentrality(40., 80.);
      c65.SetTriggerRange(3., 6.);
      c65.SetAssociatedRange(1.6, 1.8);
      Correlators.push_back(c65);

  /*      // Book histograms 
	// book(_h["0111"], 1, 1, 1);
	//(Fig.2)
	//correlations functions for charged baryons and charged mesons for 0-20% Cu-Cu collisions
	//and 40-80% Au+Au collisions  
	 book(_h["0211"], 2, 1, 1);
	 book(_h["0212"], 2, 1, 2);
	 book(_h["0213"], 2, 1, 3);
	 book(_h["0214"], 2, 1, 4);
	 book(_h["0215"], 2, 1, 5);
	 book(_h["0216"], 2, 1, 6);
	//(Fig.3)
	//ratio of baryons/mesons in 0-60% Cu+Cu collisions
	 book(_h["0311"], 3, 1, 1);
	//(Fig.4)
	//jet-like yeild as a function of trigger particles for charged mesons and charged baryons
	//in minimum bias d+Au collisions 
	 book(_h["0411"], 4, 1, 1);
	 book(_h["0412"], 4, 1, 2);
	//(Fig.4)
	//jet-like yeild as a function of trigger particles for charged mesons and charged baryons
	//in 0-60% Cu+Cu collision
	 book(_h["0511"], 5, 1, 1);
	 book(_h["0512"], 5, 1, 2);
	//(Fig.4)
	//jet-like yeild as a function of trigger particles for charged mesons and charged baryons
	//in 40-80% Au+Au collision
	 book(_h["0611"], 6, 1, 1);
	 book(_h["0612"], 6, 1, 2);
// 	//(Fig.4)
// 	//jet-like yeild as a function of trigger particles for charged mesons and charged baryons
// 	//in 0-12% Au+Au collision
// 	 book(_h["0711"], 7, 1, 1);
// 	//(Fig.5)
// 	//centrality dependence of the jet-like yield of charged mesons and charged baryons in Cu+Cu
// 	//collisions 
// 	 book(_h["0811"], 8, 1, 1);
// 	//(Fig.5)
// 	//centrality dependence of the jet-like yield of charged mesons and charged baryons in Au+Au
// 	//collisions
// 	 book(_h["0911"], 9, 1, 1);
// 	//(Fig.5)
// 	//centrality dependence of the jet-like yield of charged mesons and charged baryons in d+Au
// 	//collisions
// 	 book(_h["1011"], 10, 1, 1);
 	//(Fig.6)
  	//jet-like yeild as a function of associated trigger particles for charged mesons and charged
 	//baryons in d+Au and 0-60% Cu+Cu collisions  
 	 book(_h["1111"], 11, 1, 1);
 	 book(_h["1112"], 11, 1, 2);
 	 book(_h["1113"], 11, 1, 3);
 	 book(_h["1114"], 11, 1, 4);
*/ 

for(unsigned int i = 1; i<= Correlators.size(); i++)
      {
          book(sow[i],"sow" + to_string(i));
          book(_h["031" + to_string(i)], 3, 1, i);
          book(_h["041" + to_string(i)], 4, 1, i);
          book(_h["061" + to_string(i)], 6, 1, i);
          book(_DeltaPhi[i], "DeltaPhi" + to_string(i), 24, 0, M_PI);
          book(_DeltaPhiSub[i], "DeltaPhiSub" + to_string(i), 24, 0, M_PI);
          book(_DeltaEta[i], "DeltaEta" + to_string(i), 20, 0, 2);
      }
          nEvents.assign(Correlators.size()+1, 0); 
          nTriggers.assign(Correlators.size()+1, 0); 
      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const ChargedFinalState& cfsTrig = apply<ChargedFinalState>(event, "CFSTrig");

//==================================================
//Select the histograms accordingly to the collision system, beam energy and centrality
//WARNING: d-Au is now implemented
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
     if (beam.first.pid() == 100029053 && beam.second.pid() == 1000010020)
      {
	CollSystem = "dAu";
      }
     if (beam.first.pid() == 1000010020  && beam.second.pid() == 100029053)
      {
        CollSystem = "dAu";
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
      const CentralityProjection& centrProj = apply<CentralityProjection>(event, "V0M");
      double centr = centrProj();

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
double dEta = deltaEta(pTrig, pAssoc);

for(Correlator& corr : Correlators)
                {
                    if(!corr.CheckConditionsMaxTrigger(SysAndEnergy, centr, pTrig.pt()/GeV, pAssoc.pt()/GeV)) continue;
                    
                    double etaCorrection = EtaEffCorrection(abs(dEta), *_DeltaEta[corr.GetIndex()]);
                    
                    if(abs(dPhi) < 0.78)
                    {
                        _h["031" + to_string(corr.GetIndex())]->fill(-abs(dEta), 0.5/etaCorrection);
                    }
                    
                    if(abs(dEta) < 1.78)
                    {
                        _h["041" + to_string(corr.GetIndex())]->fill(-abs(dPhi), 0.5);
                        _h["061" + to_string(corr.GetIndex())]->fill(-abs(dPhi), 0.5/etaCorrection);
                        _DeltaPhi[corr.GetIndex()]->fill(abs(dPhi), 0.5/etaCorrection);
                        _DeltaPhiSub[corr.GetIndex()]->fill(abs(dPhi), 0.5/etaCorrection);
                    }
                    
                }//end of correlations loop
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
                _h["031" + to_string(i)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                _h["041" + to_string(i)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                _h["061" + to_string(i)]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                
                _DeltaPhi[i]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                _DeltaPhiSub[i]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                
                vector<int> n{2,3};
                SubtractBackground(*_DeltaPhi[i], *_h["061" + to_string(i)], n, 0.63, 2.51);
                SubtractBackground(*_DeltaPhi[i], *_DeltaPhiSub[i], n, 0.63, 2.51);
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
    bool fillTrigger = true;
    vector<int> nTriggers;
    vector<int> nEvents;
    vector<Correlator> Correlators;

  };

// The hook for the plugin system
RIVET_DECLARE_PLUGIN(STAR_2016_I1429700);


}



/*/ Veto event for too large centralities since those are not used in the analysis at all
      if ((centr < 0.) || (centr > 80.)){
        vetoEvent;
      }


    // loop over charged final state particles
      for(const Particle& p : cfs.particles()) {

	       // protections against mc generators decaying long-lived particles
	        if (!( p.hasAncestorWith(Cuts::pid == 310) || p.hasAncestorWith(Cuts::pid == -310)  ||     // K0s
	           p.hasAncestorWith(Cuts::pid == 130)  || p.hasAncestorWith(Cuts::pid == -130)  ||     // K0l
	           p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) ||     // Xi0
	           p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) ||     // Lambda
	           p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) ||     // Sigma+/-
	           p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) ||     // Xi-/+
	           p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334) ))    // Omega-/+
	        {
          if (centr >= 0. && centr <= 5.) {
            switch (p.pid()) {
              	case 211: // pi+
              	  {

              	    break;
              	  }
              	case -211: // pi-
              	  {

              	    break;
              	  }
              	case 2212: // proton
              	  {

              	    break;
              	  }
              	case -2212: // anti-proton
                  {

              	    break;
              	  }
              	case 321: // K+
              	  {

              	    break;
              	  }
              	case -321: // K-
              	  {

              	    break;
              	  }
              }
          }

          if (centr > 5. && centr <= 10.) {
            	switch (p.pid()) {
                      	case 211: // pi+
                      	  {

                      	    break;
                      	  }
                      	case -211: // pi-
                      	  {

                      	    break;
                      	  }
                      	case 2212: // proton
                      	  {

                      	    break;
                      	  }
                      	case -2212: // anti-proton
                          {

                      	    break;
                      	  }
                      	case 321: // K+
                      	  {

                      	    break;
                      	  }
                      	case -321: // K-
                      	  {

                      	    break;
                      	  }
                }
            }

          if (centr > 10. && centr <= 20.) {
            switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }

          if (centr > 20. && centr <= 30.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }
          if (centr > 30. && centr <= 40.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }
          if (centr > 40. && centr <= 50.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }
          if (centr > 50. && centr <= 60.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }
          if (centr > 60. && centr <= 70.) {
            switch (p.pid()) {
              case 211: // pi+
                {

                  break;
                }
              case -211: // pi-
                {

                  break;
                }
              case 2212: // proton
                {

                  break;
                }
              case -2212: // anti-proton
                {

                  break;
                }
              case 321: // K+
                {

                  break;
                }
              case -321: // K-
                {

                  break;
                }
            }
          }
          if (centr > 70. && centr <= 80.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }
          if (centr > 80. && centr <= 90.) {
              switch (p.pid()) {
                case 211: // pi+
                  {

                    break;
                  }
                case -211: // pi-
                  {

                    break;
                  }
                case 2212: // proton
                  {

                    break;
                  }
                case -2212: // anti-proton
                  {

                    break;
                  }
                case 321: // K+
                  {

                    break;
                  }
                case -321: // K-
                  {

                    break;
                  }
              }
          }


	          // particle switch
         } // primaty pi, K, p only

      } // particle loop

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double norm = sumOfWeights() *2.*M_PI;
      //scale(_h["0111"], 1./norm);

    }

    //@}


    /// name Histograms
    //{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(STAR_2016_I1429700);


}
*/
