
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
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
    
   //------------------------------------------------------------------------------------------------------------------------ 
  class Correlator {
      
    private:
    
      int _index;
      string _collSystemAndEnergy;
      pair<double,double> _centrality;
      pair<double,double> _triggerRange;
      pair<double,double> _associatedRange;
      vector<int> _pid;
      bool _noCentrality = false;
  
    public:
    
      /// Constructor
      Correlator(int index) {
          _index = index;
      }

      void SetCollSystemAndEnergy(string s){ _collSystemAndEnergy = s; }
      void SetCentrality(double cmin, double cmax){ _centrality = make_pair(cmin, cmax); }
      void SetNoCentrality(){ _noCentrality = true; }
      void SetTriggerRange(double tmin, double tmax){ _triggerRange = make_pair(tmin, tmax); }
      void SetAssociatedRange(double amin, double amax){ _associatedRange = make_pair(amin, amax); }
      void SetPID(std::initializer_list<int> pid){ _pid = pid; }
    
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
    bool CheckCentrality(double cent){ return ((cent>_centrality.first && cent<_centrality.second) || _noCentrality == true) ? true : false; }
    bool CheckTriggerRange(double tpt){ return (tpt>_triggerRange.first && tpt<_triggerRange.second) ? true : false; }
    bool CheckAssociatedRange(double apt){ return (apt>_associatedRange.first && apt<_associatedRange.second) ? true : false; }
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



  class STAR_2016_I1442357 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2016_I1442357);


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

    //Instead of directly taking the calculated Delta_eta, this method takes the center of the bin associated to the Delta_eta.
    //This is done to avoid values of zero efficiency when |Delta_eta| is close to maxDeltaEta
    double EtaEffCorrection(double deltaEta, YODA::Histo1D& hist)
    {
        double maxDeltaEta = 2.;
        
        int binEta = hist.binIndexAt(deltaEta);
        
        if(binEta < 0)
        {
            MSG_INFO("Eta is out of bounds!");
            return 0.;
        }
        
        double binCenterEta = hist.bin(binEta).xMid();
        
        return 1. - (1./maxDeltaEta)*abs(binCenterEta);
    }

    //bmin and bmax are included in the integral. Range = [bmin, bmax]
    //Give the min and max bins to calculate the integral
    double GetYieldInBinRange(YODA::Histo1D& hist, int bmin, int bmax)
    {
        double integral = 0.;
        
        if(bmin < 0 || bmax >= (int)hist.numBins())
        {
            MSG_ERROR("Out of range!");
            return 0.;
        }
        
        for(int i = bmin; i <= bmax; i++)
        {
            integral += hist.bin(i).sumW();
        }
        
        return integral;
        
    }

    //vmin is included in the integral, but vmax is not. Range = [vmin, vmax[
    //Give the min and max values to calculate the integral
    double GetYieldInUserRange(YODA::Histo1D& hist, double vmin, double vmax)
    {
        double integral = 0.;
        
        if(hist.binIndexAt(vmin) < 0 || hist.binIndexAt(vmax) < 0)
        {
            MSG_ERROR("Out of range!");
            return 0.;
        }
        
        for(int i = hist.binIndexAt(vmin); i < hist.binIndexAt(vmax); i++)
        {
            integral += hist.bin(i).sumW();
        }
        
        return integral;
        
    }
    
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
    
    bool FindHistoZT(YODA::Histo1D hist, double zT, string &histoID)
    {
                
        for(auto &bin : hist.bins())
        {
            if(zT < bin.xMin() || zT > bin.xMax()) continue;
            else
            {
                histoID += to_string(bin.xMin()) + "_" + to_string(bin.xMax());
                return true;
            }
        }
        
        return false;
        
    }


    /// Book histograms and initialise projections before the run
    void init() {             
      // the basic final-state projection: all final-state particles within the given eta acceptance

      const ChargedFinalState cfs(Cuts::abseta < 1.0 && Cuts::pT > 1*GeV);
      declare(cfs, "CFS");
      const ChargedFinalState cfsTrig(Cuts::abseta < 1.0 && Cuts::pT > 2*GeV);
      declare(cfsTrig, "CFSTrig");
      // FinalState of prompt photons and bare muons and electrons in the event
      const PromptFinalState pfsTrig(Cuts::abspid == PID::PHOTON);
      declare(pfsTrig, "PFSTrig");
      // pi0 also
      
      const PrimaryParticles primp(pdgPi0, Cuts::abseta < 1.);
      declare(primp, "PP");
        
      const PromptFinalState pfs(Cuts::abseta < 1. && Cuts::pid == 22);
      declare(pfs, "PFS");
      
      beamOpt = getOption<string>("beam","NONE");


      if(beamOpt=="PP") collSys = pp;
      else if(beamOpt=="AUAU") collSys = AuAu;
      


      // Declare centrality projection
      declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");


      //==================================================
      // Create one correlator for each set of Collisions System / Beam Energy / Centrality Interval / Trigger pT interval / Associated pT interval
      // The number used in the constructor is the y-axis index of the histograms (d00-x00-y00) of the paper
      // Ex.: Correlator c1(1); -> is the correlator for histograms _h["0311"], _h["0411"], etc
      // Ex.: Correlator c2(2); -> is the correlator for histograms _h["0312"], _h["0412"], etc
       //==================================================

    
    Correlator b1(1);
	b1.SetCollSystemAndEnergy("AuAu200GeV");
	b1.SetCentrality(0., 12.);
	b1.SetTriggerRange(1., 20.);  //pT_trig
	b1.SetAssociatedRange(1.2, 3.); // zT
	b1.SetPID(pdgPhoton);
	CorrelatorsB.push_back(b1);

	Correlator b2(2);
    b2.SetCollSystemAndEnergy("AuAu200GeV");
    b2.SetCentrality(0., 12.);
    b2.SetTriggerRange(1., 20.);  //pT_trig
    b2.SetAssociatedRange(1.2, 3.); // zT
    b2.SetPID(pdgPi0);
    CorrelatorsB.push_back(b2);

	Correlator b3(3);
    b3.SetCollSystemAndEnergy("AuAu200GeV");
    b3.SetCentrality(0., 12.);
    b3.SetTriggerRange(1., 20.);  //pT_trig
    b3.SetAssociatedRange(3., 5.); // zT
    b3.SetPID(pdgPhoton);
    CorrelatorsB.push_back(b3);

	Correlator b4(4);
    b4.SetCollSystemAndEnergy("AuAu200GeV");
    b4.SetCentrality(0., 12.);
    b4.SetTriggerRange(1., 20.);  //pT_trig
    b4.SetAssociatedRange(3., 5.); // zT
    b4.SetPID(pdgPi0);
    CorrelatorsB.push_back(b4);

	Correlator b5(5);
    b5.SetCollSystemAndEnergy("pp200GeV");
    b5.SetNoCentrality();
    b5.SetTriggerRange(1., 20.);  //pT_trig
    b5.SetAssociatedRange(1.2, 3.); // zT
    b5.SetPID(pdgPhoton);
    CorrelatorsB.push_back(b5);

	Correlator b6(6);
    b6.SetCollSystemAndEnergy("pp200GeV");
    b6.SetNoCentrality();
    b6.SetTriggerRange(1., 20.);  //pT_trig
    b6.SetAssociatedRange(1.2, 3.); // zT
    b6.SetPID(pdgPi0);
    CorrelatorsB.push_back(b6);

	Correlator b7(7);
    b7.SetCollSystemAndEnergy("pp200GeV");
    b7.SetNoCentrality();
    b7.SetTriggerRange(1., 20.);  //pT_trig
    b7.SetAssociatedRange(3., 5.); // zT
    b7.SetPID(pdgPhoton);
    CorrelatorsB.push_back(b7);

	Correlator b8(8);
    b8.SetCollSystemAndEnergy("pp200GeV");
    b8.SetNoCentrality();
    b8.SetTriggerRange(1., 20.);  //pT_trig
    b8.SetAssociatedRange(3., 5.); // zT
    b8.SetPID(pdgPi0);
    CorrelatorsB.push_back(b8);
    
    //This will book the histograms of figure 1
    for(Correlator& corr : CorrelatorsB)
	{
        book(_h["0" + to_string(corr.GetIndex()) + "11"], corr.GetIndex(), 1, 1);
        book(sow[corr.GetIndex()],"sow" + to_string(corr.GetIndex()));
        nTriggers[corr.GetIndex()] = 0;
    }
    

	Correlator c9(9);	//11,1 & 13,1 are the same
    c9.SetCollSystemAndEnergy("AuAu200GeV");
    c9.SetCentrality(0., 12.);
    c9.SetTriggerRange(1., 20.);  //pT_trig
    c9.SetAssociatedRange(1.2, 999.); // zT
    c9.SetPID(pdgPi0);
    CorrelatorsPi0.push_back(c9);

	Correlator c10(10);	//12,1 & 14,1 are the same
    c10.SetCollSystemAndEnergy("pp200GeV");
    c10.SetNoCentrality();
    c10.SetTriggerRange(1., 20.);  //pT_trig
    c10.SetAssociatedRange(1.2, 999.); // zT
    c10.SetPID(pdgPi0);
    CorrelatorsPi0.push_back(c10);
    
    for(Correlator& corr : CorrelatorsPi0)
	{
        book(_h["0" + to_string(corr.GetIndex()) + "11"], corr.GetIndex(), 1, 1);
        book(_h["0" + to_string(corr.GetIndex()+2) + "11"], corr.GetIndex(), 1, 1);
        
        string refname = mkAxisCode(corr.GetIndex(), 1, 1);
        const Histo1D& refdata = refData(refname);
        for(auto &bin : refdata.bins())
        {
            book(_DeltaPhizT["0" + to_string(corr.GetIndex()) + "11zT_" + to_string(bin.xMin()) + "_" + to_string(bin.xMax())], "DeltaPhi_0" + to_string(corr.GetIndex()) + "11zT_" + to_string(bin.xMin()) + "_" + to_string(bin.xMax()), 24, -M_PI/2., 3.*M_PI/2.);
        }
        
        book(sow[corr.GetIndex()],"sow" + to_string(corr.GetIndex()));
        nTriggers[corr.GetIndex()] = 0;
    }
	
	Correlator c11(13);
	c11.SetCollSystemAndEnergy("AuAu200GeV");
    c11.SetCentrality(0., 12.);
    c11.SetTriggerRange(1., 20.);  //pT_trig
    c11.SetAssociatedRange(1.2, 999.); // zT
    c11.SetPID(pdgPhoton);
    CorrelatorsPhoton.push_back(c11);

	Correlator c12(14);	//18,1 is the same
    c12.SetCollSystemAndEnergy("pp200GeV");
    c12.SetNoCentrality();
    c12.SetTriggerRange(1., 20.);  //pT_trig
    c12.SetAssociatedRange(1.2, 999.); // zT
    c12.SetPID(pdgPhoton);
    CorrelatorsPhoton.push_back(c12);
    
    for(Correlator& corr : CorrelatorsPhoton)
	{
        book(_h["0" + to_string(corr.GetIndex()) + "11"], corr.GetIndex(), 1, 1);
        
        string refname = mkAxisCode(corr.GetIndex(), 1, 1);
        const Histo1D& refdata = refData(refname);
        for(auto &bin : refdata.bins())
        {
            book(_DeltaPhizT["0" + to_string(corr.GetIndex()) + "11zT_" + to_string(bin.xMin()) + "_" + to_string(bin.xMax())], "DeltaPhi_0" + to_string(corr.GetIndex()) + "11zT_" + to_string(bin.xMin()) + "_" + to_string(bin.xMax()), 24, -M_PI/2., 3.*M_PI/2);
        }
        
        book(sow[corr.GetIndex()],"sow" + to_string(corr.GetIndex()));
        nTriggers[corr.GetIndex()] = 0;
    }


    /*
	for(Correlator& corr : Correlators)
	{	
		//book(sow[corr.GetIndex()],"sow" + corr.GetIndex());
        cout<<"booking for"<<corr.GetIndex()<<endl;

		if(corr.GetIndex() == 9 || corr.GetIndex() == 10)
		{
			book(_h["0" + to_string(corr.GetIndex()) + "1" + "1"], corr.GetIndex(), 1, 1);
			book(_h["0" + to_string(corr.GetIndex()+2) + "1" + "1"], corr.GetIndex()+2, 1, 1);
			book(_h["0" + to_string(corr.GetIndex()+4) + "1" + "1"], corr.GetIndex()+4, 1, 1);
		}
		
		else if(corr.GetIndex() == 11)
		{
			book(_h["0" + to_string(corr.GetIndex()+4) + "1" + "1"], corr.GetIndex()+4, 1, 1);
		}

		 else if(corr.GetIndex() == 12 || corr.GetIndex() == 13)
                {
                        book(_h["0" + to_string(corr.GetIndex()+4) + "1" + "1"], corr.GetIndex()+4, 1, 1);
			            book(_h["0" + to_string(corr.GetIndex()+6) + "1" + "1"], corr.GetIndex()+6, 1, 1);
                }

		 else if(corr.GetIndex() == 14 || corr.GetIndex() == 15)
                {
                        book(_h["0" + to_string(corr.GetIndex()+6) + "1" + "1"], corr.GetIndex()+6, 1, 1);
                }

		else
		{
            cout<<"d="<<corr.GetIndex()<<endl;
			book(_h["0" + to_string(corr.GetIndex()) + "1" + "1"], corr.GetIndex(), 1, 1);
		}

	}
	*/
/* Old Booking

      for(unsigned int i = 1; i<= Correlators.size(); i++)
      {
        book(sow[i],"sow" + to_string(i));
      //  book(_DeltaPhi[i], "DeltaPhi" + to_string(i), 24, 0, M_PI);
     //   book(_DeltaPhiSub[i], "DeltaPhiSub" + to_string(i), 24, 0, M_PI);
      

     	 for(unsigned int j = 1; j<= 21; j++)
     	 {
        	  book(_h["0" + to_string(j) + "11"], j, 1, 1); 
     	 }
      
     	 nEvents.assign(Correlators.size()+1, 0); 
     	 nTriggers.assign(Correlators.size()+1, 0); 
      } //moved from being above second for statement
    }
      // Book histograms
*/
/*
	 // *****FIGURE 1, Delta_phi(rad) // The Azimuthal correlation functions of charged hadrons per trigger,12 < pT^trig < 20 GeV/c
   ////// DO NOT NEED THIS FIGURE ///////////
   //per trigger Yields (rad^-1)
   //Data from Fig. 1a
   book(_h["0111"], 1, 1, 1); //AuAu, 1.2<pT^assoc<3 GeV/c, gamma
   book(_h["0211"], 2, 1, 1); //AuAu, 1.2<pT^assoc<3 GeV/c, pi0
   //Data from Fig. 1b
	 book(_h["0311"], 3, 1, 1); //AuAu, 3<pT^assoc<5 GeV/c, gamma
   book(_h["0411"], 4, 1, 1); //AuAu, 3<pT^assoc<5 GeV/c, pi0
   //Data from Fig. 1c
   book(_h["0511"], 5, 1, 1); //pp, 1.2<pT^assoc<3 GeV/c, gamma
   book(_h["0611"], 6, 1, 1); //pp, 1.2<pT^assoc<3 GeV/c, pi0
   //Data from Fig. 1d
   book(_h["0711"], 7, 1, 1); //pp, 3<pT^assoc<5 GeV/c, gamma
   book(_h["0811"], 8, 1, 1); //pp, 3<pT^assoc<5 GeV/c, pi0
   
   // *****FIGURE 2, zT //pi0-hadron correlations**************************************************
   //D(zT)
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   //Data from Fig. 2a
   book(_h["0911"], 9, 1, 1);  //Away-side, Au+Au
   book(_h["1011"], 10, 1, 1); //Away-side, p+p
   //Data from Fig. 2b
   book(_h["1111"], 11, 1, 1); //Near-side, Au+Au
   book(_h["1211"], 12, 1, 1); //Near-side, p+p
   
   // *****FIGURE 3, zT //Direct Photon-hadron correlations****************************************
   //D(zT)
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   book(_h["1311"], 13, 1, 1); //Away-side, AuAu 
   book(_h["1411"], 14, 1, 1); //Away-side, pp
   
   // *****FIGURE 4, zT_corr **********************************************************************
   //D(zT_Corr)
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   book(_h["1511"], 15, 1, 1); //Pi0-hadron correlations, Away-side, pp
   
   // *****FIGURE 5, zT ***************************************************************************
   //AuAu
   //I_AA
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   book(_h["1611"], 16, 1, 1); //gamma, Away-side, I_AA
   book(_h["1711"], 17, 1, 1); //pi0, Away-side, I_AA
   
   // *****FIGURE 6, zT ***************************************************************************
   //AuAu
   //Ratio
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   book(_h["1811"], 18, 1, 1); //gamma
   book(_h["1911"], 19, 1, 1); //pi0
   
   // *****FIGURE 7 *******************************************************************************
   //I_AA
   //12<pT^trig<20 GeV/c
   //1.2<pT^assoc GeV/c
   book(_h["2011"], 20, 1, 1); //gamma, Away-side, I_AA (pT^trig)
   book(_h["2111"], 21, 1, 1); //gamma, Away-side, I_AA (pT^assoc)
   
   */
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
        

        //cout<<"We are analyzing"<<endl; // NEW 7/15/20

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const PromptFinalState& pfsTrig = apply<PromptFinalState>(event, "PFSTrig");
      const ChargedFinalState& cfsTrig = apply<ChargedFinalState>(event, "CFSTrig");
      
      const PrimaryParticles& ppTrigPi0 = apply<PrimaryParticles>(event, "PP");
      const PromptFinalState& pfsTrigPhotons = apply<PromptFinalState>(event, "PFS");

      //==================================================
      // Select the histograms accordingly to the collision system, beam energy and centrality
      // WARNING: Still not implemented for d-Au
      //==================================================
      
      /*
      double nNucleons = 0.;
      string CollSystem = "Empty";
      const ParticlePair& beam = beams();
      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
          CollSystem = "AuAu";
          nNucleons = 197.;
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
      if(cmsEnergy.compare("Empty") == 0) return;
      */
      
      
      
      if(collSys == pp) SysAndEnergy = "pp200GeV";
      else if(collSys == AuAu) SysAndEnergy = "AuAu200GeV";
      


      /// DO NOT MESS WITH THIS
      // Prepare centrality projection and value
      const CentralityProjection& centrProj = apply<CentralityProjection>(event, "CMULT");
      double centr = centrProj();
      
      //cout<<"We are analyzing  2"<<endl;
      
      // Veto event for too large centralities since those are not used in the analysis at all
      if ((centr < 0.) || (centr > 12.)){
        vetoEvent;
      }
    

    double triggerptMin = 1.;
    double triggerptMax = 20.;
    double associatedptMin = 1.2;
    double associatedptMax = 999.;

    bool isVeto = true;
    
    for(Correlator& corr : CorrelatorsB)
    {
        if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
        if(!corr.CheckCentrality(centr)) continue;
        
        sow[corr.GetIndex()]->fill();
        
        isVeto = false;
    }
    
    for(Correlator& corr : CorrelatorsPi0)
    {
        if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
        if(!corr.CheckCentrality(centr)) continue;
        
        sow[corr.GetIndex()]->fill();
        
        isVeto = false;
    }
    
    for(Correlator& corr : CorrelatorsPhoton)
    {
        if(!corr.CheckCollSystemAndEnergy(SysAndEnergy)) continue;
        if(!corr.CheckCentrality(centr)) continue;
        
        sow[corr.GetIndex()]->fill();
        
        isVeto = false;
    }
    
    if(isVeto) vetoEvent;
    
    //PHOTON
    // loop over charged final state particles
    for(const Particle& pTrig : pfsTrigPhotons.particles()) {
        
        if(pTrig.pt()/GeV < triggerptMin || pTrig.pt()/GeV > triggerptMax) continue;
        
        //Check if is secondary
        if(isSecondary(pTrig)) continue;
                      
            for(Correlator& corr : CorrelatorsB)
            {
                if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
                if(!corr.CheckPID(pdgPhoton)) continue;
                nTriggers[corr.GetIndex()]++;
            }
            
            for(Correlator& corr : CorrelatorsPhoton)
            {
                if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
                nTriggers[corr.GetIndex()]++;
            }

            for(const Particle& pAssoc : cfs.particles()) {
                
                //cout<<"z_T is "<<pAssoc.pt()/pTrig.pt()<<endl;   // NEW 7/15/20
                if(pAssoc.pt()/GeV < associatedptMin || pAssoc.pt()/GeV > associatedptMax) continue;
                
                //Check if Trigger and Associated are the same particle
                if(isSameParticle(pTrig,pAssoc)) continue;
                
                //Check if is secondary
                if(isSecondary(pAssoc)) continue;
                
                double dPhi = GetDeltaPhi(pTrig, pAssoc);
                
                for(Correlator& corr : CorrelatorsB)
                {
                    if(!corr.CheckPID(pdgPhoton)) continue;
                    if(!corr.CheckAssociatedRange(pAssoc.pt()/GeV)) continue;
                
                    if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
                        
                    _h["0" + to_string(corr.GetIndex()) + "11"]->fill(dPhi);                    
                }

                for(Correlator& corr : CorrelatorsPhoton)
                {
                    if(!corr.CheckPID(pdgPhoton)) continue;
                    if(!corr.CheckAssociatedRange(pAssoc.pt()/GeV)) continue;
                    if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
                
                    double zT = pAssoc.pT()/pTrig.pT();
                    string histoID = "zT_";
                    if(!FindHistoZT(*_h["0" + to_string(corr.GetIndex()) + "11"], zT, histoID)) continue;
                    _DeltaPhizT["0" + to_string(corr.GetIndex()) + "11" + histoID]->fill(dPhi);
                    
                }
                
                
            } // end of loop over associated particles
      } // particle loop
        
    //PI0
    // loop over charged final state particles
    for(const Particle& pTrig : ppTrigPi0.particles()) {
        
        if(pTrig.pt()/GeV < triggerptMin || pTrig.pt()/GeV > triggerptMax) continue;

        //cout << "pT_trigeer: " << pTrig.pt()/GeV << endl;
        
        //Check if is secondary
        if(isSecondary(pTrig)) continue;
                      
            for(Correlator& corr : CorrelatorsB)
            {
                if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
                if(!corr.CheckPID(pdgPi0)) continue;
                nTriggers[corr.GetIndex()]++;
            }
            
            for(Correlator& corr : CorrelatorsPi0)
            {
                if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
                nTriggers[corr.GetIndex()]++;
            }
            

            for(const Particle& pAssoc : cfs.particles()) {
                
                //cout<<"z_T is "<<pAssoc.pt()/pTrig.pt()<<endl;   // NEW 7/15/20
                if(pAssoc.pt()/GeV < associatedptMin || pAssoc.pt()/GeV > associatedptMax) continue;
                
                //Check if Trigger and Associated are the same particle
                if(isSameParticle(pTrig,pAssoc)) continue;
                
                //Check if is secondary
                if(isSecondary(pAssoc)) continue;
                
                double dPhi = GetDeltaPhi(pTrig, pAssoc);
                
                for(Correlator& corr : CorrelatorsB)
                {
                    if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
                    if(!corr.CheckPID(pdgPi0)) continue;
                    if(!corr.CheckAssociatedRange(pAssoc.pt()/GeV)) continue;
                        
                    _h["0" + to_string(corr.GetIndex()) + "11"]->fill(dPhi);
                }
                
                for(Correlator& corr : CorrelatorsPi0)
                {
                    if(!corr.CheckConditions(SysAndEnergy, centr, pTrig.pt()/GeV)) continue;
                    if(!corr.CheckPID(pdgPi0)) continue;
                    if(!corr.CheckAssociatedRange(pAssoc.pt()/GeV)) continue;
                
                    double zT = pAssoc.pT()/pTrig.pT();
                    string histoID = "zT_";
                    if(!FindHistoZT(*_h["0" + to_string(corr.GetIndex()) + "11"], zT, histoID)) continue;
                    _DeltaPhizT["0" + to_string(corr.GetIndex()) + "11" + histoID]->fill(dPhi);
                }
                
                
            } // end of loop over associated particles
      } // particle loop
        
        
        
        
        
        
        
    }


    /// Normalise histograms etc., after the run
    void finalize() {
        
        for(Correlator& corr : CorrelatorsB)
        {
            cout << "Entries: " << sow[corr.GetIndex()]->numEntries() << " Triggers: " <<  nTriggers[corr.GetIndex()] << endl;
            _h["0" + to_string(corr.GetIndex()) + "11"]->scaleW(sow[corr.GetIndex()]->numEntries()/(nTriggers[corr.GetIndex()]*sow[corr.GetIndex()]->sumW()));
        }
        
        
        /*
        for(unsigned int i = 1; i <= Correlators.size(); i++)
        {
            if(nTriggers[i] > 0)
            {
                _h[to_string(i) + "11"]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                
                _DeltaPhi[i]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                _DeltaPhiSub[i]->scaleW((double)nEvents[i]/(nTriggers[i]*sow[i]->sumW()));
                
                vector<int> n{2,3};
               // SubtractBackground(*_DeltaPhi[i], *_h["061" + to_string(i)], n, 0.63, 2.51);
               // SubtractBackground(*_DeltaPhi[i], *_DeltaPhiSub[i], n, 0.63, 2.51);
            }    
        }
        */
      //normalize correlation histograms by scaling by 1.0/(Ntrig*binwidthphi*binwidtheta) in each bin BUT also be careful when rebinning.  Probably best to FIRST add histograms for correlation functions THEN normalize
      //do background subtraction ala zyam
      //calculate yields

    } //end void finalize
    
    //Histograms and variables

    map<string, Histo1DPtr> _h;
    map<int, CounterPtr> sow;
    map<string, Histo1DPtr> _DeltaPhizT;
    map<int, Histo1DPtr> _DeltaPhi;
    map<int, int> nTriggers;
    vector<Correlator> CorrelatorsPi0;
    vector<Correlator> CorrelatorsPhoton;
    vector<Correlator> CorrelatorsB;
    
    initializer_list<int> pdgPi0={111, -111};
    initializer_list<int> pdgPhoton={22};
    enum CollisionSystem {pp, AuAu};
    CollisionSystem collSys;
    string SysAndEnergy = "";
    string beamOpt;

};
  DECLARE_RIVET_PLUGIN(STAR_2016_I1442357);
}
