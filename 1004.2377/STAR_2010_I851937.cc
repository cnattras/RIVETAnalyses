// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include "../Centralities/RHICCentrality.hh"


#define _USE_MATH_DEFINES
using namespace std;
namespace Rivet {


  /// @brief Add a short analysis description here
  class STAR_2010_I851937 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2010_I851937);


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

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

        const ChargedFinalState cfs(Cuts::abseta < 0.35);
        declare(cfs, "CFS");
        
        const PrimaryParticles pp(pdgPi0, Cuts::abseta < 0.35);
        declare(pp, "PP");
        
        const PromptFinalState pfs(Cuts::abseta < 0.35 && Cuts::pid == 22);
        declare(pfs, "PFS"); }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const PrimaryParticles& ppTrigPi0 = apply<PrimaryParticles>(event, "PP");
      const PromptFinalState& pfsTrigPhotons = apply<PromptFinalState>(event, "PFS");
      //==================================================
      // Select the histograms accordingly to the collision system, beam energy and centrality
      // WARNING: Still not implemented for d-Au
      //==================================================
      double nNucleons = 0.;
      string CollSystem = "Empty";
      const ParticlePair& beam = beams();
      CollSystem = "pp";
      nNucleons = 1.;
      //if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000010020) CollSystem = "dAu";
      //if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000290630) CollSystem = "dAu";
                                                                              
      string cmsEnergy = "Empty";
      if (fuzzyEquals(sqrtS()/GeV, 200*nNucleons, 1E-3)) cmsEnergy = "200GeV";
      string SysAndEnergy = CollSystem + cmsEnergy;
	}

    /// Normalise histograms etc., after the run
    void finalize() {
	}

    std::initializer_list<int> pdgPi0 = {111, -111};  // Pion 0
    std::initializer_list<int> pdgPhoton = {22};  // Pion 0
    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}


  };


  DECLARE_RIVET_PLUGIN(STAR_2010_I851937);

}
