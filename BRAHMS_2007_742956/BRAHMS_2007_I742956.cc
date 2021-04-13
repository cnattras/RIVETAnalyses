// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "math.h"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BRAHMS_2007_I742956 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BRAHMS_2007_I742956);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
       
        const FinalState fs(Cuts::abseta < 5);//The cut is for particles we want thus <5 correponds to particles with 5 GeV or less.
        declare(fs, "fs");
        
    
        book(_h["CrsSecPIplus"], 1, 1, 1);
        book(_h["CrsSecPIminus"], 2, 1, 1);
        book(_h["CrsSecKplus"], 3, 1, 1);
        book(_h["CrsSecKminus"], 4, 1, 1);
        book(_h["CrsSecP"], 5, 1, 1);
        book(_h["CrsSecAntiP"], 6, 1, 1);

        book(_h["CrsSecPIplus2"], 7, 1, 1);
        book(_h["CrsSecPIminus2"], 8, 1, 1);
        book(_h["CrsSecKplus2"], 9, 1, 1);
        book(_h["CrsSecKminus2"], 10, 1, 1);
        book(_h["CrsSecP2"], 11, 1, 1);
        book(_h["CrsSecAntiP2"], 12, 1, 1);
        
        
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

        Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();

        for(const Particle& p : fsParticles)   //Particle Looping
        {

            double ySize = .1;
            
            if((abs(p.rapidity())<3.0)&&(abs(p.rapidity())<2.9)){ // Interested in regions near 2.95
                if(p.pid() == 211) _h["CrsSecPIplus"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize)); //Pion
                if(p.pid() == -211) _h["CrsSecPIminus"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize));
            
                if(p.pid() == 321) _h["CrsSecKplus"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize)); //Kaon
                if(p.pid() == -321) _h["CrsSecKminus"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize));
            
                if(p.pid() == 2212) _h["CrsSecP"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize)); //Pion
                if(p.pid() == -2212) _h["CrsSecAntiP"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize));
            }
            if((abs(p.rapidity())<3.35)&&(abs(p.rapidity())<3.25)){ // Interested in regions near 2.95
                if(p.pid() == 211) _h["CrsSecPIplus2"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize)); //Pion
                if(p.pid() == -211) _h["CrsSecPIminus2"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize));
<<<<<<< HEAD
            
                if(p.pid() == 321) _h["CrsSecKplus2"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize)); //Kaon
                if(p.pid() == -321) _h["CrsSecKminus2"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize));
            
=======
            
                if(p.pid() == 321) _h["CrsSecKplus2"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize)); //Kaon
                if(p.pid() == -321) _h["CrsSecKminus2"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize));
            
>>>>>>> dae383376ddc13a159f6cd722bd477cfb9cd81d7
                if(p.pid() == 2212) _h["CrsSecP2"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize)); //Pion
                if(p.pid() == -2212) _h["CrsSecAntiP2"]->fill(p.pT()/GeV,(1/p.pT())*(1/(2*pi))*(1/ySize));
            }
        
        }
        

    }


    /// Normalise histograms etc., after the run
    void finalize() {
<<<<<<< HEAD
        
        scale(_h["CrsSecPIplus"], ((crossSection()*pow(10,-6))/sumOfWeights())); // Still missing mean pT, pTbinsize and general nEvents, 1/Lum, also what does scale do?
        scale(_h["CrsSecPIminus"], ((crossSection()*pow(10,-6))/sumOfWeights()));
        scale(_h["CrsSecKplus"], ((crossSection()*pow(10,-6))/sumOfWeights()));
        scale(_h["CrsSecKminus"], ((crossSection()*pow(10,-6))/sumOfWeights()));
        scale(_h["CrsSecP"], ((crossSection()*pow(10,-6))/sumOfWeights()));
        scale(_h["CrsSecAntiP"], ((crossSection()*pow(10,-6))/sumOfWeights()));
        scale(_h["CrsSecPIplus2"], ((crossSection()*pow(10,-6))/sumOfWeights()));
        scale(_h["CrsSecPIminus2"], ((crossSection()*pow(10,-6))/sumOfWeights()));
        scale(_h["CrsSecKplus2"], ((crossSection()*pow(10,-6))/sumOfWeights()));
        scale(_h["CrsSecKminus2"], ((crossSection()*pow(10,-6))/sumOfWeights()));
        scale(_h["CrsSecP2"], ((crossSection()*pow(10,-6))/sumOfWeights()));
        scale(_h["CrsSecAntiP2"], ((crossSection()*pow(10,-6))/sumOfWeights()));
=======
//
//      normalize(_h["PtDist"]); // normalize to unity
//      normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
//
//        scale(_h["CSPiplus"], (crossSection()/picobarn/sumOfWeights())); // Still missing mean pT, pTbinsize and general nEvents, 1/Lum, also what does scale do?
//        scale(_h["CSPiminus"],(crossSection()/picobarn/sumOfWeights()));
//        scale(_h["CSKplus"], (crossSection()/picobarn/sumOfWeights()));
//        scale(_h["CSKminus"], (crossSection()/picobarn/sumOfWeights()));
//        scale(_h["CSPPt"], (crossSection()/picobarn/sumOfWeights()));
//        scale(_h["CSAntiPPt"], (crossSection()/picobarn/sumOfWeights()));
//        scale(_h["CSPiplus2"], (crossSection()/picobarn/sumOfWeights()));
//        scale(_h["CSPiminus2"], (crossSection()/picobarn/sumOfWeights()));
//        scale(_h["CSKplus2"], (crossSection()/picobarn/sumOfWeights()));
//        scale(_h["CSKminus2"], (crossSection()/picobarn/sumOfWeights()));
//        scale(_h["CSPPt2"], (crossSection()/picobarn/sumOfWeights()));
//        scale(_h["CSAntiPPt2"], (crossSection()/picobarn/sumOfWeights()));
>>>>>>> dae383376ddc13a159f6cd722bd477cfb9cd81d7
        
        
    }
 
    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}


  };


  DECLARE_RIVET_PLUGIN(BRAHMS_2007_I742956);

}

//        _h["CrsSecPIplus"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
//        _h["CrsSecPIminus"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
//        _h["CrsSecKplus"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
//        _h["CrsSecKminus"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
//        _h["CrsSecP"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
//        _h["CrsSecAntiP"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
//        _h["CrsSecPIplus2"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
//        _h["CrsSecPIminus2"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
//        _h["CrsSecKplus2"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
//        _h["CrsSecKminus2"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
//        _h["CrsSecP2"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
//        _h["CrsSecAntiP2"]->scaleW((crossSection()/picobarn)/_c["sow"]->sumW());
