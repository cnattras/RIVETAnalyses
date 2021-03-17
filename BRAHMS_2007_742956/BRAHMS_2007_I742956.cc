// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

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
//        const FinalState fs(Cuts::abseta > 5.0);
        
        const PromptFinalState pfs(Cuts::pT > 5.0);
              declare(pfs, "pfs");

          // Book histograms
          // specify custom binning
        vector<double> pTedges{0.35, 0.4, 0.45, 0.5, 0.55, .6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 3.0, 3.4, 3.8, 4.2, 5.2};

        //First rapidity y = 2.95
        book(_h["CSPi+Pt"], "CSPi+Pt", pTedges);
        book(_h["CSPi-Pt"], "CSPi-Pt", pTedges);
        book(_h["CSK+Pt"], "CSK+Pt", pTedges);
        book(_h["CSK-Pt"], "CSK-Pt", pTedges);
        book(_h["CSPPt"], "CSPPt", pTedges);
        book(_h["CSAntiPPt"], "CSAntiPPt", pTedges);
        
        book(_h["CSPi+Pt2"], "CSPi+Pt2", pTedges);
        book(_h["CSPi-Pt2"], "CSPi-Pt2", pTedges);
        book(_h["CSK+Pt2"], "CSK+Pt2", pTedges);
        book(_h["CSK-Pt2"], "CSK-Pt2", pTedges);
        book(_h["CSPPt2"], "CSPPt2", pTedges);
        book(_h["CSAntiPPt2"], "CSAntiPPt2", pTedges);
    
        
        //Second rapidity y =3.3
        book(_h["CrsSecPI+"], 1, 1, 1);
        book(_h["CrsSecPI-"], 2, 1, 1);
        book(_h["CrsSecK+"], 3, 1, 1);
        book(_h["CrsSecK-"], 4, 1, 1);
        book(_h["CrsSecP"], 5, 1, 1);
        book(_h["CrsSecAntiP"], 6, 1, 1);

        book(_h["CrsSecPI+R2"], 7, 1, 1);
        book(_h["CrsSecPI-R2"], 8, 1, 1);
        book(_h["CrsSecK+R2"], 9, 1, 1);
        book(_h["CrsSecPR2"], 10, 1, 1);
        book(_h["CrsSecAntPR2"], 11, 1, 1);
        book(_h["CrsSecK+R2"], 12, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

        Particles pfsParticles = applyProjection<FinalState>(event,"pfs").particles();

        for(const Particle& p : pfsParticles) //Particle Looping
        {
            if((abs(p.rapidity())/GeV<3.0)&&(abs(p.rapidity())/GeV<2.9)){ // Interested in regions near 2.95
                if(p.pid() == 211) _h["CSPi+Pt"]->fill(p.pT()/GeV); //Pion
                if(p.pid() == -211) _h["CSPi-Pt"]->fill(p.pT()/GeV);
                
                if(p.pid() == 321) _h["CSK+Pt"]->fill(p.pT()/GeV); //Kaon
                if(p.pid() == -321) _h["CSK-Pt"]->fill(p.pT()/GeV);
                
                if(p.pid() == 2212) _h["CSPPt"]->fill(p.pT()/GeV); //Pion
                if(p.pid() == -2212) _h["CSAntiPPt"]->fill(p.pT()/GeV);
            }
            if((abs(p.rapidity())/GeV<3.35)&&(abs(p.rapidity())/GeV<3.25)){ // Interested in regions near 2.95
                if(p.pid() == 211) _h["CSPi+Pt2"]->fill(p.pT()/GeV); //Pion
                if(p.pid() == -211) _h["CSPi-Pt2"]->fill(p.pT()/GeV);
                
                if(p.pid() == 321) _h["CSK+Pt2"]->fill(p.pT()/GeV); //Kaon
                if(p.pid() == -321) _h["CSK-Pt2"]->fill(p.pT()/GeV);
                
                if(p.pid() == 2212) _h["CSPPt2"]->fill(p.pT()/GeV); //Pion
                if(p.pid() == -2212) _h["CSAntiPPt2"]->fill(p.pT()/GeV);
            }
            
        }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
//
//      normalize(_h["XXXX"]); // normalize to unity
//      normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
//      scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

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
