// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2008_I776624 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2008_I776624);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {


      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      // // Initialise and register projections

      // // The basic final-state projection:
      // // all final-state particles within
      // // the given eta acceptance
      const UnstableParticles up(Cuts::abseta < 0.5 && Cuts::abspid == 443);
      
      // book HEPData histograms
      book(_h2D["YAA_pT_mid_020"], 1, 1, 1);  
      book(_h2D["YAA_pT_mid_2040"], 2, 1, 1);  
      book(_h2D["YAA_pT_mid_4060"], 3, 1, 1); 
      book(_h2D["YAA_pT_mid_6094"], 4, 1, 1); 
     
      // book(_h2D["YAA_pT_fwd_cent_1"], 5, 1, 1); 
      // book(_h2D["YAA_pT_fwd_cent_2"], 6, 1, 2); 
      // book(_h2D["YAA_pT_fwd_cent_3"], 7, 1, 3); 
      // book(_h2D["YAA_pT_fwd_cent_4"], 8, 1, 4); 

      // book(_h2D["ptsq_npart_mid"], 9, 1, 1); 
      // book(_h2D["ptsq_npart_fwd"], 10, 1, 1); 
   
      // book(_h2D["RAA_pT_020_mid"], 11, 1, 1); 
      // book(_h2D["RAA_pT_020_fwd"], 12, 1, 2);

      // book(_h2D["RAA_y_020"], 13, 1, 1); 
      // book(_h2D["RAA_npart_mid"], 14, 1, 1);
      // book(_h2D["RAA_npart_fwd"], 15, 1, 1);
      // book(_h2D["RAA_CNM_ratio"], 12, 1, 1);

      // // book RIVET counter histograms
      // book(_c["sow_AuAu0020"], sow_AuAu0020);
      // book(_c["sow_AuAu2040"], sow_AuAu2040);
      // book(_c["sow_AuAu4060"], sow_AuAu4060);
      // book(_c["sow_AuAu6094"], sow_AuAu6094);
      

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

     
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();
      // //const double r = rap();
      // // double absr = abs(r);    // need rapidity selection

      // if(c < 20.)   
      // 	{
      // 	  _c["sow_AuAu0020"]->fill();
      // 	}
      // else if(c >= 20. && c < 40.)
      // 	{
      // 	  _c["sow_AuAu2040"]->fill();
      // 	}
      //  else if(c >= 40. && c < 60.)
      // 	{
      // 	  _c["sow_AuAu4060"]->fill();
      // 	}
      //  else if(c >= 60. && c < 94.)
      // 	{
      // 	  _c["sow_AuAu6094"]->fill();
      // 	}
      

      Particles upParticles = applyProjection<UnstableParticles>(event,"up").particles();

      for(const Particle& p : upParticles) 
      	{
	  
	  if(c < 20. && p.pid() == 443)  // 11
	    _h2D["YAA_pT_mid_020"]->fill(p.pT()/GeV, 1.0);
	
	  if(c < 40. && c > 20. && p.pid() == 443)
	    _h2D["YAA_pT_mid_2040"]->fill(p.pT()/GeV, 1.0);

	  if(c < 60. && c > 40. && p.pid() == 443)
	    _h2D["YAA_pT_mid_4060"]->fill(p.pT()/GeV, 1.0);

	  if(c < 94. && c > 60. && p.pid() == 443)
	    _h2D["YAA_pT_mid_6094"]->fill(p.pT()/GeV, 1.0);
	
      	}// end for loop


    }


    /// Normalise histograms etc., after the run
    void finalize() {

  
      normalize(_h2D["YAA_pT_mid_020"]);
      normalize(_h2D["YAA_pT_mid_2040"]);
      normalize(_h2D["YAA_pT_mid_4060"]);
      normalize(_h2D["YAA_pT_mid_6094"]);
      
      // normalize(_h["XXXX"]); // normalize to unity
      // normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
      // scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Histo2DPtr> _h2D;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2008_I776624);

}
