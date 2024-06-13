// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2020_I1773662 : public Analysis {
    public:

      /// Constructor
      DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2020_I1773662);


      /// @name Analysis methods
      //@{

      /// Book histograms and initialise projections before the run
      void init() {
        declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

        // Initialise and register projections
        //
        // const FinalState fs(Cuts::abseta < 4.9);
        // PromptFinalState bare_leps(Cuts::abspid == PID::MUON)
        const UnstableParticles up(Cuts::abseta < 0.5 && Cuts::absrap<=2.2 && Cuts::absrap>1.2 &&  Cuts::pT > 10*GeV);// && Cuts::abspid == 443);
        //const UnstableParticles up(Cuts::abseta < 0.5 && Cuts::abspid == 10443);

        declare(up, "up");
        beamOpt = getOption<string>("beam","NONE");

        if(beamOpt=="PP510") collsys = pp510;
        else if(beamOpt=="AUAU200") collsys = AuAu200;


        //figure 5- paper
        book(_h["Fig5_Sigma"], 1, 1, 1);
        //figure 6-paper
        book(_h["Fig6_Sigma"], 2, 1, 1);
        book(_c["sow_pp"], "sow_pp");

      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {

        const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
        const double c = cent();

      
        Particles upParticles = applyProjection<UnstableParticles>(event,"up").particles();

        if(collsys==pp510)
        {
          if(c < 5.) 
          {
            _c["sow_pp"]->fill();
          }


          //double absrap = p.absrap();
          for(const Particle& p : upParticles) 
          {
            if(p.pid() == 443)// && absrap<=2.2 && absrap>1.2) 
              _h["Fig5_Sigma"]->fill(p.pT()/GeV, 1.0);
          }//end of for loop

        }
        
      }


      /// Normalise histograms etc., after the run
      void finalize() {

        //double scale = 1./(2*M_PI);// for invariant yield
        // normalize(_h["Fig5_Sigma"], crossSection()/nanobarn/sumW());
        scale(_h["Fig5_Sigma"], crossSection()/nanobarn/sumW());


     
      }

      //@}


      /// @name Histograms
      //@{
      map<string, Histo1DPtr> _h;
      map<string, Profile1DPtr> _p;
      map<string, CounterPtr> _c;
      string beamOpt = "";
      enum CollisionSystem {pp510, AuAu200};
      CollisionSystem collsys;

      //@}


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2020_I1773662);

}
