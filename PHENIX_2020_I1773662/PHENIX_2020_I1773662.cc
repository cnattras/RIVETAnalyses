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


        // // book RIVET counter histograms
        // book(_c["sow_pp"], sow_pp);
        // book(_c["sow_AuAu2040"], sow_AuAu2040);
        // book(_c["sow_AuAu4060"], sow_AuAu4060);
        // book(_c["sow_AuAu6094"], sow_AuAu6094);



        /*     // The basic final-state projection:
        // all final-state particles within
        // the given eta acceptance
        // const FinalState fs(Cuts::abseta < 4.9);

        // The final-state particles declared above are clustered using FastJet with
        // the anti-kT algorithm and a jet-radius parameter 0.4
        // muons and neutrinos are excluded from the clustering
        // FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
        // declare(jetfs, "jets");

        // FinalState of prompt photons and bare muons and electrons in the event
        //  PromptFinalState photons(Cuts::abspid == PID::PHOTON);
        // PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

        // Dress the prompt bare leptons with prompt photons within dR < 0.1,
        // and apply some fiducial cuts on the dressed leptons
        //  Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
        // DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
        // declare(dressed_leps, "leptons");

        // Missing momentum
        //declare(MissingMomentum(fs), "MET");

        // Book histograms
        // specify custom binning
        // book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
        // book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
        // book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
        // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
        // book(_h["AAAA"], 1, 1, 1);
        // book(_p["BBBB"], 2, 1, 1);
        // book(_c["CCCC"], 3, 1, 1);

*/
      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {

        const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
        const double c = cent();

        /*     
               if(c < 20.)   
               {
               i  _c["sow_AuAu0020"]->fill();
               }
               else if(c >= 20. && c < 40.)
               {
               _c["sow_AuAu2040"]->fill();
               }
        //  else if(c >= 40. && c < 60.)
        // 	{
        // 	  _c["sow_AuAu4060"]->fill();
        // 	}
        //  else if(c >= 60. && c < 94.)
        // 	{
        // 	  _c["sow_AuAu6094"]->fill();
        // 	}
        */  

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
        /*  // Retrieve dressed leptons, sorted by pT
            vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();

        // Retrieve clustered jets, sorted by pT, with a minimum pT cut
        Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

        // Remove all jets within dR < 0.2 of a dressed lepton
        idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

        // Select jets ghost-associated to B-hadrons with a certain fiducial selection
        Jets bjets = filter_select(jets, [](const Jet& jet) {
        return  jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
        });

        // Veto event if there are no b-jets
        if (bjets.empty())  vetoEvent;

        // Apply a missing-momentum cut
        if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent;

        // Fill histogram with leading b-jet pT
        _h["XXXX"]->fill(bjets[0].pT()/GeV);

*/
      }


      /// Normalise histograms etc., after the run
      void finalize() {

        //double scale = 1./(2*M_PI);// for invariant yield
        // normalize(_h["Fig5_Sigma"], crossSection()/nanobarn/sumW());
        scale(_h["Fig5_Sigma"], crossSection()/nanobarn/sumW());



        //double factor = 

        /* normalize(_c["sow_AuAu2040"]);


           normalize(_h["XXXX"]); // normalize to unity
           normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
           scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

*/
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
