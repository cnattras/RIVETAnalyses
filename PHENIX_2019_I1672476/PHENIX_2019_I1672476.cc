// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include <math.h>
#define _USE_MATH_DEFINES
namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2019_I1672476 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2019_I1672476);
    

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
      // Initialise and register projections
      
      const PromptFinalState pfs(Cuts::abseta < 0.35 && Cuts::pid == 22);
      declare(pfs, "pfs");
      const FinalState fs(Cuts::abseta < 0.35 && Cuts::abscharge > 0);
      declare(fs, "fs");

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      

      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning
      book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
      book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
      book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)

      book(_h["AuAu62_c0-20"], 1, 1, 1);
      book(_h["AuAu62_c0-86"], 1, 1, 2);
      book(_h["AuAu39_c0-86"], 2, 1, 1);
      book(_p["AuAu39_chPMult"], 3, 1, 1);
      book(_p["AuAu62_chPMult"], 4, 1, 1);
      book(_p["AuAu200_chPMult"], 4, 1, 2);      
     // book(_h["fig2-2-a"], 5, 1, 1);
     // book(_h["fig2-2-b"], 5, 1, 2);
     // book(_h["fig2-2-c"], 5, 1, 3);
      book(_h["fig3-1a"], 6, 1, 1);
      book(_h["fig3-1b"], 7, 1, 1);
      book(_h["fig3-1c"], 8, 1, 1);
      book(_h["fig3-1d"], 9, 1, 1);
      book(_h["fig3-1e"], 10, 1, 1);
      book(_h["fig3-1f"], 11, 1, 1);
      book(_h["fig3-2a"], 12, 1, 1);
      book(_h["fig3-2b"], 13, 1, 1);
      book(_h["fig3-2c"], 14, 1, 1);
      book(_h["fig3-2d"], 15, 1, 1);
      book(_h["fig3-2e"], 16, 1, 1);
      book(_h["fig3-2f"], 17, 1, 1);
      book(_h["fig4-1a"], 18, 1, 1);
      book(_h["fig4-1b"], 19, 1, 1);
      book(_h["fig4-1c"], 20, 1, 1);
     // book(_h["fig4-2-a"], 21, 1, 1);
     //
      book(sow["sow-fig1-1-a"],"sow-fig1-1-a");//these are currently unused counters
      book(sow["sow-fig1-1-b"],"sow-fig1-1-b");
      book(sow["sow-fig1-2"],"sow-fig1-2");
      book(sow["sow-fig2-1a"],"sow-fig2-1a");
      book(sow["sow-fig2-1b-a"],"sow-fig2-1b-a");
      book(sow["sow-fig2-1b-b"],"sow-fig2-1b-b");
      book(sow["sow-fig3-1a"],"sow-fig3-1a");
      book(sow["sow-fig3-1b"],"sow-fig3-1b");
      book(sow["sow-fig3-1c"],"sow-fig3-1c");
      book(sow["sow-fig3-1d"],"sow-fig3-1d");
      book(sow["sow-fig3-1e"],"sow-fig3-1e");
      book(sow["sow-fig3-1f"],"sow-fig3-1f");
      book(sow["sow-fig3-2a"],"sow-fig3-2a");
      book(sow["sow-fig3-2b"],"sow-fig3-2b");
      book(sow["sow-fig3-2c"],"sow-fig3-2c");
      book(sow["sow-fig3-2d"],"sow-fig3-2d");
      book(sow["sow-fig3-2e"],"sow-fig3-2e");
      book(sow["sow-fig3-2f"],"sow-fig3-2f");
      book(sow["sow-fig4-1a"],"sow-fig4-1a");
      book(sow["sow-fig4-1b"],"sow-fig4-1b");
      book(sow["sow-fig4-1c"],"sow-fig4-1c");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const PromptFinalState pfs = apply<PromptFinalState>(event, "pfs");
      const FinalState fs = apply<FinalState>(event, "fs");
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();
      const ParticlePair& beam = beams();
 
      int NN = 0;


      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
	 NN = 197.;
         if (fuzzyEquals(sqrtS()/GeV, 39*NN, 1E-3)) collSystem = AuAu39;
	 if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 1E-3)) collSystem = AuAu62;
      }
	 Particles photons = applyProjection<PromptFinalState>(event, "pfs").particles();
	 Particles chargedParticles = applyProjection<FinalState>(event, "fs").particles();
      if(collSystem == AuAu62)//for now I am attempting to only fill the first three histogram(AuAu62 cent0-20 and 0-86, and AuAu39 cent0-86), where pT is x-ax, and inv yield is y-ax
      {
         if((c >= 0.) && (c < 20.))
	 {
	    for(const Particle& p : photons)
	    {
		double partPt = p.pT()/GeV;
		//double pt_weight = 1./(partPt*2.*M_PI);
		_h["AuAu62_c0-20"]->fill(partPt); //I Copied this from Chrsital's cc file for PHENIX_2012_I1107625, But i think mine should have inv yield instead of pt_weight(not sure what this is)
	    }
	 }
	 else if((c >= 0.) && (c < 86. ))
	 {
	     for(const Particle& p : photons)
	     {
		 double partPt = p.pT()/GeV;
		 //double pt_weight = 1./(partPt*2.*M_PI);
		 _h["AuAu62_c0-86"]->fill(partPt);
	     }
	 }

      }

      if (collSystem == AuAu39)
      {
	  if((c >= 0.) && (c <86.))
	  {
	      for(const Particle& p : photons)
	      {
		  double partPt = p.pT()/GeV;
		  //double pt_weight = 1./(partPt*2.*M_PI);
		  _h["AuAu39_c0-86"]->fill(partPt);
	      }
	  }
      }
     

      int nCh39count = 0;
      int nCh62count = 0;
      int nCh200count = 0;
      int absEta = 0.7;
       

      //count charged particels for diff AuAu enegeries
      if (collSystem == AuAu39)
      {
	   nCh39count ++;   
      }
      else if(collSystem == AuAu62)
      {
	   nCh62count ++;
      }
      else if(collSystem == AuAu200)
      {
	   nCh200count ++;
      }

      //fill histos with cent c and charged particle count/abs val of Eta
      if (collSystem == AuAu39)
      {
	   _p["AuAu39_chPMult"]->fill(c,nCh39count/absEta);
      }  
      else if(collSystem == AuAu62)
      {
	   _p["AuAu62_chPMult"]->fill(c,nCh62count/absEta);
      }
      else if(collSystem ==AuAu200)
      {
	   _p["AuAu200_chPMult"]->fill(c,nCh200count/absEta);
      }






 

      // Retrieve dressed leptons, sorted by pT
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

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h["XXXX"]); // normalize to unity
      normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
      scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)
    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> sow;
    enum CollisionSystem {AuAu39, AuAu62, AuAu200};
    CollisionSystem collSystem;
    string beamOpt;

    //@}


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2019_I1672476);

}
