// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
//#include "Rivet/Projections/FastJets.hh"
//#include "Rivet/Projections/DressedLeptons.hh"
//#include "Rivet/Projections/MissingMomentum.hh"
//#include "Rivet/Projections/PromptFinalState.hh"
//#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <fstream>
#include <iostream>
#include <math.h>
#define _USE_MATH_DEFINES

using namespace std;

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2016_I1427723 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2016_I1427723);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // the basic final-state projection: all final-state particles within the given eta acceptance

      //const ChargedFinalState cfs(Cuts::abseta < 0.6 && Cuts::pT > 0*MeV);
      // Declare centrality projection
      const FinalState fs(Cuts::abseta < 0.6);
      declare(fs, "FS");
      const ChargedFinalState cfs(Cuts::abseta < 0.5);
      declare(cfs, "CFS");

      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");
      //this histogram is just the npart but it was interpreted by RIVET as another y value
      //book(_h["0111"], 1, 1, 1);
      //this histogram is the et
      book(_h["et"], 1, 1, 2);
      //this histogram is ET/(npart/2).  It is a bit redundant with the previous plot but often a more useful way to plot things
      book(_h["etovernpart"], 1, 1, 3);
      //this is ET/(nquark/2).  redundant with the previous histogram
      //book(_h["0114"], 1, 1, 4);
      //this is npart again
      //book(_h["0211"], 2, 1, 1);
      //this is nquark
      //book(_h["0212"], 2, 1, 2);
      //this is the area
      //book(_h["0213"], 2, 1, 3);
      //this is ET/Nch
      book(_h["etovernch"], 3, 1, 1);
      //ANTONIO - I need to add two histograms with the same centrality binning as the other histograms but to store (1) the number of events and (2) the charged particle multiplicity
      
      // the basic final-state projection: 
      // all final-state particles within 
      // the given eta acceptance

      // the final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      //FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      //declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      //PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      //PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

      // dress the prompt bare leptons with prompt photons within dR < 0.1
      // apply some fiducial cuts on the dressed leptons
      //Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
      //DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      //declare(dressed_leps, "leptons");

      // missing momentum
      //declare(MissingMomentum(fs), "MET");

      // The centrality bins upper bin edges.
      centralityBins = { 2.5,5.0,7.5,10., 12.5,15.0,17.5, 20.,22.5,25.0,27.5,30., 32.5,35.0,37.5,40., 42.5,45.0,47.5,50., 52.5,55.0,57.5,60., 62.5,65.0,67.5,70., 72.5,75.0,77.5,80. };

      npart = {396.1, 369.3, 342.4, 316.5, 292.5,
		 269.9, 248.8, 229.3, 210.8, 193.5,
		 177  , 161.9, 147.6, 134.5, 121.9,
		 110.2, 99.17, 89.06, 79.52, 71.07,
		 62.98, 55.51, 48.61, 42.56, 36.88,
		 31.78, 27.25, 23.17, 19.51, 16.42,
		 13.76, 11.41, 9.615, 8.011, 6.698,
		 5.616, 4.747, 4.004, 3.371, 2.755};
      //for (int i = 0; i < nCentralityBins; ++i) {
      //book(nEvents[centralityBins[i]], "nEvents_" + toString(i)); 
      //book(nEvents[centralityBins[i]], "nCh_" + toString(i)); 
	//nEvents[centralityBins[i]] = bookCounter("nEvents_" + toString(i));
	//nCh[centralityBins[i]] = bookCounter("nCh_" + toString(i));
      //}
      //for (int i = 0;i < nCentralityBins; ++i)
      //{
      // book(nEvents[centralityBins[i]], "nEvents_" + toString(i)); 
      //book(nEvents[centralityBins[i]], "nCh_" + toString(i)); 
	//nEvents[centralityBins[i]] = bookCounter("nEvents_" + toString(i));
	//nCh[centralityBins[i]] = bookCounter("nCh_" + toString(i));
      //}
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

      // Prepare centrality projection and value
      const CentralityProjection& centrProj = apply<CentralityProjection>(event, "V0M");
      double centr = centrProj();

      // Veto event for too large centralities since those are not used in the analysis at all
      if ((centr < 0.) || (centr > 80.)){
        vetoEvent;
      }

//       // Find the correct centrality bin and count the number of events
      //auto eventsItr = nEvents.upper_bound(centr);
       //if (eventsItr == nEvents.end()) return;
       //eventsItr->second->fill();

      //cout<<"Hi"<<endl;
      //calculate ET in the event
      int centbin = ((int)centr/2.5);
      float mynpart = 1;
      if(centbin<nCentralityBins) mynpart = npart[centbin];
      float et = 0;
      for(const Particle& p : fs.particles()) {
	et = et+ p.momentum().Et()/GeV;
      }
      float mult = 0;
      for(const Particle& p : cfs.particles()) {
	mult++;
      }
      _h["et"]->fill(centr,et);
      _h["etovernpart"]->fill(centr,et/(mynpart/2));
      _h["etovernch"]->fill(centr,et/mult);
      //ANTONIO - this is where I would fill the event and Nch histograms
      //cout<<"et "<<et<<" mult "<<mult<<" et/mult "<<et/mult<<" centrality "<<centr<<" cent bin "<<((int)centr/2.5)<<endl;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      //ANTONIO - here is where I would need to (1) scale all histograms by the number of events including the Nch histogram (2) divide etovernpart by the npart (3) divide etovernch by nch

    }

    //@}

    vector<double> centralityBins;
    int nCentralityBins = 32;
    vector<double> npart;
    /// @name Histograms
    //@{
    map<string, Profile1DPtr> _h;
    //map<string, Profile1DPtr> _p;
    //map<double, CounterPtr> nEvents;
    //map<double, CounterPtr> nEvents;
    //map<double, CounterPtr> nCh;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2016_I1427723);


}
