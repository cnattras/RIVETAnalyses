// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "../Centralities/RHICCentrality.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2024_IPAT : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2024_IPAT);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
//PAT: You will need a line like this but it should specify 
  // declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      const FinalState fs(Cuts::abseta < 0.9 && Cuts::pT > 0.15*GeV && Cuts::abscharge > 0);
      declare(fs, "fs");
      // associated particles
      const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::pT > 0.15*GeV && Cuts::abscharge > 0);
      declare(aprim, "aprim");
      //The second primary particles is so that the jet specctra goes over all particles instead of cutting out pT < 0.15GeV
      //const ALICE::PrimaryParticles aprimall(Cuts::abseta < 0.9 && Cuts::abscharge > 0 && (Cuts::abspid == Rivet::PID::PIPLUS || Cuts::abspid == Rivet::PID::KPLUS || Cuts::abspid == Rivet::PID::PROTON || Cuts::abspid == Rivet::PID::ELECTRON || Cuts::abspid == Rivet::PID::MUON));


      //Pat, change this so that it accepts charged primary particles OR photons
      //You also want to 
      //This is the jet constituents
      const ALICE::PrimaryParticles aprimall(Cuts::abseta < 0.7 && Cuts::abscharge > 0);
      declare(aprimall, "aprimall");


      // jets à la ALICE - Jet area will be available using the pseudojet
      fastjet::AreaType fjAreaType = fastjet::active_area_explicit_ghosts;
      fastjet::GhostedAreaSpec fjGhostAreaSpec = fastjet::GhostedAreaSpec(1., 1, 0.005, 1., 0.1, 1e-100);
      fjAreaDef = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
      FastJets jetfs(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.2, fjAreaDef);
      declare(jetfs, "jetsfs");

      // Book histograms

      //particle types: (pi+ + pi-), (p), (K+ + K-)
      //particle momentum: (1, 1.5), (1.5, 2), (2,3), (3,4), (4,5), (5,6), (6,10)
      // jet momentum: 20-40 GeV

      // Desired plots
      // Deta and Dphi for each particle type and particle momentum
      // Deta: 60 bins from -1.5 to 1.5
      // Dphi: 48 bins from -pi/2 to 3pi/2

      // Integrated yield for each particle type as a function of momentum
      // Yield: [(1, 1.5), (1.5, 2), (2,3), (3,4), (4,5), (5,6), (6,10)]
      // Ratio of K to pi and p to pi as a function of momentum
      // Ratio: [(1, 1.5), (1.5, 2), (2,3), (3,4), (4,5), (5,6), (6,10)]

      // specify custom binning

      // Booking proton histos
      // Momenta (1, 1.5)
      book(_histos["dphi_p_1_1.5"], "dphi_p_1_1.5", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_p_1_1.5"], "deta_p_1_1.5", 60, -1.5, 1.5);
      // Momenta (1.5, 2)
      book(_histos["dphi_p_1.5_2"], "dphi_p_1.5_2", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_p_1.5_2"], "deta_p_1.5_2", 60, -1.5, 1.5);
      // Momenta (2, 3)
      book(_histos["dphi_p_2_3"], "dphi_p_2_3", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_p_2_3"], "deta_p_2_3", 60, -1.5, 1.5);
      // Momenta (3, 4)
      book(_histos["dphi_p_3_4"], "dphi_p_3_4", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_p_3_4"], "deta_p_3_4", 60, -1.5, 1.5);
      // Momenta (4,5)
      book(_histos["dphi_p_4_5"], "dphi_p_4_5", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_p_4_5"], "deta_p_4_5", 60, -1.5, 1.5);
      // Momenta (5, 6)
      book(_histos["dphi_p_5_6"], "dphi_p_5_6", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_p_5_6"], "deta_p_5_6", 60, -1.5, 1.5);
      // Momenta (6, 10)
      book(_histos["dphi_p_6_10"], "dphi_p_6_10", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_p_6_10"], "deta_p_6_10", 60, -1.5, 1.5);
      // Yield
      book(_histos["yield_p"], "yield_p", {1, 1.5, 2, 3, 4, 5, 6, 10});


      // Booking pion histos
      // Momenta (1, 1.5)
      book(_histos["dphi_pi_1_1.5"], "dphi_pi_1_1.5", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_pi_1_1.5"], "deta_pi_1_1.5", 60, -1.5, 1.5);
      // Momenta (1.5, 2)
      book(_histos["dphi_pi_1.5_2"], "dphi_pi_1.5_2", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_pi_1.5_2"], "deta_pi_1.5_2", 60, -1.5, 1.5);
      // Momenta (2, 3)
      book(_histos["dphi_pi_2_3"], "dphi_pi_2_3", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_pi_2_3"], "deta_pi_2_3", 60, -1.5, 1.5);
      // Momenta (3, 4)
      book(_histos["dphi_pi_3_4"], "dphi_pi_3_4", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_pi_3_4"], "deta_pi_3_4", 60, -1.5, 1.5);
      // Momenta (4,5)
      book(_histos["dphi_pi_4_5"], "dphi_pi_4_5", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_pi_4_5"], "deta_pi_4_5", 60, -1.5, 1.5);
      // Momenta (5, 6)
      book(_histos["dphi_pi_5_6"], "dphi_pi_5_6", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_pi_5_6"], "deta_pi_5_6", 60, -1.5, 1.5);
      // Momenta (6, 10)
      book(_histos["dphi_pi_6_10"], "dphi_pi_6_10", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_pi_6_10"], "deta_pi_6_10", 60, -1.5, 1.5);
      // Yield
      book(_histos["yield_pi"], "yield_pi", {1, 1.5, 2, 3, 4, 5, 6, 10});


      // Booking Kaon histos
      // Momenta (1, 1.5)
      book(_histos["dphi_K_1_1.5"], "dphi_K_1_1.5", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_K_1_1.5"], "deta_K_1_1.5", 60, -1.5, 1.5);
      // Momenta (1.5, 2)
      book(_histos["dphi_K_1.5_2"], "dphi_K_1.5_2", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_K_1.5_2"], "deta_K_1.5_2", 60, -1.5, 1.5);
      // Momenta (2, 3)
      book(_histos["dphi_K_2_3"], "dphi_K_2_3", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_K_2_3"], "deta_K_2_3", 60, -1.5, 1.5);
      // Momenta (3, 4)
      book(_histos["dphi_K_3_4"], "dphi_K_3_4", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_K_3_4"], "deta_K_3_4", 60, -1.5, 1.5);
      // Momenta (4,5)
      book(_histos["dphi_K_4_5"], "dphi_K_4_5", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_K_4_5"], "deta_K_4_5", 60, -1.5, 1.5);
      // Momenta (5, 6)
      book(_histos["dphi_K_5_6"], "dphi_K_5_6", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_K_5_6"], "deta_K_5_6", 60, -1.5, 1.5);
      // Momenta (6, 10)
      book(_histos["dphi_K_6_10"], "dphi_K_6_10", 48 , - 3.14 / 2, 3 * 3.14 / 2);
      book(_histos["deta_K_6_10"], "deta_K_6_10", 60, -1.5, 1.5);
      // Yield
      book(_histos["yield_K"], "yield_K", {1, 1.5, 2, 3, 4, 5, 6, 10});

      // Book ratios
      book(_scatters["ratio_K_to_pi"], "ratio_K_to_pi", {1, 1.5, 2, 3, 4, 5, 6, 10});
      book(_scatters["ratio_p_to_pi"], "ratio_p_to_pi", {1, 1.5, 2, 3, 4, 5, 6, 10});

      // Counters (still need work)
      book(_histos["number_of_events"], "number_of_events", 1, 0.0, 1.0);
      book(_histos["number_of_jets"], "number_of_jets", 1, 0.0, 1.0);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      std::cout << "Event number: " << event.genEvent()->event_number() << std::endl;
      _histos["number_of_events"]->fill(0.5, 1.0);

      //Get final state particles, e.g. all particles "detected"
      const FinalState fs = apply<FinalState>(event, "fs");
      //Get jets
      FastJets jetsfs = apply<FastJets>(event, "jetsfs");
      //Get particles in ALICE EMCAL acceptance for jets
      const ALICE::PrimaryParticles aprimall = apply<ALICE::PrimaryParticles>(event, "aprimall");
      const Particles ALICEparticlesall = aprimall.particles();
      //Get particles in ALICE TPC acceptance for associated hadrons
      const ALICE::PrimaryParticles aprim = apply<ALICE::PrimaryParticles>(event, "aprim");
      const Particles ALICEparticles = aprim.particles();

    


      //jetsfs.calc(ALICEparticles);
      jetsfs.calc(ALICEparticlesall);

      Jets jets = jetsfs.jetsByPt(Cuts::abseta < 0.5 && Cuts::pT >= 5.*GeV && Cuts::pT <= 40.*GeV);
      std::cout << "Number of jets: " << jets.size() << std::endl;
      std::cout << "Number of particles: " << ALICEparticles.size() << std::endl;
      for(auto jet : jets)
      {
        // why is this not counting the number of events?
        _histos["number_of_jets"]->fill(0.5, 1.0);
        for(auto particle : ALICEparticles)
        {
            auto dphi = jet.phi()-particle.phi();
            auto deta = jet.eta()-particle.eta();
            // make sure dphi is in [-pi/2,3pi/2] range
            if (dphi > 3*3.14159/2) dphi -= 2*3.14159;
            if (dphi < -3*3.14159/2) dphi += 2*3.14159;
            //put logic to fill the appropriate associated hadron momentum
            // according to the definition above
            if (particle.pid() == PID::PIPLUS || particle.pid() == PID::PIMINUS)
            {
              if (particle.pt() > 1 && particle.pt() < 1.5)
              {
                _histos["dphi_pi_1_1.5"]->fill(dphi);
                _histos["deta_pi_1_1.5"]->fill(deta);
              }
              else if (particle.pt() > 1.5 && particle.pt() < 2)
              {
                _histos["dphi_pi_1.5_2"]->fill(dphi);
                _histos["deta_pi_1.5_2"]->fill(deta);
              }
              else if (particle.pt() > 2 && particle.pt() < 3)
              {
                _histos["dphi_pi_2_3"]->fill(dphi);
                _histos["deta_pi_2_3"]->fill(deta);
              }
              else if (particle.pt() > 3 && particle.pt() < 4)
              {
                _histos["dphi_pi_3_4"]->fill(dphi);
                _histos["deta_pi_3_4"]->fill(deta);
              }
              else if (particle.pt() > 4 && particle.pt() < 5)
              {
                _histos["dphi_pi_4_5"]->fill(dphi);
                _histos["deta_pi_4_5"]->fill(deta);
              }
              else if (particle.pt() > 5 && particle.pt() < 6)
              {
                _histos["dphi_pi_5_6"]->fill(dphi);
                _histos["deta_pi_5_6"]->fill(deta);
              }
              else if (particle.pt() > 6 && particle.pt() < 10)
              {
                _histos["dphi_pi_6_10"]->fill(dphi);
                _histos["deta_pi_6_10"]->fill(deta);
              }
            }
            if (particle.pid() == PID::PROTON)
            {
                if (particle.pt() > 1 && particle.pt() < 1.5)
              {
                _histos["dphi_p_1_1.5"]->fill(dphi);
                _histos["deta_p_1_1.5"]->fill(deta);
              }
              else if (particle.pt() > 1.5 && particle.pt() < 2)
              {
                _histos["dphi_p_1.5_2"]->fill(dphi);
                _histos["deta_p_1.5_2"]->fill(deta);
              }
              else if (particle.pt() > 2 && particle.pt() < 3)
              {
                _histos["dphi_p_2_3"]->fill(dphi);
                _histos["deta_p_2_3"]->fill(deta);
              }
              else if (particle.pt() > 3 && particle.pt() < 4)
              {
                _histos["dphi_p_3_4"]->fill(dphi);
                _histos["deta_p_3_4"]->fill(deta);
              }
              else if (particle.pt() > 4 && particle.pt() < 5)
              {
                _histos["dphi_p_4_5"]->fill(dphi);
                _histos["deta_p_4_5"]->fill(deta);
              }
              else if (particle.pt() > 5 && particle.pt() < 6)
              {
                _histos["dphi_p_5_6"]->fill(dphi);
                _histos["deta_p_5_6"]->fill(deta);
              }
              else if (particle.pt() > 6 && particle.pt() < 10)
              {
                _histos["dphi_p_6_10"]->fill(dphi);
                _histos["deta_p_6_10"]->fill(deta);
              }
            }
            if (particle.pid() == PID::KPLUS || particle.pid() == PID::KMINUS)
            {
                if (particle.pt() > 1 && particle.pt() < 1.5)
              {
                _histos["dphi_K_1_1.5"]->fill(dphi);
                _histos["deta_K_1_1.5"]->fill(deta);
              }
              else if (particle.pt() > 1.5 && particle.pt() < 2)
              {
                _histos["dphi_K_1.5_2"]->fill(dphi);
                _histos["deta_K_1.5_2"]->fill(deta);
              }
              else if (particle.pt() > 2 && particle.pt() < 3)
              {
                _histos["dphi_K_2_3"]->fill(dphi);
                _histos["deta_K_2_3"]->fill(deta);
              }
              else if (particle.pt() > 3 && particle.pt() < 4)
              {
                _histos["dphi_K_3_4"]->fill(dphi);
                _histos["deta_K_3_4"]->fill(deta);
              }
              else if (particle.pt() > 4 && particle.pt() < 5)
              {
                _histos["dphi_K_4_5"]->fill(dphi);
                _histos["deta_K_4_5"]->fill(deta);
              }
              else if (particle.pt() > 5 && particle.pt() < 6)
              {
                _histos["dphi_K_5_6"]->fill(dphi);
                _histos["deta_K_5_6"]->fill(deta);
              }
              else if (particle.pt() > 6 && particle.pt() < 10)
              {
                _histos["dphi_K_6_10"]->fill(dphi);
                _histos["deta_K_6_10"]->fill(deta);
              }
            }
           
        }
        
      }
      if(jets.size() == 0) vetoEvent;

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      float numJets = _histos["number_of_jets"]->bin(0).numEntries();
      // why is number of jets = 0?
      std::cout << "Number of jets: " << numJets << std::endl;
      // extend for all histograms
      /*_histos["dphi_pi"]->scaleW(1/(numJets));
      _histos["deta_pi"]->scaleW(1/(numJets));
      _histos["dphi_p"]->scaleW(1/(numJets));
      _histos["deta_p"]->scaleW(1/(numJets));
      _histos["dphi_k"]->scaleW(1/(numJets));
      _histos["deta_k"]->scaleW(1/(numJets));*/
      // determine the minimum value of the dphi histograms and subtract that value from all bins
// Loop over all histograms for dphi, find minimum value for dphi
    for (const auto& hist : _histos) {
        if (hist.first.find("dphi_") != std::string::npos) { // Check if histogram is for dphi_pi
            // Find the minimum value in the histogram bins
            double min_value = hist.second->bin(0).height(); // Initialize with the height of the first bin
            for (size_t i = 1; i < hist.second->numBins(); ++i) {
                double bin_value = hist.second->bin(i).height();
                if (bin_value < min_value) {
                    min_value = bin_value;
                }
    //subract that value from all bins 
            
          
       }
    }
  }
      // Integrate the dphi histograms from -pi/2 to pi/2 and from pi/2 to 3pi/2, for the near-side and away-side respectively for each momentum bin and particle type and fill yield histograms

      // divide the yield histograms for K/pi and p/pi on the near-side and away-side


      //divide(_histos["dphi_k"], _histos["dphi_pi"], _scatters["dphi_ktopi"]);
// proton normalization
      _histos["dphi_p_1_1.5"]->scaleW(1/(numJets));
      _histos["deta_p_1_1.5"]->scaleW(1/(numJets));
      _histos["dphi_p_1.5_2"]->scaleW(1/(numJets));
      _histos["deta_p_1.5_2"]->scaleW(1/(numJets));
      _histos["dphi_p_2_3"]->scaleW(1/(numJets));
      _histos["deta_p_2_3"]->scaleW(1/(numJets));
      _histos["dphi_p_3_4"]->scaleW(1/(numJets));
      _histos["deta_p_3_4"]->scaleW(1/(numJets));
      _histos["dphi_p_4_5"]->scaleW(1/(numJets));
      _histos["deta_p_4_5"]->scaleW(1/(numJets));
      _histos["dphi_p_5_6"]->scaleW(1/(numJets));
      _histos["deta_p_5_6"]->scaleW(1/(numJets));
      _histos["dphi_p_6_10"]->scaleW(1/(numJets));
      _histos["deta_p_6_10"]->scaleW(1/(numJets));

// pion normalization
      _histos["dphi_pi_1_1.5"]->scaleW(1/(numJets));
      _histos["deta_pi_1_1.5"]->scaleW(1/(numJets));
      _histos["dphi_pi_1.5_2"]->scaleW(1/(numJets));
      _histos["deta_pi_1.5_2"]->scaleW(1/(numJets));
      _histos["dphi_pi_2_3"]->scaleW(1/(numJets));
      _histos["deta_pi_2_3"]->scaleW(1/(numJets));
      _histos["dphi_pi_3_4"]->scaleW(1/(numJets));
      _histos["deta_pi_3_4"]->scaleW(1/(numJets));
      _histos["dphi_pi_4_5"]->scaleW(1/(numJets));
      _histos["deta_pi_4_5"]->scaleW(1/(numJets));
      _histos["dphi_pi_5_6"]->scaleW(1/(numJets));
      _histos["deta_pi_5_6"]->scaleW(1/(numJets));
      _histos["dphi_pi_6_10"]->scaleW(1/(numJets));
      _histos["deta_pi_6_10"]->scaleW(1/(numJets));

// kaon normalization
      _histos["dphi_K_1_1.5"]->scaleW(1/(numJets));
      _histos["deta_K_1_1.5"]->scaleW(1/(numJets));
      _histos["dphi_K_1.5_2"]->scaleW(1/(numJets));
      _histos["deta_K_1.5_2"]->scaleW(1/(numJets));
      _histos["dphi_K_2_3"]->scaleW(1/(numJets));
      _histos["deta_K_2_3"]->scaleW(1/(numJets));
      _histos["dphi_K_3_4"]->scaleW(1/(numJets));
      _histos["deta_K_3_4"]->scaleW(1/(numJets));
      _histos["dphi_K_4_5"]->scaleW(1/(numJets));
      _histos["deta_K_4_5"]->scaleW(1/(numJets));
      _histos["dphi_K_5_6"]->scaleW(1/(numJets));
      _histos["deta_K_5_6"]->scaleW(1/(numJets));
      _histos["dphi_K_6_10"]->scaleW(1/(numJets));
      _histos["deta_K_6_10"]->scaleW(1/(numJets));
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _histos;
    map<string, Profile1DPtr> _profiles;
    map<string, CounterPtr> _counters;
    map<string, Scatter2DPtr> _scatters;
    /// @}
    fastjet::AreaDefinition *fjAreaDef;


  };


  RIVET_DECLARE_PLUGIN(ALICE_2024_IPAT);

}
