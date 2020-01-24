// -*- C++ -*-
#include "Rivet/Analysis.hh"
//#include "Rivet/Particle.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
//#include "Rivet/Projections/LeadingParticlesFinalState"
#include "Rivet/Math/MathUtils.hh"
//#include "Rivet/HeavyIonAnalysis.hh"
//#include "Rivet/Projections/ALICEToolsHI.hh"

#include "iostream"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2019_JETV2 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2019_JETV2);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      // the basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs((Cuts::pT > 150*MeV) && (Cuts::abseta < 0.9));
      declare (fs, "FS");
      // charged final state partcles
      ChargedFinalState cfs((Cuts::pT > 150*MeV) && (Cuts::abseta < 0.9));
      declare(cfs, "tracks");
      // the final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetcfs(cfs, FastJets::ANTIKT, 0.2, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetcfs, "jets");

      // Book histograms
      book(_h["num_part"], "N_Particles", 500, 0, 500);
      book(_h["num_jet"], "N_Jets", 50, 0, 50);
      book(_h["part_pt"], "Particle_Momentum_Spectra", 150, 0, 150);
      book(_h["jet_pt"], "Jet_Momentum_Spectra", 15, 0, 150);
      book(_h["part_phi"], "Particle_phi", 100, 0, 2*pi);
      book(_h["jet_phi"], "Jet_phi", 100, 0, 2*pi);
      book(_h["jet_phi_00"], "Jet_phi_00", 100, 0, 2*pi);
      book(_h["jet_phi_01"], "Jet_phi_01", 100, 0, 2*pi);
      book(_h["jet_phi_02"], "Jet_phi_02", 100, 0, 2*pi);
      book(_h["jet_phi_03"], "Jet_phi_03", 100, 0, 2*pi);
      book(_h["jet_phi_04"], "Jet_phi_04", 100, 0, 2*pi);
      book(_h["jet_phi_05"], "Jet_phi_05", 100, 0, 2*pi);
      book(_h["jet_phi_06"], "Jet_phi_06", 100, 0, 2*pi);
      book(_h["jet_phi_07"], "Jet_phi_07", 100, 0, 2*pi);
      book(_h["jet_phi_08"], "Jet_phi_08", 100, 0, 2*pi);
      book(_h["jet_phi_09"], "Jet_phi_09", 100, 0, 2*pi);
      book(_h["jet_phi_10"], "Jet_phi_10", 100, 0, 2*pi);
      book(_h["jet_phi_11"], "Jet_phi_11", 100, 0, 2*pi);
      book(_h["jet_phi_12"], "Jet_phi_12", 100, 0, 2*pi);
      book(_h["jet_phi_13"], "Jet_phi_13", 100, 0, 2*pi);
      book(_h["jet_phi_14"], "Jet_phi_14", 100, 0, 2*pi);
      book(_h["part_eta"], "Particle_eta", 100, -3, 3);
      book(_h["jet_eta"], "Jet_eta", 100, -3, 3);
      book(_h["ep"], "Event_plane_angle", 100, -1*pi, pi);
      //book(_h["ep1"], "Event_plane_angle_1", 360, 0, 2*pi);
      //book(_h["ep2"], "Event_plane_angle_2", 360, 0, 2*pi);
      //book(_h["ep3"], "Event_plane_angle_3", 360, 0, 2*pi);
      book(_h["part_in"], "In_plane_particles", 20, 0, 100);
      book(_h["part_out"], "Out_of_plane_particles", 20, 0, 100);
      book(_h["jet_in"], "In_plane_jets", 15, 0, 150);
      book(_h["jet_out"], "Out_of_plane_jets", 15, 0, 150);

      // Booking Scatters
      book(_h_jet_ratio, "Jet_ratio");
      book(_h_jet_v2, "Jet_v2");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Using jewel? Set this to true to force EP angle = 0
      bool zero_ep = 1;

      // centrality check
      //const double cent = centrality(event, "ImpactParameterMethod");
      //if((c < 0.) || (c > 100.)) vetoEvent;

      // particles - calculating q vectors
      double Qx = 0, Qy = 0;
      int num_part = 0;
      const Particles& trks = apply<ChargedFinalState>(event, "tracks").particles();
      for (const Particle& p : trks) {
        _h["part_pt"]->fill(p.pT());
        _h["part_phi"]->fill(p.phi());
        _h["part_eta"]->fill(p.eta());
        Qx += cos(2*p.phi());
        Qy += sin(2*p.phi());
        num_part++;
      }
      _h["num_part"]->fill(num_part);

      // event plane
      double ep_angle = 0.0;
      if (!zero_ep) {
        if (Qx == 0) vetoEvent;
        else {
          //ep_angle = 0.5*atan(Qy/Qx);
          ep_angle = mapAngle0ToPi(0.5*atan(Qy/Qx));
        }
      }
      _h["ep"]->fill(ep_angle);

      // retrieve clustered jets, sorted by pT, with a minimum pT cut
      int num_jet = 0;
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 10*GeV);
      for(const Jet& j : jets) {
        Particles constituents = j.particles();
        if (constituents[0].pT() >= 3*GeV) {
          _h["jet_pt"]->fill(j.pT());
          _h["jet_phi"]->fill(j.phi());
          _h["jet_eta"]->fill(j.eta());
          if (inRange(j.pT(), 0, 10)) _h["jet_phi_00"]->fill(j.phi());
          if (inRange(j.pT(), 10, 20)) _h["jet_phi_01"]->fill(j.phi());
          if (inRange(j.pT(), 20, 30)) _h["jet_phi_02"]->fill(j.phi());
          if (inRange(j.pT(), 30, 40)) _h["jet_phi_03"]->fill(j.phi());
          if (inRange(j.pT(), 40, 50)) _h["jet_phi_04"]->fill(j.phi());
          if (inRange(j.pT(), 50, 60)) _h["jet_phi_05"]->fill(j.phi());
          if (inRange(j.pT(), 60, 70)) _h["jet_phi_06"]->fill(j.phi());
          if (inRange(j.pT(), 70, 80)) _h["jet_phi_07"]->fill(j.phi());
          if (inRange(j.pT(), 80, 90)) _h["jet_phi_08"]->fill(j.phi());
          if (inRange(j.pT(), 90, 100)) _h["jet_phi_09"]->fill(j.phi());
          if (inRange(j.pT(), 100, 110)) _h["jet_phi_10"]->fill(j.phi());
          if (inRange(j.pT(), 110, 120)) _h["jet_phi_11"]->fill(j.phi());
          if (inRange(j.pT(), 120, 130)) _h["jet_phi_12"]->fill(j.phi());
          if (inRange(j.pT(), 130, 140)) _h["jet_phi_13"]->fill(j.phi());
          if (inRange(j.pT(), 140, 150)) _h["jet_phi_14"]->fill(j.phi());

          if (!zero_ep) {
            // counting in-plane vs out-of-plane jets
            // 1st case if EP angle is between 0 and pi/4
            if (inRange(ep_angle, 0, pi/4)) {
              // angular ranges for in-plane jets
              if (inRange(j.phi(), ep_angle - pi/4, ep_angle + pi/4) ||
                  inRange(j.phi(), ep_angle + 3*pi/4, ep_angle + 5*pi/4) ||
                  inRange(j.phi(), ep_angle + 7*pi/4, ep_angle + 9*pi/4)) {
                _h["jet_in"]->fill(j.pT());
                _h["jet_dif"]->fill(j.pT());
              }
              // angular ranges for out-of-plane jets
              if (inRange(j.phi(), ep_angle + pi/4, ep_angle + 3*pi/4) ||
                  inRange(j.phi(), ep_angle + 5*pi/4, ep_angle + 7*pi/4)) {
                _h["jet_out"]->fill(j.pT());
              }
            }
            // 2nd case if EP angle is > pi/4
            if (inRange(ep_angle, pi/4, pi/2)) {
              //angular ranges for in-plane jets
              if (inRange(j.phi(), ep_angle - pi/4, ep_angle + pi/4) ||
                  inRange(j.phi(), ep_angle + 3*pi/4, ep_angle + 5*pi/4) ||
                  inRange(j.phi(), ep_angle + 7*pi/4, 2*pi)) {
                _h["jet_in"]->fill(j.pT());
              }
              // angular ranges for out-of-plane jets
              if (inRange(j.phi(), 0, ep_angle - pi/4) ||
                  inRange(j.phi(), ep_angle + pi/4, ep_angle + 3*pi/4) ||
                  inRange(j.phi(), ep_angle + 5*pi/4, ep_angle + 7*pi/4)) {
                _h["jet_out"]->fill(j.pT());
              }
            }
          }

          if (zero_ep) {
            // counting in-plane vs out-of-plane jets
            // angular ranges for in-plane jets
            if (inRange(j.phi(), 0, pi/4) ||
                inRange(j.phi(), 3*pi/4, 5*pi/4) ||
                inRange(j.phi(), 7*pi/4, 2*pi)) {
              _h["jet_in"]->fill(j.pT());
            }
            // angular ranges for out-of-plane jets
            if (inRange(j.phi(), pi/4, 3*pi/4) ||
                inRange(j.phi(), 5*pi/4, 7*pi/4)) {
              _h["jet_out"]->fill(j.pT());
            }
          }
          // counting jets per event
          num_jet++;
        }
        _h["num_jet"]->fill(num_jet);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // calculating jet v2
      int n_bins = _h["jet_in"]->numBins();
      cout << "n_bins = " << n_bins << endl;
      for (int i = 0; i < n_bins; i++) {
        double n_in = 0, n_out = 0;
        double jet_v2 = 0, ep_res = 1.0, jet_rat = 0.0;
        cout << "Bin number " << i << endl;
        n_in = _h["jet_in"]->bin(i).sumW();
        n_out = _h["jet_out"]->bin(i).sumW();
        cout << "n_in = " << n_in << "; n_out = " << n_out << endl;
        if (n_out == 0) jet_rat = 0.0;
        else jet_rat = (double) n_in / (double) n_out;
        cout << "Jet ratio (n_in/n_out) = " << jet_rat << endl;
        if (n_in + n_out == 0) jet_v2 = 0.0;
        else jet_v2 = (pi/4)*(1/ep_res)*(n_in - n_out)/(n_in + n_out);
        cout << "jet v_2 = " << jet_v2 << endl;
      }

      divide(_h["jet_in"], _h["jet_out"], _h_jet_ratio);
      divide(*_h["jet_in"] - *_h["jet_out"], *_h["jet_in"] + *_h["jet_out"], _h_jet_v2);
      _h_jet_v2->scale(1.0, pi/4.0);
    }

  private:
    Scatter2DPtr _h_jet_ratio;
    Scatter2DPtr _h_jet_v2;

    //@}

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2019_JETV2);

}
