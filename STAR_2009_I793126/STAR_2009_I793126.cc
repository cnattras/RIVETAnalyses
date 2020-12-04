// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class STAR_2009_I793126 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2009_I793126);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      //const FinalState fs(Cuts::abseta < 0.5 && Cuts::pT > 0.15*GeV && Cuts::abscharge > 0);
      const FinalState fs(Cuts::absrap < 0.1 && Cuts::pT > 0.2*GeV && Cuts::abscharge > 0);
      declare(fs, "fs");

      beamOpt = getOption<string>("beam", "NONE");
      if (beamOpt == "DAU200") collSys = dAu200;
      if (beamOpt == "AUAU62") collSys = AuAu62;
      if (beamOpt == "AUAU130") collSys = AuAu130;

      //Event counters
      book(_c["sow_dAu200"], "sow_dAu200");
      book(_c["sow_AuAu62"], "sow_AuAu62");
      book(_c["sow_AuAu130"], "sow_AuAu130");

      //Booking histograms for figures in paper
      /*//Figure 1 AuAu200
      book(_h["Figure_1_AuAu200"], 1, 1, 1);
      //Figure 1 AuAu62
      book(_h["Figure_1_AuAu62"], 2, 1, 1);
      //Figure 2a
      book(_h["Figure_2a"], 3, 1, 1);
      //Figure 2b
      book(_h["Figure_2b_1"], 4, 1, 1);
      book(_h["Figure_2b_2"], 4, 1, 2);
      book(_h["Figure_2b_3"], 4, 1, 3);
      book(_h["Figure_2b_4"], 4, 1, 4);
      //Figure 3
      book(_h["Figure_3"], 5, 1, 1);
      //Figure 4a AuAu200
      book(_h["Figure_4a_AuAu200_1"], 6, 1, 1);
      book(_h["Figure_4a_AuAu200_2"], 6, 1, 2);
      //Figure 4a AuAu62
      book(_h["Figure_4a_AuAu62_1"], 7, 1, 1);
      book(_h["Figure_4a_AuAu62_2"], 7, 1, 2);
      //Figure 4a pp
      book(_h["Figure_4a_pp"], 8, 1, 1);
      //Figure 4b AuAu200
      book(_h["Figure_4b_AuAu200_1"], 9, 1, 1);
      book(_h["Figure_4b_AuAu200_2"], 9, 1, 2);
      //Figure 4b AuAu62
      book(_h["Figure_4b_AuAu62_1"], 10, 1, 1);
      book(_h["Figure_4b_AuAu62_2"], 10, 1, 2);
      //Figure 4b pp
      book(_h["Figure_4b_pp"], 11, 1, 1);
      //Figure 4b AuAu200
      //Figure 4b AuAu62
      //Figure 10
      book(_h["Figure_10_1"], 12, 1, 1);
      book(_h["Figure_10_2"], 12, 1, 2);
      book(_h["Figure_10_3"], 12, 1, 3);
      book(_h["Figure_10_4"], 12, 1, 4);
      book(_h["Figure_10_5"], 12, 1, 5);
      book(_h["Figure_10_6"], 12, 1, 6);
      //Figure 11a
      book(_h["Figure_11a_1"], 13, 1, 1);
      book(_h["Figure_11a_2"], 13, 1, 2);
      //Figure 11b
      book(_h["Figure_11b_1"], 14, 1, 1);
      book(_h["Figure_11b_2"], 14, 1, 2);
      //Figure 12 dAu
      book(_h2D["Figure_12_dAu"], 15, 1, 1);
      //Figure 12 pp
      book(_h2D["Figure_12_pp"], 16, 1, 1);
      //Figure 13a
      book(_h2D["Figure_13a_1"], 17, 1, 1);
      book(_h2D["Figure_13a_2"], 17, 1, 2);
      book(_h2D["Figure_13a_3"], 17, 1, 3);
      //Figure 13b
      book(_h2D["Figure_13b_1"], 18, 1, 1);
      book(_h2D["Figure_13b_2"], 18, 1, 2);
      book(_h2D["Figure_13b_3"], 18, 1, 3);
      //Figure 14a pion and kaon
      book(_h2D["Figure_14a_pion_and_kaon_1"], 19, 1, 1);
      book(_h2D["Figure_14a_pion_and_kaon_2"], 19, 1, 2);
      //Figure 14a proton
      book(_h2D["Figure_14a_proton"], 20, 1, 1);
      //Figure 14b pion and kaon
      book(_h2D["Figure_14b_pion_and_kaon_1"], 21, 1, 1);
      book(_h2D["Figure_14b_pion_and_kaon_2"], 21, 1, 2);
      //Figure 14b proton
      book(_h2D["Figure_14b_proton"], 22, 1, 1);
      //Figure 15a
      book(_h["Figure_15a_1"], 23, 1, 1);
      book(_h["Figure_15a_2"], 23, 1, 2);
      book(_h["Figure_15a_3"], 23, 1, 3);
      book(_h["Figure_15a_4"], 23, 1, 4);
      //Figure 15a Bkgd
      book(_h["Figure_15a_Bkgd"], 24, 1, 1);
      //Figure 15b
      book(_h["Figure_15b_1"], 25, 1, 1);
      book(_h["Figure_15b_2"], 25, 1, 2);
      book(_h["Figure_15b_3"], 25, 1, 3);
      book(_h["Figure_15b_4"], 25, 1, 4);
      //Figure 15b Bkgd
      book(_h["Figure_15b_Bkgd"], 26, 1, 1);
      //Figure 15c
      book(_h["Figure_15c_1"], 27, 1, 1);
      book(_h["Figure_15c_2"], 27, 1, 2);
      book(_h["Figure_15c_3"], 27, 1, 3);
      book(_h["Figure_15c_4"], 27, 1, 4);
      //Figure 15c Bkgd
      book(_h["Figure_15c_Bkgd"], 28, 1, 1);
      //Figure 15d
      book(_h["Figure_15d_1"], 29, 1, 1);
      book(_h["Figure_15d_2"], 29, 1, 2);
      book(_h["Figure_15d_3"], 29, 1, 3);
      book(_h["Figure_15d_4"], 29, 1, 4);
      //Figure 15d Bkgd
      book(_h["Figure_15d_Bkgd"], 30, 1, 1);
      //Figure 16 all and weak-decay bkgd
      book(_h2D["Figure_16_allandweakdecay_bkgd_1"], 31, 1, 1);
      book(_h2D["Figure_16_allandweakdecay_bkgd_2"], 31, 1, 2);
      //Figure 16 muon contamination
      book(_h2D["Figure_16_muon_contamination"], 32, 1, 1);
      //Figure 17 AuAu dE/dx
      book(_h2D["Figure17_AuAu_dEdx"], 33, 1, 1);
      //Figure 17 AuAu Blast-wave fit
      book(_h["Figure_17_AuAu_Blastwave_fit"], 34, 1, 1);
      //Figure 17 AuAu p_T-Gaussian fit
      book(_h["Figure_17_AuAu_pT_Gaussian_fit"], 35, 1, 1);
      //Figure 17 AuAu p_T-exponential fit
      book(_h["Figure_17_AuAu_pT_exponential_fit"], 36, 1, 1);
      //Figure 17 AuAu TOF data
      book(_h2D["Figure_17_AuAu_TOF_data"], 37, 1, 1);
      //Figure 17 dAu dE/dx
      book(_h2D["Figure_17_dAu_dEdx"], 38, 1, 1);
      //Figure 17 dAu Blast-wave fit
      book(_h["Figure_17_dAu_pT_Blastwave_fit"], 39, 1, 1);
      //Figure 17 dAu p_T-Gaussian fit
      book(_h["Figure_17_dAu_pT_Gaussian_fit"], 40, 1, 1);
      //Figure 17 dAu p_T-exponential fit
      book(_h["Figure_17_dAu_pT_exponential_fit"], 41, 1, 1);
      //Figure 17 dAu TOF data
      book(_h2D["Figure_17_dAu_TOF_data"], 42, 1, 1);*/
      //Figure 18 kaon
      book(_h["Figure_18_kaon_1"], 43, 1, 1);
      book(_h["Figure_18_kaon_2"], 43, 1, 2);
      book(_h["Figure_18_kaon_3"], 43, 1, 3);
      book(_h["Figure_18_kaon_4"], 43, 1, 4);
      book(_h["Figure_18_kaon_5"], 43, 1, 5);
      book(_h["Figure_18_kaon_6"], 43, 1, 6);
      book(_h["Figure_18_kaon_7"], 43, 1, 7);
      book(_h["Figure_18_kaon_8"], 43, 1, 8);
      //Figure 18 pion
      book(_h["Figure_18_pion_1"], 44, 1, 1);
      book(_h["Figure_18_pion_2"], 44, 1, 2);
      book(_h["Figure_18_pion_3"], 44, 1, 3);
      book(_h["Figure_18_pion_4"], 44, 1, 4);
      book(_h["Figure_18_pion_5"], 44, 1, 5);
      book(_h["Figure_18_pion_6"], 44, 1, 6);
      book(_h["Figure_18_pion_7"], 44, 1, 7);
      book(_h["Figure_18_pion_8"], 44, 1, 8);
      //Figure 18 proton
      book(_h["Figure_18_proton_1"], 45, 1, 1);
      book(_h["Figure_18_proton_2"], 45, 1, 2);
      book(_h["Figure_18_proton_3"], 45, 1, 3);
      book(_h["Figure_18_proton_4"], 45, 1, 4);
      book(_h["Figure_18_proton_5"], 45, 1, 5);
      book(_h["Figure_18_proton_6"], 45, 1, 6);
      book(_h["Figure_18_proton_7"], 45, 1, 7);
      book(_h["Figure_18_proton_8"], 45, 1, 8);
      //Figure 19 kaon
      book(_h["Figure_19_kaon_1"], 46, 1, 1);
      book(_h["Figure_19_kaon_2"], 46, 1, 2);
      book(_h["Figure_19_kaon_3"], 46, 1, 3);
      book(_h["Figure_19_kaon_4"], 46, 1, 4);
      book(_h["Figure_19_kaon_5"], 46, 1, 5);
      book(_h["Figure_19_kaon_6"], 46, 1, 6);
      book(_h["Figure_19_kaon_7"], 46, 1, 7);
      book(_h["Figure_19_kaon_8"], 46, 1, 8);
      book(_h["Figure_19_kaon_9"], 46, 1, 9);
      book(_h["Figure_19_kaon_10"], 46, 1, 10);
      book(_h["Figure_19_kaon_11"], 46, 1, 11);
      book(_h["Figure_19_kaon_12"], 46, 1, 12);
      book(_h["Figure_19_kaon_13"], 46, 1, 13);
      book(_h["Figure_19_kaon_14"], 46, 1, 14);
      book(_h["Figure_19_kaon_15"], 46, 1, 15);
      book(_h["Figure_19_kaon_16"], 46, 1, 16);
      book(_h["Figure_19_kaon_17"], 46, 1, 17);
      book(_h["Figure_19_kaon_18"], 46, 1, 18);
      //Figure 19 pion
      book(_h["Figure_19_pion_1"], 47, 1, 1);
      book(_h["Figure_19_pion_2"], 47, 1, 2);
      book(_h["Figure_19_pion_3"], 47, 1, 3);
      book(_h["Figure_19_pion_4"], 47, 1, 4);
      book(_h["Figure_19_pion_5"], 47, 1, 5);
      book(_h["Figure_19_pion_6"], 47, 1, 6);
      book(_h["Figure_19_pion_7"], 47, 1, 7);
      book(_h["Figure_19_pion_8"], 47, 1, 8);
      book(_h["Figure_19_pion_9"], 47, 1, 9);
      book(_h["Figure_19_pion_10"], 47, 1, 10);
      book(_h["Figure_19_pion_11"], 47, 1, 11);
      book(_h["Figure_19_pion_12"], 47, 1, 12);
      book(_h["Figure_19_pion_13"], 47, 1, 13);
      book(_h["Figure_19_pion_14"], 47, 1, 14);
      book(_h["Figure_19_pion_15"], 47, 1, 15);
      book(_h["Figure_19_pion_16"], 47, 1, 16);
      book(_h["Figure_19_pion_17"], 47, 1, 17);
      book(_h["Figure_19_pion_18"], 47, 1, 18);
      //Figure 19 proton
      book(_h["Figure_19_proton_1"], 48, 1, 1);
      book(_h["Figure_19_proton_2"], 48, 1, 2);
      book(_h["Figure_19_proton_3"], 48, 1, 3);
      book(_h["Figure_19_proton_4"], 48, 1, 4);
      book(_h["Figure_19_proton_5"], 48, 1, 5);
      book(_h["Figure_19_proton_6"], 48, 1, 6);
      book(_h["Figure_19_proton_7"], 48, 1, 7);
      book(_h["Figure_19_proton_8"], 48, 1, 8);
      book(_h["Figure_19_proton_9"], 48, 1, 9);
      book(_h["Figure_19_proton_10"], 48, 1, 10);
      book(_h["Figure_19_proton_11"], 48, 1, 11);
      book(_h["Figure_19_proton_12"], 48, 1, 12);
      book(_h["Figure_19_proton_13"], 48, 1, 13);
      book(_h["Figure_19_proton_14"], 48, 1, 14);
      book(_h["Figure_19_proton_15"], 48, 1, 15);
      book(_h["Figure_19_proton_16"], 48, 1, 16);
      book(_h["Figure_19_proton_17"], 48, 1, 17);
      book(_h["Figure_19_proton_18"], 48, 1, 18);
      //Figure 20
      book(_h["Figure_20_1"], 49, 1, 1);
      book(_h["Figure_20_2"], 49, 1, 2);
      book(_h["Figure_20_3"], 49, 1, 3);
      book(_h["Figure_20_4"], 49, 1, 4);
      book(_h["Figure_20_5"], 49, 1, 5);
      book(_h["Figure_20_6"], 49, 1, 6);
      book(_h["Figure_20_7"], 49, 1, 7);
      book(_h["Figure_20_8"], 49, 1, 8);
      book(_h["Figure_20_9"], 49, 1, 9);
      book(_h["Figure_20_10"], 49, 1, 10);
      book(_h["Figure_20_11"], 49, 1, 11);
      book(_h["Figure_20_12"], 49, 1, 12);
      book(_h["Figure_20_13"], 49, 1, 13);
      book(_h["Figure_20_14"], 49, 1, 14);
      book(_h["Figure_20_15"], 49, 1, 15);
      book(_h["Figure_20_16"], 49, 1, 16);
      /*//Figure 24 Au+Au 62.4 GeV
      book(_h["Figure_24_AuAu_62_point4_GeV_1"], 50, 1, 1);
      book(_h["Figure_24_AuAu_62_point4_GeV_2"], 50, 1, 2);
      //Figure 24 Au+Au 130 GeV
      book(_h["Figure_24_AuAu_130_GeV_1"], 51, 1, 1);
      book(_h["Figure_24_AuAu_130_GeV_2"], 51, 1, 2);
      //Figure 24 Au+Au 200 GeV
      book(_h["Figure_24_AuAu_200_GeV_1"], 52, 1, 1);
      book(_h["Figure_24_AuAu_200_GeV_2"], 52, 1, 2);
      //Figure 30a
      book(_h["Figure_30a"], 53, 1, 1);
      //Figure 30b
      book(_h["Figure_30b"], 54, 1, 1);
      //Figure 32 AGS E859 Si+Al 5.4 GeV
      book(_h["Figure_32_AGS_E859_SiAl_5_point4_GeV"], 55, 1, 1);
      //Figure 32 AGS E866 Au+Au 4.7 GeV
      book(_h["Figure_32_AGS_AuAu_5_point7_GeV"], 56, 1, 1);
      //Figure 32 SPS NA49 Pb+Pb 17.3 GeV
      book(_h["Figure_32_SPS_NA49_PbPb_17_point3_GeV"], 57, 1, 1);
      //Figure 32 SPS NA49 Pb+Pb energy scan
      book(_h["Figure_32_SPS_NA49_PbPb_beam_energy_scan"], 58, 1, 1);
      //Figure 32 SPS NA49 S+S 20 GeV
      book(_h["Figure_32_SPS_NA49_SS_20_GeV"], 59, 1, 1);
      //Figure 32 SPS NA49 C+C/Si+Si 17.3 GeV
      book(_h["Figure_32_SPS_NA49_CC_SiSi_17_point3_GeVv"], 60, 1, 1);
      //Figure 32 STAR Au+Au 62.4 GeV
      book(_h["Figure_32_AuAu_62_point4_GeV"], 61, 1, 1);
      //Figure 32 STAR Au+Au 130 GeV
      book(_h["Figure_32_STAR_AuAu_130_GeV"], 62, 1, 1);
      //Figure 32 STAR Au-Au 200 GeV
      book(_h["Figure_32_AuAu_200_GeV"], 63, 1, 1);
      //Figure 33 AGS E859 Si+Al 5.4 GeV
      book(_h["Figure_33_AGS_E859_SiAl_5_point4_GeV"], 64, 1, 1);
      //Figure 33 AGS E866 Au+Au 4.7 GeV
      book(_h["Figure_33_AGS_E866_AuAu_4_point7_GeV"], 65, 1, 1);
      //Figure 33 SPS NA49 Pb+Pb 17.3 GeV
      book(_h["Figure_33_SPS_NA49_PbPb_17_point3_GeV"], 66, 1, 1);
      //Figure 33 SPS NA49 Pb+Pb energy scan
      book(_h["Figure_33_SPS_NA49_PbPb_energy_scan"], 67, 1, 1);
      //Figure 33 SPS NA49 S+S 20 GeV
      book(_h["Figure_33_SPS_NA49_SS_20_GeV"], 68, 1, 1);
      //Figure 33 SPS NA49 C+C/Si+Si 17.3 GeV
      book(_h["Figure_33_SPS_NA49_CC_SiSi_17_point3_GeV"], 69, 1, 1);
      //Figure 33 STAR Au+Au 62.4 GeV
      book(_h["Figure_33_STAR_AuAu_62_point4_GeV"], 70, 1, 1);
      //Figure 33 STAR Au+Au 130 GeV
      book(_h["Figure_33_STAR_AuAu_130_GeV"], 71, 1, 1);
      //Figure 33 STAR Au+Au 200 GeV
      book(_h["Figure_33_STAR_AuAu_200_GeV"], 72, 1, 1);
      //Figure 38 Becattini et al.
      book(_h["Figure_38_Becattini_et_al"], 73, 1, 1);
      //Figure 38 Andronic et al.
      book(_h["Figure_38_Andronic_et_al"], 74, 1, 1);
      //Figure 38 SIS
      book(_h["Figure_38_SIS"], 75, 1, 1);
      //Figure 38 AGS
      book(_h["Figure_38_AGS"], 76, 1, 1);
      //Figure 38 SPS
      book(_h["Figure_38_SPS"], 77, 1, 1);
      //Figure 38 STAR
      book(_h["Figure_38_STAR"], 78, 1, 1);
      //Figure 39 Elementary collisions Becattini et al.
      book(_h["Figure_39_Elementary_collisions_Becattini_et_al"], 79, 1, 1);
      //Figure 39 Becattini et al.
      book(_h["Figure_39_Becattini_et_al"], 80, 1, 1);
      //Figure 39 Andronic et al.
      book(_h["Figure_39_Andronic_et_al"], 81, 1, 1);
      //Figure 39 SIS
      book(_h["Figure_39_SIS"], 82, 1, 1);
      //Figure 39 AGS
      book(_h["Figure_39_AGS"], 83, 1, 1);
      //Figure 39 SPS
      book(_h["Figure_39_SPS"], 84, 1, 1);
      //Figure 39 STAR pp
      book(_h["Figure_39_STAR_pp"], 85, 1, 1);
      //Figure 39 STAR Tchem
      book(_h["Figure_39_STAR_Tchem"], 86, 1, 1);
      //Figure 39 EOS
      book(_h["Figure_39_EOS"], 87, 1, 1);
      //Figure 39 FOPI
      book(_h["Figure_39_FOPI"], 88, 1, 1);
      //Figure 39 E866
      book(_h["Figure_39_E866"], 89, 1, 1);
      //Figure 39 NA49
      book(_h["Figure_39_NA49"], 90, 1, 1);
      //Figure 39 STAR Tkin
      book(_h["Figure_39_STAR_Tkin"], 91, 1, 1);
      //Figure 40 FOPI
      book(_h["Figure_40_FOPI"], 92, 1, 1);
      //Figure 40 EOS
      book(_h["Figure_40_EOS"], 93, 1, 1);
      //Figure 40 E866
      book(_h["Figure_40_E866"], 94, 1, 1);
      //Figure 40 NA49
      book(_h["Figure_40_NA49"], 95, 1, 1);
      //Figure 40 STAR
      book(_h["Figure_40_STAR"], 96, 1, 1);
      //Figure 41 mu Andronic et al.
      //book(_h["Figure_41_mu_Andronic_et_al"], 97, 1, 1);
      //Figure 41 mu SIS
      //book(_h["Figure_41_mu_SIS"], 98, 1, 1);
      //Figure 41 mu AGS 4.8 GeV
      //book(_h["Figure_41_mu_AGS_4_point8_GeV"], 99, 1, 1);
      //Figure 41 mu SPS
      //book(_h["Figure_41_mu_SPS"], 100, 1, 1);
      //Figure 41 mu STAR Au+Au
      //book(_h["Figure_41_mu_STAR_AuAu"], 101, 1, 1);
      //Figure 41 mu STAR pp 200 GeV
      //book(_h["Figure_41_mu_STAR_pp_200_GeV"], 102, 1, 1);
      //Figure 41 Tch Andronic et al.
      //book(_h["Figure_41_Tch_Andronic_et_al"], 103, 1, 1);
      //Figure 41 Tch SIS
      //book(_h["Figure_41_Tch_SIS"], 104, 1, 1);
      //Figure 41 Tch AGS 4.8 GeV
      //book(_h2D["Figure_41_Tch_AGS_4_point8_GeV"], 105, 1, 1);
      //Figure 41 Tch SPS
      //book(_h["Figure_41_Tch_SPS"], 106, 1, 1);
      //Figure 41 Tch STAR Au+Au
      //book(_h["Figure_41_Tch_STAR_Au+Au"], 107, 1, 1);
      //Figure 41 Tch STAR pp 200 GeV
      //book(_h["Figure_41_Tch_STAR_pp_200_GeV"], 108, 1, 1);
      //Figure 42 Optical Glauber
      book(_h["Figure_42_Optical_Glauber"], 109, 1, 1);
      //Figure 42 MC Glauber
      book(_h["Figure_42_MC_Glauber"], 110, 1, 1);*/

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();

      Particles fsParticles = applyProjection<FinalState>(event, "fs").particles();

      //Sorting events by energy
      //Fig 18 dAu spectra @ 200 Gev
      //Sorting particles by type
      if (collSys == dAu200)
      {
	       _c["sow_dAu200"]->fill();

        for (const Particle& p : fsParticles) {
          if (p.pid() == 321) { //kaon+ (KPLUS Pdgid = 321)
            if (c < 20.0) _h["Figure_18_kaon_1"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 0-20% centrality
            if (c > 20.0 && c < 40.0) _h["Figure_18_kaon_2"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-40% centrality
            _h["Figure_18_kaon_3"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // Min Bias
            if (c > 40.0) _h["Figure_18_kaon_4"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-100% centrality
          }
          if (p.pid() == -321) { //kaon- (KMINUS Pdgid = -321)
            if (c < 20.0) _h["Figure_18_kaon_5"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV));
            if (c > 20.0 && c < 40.0) _h["Figure_18_kaon_6"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-40% centrality
            _h["Figure_18_kaon_7"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // Min Bias
            if (c > 40.0) _h["Figure_18_kaon_8"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-100% centrality
          }
          if (p.pid() == 211 && absValue(p.parents()[0]) != 3122 && absValue(p.parents()[0]) != 310 && absValue(p.parents()[0]) != 13) { //pion+ (PIPLUS Pdgid = 211)
            if (c < 20.0) _h["Figure_18_pion_1"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 0-20% centrality
            if (c > 20.0 && c < 40.0) _h["Figure_18_pion_2"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-40% centrality
            _h["Figure_18_pion_3"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // Min Bias
            if (c > 40.0) _h["Figure_18_pion_4"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-100% centrality
          }
          if (p.pid() == -211 && absValue(p.parents()[0]) != 3122 && absValue(p.parents()[0]) != 310 && absValue(p.parents()[0]) != 13) { //pion- (PIMINUS Pdgis = -211)
            if (c < 20.0) _h["Figure_18_pion_5"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV));
            if (c > 20.0 && c < 40.0) _h["Figure_18_pion_6"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-40% centrality
            _h["Figure_18_pion_7"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // Min Bias
            if (c > 40.0) _h["Figure_18_pion_8"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-100% centrality
          }
          if (p.pid() == 2212) { //proton+ (PROTON Pdgid = 2212)
            if (c < 20.0) _h["Figure_18_proton_1"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 0-20% centrality
            if (c > 20.0 && c < 40.0) _h["Figure_18_proton_2"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-40% centrality
            _h["Figure_18_proton_3"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // Min Bias
            if (c > 40.0) _h["Figure_18_proton_4"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-100% centrality
          }
          if (p.pid() == -2212) { //proton- (ANTIPROTON Pdgid = -2212)
            if (c < 20.0) _h["Figure_18_proton_5"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV));
            if (c > 20.0 && c < 40.0) _h["Figure_18_proton_6"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-40% centrality
            _h["Figure_18_proton_7"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // Min Bias
            if (c > 40.0) _h["Figure_18_proton_8"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-100% centrality
          }
        }
      }

      //Figure 19 AuAu @ 62.4 GeV
      if (collSys == AuAu62)
      {
        _c["sow_AuAu62"]->fill();

        for (const Particle& p : fsParticles) {
          if (p.pid() == 321) { //kaon+ (KPLUS Pdgid = 321)
            if (c < 5.0) _h["Figure_19_kaon_1"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 0-5% centrality
            if (c > 5.0 && c < 10.0) _h["Figure_19_kaon_2"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 5-10% centrality
            if (c > 10.0 && c < 20.0) _h["Figure_19_kaon_3"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 10-20% centrality
            if (c > 20.0 && c < 30.0) _h["Figure_19_kaon_4"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-30% centrality
            if (c > 30.0 && c < 40.0) _h["Figure_19_kaon_5"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 30-40% centrality
            if (c > 40.0 && c < 50.0) _h["Figure_19_kaon_6"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-50% centrality
            if (c > 50.0 && c < 60.0) _h["Figure_19_kaon_7"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 50-60% centrality
            if (c > 60.0 && c < 70.0) _h["Figure_19_kaon_8"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 60-70% centrality
            if (c > 70.0 && c < 80.0) _h["Figure_19_kaon_9"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 70-80% centrality
          }
          if (p.pid() == -321) { //kaon- (KMINUS Pdgid = -321)
            if (c < 5.0) _h["Figure_19_kaon_10"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 0-5% centrality
            if (c > 5.0 && c < 10.0) _h["Figure_19_kaon_11"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 5-10% centrality
            if (c > 10.0 && c < 20.0) _h["Figure_19_kaon_12"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 10-20% centrality
            if (c > 20.0 && c < 30.0) _h["Figure_19_kaon_13"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-30% centrality
            if (c > 30.0 && c < 40.0) _h["Figure_19_kaon_14"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 30-40% centrality
            if (c > 40.0 && c < 50.0) _h["Figure_19_kaon_15"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-50% centrality
            if (c > 50.0 && c < 60.0) _h["Figure_19_kaon_16"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 50-60% centrality
            if (c > 60.0 && c < 70.0) _h["Figure_19_kaon_17"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 60-70% centrality
            if (c > 70.0 && c < 80.0) _h["Figure_19_kaon_18"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 70-80% centrality
          }
          if (p.pid() == 211 && absValue(p.parents()[0]) != 3122 && absValue(p.parents()[0]) != 310 && absValue(p.parents()[0]) != 13) { //pion+ (PIPLUS Pdgid = 211)
            if (c < 5.0) _h["Figure_19_pion_1"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 0-5% centrality
            if (c > 5.0 && c < 10.0) _h["Figure_19_pion_2"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 5-10% centrality
            if (c > 10.0 && c < 20.0) _h["Figure_19_pion_3"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 10-20% centrality
            if (c > 20.0 && c < 30.0) _h["Figure_19_pion_4"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-30% centrality
            if (c > 30.0 && c < 40.0) _h["Figure_19_pion_5"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 30-40% centrality
            if (c > 40.0 && c < 50.0) _h["Figure_19_pion_6"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-50% centrality
            if (c > 50.0 && c < 60.0) _h["Figure_19_pion_7"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 50-60% centrality
            if (c > 60.0 && c < 70.0) _h["Figure_19_pion_8"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 60-70% centrality
            if (c > 70.0 && c < 80.0) _h["Figure_19_pion_9"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 70-80% centrality
          }
          if (p.pid() == -211 && absValue(p.parents()[0]) != 3122 && absValue(p.parents()[0]) != 310 && absValue(p.parents()[0]) != 13) { //pion- (PIMINUS Pdgis = -211)
            if (c < 5.0) _h["Figure_19_pion_10"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 0-5% centrality
            if (c > 5.0 && c < 10.0) _h["Figure_19_pion_11"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 5-10% centrality
            if (c > 10.0 && c < 20.0) _h["Figure_19_pion_12"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 10-20% centrality
            if (c > 20.0 && c < 30.0) _h["Figure_19_pion_13"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-30% centrality
            if (c > 30.0 && c < 40.0) _h["Figure_19_pion_14"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 30-40% centrality
            if (c > 40.0 && c < 50.0) _h["Figure_19_pion_15"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-50% centrality
            if (c > 50.0 && c < 60.0) _h["Figure_19_pion_16"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 50-60% centrality
            if (c > 60.0 && c < 70.0) _h["Figure_19_pion_17"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 60-70% centrality
            if (c > 70.0 && c < 80.0) _h["Figure_19_pion_18"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 70-80% centrality
          }
          if (p.pid() == 2212) { //proton+ (PROTON Pdgid = 2212)
            if (c < 5.0) _h["Figure_19_proton_1"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 0-5% centrality
            if (c > 5.0 && c < 10.0) _h["Figure_19_proton_2"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 5-10% centrality
            if (c > 10.0 && c < 20.0) _h["Figure_19_proton_3"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 10-20% centrality
            if (c > 20.0 && c < 30.0) _h["Figure_19_proton_4"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-30% centrality
            if (c > 30.0 && c < 40.0) _h["Figure_19_proton_5"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 30-40% centrality
            if (c > 40.0 && c < 50.0) _h["Figure_19_proton_6"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-50% centrality
            if (c > 50.0 && c < 60.0) _h["Figure_19_proton_7"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 50-60% centrality
            if (c > 60.0 && c < 70.0) _h["Figure_19_proton_8"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 60-70% centrality
            if (c > 70.0 && c < 80.0) _h["Figure_19_proton_9"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 70-80% centrality
          }
          if (p.pid() == -2212) { //proton- (ANTIPROTON Pdgid = -2212)
            if (c < 5.0) _h["Figure_19_proton_10"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 0-5% centrality
            if (c > 5.0 && c < 10.0) _h["Figure_19_proton_11"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 5-10% centrality
            if (c > 10.0 && c < 20.0) _h["Figure_19_proton_12"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 10-20% centrality
            if (c > 20.0 && c < 30.0) _h["Figure_19_proton_13"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 20-30% centrality
            if (c > 30.0 && c < 40.0) _h["Figure_19_proton_14"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 30-40% centrality
            if (c > 40.0 && c < 50.0) _h["Figure_19_proton_15"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 40-50% centrality
            if (c > 50.0 && c < 60.0) _h["Figure_19_proton_16"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 50-60% centrality
            if (c > 60.0 && c < 70.0) _h["Figure_19_proton_17"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 60-70% centrality
            if (c > 70.0 && c < 80.0) _h["Figure_19_proton_18"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 70-80% centrality
          }
        }
      }

      //Figure 20 AuAu @ 130 GeV
      if (collSys == AuAu130)
      {
        _c["sow_AuAu130"]->fill();

        for (const Particle& p : fsParticles) {
          if (p.pid() == 211 && absValue(p.parents()[0]) != 3122 && absValue(p.parents()[0]) != 310 && absValue(p.parents()[0]) != 13) { //pion+ (PIPLUS Pdgid = 211)
            if (c < 6.0) _h["Figure_20_1"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 0-6% centrality
            if (c > 6.0 && c < 11.0) _h["Figure_20_2"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 6-11% centrality
            if (c > 11.0 && c < 18.0) _h["Figure_20_3"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 11-18% centrality
            if (c > 18.0 && c < 26.0) _h["Figure_20_4"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 18-26% centrality
            if (c > 26.0 && c < 34.0) _h["Figure_20_5"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 26-34% centrality
            if (c > 34.0 && c < 45.0) _h["Figure_20_6"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 34-45% centrality
            if (c > 45.0 && c < 58.0) _h["Figure_20_7"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 45-58% centrality
            if (c > 58.0 && c < 85.0) _h["Figure_20_8"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 58-85% centrality
          }
          if (p.pid() == -211 && absValue(p.parents()[0]) != 3122 && absValue(p.parents()[0]) != 310 && absValue(p.parents()[0]) != 13) { //pion- (PIMINUS Pdgis = -211)
            if (c < 6.0) _h["Figure_20_9"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 0-6% centrality
            if (c > 6.0 && c < 11.0) _h["Figure_20_10"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 6-11% centrality
            if (c > 11.0 && c < 18.0) _h["Figure_20_11"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 11-18% centrality
            if (c > 18.0 && c < 26.0) _h["Figure_20_12"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 18-26% centrality
            if (c > 26.0 && c < 34.0) _h["Figure_20_13"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 26-34% centrality
            if (c > 34.0 && c < 45.0) _h["Figure_20_14"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 34-45% centrality
            if (c > 45.0 && c < 58.0) _h["Figure_20_15"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 45-58% centrality
            if (c > 58.0 && c < 85.0) _h["Figure_20_16"]->fill(p.pT()/GeV, 1.0/(2.0*3.14*p.pT()/GeV)); // 58-85% centrality
          }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h["XXXX"]); // normalize to unity
      //normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
      //scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

      //_h["AAAA"]->scaleW(1./_c["CCCC"]->sumW());

      if (collSys == dAu200)
      {
        _h["Figure_18_kaon_1"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_kaon_2"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_kaon_3"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
	      _h["Figure_18_kaon_4"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
	      _h["Figure_18_kaon_5"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
	      _h["Figure_18_kaon_6"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_kaon_7"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
	      _h["Figure_18_kaon_8"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_pion_1"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_pion_2"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_pion_3"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_pion_4"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_pion_5"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_pion_6"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_pion_7"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_pion_8"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_proton_1"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_proton_2"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_proton_3"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_proton_4"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_proton_5"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_proton_6"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_proton_7"]->scaleW(1.0/_c["sow_dAu200"]->sumW());
        _h["Figure_18_proton_8"]->scaleW(1.0/_c["sow_dAu200"]->sumW());

        /*normalize(_h["Figure_18_kaon_1"]);
        normalize(_h["Figure_18_kaon_2"]);
        normalize(_h["Figure_18_kaon_3"]);
        normalize(_h["Figure_18_kaon_4"]);
        normalize(_h["Figure_18_kaon_5"]);
        normalize(_h["Figure_18_kaon_6"]);
        normalize(_h["Figure_18_kaon_7"]);
        normalize(_h["Figure_18_kaon_8"]);
        normalize(_h["Figure_18_pion_1"]);
        normalize(_h["Figure_18_pion_2"]);
        normalize(_h["Figure_18_pion_3"]);
        normalize(_h["Figure_18_pion_4"]);
        normalize(_h["Figure_18_pion_5"]);
        normalize(_h["Figure_18_pion_6"]);
        normalize(_h["Figure_18_pion_7"]);
        normalize(_h["Figure_18_pion_8"]);
        normalize(_h["Figure_18_proton_1"]);
        normalize(_h["Figure_18_proton_2"]);
        normalize(_h["Figure_18_proton_3"]);
        normalize(_h["Figure_18_proton_4"]);
        normalize(_h["Figure_18_proton_5"]);
        normalize(_h["Figure_18_proton_6"]);
        normalize(_h["Figure_18_proton_7"]);
        normalize(_h["Figure_18_proton_8"]);

        scale(_h["Figure_18_kaon_1"], 4.0);
        scale(_h["Figure_18_kaon_2"], 2.0);
        scale(_h["Figure_18_kaon_3"], 1.0);
        scale(_h["Figure_18_kaon_4"], 0.5);
        scale(_h["Figure_18_kaon_5"], 4.0);
        scale(_h["Figure_18_kaon_6"], 2.0);
        scale(_h["Figure_18_kaon_7"], 1.0);
        scale(_h["Figure_18_kaon_8"], 0.5);
        scale(_h["Figure_18_pion_1"], 4.0);
        scale(_h["Figure_18_pion_2"], 2.0);
        scale(_h["Figure_18_pion_3"], 1.0);
        scale(_h["Figure_18_pion_4"], 0.5);
        scale(_h["Figure_18_pion_5"], 4.0);
        scale(_h["Figure_18_pion_6"], 2.0);
        scale(_h["Figure_18_pion_7"], 1.0);
        scale(_h["Figure_18_pion_8"], 0.5);
        scale(_h["Figure_18_proton_1"], 4.0);
        scale(_h["Figure_18_proton_2"], 2.0);
        scale(_h["Figure_18_proton_3"], 1.0);
        scale(_h["Figure_18_proton_4"], 0.5);
        scale(_h["Figure_18_proton_5"], 4.0);
        scale(_h["Figure_18_proton_6"], 2.0);
        scale(_h["Figure_18_proton_7"], 1.0);
        scale(_h["Figure_18_proton_8"], 0.5);*/
      }

      if (collSys == AuAu62)
      {
        _h["Figure_19_kaon_1"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_kaon_2"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_kaon_3"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_kaon_4"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_kaon_5"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_kaon_6"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_kaon_7"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_kaon_8"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_kaon_9"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_kaon_10"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_kaon_11"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_kaon_12"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_kaon_13"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_kaon_14"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_kaon_15"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_kaon_16"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_kaon_17"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_kaon_18"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_pion_1"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_pion_2"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_pion_3"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_pion_4"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_pion_5"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_pion_6"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_pion_7"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_pion_8"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_pion_9"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_pion_10"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_pion_11"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_pion_12"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_pion_13"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_pion_14"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_pion_15"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_pion_16"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_pion_17"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
	      _h["Figure_19_pion_18"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_1"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_2"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_3"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_4"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_5"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_6"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_7"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_8"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_9"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_10"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_11"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_12"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_13"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_14"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_15"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_16"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_17"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());
        _h["Figure_19_proton_18"]->scaleW(1.0/_c["sow_AuAu62"]->sumW());

        /*normalize(_h["Figure_19_kaon_1"]);
        normalize(_h["Figure_19_kaon_2"]);
        normalize(_h["Figure_19_kaon_3"]);
        normalize(_h["Figure_19_kaon_4"]);
        normalize(_h["Figure_19_kaon_5"]);
        normalize(_h["Figure_19_kaon_6"]);
        normalize(_h["Figure_19_kaon_7"]);
        normalize(_h["Figure_19_kaon_8"]);
        normalize(_h["Figure_19_kaon_9"]);
        normalize(_h["Figure_19_kaon_10"]);
        normalize(_h["Figure_19_kaon_11"]);
        normalize(_h["Figure_19_kaon_12"]);
        normalize(_h["Figure_19_kaon_13"]);
        normalize(_h["Figure_19_kaon_14"]);
        normalize(_h["Figure_19_kaon_15"]);
        normalize(_h["Figure_19_kaon_16"]);
        normalize(_h["Figure_19_kaon_17"]);
        normalize(_h["Figure_19_kaon_18"]);
        normalize(_h["Figure_19_pion_1"]);
        normalize(_h["Figure_19_pion_2"]);
        normalize(_h["Figure_19_pion_3"]);
        normalize(_h["Figure_19_pion_4"]);
        normalize(_h["Figure_19_pion_5"]);
        normalize(_h["Figure_19_pion_6"]);
        normalize(_h["Figure_19_pion_7"]);
        normalize(_h["Figure_19_pion_8"]);
        normalize(_h["Figure_19_pion_9"]);
        normalize(_h["Figure_19_pion_10"]);
        normalize(_h["Figure_19_pion_11"]);
        normalize(_h["Figure_19_pion_12"]);
        normalize(_h["Figure_19_pion_13"]);
        normalize(_h["Figure_19_pion_14"]);
        normalize(_h["Figure_19_pion_15"]);
        normalize(_h["Figure_19_pion_16"]);
        normalize(_h["Figure_19_pion_17"]);
        normalize(_h["Figure_19_pion_18"]);
        normalize(_h["Figure_19_proton_1"]);
        normalize(_h["Figure_19_proton_2"]);
        normalize(_h["Figure_19_proton_3"]);
        normalize(_h["Figure_19_proton_4"]);
        normalize(_h["Figure_19_proton_5"]);
        normalize(_h["Figure_19_proton_6"]);
        normalize(_h["Figure_19_proton_7"]);
        normalize(_h["Figure_19_proton_8"]);
        normalize(_h["Figure_19_proton_9"]);
        normalize(_h["Figure_19_proton_10"]);
        normalize(_h["Figure_19_proton_11"]);
        normalize(_h["Figure_19_proton_12"]);
        normalize(_h["Figure_19_proton_13"]);
        normalize(_h["Figure_19_proton_14"]);
        normalize(_h["Figure_19_proton_15"]);
        normalize(_h["Figure_19_proton_16"]);
        normalize(_h["Figure_19_proton_17"]);
        normalize(_h["Figure_19_proton_18"]);*/
      }

      if (collSys == AuAu130)
      {
        _h["Figure_20_1"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
        _h["Figure_20_2"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
        _h["Figure_20_3"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
	      _h["Figure_20_4"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
	      _h["Figure_20_5"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
	      _h["Figure_20_6"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
        _h["Figure_20_7"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
	      _h["Figure_20_8"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
        _h["Figure_20_9"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
        _h["Figure_20_10"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
        _h["Figure_20_12"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
        _h["Figure_20_13"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
        _h["Figure_20_14"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
        _h["Figure_20_15"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());
        _h["Figure_20_16"]->scaleW(1.0/_c["sow_AuAu130"]->sumW());

        /*normalize(_h["Figure_20_1"]);
        normalize(_h["Figure_20_2"]);
        normalize(_h["Figure_20_3"]);
        normalize(_h["Figure_20_4"]);
        normalize(_h["Figure_20_5"]);
        normalize(_h["Figure_20_6"]);
        normalize(_h["Figure_20_7"]);
        normalize(_h["Figure_20_8"]);
        normalize(_h["Figure_20_9"]);
        normalize(_h["Figure_20_10"]);
        normalize(_h["Figure_20_11"]);
        normalize(_h["Figure_20_12"]);
        normalize(_h["Figure_20_13"]);
        normalize(_h["Figure_20_14"]);
        normalize(_h["Figure_20_15"]);
        normalize(_h["Figure_20_16"]);*/
      }

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Histo2DPtr> _h2D;
    string beamOpt = "";
    //string beam = "";
    enum CollisionSystem {dAu200, AuAu62, AuAu130};
    CollisionSystem collSys;
    //@}


  };


  DECLARE_RIVET_PLUGIN(STAR_2009_I793126);

}
