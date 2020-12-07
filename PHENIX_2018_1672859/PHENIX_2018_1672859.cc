// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "../Centralities/RHICCentrality.hh"

const double TAB[6] = {2.54, 8.8, 6.0, 7.5, 3.1, 1.0};

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2018_1672859 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2018_1672859);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init()
    {
      //~ eta pdg id: 221
      //~ pi0 pdg id: 111
      
      beamOpt = getOption<string>("beam", "NONE");
      declareCentrality(RHICCentrality("PHENIX"),"RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      // Initialise and register projections
      // The basic final-state projection: all final-state particles within the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);
      declare(fs, "fs");
      
      const UnstableParticles up(Cuts::abseta < 4.9);
      declare(UnstableParticles(), "UFS");
      // declare(up, "up");

      // Book histograms specify custom binning
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
    
    book(_h["pi0_minbias"], 1, 1, 1);
      // book(_h["pi0_0010"],  2, 1, 1);
      // book(_h["pi0_1020"],  3, 1, 1);
      // book(_h["pi0_0020"],  4, 1, 1);
      // book(_h["pi0_2040"],  5, 1, 1);
      // book(_h["pi0_4060"],  6, 1, 1);
      // book(_h["pi0_6090"],  7, 1, 1);
      
      book(_h["eta_minbias"],  8, 1, 1);
      // book(_h["eta_0020"],  9, 1, 1);
      // book(_h["eta_2040"], 10, 1, 1);
      // book(_h["eta_4060"], 11, 1, 1);
      // book(_h["eta_6090"], 12, 1, 1);
      
      string refname13 = mkAxisCode(13, 1, 1);
      const Scatter2D& refdata13 = refData(refname13);
      book(_h["eta_mb"], refname13 + "_eta_mb", refdata13);
      book(_h["pi0_mb"], refname13 + "_pi0_mb", refdata13);
      book(_s["eta_over_pi0_mb"], refname13);
      
      string refname18 = mkAxisCode(18, 1, 1);
      const Scatter2D& refdata18 = refData(refname18);
      book(_h["CuAu_minbias_pi0"], refname18 + "_CuAu_minbias", refdata18);
      book(_h["pp_pi0"], refname18 + "_pp", refdata18);
      book(_s["R_AB_pi0"], refname18);

      string refname23 = mkAxisCode(23, 1, 1);
      const Scatter2D& refdata23 = refData(refname23);
      book(_h["CuAu_minbias_eta"], refname23 + "_CuAu_minbias", refdata23);
      book(_h["pp_eta"], refname23 + "_pp", refdata23);
      book(_s["R_AB_eta"], refname23);
      
      /* // Histos for other plots
      // string refname14 = mkAxisCode(14, 1, 1);
      // const Scatter2D& refdata14 = refData(refname14);
      // book(_h["eta_0p"], refname14 + "_eta_0p", refdata14);
      // book(_h["pi0_0p"], refname14 + "_pi0_0p", refdata14);
      // book(_s["eta_over_pi0_0p"], refname14);
      
      // string refname15 = mkAxisCode(15, 1, 1);
      // const Scatter2D& refdata15 = refData(refname15);
      // book(_h["eta_20p"], refname15 + "_eta_20p", refdata15);
      // book(_h["pi0_20p"], refname15 + "_pi0_20p", refdata15);
      // book(_s["eta_over_pi0_20p"], refname15);
      
      // string refname16 = mkAxisCode(16, 1, 1);
      // const Scatter2D& refdata16 = refData(refname16);
      // book(_h["eta_40p"], refname16 + "_eta_40p", refdata16);
      // book(_h["pi0_40p"], refname16 + "_pi0_40p", refdata16);
      // book(_s["eta_over_pi0_40p"], refname16);
      
      // string refname17 = mkAxisCode(17, 1, 1);
      // const Scatter2D& refdata17 = refData(refname17);
      // book(_h["eta_60p"], refname17 + "_eta_60p", refdata17);
      // book(_h["pi0_60p"], refname17 + "_pi0_60p", refdata17);
      // book(_s["eta_over_pi0_60p"], refname17);
      */
    }


    /// Perform the per-event analysis
    void analyze(const Event& event)
    {
      const UnstableParticles& usp = applyProjection<UnstableParticles>(event, "UFS");
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();
      
      for(const Particle & p : usp.particles())
      {
        if(beamOpt == "CUAU200")
        {
          if(p.pid() == PID::PI0)
          {
            _h["pi0_minbias"]->fill(p.pT()/GeV);
            _h["pi0_mb"]->fill(p.pT()/GeV);
            // _h["CuAu_minbias_pi0"]->fill(p.pT()/GeV);
            
            // if(c <= 10.) _h["pi0_0010"]->fill(p.pT()/GeV);
            // if(c > 10 && c < 20.) _h["pi0_1020"]->fill(p.pT()/GeV);
            // if(c > 00 && c < 20.) _h["pi0_0020"]->fill(p.pT()/GeV);
            // if(c > 20 && c < 40.) _h["pi0_2040"]->fill(p.pT()/GeV);
            // if(c > 40 && c < 60.) _h["pi0_4060"]->fill(p.pT()/GeV);
            // if(c > 60 && c < 90.) _h["pi0_6090"]->fill(p.pT()/GeV);
          }
          else if(p.pid() == PID::ETA)
          {
            _h["eta_minbias"]->fill(p.pT()/GeV);
            _h["eta_mb"]->fill(p.pT()/GeV);
            // _h["CuAu_minbias_eta"]->fill(p.pT()/GeV);
            
            // if(c > 00 && c < 20.) _h["eta_0020"]->fill(p.pT()/GeV);
            // if(c > 20 && c < 40.) _h["eta_2040"]->fill(p.pT()/GeV);
            // if(c > 40 && c < 60.) _h["eta_4060"]->fill(p.pT()/GeV);
            // if(c > 60 && c < 90.) _h["eta_6090"]->fill(p.pT()/GeV);
          }
        }        
        else if(beamOpt == "PP200")
        {          
          if(p.pid() == PID::PI0)
          {
            _h["pp_pi0"]->fill(p.pT()/GeV);
            // if(c <= 10.) _h["pi0_0010"]->fill(p.pT()/GeV);
            // if(c > 10 && c < 20.) _h["pi0_1020"]->fill(p.pT()/GeV);
            // if(c > 00 && c < 20.) _h["pi0_0020"]->fill(p.pT()/GeV);
            // if(c > 20 && c < 40.) _h["pi0_2040"]->fill(p.pT()/GeV);
            // if(c > 40 && c < 60.) _h["pi0_4060"]->fill(p.pT()/GeV);
            // if(c > 60 && c < 90.) _h["pi0_6090"]->fill(p.pT()/GeV);
          }
          else
          {
            _h["pp_eta"]->fill(p.pT()/GeV);            
            // if(c > 00 && c < 20.) _h["eta_0020"]->fill(p.pT()/GeV);
            // if(c > 20 && c < 40.) _h["eta_2040"]->fill(p.pT()/GeV);
            // if(c > 40 && c < 60.) _h["eta_4060"]->fill(p.pT()/GeV);
            // if(c > 60 && c < 90.) _h["eta_6090"]->fill(p.pT()/GeV);
          }
        }
      }
      
    }


    /// Normalise histograms etc., after the run
    void finalize()
    {
      normalize(_h["pi0_minbias"]); // normalize to unity
      // normalize(_h["pi0_0010"]); // normalize to unity
      // normalize(_h["pi0_1020"]); // normalize to unity
      // normalize(_h["pi0_0020"]); // normalize to unity
      // normalize(_h["pi0_2040"]); // normalize to unity
      // normalize(_h["pi0_4060"]); // normalize to unity
      // normalize(_h["pi0_6090"]); // normalize to unity
      
      normalize(_h["eta_minbias"]); // normalize to unity
      // normalize(_h["eta_0020"]); // normalize to unity
      // normalize(_h["eta_2040"]); // normalize to unity
      // normalize(_h["eta_4060"]); // normalize to unity
      // normalize(_h["eta_6090"]); // normalize to unity
      
      divide(_h["eta_mb"],_h["pi0_mb"],_s["eta_over_pi0_mb"]);
      divide(_h["CuAu_minbias_pi0"],_h["pp_pi0"],_s["R_AB_pi0"]);
      divide(_h["CuAu_minbias_eta"],_h["pp_eta"],_s["R_AB_eta"]);
      
      bool pi0_pp_available = false;
      bool pi0_CuAu_available = false;
      bool eta_pp_available = false;
      bool eta_CuAu_available = false;

      for (auto element : _c)
      {
        string name = element.second->name();
        if (name.find("pi0_minbias") != std::string::npos)
        {
          if (element.second->sumW()>0) pi0_CuAu_available = true;
          else
          {
            pi0_CuAu_available=false;
            break;
          }
        }
        else if (name.find("pp_pi0") != std::string::npos)
        {
          if (element.second->sumW()>0) pi0_pp_available=true;
          else
          {
            pi0_pp_available=false;
            break;
          }
        }
        else if (name.find("eta_minbias") != std::string::npos)
        {
          if (element.second->sumW()>0) eta_pp_available=true;
          else
          {
            eta_pp_available=false;
            break;
          }
        }
        else if (name.find("pp_eta") != std::string::npos)
        {
          if (element.second->sumW()>0) eta_pp_available=true;
          else
          {
            eta_pp_available=false;
            break;
          }
        }
      }
      if((!pi0_pp_available) || (!pi0_CuAu_available) || (!eta_pp_available) || (!eta_CuAu_available)) return;
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;
    string beamOpt = "";
    ///@}


  };

  DECLARE_RIVET_PLUGIN(PHENIX_2018_1672859);

}
