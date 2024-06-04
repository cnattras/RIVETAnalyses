// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/HadronicFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2001_I562409 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2001_I562409);


    /// @name Analysis methods
    ///@{

    //create binShift function
    void binShift(YODA::Histo1D& histogram) {
        std::vector<YODA::HistoBin1D> binlist = histogram.bins();
        int n = 0;
        for (YODA::HistoBin1D bins : binlist) {
            double p_high = bins.xMax();
            double p_low = bins.xMin();
            //Now calculate f_corr
            if (bins.xMin() == binlist[0].xMin()) { //Check if we are working with first bin
                float b = 1 / (p_high - p_low) * log(binlist[0].height()/binlist[1].height());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
                n += 1;
            } else if (bins.xMin() == binlist.back().xMin()){ //Check if we are working with last bin
                float b = 1 / (p_high - p_low) * log(binlist[binlist.size()-2].height() / binlist.back().height());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
            } else { //Check if we are working with any middle bin
                float b = 1 / (p_high - p_low) * log(binlist[n-1].height() / binlist[n+1].height());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
                n += 1;
            }
        }
    }

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const HadronicFinalState hfs(Cuts::abseta < 0.5 && Cuts::abscharge > 0);
      declare(hfs, "hfs");

      const UnstableParticles pi0UP(Cuts::abseta < 0.5 && Cuts::abspid == 111);
      declare(pi0UP,"pi0UP");

      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_h["ChHadronsCent0_10"], 5, 1, 1);
      book(_h["ChHadronsCent60_80"], 6, 1, 1);
      book(_c["Cent0_10"], "_Cent0_10");
      book(_c["Cent60_80"], "_Cent60_80");
      book(_c["pp_130"], "_Sow_pp");
      book(_h["Pi0PbScCent0_10"], 1, 1, 1);
      book(_h["Pi0PbScCent60_80"], 2, 1, 1);
      book(_h["Pi0PbGlCent0_10"], 3, 1, 1);
      book(_h["Pi0PbGlCent60_80"], 4, 1, 1);
      string refnameRaa1  = mkAxisCode(7,1,1);
      const Scatter2D& refdataRaa1 =refData(refnameRaa1);
      book(_h["c10Pt_AuAu130_t7"], "_" + refnameRaa1 + "_AuAu130", refdataRaa1);
      book(_h["c10Pt130_pp_t7"], "_" + refnameRaa1 + "_pp130", refdataRaa1);
      book(_s["Raa_c010_AuAu130_t7"], refnameRaa1);
      string refnameRaa2  = mkAxisCode(8,1,1);
      const Scatter2D& refdataRaa2 =refData(refnameRaa2);
      book(_h["c10Pt_AuAu130_t8"], "_" + refnameRaa2 + "_AuAu130", refdataRaa2);
      book(_h["c10Pt130_pp_t8"], "_" + refnameRaa2 + "_pp130", refdataRaa2);
      book(_s["Raa_c010_AuAu130_t8"], refnameRaa2);
      string refnameRaa3  = mkAxisCode(9,1,1);
      const Scatter2D& refdataRaa3 =refData(refnameRaa3);
      book(_h["c10Pt_AuAu130_t9"], "_" + refnameRaa3 + "_AuAu130", refdataRaa3);
      book(_h["c10Pt130_pp_t9"], "_" + refnameRaa3 + "_pp130", refdataRaa3);
      book(_s["Raa_c010_AuAu130_t9"], refnameRaa3);
      
      string refnameRCP1  = mkAxisCode(10,1,1);
      const Scatter2D& refdataRCP1 =refData(refnameRCP1);
      book(_h["c10Pt_AuAu130_t10"], "_" + refnameRCP1 + "_cAuAu130", refdataRCP1);
      book(_h["p80Pt130_AuAu130_t10"], "_" + refnameRCP1 + "_pAuAu130", refdataRCP1);
      book(_s["RCP_AuAu130_t10"], refnameRCP1);

      string refnameRCP2  = mkAxisCode(11,1,1);
      const Scatter2D& refdataRCP2 =refData(refnameRCP2);
      book(_h["c10Pt_AuAu130_t11"], "_" + refnameRCP2 + "_cAuAu130", refdataRCP2);
      book(_h["p80Pt130_AuAu130_t11"], "_" + refnameRCP2 + "_pAuAu130", refdataRCP2);
      book(_s["RCP_AuAu130_t11"], refnameRCP2);

      string refnameRCP3  = mkAxisCode(12,1,1);
      const Scatter2D& refdataRCP3 =refData(refnameRCP3);
      book(_h["c10Pt_AuAu130_t12"], "_" + refnameRCP3 + "_cAuAu130", refdataRCP3);
      book(_h["p80Pt130_AuAu130_t12"], "_" + refnameRCP3 + "_pAuAu130", refdataRCP3);
      book(_s["RCP_AuAu130_t12"], refnameRCP3);


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
      const double cent = centProj();
      
      const HadronicFinalState hfs = apply<HadronicFinalState>(event, "hfs");
      const Particles hfsParticles = hfs.particles();

      const UnstableParticles pi0UP = apply<UnstableParticles>(event, "pi0UP");
      const Particles pi0UPParticles = pi0UP.particles();

      double inv2PI = 1./(2.*M_PI);
      const ParticlePair& beam = beams();
      string collSys;

      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
    {
    collSys = "AuAu";
}
if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
    {
    collSys = "pp";
    }
            if (collSys == "pp")
            {
              _c["pp_130"]->fill();
              for(auto p : hfsParticles) //loop over charged hadrons
              {
                      _h["c10Pt130_pp_t9"]->fill(p.pT()/GeV); 
          
              }
              for(auto p : pi0UPParticles) //loop over pi0s
              {
		      _h["c10Pt130_pp_t7"]->fill(p.pT()/GeV);
          _h["c10Pt130_pp_t8"]->fill(p.pT()/GeV);
         
              }
            }
          else if (collSys == "AuAu")
          {
             if(cent < 10.) //Check centrality of the event
      {
              _c["Cent0_10"]->fill(); //fill counter for 0-10% most central events
              for(auto p : hfsParticles) //loop over charged hadrons
              {
                      _h["ChHadronsCent0_10"]->fill(p.pT()/GeV, inv2PI/(2.*p.pT()/GeV)); //additional 1/2 factor to take into account h^(+)+h^(-)/2
                      _h["c10Pt_AuAu130_t9"]->fill(p.pT()/GeV);
                      _h["c10Pt_AuAu130_t12"]->fill(p.pT()/GeV);
              }
              
              for(auto p : pi0UPParticles) //loop over pi0s
              {
		      _h["Pi0PbScCent0_10"]->fill(p.pT()/GeV, inv2PI/(p.pT()/GeV));
          _h["Pi0PbGlCent0_10"]->fill(p.pT()/GeV, inv2PI/(p.pT()/GeV));
          _h["c10Pt_AuAu130_t7"]->fill(p.pT()/GeV);
          _h["c10Pt_AuAu130_t8"]->fill(p.pT()/GeV);
          _h["c10Pt_AuAu130_t10"]->fill(p.pT()/GeV);
          _h["c10Pt_AuAu130_t11"]->fill(p.pT()/GeV);
              }
      }
      else if(cent >= 60. && cent < 80.) //Check centrality of the event
      {
              _c["Cent60_80"]->fill(); //fill counter for 60-80% most central events
              for(auto p : hfsParticles) //loop over charged hadrons
              {
                      _h["ChHadronsCent60_80"]->fill(p.pT()/GeV, inv2PI/(2.*p.pT()/GeV)); //additional 1/2 factor to take into account h^(+)+h^(-)/2
                    _h["p80Pt130_AuAu130_t12"]->fill(p.pT()/GeV);
              }

              for(auto p : pi0UPParticles) //loop over pi0s
              {
                      _h["Pi0PbScCent60_80"]->fill(p.pT()/GeV, inv2PI/(p.pT()/GeV));
		      _h["Pi0PbGlCent60_80"]->fill(p.pT()/GeV, inv2PI/(p.pT()/GeV)); 
          _h["p80Pt130_AuAu130_t10"]->fill(p.pT()/GeV);
          _h["p80Pt130_AuAu130_t11"]->fill(p.pT()/GeV);
              }
      }
          }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //Divide histograms by number of events
      _h["ChHadronsCent0_10"]->scaleW(1./_c["Cent0_10"]->sumW());
      _h["ChHadronsCent60_80"]->scaleW(1./_c["Cent60_80"]->sumW());
      _h["Pi0PbScCent0_10"]->scaleW(1./_c["Cent0_10"]->sumW());
      _h["Pi0PbScCent60_80"]->scaleW(1./_c["Cent60_80"]->sumW());
      _h["Pi0PbGlCent0_10"]->scaleW(1./_c["Cent0_10"]->sumW());
      _h["Pi0PbGlCent60_80"]->scaleW(1./_c["Cent0_10"]->sumW());
     _h["c10Pt_AuAu130_t7"]->scaleW(1./_c["Cent0_10"]->sumW());
     _h["c10Pt_AuAu130_t8"]->scaleW(1./_c["Cent0_10"]->sumW());
     _h["c10Pt_AuAu130_t9"]->scaleW(1./_c["Cent0_10"]->sumW());
     _h["c10Pt_AuAu130_t10"]->scaleW(1./_c["Cent0_10"]->sumW());
     _h["c10Pt_AuAu130_t11"]->scaleW(1./_c["Cent0_10"]->sumW());
     _h["c10Pt_AuAu130_t12"]->scaleW(1./_c["Cent0_10"]->sumW());
     _h["c10Pt130_pp_t7"]->scaleW(1./_c["pp_130"]->sumW());
     _h["c10Pt130_pp_t8"]->scaleW(1./_c["pp_130"]->sumW());
     _h["c10Pt130_pp_t9"]->scaleW(1./_c["pp_130"]->sumW());
     _h["p80Pt130_AuAu130_t10"]->scaleW(1./_c["Cent60_80"]->sumW());
     _h["p80Pt130_AuAu130_t11"]->scaleW(1./_c["Cent60_80"]->sumW());
     _h["p80Pt130_AuAu130_t12"]->scaleW(1./_c["Cent60_80"]->sumW());
      divide(_h["c10Pt_AuAu130_t7"], _h["c10Pt130_pp_t7"],_s["Raa_c010_AuAu130_t7"]);
      divide(_h["c10Pt_AuAu130_t8"], _h["c10Pt130_pp_t8"],_s["Raa_c010_AuAu130_t8"]);
      divide(_h["c10Pt_AuAu130_t9"], _h["c10Pt130_pp_t9"],_s["Raa_c010_AuAu130_t9"]);
      divide(_h["c10Pt_AuAu130_t10"], _h["p80Pt130_AuAu130_t10"], _s["RCP_AuAu130_t10"]);
      divide(_h["c10Pt_AuAu130_t11"], _h["p80Pt130_AuAu130_t11"], _s["RCP_AuAu130_t11"]);
      divide(_h["c10Pt_AuAu130_t12"], _h["p80Pt130_AuAu130_t12"], _s["RCP_AuAu130_t12"]);
     _s["Raa_c010_AuAu130_t7"]->scaleY(1./905.);
     _s["Raa_c010_AuAu130_t8"]->scaleY(1./905.);
     _s["Raa_c010_AuAu130_t9"]->scaleY(1./905.);
     _s["RCP_AuAu130_t10"]->scaleY(1./45.);
     _s["RCP_AuAu130_t11"]->scaleY(1./45.);
     _s["RCP_AuAu130_t12"]->scaleY(1./45.);
    }

      ///@}
      //@name Histograms
      ///@{
      map<string, Histo1DPtr> _h;
      map<string, Profile1DPtr> _p;
      map<string, CounterPtr> _c;
      map<string, Scatter2DPtr> _s;
      ///@}
  };


  RIVET_DECLARE_PLUGIN(PHENIX_2001_I562409);

}
