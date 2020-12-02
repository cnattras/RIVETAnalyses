// -*- C++ -*-
#include <boost/math/special_functions/bessel.hpp> 

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

#include "../Centralities/RHICCentrality.hh"

const int NCEN = 6; // Centrality : 0-10, 10-20, 20-30, 30-40, 40-50, 50-60\%
const int NPTB = 2; // Pt : 0.75 - 1.0, 1.75 - 2.0 GeV/$c$
const int NHAR = 3; // Harmonics : $v_2$, $v_3$, $v_4$ 
const int NSBE = 3; // EP Subevent, 1 < |eta| < 2.8, -2.8 < eta < -1, 1 < eta < 2.8

namespace Rivet {


  /// @brief Azimuthal anisotropy for charged hadrons in Au+Au collisions at $\sqrt{s_{NN}}$=200 GeV.
  class PHENIX_2011_I900703 : public Analysis {
    public:

      /// Constructor
      DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2011_I900703);

      /// @name Analysis methods
      //@{
      double resEventPlane(double chi) {
        //Calculates the event plane resolution as a function of chi

        double pi  = acos(-1.);
        double con = sqrt(pi/2.)/2.;
        double arg = chi * chi / 4.;

        double res = con * chi * exp(-arg) * (cyl_bessel_i(0, arg) + cyl_bessel_i(1, arg));
        return res;
      }

      double chi(double res) {
        // Calculates chi from the event plane resolution
        double chi = 2.0;
        double delta = 1.0;

        for (int i = 0; i < 200; i++) {
          chi = (resEventPlane(chi) < res) ? chi + delta : chi - delta;
          delta = delta / 2.;
        }

        return chi;
      }

      // Histogram Names
      const char *name_vn_pt[NCEN][NHAR] =
      {
        // $v_2$        $v_3$          $v_4$ 
        {"d01-x01-y01", "d01-x01-y02", "d01-x01-y03"}, // 0-10%
        {"d02-x01-y01", "d02-x01-y02", "d02-x01-y03"}, // 10-20%
        {"d03-x01-y01", "d03-x01-y02", "d03-x01-y03"}, // 20-30%
        {"d04-x01-y01", "d04-x01-y02", "d04-x01-y03"}, // 30-40%
        {"d05-x01-y01", "d05-x01-y02", "d05-x01-y03"}, // 40-50%
        {"d06-x01-y01", "d06-x01-y02", ""}             // 50-60%, No $v_4$ measurment in 50-60%
      };

      const char *name_vn_cen[NPTB][NHAR] =
      {
        // $v_2$        $v_3$          $v_4$ 
        {"d08-x01-y01", "d08-x01-y03", "d09-x01-y01"}, // pT = 0.75 - 1 GeV/c
        {"d08-x01-y02", "d08-x01-y04", "d09-x01-y02"}  // pT = 1.75 - 2 GeV/c
      };

      const char *name_epcor[NHAR] =
      {
        // Psi2    Psi3       Psi4 
        "Psi2Cor", "Psi3Cor", "Psi4Cor"
      };

      /// Book histograms and initialise projections before the run
      void init() {

        // Initialise and register projections
        // Declare beamOpt
        beamOpt = getOption<string>("beam", "NONE");

        // Declare Centrality
        declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

        // Declare Particles in this RIVET Analysis
        // Charged Particles for EP determination and Vn measurements
        const FinalState fs(Cuts::abseta < 2.8);
        declare(fs,"fs");

        // Book Histograms
        for(int icen = 0; icen<NCEN; icen++){
          for(int ihar = 0; ihar<NHAR; ihar++){
            if(icen==5 && ihar==2) continue; // No $v_4$ measurment in 50-60%
            book(_p_vn_pt[name_vn_pt[icen][ihar]], icen+1, 1, ihar+1);
          }
        }

        for(int iptb = 0; iptb<NPTB; iptb++){
          for(int ihar = 0; ihar<NHAR; ihar++){

            int dd = 8;
            if(ihar==2) dd = 9;

            int yy = iptb + 1;
            if(ihar==1) yy = iptb + 1 + 2*ihar;
            book(_p_vn_cen[name_vn_cen[iptb][ihar]], dd, 1, yy);
          }
        }

        for(int ihar = 0; ihar<NHAR; ihar++){
          book(_p_epcor_cen[name_epcor[ihar]], name_epcor[ihar], NCEN, 0., 60.);
        }

      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {

        // A. Get Centrality
        const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
        const double cent = centProj();
        const int icen = (int) cent / 10.; // Categorize centrality in 10% step.
        if (icen < 0 || icen >= 6) vetoEvent;

        // B. Calculate Event Plane
        double Qx[NHAR][NSBE]={0.}; // Q-vector in x
        double Qy[NHAR][NSBE]={0.}; // Q-vector in y
        double ep[NHAR][NSBE]={0.}; // Event Plane

        // Fill Q-vectors
        Particles particles = applyProjection<FinalState>(event,"fs").particles();
        for (const Particle& p :particles) {
          double eta = p.pseudorapidity();
          double phi = p.phi();
          for(int ihar = 0; ihar<NHAR; ihar++){
            double qx = cos( (ihar + 2) * phi);
            double qy = sin( (ihar + 2) * phi);
            if( fabs(eta) > 1 && fabs(eta) < 2.8) {
              Qx[ihar][0] += qx;
              Qy[ihar][0] += qy;
            }
            if( eta > -2.8 && eta < -1.0) {
              Qx[ihar][1] += qx;
              Qy[ihar][1] += qy;
            }
            else if( eta > 1.0 && eta < 2.8) {
              Qx[ihar][2] += qx;
              Qy[ihar][2] += qy;
            }
          }
        }

        // Determine Evnet Plane
        for(int ihar = 0; ihar<NHAR; ihar++){
          for(int isbe = 0; isbe<NSBE; isbe++){
            ep[ihar][isbe] = atan2(Qy[ihar][isbe], Qx[ihar][isbe]) / (ihar + 2);
          }
        }

        //Measure Event Plane Correlations
        for(int ihar = 0; ihar<NHAR; ihar++){
          _p_epcor_cen[name_epcor[ihar]]->fill(cent, cos( (ihar + 2) * (ep[ihar][1] - ep[ihar][2])));
        }

        //C. Fill vn vs pT and centrality
        // vn vs pT and centrality
        for (const Particle& p :particles) {
          const double pT = p.pt() / GeV;
          const double eta = p.pseudorapidity();
          const double phi = p.phi();
          const int abscharge = p.abscharge();

          if(fabs(eta) > 0.35) continue; // PHENIX acceptance
          if(abscharge < 1) continue; // Select charged particles
          for(int ihar = 0; ihar<NHAR; ihar++){
            if(icen==5 && ihar==2) continue; // No $v_4$ measurment in 50-60%
            _p_vn_pt[name_vn_pt[icen][ihar]]->fill(pT, cos( (ihar + 2.) * (phi - ep[ihar][0])));
          }

          int iptb = -1;
          if(pT > 0.75 && pT < 1.0 ) iptb = 0;
          else if(pT > 1.75 && pT < 2.0 ) iptb = 1;
          if(iptb>-1){
            for(int ihar = 0; ihar<NHAR; ihar++){
              _p_vn_cen[name_vn_cen[iptb][ihar]]->fill(cent, cos( (ihar + 2.) * (phi - ep[ihar][0])));
            }
          }
        }

      }

      /// Normalise histograms etc., after the run
      void finalize() {

        //D. Scale vn with the EP resolutions
        for(int icen = 0; icen<NCEN; icen++){ 
          for(int ihar = 0; ihar<NHAR; ihar++){
            if(icen==5 && ihar==2) continue; // No $v_4$ measurment in 50-60%
            double epcor = _p_epcor_cen[name_epcor[ihar]]-> bin(icen).mean();
            double chi_sub = chi(sqrt(epcor));
            double res = resEventPlane(sqrt(2.0)*chi_sub);
            _p_vn_pt[name_vn_pt[icen][ihar]] -> scaleY(1./res);  
          }
        }

        for(int iptb = 0; iptb<NPTB; iptb++){
          for(int ihar = 0; ihar<NHAR; ihar++){
            for(int icen = 0; icen<NCEN; icen++){ 
              if(icen==5 && ihar==2) continue; // No $v_4$ measurment in 50-60%
              double epcor = _p_epcor_cen[name_epcor[ihar]]-> bin(icen).mean();
              double chi_sub = chi(sqrt(epcor));
              double res = resEventPlane(sqrt(2.0)*chi_sub);
              _p_vn_cen[name_vn_cen[iptb][ihar]] -> bin(icen).scaleY(1./res);  
            }
          }
        }

      }

      //@}

      map<string, Profile1DPtr> _p_vn_pt;
      map<string, Profile1DPtr> _p_vn_cen;
      map<string, Profile1DPtr> _p_epcor_cen;
      string beamOpt = "";

      /// @name Histograms
      //@{
      //@}


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2011_I900703);

}
