// -*- C++ -*-
#include "Rivet/Analysis.hh"

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/AliceCommon.hh"

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/UnstableParticles.hh"

#include "../Centralities/RHICCentrality.hh" //external header for Centrality calculation
#include <math.h>
#define _USE_MATH_DEFINES


namespace Rivet {


  class PHENIX_2008_I778168 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2008_I778168);
    
    //create binShift function
    void binShift(YODA::Histo1D& histogram) {
        const auto& binlist = histogram.bins();
        const auto& lastBin = histogram.bin(histogram.numBins());
        int n = 0;
        for (auto& bins : binlist) {
            double p_high = bins.xMax();
            double p_low = bins.xMin();
            //Now calculate f_corr
            if (bins.xMin() == binlist[0].xMin()) { //Check if we are working with first bin
                float b = 1 / (p_high - p_low) * log(binlist[0].sumW()/binlist[1].sumW());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
                n += 1;
            } else if (bins.xMin() == lastBin.xMin()){ //Check if we are working with last bin
                float b = 1 / (p_high - p_low) * log(binlist[binlist.size()-2].sumW() / lastBin.sumW());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
            } else { //Check if we are working with any middle bin
                float b = 1 / (p_high - p_low) * log(binlist[n-1].sumW() / binlist[n+1].sumW());
                float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
                histogram.bin(n).scaleW(f_corr);
                n += 1;
            }
        }
    }

    void init() {

      //const PrimaryParticles fs(pdgIds, Cuts::abseta < 1. && Cuts::abscharge == 0);
      //declare(fs, "fs");

      const UnstableParticles up(Cuts::abseta < 0.35 && Cuts::pid == 111);
      declare(up, "up");

      beamOpt = getOption<string>("beam","NONE");
      
      const ParticlePair& beam = beams();


        if (beamOpt == "NONE") {

          if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
          {
            float NN = 197.;
            if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = AuAu;
          }
          if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
          {
            float NN = 1.;
            if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = pp;
          }
        }

        else if (beamOpt == "PP200") collSys = pp;
	    else if (beamOpt == "AUAU200") collSys = AuAu;

      //if(beamOpt=="PP") collSys = pp;
      //else if(beamOpt=="AUAU") collSys = AuAu;

    declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");


 //Yields_________________
       book(_h["ptyieldsc5"], 2, 1, 3);
       book(_h["ptyieldsc10"], 1, 1, 1);
       book(_h["ptyieldsc20"], 1, 1, 2);
       book(_h["ptyieldsc30"], 2, 1, 1);
       book(_h["ptyieldsc40"], 2, 1, 2);
       book(_h["ptyieldsc50"], 3, 1, 1);
       book(_h["ptyieldsc60"], 3, 1, 2);
       book(_h["ptyieldsc70"], 4, 1, 1);
       book(_h["ptyieldsc80"], 5, 1, 1);
       book(_h["ptyieldsc92"], 5, 1, 2);
       book(_h["ptyieldscall"], 1, 1, 3);

       book(_c["sow_AuAu_c5"],"_sow_AuAu_c5");
       book(_c["sow_AuAu_c10"],"_sow_AuAu_c10");
       book(_c["sow_AuAu_c20"],"_sow_AuAu_c20");
       book(_c["sow_AuAu_c30"],"_sow_AuAu_c30");
       book(_c["sow_AuAu_c40"],"_sow_AuAu_c40");
       book(_c["sow_AuAu_c50"],"_sow_AuAu_c50");
       book(_c["sow_AuAu_c60"],"_sow_AuAu_c60");
       book(_c["sow_AuAu_c70"],"_sow_AuAu_c70");
       book(_c["sow_AuAu_c80"],"_sow_AuAu_c80");
       book(_c["sow_AuAu_c92"],"_sow_AuAu_c92");
       book(_c["sow_AuAu_call"],"_sow_AuAu_call");

 //RAA _______________________________
       string refnameRaa1 = mkAxisCode(11,1,1);
             const Estimate1D& refdataRaa1 =refData(refnameRaa1);
       book(_h["c5Pt_AuAu"],"_" + refnameRaa1 + "_AuAu", refdataRaa1);
       book(_h["c5Pt_pp"], "_" +refnameRaa1 + "_pp", refdataRaa1);
       book(hRaa["Raa_c05_AuAu"], refnameRaa1);

       string refnameRaa2 = mkAxisCode(6,1,1);
             const Estimate1D& refdataRaa2 =refData(refnameRaa2);
       book(_h["c10Pt_AuAu"],"_" + refnameRaa2 + "_AuAu", refdataRaa2);
       book(_h["c10Pt_pp"],"_" + refnameRaa2 + "_pp", refdataRaa2);
       book(hRaa["Raa_c010_AuAu"], refnameRaa2);

       string refnameRaa3 = mkAxisCode(6,1,2);
             const Estimate1D& refdataRaa3 =refData(refnameRaa3);
       book(_h["c20Pt_AuAu"],"_" + refnameRaa3 + "_AuAu", refdataRaa3);
       book(_h["c20Pt_pp"],  "_" +refnameRaa3 + "_pp", refdataRaa3);
       book(hRaa["Raa_c1020_AuAu"],"_" + refnameRaa3);

       string refnameRaa4 = mkAxisCode(7,1,1);
             const Estimate1D& refdataRaa4 =refData(refnameRaa4);
       book(_h["c30Pt_AuAu"],"_" + refnameRaa4 + "_AuAu", refdataRaa4);
       book(_h["c30Pt_pp"],"_" + refnameRaa4 + "_pp", refdataRaa4);
       book(hRaa["Raa_c2030_AuAu"], refnameRaa4);

       string refnameRaa5 = mkAxisCode(7,1,2);
             const Estimate1D& refdataRaa5 =refData(refnameRaa5);
       book(_h["c40Pt_AuAu"], "_" + refnameRaa5 + "_AuAu", refdataRaa5);
       book(_h["c40Pt_pp"], "_" + refnameRaa5 + "_pp", refdataRaa5);
       book(hRaa["Raa_c3040_AuAu"], refnameRaa5);

       string refnameRaa6 = mkAxisCode(8,1,1);
             const Estimate1D& refdataRaa6 =refData(refnameRaa6);
       book(_h["c50Pt_AuAu"], "_" + refnameRaa6 + "_AuAu", refdataRaa6);
       book(_h["c50Pt_pp"], "_" + refnameRaa6 + "_pp", refdataRaa6);
       book(hRaa["Raa_c4050_AuAu"], refnameRaa6);

       string refnameRaa7 = mkAxisCode(8,1,2);
             const Estimate1D& refdataRaa7 =refData(refnameRaa7);
       book(_h["c60Pt_AuAu"], "_" + refnameRaa7 + "_AuAu", refdataRaa7);
       book(_h["c60Pt_pp"], "_" + refnameRaa7 + "_pp", refdataRaa7);
       book(hRaa["Raa_c5060_AuAu"], refnameRaa7);

       string refnameRaa8 = mkAxisCode(9,1,1);
             const Estimate1D& refdataRaa8 =refData(refnameRaa8);
       book(_h["c70Pt_AuAu"], "_" + refnameRaa8 + "_AuAu", refdataRaa8);
       book(_h["c70Pt_pp"], "_" + refnameRaa8 + "_pp", refdataRaa8);
       book(hRaa["Raa_c6070_AuAu"], refnameRaa8);

       string refnameRaa9 = mkAxisCode(10,1,1);
             const Estimate1D& refdataRaa9 =refData(refnameRaa9);
       book(_h["c80Pt_AuAu"], "_" + refnameRaa9 + "_AuAu", refdataRaa9);
       book(_h["c80Pt_pp"], "_" + refnameRaa9 + "_pp", refdataRaa9);
       book(hRaa["Raa_c7080_AuAu"], refnameRaa9);

       string refnameRaa10 = mkAxisCode(10,1,2);
             const Estimate1D& refdataRaa10 =refData(refnameRaa10);
       book(_h["c92Pt_AuAu"], "_" + refnameRaa10 + "_AuAu", refdataRaa10);
       book(_h["c92Pt_pp"], "_" + refnameRaa10 + "_pp", refdataRaa10);
       book(hRaa["Raa_c8092_AuAu"], refnameRaa10);

       string refnameRaa11 = mkAxisCode(6,1,3);
             const Estimate1D& refdataRaa11 =refData(refnameRaa11);
       book(_h["callPt_AuAu"], "_" + refnameRaa11 + "_AuAu", refdataRaa11);
       book(_h["callPt_pp"], "_" + refnameRaa11 + "_pp", refdataRaa11);
       book(hRaa["Raa_minbias_AuAu"], refnameRaa11);


       book(_c["sow_pp"],"_sow_pp");

 //Centrality vs RAA
       /*
       string refnameRaa13 = mkAxisCode(11,1,1);
       const Estimate1D& refdataRaa13 =refData(refnameRaa13);
       book(hRaaNpart["Raa_pt520c05_AuAu"], refnameRaa13, refdataRaa13);

      string refnameRaa14 = mkAxisCode(12,1,1);
       const Estimate1D& refdataRaa14 =refData(refnameRaa14);
       book(hRaaNpart["Raa_pt520c010_AuAu"], refnameRaa14, refdataRaa14);

      string refnameRaa15 = mkAxisCode(13,1,1);
       const Estimate1D& refdataRaa15 =refData(refnameRaa15);
       book(hRaaNpart["Raa_pt1020c05_AuAu"], refnameRaa15, refdataRaa15);

      string refnameRaa16 = mkAxisCode(14,1,1);
       const Estimate1D& refdataRaa16 =refData(refnameRaa16);
       book(hRaaNpart["Raa_pt1020c010_AuAu"], refnameRaa16, refdataRaa16);
       */

// // Mani: Is this correct?? Is there no need for 'all' centralities?
      /*
       centBins.insert(pair<string, int>("Raa_c010_AuAu",0));
       centBins.insert(pair<string, int>("Raa_c1020_AuAu",1));
       centBins.insert(pair<string, int>("Raa_c2030_AuAu",2));
       centBins.insert(pair<string, int>("Raa_c3040_AuAu",3));
       centBins.insert(pair<string, int>("Raa_c4050_AuAu",4));
       centBins.insert(pair<string, int>("Raa_c5060_AuAu",5));
       centBins.insert(pair<string, int>("Raa_c6070_AuAu",6));
       centBins.insert(pair<string, int>("Raa_c7080_AuAu",7));
       centBins.insert(pair<string, int>("Raa_c8092_AuAu",8));
       centBins.insert(pair<string, int>("Raa_call_AuAu",9));
       */
    }


    void analyze(const Event& event) {

        Particles neutralParticles = apply<UnstableParticles>(event,"up").particles();

        /*const ParticlePair& beam = beams();
        string CollSystem = "Empty";

        if (beamOpt == "NONE") {

          if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
          {
                  CollSystem = "AuAu";
          }
          if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
          {
                  CollSystem = "pp";
          }
        }

        else if (beamOpt == "PP200") CollSystem = "pp";
	else if (beamOpt == "AUAU200") CollSystem = "AuAu";*/

        if(collSys == pp)
        {
                 _c["sow_pp"]->fill();
                 for(Particle p : neutralParticles)
                 {
                     _h["c5Pt_pp"]->fill(p.pT()/GeV);
                     _h["c10Pt_pp"]->fill(p.pT()/GeV);
                     _h["c20Pt_pp"]->fill(p.pT()/GeV);
                     _h["c30Pt_pp"]->fill(p.pT()/GeV);
                     _h["c40Pt_pp"]->fill(p.pT()/GeV);
                     _h["c50Pt_pp"]->fill(p.pT()/GeV);
                     _h["c60Pt_pp"]->fill(p.pT()/GeV);
                     _h["c70Pt_pp"]->fill(p.pT()/GeV);
                     _h["c80Pt_pp"]->fill(p.pT()/GeV);
                     _h["c92Pt_pp"]->fill(p.pT()/GeV);
                     _h["callPt_pp"]->fill(p.pT()/GeV);
                 }
                 return;
        }

        const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
        const double c = cent();

        if ((c < 0.) || (c >= 92.)) vetoEvent;

        if (collSys == AuAu)
        {
            _c["sow_AuAu_call"]->fill();

            if((c >= 0.) && (c < 5.))
            {
                _c["sow_AuAu_c5"]->fill();
                for(const Particle& p : neutralParticles)
                {
                    double partPt = p.pT()/GeV;
                    double pt_weight = 1./(partPt*2.*M_PI);

                    _h["ptyieldsc5"]->fill(partPt, pt_weight);
                    _h["c5Pt_AuAu"]->fill(p.pT()/GeV);
                }
            }
            else if((c >= 0.) && (c < 10.))
            {
                _c["sow_AuAu_c10"]->fill();
                for(const Particle& p : neutralParticles)
                {
                    double partPt = p.pT()/GeV;
                    double pt_weight = 1./(partPt*2.*M_PI);
                    _h["ptyieldsc10"]->fill(partPt, pt_weight);
                    _h["c10Pt_AuAu"]->fill(p.pT()/GeV);
                    _h["callPt_AuAu"]->fill(p.pT()/GeV);
                    _h["ptyieldscall"]->fill(p.pT()/GeV, pt_weight);
                }
            }
            else if((c >= 10.) && (c < 20.))
            {
                _c["sow_AuAu_c20"]->fill();
                for(const Particle& p : neutralParticles)
                {
                    double partPt = p.pT()/GeV;
                    double pt_weight = 1./(partPt*2.*M_PI);
                    _h["ptyieldsc20"]->fill(partPt, pt_weight);
                    _h["c20Pt_AuAu"]->fill(p.pT()/GeV);
                    _h["callPt_AuAu"]->fill(p.pT()/GeV);
                    _h["ptyieldscall"]->fill(p.pT()/GeV, pt_weight);
                }
            }
            else if((c >= 20.) && (c < 30.))
            {
                _c["sow_AuAu_c30"]->fill();
                for(const Particle& p : neutralParticles)
                {
                    double partPt = p.pT()/GeV;
                    double pt_weight = 1./(partPt*2.*M_PI);
                    _h["ptyieldsc30"]->fill(partPt, pt_weight);
                    _h["c30Pt_AuAu"]->fill(p.pT()/GeV);
                    _h["callPt_AuAu"]->fill(p.pT()/GeV);
                    _h["ptyieldscall"]->fill(p.pT()/GeV, pt_weight);
                }
            }
            else if((c >= 30.) && (c < 40.))
            {
                _c["sow_AuAu_c40"]->fill();
                for(const Particle& p : neutralParticles)
                {
                    double partPt = p.pT()/GeV;
                    double pt_weight = 1./(partPt*2.*M_PI);
                    _h["ptyieldsc40"]->fill(partPt, pt_weight);
                    _h["c40Pt_AuAu"]->fill(p.pT()/GeV);
                    _h["callPt_AuAu"]->fill(p.pT()/GeV);
                    _h["ptyieldscall"]->fill(p.pT()/GeV, pt_weight);
                }
            }
            else if((c >= 40.) && (c < 50.))
            {
                _c["sow_AuAu_c50"]->fill();
                for(const Particle& p : neutralParticles)
                {
                    double partPt = p.pT()/GeV;
                    double pt_weight = 1./(partPt*2.*M_PI);
                    _h["ptyieldsc50"]->fill(partPt, pt_weight);
                    _h["c50Pt_AuAu"]->fill(p.pT()/GeV);
                    _h["callPt_AuAu"]->fill(p.pT()/GeV);
                    _h["ptyieldscall"]->fill(p.pT()/GeV, pt_weight);
                }
            }
            else if((c >= 50.) && (c < 60.))
            {
                _c["sow_AuAu_c60"]->fill();
                for(const Particle& p : neutralParticles)
                {
                    double partPt = p.pT()/GeV;
                    double pt_weight = 1./(partPt*2.*M_PI);
                    _h["ptyieldsc60"]->fill(partPt, pt_weight);
                    _h["c60Pt_AuAu"]->fill(p.pT()/GeV);
                    _h["callPt_AuAu"]->fill(p.pT()/GeV);
                    _h["ptyieldscall"]->fill(p.pT()/GeV, pt_weight);
                }
            }
            else if((c >= 60.) && (c < 70.))
            {
                _c["sow_AuAu_c70"]->fill();
                for(const Particle& p : neutralParticles)
                {
                    double partPt = p.pT()/GeV;
                    double pt_weight = 1./(partPt*2.*M_PI);
                    _h["ptyieldsc70"]->fill(partPt, pt_weight);
                    _h["c70Pt_AuAu"]->fill(p.pT()/GeV);
                    _h["callPt_AuAu"]->fill(p.pT()/GeV);
                    _h["ptyieldscall"]->fill(p.pT()/GeV, pt_weight);
                }
            }
            else if((c >= 70.) && (c < 80.))
            {
                _c["sow_AuAu_c80"]->fill();
                for(const Particle& p : neutralParticles)
                {
                    double partPt = p.pT()/GeV;
                    double pt_weight = 1./(partPt*2.*M_PI);
                    _h["ptyieldsc80"]->fill(partPt, pt_weight);
                    _h["c80Pt_AuAu"]->fill(p.pT()/GeV);
                    _h["callPt_AuAu"]->fill(p.pT()/GeV);
                    _h["ptyieldscall"]->fill(p.pT()/GeV, pt_weight);
                }
            }
 	        else if((c >= 80.) && (c < 92.))
            {
                _c["sow_AuAu_c92"]->fill();
                for(const Particle& p : neutralParticles)
                {
                    double partPt = p.pT()/GeV;
                    double pt_weight = 1./(partPt*2.*M_PI);
                    _h["ptyieldsc92"]->fill(partPt, pt_weight);
                    _h["c92Pt_AuAu"]->fill(p.pT()/GeV);
                    _h["callPt_AuAu"]->fill(p.pT()/GeV);
                    _h["ptyieldscall"]->fill(p.pT()/GeV, pt_weight);
                }
            }

 	  }
    }

    void finalize() {

// Yields, don't forget the Yields ya wanker.

// R_AA for 0_5% Au+Au
    // Normalize Au+Au yield
        binShift(*_h["c5Pt_AuAu"]);
        if (_c["sow_AuAu_c5"]->sumW() > 0) {
        _h["c5Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c5"]->sumW());
        }
    // Normalize p+p yield
        binShift(*_h["c5Pt_pp"]);
        if (_c["sow_pp"]->sumW() > 0) {
        _h["c5Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
        }
    // Only divide if BOTH are valid
        if (_c["sow_AuAu_c5"]->sumW() > 0 && _c["sow_pp"]->sumW() > 0) {
        divide(_h["c5Pt_AuAu"],_h["c5Pt_pp"],hRaa["Raa_c05_AuAu"]);
    // Scale by <T_AA>
        hRaa["Raa_c05_AuAu"]->scale(1./1065.4);
        }

// R_AA for 0_10% Au+Au
    // Normalize Au+Au yield
        binShift(*_h["c10Pt_AuAu"]);
        if (_c["sow_AuAu_c10"]->sumW() > 0) {
        _h["c10Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c10"]->sumW());
        }
    // Normalize p+p yield
        binShift(*_h["c10Pt_pp"]);
        if (_c["sow_pp"]->sumW() > 0) {
        _h["c10Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
        }
    // Only divide if BOTH are valid
        if (_c["sow_AuAu_c10"]->sumW() > 0 && _c["sow_pp"]->sumW() > 0) {
        divide(_h["c10Pt_AuAu"],_h["c10Pt_pp"],hRaa["Raa_c010_AuAu"]);
    // Scale by <T_AA>
        hRaa["Raa_c010_AuAu"]->scale(1./955.4);
        }

// R_AA for 10_20% Au+Au
    // Normalize Au+Au yield
        binShift(*_h["c20Pt_AuAu"]);
        if (_c["sow_AuAu_c20"]->sumW() > 0) {
        _h["c20Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c20"]->sumW());
        }
    // Normalize p+p yield
        binShift(*_h["c20Pt_pp"]);
        if (_c["sow_pp"]->sumW() > 0) {
        _h["c20Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
        }
    // Only divide if BOTH are valid
        if (_c["sow_AuAu_c20"]->sumW() > 0 && _c["sow_pp"]->sumW() > 0) {
        divide(_h["c20Pt_AuAu"],_h["c20Pt_pp"],hRaa["Raa_c1020_AuAu"]);
    // Scale by <T_AA>
        hRaa["Raa_c1020_AuAu"]->scale(1./602.6);
        }

// R_AA for 20_30% Au+Au
    // Normalize Au+Au yield
        binShift(*_h["c30Pt_AuAu"]);
        if (_c["sow_AuAu_c30"]->sumW() > 0) {
        _h["c30Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c30"]->sumW());
        }
    // Normalize p+p yield
        binShift(*_h["c30Pt_pp"]);
        if (_c["sow_pp"]->sumW() > 0) {
        _h["c30Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
        }
    // Only divide if BOTH are valid
        if (_c["sow_AuAu_c30"]->sumW() > 0 && _c["sow_pp"]->sumW() > 0) {
        divide(_h["c30Pt_AuAu"],_h["c30Pt_pp"],hRaa["Raa_c2030_AuAu"]);
    // Scale by <T_AA>
        hRaa["Raa_c2030_AuAu"]->scale(1./373.8);
        }

// R_AA for 30_40% Au+Au
    // Normalize Au+Au yield
        binShift(*_h["c40Pt_AuAu"]);
        if (_c["sow_AuAu_c40"]->sumW() > 0) {
        _h["c40Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c40"]->sumW());
        }
    // Normalize p+p yield
        binShift(*_h["c40Pt_pp"]);
        if (_c["sow_pp"]->sumW() > 0) {
        _h["c40Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
        }
    // Only divide if BOTH are valid
        if (_c["sow_AuAu_c40"]->sumW() > 0 && _c["sow_pp"]->sumW() > 0) {
        divide(_h["c40Pt_AuAu"],_h["c40Pt_pp"],hRaa["Raa_c3040_AuAu"]);
    // Scale by <T_AA>
        hRaa["Raa_c3040_AuAu"]->scale(1./219.8);
        }

// R_AA for 40_50% Au+Au
    // Normalize Au+Au yield
        binShift(*_h["c50Pt_AuAu"]);
        if (_c["sow_AuAu_c50"]->sumW() > 0) {
        _h["c50Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c50"]->sumW());
        }
    // Normalize p+p yield
        binShift(*_h["c50Pt_pp"]);
        if (_c["sow_pp"]->sumW() > 0) {
        _h["c50Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
        }
    // Only divide if BOTH are valid
        if (_c["sow_AuAu_c50"]->sumW() > 0 && _c["sow_pp"]->sumW() > 0) {
        divide(_h["c50Pt_AuAu"],_h["c50Pt_pp"],hRaa["Raa_c4050_AuAu"]);
    // Scale by <T_AA>
        hRaa["Raa_c4050_AuAu"]->scale(1./120.3);
        }

// R_AA for 50_60% Au+Au
    // Normalize Au+Au yield
        binShift(*_h["c60Pt_AuAu"]);
        if (_c["sow_AuAu_c60"]->sumW() > 0) {
        _h["c60Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c60"]->sumW());
        }
    // Normalize p+p yield
        binShift(*_h["c60Pt_pp"]);
        if (_c["sow_pp"]->sumW() > 0) {
        _h["c60Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
        }
    // Only divide if BOTH are valid
        if (_c["sow_AuAu_c60"]->sumW() > 0 && _c["sow_pp"]->sumW() > 0) {
        divide(_h["c60Pt_AuAu"],_h["c60Pt_pp"],hRaa["Raa_c5060_AuAu"]);
    // Scale by <T_AA>
        hRaa["Raa_c5060_AuAu"]->scale(1./61.0);
        }

// R_AA for 60_70% Au+Au
    // Normalize Au+Au yield
        binShift(*_h["c70Pt_AuAu"]);
        if (_c["sow_AuAu_c70"]->sumW() > 0) {
        _h["c70Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c70"]->sumW());
        }
    // Normalize p+p yield
        binShift(*_h["c70Pt_pp"]);
        if (_c["sow_pp"]->sumW() > 0) {
        _h["c70Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
        }
    // Only divide if BOTH are valid
        if (_c["sow_AuAu_c70"]->sumW() > 0 && _c["sow_pp"]->sumW() > 0) {
        divide(_h["c70Pt_AuAu"],_h["c70Pt_pp"],hRaa["Raa_c6070_AuAu"]);
    // Scale by <T_AA>
        hRaa["Raa_c6070_AuAu"]->scale(1./28.5);
        }

// R_AA for 70_80% Au+Au
    // Normalize Au+Au yield
        binShift(*_h["c80Pt_AuAu"]);
        if (_c["sow_AuAu_c80"]->sumW() > 0) {
        _h["c80Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c80"]->sumW());
        }
    // Normalize p+p yield
        binShift(*_h["c80Pt_pp"]);
        if (_c["sow_pp"]->sumW() > 0) {
        _h["c80Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
        }
    // Only divide if BOTH are valid
        if (_c["sow_AuAu_c80"]->sumW() > 0 && _c["sow_pp"]->sumW() > 0) {
        divide(_h["c80Pt_AuAu"],_h["c80Pt_pp"],hRaa["Raa_c7080_AuAu"]);
    // Scale by <T_AA>
        hRaa["Raa_c7080_AuAu"]->scale(1./12.4);
        }

// R_AA for 80_92% Au+Au
    // Normalize Au+Au yield
        binShift(*_h["c92Pt_AuAu"]);
        if (_c["sow_AuAu_c92"]->sumW() > 0) {
        _h["c92Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c92"]->sumW());
        }
    // Normalize p+p yield
        binShift(*_h["c92Pt_pp"]);
        if (_c["sow_pp"]->sumW() > 0) {
        _h["c92Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
        }
    // Only divide if BOTH are valid
        if (_c["sow_AuAu_c92"]->sumW() > 0 && _c["sow_pp"]->sumW() > 0) {
        divide(_h["c92Pt_AuAu"],_h["c92Pt_pp"],hRaa["Raa_c8092_AuAu"]);
    // Scale by <T_AA>
        hRaa["Raa_c8092_AuAu"]->scale(1./4.9);
        }

// R_AA for 0_92% min_bias Au+Au
    // Normalize Au+Au yield
        binShift(*_h["callPt_AuAu"]);
        if (_c["sow_AuAu_call"]->sumW() > 0) {
        _h["callPt_AuAu"]->scaleW(1./_c["sow_AuAu_call"]->sumW());
        }
    // Normalize p+p yield
        binShift(*_h["callPt_pp"]);
        if (_c["sow_pp"]->sumW() > 0) {
        _h["callPt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
        }
    // Only divide if BOTH are valid
        if (_c["sow_AuAu_call"]->sumW() > 0 && _c["sow_pp"]->sumW() > 0) {
        divide(_h["callPt_AuAu"],_h["callPt_pp"],hRaa["Raa_minbias_AuAu"]);
    // Scale by <T_AA>
        hRaa["Raa_minbias_AuAu"]->scale(1./257.8);
        }


// // Scale???
// // Directly copy/pasted?
       /*
        * ANTONIO'S COMMENT: I think it is better to scale the pp (or AuAu) instead of scaling the Raa plot ****************
        * We can do this when we have the Ncoll
       hRaa["Raa_c05_AuAu"]->scale(1./777.2);
       hRaa["Raa_c010_AuAu"]->scale(1./777.2);
       hRaa["Raa_c1020_AuAu"]->scale(1./496.7);
       hRaa["Raa_c2030_AuAu"]->scale(1./253.6);
       hRaa["Raa_c3040_AuAu"]->scale(1./81.81);
       hRaa["Raa_c4050_AuAu"]->scale(1./13.88);
       hRaa["Raa_c5060_AuAu"]->scale(1./13.88);
       hRaa["Raa_c6070_AuAu"]->scale(1./13.88);
       hRaa["Raa_c7080_AuAu"]->scale(1./13.88);
       hRaa["Raa_c8092_AuAu"]->scale(1./13.88);
       hRaa["Raa_minbias_AuAu"]->scale(1./401.6);
*/
    //Centrality vs RAA_______________________________

       /*
       for(auto element : hRaa)
       {
           string name = element.first;
           if(name.find("Raa_minbias_AuAu") != std::string::npos) continue;

           YODA::Estimate1D yodaRaa = *element.second;

           double averageRaaPt520 = 0.;
           double averageRaaPt1020 = 0.;
           double averageRaaPt520Error = 0.;
           double averageRaaPt1020Error = 0.;
           int nbinsPt520 = 0;
           int nbinsPt1020 = 0;




           for(auto &point : yodaRaa.points())
           {

               if(point.x() > 5. && point.x() < 20.)
               {
                   if(!isnan(point.y()))
                   {
                       averageRaaPt520 += point.y();
                       averageRaaPt520Error += point.yErrPlus()*point.yErrPlus();
                       nbinsPt520++;
                   }

               }
               if(point.x() > 10. && point.x() < 20)
               {
                   if(!isnan(point.y()))
                   {
                       averageRaaPt1020 += point.y();
                       averageRaaPt1020Error += point.yErrPlus()*point.yErrPlus();
                       nbinsPt1020++;
                   }
               }
           }


           if(name.find("39") != std::string::npos)
           {
               hRaaNpart["Raa_pt46_AuAu39"]->point(centBins[name]).setY(averageRaaPt46/nbinsPt46);
               hRaaNpart["Raa_pt46_AuAu39"]->point(centBins[name]).setYErr(sqrt(averageRaaPt46Error)/nbinsPt46);
               hRaaNpart["Raa_pt610_AuAu39"]->point(centBins[name]).setY(averageRaaPt610/nbinsPt610);
               hRaaNpart["Raa_pt610_AuAu39"]->point(centBins[name]).setYErr(sqrt(averageRaaPt610Error)/nbinsPt610);
           }
           else if(name.find("62") != std::string::npos)
           {
               hRaaNpart["Raa_pt46_AuAu62"]->point(centBins[name]).setY(averageRaaPt46/nbinsPt46);
               hRaaNpart["Raa_pt46_AuAu62"]->point(centBins[name]).setYErr(sqrt(averageRaaPt46Error)/nbinsPt46);
              hRaaNpart["Raa_pt610_AuAu62"]->point(centBins[name]).setY(averageRaaPt610/nbinsPt610);
               hRaaNpart["Raa_pt610_AuAu62"]->point(centBins[name]).setYErr(sqrt(averageRaaPt610Error)/nbinsPt610);
           }

       }*/

    }

    map<string, Histo1DPtr> _h;

    map<string, Estimate1DPtr> hRaa;
    map<string, CounterPtr> _c;

     map<string, Estimate1DPtr> hRaaNpart;
     //map<string, int> centBins;
    string beamOpt;
    enum CollisionSystem {pp, AuAu};
    CollisionSystem collSys;

  };



  RIVET_DECLARE_PLUGIN(PHENIX_2008_I778168);

}
