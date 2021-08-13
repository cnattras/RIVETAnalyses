// -*- C++ -*-
#include "Rivet/Analysis.hh"

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/AliceCommon.hh"

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
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
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2008_I778168);


    void init() {

      std::initializer_list<int> pdgIds = {111};  // Pion 0

      //const PrimaryParticles fs(pdgIds, Cuts::abseta < 1. && Cuts::abscharge == 0);
      //declare(fs, "fs");

      const UnstableParticles up(Cuts::abseta < 0.35 && Cuts::pid == 111);
      declare(up, "up");

      beamOpt = getOption<string>("beam","NONE");


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
       string refnameRaa1 = mkAxisCode(7,1,3);
             const Scatter2D& refdataRaa1 =refData(refnameRaa1);
       book(_h["c5Pt_AuAu"], refnameRaa1 + "_AuAu", refdataRaa1);
       book(_h["c5Pt_pp"], refnameRaa1 + "_pp", refdataRaa1);
       book(hRaa["Raa_c05_AuAu"], refnameRaa1);

       string refnameRaa2 = mkAxisCode(6,1,1);
             const Scatter2D& refdataRaa2 =refData(refnameRaa2);
       book(_h["c10Pt_AuAu"], refnameRaa2 + "_AuAu", refdataRaa2);
       book(_h["c10Pt_pp"], refnameRaa2 + "_pp", refdataRaa2);
       book(hRaa["Raa_c010_AuAu"], refnameRaa2);

       string refnameRaa3 = mkAxisCode(6,1,2);
             const Scatter2D& refdataRaa3 =refData(refnameRaa3);
       book(_h["c20Pt_AuAu"], refnameRaa3 + "_AuAu", refdataRaa3);
       book(_h["c20Pt_pp"], refnameRaa3 + "_pp", refdataRaa3);
       book(hRaa["Raa_c1020_AuAu"], refnameRaa3);

       string refnameRaa4 = mkAxisCode(7,1,1);
             const Scatter2D& refdataRaa4 =refData(refnameRaa4);
       book(_h["c30Pt_AuAu"], refnameRaa4 + "_AuAu", refdataRaa4);
       book(_h["c30Pt_pp"], refnameRaa4 + "_pp", refdataRaa4);
       book(hRaa["Raa_c2030_AuAu"], refnameRaa4);

       string refnameRaa5 = mkAxisCode(7,1,2);
             const Scatter2D& refdataRaa5 =refData(refnameRaa5);
       book(_h["c40Pt_AuAu"], refnameRaa5 + "_AuAu", refdataRaa5);
       book(_h["c40Pt_pp"], refnameRaa5 + "_pp", refdataRaa5);
       book(hRaa["Raa_c3040_AuAu"], refnameRaa5);

       string refnameRaa6 = mkAxisCode(8,1,1);
             const Scatter2D& refdataRaa6 =refData(refnameRaa6);
       book(_h["c50Pt_AuAu"], refnameRaa6 + "_AuAu", refdataRaa6);
       book(_h["c50Pt_pp"], refnameRaa6 + "_pp", refdataRaa6);
       book(hRaa["Raa_c4050_AuAu"], refnameRaa6);

       string refnameRaa7 = mkAxisCode(8,1,2);
             const Scatter2D& refdataRaa7 =refData(refnameRaa7);
       book(_h["c60Pt_AuAu"], refnameRaa7 + "_AuAu", refdataRaa7);
       book(_h["c60Pt_pp"], refnameRaa7 + "_pp", refdataRaa7);
       book(hRaa["Raa_c5060_AuAu"], refnameRaa7);

       string refnameRaa8 = mkAxisCode(9,1,1);
             const Scatter2D& refdataRaa8 =refData(refnameRaa8);
       book(_h["c70Pt_AuAu"], refnameRaa8 + "_AuAu", refdataRaa8);
       book(_h["c70Pt_pp"], refnameRaa8 + "_pp", refdataRaa8);
       book(hRaa["Raa_c6070_AuAu"], refnameRaa8);

       string refnameRaa9 = mkAxisCode(10,1,1);
             const Scatter2D& refdataRaa9 =refData(refnameRaa9);
       book(_h["c80Pt_AuAu"], refnameRaa9 + "_AuAu", refdataRaa9);
       book(_h["c80Pt_pp"], refnameRaa9 + "_pp", refdataRaa9);
       book(hRaa["Raa_c7080_AuAu"], refnameRaa9);

       string refnameRaa10 = mkAxisCode(10,1,2);
             const Scatter2D& refdataRaa10 =refData(refnameRaa10);
       book(_h["c92Pt_AuAu"], refnameRaa10 + "_AuAu", refdataRaa10);
       book(_h["c92Pt_pp"], refnameRaa10 + "_pp", refdataRaa10);
       book(hRaa["Raa_c8092_AuAu"], refnameRaa10);

       string refnameRaa11 = mkAxisCode(6,1,3);
             const Scatter2D& refdataRaa11 =refData(refnameRaa11);
       book(_h["callPt_AuAu"], refnameRaa11 + "_AuAu", refdataRaa11);
       book(_h["callPt_pp"], refnameRaa11 + "_pp", refdataRaa11);
       book(hRaa["Raa_minbias_AuAu"], refnameRaa11);


       book(_c["sow_pp"],"_sow_pp");

 //Centrality vs RAA
       /*
       string refnameRaa13 = mkAxisCode(11,1,1);
       const Scatter2D& refdataRaa13 =refData(refnameRaa13);
       book(hRaaNpart["Raa_pt520c05_AuAu"], refnameRaa13, refdataRaa13);

      string refnameRaa14 = mkAxisCode(12,1,1);
       const Scatter2D& refdataRaa14 =refData(refnameRaa14);
       book(hRaaNpart["Raa_pt520c010_AuAu"], refnameRaa14, refdataRaa14);

      string refnameRaa15 = mkAxisCode(13,1,1);
       const Scatter2D& refdataRaa15 =refData(refnameRaa15);
       book(hRaaNpart["Raa_pt1020c05_AuAu"], refnameRaa15, refdataRaa15);

      string refnameRaa16 = mkAxisCode(14,1,1);
       const Scatter2D& refdataRaa16 =refData(refnameRaa16);
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

        Particles neutralParticles = applyProjection<UnstableParticles>(event,"up").particles();

        const ParticlePair& beam = beams();
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
	else if (beamOpt == "AUAU200") CollSystem = "AuAu";

        if(CollSystem == "pp")
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

        if (CollSystem == "AuAu")
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

 //Yields_________________
       _h["ptyieldsc5"]->scaleW(1./_c["sow_AuAu_c5"]->sumW());
       _h["ptyieldsc10"]->scaleW(1./_c["sow_AuAu_c10"]->sumW());
       _h["ptyieldsc20"]->scaleW(1./_c["sow_AuAu_c20"]->sumW());
       _h["ptyieldsc30"]->scaleW(1./_c["sow_AuAu_c30"]->sumW());
       _h["ptyieldsc40"]->scaleW(1./_c["sow_AuAu_c40"]->sumW());
       _h["ptyieldsc50"]->scaleW(1./_c["sow_AuAu_c50"]->sumW());
       _h["ptyieldsc60"]->scaleW(1./_c["sow_AuAu_c60"]->sumW());
       _h["ptyieldsc70"]->scaleW(1./_c["sow_AuAu_c70"]->sumW());
       _h["ptyieldsc80"]->scaleW(1./_c["sow_AuAu_c80"]->sumW());
       _h["ptyieldsc92"]->scaleW(1./_c["sow_AuAu_c92"]->sumW());
       _h["ptyieldscall"]->scaleW(1./_c["sow_AuAu_call"]->sumW());

 //RAA _______________________________
       _h["c5Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c5"]->sumW());
       _h["c10Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c10"]->sumW());
       _h["c20Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c20"]->sumW());
       _h["c30Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c30"]->sumW());
       _h["c40Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c40"]->sumW());
       _h["c50Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c50"]->sumW());
       _h["c60Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c60"]->sumW());
       _h["c70Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c70"]->sumW());
       _h["c80Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c80"]->sumW());
       _h["c92Pt_AuAu"]->scaleW(1./_c["sow_AuAu_c92"]->sumW());
       _h["callPt_AuAu"]->scaleW(1./_c["sow_AuAu_call"]->sumW());

       _h["c5Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
       _h["c10Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
       _h["c20Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
       _h["c30Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
       _h["c40Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
       _h["c50Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
       _h["c60Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
       _h["c70Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
       _h["c80Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
       _h["c92Pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());
       _h["callPt_pp"]->scaleW(1./_c["sow_pp"]->sumW());

       divide(_h["c5Pt_AuAu"],_h["c5Pt_pp"],hRaa["Raa_c05_AuAu"]);
       divide(_h["c10Pt_AuAu"],_h["c10Pt_pp"],hRaa["Raa_c010_AuAu"]);
       divide(_h["c20Pt_AuAu"],_h["c20Pt_pp"],hRaa["Raa_c1020_AuAu"]);
       divide(_h["c30Pt_AuAu"],_h["c30Pt_pp"],hRaa["Raa_c2030_AuAu"]);
       divide(_h["c40Pt_AuAu"],_h["c40Pt_pp"],hRaa["Raa_c3040_AuAu"]);
       divide(_h["c50Pt_AuAu"],_h["c50Pt_pp"],hRaa["Raa_c4050_AuAu"]);
       divide(_h["c60Pt_AuAu"],_h["c60Pt_pp"],hRaa["Raa_c5060_AuAu"]);
       divide(_h["c70Pt_AuAu"],_h["c70Pt_pp"],hRaa["Raa_c6070_AuAu"]);
       divide(_h["c80Pt_AuAu"],_h["c80Pt_pp"],hRaa["Raa_c7080_AuAu"]);
       divide(_h["c92Pt_AuAu"],_h["c92Pt_pp"],hRaa["Raa_c8092_AuAu"]);
       divide(_h["callPt_AuAu"],_h["callPt_pp"],hRaa["Raa_minbias_AuAu"]);

       hRaa["Raa_c05_AuAu"]->scaleY(1./1065.4);
       hRaa["Raa_c010_AuAu"]->scaleY(1./955.4);
       hRaa["Raa_c1020_AuAu"]->scaleY(1./602.6);
       hRaa["Raa_c2030_AuAu"]->scaleY(1./373.8);
       hRaa["Raa_c3040_AuAu"]->scaleY(1./219.8);
       hRaa["Raa_c4050_AuAu"]->scaleY(1./120.3);
       hRaa["Raa_c5060_AuAu"]->scaleY(1./61.0);
       hRaa["Raa_c6070_AuAu"]->scaleY(1./28.5);
       hRaa["Raa_c7080_AuAu"]->scaleY(1./12.4);
       hRaa["Raa_c8092_AuAu"]->scaleY(1./4.9);
       hRaa["Raa_minbias_AuAu"]->scaleY(1./257.8);

// // Scale???
// // Directly copy/pasted?
       /*
        * ANTONIO'S COMMENT: I think it is better to scale the pp (or AuAu) instead of scaling the Raa plot ****************
        * We can do this when we have the Ncoll
       hRaa["Raa_c05_AuAu"]->scaleY(1./777.2);
       hRaa["Raa_c010_AuAu"]->scaleY(1./777.2);
       hRaa["Raa_c1020_AuAu"]->scaleY(1./496.7);
       hRaa["Raa_c2030_AuAu"]->scaleY(1./253.6);
       hRaa["Raa_c3040_AuAu"]->scaleY(1./81.81);
       hRaa["Raa_c4050_AuAu"]->scaleY(1./13.88);
       hRaa["Raa_c5060_AuAu"]->scaleY(1./13.88);
       hRaa["Raa_c6070_AuAu"]->scaleY(1./13.88);
       hRaa["Raa_c7080_AuAu"]->scaleY(1./13.88);
       hRaa["Raa_c8092_AuAu"]->scaleY(1./13.88);
       hRaa["Raa_minbias_AuAu"]->scaleY(1./401.6);
*/
    //Centrality vs RAA_______________________________

       /*
       for(auto element : hRaa)
       {
           string name = element.first;
           if(name.find("Raa_minbias_AuAu") != std::string::npos) continue;

           YODA::Scatter2D yodaRaa = *element.second;

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

    map<string, Scatter2DPtr> hRaa;
    map<string, CounterPtr> _c;

     map<string, Scatter2DPtr> hRaaNpart;
     //map<string, int> centBins;
    string beamOpt = "NONE";
    //enum CollisionSystem {pp, AuAu};
    //CollisionSystem collSys;

  };



  DECLARE_RIVET_PLUGIN(PHENIX_2008_I778168);

}
