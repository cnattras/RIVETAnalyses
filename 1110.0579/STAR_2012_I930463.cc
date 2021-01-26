// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "../Centralities/RHICCentrality.hh" //external header for Centrality calculation
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES

namespace Rivet {


  /// @brief Add a short analysis description here
  class STAR_2012_I930463 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2012_I930463);

  bool getBinCenter(YODA::Histo1D hist, double pT, double &binCenter)
  {
      if(pT > hist.xMin() && pT < hist.xMax())
      {
          binCenter = hist.bin(hist.binIndexAt(pT)).xMid();
          return true;
      }
      else return false;
  }


	void init() {
		std::initializer_list<int> pdgIds = { 321, 211, 2212};  // pi+ 211  K+ 321   proton 2212	K0S 310		rho0 113

		//charged particles
		const PrimaryParticles cp(pdgIds, Cuts::absrap < 0.5 && Cuts::pT > 3.5*GeV && Cuts::abscharge > 0);
		declare(cp, "cp");

		//neutral particles (changed fs to np)
		//const PrimaryParticles np(pdgIds, Cuts::absrap < 0.5 && Cuts::abscharge == 0);
		//declare(np, "np");

    const UnstableParticles np(Cuts::absrap < 0.5 && Cuts::pT > 3.5*GeV && (Cuts::abspid == 310 || Cuts::abspid == 113) );
    declare(np, "np");


		declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");


		book(sow["sow_pp"], "sow_pp");
		book(sow["sow_AuAuc12"], "sow_AuAuc12");

		//Figure 1 Yield pp
		book(hPionPosPt["ptyieldspp"], 1, 1, 1);
		book(hPionNegPt["ptyieldspp"], 1, 1, 2);
		book(hKaonPosPt["ptyieldspp"], 1, 1, 3);
		book(hKaonNegPt["ptyieldspp"], 1, 1, 4);
		book(hProtPosPt["ptyieldspp"], 1, 1, 5);
		book(hProtNegPt["ptyieldspp"], 1, 1, 6);
		book(hKaon0SPt["ptyieldspp"], 2, 1, 1);
		book(hRho0Pt["ptyieldspp"], 3, 1, 1);

		//Figure 1 Yield AUAU
		book(hKpPosPt["ptyieldsAuAuc12"], 4, 1, 1);
		book(hKpNegPt["ptyieldsAuAuc12"], 4, 1, 2);
		book(hRho0Pt["ptyieldsAuAuc12"], 5, 1, 1);
		book(hKaon0SPt["ptyieldsAuAuc12"], 6, 1, 1);

		//Figure 2 Yield Ratio pp
		string refname1 = mkAxisCode(7, 1, 1);
		const Scatter2D& refdata1 = refData(refname1);
		book(hPionNegPt["pp1"], refname1 + "_PionNeg", refdata1);
		book(hPionPosPt["pp1"], refname1 + "_PionPos", refdata1);
		book(RatioPion["pp"], refname1);

		string refname2 = mkAxisCode(7, 1, 2);
		const Scatter2D& refdata2 = refData(refname2);
		book(hProtNegPt["pp1"], refname2 + "_ProtNeg", refdata2);
		book(hProtPosPt["pp1"], refname2 + "_ProtPos", refdata2);
		book(RatioProt["pp"], refname2);

		string refname3 = mkAxisCode(7, 1, 3);
		const Scatter2D& refdata3 = refData(refname3);
		book(hKaonNegPt["pp"], refname3 + "_KaonNeg", refdata3);
		book(hKaonPosPt["pp"], refname3 + "_KaonPos", refdata3);
		book(RatioKaon["pp"], refname3);

		string refname4 = mkAxisCode(7, 1, 4);
		const Scatter2D& refdata4 = refData(refname4);
		book(hProtPosPt["pp2"], refname4 + "_ProtPos", refdata4);
		book(hPionPosPt["pp2"], refname4 + "_PionPos", refdata4);
		book(Ratioppipos["pp"], refname4);

		string refname5 = mkAxisCode(7, 1, 5);
		const Scatter2D& refdata5 = refData(refname5);
		book(hProtNegPt["pp2"], refname5 + "_ProtNeg", refdata5);
		book(hPionNegPt["pp2"], refname5 + "_PionNeg", refdata5);
		book(Ratioppineg["pp"], refname5);

		string refname6 = mkAxisCode(7, 1, 6);
		const Scatter2D& refdata6 = refData(refname6);
		book(hKaonPt["pp"], refname6 + "_Kaon", refdata6);
		book(hPionPt["pp1"], refname6 + "_Pion", refdata6);
		book(RatioKpi["pp"], refname6);

		string refname7 = mkAxisCode(8, 1, 1);
		const Scatter2D& refdata7 = refData(refname7);
		book(hKaon0SPt["pp"], refname7 + "_Kaon0S", refdata7);
		book(hPionPt["pp2"], refname7 + "_Pion", refdata7);
		book(RatioK0spi["pp"], refname7);

		//Figure 2 Yield Ratio AUAU
		string refname8 = mkAxisCode(9, 1, 1);
		const Scatter2D& refdata8 = refData(refname8);
		book(hProtPosPt["AuAuc12"], refname8 + "_ProtPos", refdata8);
		book(hPionPosPt["AuAuc12"], refname8 + "_PionPos", refdata8);
		book(Ratioppipos["AuAuc12"], refname8);

		string refname9 = mkAxisCode(9, 1, 2);
		const Scatter2D& refdata9 = refData(refname9);
		book(hProtNegPt["AuAuc12"], refname9 + "_ProtNeg", refdata9);
		book(hPionNegPt["AuAuc12"], refname9 + "_PionNeg", refdata9);
		book(Ratioppineg["AuAuc12"], refname9);

		//Figure 3 RAA

		string refname10 = mkAxisCode(10, 1, 1);
		const Scatter2D& refdata10 = refData(refname10);
		book(hPionPt["Raa_c12_AuAu"], refname10 + "_AuAu", refdata10);
		book(hPionPt["Raa_c12_pp"], refname10 + "_pp", refdata10);
		book(hRaa["pi_c12_AuAu"], refname10);

		string refname11 = mkAxisCode(11, 1, 1);
		const Scatter2D& refdata11 = refData(refname11);
		book(hKpPt["Raa_c12_AuAu"], refname11 + "_AuAu", refdata11);
		book(hKpPt["Raa_c12_pp"], refname11 + "_pp", refdata11);
		book(hRaa["Kp_c12_AuAu"], refname11);

		string refname12 = mkAxisCode(12, 1, 1);
		const Scatter2D& refdata12 = refData(refname12);
		book(hKaon0SPt["Raa_c12_AuAu"], refname12 + "_AuAu", refdata12);
		book(hKaon0SPt["Raa_c12_pp"], refname12 + "_pp", refdata12);
		book(hRaa["K0S_c12_AuAu"], refname12);

		string refname13 = mkAxisCode(13, 1, 1);
		const Scatter2D& refdata13 = refData(refname13);
		book(hRho0Pt["Raa_c12_AuAu"], refname13 + "_AuAu", refdata13);
		book(hRho0Pt["Raa_c12_pp"], refname13 + "_pp", refdata13);
		book(hRaa["Rho0_c12_AuAu"], refname13);

		//Figure 3 RAA Ratio ?? Best way to do Raa ratio?


    }



    void analyze(const Event& event) {

      const ParticlePair& beam = beams();

      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970) collSys = AuAu200;
      else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = pp;

      Particles chargedParticles = applyProjection<PrimaryParticles>(event, "cp").particles();
  		Particles neutralParticles = applyProjection<UnstableParticles>(event, "np").particles();

      if (collSys == pp)
      {
          sow["sow_pp"]->fill();
          for (Particle p : chargedParticles)
          {
              double partPt = p.pT() / GeV;
              double pt_weight = 1. / (2. * M_PI);
              double binCenter = 0.;

              switch (p.pid()) {

                  case 211: // pi+
				          {
                      if(getBinCenter(*hPionPosPt["ptyieldspp"], partPt, binCenter))
                      {
                          pt_weight /= binCenter;
                          hPionPosPt["ptyieldspp"]->fill(partPt, pt_weight);
                      }
                      hPionPosPt["pp1"]->fill(partPt);
                      hPionPosPt["pp2"]->fill(partPt);
                      hPionPt["pp1"]->fill(partPt);
                      hPionPt["pp2"]->fill(partPt);
                      hPionPt["Raa_c12_pp"]->fill(partPt);
                      break;
                  }
                  case -211: // pi-
                  {
                      if(getBinCenter(*hPionNegPt["ptyieldspp"], partPt, binCenter))
                      {
                          pt_weight /= binCenter;
                          hPionNegPt["ptyieldspp"]->fill(partPt, pt_weight);
                      }
                      hPionNegPt["pp1"]->fill(partPt);
                      hPionNegPt["pp2"]->fill(partPt);
                      hPionPt["pp1"]->fill(partPt);
                      hPionPt["pp2"]->fill(partPt);
                      hPionPt["Raa_c12_pp"]->fill(partPt);
                      break;
                  }
                  case  321: // K+
                  {
                      if(getBinCenter(*hKaonPosPt["ptyieldspp"], partPt, binCenter))
                      {
                          pt_weight /= binCenter;
                          hKaonPosPt["ptyieldspp"]->fill(partPt, pt_weight);
                      }
                      hKaonPosPt["pp"]->fill(partPt);
                      hKaonPt["pp"]->fill(partPt);
                      break;
                  }
                  case  -321: // K-
                  {
                      if(getBinCenter(*hKaonNegPt["ptyieldspp"], partPt, binCenter))
                      {
                          pt_weight /= binCenter;
                          hKaonNegPt["ptyieldspp"]->fill(partPt, pt_weight);
                      }
                      hKaonNegPt["pp"]->fill(partPt);
                      hKaonPt["pp"]->fill(partPt);
                      hKpPt["Raa_c12_pp"]->fill(partPt);
                      break;
                  }
                  case 2212: // proton
                  {
                      if(getBinCenter(*hProtPosPt["ptyieldspp"], partPt, binCenter))
                      {
                          pt_weight /= binCenter;
                          hProtPosPt["ptyieldspp"]->fill(partPt, pt_weight);
                      }
                      hProtPosPt["pp1"]->fill(partPt);
                      hProtPosPt["pp2"]->fill(partPt);
                      hKpPt["Raa_c12_pp"]->fill(partPt);
                      break;
                  }
                  case -2212: // anti-proton
                  {
                      if(getBinCenter(*hProtNegPt["ptyieldspp"], partPt, binCenter))
                      {
                          pt_weight /= binCenter;
                          hProtNegPt["ptyieldspp"]->fill(partPt, pt_weight);
                      }
                      hProtNegPt["pp1"]->fill(partPt);
                      hProtNegPt["pp2"]->fill(partPt);
                      break;
                  }
              }
          }

          for (Particle p : neutralParticles)
          {
              double partPt = p.pT() / GeV;
              double pt_weight = 1. / (2. * M_PI);
              double binCenter = 0.;

              switch (p.pid()) {
                  case 310: // K0S
                  {
                      if(getBinCenter(*hKaon0SPt["ptyieldspp"], partPt, binCenter))
                      {
                          pt_weight /= binCenter;
                          hKaon0SPt["ptyieldspp"]->fill(partPt, pt_weight);
                      }
                      hKaon0SPt["pp"]->fill(partPt);
                      hKaon0SPt["Raa_c12_pp"]->fill(partPt);
                      break;
                  }
                  case 113: // rho0
                  {
                      if(getBinCenter(*hRho0Pt["ptyieldspp"], partPt, binCenter))
                      {
                          pt_weight /= binCenter;
                          hRho0Pt["ptyieldspp"]->fill(partPt, pt_weight);
                      }
                      hRho0Pt["Raa_c12_pp"]->fill(partPt);
                      break;
                  }
              }
          }

          return;
		}


		const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
		const double c = cent();

		if (collSys == AuAu200)

		{
			if ((c < 0.) || (c > 12.)) vetoEvent;
			sow["sow_AuAuc12"]->fill();
			for (Particle p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (2. * M_PI);
        double binCenter = 0.;

				switch (p.pid()) {
				case 211: // pi+
				{
					hPionPosPt["AuAuc12"]->fill(partPt);
					hPionPt["Raa_c12_AuAu"]->fill(partPt);
					break;
				}
				case -211: // pi-
				{
					hPionNegPt["AuAuc12"]->fill(partPt);
					hPionPt["Raa_c12_AuAu"]->fill(partPt);
					break;
				}
				case  321: // K+
				{
          if(getBinCenter(*hKpPosPt["ptyieldsAuAuc12"], partPt, binCenter))
          {
              pt_weight /= binCenter;
              hKpPosPt["ptyieldsAuAuc12"]->fill(partPt, pt_weight);
          }
					hKpPt["Raa_c12_AuAu"]->fill(partPt);
					break;
				}
				case  -321: // K-
				{
          if(getBinCenter(*hKpNegPt["ptyieldsAuAuc12"], partPt, binCenter))
          {
              pt_weight /= binCenter;
              hKpNegPt["ptyieldsAuAuc12"]->fill(partPt, pt_weight);
          }
					hKpPt["Raa_c12_AuAu"]->fill(partPt);
					break;
				}
				case 2212: // proton
				{
          if(getBinCenter(*hKpPosPt["ptyieldsAuAuc12"], partPt, binCenter))
          {
              pt_weight /= binCenter;
              hKpPosPt["ptyieldsAuAuc12"]->fill(partPt, pt_weight);
          }
					hProtPosPt["AuAuc12"]->fill(partPt);
					hKpPt["Raa_c12_AuAu"]->fill(partPt);
					break;
				}
				case -2212: // anti-proton
				{
          if(getBinCenter(*hKpNegPt["ptyieldsAuAuc12"], partPt, binCenter))
          {
              pt_weight /= binCenter;
              hKpNegPt["ptyieldsAuAuc12"]->fill(partPt, pt_weight);
          }
					hProtNegPt["AuAuc12"]->fill(partPt);
					hKpPt["Raa_c12_AuAu"]->fill(partPt);
					break;
				}
				}
			}

			for (Particle p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (2. * M_PI);
        double binCenter = 0.;

				switch (p.pid()) {
				case 310: // K0S
				{
          if(getBinCenter(*hKaon0SPt["ptyieldsAuAuc12"], partPt, binCenter))
          {
              pt_weight /= binCenter;
              hKaon0SPt["ptyieldsAuAuc12"]->fill(partPt, pt_weight);
          }
					hKaon0SPt["Raa_c12_AuAu"]->fill(partPt);
					break;
				}
				case 113: // rho0
				{
          if(getBinCenter(*hRho0Pt["ptyieldsAuAuc12"], partPt, binCenter))
          {
              pt_weight /= binCenter;
              hRho0Pt["ptyieldsAuAuc12"]->fill(partPt, pt_weight);
          }
					hRho0Pt["Raa_c12_AuAu"]->fill(partPt);
					break;
				}
				}
			}

			return;
		}


    }


    void finalize() {
		bool AuAu200_available = false;
		bool pp_available = false;

    if(sow["sow_pp"]->sumW() > 0) pp_available = true;
    if(sow["sow_AuAuc12"]->sumW() > 0) AuAu200_available = true;

		if (!(AuAu200_available && pp_available)) return;


		//Figure 1 Yield pp
		hPionPosPt["ptyieldspp"]->scaleW(1. / sow["sow_pp"]->sumW());
		hPionNegPt["ptyieldspp"]->scaleW(1. / sow["sow_pp"]->sumW());
		hKaonPosPt["ptyieldspp"]->scaleW(1. / sow["sow_pp"]->sumW());
		hKaonNegPt["ptyieldspp"]->scaleW(1. / sow["sow_pp"]->sumW());
		hProtPosPt["ptyieldspp"]->scaleW(1. / sow["sow_pp"]->sumW());
		hProtNegPt["ptyieldspp"]->scaleW(1. / sow["sow_pp"]->sumW());
		hKaon0SPt["ptyieldspp"]->scaleW(1. / sow["sow_pp"]->sumW());
		hRho0Pt["ptyieldspp"]->scaleW(1. / sow["sow_pp"]->sumW());

		//Figure 1 Yield AUAU
		hKpPosPt["ptyieldsAuAuc12"]->scaleW(1. / sow["sow_AuAuc12"]->sumW());
		hKpNegPt["ptyieldsAuAuc12"]->scaleW(1. / sow["sow_AuAuc12"]->sumW());
		hRho0Pt["ptyieldsAuAuc12"]->scaleW(1. / sow["sow_AuAuc12"]->sumW());
		hKaon0SPt["ptyieldsAuAuc12"]->scaleW(1. / sow["sow_AuAuc12"]->sumW());

		//Figure 2 Yield Ratio pp
		divide(hPionNegPt["pp1"], hPionPosPt["pp1"], RatioPion["pp"]);
		divide(hProtNegPt["pp1"], hProtPosPt["pp1"], RatioProt["pp"]);
		divide(hKaonNegPt["pp"], hKaonPosPt["pp"], RatioKaon["pp"]);
		divide(hProtPosPt["pp2"], hPionPosPt["pp2"], Ratioppipos["pp"]);
		divide(hProtNegPt["pp2"], hPionNegPt["pp2"], Ratioppineg["pp"]);
		divide(hKaonPt["pp"], hPionPt["pp1"], RatioKpi["pp"]);
		divide(hKaon0SPt["pp"], hPionPt["pp2"], RatioK0spi["pp"]);

		//Figure 2 Yield Ratio AUAU
		divide(hProtPosPt["AuAuc12"], hPionPosPt["AuAuc12"], Ratioppipos["AuAuc12"]);
		divide(hProtNegPt["AuAuc12"], hPionNegPt["AuAuc12"], Ratioppineg["AuAuc12"]);

		//Figure 3 RAA
		hPionPt["Raa_c12_AuAu"]->scaleW(1. / sow["sow_AuAuc12"]->sumW());
		hPionPt["Raa_c12_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hPionPt["Raa_c12_AuAu"], hPionPt["Raa_c12_pp"], hRaa["pi_c12_AuAu"]);
		//hRaa["pi_c12_AuAu"]->scaleY(1. / 960.2);

		hKpPt["Raa_c12_AuAu"]->scaleW(1. / sow["sow_AuAuc12"]->sumW());
		hKpPt["Raa_c12_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hKpPt["Raa_c12_AuAu"], hKpPt["Raa_c12_pp"], hRaa["Kp_c12_AuAu"]);
		//hRaa["Kp_c12_AuAu"]->scaleY(1. / 960.2);

		hKaon0SPt["Raa_c12_AuAu"]->scaleW(1. / sow["sow_AuAuc12"]->sumW());
		hKaon0SPt["Raa_c12_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hKaon0SPt["Raa_c12_AuAu"], hKaon0SPt["Raa_c12_pp"], hRaa["K0S_c12_AuAu"]);
		//hRaa["K0S_c12_AuAu"]->scaleY(1. / 960.2);

		hRho0Pt["Raa_c12_AuAu"]->scaleW(1. / sow["sow_AuAuc12"]->sumW());
		hRho0Pt["Raa_c12_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
		divide(hRho0Pt["Raa_c12_AuAu"], hRho0Pt["Raa_c12_pp"], hRaa["Rho0_c12_AuAu"]);
		//hRaa["Rho0_c13_AuAu"]->scaleY(1. / 960.2);

	    	//Figure 3 RAA Ratio_____________Need to implement


    }


	map<string, Histo1DPtr> hKaonNegPt;
	map<string, Histo1DPtr> hKaonPosPt;
	map<string, Histo1DPtr> hPionNegPt;
	map<string, Histo1DPtr> hPionPosPt;
	map<string, Histo1DPtr> hProtNegPt;
	map<string, Histo1DPtr> hProtPosPt;
	map<string, Histo1DPtr> hKaon0SPt;
	map<string, Histo1DPtr> hRho0Pt;
	map<string, Histo1DPtr> hKpNegPt;
	map<string, Histo1DPtr> hKpPosPt;
	map<string, Histo1DPtr> hKpPt;


	map<string, Histo1DPtr> hKaonPt;
	map<string, Histo1DPtr> hPionPt;
	map<string, Histo1DPtr> hProtPt;

	map<string, Scatter2DPtr> RatioPion;
	map<string, Scatter2DPtr> RatioProt;
	map<string, Scatter2DPtr> RatioKaon;
	map<string, Scatter2DPtr> Ratioppipos;
	map<string, Scatter2DPtr> Ratioppineg;
  map<string, Scatter2DPtr> RatioKpi;
  map<string, Scatter2DPtr> RatioK0spi;

	map<string, Scatter2DPtr> hRaa;


	map<string, CounterPtr> sow;
	string beamOpt;
	enum CollisionSystem { pp, AuAu200 };
	CollisionSystem collSys;

  };


  DECLARE_RIVET_PLUGIN(STAR_2012_I930463);

}
