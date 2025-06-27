// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "../Centralities/RHICCentrality.hh" //external header for Centrality calculation
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES

namespace Rivet {

	class STAR_2019_I1711377 : public Analysis {
	public:

		RIVET_DEFAULT_ANALYSIS_CTOR(STAR_2019_I1711377);


		void init() {
			std::initializer_list<int> pdgIds = { 421, -421 }; //D0 and D0bar

			const PrimaryParticles fs(pdgIds, Cuts::absrap < 1 && Cuts::abscharge == 0);
			declare(fs, "fs");

			beamOpt = getOption<string>("beam", "NONE");

			if (beamOpt == "PP") collSys = pp;
			else if (beamOpt == "AUAU200") collSys = AuAu200;


			if (!(collSys == pp)) declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");


			book(sow["sow_pp"], "sow_pp");

			for (int i = 0, N = AUAUCentralityBins.size(); i < N; ++i)
			{

				//yields (fig 22)_________________
				book(hD0D0barPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 1, 1, 1 + i);

				//yields (fig 33)_________________
				book(hD0Pt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 6, 1, 1 + i);
				book(hD0barPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 6, 1, 5 + i);

				book(sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])], "sow_AUAUc" + std::to_string(AUAUCentralityBins[i]));

				//ratios (fig 34)_________________
				string refname1 = mkAxisCode(7, 1, 1 + i);
				const Estimate1D& refdata1 = refData(refname1);
				book(r_D0bar["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname1 + "_D0bar", refdata1);
				book(r_D0["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname1 + "_D0", refdata1);
				book(RatioD0barD0["AuAuc" + std::to_string(AUAUCentralityBins[i])], refname1);

			}

			//yields (fig 22)_________________
			book(hD0D0barPt["ptyieldsAuAuc1040"], 1, 1, 6);
			book(hD0D0barPt["ptyieldsAuAuc4080"], 1, 1, 7);

			//RAA (fig 31)_________________
			string refname2 = mkAxisCode(5, 1, 1);
			const Estimate1D& refdata2 = refData(refname2);
			book(hD0Pt["Raa_c10_AuAu"], refname2 + "_AuAu", refdata2);
			book(hD0Pt["Raa_c10_pp"], refname2 + "_pp", refdata2);
			book(hRaa["D0_c10_AuAu"], refname2);

			string refname3 = mkAxisCode(5, 1, 2);
			const Estimate1D& refdata3 = refData(refname3);
			book(hD0Pt["Raa_c40_AuAu"], refname3 + "_AuAu", refdata3);
			book(hD0Pt["Raa_c40_pp"], refname3 + "_pp", refdata3);
			book(hRaa["D0_c40_AuAu"], refname3);

			book(sow["sow_AUAUc1040"], "sow_AUAUc1040");

			string refname4 = mkAxisCode(5, 1, 3);
			const Estimate1D& refdata4 = refData(refname4);
			book(hD0Pt["Raa_c80_AuAu"], refname4 + "_AuAu", refdata4);
			book(hD0Pt["Raa_c80_pp"], refname4 + "_pp", refdata4);
			book(hRaa["D0_c80_AuAu"], refname4);

			book(sow["sow_AUAUc4080"], "sow_AUAUc4080");

			//RCP ref 6080 (figs 29 & 35)_________________
			string refname5 = mkAxisCode(8, 1, 1);
			const Estimate1D& refdata5 = refData(refname5);
			book(hD0Pt["AuAuc0010a"], refname5 + "_D0", refdata5);
			book(hD0Pt["AuAuc6080a"], refname5 + "_D0", refdata5);
			book(hRcp["D0_c00106080_AuAu"], refname5);

			string refname6 = mkAxisCode(8, 1, 2);
			const Estimate1D& refdata6 = refData(refname6);
			book(hD0Pt["AuAuc1020a"], refname6 + "_D0", refdata6);
			book(hD0Pt["AuAuc6080b"], refname6 + "_D0", refdata6);
			book(hRcp["D0_c10206080_AuAu"], refname6);

			string refname7 = mkAxisCode(8, 1, 3);
			const Estimate1D& refdata7 = refData(refname7);
			book(hD0Pt["AuAuc2040a"], refname7 + "_D0", refdata7);
			book(hD0Pt["AuAuc6080c"], refname7 + "_D0", refdata7);
			book(hRcp["D0_c20406080_AuAu"], refname7);

			string refname8 = mkAxisCode(8, 1, 4);
			const Estimate1D& refdata8 = refData(refname8);
			book(hD0Pt["AuAuc4060a"], refname8 + "_D0", refdata8);
			book(hD0Pt["AuAuc6080d"], refname8 + "_D0", refdata8);
			book(hRcp["D0_c40606080_AuAu"], refname8);

			//RCP ref 4060 (figs 30 & 36)_________________
			string refname9 = mkAxisCode(9, 1, 1);
			const Estimate1D& refdata9 = refData(refname9);
			book(hD0Pt["AuAuc0010b"], refname9 + "_D0", refdata9);
			book(hD0Pt["AuAuc4060b"], refname9 + "_D0", refdata9);
			book(hRcp["D0_c00104060_AuAu"], refname9);

			string refname10 = mkAxisCode(9, 1, 2);
			const Estimate1D& refdata10 = refData(refname10);
			book(hD0Pt["AuAuc1020b"], refname10 + "_D0", refdata10);
			book(hD0Pt["AuAuc4060c"], refname10 + "_D0", refdata10);
			book(hRcp["D0_c10204060_AuAu"], refname10);

			string refname11 = mkAxisCode(9, 1, 3);
			const Estimate1D& refdata11 = refData(refname11);
			book(hD0Pt["AuAuc2040b"], refname11 + "_D0", refdata11);
			book(hD0Pt["AuAuc4060d"], refname11 + "_D0", refdata11);
			book(hRcp["D0_c20404060_AuAu"], refname11);

			book(hD0Pt["ptyieldspp"], 2, 1, 1);


			//NEED X-SECTIONS (3, 1, 1) (3, 1, 2) (4, 1, 1)


		}


		void analyze(const Event& event) {
			Particles neutralParticles = apply<PrimaryParticles>(event, "fs").particles();

			if (collSys == pp)
			{
				sow["sow_pp"]->fill();
				for (Particle p : neutralParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					if (p.pid() == 421) // D0
					{
						hD0Pt["ptyieldspp"]->fill(partPt, pt_weight);
						hD0Pt["Raa_c10_pp"]->fill(partPt);
						hD0Pt["Raa_c40_pp"]->fill(partPt);
						hD0Pt["Raa_c80_pp"]->fill(partPt);

					}
				}
				return;
			}

			const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
			const double c = cent();

			if (collSys == AuAu200)
			{
				if ((c < 0.) || (c > 80.)) vetoEvent;



				if ((c >= 0.) && (c < 10.))
				{
					sow["sow_AUAUc10"]->fill();
					for (const Particle& p : neutralParticles)
					{
						double partPt = p.pT() / GeV;
						double pt_weight = 1. / (partPt * 2. * M_PI);

						switch (p.pid()) {
						case 421: // D0
						{
							hD0D0barPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
							hD0D0barPt["ptyieldsAuAuc1040"]->fill(partPt, pt_weight);
							hD0Pt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
							r_D0["AuAuc10"]->fill(partPt);
							hD0Pt["Raa_c10_AuAu"]->fill(partPt);
							hD0Pt["AuAuc0010a"]->fill(partPt);
							hD0Pt["AuAuc0010b"]->fill(partPt);

							break;
						}
						case -421: // D0bar
						{
							hD0D0barPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
							hD0D0barPt["ptyieldsAuAuc1040"]->fill(partPt, pt_weight);
							hD0barPt["ptyieldsAuAuc10"]->fill(partPt, pt_weight);
							r_D0bar["AuAuc10"]->fill(partPt);

							break;
						}
						}
					}
				}

				else if ((c >= 10.) && (c < 20.))
				{
					sow["sow_AUAUc20"]->fill();
					sow["sow_AUAUc1040"]->fill();
					for (const Particle& p : neutralParticles)
					{
						double partPt = p.pT() / GeV;
						double pt_weight = 1. / (partPt * 2. * M_PI);

						switch (p.pid()) {
						case 421: // D0
						{
							hD0D0barPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
							hD0D0barPt["ptyieldsAuAuc1040"]->fill(partPt, pt_weight);
							hD0Pt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
							r_D0["AuAuc20"]->fill(partPt);
							hD0Pt["Raa_c40_AuAu"]->fill(partPt);
							hD0Pt["AuAuc1020a"]->fill(partPt);
							hD0Pt["AuAuc1020b"]->fill(partPt);
							

							break;
						}
						case -421: // D0bar
						{
							hD0D0barPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
							hD0D0barPt["ptyieldsAuAuc1040"]->fill(partPt, pt_weight);
							hD0barPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
							r_D0bar["AuAuc20"]->fill(partPt);

							break;
						}
						}
					}
				}

				else if ((c >= 20.) && (c < 40.))
				{
					sow["sow_AUAUc40"]->fill();
					sow["sow_AUAUc1040"]->fill();
					for (const Particle& p : neutralParticles)
					{
						double partPt = p.pT() / GeV;
						double pt_weight = 1. / (partPt * 2. * M_PI);

						switch (p.pid()) {
						case 421: // D0
						{
							hD0D0barPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
							hD0D0barPt["ptyieldsAuAuc1040"]->fill(partPt, pt_weight);
							hD0Pt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
							r_D0["AuAuc40"]->fill(partPt);
							hD0Pt["Raa_c40_AuAu"]->fill(partPt);
							hD0Pt["AuAuc2040a"]->fill(partPt);
							hD0Pt["AuAuc2040b"]->fill(partPt);

							break;
						}
						case -421: // D0bar
						{
							hD0D0barPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
							hD0D0barPt["ptyieldsAuAuc1040"]->fill(partPt, pt_weight);
							hD0barPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
							r_D0bar["AuAuc40"]->fill(partPt);

							break;
						}
						}
					}
				}
				else if ((c >= 40.) && (c < 60.))
				{
					sow["sow_AUAUc60"]->fill();
					sow["sow_AUAUc4080"]->fill();
					for (const Particle& p : neutralParticles)
					{
						double partPt = p.pT() / GeV;
						double pt_weight = 1. / (partPt * 2. * M_PI);

						switch (p.pid()) {
						case 421: // D0
						{
							hD0D0barPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
							hD0D0barPt["ptyieldsAuAuc4080"]->fill(partPt, pt_weight);
							hD0Pt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
							r_D0["AuAuc60"]->fill(partPt);
							hD0Pt["Raa_c80_AuAu"]->fill(partPt);
							hD0Pt["AuAuc4060a"]->fill(partPt);
							hD0Pt["AuAuc4060b"]->fill(partPt);
							hD0Pt["AuAuc4060c"]->fill(partPt);
							hD0Pt["AuAuc4060d"]->fill(partPt);

							break;
						}
						case -421: // D0bar
						{
							hD0D0barPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
							hD0D0barPt["ptyieldsAuAuc4080"]->fill(partPt, pt_weight);
							hD0barPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
							r_D0bar["AuAuc60"]->fill(partPt);

							break;
						}
						}
					}
				}

				else if ((c >= 60.) && (c < 80.))
				{
					sow["sow_AUAUc80"]->fill();
					sow["sow_AUAUc4080"]->fill();
					for (const Particle& p : neutralParticles)
					{
						double partPt = p.pT() / GeV;
						double pt_weight = 1. / (partPt * 2. * M_PI);

						switch (p.pid()) {
						case 421: // D0
						{
							hD0D0barPt["ptyieldsAuAuc80"]->fill(partPt, pt_weight);
							hD0D0barPt["ptyieldsAuAuc4080"]->fill(partPt, pt_weight);
							hD0Pt["ptyieldsAuAuc80"]->fill(partPt, pt_weight);
							r_D0["AuAuc80"]->fill(partPt);
							hD0Pt["Raa_c80_AuAu"]->fill(partPt);
							hD0Pt["AuAuc6080a"]->fill(partPt);
							hD0Pt["AuAuc6080b"]->fill(partPt);
							hD0Pt["AuAuc6080c"]->fill(partPt);
							hD0Pt["AuAuc6080d"]->fill(partPt);


							break;
						}
						case -421: // D0bar
						{
							hD0D0barPt["ptyieldsAuAuc80"]->fill(partPt, pt_weight);
							hD0D0barPt["ptyieldsAuAuc4080"]->fill(partPt, pt_weight);
							hD0barPt["ptyieldsAuAuc80"]->fill(partPt, pt_weight);
							r_D0bar["AuAuc80"]->fill(partPt);

							break;
						}
						}
					}
				}
				return;
			}
		}

		void finalize() {
			bool AuAu200_available = false;
			bool pp_available = false;

			for (auto element : hD0Pt)
			{
				string name = element.second->name();
				if (name.find("AuAu") != std::string::npos)
				{
					if (element.second->numEntries() > 0) AuAu200_available = true;
					else
					{
						AuAu200_available = false;
						break;
					}
				}

				else if (name.find("pp") != std::string::npos)
				{
					if (element.second->numEntries() > 0) pp_available = true;
					else
					{
						pp_available = false;
						break;
					}
				}
			}
			for (auto element : hD0barPt)
			{
				string name = element.second->name();
				if (name.find("AuAu") != std::string::npos)
				{
					if (element.second->numEntries() > 0) AuAu200_available = true;
					else
					{
						AuAu200_available = false;
						break;
					}
				}

				else if (name.find("pp") != std::string::npos)
				{
					if (element.second->numEntries() > 0) pp_available = true;
					else
					{
						pp_available = false;
						break;
					}
				}
			}


			if (!(AuAu200_available && pp_available)) return;

			for (int i = 0, N = AUAUCentralityBins.size(); i < N; ++i)
			{
				//yields _________________
				hD0D0barPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / (sow["sow_AUAUc" +  std::to_string(AUAUCentralityBins[i])])->sumW()*2);
				hD0Pt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
				hD0barPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());

				//Ratio of yields 
				divide(r_D0bar["AuAuc" + std::to_string(AUAUCentralityBins[i])], r_D0["AuAuc" + std::to_string(AUAUCentralityBins[i])], RatioD0barD0["AuAuc" + std::to_string(AUAUCentralityBins[i])]);
			}
			
			//yields _________________
			hD0D0barPt["ptyieldsAuAuc1040"]->scaleW(1. / (sow["sow_AUAUc1040"])->sumW()*2);
			hD0D0barPt["ptyieldsAuAuc4080"]->scaleW(1. / (sow["sow_AUAUc4080"])->sumW()*2);

			//RAA
			hD0Pt["Raa_c10_AuAu"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
			hD0Pt["Raa_c10_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
			divide(hD0Pt["Raa_c10_AuAu"], hD0Pt["Raa_c10_pp"], hRaa["D0_c10_AuAu"]);

			hD0Pt["Raa_c40_AuAu"]->scaleW(1. / sow["sow_AUAUc1040"]->sumW());
			hD0Pt["Raa_c40_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
			divide(hD0Pt["Raa_c40_AuAu"], hD0Pt["Raa_c40_pp"], hRaa["D0_c40_AuAu"]);

			hD0Pt["Raa_c80_AuAu"]->scaleW(1. / sow["sow_AUAUc4080"]->sumW());
			hD0Pt["Raa_c80_pp"]->scaleW(1. / sow["sow_pp"]->sumW());
			divide(hD0Pt["Raa_c80_AuAu"], hD0Pt["Raa_c80_pp"], hRaa["D0_c80_AuAu"]);

			//RCP ref 6080 _________________
			hD0Pt["AuAuc0010a"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
			hD0Pt["AuAuc6080a"]->scaleW(1. / sow["sow_AUAUc80"]->sumW());
			divide(hD0Pt["AuAuc0010a"], hD0Pt["AuAuc6080a"], hRcp["D0_c00106080_AuAu"]);

			hD0Pt["AuAuc1020a"]->scaleW(1. / sow["sow_AUAUc20"]->sumW());
			hD0Pt["AuAuc6080b"]->scaleW(1. / sow["sow_AUAUc80"]->sumW());
			divide(hD0Pt["AuAuc1020a"], hD0Pt["AuAuc6080b"], hRcp["D0_c10206080_AuAu"]);

			hD0Pt["AuAuc2040a"]->scaleW(1. / sow["sow_AUAUc40"]->sumW());
			hD0Pt["AuAuc6080c"]->scaleW(1. / sow["sow_AUAUc80"]->sumW());
			divide(hD0Pt["AuAuc2040a"], hD0Pt["AuAuc6080c"], hRcp["D0_c20406080_AuAu"]);

			hD0Pt["AuAuc4060a"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
			hD0Pt["AuAuc6080d"]->scaleW(1. / sow["sow_AUAUc80"]->sumW());
			divide(hD0Pt["AuAuc4060a"], hD0Pt["AuAuc6080d"], hRcp["D0_c40606080_AuAu"]);
			
			//RCP ref 4060 _________________
			hD0Pt["AuAuc0010b"]->scaleW(1. / sow["sow_AUAUc10"]->sumW());
			hD0Pt["AuAuc4060b"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
			divide(hD0Pt["AuAuc0010b"], hD0Pt["AuAuc4060b"], hRcp["D0_c00104060_AuAu"]);

			hD0Pt["AuAuc1020b"]->scaleW(1. / sow["sow_AUAUc20"]->sumW());
			hD0Pt["AuAuc4060c"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
			divide(hD0Pt["AuAuc1020b"], hD0Pt["AuAuc4060c"], hRcp["D0_c10204060_AuAu"]);

			hD0Pt["AuAuc2040b"]->scaleW(1. / sow["sow_AUAUc40"]->sumW());
			hD0Pt["AuAuc4060d"]->scaleW(1. / sow["sow_AUAUc60"]->sumW());
			divide(hD0Pt["AuAuc2040b"], hD0Pt["AuAuc4060d"], hRcp["D0_c20404060_AuAu"]);


			hD0Pt["ptyieldspp"]->scaleW(1. / sow["sow_pp"]->sumW());

			//scaling for RAA 
			hRaa["D0_c10_AuAu"]->scale(1. / 938.8);
			//hRaa["D0_c40_AuAu"]->scale(1. / 609.5); \\ 10-20: 579.9 and 20-40: 288.3
			//hRaa["D0_c80_AuAu"]->scale(1. / 300.8); \\ 40-60: 91.3 and 60-80: 21.3

		}

		map<string, Histo1DPtr> hD0D0barPt;
		map<string, Histo1DPtr> hD0Pt;
		map<string, Histo1DPtr> hD0barPt;

		map<string, Histo1DPtr> r_D0bar;
		map<string, Histo1DPtr> r_D0;
		map<string, Estimate1DPtr> RatioD0barD0;


		map<string, Estimate1DPtr> hRaa;
		map<string, Estimate1DPtr> hRcp;
		map<string, CounterPtr> sow;
		string beamOpt;
		enum CollisionSystem { pp, AuAu200 };
		CollisionSystem collSys;
		vector<int> AUAUCentralityBins{ 10, 20, 40, 60, 80 };


  };


  RIVET_DECLARE_PLUGIN(STAR_2019_I1711377);

}
