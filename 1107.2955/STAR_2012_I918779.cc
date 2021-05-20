// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
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


  /// @brief Add a short analysis description here
  class STAR_2012_I918779 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2012_I918779);

	void init() {
		std::initializer_list<int> pdgIds = { 310, 3122, -3122, 3312, -3312, 3334, -3334 };  // 310 K0S , 3122 LAMBDA, -3122 ANTILAMBDA, 3312 XI-, -3312 ANTIXI+, 3334 OMEGA-,	-3334 ANTIOMEGA+	 



		//charged particles
		const PrimaryParticles cp(pdgIds, Cuts::absrap < 0.5 && Cuts::abscharge > 0);
		declare(cp, "cp");

		//neutral particles 
		const PrimaryParticles np(pdgIds, Cuts::absrap < 0.5 && Cuts::abscharge == 0);
		declare(np, "np");

		beamOpt = getOption<string>("beam", "NONE");

		if (beamOpt == "PP") collSys = pp;
		else if (beamOpt == "CUCU200") collSys = CuCu200;
		else if (beamOpt == "AUAU200") collSys = AuAu200;


		if (!(collSys == pp)) declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");


		book(sow["sow_pp"], "sow_pp");

		for (int i = 0, N = CUCUCentralityBins.size(); i < N; ++i)
		{

			//Figure 1 Yields _________________
			book(hKaon0SPt["ptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])], 1, 1, 1 + i);
			book(hLambdaPt["LambdaptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])], 2, 1, 1 + i);
			book(hLambdaPt["LambdabarptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])], 2, 1, 6 + i);
			book(hXiPt["XiptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])], 3, 1, 1 + i);
			book(hXiPt["XibarptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])], 3, 1, 6 + i);
			book(hOmegaPt["ptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])], 4, 1, 1 + i);

			book(hKaon0SPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 5, 1, 1 + i);
			book(hLambdaPt["LambdaptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 6, 1, 1 + i);
			book(hLambdaPt["LambdabarptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])], 6, 1, 6 + i);

			book(sow["sow_CUCUc" + std::to_string(CUCUCentralityBins[i])], "sow_CUCUc" + std::to_string(CUCUCentralityBins[i]));
			book(sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])], "sow_AUAUc" + std::to_string(AUAUCentralityBins[i]));
		}
			

		//Figure 2 spectra/Npart _______________

		book(hKaon0SPt["spectraCuCuc10"], 7, 1, 1);
		book(hKaon0SPt["spectraAuAuc40"], 7, 1, 2);

		book(hLambdaPt["spectraCuCuc10"], 7, 1, 3);
		book(hLambdaPt["spectraAuAuc40"], 7, 1, 4);

		book(hXiPt["spectraCuCuc10"], 8, 1, 1);
		book(hXiPt["spectraAuAuc40"], 9, 1, 1);
		
		book(hOmegaPt["spectraCuCuc10"], 10, 1, 1);
		book(hOmegaPt["spectraAuAuc40"], 11, 1, 1);


		//Figure 3 Strangeness Enhancement vs Centrality 12/1/1-5 and 13/1/1-5



	}



	void analyze(const Event& event) {
		Particles chargedParticles = applyProjection<PrimaryParticles>(event, "cp").particles();
		Particles neutralParticles = applyProjection<PrimaryParticles>(event, "np").particles();


		//pp for figure 3
		if (collSys == pp)
		{
			sow["sow_pp"]->fill();
			for (Particle p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {

				case 3312: //Xi
				{
					
					break;
				}
				case  -3312: //Xibar
				{
					
					break;
				}
				case 3334: // Omega
				{
					
					break;
				}
				case -3334: // Omegabar
				{
					
					break;
				}
				}
			}

			for (Particle p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {

				case 310: // K0S
				{
					
					break;
				}
				case 3122: // lambda
				{
					
					break;
				}
				case -3122: // lambdabar
				{
					
					break;
				}
				}
			}

			return;
		}

		const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
		const double c = cent();

		if (collSys == CuCu200)

		{
			if ((c < 0.) || (c > 60.)) vetoEvent;


			if ((c >= 0.) && (c < 10.))
			{
				sow["sow_CuCuc10"]->fill();
				for (Particle p : chargedParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 3312: //Xi
					{
						hXiPt["XiptyieldsCuCuc10"]->fill(partPt, pt_weight);
						hXiPt["spectraCuCuc10"]->fill(partPt, pt_weight);
						break;
					}
					case  -3312: //Xibar
					{
						hXiPt["XibarptyieldsCuCuc10"]->fill(partPt, pt_weight);
						hXiPt["spectraCuCuc10"]->fill(partPt, pt_weight);
						break;
					}
					case 3334: // Omega
					{
						hOmegaPt["ptyieldsCuCuc10"]->fill(partPt, pt_weight);
						hOmegaPt["spectraCuCuc10"]->fill(partPt, pt_weight);
						break;
					}
					case -3334: // Omegabar
					{
						hOmegaPt["ptyieldsCuCuc10"]->fill(partPt, pt_weight);
						hOmegaPt["spectraCuCuc10"]->fill(partPt, pt_weight);
						break;
					}
					}
				}

				for (Particle p : neutralParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 310: // K0S
					{
						hKaon0SPt["ptyieldsCuCuc10"]->fill(partPt, pt_weight);
						hKaon0SPt["spectraCuCuc10"]->fill(partPt, pt_weight);
						break;
					}
					case 3122: // lambda
					{
						hLambdaPt["LambdaptyieldsCuCuc10"]->fill(partPt, pt_weight);
						hLambdaPt["spectraCuCuc10"]->fill(partPt, pt_weight);
						break;
					}
					case -3122: // lambdabar
					{
						hLambdaPt["LambdabarptyieldsCuCuc10"]->fill(partPt, pt_weight);
						hLambdaPt["spectraCuCuc10"]->fill(partPt, pt_weight);
						break;
					}
					}
				}
			}
			else if ((c >= 10.) && (c < 20.))
			{
				sow["sow_CuCuc20"]->fill();


				for (Particle p : chargedParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 3312: //Xi
					{
						hXiPt["XiptyieldsCuCuc20"]->fill(partPt, pt_weight);
						break;
					}
					case  -3312: //Xibar
					{
						hXiPt["XibarptyieldsCuCuc20"]->fill(partPt, pt_weight);
						break;
					}
					case 3334: // Omega
					{
						hOmegaPt["ptyieldsCuCuc20"]->fill(partPt, pt_weight);
						break;
					}
					case -3334: // Omegabar
					{
						hOmegaPt["ptyieldsCuCuc20"]->fill(partPt, pt_weight);
						break;
					}
					}
				}

				for (Particle p : neutralParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 310: // K0S
					{
						hKaon0SPt["ptyieldsCuCuc20"]->fill(partPt, pt_weight);
						break;
					}
					case 3122: // lambda
					{
						hLambdaPt["LambdaptyieldsCuCuc20"]->fill(partPt, pt_weight);
						break;
					}
					case -3122: // lambdabar
					{
						hLambdaPt["LambdabarptyieldsCuCuc20"]->fill(partPt, pt_weight);
						break;
					}
					}
				}
			}

			else if ((c >= 20.) && (c < 30.))
			{
			sow["sow_CuCuc30"]->fill();


			for (Particle p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {

				case 3312: //Xi
				{
					hXiPt["XiptyieldsCuCuc30"]->fill(partPt, pt_weight);
					break;
				}
				case  -3312: //Xibar
				{
					hXiPt["XibarptyieldsCuCuc30"]->fill(partPt, pt_weight);
					break;
				}
				case 3334: // Omega
				{
					hOmegaPt["ptyieldsCuCuc30"]->fill(partPt, pt_weight);
					break;
				}
				case -3334: // Omegabar
				{
					hOmegaPt["ptyieldsCuCuc30"]->fill(partPt, pt_weight);
					break;
				}
				}
			}

			for (Particle p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {

				case 310: // K0S
				{
					hKaon0SPt["ptyieldsCuCuc30"]->fill(partPt, pt_weight);
					break;
				}
				case 3122: // lambda
				{
					hLambdaPt["LambdaptyieldsCuCuc30"]->fill(partPt, pt_weight);
					break;
				}
				case -3122: // lambdabar
				{
					hLambdaPt["LambdabarptyieldsCuCuc30"]->fill(partPt, pt_weight);
					break;
				}
				}
			}
			}

			else if ((c >= 30.) && (c < 40.))
			{
			sow["sow_CuCuc40"]->fill();


			for (Particle p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {

				case 3312: //Xi
				{
					hXiPt["XiptyieldsCuCuc40"]->fill(partPt, pt_weight);
					break;
				}
				case  -3312: //Xibar
				{
					hXiPt["XibarptyieldsCuCuc40"]->fill(partPt, pt_weight);
					break;
				}
				case 3334: // Omega
				{
					hOmegaPt["ptyieldsCuCuc40"]->fill(partPt, pt_weight);
					break;
				}
				case -3334: // Omegabar
				{
					hOmegaPt["ptyieldsCuCuc40"]->fill(partPt, pt_weight);
					break;
				}
				}
			}

			for (Particle p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {

				case 310: // K0S
				{
					hKaon0SPt["ptyieldsCuCuc40"]->fill(partPt, pt_weight);
					break;
				}
				case 3122: // lambda
				{
					hLambdaPt["LambdaptyieldsCuCuc40"]->fill(partPt, pt_weight);
					break;
				}
				case -3122: // lambdabar
				{
					hLambdaPt["LambdabarptyieldsCuCuc40"]->fill(partPt, pt_weight);
					break;
				}
				}
			}
			}

			else if ((c >= 40.) && (c < 60.))
			{
			sow["sow_CuCuc60"]->fill();


			for (Particle p : chargedParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {

				case 3312: //Xi
				{
					hXiPt["XiptyieldsCuCuc60"]->fill(partPt, pt_weight);
					break;
				}
				case  -3312: //Xibar
				{
					hXiPt["XibarptyieldsCuCuc60"]->fill(partPt, pt_weight);
					break;
				}
				case 3334: // Omega
				{
					hOmegaPt["ptyieldsCuCuc60"]->fill(partPt, pt_weight);
					break;
				}
				case -3334: // Omegabar
				{
					hOmegaPt["ptyieldsCuCuc60"]->fill(partPt, pt_weight);
					break;
				}
				}
			}

			for (Particle p : neutralParticles)
			{
				double partPt = p.pT() / GeV;
				double pt_weight = 1. / (partPt * 2. * M_PI);

				switch (p.pid()) {

				case 310: // K0S
				{
					hKaon0SPt["ptyieldsCuCuc60"]->fill(partPt, pt_weight);
					break;
				}
				case 3122: // lambda
				{
					hLambdaPt["LambdaptyieldsCuCuc60"]->fill(partPt, pt_weight);
					break;
				}
				case -3122: // lambdabar
				{
					hLambdaPt["LambdabarptyieldsCuCuc60"]->fill(partPt, pt_weight);
					break;
				}
				}
			}
			}
			return;
		}

		if (collSys == AuAu200)

		{
			if ((c < 0.) || (c > 80.)) vetoEvent;


			if ((c >= 0.) && (c < 5.))
			{
				sow["sow_AuAuc5"]->fill();
				for (Particle p : chargedParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 3312: //Xi
					{
						
						break;
					}
					case  -3312: //Xibar
					{
						
						break;
					}
					case 3334: // Omega
					{
						
						break;
					}
					case -3334: // Omegabar
					{
						
						break;
					}
					}
				}

				for (Particle p : neutralParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 310: // K0S
					{
						hKaon0SPt["ptyieldsAuAuc5"]->fill(partPt, pt_weight);
						break;
					}
					case 3122: // lambda
					{
						hLambdaPt["LambdaptyieldsAuAuc5"]->fill(partPt, pt_weight);
						break;
					}
					case -3122: // lambdabar
					{
						hLambdaPt["LambdabarptyieldsAuAuc5"]->fill(partPt, pt_weight);
						break;
					}
					}
				}
			}
			else if ((c >= 10.) && (c < 20.))
			{
				sow["sow_AuAuc20"]->fill();


				for (Particle p : chargedParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 3312: //Xi
					{
						
						break;
					}
					case  -3312: //Xibar
					{
						
						break;
					}
					case 3334: // Omega
					{
					
						break;
					}
					case -3334: // Omegabar
					{
					
						break;
					}
					}
				}

				for (Particle p : neutralParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 310: // K0S
					{
						hKaon0SPt["ptyieldsAuAuc20"]->fill(partPt, pt_weight);
						break;
					}
					case 3122: // lambda
					{
						hLambdaPt["LambdaptyieldsAuAuc20"]->fill(partPt, pt_weight);
						break;
					}
					case -3122: // lambdabar
					{
						hLambdaPt["LambdabarptyieldsAuAuc20"]->fill(partPt, pt_weight);
						break;
					}
					}
				}
			}

			else if ((c >= 20.) && (c < 40.))
			{
				sow["sow_AuAuc40"]->fill();


				for (Particle p : chargedParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 3312: //Xi
					{
						hXiPt["spectraAuAuc40"]->fill(partPt, pt_weight);
						break;
					}
					case  -3312: //Xibar
					{
						hXiPt["spectraAuAuc40"]->fill(partPt, pt_weight);
						break;
					}
					case 3334: // Omega
					{
						hOmegaPt["spectraAuAuc40"]->fill(partPt, pt_weight);
						break;
					}
					case -3334: // Omegabar
					{
						hOmegaPt["spectraAuAuc40"]->fill(partPt, pt_weight);
						break;
					}
					}
				}

				for (Particle p : neutralParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 310: // K0S
					{
						hKaon0SPt["ptyieldsAuAuc40"]->fill(partPt, pt_weight);
						hKaon0SPt["spectraAuAuc40"]->fill(partPt, pt_weight);
						break;
					}
					case 3122: // lambda
					{
						hLambdaPt["LambdaptyieldsAuAuc40"]->fill(partPt, pt_weight);
						hLambdaPt["spectraAuAuc40"]->fill(partPt, pt_weight);
						break;
					}
					case -3122: // lambdabar
					{
						hLambdaPt["LambdabarptyieldsAuAuc40"]->fill(partPt, pt_weight);
						hLambdaPt["spectraAuAuc40"]->fill(partPt, pt_weight);
						break;
					}
					}
				}
			}

			else if ((c >= 40.) && (c < 60.))
			{
				sow["sow_AuAuc60"]->fill();


				for (Particle p : chargedParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 3312: //Xi
					{
						
						break;
					}
					case  -3312: //Xibar
					{
						
						break;
					}
					case 3334: // Omega
					{
						
						break;
					}
					case -3334: // Omegabar
					{
						
						break;
					}
					}
				}

				for (Particle p : neutralParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 310: // K0S
					{
						hKaon0SPt["ptyieldsAuAuc60"]->fill(partPt, pt_weight);
						break;
					}
					case 3122: // lambda
					{
						hLambdaPt["LambdaptyieldsAuAuc60"]->fill(partPt, pt_weight);
						break;
					}
					case -3122: // lambdabar
					{
						hLambdaPt["LambdabarptyieldsAuAuc60"]->fill(partPt, pt_weight);
						break;
					}
					}
				}
			}

			else if ((c >= 60.) && (c < 80.))
			{
				sow["sow_AuAuc80"]->fill();


				for (Particle p : chargedParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 3312: //Xi
					{
						
						break;
					}
					case  -3312: //Xibar
					{
						
						break;
					}
					case 3334: // Omega
					{
						
						break;
					}
					case -3334: // Omegabar
					{
						
						break;
					}
					}
				}

				for (Particle p : neutralParticles)
				{
					double partPt = p.pT() / GeV;
					double pt_weight = 1. / (partPt * 2. * M_PI);

					switch (p.pid()) {

					case 310: // K0S
					{
						hKaon0SPt["ptyieldsAuAuc80"]->fill(partPt, pt_weight);
						break;
					}
					case 3122: // lambda
					{
						hLambdaPt["LambdaptyieldsAuAuc80"]->fill(partPt, pt_weight);
						break;
					}
					case -3122: // lambdabar
					{
						hLambdaPt["LambdabarptyieldsAuAuc80"]->fill(partPt, pt_weight);
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
		bool CuCu200_available = false;
		bool pp_available = false;

		for (auto element : hKaon0SPt)
		{
			string name = element.second->name();
			if (name.find("CuCu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) CuCu200_available = true;
				else
				{
					CuCu200_available = false;
					break;
				}
			}
			else if (name.find("AuAu") != std::string::npos)
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

		for (auto element : hLambdaPt)
		{
			string name = element.second->name();
			if (name.find("CuCu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) CuCu200_available = true;
				else
				{
					CuCu200_available = false;
					break;
				}
			}
			else if (name.find("AuAu") != std::string::npos)
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

		for (auto element : hXiPt)
		{
			string name = element.second->name();
			if (name.find("CuCu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) CuCu200_available = true;
				else
				{
					CuCu200_available = false;
					break;
				}
			}
			else if (name.find("AuAu") != std::string::npos)
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

		for (auto element : hOmegaPt)
		{
			string name = element.second->name();
			if (name.find("CuCu") != std::string::npos)
			{
				if (element.second->numEntries() > 0) CuCu200_available = true;
				else
				{
					CuCu200_available = false;
					break;
				}
			}
			else if (name.find("AuAu") != std::string::npos)
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
		

		if (!(CuCu200_available && AuAu200_available && pp_available)) return;


		//Figure 1 Yields _________________
		for (int i = 0, N = CUCUCentralityBins.size(); i < N; ++i)
		{
			hKaon0SPt["ptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])]->scaleW(1. / sow["sow_CUCUc" + std::to_string(CUCUCentralityBins[i])]->sumW());
			hLambdaPt["LambdaptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])]->scaleW(1. / sow["sow_CUCUc" + std::to_string(CUCUCentralityBins[i])]->sumW());
			hLambdaPt["LambdabarptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])]->scaleW(1. / sow["sow_CUCUc" + std::to_string(CUCUCentralityBins[i])]->sumW());
			hXiPt["XiptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])]->scaleW(1. / sow["sow_CUCUc" + std::to_string(CUCUCentralityBins[i])]->sumW());
			hXiPt["XibarptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])]->scaleW(1. / sow["sow_CUCUc" + std::to_string(CUCUCentralityBins[i])]->sumW());
			hOmegaPt["ptyieldsCuCuc" + std::to_string(CUCUCentralityBins[i])]->scaleW(1. / sow["sow_CUCUc" + std::to_string(CUCUCentralityBins[i])]->sumW());

			hKaon0SPt["ptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
			hLambdaPt["LambdaptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
			hLambdaPt["LambdabarptyieldsAuAuc" + std::to_string(AUAUCentralityBins[i])]->scaleW(1. / sow["sow_AUAUc" + std::to_string(AUAUCentralityBins[i])]->sumW());
		}
		
		//Figure 2 spectra/Npart _______________
		hKaon0SPt["spectraCuCuc10"]->scaleW(1. / 99);
		hKaon0SPt["spectraAuAuc40"]->scaleW(1. / 141);

		hLambdaPt["spectraCuCuc10"]->scaleW(1. / 99);
		hLambdaPt["spectraAuAuc40"]->scaleW(1. / 141);

		hXiPt["spectraCuCuc10"]->scaleW(1. / 99);
		hXiPt["spectraAuAuc40"]->scaleW(1. / 141);

		hOmegaPt["spectraCuCuc10"]->scaleW(1. / 99);
		hOmegaPt["spectraAuAuc40"]->scaleW(1. / 141);





	}


	map<string, Histo1DPtr> hKaon0SPt;
	map<string, Histo1DPtr> hLambdaPt;
	map<string, Histo1DPtr> hXiPt;
	map<string, Histo1DPtr> hOmegaPt;

	map<string, CounterPtr> sow;
	string beamOpt;
	enum CollisionSystem { pp, CuCu200, AuAu200 };
	CollisionSystem collSys;
	vector<int> CUCUCentralityBins{ 10, 20, 30, 40, 60 };
	vector<int> AUAUCentralityBins{ 5, 20, 40, 60, 80 };

  };


  DECLARE_RIVET_PLUGIN(STAR_2012_I918779);

}
