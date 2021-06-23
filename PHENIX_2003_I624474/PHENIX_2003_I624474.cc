// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include <math.h>
#define _USE_MATH_DEFINES
namespace Rivet {


	/// @brief Add a short analysis description here
	class PHENIX_2003_I624474 : public Analysis {
	public:

	/// Constructor
	DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2003_I624474);

	/// Book histograms and initialise projections before the run
	void init() {
	
	// Particles: pi+, pi-, K+, K-, p, pbar, pi0
	// Initialise and register projections
      
	// The basic final-state projection:
	// all final-state particles within
	// the given eta acceptance
	std::intializer_list<int> pdgIds ={211, 321, 2212, -211, -321, -2212};

	const PrimaryParticles cp(pdgIds, Cuts::abseta < .35 && Cuts::pT > .5 && Cuts::abscharge > 0);	//fix pT//
	declare(cp, "cp");

	const UnstableParticles np(Cuts::abseta < .35 && Cuts::abspid == 111);
	declare(np, "np");

	declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

	beamOpt = getOption<string>("beam", "NONE");
	if (beamOpt == "AUAU200") collSys = AUAU200;

	// Book histograms
	// specify custom binning
      

	//___PiPlus yields___

	//Invariant yield of piplus minimum bias
	book(hAUAU_Yields["PiplusMin"],1,1,1);

	//Invariant yield of piplus 0-5%
	book(hAUAU_Yields["Piplus5"],1,1,2);

	//Invariant yield of piplus 5-10%
	book(hAUAU_Yields["Piplus10"],1,1,3);

	//Invariant yield of piplus 10-15% 
	book(hAUAU_Yields["Piplus15"],1,1,4);

	//Invariant yield of piplus 15-20%
	book(hAUAU_Yields["Piplus20"],1,1,5);

	//Invariant yield of piplus 20-30%
	book(hAUAU_Yields["Piplus30"],1,1,6);

	//Invariant yield of piplus 30-40%
	book(hAUAU_Yields["Piplus40"],1,1,7);

	//Invariant yield of piplus 40-50%
	book(hAUAU_Yields["Piplus50"],1,1,8);

	//Invariant yield of piplus 50-60%
	book(hAUAU_Yields["Piplus60"],1,1,9);

	//Invariant yield of piplus 60-70%
	book(hAUAU_Yields["Piplus70"],1,1,10);

	//Invariant yield of piplus 70-80%
	book(hAUAU_Yields["Piplus80"],1,1,11);

	//Invariant yield of piplus 80-92%
	book(hAUAU_Yields["Piplus92"],1,1,12);

	//Invariant yield of piplus 60-92%
	book(hAUAU_Yields["Piplus60_92"],1,1,13);


	//___PiMinus yields___

	//Invariant yield of piminus minimum bias
	book(hAUAU_Yields["PiminusMin"],2,1,1);

	//Invariant yield of piminus 0-5%
	book(hAUAU_Yields["Piminus5"],2,1,2);

	//Invariant yield of piminus 5-10%
	book(hAUAU_Yields["Piminus10"],2,1,3);

	//Invariant yield of piminus 10-15% 
	book(hAUAU_Yields["Piminus15"],2,1,4);

	//Invariant yield of piminus 15-20%
	book(hAUAU_Yields["Piminus20"],2,1,5);

	//Invariant yield of piminus 20-30%
	book(hAUAU_Yields["Piminus30"],2,1,6);

	//Invariant yield of piminus 30-40%
	book(hAUAU_Yields["Piminus40"],2,1,7);

	//Invariant yield of piminus 40-50%
	book(hAUAU_Yields["Piminus50"],2,1,8);

	//Invariant yield of piminus 50-60%
	book(hAUAU_Yields["Piminus60"],2,1,9);

	//Invariant yield of piminus 60-70%
	book(hAUAU_Yields["Piminus70"],2,1,10);

	//Invariant yield of piminus 70-80%
	book(hAUAU_Yields["Piminus80"],2,1,11);

	//Invariant yield of piminus 80-92%
	book(hAUAU_Yields["Piminus92"],2,1,12);

	//Invariant yield of piminus 60-92%
	book(hAUAU_Yields["Piminus60_92"],2,1,13);


	//___KPlus yields___

	//Invariant yield of Kplus minimum bias
	book(hAUAU_Yields["KplusMin"],3,1,1);

	//Invariant yield of Kplus 0-5%
	book(hAUAU_Yields["Kplus5"],3,1,2);

	//Invariant yield of Kplus 5-10%
	book(hAUAU_Yields["Kplus10"],3,1,3);

	//Invariant yield of Kplus 10-15% 
	book(hAUAU_Yields["Kplus15"],3,1,4);

	//Invariant yield of Kplus 15-20%
	book(hAUAU_Yields["Kplus20"],3,1,5);

	//Invariant yield of Kplus 20-30%
	book(hAUAU_Yields["Kplus30"],3,1,6);

	//Invariant yield of Kplus 30-40%
	book(hAUAU_Yields["Kplus40"],3,1,7);

	//Invariant yield of Kplus 40-50%
	book(hAUAU_Yields["Kplus50"],3,1,8);

	//Invariant yield of Kplus 50-60%
	book(hAUAU_Yields["Kplus60"],3,1,9);

	//Invariant yield of Kplus 60-70%
	book(hAUAU_Yields["Kplus70"],3,1,10);

	//Invariant yield of Kplus 70-80%
	book(hAUAU_Yields["Kplus80"],3,1,11);

	//Invariant yield of Kplus 80-92%
	book(hAUAU_Yields["Kplus92"],3,1,12);

	//Invariant yield of Kplus 60-92%
	book(hAUAU_Yields["Kplus60_92"],3,1,13);


	//___KMinus yields___

	//Invariant yield of Kminus minimum bias
	book(hAUAU_Yields["KminusMin"],4,1,1);

	//Invariant yield of Kminus 0-5%
	book(hAUAU_Yields["Kminus5"],4,1,2);

	//Invariant yield of Kminus 5-10%
	book(hAUAU_Yields["Kminus10"],4,1,3);

	//Invariant yield of Kminus 10-15% 
	book(hAUAU_Yields["Kminus15"],4,1,4);

	//Invariant yield of Kminus 15-20%
	book(hAUAU_Yields["Kminus20"],4,1,5);

	//Invariant yield of Kminus 20-30%
	book(hAUAU_Yields["Kminus30"],4,1,6);

	//Invariant yield of Kminus 30-40%
	book(hAUAU_Yields["Kminus40"],4,1,7);

	//Invariant yield of Kminus 40-50%
	book(hAUAU_Yields["Kminus50"],4,1,8);

	//Invariant yield of Kminus 50-60%
	book(hAUAU_Yields["Kminus60"],4,1,9);

	//Invariant yield of Kminus 60-70%
	book(hAUAU_Yields["Kminus70"],4,1,10);

	//Invariant yield of Kminus 70-80%
	book(hAUAU_Yields["Kminus80"],4,1,11);

	//Invariant yield of Kminus 80-92%
	book(hAUAU_Yields["Kminus92"],4,1,12);

	//Invariant yield of Kminus 60-92%
	book(hAUAU_Yields["Kminus60_92"],4,1,13);


	//___Proton yields___

	//Invariant yield of protons minimum bias
	book(hAUAU_Yields["ProtonsMin"],5,1,1);

	//Invariant yield of protons 0-5%
	book(hAUAU_Yields["Protons5"],5,1,2);

	//Invariant yield of protons 5-10%
	book(hAUAU_Yields["Protons10"],5,1,3);

	//Invariant yield of protons 10-15% 
	book(hAUAU_Yields["Protons15"],5,1,4);

	//Invariant yield of protons 15-20%
	book(hAUAU_Yields["Protons20"],5,1,5);

	//Invariant yield of protons 20-30%
	book(hAUAU_Yields["Protons30"],5,1,6);

	//Invariant yield of protons 30-40%
	book(hAUAU_Yields["Protons40"],5,1,7);

	//Invariant yield of protons 40-50%
	book(hAUAU_Yields["Protons50"],5,1,8);

	//Invariant yield of protons 50-60%
	book(hAUAU_Yields["Protons60"],5,1,9);

	//Invariant yield of protons 60-70%
	book(hAUAU_Yields["Protons70"],5,1,10);

	//Invariant yield of protons 70-80%
	book(hAUAU_Yields["Protons80"],5,1,11);

	//Invariant yield of protons 80-92%
	book(hAUAU_Yields["Protons92"],5,1,12);

	//Invariant yield of protons 60-92%
	book(hAUAU_Yields["Protons60_92"],5,1,13);


	//___Pbar yields___

	//Invariant yield of pbar minimum bias
	book(hAUAU_Yields["PbarMin"],6,1,1);

	//Invariant yield of pbar 0-5%
	book(hAUAU_Yields["Pbar5"],6,1,2);

	//Invariant yield of pbar 5-10%
	book(hAUAU_Yields["Pbar10"],6,1,3);

	//Invariant yield of pbar 10-15% 
	book(hAUAU_Yields["Pbar15"],6,1,4);

	//Invariant yield of pbar 15-20%
	book(hAUAU_Yields["Pbar20"],6,1,5);

	//Invariant yield of pbar 20-30%
	book(hAUAU_Yields["Pbar30"],6,1,6);

	//Invariant yield of pbar 30-40%
	book(hAUAU_Yields["Pbar40"],6,1,7);

	//Invariant yield of pbar 40-50%
	book(hAUAU_Yields["Pbar50"],6,1,8);

	//Invariant yield of pbar 50-60%
	book(hAUAU_Yields["Pbar60"],6,1,9);

	//Invariant yield of pbar 60-70%
	book(hAUAU_Yields["Pbar70"],6,1,10);

	//Invariant yield of pbar 70-80%
	book(hAUAU_Yields["Pbar80"],6,1,11);

	//Invariant yield of pbar 80-92%
	book(hAUAU_Yields["Pbar92"],7,1,1);

	//Invariant yield of pbar 60-92%
	book(hAUAU_Yields["Pbar60_92"],6,1,12);



	//____Counters____

	book(sow["sow_AUAUmin"],"sow_AUAUmin");
	book(sow["sow_AUAU5"],"sow_AUAU5");
	book(sow["sow_AUAU10"],"sow_AUAU10");
	book(sow["sow_AUAU15"],"sow_AUAU15");
	book(sow["sow_AUAU20"],"sow_AUAU20");
	book(sow["sow_AUAU30"],"sow_AUAU30");
	book(sow["sow_AUAU40"],"sow_AUAU40");
	book(sow["sow_AUAU50"],"sow_AUAU50");
	book(sow["sow_AUAU60"],"sow_AUAU60");
	book(sow["sow_AUAU70"],"sow_AUAU70");
	book(sow["sow_AUAU80"],"sow_AUAU80");
	book(sow["sow_AUAU92"],"sow_AUAU92");
	book(sow["sow_AUAU60_92"],"sow_AUAU60_92");
	book(sow["sow_AUAU0_10"],"sow_AUAU0_10");
	book(sow["sow_AUAUall"],"sow_AUAUall");

	//____Rcp____
	
	string refnameRcpPi = mkAxiscode(20,1,1);
	book(hRcp["Pions"], refnameRcpPi);

	string refnameRcpK = mkAxiscode(21,1,1);
	book(hRcp["Kaons"], refnameRcpK);

	string refnameRcpP = mkAxiscode(22,1,1);
	book(hRcp["Pbar+P"], refnameRcpP);

	string refnameRcpPi0 = mkAxiscode(23,1,1);
	book(hRcp["Pi0"], refnameRcpPi0);
	
	}


    /// Perform the per-event analysis
	void analyze(const Event& event) {

		const PrimaryParticles cp = apply<PrimaryParticles>(event, "cp");
		const UnstableParticles np = apply<UnstableParticles>(event, "np");
		const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
		const double c = cent();
		const ParticlePair& beam = beams();
		const Particles chargedParticles = cp.particles();

		if ((c < 0.) || (c > 92.2)) vetoEvent;

		sow["sow_AUAUall"]->fill();
			
		if ((c >= 0.) && (c < 5.)) {
				
			sow["sow_AUAU5"]->fill();
			sow["sow_AUAU0_10"]->fill();

			for (const Particle& p : chargedParticles) {
		
				double partPt = p.pT()/GeV;
				double pt_weight = 1./(partPt*2.*M_PI);
					
				switch (p.pid() {
		
					case 211:	//pi+
						
						hAUAU_Yields["Piplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piplus5"]->fill(partPt);
			
					case -211:	//pi-
						
						hAUAU_Yields["Piminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Piminus5"]->fill(partPt);
					
					case 321:	//K+
						
						hAUAU_Yields["Kplusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kplus5"]->fill(partPt);			

					case -321:	//K-
						
						hAUAU_Yields["Kminusmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Kminus5"]->fill(partPt);
					
					case 2212:	//proton
						
						hAUAU_Yields["Protonsmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Protons5"]->fill(partPt);			

					case -2212:	//antiproton
						
						hAUAU_Yields["Pbarmin"]->fill(partPt, pt_weight);
						hAUAU_Yields["Pbar5"]->fill(partPt);
				}
			}
		}
	}


	void finalize() {

		//____Yields____
		hAUAU_Yields["Piplusmin"]->scaleW(1./sow["AUAU_min"]->sumW());	//minimum bias centrality
		hAUAU_Yields["Piminusmin"]->scaleW(1./sow["AUAU_min"]->sumW());
		hAUAU_Yields["Kplusmin"]->scaleW(1./sow["AUAU_min"]->sumW());
		hAUAU_Yields["Kminusmin"]->scaleW(1./sow["AUAU_min"]->sumW());
		hAUAU_Yields["Protonsmin"]->scaleW(1./sow["AUAU_min"]->sumW());
		hAUAU_Yields["Pbarmin"]->scaleW(1./sow["AUAU_min"]->sumW());

		hAUAU_Yields["Piplus5"]->scaleW(1./sow["AUAU_5"]->sumW());	//0-5% centrality
		hAUAU_Yields["Piminus5"]->scaleW(1./sow["AUAU_5"]->sumW());
		hAUAU_Yields["Kplus5"]->scaleW(1./sow["AUAU_5"]->sumW());
		hAUAU_Yields["Kminus5"]->scaleW(1./sow["AUAU_5"]->sumW());
		hAUAU_Yields["Protons5"]->scaleW(1./sow["AUAU_5"]->sumW());
		hAUAU_Yields["Pbar5"]->scaleW(1./sow["AUAU_5"]->sumW());

		hAUAU_Yields["Piplus10"]->scaleW(1./sow["AUAU_10"]->sumW());	//5-10% centrality
		hAUAU_Yields["Piminus10"]->scaleW(1./sow["AUAU_10"]->sumW());
		hAUAU_Yields["Kplus10"]->scaleW(1./sow["AUAU_10"]->sumW());
		hAUAU_Yields["Kminus10"]->scaleW(1./sow["AUAU_10"]->sumW());
		hAUAU_Yields["Protons10"]->scaleW(1./sow["AUAU_10"]->sumW());
		hAUAU_Yields["Pbar10"]->scaleW(1./sow["AUAU_10"]->sumW());

		hAUAU_Yields["Piplus15"]->scaleW(1./sow["AUAU_15"]->sumW());	//10-15% centrality
		hAUAU_Yields["Piminus15"]->scaleW(1./sow["AUAU_15"]->sumW());
		hAUAU_Yields["Kplus15"]->scaleW(1./sow["AUAU_15"]->sumW());
		hAUAU_Yields["Kminus15"]->scaleW(1./sow["AUAU_15"]->sumW());
		hAUAU_Yields["Protons15"]->scaleW(1./sow["AUAU_15"]->sumW());
		hAUAU_Yields["Pbar15"]->scaleW(1./sow["AUAU_15"]->sumW());

		hAUAU_Yields["Piplus20"]->scaleW(1./sow["AUAU_20"]->sumW());	//15-20% centrality
		hAUAU_Yields["Piminus20"]->scaleW(1./sow["AUAU_20"]->sumW());
		hAUAU_Yields["Kplus20"]->scaleW(1./sow["AUAU_20"]->sumW());
		hAUAU_Yields["Kminus20"]->scaleW(1./sow["AUAU_20"]->sumW());
		hAUAU_Yields["Protons20"]->scaleW(1./sow["AUAU_20"]->sumW());
		hAUAU_Yields["Pbar20"]->scaleW(1./sow["AUAU_20"]->sumW());

		hAUAU_Yields["Piplus30"]->scaleW(1./sow["AUAU_30"]->sumW());	//20-30% centrality
		hAUAU_Yields["Piminus30"]->scaleW(1./sow["AUAU_30"]->sumW());
		hAUAU_Yields["Kplus30"]->scaleW(1./sow["AUAU_30"]->sumW());
		hAUAU_Yields["Kminus30"]->scaleW(1./sow["AUAU_30"]->sumW());
		hAUAU_Yields["Protons30"]->scaleW(1./sow["AUAU_30"]->sumW());
		hAUAU_Yields["Pbar30"]->scaleW(1./sow["AUAU_30"]->sumW());

		hAUAU_Yields["Piplus40"]->scaleW(1./sow["AUAU_40"]->sumW());	//30-40% centrality
		hAUAU_Yields["Piminus40"]->scaleW(1./sow["AUAU_40"]->sumW());
		hAUAU_Yields["Kplus40"]->scaleW(1./sow["AUAU_40"]->sumW());
		hAUAU_Yields["Kminus40"]->scaleW(1./sow["AUAU_40"]->sumW());
		hAUAU_Yields["Protons40"]->scaleW(1./sow["AUAU_40"]->sumW());
		hAUAU_Yields["Pbar40"]->scaleW(1./sow["AUAU_40"]->sumW());

		hAUAU_Yields["Piplus50"]->scaleW(1./sow["AUAU_50"]->sumW());	//40-50% centrality
		hAUAU_Yields["Piminus50"]->scaleW(1./sow["AUAU_50"]->sumW());
		hAUAU_Yields["Kplus50"]->scaleW(1./sow["AUAU_50"]->sumW());
		hAUAU_Yields["Kminus50"]->scaleW(1./sow["AUAU_50"]->sumW());
		hAUAU_Yields["Protons50"]->scaleW(1./sow["AUAU_50"]->sumW());
		hAUAU_Yields["Pbar50"]->scaleW(1./sow["AUAU_50"]->sumW());

		hAUAU_Yields["Piplus60"]->scaleW(1./sow["AUAU_60"]->sumW());	//50-60% centrality
		hAUAU_Yields["Piminus60"]->scaleW(1./sow["AUAU_60"]->sumW());
		hAUAU_Yields["Kplus60"]->scaleW(1./sow["AUAU_60"]->sumW());
		hAUAU_Yields["Kminus60"]->scaleW(1./sow["AUAU_60"]->sumW());
		hAUAU_Yields["Protons60"]->scaleW(1./sow["AUAU_60"]->sumW());
		hAUAU_Yields["Pbar60"]->scaleW(1./sow["AUAU_60"]->sumW());

		hAUAU_Yields["Piplus70"]->scaleW(1./sow["AUAU_70"]->sumW());	//60-70% centrality
		hAUAU_Yields["Piminus70"]->scaleW(1./sow["AUAU_70"]->sumW());
		hAUAU_Yields["Kplus70"]->scaleW(1./sow["AUAU_70"]->sumW());
		hAUAU_Yields["Kminus70"]->scaleW(1./sow["AUAU_70"]->sumW());
		hAUAU_Yields["Protons70"]->scaleW(1./sow["AUAU_70"]->sumW());
		hAUAU_Yields["Pbar70"]->scaleW(1./sow["AUAU_70"]->sumW());

		hAUAU_Yields["Piplus80"]->scaleW(1./sow["AUAU_80"]->sumW());	//70-80% centrality
		hAUAU_Yields["Piminus80"]->scaleW(1./sow["AUAU_80"]->sumW());
		hAUAU_Yields["Kplus80"]->scaleW(1./sow["AUAU_80"]->sumW());
		hAUAU_Yields["Kminus80"]->scaleW(1./sow["AUAU_80"]->sumW());
		hAUAU_Yields["Protons80"]->scaleW(1./sow["AUAU_80"]->sumW());
		hAUAU_Yields["Pbar80"]->scaleW(1./sow["AUAU_80"]->sumW());

		hAUAU_Yields["Piplus92"]->scaleW(1./sow["AUAU_92"]->sumW());	//80-92% centrality
		hAUAU_Yields["Piminus92"]->scaleW(1./sow["AUAU_92"]->sumW());
		hAUAU_Yields["Kplus92"]->scaleW(1./sow["AUAU_92"]->sumW());
		hAUAU_Yields["Kminus92"]->scaleW(1./sow["AUAU_92"]->sumW());
		hAUAU_Yields["Protons92"]->scaleW(1./sow["AUAU_92"]->sumW());
		hAUAU_Yields["Pbar92"]->scaleW(1./sow["AUAU_92"]->sumW());

		hAUAU_Yields["Piplus60_92"]->scaleW(1./sow["AUAU_60_92"]->sumW());	//60-92% centrality
		hAUAU_Yields["Piminus60_92"]->scaleW(1./sow["AUAU_60_92"]->sumW());
		hAUAU_Yields["Kplus60_92"]->scaleW(1./sow["AUAU_60_92"]->sumW());
		hAUAU_Yields["Kminus60_92"]->scaleW(1./sow["AUAU_60_92"]->sumW());
		hAUAU_Yields["Protons60_92"]->scaleW(1./sow["AUAU_60_92"]->sumW());
		hAUAU_Yields["Pbar60_92"]->scaleW(1./sow["AUAU_60_92"]->sumW());


		//____Rcp____

	};

	map<string, Histo1DPtr> hAUAU_Yields;
	map<string, CounterPtr> sow;

}
