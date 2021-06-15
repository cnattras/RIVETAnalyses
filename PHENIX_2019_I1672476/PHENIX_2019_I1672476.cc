// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include <math.h>
#define _USE_MATH_DEFINES
namespace Rivet {

/// @brief Add a short analysis description here
class PHENIX_2019_I1672476 : public Analysis {
public:

	/// Constructor
	DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2019_I1672476);
    

	/// Book histograms and initialise projections before the run
	void init() {
		declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
      
		const PromptFinalState pfs(Cuts::abseta < 0.35 && Cuts::pid == 22);
		declare(pfs, "pfs");
		const FinalState fs(Cuts::abseta < 0.35 && Cuts::abscharge > 0);
		declare(fs, "fs");
      

		book(_h["AuAu62_c0-20"], 1, 1, 1);
		book(_h["AuAu62_c0-86"], 1, 1, 2);
		book(_h["AuAu39_c0-86"], 2, 1, 1);
		book(_p["AuAu39_chPMult"], 3, 1, 1);
		book(_p["AuAu62_chPMult"], 4, 1, 1);
		book(_p["AuAu200_chPMult"], 4, 1, 2);      
		// book(_h["fig2-2-a"], 5, 1, 1);
		// book(_h["fig2-2-b"], 5, 1, 2);
		// book(_h["fig2-2-c"], 5, 1, 3);
		book(_h["AuAu62_chPMult_c0-40"], 6, 1, 1);
		book(_h["AuAu62_chPMult_c0-86"], 7, 1, 1);
		book(_h["AuAu39_chPMult_c0-86"], 8, 1, 1);
		book(_h["CuCu200_chPMult_c0-40"], 9, 1, 1);
		book(_h["CuCu200_chPMult_c0-40"], 10, 1, 1);
		book(_h["pp200_chPMult"], 11, 1, 1);
		book(_h["fig3-2a"], 12, 1, 1);
		book(_h["fig3-2b"], 13, 1, 1);
		book(_h["fig3-2c"], 14, 1, 1);
		book(_h["fig3-2d"], 15, 1, 1);
		book(_h["fig3-2e"], 16, 1, 1);
		book(_h["fig3-2f"], 17, 1, 1);
		book(_h["fig4-1a"], 18, 1, 1);
		book(_h["fig4-1b"], 19, 1, 1);
		book(_h["fig4-1c"], 20, 1, 1);
		// book(_h["fig4-2-a"], 21, 1, 1);
		//
		/*
		book(sow["sow-fig1-1-a"],"sow-fig1-1-a");//these are currently unused counters
		book(sow["sow-fig1-1-b"],"sow-fig1-1-b");
		book(sow["sow-fig1-2"],"sow-fig1-2");
		book(sow["sow-fig2-1a"],"sow-fig2-1a");
		book(sow["sow-fig2-1b-a"],"sow-fig2-1b-a");
		book(sow["sow-fig2-1b-b"],"sow-fig2-1b-b");
		book(sow["sow-fig3-1a"],"sow-fig3-1a");
		book(sow["sow-fig3-1b"],"sow-fig3-1b");
		book(sow["sow-fig3-1c"],"sow-fig3-1c");
		book(sow["sow-fig3-1d"],"sow-fig3-1d");
		book(sow["sow-fig3-1e"],"sow-fig3-1e");
		book(sow["sow-fig3-1f"],"sow-fig3-1f");
		book(sow["sow-fig3-2a"],"sow-fig3-2a");
		book(sow["sow-fig3-2b"],"sow-fig3-2b");
		book(sow["sow-fig3-2c"],"sow-fig3-2c");
		book(sow["sow-fig3-2d"],"sow-fig3-2d");
		book(sow["sow-fig3-2e"],"sow-fig3-2e");
		book(sow["sow-fig3-2f"],"sow-fig3-2f");
		book(sow["sow-fig4-1a"],"sow-fig4-1a");
		book(sow["sow-fig4-1b"],"sow-fig4-1b");
		book(sow["sow-fig4-1c"],"sow-fig4-1c");
		*/
}


	void analyze(const Event& event) {

		const PromptFinalState pfs = apply<PromptFinalState>(event, "pfs");
		const FinalState fs = apply<FinalState>(event, "fs");
		const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
		const double c = cent();
		const ParticlePair& beam = beams();
 
		int NN = 0;


		if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
		{
			NN = 197.;

			if (fuzzyEquals(sqrtS()/GeV, 39*NN, 1E-3)) collSystem = AuAu39;
			if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 1E-3)) collSystem = AuAu62;
			if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSystem = AuAu200;
		}

		if (beam.first.pid() == 1000290640 && beam.second.pid() == 1000290640)
		{
			NN = 64.;

			if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSystem = CuCu200;
		}
		
		if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
		{
			if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSystem = pp200;
		}

	

		Particles photons = applyProjection<PromptFinalState>(event, "pfs").particles();
		Particles chargedParticles = applyProjection<FinalState>(event, "fs").particles();



	if(collSystem == AuAu62)
	{
		if((c >= 0.) && (c < 20.))
		{
			for(const Particle& p : photons)
			{
				double partPt = p.pT()/GeV;
				double pt_weight = 1./(partPt*2.*M_PI);
				_h["AuAu62_c0-20"]->fill(partPt, pt_weight); 
				_h["AuAu62_c0-86"]->fill(partPt, pt_weight);
			}
		}

		else if((c >= 0.) && (c < 86. ))
		{
			for(const Particle& p : photons)
			{
				double partPt = p.pT()/GeV;
				double pt_weight = 1./(partPt*2.*M_PI);
				_h["AuAu62_c0-86"]->fill(partPt, pt_weight);
			}
		}

	}

	if (collSystem == AuAu39)
	{
		if((c >= 0.) && (c <86.))
		{
			for(const Particle& p : photons)
			{
				double partPt = p.pT()/GeV;
				double pt_weight = 1./(partPt*2.*M_PI);
				_h["AuAu39_c0-86"]->fill(partPt, pt_weight);
			}
		}
	}
     

	int nAuAu39 = 0;
	int nAuAu62 = 0;
	int nAuAu200 = 0;
	int nCuCu200 = 0;
	int npp200 = 0;
	int absEta = 0.7;
       

	//count charged particles for diff AuAu enegeries
	for(const Particle& p : chargedParticles)
	{
		if (collSystem == AuAu39)
		{
			nAuAu39 ++;   
		}

		else if(collSystem == AuAu62)
		{
			nAuAu62 ++;
		}

		else if(collSystem == AuAu200)
		{
			nAuAu200 ++;
		}

		else if(collSystem == CuCu200)
		{
			nCuCu200 ++;
		}

		else if(collSystem == pp200)
		{
			npp200 ++;
		}
	}
      

	//fill histos with cent c and charged particle count/abs val of Eta
	if (collSystem == AuAu39)
	{
		_p["AuAu39_chPMult"]->fill(c,nAuAu39/absEta);

		if((c >= 0.) && (c < 86.))
		{
			_p["AuAu39_chPMult_c0-86"]->fill(c,nAuAu39/absEta);
		}
	}

	else if(collSystem == AuAu62)
	{
		_p["AuAu62_chPMult"]->fill(c,nAuAu62/absEta);

		if((c >= 0.) && (c < 40.))
		{
			_p["AuAu62_chPMult_c0-40"]->fill(c,nAuAu62/absEta);
		}

		if((c >= 0.) && (c < 86.))
		{
			_p["AuAu62_chPMult_c0-86"]->fill(c,nAuAu62/absEta);
		}
	}

	else if(collSystem == AuAu200)
	{
		_p["AuAu200_chPMult"]->fill(c,nAuAu200/absEta);
	}
	
	else if(collSystem == CuCu200)
	{
		if((c >= 0.) && (c < 40.))
		{
			_p["CuCu200_chPMult_c0-40"]->fill(c,nCuCu200/absEta);
		}
		
		if((c >= 0.) && (c < 94.))
		{
			_p["CuCu200_chPMult_c0-94"]->fill(c,nCuCu200/absEta);

		}
	}
	
	else if(collSystem == pp200)
	{
		_p["pp200_chPMult"]->fill(c,npp200/absEta);
	}





 


}


	/// Normalise histograms etc., after the run
	void finalize() {

		normalize(_h["XXXX"]); // normalize to unity
		normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
		scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)
	}



	/// @name Histograms
	//@{
	map<string, Histo1DPtr> _h;
	map<string, Profile1DPtr> _p;
	map<string, CounterPtr> sow;
	enum CollisionSystem {AuAu39, AuAu62, AuAu200, CuCu200, pp200};
	CollisionSystem collSystem;
	string beamOpt;




	};


	DECLARE_RIVET_PLUGIN(PHENIX_2019_I1672476);

}
