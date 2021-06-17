// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
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
		book(_p["AuAu39_chPMult_c0-55"], 3, 1, 1);
		book(_p["AuAu62_chPMult_c0-60"], 4, 1, 1);
		book(_p["AuAu200_chPMult_c0-60"], 4, 1, 2);      
		// book(_h["fig2-2-a"], 5, 1, 1);
		// book(_h["fig2-2-b"], 5, 1, 2);
		// book(_h["fig2-2-c"], 5, 1, 3);
		book(_p["AuAu62_chPMult_c0-40"], 6, 1, 1);
		book(_p["AuAu62_chPMult_c0-86"], 7, 1, 1);
		book(_p["AuAu39_chPMult_c0-86"], 8, 1, 1);
		book(_p["CuCu200_chPMult_c0-40"], 9, 1, 1);
		book(_p["CuCu200_chPMult_c0-94"], 10, 1, 1);
		book(_p["pp200_chPMult"], 11, 1, 1);
		book(_h["AuAu62_phtYld_c0-40"], 12, 1, 1);
		book(_h["AuAu62_phtYld_c0-86"], 13, 1, 1);
		book(_h["AuAu39_phtYld_c0-86"], 14, 1, 1);
		book(_h["CuCu200_phtYld_c0-40"], 15, 1, 1);
		book(_h["CuCu200_phtYld_c0-94"], 16, 1, 1);
		book(_h["pp200_phtYld"], 17, 1, 1);
		book(_p["AuAu200_chPMult_c0-92"], 18, 1, 1);
		//book(_h["pp200_chPMult_c0-100"], 19, 1, 1); exact same as 3-1f
		book(_p["pp62_chPMult"], 20, 1, 1);
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
			if (fuzzyEquals(sqrtS()/GeV, 62*NN, 1E-3)) collSystem = pp62;
		}

	

		const Particles photons = pfs.particles(Cuts::pT > 1.*GeV && Cuts::pT < 5.*GeV);
		const Particles chargedParticles = fs.particles();



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
     

	int absEta = 0.7;

      

	//fill histos with cent c and charged particle count/abs val of Eta
	if (collSystem == AuAu39)
	{
		if((c >= 0.) && (c < 55.))
		{
			_p["AuAu39_chPMult_c0-55"]->fill(c,chargedParticles.size()/absEta);
		}

		if((c >= 0.) && (c < 86.))
		{
			_p["AuAu39_chPMult_c0-86"]->fill(c,chargedParticles.size()/absEta);
			_p["AuAu39_phtYld_c0-86"]->fill(c,photons.size()/absEta);
		}
	}

	else if(collSystem == AuAu62)
	{
		if((c >= 0.) && (c < 60.))
		{
			_p["AuAu62_chPMult_c0-60"]->fill(c,chargedParticles.size()/absEta);
		}

		if((c >= 0.) && (c < 40.))
		{
			_p["AuAu62_chPMult_c0-40"]->fill(c,chargedParticles.size()/absEta);
			_p["AuAu62_phtYld_c0-40"]->fill(c,photons.size()/absEta);
		}

		if((c >= 0.) && (c < 86.))
		{
			_p["AuAu62_chPMult_c0-86"]->fill(c,chargedParticles.size()/absEta);
			_p["AuAu62_phtYld_c0-86"]->fill(c,photons.size()/absEta);
		}
	}

	else if(collSystem == AuAu200)
	{
		if((c >= 0.) && (c < 60.))
		{
			_p["AuAu200_chPMult_c0-60"]->fill(c,chargedParticles.size()/absEta);
		
		}
		
		if((c >= 0.) && (c < 92.))
		{
			_p["AuAu200_chPMult_c0-92"]->fill(c,chargedParticles.size()/absEta);
		}
	}
	
	else if(collSystem == CuCu200)
	{
		if((c >= 0.) && (c < 40.))
		{
			_p["CuCu200_chPMult_c0-40"]->fill(c,chargedParticles.size()/absEta);
			_p["CuCu200_phtYld_c0-40"]->fill(c,photons.size()/absEta);
		}
		
		if((c >= 0.) && (c < 94.))
		{
			_p["CuCu200_chPMult_c0-94"]->fill(c,chargedParticles.size()/absEta);
			_p["CuCu200_phtYld_c0-40"]->fill(c,photons.size()/absEta);

		}
	}
	
	else if(collSystem == pp200)
	{
		_p["pp200_chPMult"]->fill(c,chargedParticles.size()/absEta);
		_p["pp200_phtYld"]->fill(c,photons.size()/absEta);
	}

	else if(collSystem == pp62)
	{
		_p["pp62_chPMult"]->fill(c,chargedParticles.size()/absEta);
	}



 


	}


	/// Normalise histograms etc., after the run
	void finalize() {

		normalize(_h["XXXX"]); // normalize to unity
		normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
		scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)
	}



	/// @name Histograms
	
	map<string, Histo1DPtr> _h;
	map<string, Profile1DPtr> _p;
	map<string, CounterPtr> sow;
	enum CollisionSystem {AuAu39, AuAu62, AuAu200, CuCu200, pp200, pp62};
	CollisionSystem collSystem;
	string beamOpt;
	

	

	};


	DECLARE_RIVET_PLUGIN(PHENIX_2019_I1672476);

}
