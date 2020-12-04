// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "RHICCentrality.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class STAR_2016_I1420183 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2016_I1420183);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

		//declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

		// const FinalState fs(Cuts::abseta < 0.5 && Cuts::pT > 0.15*GeV);
		// declare(fs, "fs");

		// const PromptFinalState pfs(Cuts::abseta < 0.5);
		// declare(pfs, "pfs");

		// const UnstableParticles up(Cuts::abseta < 0.5 && Cuts::abspid == 421);
		// declare(up, "up");

		// book(_h["AAAA"], 1, 1, 1);
		// book(_c["sow_AuAu0010"], "sow_AuAu0010");
		// book(_c["sow_AuAu3050"], "sow_AuAu3050");
		// book(_c["sow_pp"], "sow_pp");
		
		// book(_p["XXXX"], 5, 1 ,1);
		// book(_p["YYYY"], "YYYY", 2, 0.0, 4.0);

		// book(_c["CCCC"], "CCCC");

		
		// string refname = mkAxisCode(5, 1, 1);
		// const Scatter2D& refdata = refData(refname);
		// book(_h["Kaon"], refname + "_Kaon", refdata);
		// book(_h["Pion"], refname + "_Pion", refdata);
		// book(_s["KaonOverPion"], refname);

		declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

		const FinalState fs(Cuts::abseta < 1.0 && Cuts::pT > 0.15*GeV);
		declare(fs, "fs");

		const PromptFinalState pfs(Cuts::absrap < 1.0);
		declare(pfs, "pfs");

		const UnstableParticles up(Cuts::absrap < 1.0);
		declare(up, "up");
		
		book(_h["particles"], 6, 1, 1);
		book(_h["Jpsi"], 5, 1, 1);
		book(_c["eventW"], "eventW");
		
		beamOpt = getOption<string>("beam", "NONE");

		string refname = mkAxisCode(1, 1, 1);
		const Scatter2D& refdata = refData(refname);
		book(_h["dAu"], refname + "_dAu", refdata);
		book(_h["pp"], refname + "_pp", refdata);
		book(_s["Raa"], refname);

		book(_c["dAu"], "sow_dAu");
		book(_c["pp"], "sow_pp");
		
		book(_c["sow_dAu0010"], "sow_dAu0010");


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
		
		const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
		const double c = cent();

		// if(c < 10.)
		// {
		// _c["sow_AuAu0010"]->fill();
		// }
		// else if(c >= 30. && c < 50.)
		// {
		// _c["sow_AuAu3050"]->fill();
		// }
		// if(c > 50.) vetoEvent;

		// Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();

		// for(const Particle& p : fsParticles) 
		// {
			// if(c < 10. && p.pid() == 321) _h["AAAA"]->fill(p.pT()/GeV);
		// }
		
		_c["eventW"]->fill();
		
		Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();
		Particles pfsParticles = applyProjection<PromptFinalState>(event,"pfs").particles();
		
		for(const Particle& fs : fsParticles) 
		{
			_h["particles"]->fill(fs.pT()/GeV);
		}

		Particles upParticles = applyProjection<UnstableParticles>(event,"up").particles();

		for(const Particle& p : upParticles) 
		{
			if(p.pid() == 421) _h["Jpsi"]->fill(p.pT()/GeV);
		}
		
		if(beamOpt=="pp")
		{
			_c["sow_pp"]->fill();
			for(const Particle& p : upParticles) 
			{
				if(p.pid() == 421) _h["pp"]->fill(p.pT()/GeV);
			}
			return;
		}
		else if(beamOpt=="dAu")
		{
			//const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
			//const double c = cent();

			//if(c > 10) vetoEvent;
			_c["sow_dAu0010"]->fill();
			
			for(const Particle& p : upParticles) 
			{
				if(p.pid() == 421) _h["dAu"]->fill(p.pT()/GeV);
			}
		}



    }


    /// Normalise histograms etc., after the run
    void finalize() {

		// double scale = 1./(2*M_PI);
		// _h["AAAA"]->scaleW(scale/_c["CCCC"]->sumW());
		// _h["Kaon"]->scaleW(1./_c["CCCC"]->sumW());
		// _h["Pion"]->scaleW(1./_c["CCCC"]->sumW());
		// divide(_h["Kaon"], _h["Pion"], _s["KaonOverPion"]);
		
		
		bool pp_available = false;
		bool dAu_available = false;
					
		for (auto element : _c)
		{
		string name = element.second->name();
			if (name.find("dAu") != std::string::npos)
			{
				  if (element.second->sumW()>0) dAu_available=true;
				  else
				  {
						dAu_available=false;
						break;
				  }
			}
			else if (name.find("pp") != std::string::npos)
			{
				  if (element.second->sumW()>0) pp_available=true;
				  else
				  {
						pp_available=false;
						break;
				  }
			}
		}
		if((!pp_available) || (!dAu_available)) return;	
		
		
		
		double scale = 1./(2*M_PI);
		_h["Jpsi"]->scaleW(scale/_c["eventW"]->sumW());
		//_h["Jpsi"]->scaleW(scale);
		
		
    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
	map<string, Scatter2DPtr> _s;
	string beamOpt = "";
    //@}

  };


  DECLARE_RIVET_PLUGIN(STAR_2016_I1420183);

}
