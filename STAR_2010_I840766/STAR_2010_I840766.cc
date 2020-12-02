// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Centralities/RHICCentrality.hh"


// data: /phenix/scratch/cen/sampleHepData/ 
namespace Rivet {


	/// @brief Add a short analysis description here
	class STAR_2010_I840766 : public Analysis {
	public:

		/// Constructor
		DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2010_I840766);


		/// @name Analysis methods
		///@{

		/// Book histograms and initialise projections before the run
		void init() {
			beamOpt = getOption<string>("beam", "NONE");
			declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

			const FinalState fs(Cuts::abseta < 4.9);
			declare(fs, "fs");

			const PromptFinalState pfs(Cuts::abseta < 0.5);
			declare(pfs, "pfs");

			const UnstableParticles up(Cuts::abseta < 0.5 && Cuts::abspid == 111); // 111 = pi0, 221 = eta, 22 =? direct photon
			declare(up, "up");

			//book(_h["Table1"], 1, 1, 2);

			string refname = mkAxisCode(1, 1, 2);
			const Scatter2D& refdata = refData(refname);
			book(_h["AuAu"], refname + "_AuAu", refdata);
			book(_h["pp"], refname + "_pp", refdata);
			//book(_s["Raa"], refname);



			book(_h["Table2"], 2, 1, 2);
			book(_h["Table3"], 3, 1, 2);
			book(_h["Table4"], 4, 1, 2);
			book(_h["Table5"], 5, 1, 2);
			book(_h["Table6"], 6, 1, 2);
			book(_h["Table7"], 7, 1, 2);
			book(_h["Table8"], 8, 1, 2);
			book(_h["Table9"], 9, 1, 2);
			book(_h["Table10"], 10, 1, 2);
			book(_h["Table11"], 11, 1, 2);

			book(_c["sow_AuAu"], "sow_AuAu");
			book(_c["sow_pp"], "sow_pp");


		}


		/// Perform the per-event analysis
		void analyze(const Event& event) {

			Particles unParticles = applyProjection<UnstableParticles>(event,"up").particles();
			if(beamOpt=="pp")
			{
				_c["sow_pp"]->fill();
				for(const Particle& p : unParticles) 
				{
					if(p.pid() == 111) _h["pp"]->fill(p.pT()/GeV);
				}
				return;
			}
			else if(beamOpt=="AuAu")
			{
				_c["sow_AuAu"]->fill();
				for(const Particle& p : unParticles) 
				{
						if(p.pid() == 111) _h["AuAu"]->fill(p.pT()/GeV);
				}
			}
		}


		/// Normalise histograms etc., after the run
		void finalize() {

			bool pp_available = false;
		   bool AuAu_available = false;
		        
      	for (auto element : _c)
      	{
				string name = element.second->name();
       		//cout<<"Name is "<<name<<endl;
        		if (name.find("AuAu") != std::string::npos)
        		{
       			if (element.second->sumW()>0) AuAu_available=true;
       			else
       			{
         				AuAu_available=false;
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
			if((!pp_available) || (!AuAu_available)) return;
        	//cout<<"PP sum: "<<_c["sow_pp"]->sumW()<<endl;
			if(_c["sow_AuAu"]->sumW() > 0)
				_h["AuAu"]->scaleW(1./_c["sow_AuAu"]->sumW());
			if(_c["sow_pp"]->sumW() > 0)
				_h["pp"]->scaleW(1./_c["sow_pp"]->sumW());

		}

		///@}


		/// @name Histograms
		///@{
		map<string, Histo1DPtr> _h;
		map<string, Profile1DPtr> _p;
		map<string, CounterPtr> _c;
		map<string, Scatter2DPtr> _s;
		string beamOpt = "";
		///@}


	};


	DECLARE_RIVET_PLUGIN(STAR_2010_I840766);

}
