// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "RHICCentrality.hh"
#include <math.h>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES




// data: /phenix/scratch/cen/sampleHepData/ 
namespace Rivet {


	/// @brief Add a short analysis description here
	class STAR_2010_I840766 : public Analysis {
	public:

		/// Constructor
		RIVET_DEFAULT_ANALYSIS_CTOR(STAR_2010_I840766);


		/// @name Analysis methods
		///@{

		/// Book histograms and initialise projections before the run
		void init() {
			beamOpt = getOption<string>("beam", "NONE");
			declareCentrality(RHICCentrality("STAR"), "RHIC_2019_CentralityCalibration:exp=STAR", "CMULT", "CMULT");

			const FinalState fs(Cuts::abseta < 1.0);
			declare(fs, "fs");

			const PromptFinalState pfs(Cuts::abseta < 1.0);
			declare(pfs, "pfs");

			const UnstableParticles up(Cuts::abseta < 1.0); 
			declare(up, "up");

			book(_h["Table1"], 1, 1, 2);
			book(_h["Table2"], 2, 1, 2);

			//book(_h["Table3"], 3, 1, 2);
			string refname = mkAxisCode(3, 1, 2);
			const Estimate1D& EtaPi_pp = refData(refname);
			book(_h["Table3_eta"], refname + "_eta", EtaPi_pp);
			book(_h["Table3_pi"], refname + "_pi", EtaPi_pp);
			book(_s["Table3"], refname);

			//book(_h["Table4"], 4, 1, 2);
         refname = mkAxisCode(4, 1, 2);
         const Estimate1D& EtaPi_dAu = refData(refname);
         book(_h["Table4_eta"], refname + "_eta", EtaPi_dAu);
         book(_h["Table4_pi"], refname + "_pi", EtaPi_dAu);
         book(_s["Table4"], refname);


			//book(_h["Table5"], 5, 1, 2);
         refname = mkAxisCode(5, 1, 2); // pi0
         const Estimate1D& R_dAu_pi = refData(refname);
         book(_h["Table5_pp"], refname + "_pp", R_dAu_pi);
         book(_h["Table5_dAu"], refname + "_dAu", R_dAu_pi);
         book(_s["Table5"], refname);

			//book(_h["Table6"], 6, 1, 2);
         refname = mkAxisCode(6, 1, 2); // eta
         const Estimate1D& R_dAu_eta = refData(refname);
         book(_h["Table6_pp"], refname + "_pp", R_dAu_eta);
         book(_h["Table6_dAu"], refname + "_dAu", R_dAu_eta);
         book(_s["Table6"], refname);

			//book(_h["Table7"], 7, 1, 2);
         refname = mkAxisCode(7, 1, 2);
         const Estimate1D& R_CP_pi = refData(refname);
         book(_h["Table7_Central"], refname + "_C", R_CP_pi);
         book(_h["Table7_Peripheral"], refname + "_P", R_CP_pi);
         book(_s["Table7"], refname);



			//book(_h["Table8"], 8, 1, 2);
         refname = mkAxisCode(8, 1, 2);
         const Estimate1D& R_gamma_pp = refData(refname);
         book(_h["Gamma_inclusive_pp"], refname + "_inclusive_pp", R_gamma_pp);
         book(_h["Gamma_decay_pp"], refname + "_decay_pp", R_gamma_pp);
         book(_h["Gamma_prompt_pp"], refname + "_prompt_pp", R_gamma_pp);
         book(_s["Table8"], refname);

         //book(_h["Table9"], 9, 1, 2);
         refname = mkAxisCode(9, 1, 2);
         const Estimate1D& R_gamma_dAu = refData(refname);
         book(_h["Gamma_inclusive_dAu"], refname + "_inclusive_dAu", R_gamma_dAu);
         book(_h["Gamma_decay_dAu"], refname + "_decay_dAu", R_gamma_dAu);
         book(_h["Gamma_prompt_dAu"], refname + "_prompt_dAu", R_gamma_dAu);
         book(_s["Table9"], refname);


			book(_h["Table10"], 10, 1, 2);
			book(_h["Table11"], 11, 1, 2);

			book(_c["sow_dAu"], "sow_dAu");
			book(_c["sow_pp"], "sow_pp");


		}


		/// Perform the per-event analysis
		void analyze(const Event& event) {

			Particles upParticles = apply<UnstableParticles>(event,"up").particles();
         Particles promtParticles = apply<PromptFinalState>(event,"pfs").particles();
         Particles fsParticles = apply<FinalState>(event,"pfs").particles();

			if(beamOpt=="pp")
			{
				_c["sow_pp"]->fill();
				for(const Particle& p : upParticles) 
				{
					if(p.pid() == 111) // 111 = pi0, 221 = eta, 22 =? direct photon
               {
                  _h["Table1"]->fill(p.pT()/GeV); // cross section for pi0
                  _h["Table3_pi"]->fill(p.pT()/GeV);
                  _h["Table5_pp"]->fill(p.pT()/GeV);
               }  

               if(p.pid() == 221) // 111 = pi0, 221 = eta, 22 =? direct photon
               {
                  _h["Table3_eta"]->fill(p.pT()/GeV);
                  _h["Table6_pp"]->fill(p.pT()/GeV);
               }
            }
            for(const Particle& p : promtParticles) 
            {
               if(p.pid() == 22)
               {
                  _h["Table10"]->fill(p.pT()/GeV); // cross section for direct gamma
                  _h["Gamma_prompt_pp"]->fill(p.pT()/GeV);
               }

				}
            for(const Particle& p : fsParticles) 
            {
               if(p.pid() == 22)
               {
                  _h["Gamma_inclusive_pp"]->fill(p.pT()/GeV);
               }

            }
				return;
			}
			else if(beamOpt=="dAu")
			{
            const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
            const double c = cent();

				_c["sow_dAu"]->fill();
				for(const Particle& p : upParticles) 
				{
					if(p.pid() == 111)
               {   
                  _h["Table2"]->fill(p.pT()/GeV); // cross section for pi0
                  _h["Table4_pi"]->fill(p.pT()/GeV);
                  _h["Table5_dAu"]->fill(p.pT()/GeV);
                  if(c < 20)
                     _h["Table7_Central"]->fill(p.pT()/GeV);
                  else if(c > 40)
                     _h["Table7_Peripheral"]->fill(p.pT()/GeV);
               }
                  
               if(p.pid() == 221) // 111 = pi0, 221 = eta, 22 =? direct photon
               {
                  _h["Table4_eta"]->fill(p.pT()/GeV);
                  _h["Table6_dAu"]->fill(p.pT()/GeV);
               }
            }
            for(const Particle& p : promtParticles) 
            {
               if(p.pid() == 22)
               {
                  _h["Table11"]->fill(p.pT()/GeV); // cross section for direct gamma
                  _h["Gamma_prompt_dAu"]->fill(p.pT()/GeV);
               }
               
				}
            for(const Particle& p : fsParticles) 
            {
               if(p.pid() == 22)
               {
                  _h["Gamma_inclusive_dAu"]->fill(p.pT()/GeV);
               }

            }
			}
		}


		/// Normalise histograms etc., after the run
		void finalize() {

         bool pp_available = false;
         bool dAu_available = false;
		        
      	for (auto element : _c)
      	{
				string name = element.second->name();
       		//cout<<"Name is "<<name<<endl;
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
        	//cout<<"PP sum: "<<_c["sow_pp"]->sumW()<<endl;

			if(_c["sow_pp"]->sumW() > 0)
			{
            _h["Table1"]->scaleW(1./_c["sow_pp"]->sumW());
            _h["Table10"]->scaleW(1./_c["sow_pp"]->sumW());
         }	
			if(_c["sow_dAu"]->sumW() > 0)
         {
				_h["Table2"]->scaleW(1./_c["sow_dAu"]->sumW());
            _h["Table11"]->scaleW(1./_c["sow_dAu"]->sumW());
         }

         divide(_h["Table3_eta"], _h["Table3_pi"], _s["Table3"]);
         divide(_h["Table4_eta"], _h["Table4_pi"], _s["Table4"]);

         divide(_h["Table5_dAu"], _h["Table5_pp"], _s["Table5"]);
         divide(_h["Table6_pp"], _h["Table6_dAu"], _s["Table6"]);
         divide(_h["Table7_Central"], _h["Table7_Peripheral"], _s["Table7"]);

         *_h["Gamma_decay_pp"] = *_h["Gamma_inclusive_pp"] - *_h["Gamma_prompt_pp"];
         *_h["Gamma_decay_dAu"] = *_h["Gamma_inclusive_dAu"] - *_h["Gamma_prompt_dAu"];

         divide(_h["Gamma_inclusive_pp"], _h["Gamma_decay_pp"], _s["Table8"]);
         divide(_h["Gamma_inclusive_dAu"], _h["Gamma_decay_dAu"], _s["Table9"]);



		}

		///@}


		/// @name Histograms
		///@{
		map<string, Histo1DPtr> _h;
		map<string, Profile1DPtr> _p;
		map<string, CounterPtr> _c;
		map<string, Estimate1DPtr> _s;
		string beamOpt = "";
		///@}


	};


	RIVET_DECLARE_PLUGIN(STAR_2010_I840766);

}
