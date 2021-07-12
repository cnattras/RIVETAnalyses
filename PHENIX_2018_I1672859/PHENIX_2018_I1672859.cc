// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "../Centralities/RHICCentrality.hh"
#include <math.h>
#define _USE_MATH_DEFINES
const double TAB[6] = {2.54, 8.8, 6.0, 7.5, 3.1, 1.0};

namespace Rivet {


/// @brief Add a short analysis description here
class PHENIX_2018_I1672859 : public Analysis {
public:

	/// Constructor
	DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2018_I1672859);

	/// Book histograms and initialise projections before the run
	void init()
	{
		//~ eta pdg id: 221
		//~ pi0 pdg id: 111
		beamOpt = getOption<string>("beam", "NONE");      
		declareCentrality(RHICCentrality("PHENIX"),"RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

		// Initialise and register projections
		// The basic final-state projection: all final-state particles within the given eta acceptance
		const FinalState fs(Cuts::abseta < 0.35);
		declare(fs, "fs");
      
		const UnstableParticles usp(Cuts::abseta < 0.35);
		declare(usp, "UFS");
		// declare(up, "up");

		// Book histograms specify custom binning
		// take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
		book(_h["pi0_minbias"], 1, 1, 1);
		book(_h["pi0_0010"],  2, 1, 1);
		book(_h["pi0_1020"],  3, 1, 1);
		book(_h["pi0_0020"],  4, 1, 1);
		book(_h["pi0_2040"],  5, 1, 1);
		book(_h["pi0_4060"],  6, 1, 1);
		book(_h["pi0_6090"],  7, 1, 1); //need to implement
      
		book(_h["eta_minbias"],  8, 1, 1);
		book(_h["eta_0020"],  9, 1, 1);
		book(_h["eta_2040"], 10, 1, 1);
		book(_h["eta_4060"], 11, 1, 1);
		book(_h["eta_6090"], 12, 1, 1); //need to implement
      
		//RAB plots
      
		string refname18 = mkAxisCode(18, 1, 1);
		const Scatter2D& refdata18 = refData(refname18);
		book(_h["CuAu_minbias_pi0"], refname18 + "_CuAu_minbias", refdata18);
		book(_h["pp_pi0"], refname18 + "_pp", refdata18);
		book(_s["R_AB_pi0"], refname18);

		string refname23 = mkAxisCode(23, 1, 1);
		const Scatter2D& refdata23 = refData(refname23);
		book(_h["CuAu_minbias_eta"], refname23 + "_CuAu_minbias", refdata23);
		book(_h["pp_eta"], refname23 + "_pp", refdata23);
		book(_s["R_AB_eta"], refname23);
      
		// Histos for ratio plots
 
		string refname13 = mkAxisCode(13, 1, 1);
 		const Scatter2D& refdata13 = refData(refname13);
		book(_h["eta_mb"], refname13 + "_eta_mb", refdata13);
		book(_h["pi0_mb"], refname13 + "_pi0_mb", refdata13);
		book(_s["eta_over_pi0_mb"], refname13); 

		string refname14 = mkAxisCode(14, 1, 1);
		const Scatter2D& refdata14 = refData(refname14);
		book(_h["eta_0p"], refname14 + "_eta_0p", refdata14);
		book(_h["pi0_0p"], refname14 + "_pi0_0p", refdata14);
		book(_s["eta_over_pi0_0p"], refname14);
      
		string refname15 = mkAxisCode(15, 1, 1);
		const Scatter2D& refdata15 = refData(refname15);
		book(_h["eta_20p"], refname15 + "_eta_20p", refdata15);
		book(_h["pi0_20p"], refname15 + "_pi0_20p", refdata15);
		book(_s["eta_over_pi0_20p"], refname15);
      
		string refname16 = mkAxisCode(16, 1, 1);
		const Scatter2D& refdata16 = refData(refname16);
		book(_h["eta_40p"], refname16 + "_eta_40p", refdata16);
		book(_h["pi0_40p"], refname16 + "_pi0_40p", refdata16);
		book(_s["eta_over_pi0_40p"], refname16);
      
		string refname17 = mkAxisCode(17, 1, 1);
		const Scatter2D& refdata17 = refData(refname17);
		book(_h["eta_60p"], refname17 + "_eta_60p", refdata17);
		book(_h["pi0_60p"], refname17 + "_pi0_60p", refdata17);
		book(_s["eta_over_pi0_60p"], refname17);


		book(sow["pp_pi0_N"],"pp_pi0_N");
		book(sow["pi0_0010_N"],"pi0_0010_N");
		book(sow["pi0_1020_N"],"pi0_1020_N");
    		book(sow["pi0_0020_N"],"pi0_0020_N");
		book(sow["pi0_2040_N"],"pi0_2040_N");
		book(sow["pi0_4060_N"],"pi0_4060_N");
		book(sow["pi0_6090_N"],"pi0_6090_N");
		book(sow["pi0_minbias_N"],"pi0_minbias_N");
		book(sow["pi0_mb_N"],"pi0_mb_N");
		book(sow["pi0_0p_N"],"pi0_0p_N");
		book(sow["pi0_20p_N"],"pi0_20p_N");
		book(sow["pi0_40p_N"],"pi0_40p_N");
		book(sow["pi0_60p_N"],"pi0_60p_N");
		book(sow["CuAu_minbias_pi0_N"],"CuAu_minbias_pi0_N");

		book(sow["pp_eta_N"],"pp_eta_N");
		book(sow["eta_0020_N"],"eta_0020_N");
		book(sow["eta_2040_N"],"eta_2040_N");
		book(sow["eta_4060_N"],"eta_4060_N");
		book(sow["eta_6090_N"],"eta_6090_N");
		book(sow["eta_minbias_N"],"eta_minbias_N");
		book(sow["eta_mb_N"],"eta_mb_N");
		book(sow["eta_0p_N"],"eta_0p_N");
		book(sow["eta_20p_N"],"eta_20p_N");
		book(sow["eta_40p_N"],"eta_40p_N");
		book(sow["eta_60p_N"],"eta_60p_N");
		book(sow["CuAu_minbias_eta_N"],"CuAu_minbias_eta_N");

 
 
	}


	/// Perform the per-event analysis
	void analyze(const Event& event)
	{
		const UnstableParticles& usp = applyProjection<UnstableParticles>(event, "UFS");
		const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
		const double c = cent();
      
		for(const Particle & p : usp.particles())
		{
			if(beamOpt == "CUAU200") 
			{
				if(c > 90.) vetoEvent;
				
				if(p.pid() == PID::PI0)
				{
            
					if(c <= 10.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["pi0_0010"]->fill(partPt, ptWeight);
						sow["pi0_0010_N"]->fill();		
					}
					if(c > 10 && c < 20.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["pi0_1020"]->fill(partPt, ptWeight);
						sow["pi0_1020_N"]->fill();
					}
					if(c > 0 && c < 20.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["pi0_0020"]->fill(partPt, ptWeight);
						sow["pi0_0020_N"]->fill();
						_h["pi0_0p"]->fill(partPt, ptWeight);
						sow["pi0_0p_N"]->fill();


					}
					if(c > 20 && c < 40.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["pi0_2040"]->fill(partPt, ptWeight);
						sow["pi0_2040_N"]->fill();
						_h["pi0_20p"]->fill(partPt, ptWeight);
						sow["pi0_20p_N"]->fill();


					}
					if(c > 40 && c < 60.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);			
						_h["pi0_4060"]->fill(partPt, ptWeight);
						sow["pi0_4060_N"]->fill();
						_h["pi0_40p"]->fill(partPt, ptWeight);
						sow["pi0_40p_N"]->fill();

					}
					if(c > 60 && c < 90.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);
						_h["pi0_6090"]->fill(partPt, ptWeight);
						sow["pi0_6090_N"]->fill();
						_h["pi0_60p"]->fill(partPt, ptWeight);
						sow["pi0_60p_N"]->fill();

					}

					double partPt = p.pT()/GeV;
					double ptWeight = 1./(partPt*2.*M_PI);
					_h["pi0_minbias"]->fill(partPt, ptWeight);
					_h["pi0_mb"]->fill(partPt, ptWeight);
					_h["CuAu_minbias_pi0"]->fill(partPt, ptWeight);
					sow["pi0_minbias_N"]->fill();
					sow["pi0_mb_N"]->fill();
					sow["CuAu_minbias_pi0_N"]->fill();
				}
          
				else if(p.pid() == PID::ETA)
				{
            
					if(c > 0 && c < 20.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["eta_0020"]->fill(partPt, ptWeight);
						sow["eta_0020_N"]->fill();
						_h["eta_0p"]->fill(partPt, ptWeight);
						sow["eta_0p_N"]->fill();	
					}
					if(c > 20 && c < 40.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["eta_2040"]->fill(partPt, ptWeight);
						sow["eta_2040_N"]->fill();
						_h["eta_20p"]->fill(partPt, ptWeight);	
						sow["eta_20p_N"]->fill();

					}
					if(c > 40 && c < 60.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["eta_4060"]->fill(partPt, ptWeight);
						sow["eta_4060_N"]->fill();
						_h["eta_40p"]->fill(partPt, ptWeight);
						sow["eta_40p_N"]->fill();

					}
					if(c > 60 && c < 90.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);
						_h["eta_6090"]->fill(partPt, ptWeight);
						sow["eta_6090_N"]->fill();
						_h["eta_60p"]->fill(partPt, ptWeight);
						sow["eta_60p_N"]->fill();

					}

					double partPt = p.pT()/GeV;
					double ptWeight = 1./(partPt*2.*M_PI);	
					_h["eta_minbias"]->fill(partPt, ptWeight);
					_h["eta_mb"]->fill(partPt, ptWeight);	
					
					_h["CuAu_minbias_eta"]->fill(partPt, ptWeight);
					sow["eta_minbias_N"]->fill();
					sow["eta_mb_N"]->fill();
					sow["CuAu_minbias_eta_N"]->fill();
					
				}
			}        
			else if(beamOpt == "PP200")
			{          
				if(p.pid() == PID::PI0)
				{
					_h["pp_pi0"]->fill(p.pT()/GeV);
					sow["pp_pi0_N"]->fill();
					_h["pi0_0010"]->fill(p.pT()/GeV); //since this is pp do we need to fill for all centralities or non since it does not have centrality, in this case just pp_pi0
					_h["pi0_1020"]->fill(p.pT()/GeV);
					_h["pi0_0020"]->fill(p.pT()/GeV);
					_h["pi0_2040"]->fill(p.pT()/GeV);
					_h["pi0_4060"]->fill(p.pT()/GeV);
				}
				else
				{
					_h["pp_eta"]->fill(p.pT()/GeV);
					sow["pp_eta_N"]->fill();            
				//	if(c > 0 && c < 20.) _h["eta_0020"]->fill(p.pT()/GeV);
				//	if(c > 20 && c < 40.) _h["eta_2040"]->fill(p.pT()/GeV);
				//	if(c > 40 && c < 60.) _h["eta_4060"]->fill(p.pT()/GeV);
				}
			}
		}
      
	}


	/// Normalise histograms etc., after the run
	void finalize()
	{
/*		normalize(_h["pi0_minbias"]); // normalize to unity
		normalize(_h["pi0_0010"]); // normalize to unity
		normalize(_h["pi0_1020"]); // normalize to unity
		normalize(_h["pi0_0020"]); // normalize to unity
		normalize(_h["pi0_2040"]); // normalize to unity
		normalize(_h["pi0_4060"]); // normalize to unity
		
      
		normalize(_h["eta_minbias"]); // normalize to unity
		normalize(_h["eta_0020"]); // normalize to unity
		normalize(_h["eta_2040"]); // normalize to unity
		normalize(_h["eta_4060"]); // normalize to unity
*/

		//dividing yields by # of events N
	
		_h["pp_pi0"]->scaleW(1./(sow["pp_pi0_N"]->sumW()));	
		_h["pi0_0010"]->scaleW(1./(sow["pi0_0010_N"]->sumW()));	
		_h["pi0_1020"]->scaleW(1./(sow["pi0_1020_N"]->sumW()));
		_h["pi0_0020"]->scaleW(1./(sow["pi0_0020_N"]->sumW()));	
		_h["pi0_2040"]->scaleW(1./(sow["pi0_2040_N"]->sumW()));	
		_h["pi0_4060"]->scaleW(1./(sow["pi0_4060_N"]->sumW()));
		_h["pi0_6090"]->scaleW(1./(sow["pi0_6090_N"]->sumW()));	
		_h["pi0_minbias"]->scaleW(1./(sow["pi0_minbias_N"]->sumW()));
	
		_h["pi0_mb"]->scaleW(1./(sow["pi0_mb_N"]->sumW()));
		_h["pi0_0p"]->scaleW(1./(sow["pi0_0p_N"]->sumW()));
		_h["pi0_20p"]->scaleW(1./(sow["pi0_20p_N"]->sumW()));
		_h["pi0_40p"]->scaleW(1./(sow["pi0_40p_N"]->sumW()));
		_h["pi0_60p"]->scaleW(1./(sow["pi0_60p_N"]->sumW()));			
		_h["CuAu_minbias_pi0"]->scaleW(1./(sow["CuAu_minbias_pi0_N"]->sumW()));	


		_h["pp_eta"]->scaleW(1./(sow["pp_eta_N"]->sumW()));
		_h["eta_0020"]->scaleW(1./(sow["eta_0020_N"]->sumW()));		
		_h["eta_2040"]->scaleW(1./(sow["eta_2040_N"]->sumW()));	
		_h["eta_4060"]->scaleW(1./(sow["eta_4060_N"]->sumW()));
		_h["eta_6090"]->scaleW(1./(sow["eta_6090_N"]->sumW()));		
		_h["eta_minbias"]->scaleW(1./(sow["eta_minbias_N"]->sumW()));
	
		_h["eta_mb"]->scaleW(1./(sow["eta_mb_N"]->sumW()));
		_h["eta_0p"]->scaleW(1./(sow["eta_0p_N"]->sumW()));
		_h["eta_20p"]->scaleW(1./(sow["eta_20p_N"]->sumW()));
		_h["eta_40p"]->scaleW(1./(sow["eta_40p_N"]->sumW()));
		_h["eta_60p"]->scaleW(1./(sow["eta_60p_N"]->sumW()));			
		_h["CuAu_minbias_eta"]->scaleW(1./(sow["CuAu_minbias_eta_N"]->sumW()));	
			

		//ratio plots
      
		divide(_h["eta_mb"],_h["pi0_mb"],_s["eta_over_pi0_mb"]);
		divide(_h["eta_0p"],_h["pi0_0p"],_s["eta_over_pi0_0p"]);
		divide(_h["eta_20p"],_h["pi0_20p"],_s["eta_over_pi0_20p"]);
		divide(_h["eta_40p"],_h["pi0_40p"],_s["eta_over_pi0_40p"]);
		divide(_h["eta_60p"],_h["pi0_60p"],_s["eta_over_pi0_60p"]);


		//RAB plots

		divide(_h["CuAu_minbias_pi0"],_h["pp_pi0"],_s["R_AB_pi0"]);
		divide(_h["CuAu_minbias_eta"],_h["pp_eta"],_s["R_AB_eta"]);

      
		bool pi0_pp_available = false;
		bool pi0_CuAu_available = false;
		bool eta_pp_available = false;
		bool eta_CuAu_available = false;












		for (auto element : sow)
		{
			string name = element.second->name();
			if (name.find("pi0_minbias") != std::string::npos)
			{
				if (element.second->sumW()>0) pi0_CuAu_available = true;
				else
				{
					pi0_CuAu_available=false;
					break;
				}
			}
			else if (name.find("pp_pi0") != std::string::npos)
			{
				if (element.second->sumW()>0) pi0_pp_available=true;
				else
				{
					pi0_pp_available=false;
					break;
				}
			}
			else if (name.find("eta_minbias") != std::string::npos)
			{
				if (element.second->sumW()>0) eta_pp_available=true;
				else
				{
					eta_pp_available=false;
					break;
				}
			}
			else if (name.find("pp_eta") != std::string::npos)
			{
				if (element.second->sumW()>0) eta_pp_available=true;
				else
				{
					eta_pp_available=false;
					break;
				}
			}
		}
	if((!pi0_pp_available) || (!pi0_CuAu_available) || (!eta_pp_available) || (!eta_CuAu_available)) return;
	}



	/// @name Histograms
	///@{
	map<string, Histo1DPtr> _h;
	map<string, Profile1DPtr> _p;
	map<string, CounterPtr> sow;
	map<string, Scatter2DPtr> _s;
	string beamOpt = "";
	///@}


};

	DECLARE_RIVET_PLUGIN(PHENIX_2018_I1672859);

}
