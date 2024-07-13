// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "../Centralities/RHICCentrality.hh"
#include <math.h>
#define _USE_MATH_DEFINES
const double TAB[9] = {2.54, 8.8, 6.0, 7.5, 3.1, 1.0, 0.179, 8.9, 2.86};

//TAB/TAA  values order: CuAuMinbias, CuAu0-10, CuAu10-20, CuAu0-20, CuAu20-40, CuAu40-60, CuAu60-90, AuAu20-30, AuAu40-50. All in units of mb^-1.
//still missing CuCu0-10

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
		book(_h["pi0_6090"],  7, 1, 1);
      
		book(_h["eta_minbias"],  8, 1, 1);
		book(_h["eta_0020"],  9, 1, 1);
		book(_h["eta_2040"], 10, 1, 1);
		book(_h["eta_4060"], 11, 1, 1);
		book(_h["eta_6090"], 12, 1, 1);		//60-90 not shown in paper, but we have the data
      
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

		//RAB plots
      
		string refname18 = mkAxisCode(18, 1, 1);
		const Scatter2D& refdata18 = refData(refname18);
		book(_h["CuAu_minbias_pi0"], refname18 + "_CuAu_minbias", refdata18);
		book(_h["pp_pi0_minbias"], refname18 + "_pp_minbias", refdata18);
		book(_s["R_AB_pi0_minbias"], refname18);

		string refname19 = mkAxisCode(19, 1, 1);
		const Scatter2D& refdata19 = refData(refname19);
		book(_h["CuAu_0020_pi0"], refname19 + "_CuAu_0020", refdata19);
		book(_h["pp_pi0_0020"], refname19 + "_pp", refdata19);
		book(_s["R_AB_pi0_0020"], refname19);

		string refname20 = mkAxisCode(20, 1, 1);
		const Scatter2D& refdata20 = refData(refname20);
		book(_h["CuAu_2040_pi0"], refname20 + "_CuAu_2040", refdata20);
		book(_h["pp_pi0_2040"], refname20 + "_pp", refdata20);
		book(_s["R_AB_pi0_2040"], refname20);

		string refname21 = mkAxisCode(21, 1, 1);
		const Scatter2D& refdata21 = refData(refname21);
		book(_h["CuAu_4060_pi0"], refname21 + "_CuAu_4060", refdata21);
		book(_h["pp_pi0_4060"], refname21 + "_pp", refdata21);
		book(_s["R_AB_pi0_4060"], refname21);

		string refname22 = mkAxisCode(22, 1, 1);
		const Scatter2D& refdata22 = refData(refname22);
		book(_h["CuAu_6090_pi0"], refname22 + "_CuAu_6090", refdata22);
		book(_h["pp_pi0_6090"], refname22 + "_pp", refdata22);
		book(_s["R_AB_pi0_6090"], refname22);

		string refname23 = mkAxisCode(23, 1, 1);
		const Scatter2D& refdata23 = refData(refname23);
		book(_h["CuAu_minbias_eta"], refname23 + "_CuAu_minbias", refdata23);
		book(_h["pp_eta_minbias"], refname23 + "_pp", refdata23);
		book(_s["R_AB_eta_minbias"], refname23);

		string refname24 = mkAxisCode(24, 1, 1);
		const Scatter2D& refdata24 = refData(refname24);
		book(_h["CuAu_0020_eta"], refname24 + "_CuAu_0020", refdata24);
		book(_h["pp_eta_0020"], refname24 + "_pp", refdata24);
		book(_s["R_AB_eta_0020"], refname24);

		string refname25 = mkAxisCode(25, 1, 1);
		const Scatter2D& refdata25 = refData(refname25);
		book(_h["CuAu_2040_eta"], refname25 + "_CuAu_2040", refdata25);
		book(_h["pp_eta_2040"], refname25 + "_pp", refdata25);
		book(_s["R_AB_eta_2040"], refname25);

		string refname26 = mkAxisCode(26, 1, 1);
		const Scatter2D& refdata26 = refData(refname26);
		book(_h["CuAu_4060_eta"], refname26 + "_CuAu_4060", refdata26);
		book(_h["pp_eta_4060"], refname26 + "_pp", refdata26);
		book(_s["R_AB_eta_4060"], refname26);

		string refname27 = mkAxisCode(27, 1, 1);
		const Scatter2D& refdata27 = refData(refname27);
		book(_h["CuAu_6090_eta"], refname27 + "_CuAu_6090", refdata27);
		book(_h["pp_eta_6090"], refname27 + "_pp", refdata27);
		book(_s["R_AB_eta_6090"], refname27);

		string refname42 = mkAxisCode(42, 1, 1);
		const Scatter2D& refdata42 = refData(refname42);
		book(_h["AuAu_2030_pi0"], refname42 + "_AuAu_2030", refdata42);
		book(_h["pp_pi0_2030"], refname42 + "_pp", refdata42);
		book(_s["R_AB_pi0_2030_AuAu"], refname42);

		string refname43 = mkAxisCode(43, 1, 1);
		const Scatter2D& refdata43 = refData(refname43);
		book(_h["CuCu_0010_pi0"], refname43 + "_CuCu_0010", refdata43);
		book(_h["pp_pi0_0010"], refname43 + "_pp", refdata43);
		book(_s["R_AB_pi0_0010_CuCu"], refname43);

		string refname44 = mkAxisCode(44, 1, 1);
		const Scatter2D& refdata44 = refData(refname44);
		book(_h["AuAu_4050_pi0"], refname44 + "_AuAu_4050", refdata44);
		book(_h["pp_pi0_4050"], refname44 + "_pp", refdata44);
		book(_s["R_AB_pi0_4050_AuAu"], refname44);

		string refname45 = mkAxisCode(45, 1, 1);
		const Scatter2D& refdata45 = refData(refname45);
		book(_h["CuAu_pi0_int5"], refname45 + "_CuAu_int5", refdata45);
		book(_h["pp_pi0_int5_CuAu"], refname45 + "pp", refdata45);
		book(_s["R_AB_pi0_int5_CuAu"], refname45);

		string refname46 = mkAxisCode(46, 1, 1);
		const Scatter2D& refdata46 = refData(refname46);
		book(_h["AuAu_pi0_int5"], refname46 + "_AuAu_int5", refdata46);
		book(_h["pp_pi0_int5_AuAu"], refname46 + "pp", refdata46);
		book(_s["R_AB_pi0_int5_AuAu"], refname46);

		string refname47 = mkAxisCode(47, 1, 1);
		const Scatter2D& refdata47 = refData(refname47);
		book(_h["CuAu_pi0_int10"], refname47 + "_CuAu_int10", refdata47);
		book(_h["pp_pi0_int10_CuAu"], refname47 + "pp", refdata47);
		book(_s["R_AB_pi0_int10_CuAu"], refname47);

		string refname48 = mkAxisCode(48, 1, 1);
		const Scatter2D& refdata48 = refData(refname48);
		book(_h["AuAu_pi0_int10"], refname48 + "_AuAu_int10", refdata48);
		book(_h["pp_pi0_int10_AuAu"], refname48 + "pp", refdata48);
		book(_s["R_AB_pi0_int10_AuAu"], refname48);



		//Counters for normalization

		book(sow["Nevent_minbias_CuAu"],"_Nevent_minbias_CuAu"); 
		book(sow["Nevent_0010_CuAu"],"_Nevent_0010_CuAu"); 
		book(sow["Nevent_1020_CuAu"],"_Nevent_1020_CuAu");
		book(sow["Nevent_0020_CuAu"],"_Nevent_0020_CuAu");   
		book(sow["Nevent_2040_CuAu"],"_Nevent_2040_CuAu");
		book(sow["Nevent_4060_CuAu"],"_Nevent_4060_CuAu"); 
 		book(sow["Nevent_6090_CuAu"],"_Nevent_6090_CuAu");
		//book(sow["Nevent_int5_CuAu"],"Nevent_int5_CuAu");
		//book(sow["Nevent_int10_CuAu"],"Nevent_int10_CuAu");

		book(sow["Nevent_PP"],"_Nevent_PP");
		//book(sow["Nevent_int5_PP"],"Nevent_int5_PP");
		//book(sow["Nevent_int10_PP"],"Nevent_int10_PP");
	
		book(sow["Nevent_2030_AuAu"],"_Nevent_2030_AuAu");
		book(sow["Nevent_4050_AuAu"],"_Nevent_4050_AuAu");
		book(sow["Nevent_minbias_AuAu"],"_Nevent_minbias_AuAu");
		//book(sow["Nevent_int5_AuAu"],"Nevent_int5_AuAu");
		//book(sow["Nevent_int10_AuAu"],"Nevent_int10_AuAu");

		book(sow["Nevent_0010_CuCu"],"_Nevent_0010_CuCu");
	
	}


	/// Perform the per-event analysis
	void analyze(const Event& event)
	{

		const ParticlePair& beam = beams();
		const UnstableParticles& usp = applyProjection<UnstableParticles>(event, "UFS");
		const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
		const double c = cent();
      
		if (beamOpt == "NONE") {
		
			if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
			{
				float NN = 197.;
            	if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSys = AuAu200;
			}

			else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) collSys = pp;
			{
				float NN = 1.;
            	if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSys = pp;
			}

			else if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000290630)
			{
				float NN = 63.;
				if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSys = CuCu200;
			}

			else if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000791970)
			{
				if (fuzzyEquals(sqrtS()/GeV, 200*sqrt(63*197), 1E-3)) collSys = CuAu200;
			}
		}

		else if (beamOpt == "AUAU200") collSys = AuAu200;
		
		else if (beamOpt == "PP200") collSys = pp;

		else if (beamOpt == "CUCU200") collSys = CuCu200;
		
		else if (beamOpt == "CUAU200") collSys = CuAu200;

		for(const Particle & p : usp.particles())
		{
			if(collSys == CuAu200) 
			{
				if(c > 90.) vetoEvent;
				
				if(p.pid() == PID::PI0)
				{
            
					if(c <= 10.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["pi0_0010"]->fill(partPt, ptWeight);	
					}
					if(c > 10 && c < 20.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["pi0_1020"]->fill(partPt, ptWeight);
					}
					if(c > 0 && c < 20.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["pi0_0020"]->fill(partPt, ptWeight);
						_h["pi0_0p"]->fill(partPt, ptWeight);
						_h["CuAu_0020_pi0"]->fill(partPt);

					}
					if(c > 20 && c < 40.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["pi0_2040"]->fill(partPt, ptWeight);
						_h["pi0_20p"]->fill(partPt, ptWeight);
						_h["CuAu_2040_pi0"]->fill(partPt);
							
					}
					if(c > 40 && c < 60.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);			
						_h["pi0_4060"]->fill(partPt, ptWeight);
						_h["pi0_40p"]->fill(partPt, ptWeight);
						_h["CuAu_4060_pi0"]->fill(partPt);


					}
					if(c > 60 && c < 90.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);
						_h["pi0_6090"]->fill(partPt, ptWeight);
						_h["pi0_60p"]->fill(partPt, ptWeight);
						_h["CuAu_6090_pi0"]->fill(partPt);


					}

					double partPt = p.pT()/GeV;
					double ptWeight = 1./(partPt*2.*M_PI);
					_h["pi0_minbias"]->fill(partPt, ptWeight);
					_h["pi0_mb"]->fill(partPt, ptWeight);
					_h["CuAu_minbias_pi0"]->fill(partPt);


					if(p.pT()/GeV > 5) _h["CuAu_pi0_int5"]->fill(c);
					
					if(p.pT()/GeV > 10) _h["CuAu_pi0_int10"]->fill(c);					
				}
          
				else if(p.pid() == PID::ETA)
				{
            
					if(c > 0 && c < 20.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["eta_0020"]->fill(partPt, ptWeight);
						_h["eta_0p"]->fill(partPt, ptWeight);
						_h["CuAu_0020_eta"]->fill(partPt);
	
					}
					if(c > 20 && c < 40.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["eta_2040"]->fill(partPt, ptWeight);
						_h["eta_20p"]->fill(partPt, ptWeight);	
						_h["CuAu_2040_eta"]->fill(partPt);

					}
					if(c > 40 && c < 60.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);	
						_h["eta_4060"]->fill(partPt, ptWeight);
						_h["eta_40p"]->fill(partPt, ptWeight);
						_h["CuAu_4060_eta"]->fill(partPt);

					}
					if(c > 60 && c < 90.)
					{
						double partPt = p.pT()/GeV;
						double ptWeight = 1./(partPt*2.*M_PI);
						_h["eta_6090"]->fill(partPt, ptWeight);
						_h["eta_60p"]->fill(partPt, ptWeight);
						_h["CuAu_6090_eta"]->fill(partPt);

					}

					double partPt = p.pT()/GeV;
					double ptWeight = 1./(partPt*2.*M_PI);	
					_h["eta_minbias"]->fill(partPt, ptWeight);
					_h["eta_mb"]->fill(partPt, ptWeight);	
					_h["CuAu_minbias_eta"]->fill(partPt);
					
				}
			}        
			else if(collSys == pp)
			{          
				if(p.pid() == PID::PI0)
				{
					
					_h["pp_pi0_minbias"]->fill(p.pT()/GeV);
					
					_h["pp_pi0_0020"]->fill(p.pT()/GeV);

					_h["pp_pi0_2040"]->fill(p.pT()/GeV);

					_h["pp_pi0_4060"]->fill(p.pT()/GeV);

					_h["pp_pi0_6090"]->fill(p.pT()/GeV);

					_h["pp_pi0_0010"]->fill(p.pT()/GeV);

					_h["pp_pi0_2030"]->fill(p.pT()/GeV);

					_h["pp_pi0_4050"]->fill(p.pT()/GeV);
			
					if(p.pT()/GeV > 5) 
					{
						_h["pp_pi0_int5_CuAu"]->fill(c);
					
						_h["pp_pi0_int5_AuAu"]->fill(c);

					}
					if(p.pT()/GeV > 10)
					{
						_h["pp_pi0_int10_CuAu"]->fill(c);

						_h["pp_pi0_int10_AuAu"]->fill(c);

					}


				}
				else if(p.pid() == PID::ETA)
				{
					_h["pp_eta_minbias"]->fill(p.pT()/GeV);

					_h["pp_eta_0020"]->fill(p.pT()/GeV);
           
					_h["pp_eta_2040"]->fill(p.pT()/GeV);
 
					_h["pp_eta_4060"]->fill(p.pT()/GeV);

					_h["pp_eta_6090"]->fill(p.pT()/GeV);
					

				}
			}
			else if(collSys == AuAu200)
			{
				if(p.pid() == PID::PI0)
				{
					if(c > 20 && c < 30.) _h["AuAu_2030_pi0"]->fill(p.pT()/GeV);
					
					if(c > 40 && c < 50.) _h["AuAu_4050_pi0"]->fill(p.pT()/GeV);

					if(p.pT()/GeV > 5) _h["AuAu_pi0_int5"]->fill(c);
					
					if(p.pT()/GeV > 10) _h["AuAu_pi0_int10"]->fill(c);						

				}
			}
			else if(collSys == CuCu200)
			{
				if(p.pid() == PID::PI0)
				{
					if(c > 0 && c < 10.) _h["CuCu_0010_pi0"]->fill(p.pT()/GeV);


				}
			}

		}

		//fill counters outside of the particle loop

		if(collSys == CuAu200)
		{
	//		if(p.pT()/GeV > 5) sow["Nevent_int5_CuAu"]->fill();

	//		if(p.pT()/GeV > 10) sow["Nevent_int10_CuAu"]->fill();
	
			sow["Nevent_minbias_CuAu"]->fill();

			if(c > 0 && c < 10.) sow["Nevent_0010_CuAu"]->fill();
			
			if(c > 10 && c < 20.) sow["Nevent_1020_CuAu"]->fill();

			if(c > 0 && c < 20.) sow["Nevent_0020_CuAu"]->fill();

			if(c > 20 && c < 40.) sow["Nevent_2040_CuAu"]->fill();

			if(c > 40 && c < 60.) sow["Nevent_4060_CuAu"]->fill();

			if(c > 60 && c < 90.) sow["Nevent_6090_CuAu"]->fill();

		}
		else if(collSys == pp)
		{
			sow["Nevent_PP"]->fill();
				
	//		if(p.pT()/GeV > 5) sow["Nevent_int5_PP"]->fill();

	//		if(p.pT()/GeV > 10) sow["Nevent_int10_PP"]->fill();
		}
		else if(collSys == AuAu200)
		{
	//		if(p.pT()/GeV > 5) sow["Nevent_int5_AuAu"]->fill();

	//		if(p.pT()/GeV > 10) sow["Nevent_int10_AuAu"]->fill();
	//
			sow["Nevent_minbias_AuAu"]->fill();

			if(c > 20 && c < 30.) sow["Nevent_2030_AuAu"]->fill();	
      			
			if(c > 40 && c < 50.) sow["Nevent_4050_AuAu"]->fill();
		}
		else if(collSys == CuCu200)
		{
			if(c > 0 && c < 10.) sow["Nevent_0010_CuCu"]->fill();
		}
	}


	/// Normalise histograms etc., after the run
	void finalize()
	{

		//dividing yields by # of events N

			//pi 0 normalizations			
		
		_h["pi0_0010"]->scaleW(1./(sow["Nevent_0010_CuAu"]->sumW()));	
		_h["pi0_1020"]->scaleW(1./(sow["Nevent_1020_CuAu"]->sumW()));
		_h["pi0_0020"]->scaleW(1./(sow["Nevent_0020_CuAu"]->sumW()));	
		_h["pi0_2040"]->scaleW(1./(sow["Nevent_2040_CuAu"]->sumW()));	
		_h["pi0_4060"]->scaleW(1./(sow["Nevent_4060_CuAu"]->sumW()));
		_h["pi0_6090"]->scaleW(1./(sow["Nevent_6090_CuAu"]->sumW()));	
		_h["pi0_minbias"]->scaleW(1./(sow["Nevent_minbias_CuAu"]->sumW()));
	
		_h["pi0_mb"]->scaleW(1./(sow["Nevent_minbias_CuAu"]->sumW()));
		_h["pi0_0p"]->scaleW(1./(sow["Nevent_0020_CuAu"]->sumW()));
		_h["pi0_20p"]->scaleW(1./(sow["Nevent_2040_CuAu"]->sumW()));
		_h["pi0_40p"]->scaleW(1./(sow["Nevent_4060_CuAu"]->sumW()));
		_h["pi0_60p"]->scaleW(1./(sow["Nevent_6090_CuAu"]->sumW()));
			
		_h["CuAu_minbias_pi0"]->scaleW(1./(sow["Nevent_minbias_CuAu"]->sumW()));	
		_h["CuAu_0020_pi0"]->scaleW(1./(sow["Nevent_0020_CuAu"]->sumW()));
		_h["CuAu_2040_pi0"]->scaleW(1./(sow["Nevent_2040_CuAu"]->sumW()));
		_h["CuAu_4060_pi0"]->scaleW(1./(sow["Nevent_4060_CuAu"]->sumW()));
		_h["CuAu_6090_pi0"]->scaleW(1./(sow["Nevent_6090_CuAu"]->sumW()));
		
		_h["AuAu_2030_pi0"]->scaleW(1./(sow["Nevent_2030_AuAu"]->sumW()));
		_h["AuAu_4050_pi0"]->scaleW(1./(sow["Nevent_4050_AuAu"]->sumW()));
		_h["CuCu_0010_pi0"]->scaleW(1./(sow["Nevent_0010_CuCu"]->sumW()));

		_h["CuAu_pi0_int5"]->scaleW(1./(sow["Nevent_minbias_CuAu"]->sumW()));
		_h["CuAu_pi0_int10"]->scaleW(1./(sow["Nevent_minbias_CuAu"]->sumW()));
		_h["AuAu_pi0_int5"]->scaleW(1./(sow["Nevent_minbias_AuAu"]->sumW()));
		_h["AuAu_pi0_int10"]->scaleW(1./(sow["Nevent_minbias_AuAu"]->sumW()));	

			//eta normalizations

		_h["eta_0020"]->scaleW(1./(sow["Nevent_0020_CuAu"]->sumW()));		
		_h["eta_2040"]->scaleW(1./(sow["Nevent_2040_CuAu"]->sumW()));	
		_h["eta_4060"]->scaleW(1./(sow["Nevent_4060_CuAu"]->sumW()));
		_h["eta_6090"]->scaleW(1./(sow["Nevent_6090_CuAu"]->sumW()));		
		_h["eta_minbias"]->scaleW(1./(sow["Nevent_minbias_CuAu"]->sumW()));
	
		_h["eta_mb"]->scaleW(1./(sow["Nevent_minbias_CuAu"]->sumW()));
		_h["eta_0p"]->scaleW(1./(sow["Nevent_0020_CuAu"]->sumW()));
		_h["eta_20p"]->scaleW(1./(sow["Nevent_2040_CuAu"]->sumW()));
		_h["eta_40p"]->scaleW(1./(sow["Nevent_4060_CuAu"]->sumW()));
		_h["eta_60p"]->scaleW(1./(sow["Nevent_6090_CuAu"]->sumW()));
			
		_h["CuAu_minbias_eta"]->scaleW(1./(sow["Nevent_minbias_CuAu"]->sumW()));	
		_h["CuAu_0020_eta"]->scaleW(1./(sow["Nevent_0020_CuAu"]->sumW()));
		_h["CuAu_2040_eta"]->scaleW(1./(sow["Nevent_2040_CuAu"]->sumW()));
		_h["CuAu_4060_eta"]->scaleW(1./(sow["Nevent_4060_CuAu"]->sumW()));
		_h["CuAu_6090_eta"]->scaleW(1./(sow["Nevent_6090_CuAu"]->sumW()));

			//pp normalizations

		_h["pp_pi0_minbias"]->scaleW(1./(sow["Nevent_PP"]->sumW()));
		_h["pp_pi0_0020"]->scaleW(1./(sow["Nevent_PP"]->sumW()));			
		_h["pp_pi0_2040"]->scaleW(1./(sow["Nevent_PP"]->sumW()));	
		_h["pp_pi0_4060"]->scaleW(1./(sow["Nevent_PP"]->sumW()));	
		_h["pp_pi0_6090"]->scaleW(1./(sow["Nevent_PP"]->sumW()));

		_h["pp_eta_minbias"]->scaleW(1./(sow["Nevent_PP"]->sumW()));	
		_h["pp_eta_0020"]->scaleW(1./(sow["Nevent_PP"]->sumW()));		
		_h["pp_eta_2040"]->scaleW(1./(sow["Nevent_PP"]->sumW()));	
		_h["pp_eta_4060"]->scaleW(1./(sow["Nevent_PP"]->sumW()));	
		_h["pp_eta_6090"]->scaleW(1./(sow["Nevent_PP"]->sumW()));

		_h["pp_pi0_0010"]->scaleW(1./(sow["Nevent_PP"]->sumW()));	
		_h["pp_pi0_2030"]->scaleW(1./(sow["Nevent_PP"]->sumW()));	
		_h["pp_pi0_4050"]->scaleW(1./(sow["Nevent_PP"]->sumW()));

		_h["pp_pi0_int5_CuAu"]->scaleW(1./(sow["Nevent_PP"]->sumW()));
		_h["pp_pi0_int10_CuAu"]->scaleW(1./(sow["Nevent_PP"]->sumW()));
		_h["pp_pi0_int5_AuAu"]->scaleW(1./(sow["Nevent_PP"]->sumW()));
		_h["pp_pi0_int10_AuAu"]->scaleW(1./(sow["Nevent_PP"]->sumW()));	

		//ratio plots
      
		divide(_h["eta_mb"],_h["pi0_mb"],_s["eta_over_pi0_mb"]);
		divide(_h["eta_0p"],_h["pi0_0p"],_s["eta_over_pi0_0p"]);
		divide(_h["eta_20p"],_h["pi0_20p"],_s["eta_over_pi0_20p"]);
		divide(_h["eta_40p"],_h["pi0_40p"],_s["eta_over_pi0_40p"]);
		divide(_h["eta_60p"],_h["pi0_60p"],_s["eta_over_pi0_60p"]);

		//RAB plots

		divide(_h["CuAu_minbias_pi0"],_h["pp_pi0_minbias"],_s["R_AB_pi0_minbias"]);
		divide(_h["CuAu_0020_pi0"],_h["pp_pi0_0020"],_s["R_AB_pi0_0020"]);
		divide(_h["CuAu_2040_pi0"],_h["pp_pi0_2040"],_s["R_AB_pi0_2040"]);
		divide(_h["CuAu_4060_pi0"],_h["pp_pi0_4060"],_s["R_AB_pi0_4060"]);
		divide(_h["CuAu_6090_pi0"],_h["pp_pi0_6090"],_s["R_AB_pi0_6090"]);

		divide(_h["CuAu_minbias_eta"],_h["pp_eta_minbias"],_s["R_AB_eta_minbias"]);
		divide(_h["CuAu_0020_eta"],_h["pp_eta_0020"],_s["R_AB_eta_0020"]);
		divide(_h["CuAu_2040_eta"],_h["pp_eta_2040"],_s["R_AB_eta_2040"]);
		divide(_h["CuAu_4060_eta"],_h["pp_eta_4060"],_s["R_AB_eta_4060"]);
		divide(_h["CuAu_6090_eta"],_h["pp_eta_6090"],_s["R_AB_eta_6090"]);

		divide(_h["CuCu_0010_pi0"],_h["pp_pi0_0010"],_s["R_AB_pi0_0010_CuCu"]);
		divide(_h["AuAu_2030_pi0"],_h["pp_pi0_2030"],_s["R_AB_pi0_2030_AuAu"]);
		divide(_h["AuAu_4050_pi0"],_h["pp_pi0_4050"],_s["R_AB_pi0_4050_AuAu"]);

		divide(_h["CuAu_pi0_int5"],_h["pp_pi0_int5_CuAu"],_s["R_AB_pi0_int5_CuAu"]);
		divide(_h["CuAu_pi0_int10"],_h["pp_pi0_int10_CuAu"],_s["R_AB_pi0_int10_CuAu"]);
		divide(_h["AuAu_pi0_int5"],_h["pp_pi0_int5_AuAu"],_s["R_AB_pi0_int5_AuAu"]);
		divide(_h["AuAu_pi0_int10"],_h["pp_pi0_int10_AuAu"],_s["R_AB_pi0_int10_AuAu"]);


		//RAB normalization by TAb Nuclear overlap 

		_s["R_AB_pi0_minbias"]->scaleY(1./TAB[0]);
		_s["R_AB_pi0_0020"]->scaleY(1./TAB[3]);
		_s["R_AB_pi0_2040"]->scaleY(1./TAB[4]);
		_s["R_AB_pi0_4060"]->scaleY(1./TAB[5]);
		_s["R_AB_pi0_6090"]->scaleY(1./TAB[6]);


		_s["R_AB_eta_minbias"]->scaleY(1./TAB[0]);
		_s["R_AB_eta_0020"]->scaleY(1./TAB[3]);
		_s["R_AB_eta_2040"]->scaleY(1./TAB[4]);
		_s["R_AB_eta_4060"]->scaleY(1./TAB[5]);
		_s["R_AB_eta_4060"]->scaleY(1./TAB[6]);

		_s["R_AB_pi0_2030_AuAu"]->scaleY(1./TAB[7]);
		_s["R_AB_pi0_4050_AuAu"]->scaleY(1./TAB[8]);

/*

		normalize(_h["pi0_minbias"]); // normalize to unity		//not sure if this is needed rn
		normalize(_h["pi0_0010"]); // normalize to unity
		normalize(_h["pi0_1020"]); // normalize to unity
		normalize(_h["pi0_0020"]); // normalize to unity
		normalize(_h["pi0_2040"]); // normalize to unity
		normalize(_h["pi0_4060"]); // normalize to unity
		normalize(_h["pi0_6090"]);	
      
		normalize(_h["eta_minbias"]); // normalize to unity
		normalize(_h["eta_0020"]); // normalize to unity
		normalize(_h["eta_2040"]); // normalize to unity
		normalize(_h["eta_4060"]); // normalize to unity
		normalize(_h["eta_6090"]);

*/


		//checking beam type used and which histos will be empty
      
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
	string beamOpt = "NONE";
	enum CollisionSystem {pp, AuAu200, CuCu200, CuAu200};
	CollisionSystem collSys;
	///@}


};

	DECLARE_RIVET_PLUGIN(PHENIX_2018_I1672859);

}
