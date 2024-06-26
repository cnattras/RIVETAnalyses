// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
//#include <math.h>
//#define _USE_MATH_DEFINES

namespace Rivet {


  class PHENIX_2008_I776624 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2008_I776624);

    /////////////////////////////////////////////////////////////////////////
    //create binShift function
    void binShift(YODA::Histo1D& histogram) {
    std::vector<YODA::HistoBin1D> binlist = histogram.bins();
    int n = 0;
    for (YODA::HistoBin1D bins : binlist) {
        double p_high = bins.xMax();
        double p_low = bins.xMin();
        //Now calculate f_corr
        if (bins.xMin() == binlist[0].xMin()) { //Check if we are working with first bin
            float b = 1 / (p_high - p_low) * log(binlist[0].height()/binlist[1].height());
            float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
            histogram.bin(n).scaleW(f_corr);
            n += 1;
        } else if (bins.xMin() == binlist.back().xMin()){ //Check if we are working with last bin
            float b = 1 / (p_high - p_low) * log(binlist[binlist.size()-2].height() / binlist.back().height());
            float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
            histogram.bin(n).scaleW(f_corr);
        } else { //Check if we are working with any middle bin
            float b = 1 / (p_high - p_low) * log(binlist[n-1].height() / binlist[n+1].height());
            float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
            histogram.bin(n).scaleW(f_corr);
            n += 1;
        }
    }
}

    void init() {

      const UnstableParticles up_fwd(Cuts::absrap < 2.2 && Cuts::absrap > 1.2 && Cuts::abspid == 443 && Cuts::abscharge == 0);  // fwd rap
      declare(up_fwd,"up_fwd");
      const UnstableParticles up_mid(Cuts::absrap < 0.35 && Cuts::abspid == 443 && Cuts::abscharge == 0); // mid rap
      declare(up_mid,"up_mid");
      const UnstableParticles up_all(Cuts::abspid == 443 && Cuts::abscharge == 0); // no rapiidty cuts here
      declare(up_mid,"up_all");

      beamOpt = getOption<string>("beam","NONE");
      const ParticlePair& beam = beams();

      if (beamOpt == "NONE") {    
      if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000290630) CollSys = CuCu200;
      else if (beam.first.pid() == 2212 && beam.second.pid() == 2212) CollSys = pp200;
      }
      if(beamOpt=="PP200") CollSys = pp200;
      else if(beamOpt=="CUCU200") CollSys = CuCu200;
     
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      //declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      /////////////////////////////////////////////////////////////
     //  FIGURE 1 (Tables 1-8) - no pp histograms needed
      book(_h_1D["YAA_pT_mid_020"], 1, 1, 1);  
      book(_h_1D["YAA_pT_mid_2040"], 2, 1, 1);  
      book(_h_1D["YAA_pT_mid_4060"], 3, 1, 1); 
      book(_h_1D["YAA_pT_mid_6094"], 4, 1, 1);  // up to 10 GeV at mid rapidity
     
      book(_h_1D["YAA_pT_fwd_020"], 5, 1, 1); 
      book(_h_1D["YAA_pT_fwd_2040"], 6, 1, 1); 
      book(_h_1D["YAA_pT_fwd_4060"], 7, 1, 1); 
      book(_h_1D["YAA_pT_fwd_6094"], 8, 1, 1);   // up to 6 GeV at fwd

     // counter histograms
      book(_c["c_YAA_pT_mid_020"], "_c_YAA_pT_mid_020");  
      book(_c["c_YAA_pT_mid_2040"], "_c_YAA_pT_mid_2040");  
      book(_c["c_YAA_pT_mid_4060"], "_c_YAA_pT_mid_4060");  
      book(_c["c_YAA_pT_mid_6094"], "_c_YAA_pT_mid_6094");  
     
      book(_c["c_YAA_pT_fwd_020"], "_c_YAA_pT_fwd_020");
      book(_c["c_YAA_pT_fwd_2040"], "_c_YAA_pT_fwd_2040");
      book(_c["c_YAA_pT_fwd_4060"], "_c_YAA_pT_fwd_4060");
      book(_c["c_YAA_pT_fwd_6094"], "_c_YAA_pT_fwd_6094");

      book(_c["c_pp"], "_c_pp");
   
     //  /////////////////////////////////////////////////////////////
   //    // FIGURE 2
   //    These are functions of N_{part}
   //    // x-axis for CuCu ptsq Table 9 as function of Npart - up to 5 GeV
   //    vector<double> centBins3{0.0, 20.0, 40.0, 60.0};   // no 60-94 provided for mid rapidity
   //    book(_h1D_npart3["ptsq_mid_cent"], "ptsq_mid_cent", centBins3);  
   //    // x-axis for CuCu ptsq Table 10
   //    vector<double> centBins4{0.0, 20.0, 40.0, 60.0, 94.0};
   //    book(_h1D_npart4["ptsq_fwd_cent"], "ptsq_fwd_cent", centBins4);
   
   // // 2D scatter connected to the data
   //    book( _h2D_ptsq_npart_mid["ptsq_npart_mid"], 9,1,1);  
   //    book( _h2D_ptsq_npart_fwd["ptsq_npart_fwd"], 10,1,1);  
   //    /////////////////////////////////////////////////////////////////////////////////////////

      // FIGURE 3 - Only 0-20% Events shown in figures
    
      string refnameRaa1 = mkAxisCode(11,1,1);
      const Scatter2D& refdataRaa1 =refData(refnameRaa1);
      book(_h_RAA_1D["020_pT_mid_CuCu"], refnameRaa1 + "_CuCu200", refdataRaa1);
      book(_h_RAA_1D["pT_mid_pp"], refnameRaa1 + "_pp200", refdataRaa1);
      book(_h2D_RAA["RAA_pT_mid_020"], refnameRaa1);

      string refnameRaa2 = mkAxisCode(12,1,1);
      const Scatter2D& refdataRaa2 =refData(refnameRaa2);
      book(_h_RAA_1D["020_pT_fwd_CuCu"], refnameRaa2 + "_CuCu200", refdataRaa2);
      book(_h_RAA_1D["pT_fwd_pp"], refnameRaa2 + "_pp200", refdataRaa2);
      book(_h2D_RAA["RAA_pT_fwd_020"], refnameRaa2);

      centBins1.insert(pair<string, int>("RAA_pT_mid_020",0));
      centBins1.insert(pair<string, int>("RAA_pT_fwd_020",1));
 
      book(_c["c_CuCu_rap_all_020"], "_c_CuCu_rap_all_020");

      string refnameRaa3 = mkAxisCode(13,1,1);
      const Scatter2D& refdataRaa3 =refData(refnameRaa3);
      book(_h_RAA_1D_rap["020_rap_all_CuCu"], refnameRaa3 + "_CuCu200", refdataRaa3);
      book(_h_RAA_1D_rap["rap_all_pp"], refnameRaa3 + "_pp200", refdataRaa3);
      book(_h2D_RAA_rap["RAA_rap_all_020"], refnameRaa3);

      centBins2.insert(pair<string, int>("RAA_rap_all_020",0));

      // 	/
	////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////

      // FIGURE 4
   
      // book(_c["c_CuCu_010"], "c_CuCu_010");
      // book(_c["c_CuCu_1020"], "c_CuCu_1020");
      // book(_c["c_CuCu_2030"], "c_CuCu_2030");
      // book(_c["c_CuCu_3040"], "c_CuCu_3040");
      // book(_c["c_CuCu_4050"], "c_CuCu_4050");
      // book(_c["c_CuCu_5060"], "c_CuCu_5060");
      // book(_c["c_CuCu_6094"], "c_CuCu_6094");
      // book(_c["c_CuCu_6070"], "c_CuCu_6070");
      // book(_c["c_CuCu_7094"], "c_CuCu_7094");


      // vector<double> pTBins{0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0};  
      // book(_h_RAA_1D_cent_pT["pT_cent_mid_pp"], "pT_cent_mid_pp", pTBins);
      // book(_h_RAA_1D_cent_pT["pT_cent_fwd_pp"], "pT_cent_fwd_pp", pTBins);

      // book(_h_RAA_1D_cent_pT["010_pT_mid_CuCu"], "010_pT_mid_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["1020_pT_mid_CuCu"], "1020_pT_mid_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["2030_pT_mid_CuCu"], "2030_pT_mid_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["3040_pT_mid_CuCu"], "3040_pT_mid_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["4050_pT_mid_CuCu"], "4050_pT_mid_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["5060_pT_mid_CuCu"], "5060_pT_mid_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["6094_pT_mid_CuCu"], "6094_pT_mid_CuCu", pTBins);

      // book(_h_RAA_1D_cent_pT["010_pT_fwd_CuCu"], "010_pT_fwd_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["1020_pT_fwd_CuCu"], "1020_pT_fwd_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["2030_pT_fwd_CuCu"], "2030_pT_fwd_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["3040_pT_fwd_CuCu"], "3040_pT_fwd_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["4050_pT_fwd_CuCu"], "4050_pT_fwd_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["5060_pT_fwd_CuCu"], "5060_pT_fwd_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["6070_pT_fwd_CuCu"], "6070_pT_fwd_CuCu", pTBins);
      // book(_h_RAA_1D_cent_pT["7094_pT_fwd_CuCu"], "7094_pT_fwd_CuCu", pTBins);

      // // for RAA vs pT per centrality slice, before summing over pT
      // book(_h2D_RAA_cent_pT["RAA_pT_mid_010"], "RAA_pT_mid_010");
      // book(_h2D_RAA_cent_pT["RAA_pT_mid_1020"], "RAA_pT_mid_1020");
      // book(_h2D_RAA_cent_pT["RAA_pT_mid_2030"], "RAA_pT_mid_2030");
      // book(_h2D_RAA_cent_pT["RAA_pT_mid_3040"], "RAA_pT_mid_3040");
      // book(_h2D_RAA_cent_pT["RAA_pT_mid_4050"], "RAA_pT_mid_4050");
      // book(_h2D_RAA_cent_pT["RAA_pT_mid_5060"], "RAA_pT_mid_5060");
      // book(_h2D_RAA_cent_pT["RAA_pT_mid_6094"], "RAA_pT_mid_6094");
    
      // book(_h2D_RAA_cent_pT["RAA_pT_fwd_010"], "RAA_pT_fwd_010");
      // book(_h2D_RAA_cent_pT["RAA_pT_fwd_1020"], "RAA_pT_fwd_1020");
      // book(_h2D_RAA_cent_pT["RAA_pT_fwd_2030"], "RAA_pT_fwd_2030");
      // book(_h2D_RAA_cent_pT["RAA_pT_fwd_3040"], "RAA_pT_fwd_3040");
      // book(_h2D_RAA_cent_pT["RAA_pT_fwd_4050"], "RAA_pT_fwd_4050");
      // book(_h2D_RAA_cent_pT["RAA_pT_fwd_5060"], "RAA_pT_fwd_5060");
      // book(_h2D_RAA_cent_pT["RAA_pT_fwd_6070"], "RAA_pT_fwd_6070");
      // book(_h2D_RAA_cent_pT["RAA_pT_fwd_7094"], "RAA_pT_fwd_7094");

      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////

      // for pT integrated results using SetPoint
      // string refnameRaa4 = mkAxisCode(14,1,1);
      // const Scatter2D& refdataRaa4 =refData(refnameRaa4);
      // book(_h_RAA_1D_cent["Cent_mid_CuCu"], refnameRaa4 + "_CuCu200", refdataRaa4);
      // book(_h_RAA_1D_cent["mid_pp"], refnameRaa4 + "_pp200", refdataRaa4);
      // book(_h2D_RAA_mid_Npart["RAA_mid_Npart"], refnameRaa4);

      // string refnameRaa5 = mkAxisCode(15,1,1);
      // const Scatter2D& refdataRaa5 =refData(refnameRaa5);
      // book(_h_RAA_1D_cent["Cent_fwd_CuCu"], refnameRaa5 + "_CuCu200", refdataRaa5);
      // book(_h_RAA_1D_cent["fwd_pp"], refnameRaa5 + "_pp200", refdataRaa5);
      // book(_h2D_RAA_fwd_Npart["RAA_fwd_Npart"], refnameRaa5);

      // centBins3.insert(pair<string, int>("RAA_mid_Npart",0));
      // centBins4.insert(pair<string, int>("RAA_fwd_Npart",0));

 
    }

    /////////////////////////////////////////////////////////////////////////
    void analyze(const Event& event) {

     
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();
    
      Particles FwdParticles = applyProjection<UnstableParticles>(event,"up_fwd").particles();  // muons
      Particles MidParticles = applyProjection<UnstableParticles>(event,"up_mid").particles();  // electrons
      Particles AllParticles = applyProjection<UnstableParticles>(event,"up_all").particles();  // electrons

      if(CollSys==CuCu200)
      	{
      	  // CuCu event counters 
      	  /////////////////////////////////////////////////////////
      	  if((c > 0.) && (c <= 20.))
      	    {
	      _c["c_YAA_pT_fwd_020"]->fill();
	      _c["c_YAA_pT_mid_020"]->fill();
	      _c["c_CuCu_rap_all_020"]->fill();
      	    }
	  else if((c > 20.) && (c <= 40.))
      	    {
	      _c["c_YAA_pT_fwd_2040"]->fill();
	      _c["c_YAA_pT_mid_2040"]->fill();
      	    }
	  else if((c > 40.) && (c <= 60.))
	    {
	      _c["c_YAA_pT_fwd_4060"]->fill();
	      _c["c_YAA_pT_mid_4060"]->fill();
	    }
	  else if((c > 60.) && (c <= 94.))
	    {
	      _c["c_YAA_pT_fwd_6094"]->fill();
	      _c["c_YAA_pT_mid_6094"]->fill();
	    }
      	  if(c > 94.) vetoEvent;
      	  /////////////////////////////////////////////////////////
	  
	  // if((c > 0.) && (c <= 10.))
	  //   _c["c_CuCu_010"]->fill();
	  // if((c > 10.) && (c <= 20.))
	  //   _c["c_CuCu_1020"]->fill();
	  // if((c > 20.) && (c <= 30.))
	  //   _c["c_CuCu_2030"]->fill();
	  // if((c > 30.) && (c <= 40.))
	  //   _c["c_CuCu_3040"]->fill();
	  // if((c > 40.) && (c <= 50.))
	  //   _c["c_CuCu_4050"]->fill();
	  // if((c > 50.) && (c <= 60.))
	  //   _c["c_CuCu_5060"]->fill();
	  // if((c > 60.) && (c <= 70.))
	  //   _c["c_CuCu_6070"]->fill();
	  // if((c > 60.) && (c <= 94.))
	  //   _c["c_CuCu_6094"]->fill();
	  // if((c > 70.) && (c <= 94.))
	  //   _c["c_CuCu_7094"]->fill();
	
  
      	  // // CuCu fwd J/psi
      	  // /////////////////////////////////////////////////////////
	  for(const Particle& p : FwdParticles)
      	    {

	      double jpsi_pT = p.pT()/GeV;
	      double pt_weight = 1./(jpsi_pT*2.*M_PI);
	    
      	      if(p.pid() == 443)
      	  	{
      	  	  if( (c > 0.) && (c <= 20.))
      	  	    {
      	  	      if(jpsi_pT < 6.)
			{
			  _h_1D["YAA_pT_fwd_020"]->fill(jpsi_pT, pt_weight); // Fig 1
			  _h_RAA_1D["020_pT_fwd_CuCu"]->fill(jpsi_pT);  // Fig 3
			}
		    }
      	  	  else if((c > 20.) && (c <= 40.))
      	  	    {
      	  	      if(jpsi_pT < 6.)
      	  		_h_1D["YAA_pT_fwd_2040"]->fill(jpsi_pT, pt_weight);
		    }
      	  	  else if((c > 40.) && (c <= 60.))
      	  	    {
      	  	      if(jpsi_pT < 6.)
      	  		_h_1D["YAA_pT_fwd_4060"]->fill(jpsi_pT, pt_weight);
		    }
      	  	  else if((c > 60.) && (c <= 94.))
      	  	    {
      	  	      if(jpsi_pT < 6.)
      	  		_h_1D["YAA_pT_fwd_6094"]->fill(jpsi_pT, pt_weight);
      	  	    }
      	  	  else if(c > 94.) vetoEvent;
      	  	} // 443
      	    } // for

      	  // // CuCu mid J/psi
      	  // /////////////////////////////////////////////////////////
      	  for(const Particle& p : MidParticles)
      	    {
      	      if(p.pid() == 443)
      	  	{

		  double jpsi_pT = p.pT()/GeV;
		  double pt_weight = 1./(jpsi_pT*2.*M_PI);

      	  	  if( (c > 0.) && (c <= 20.))
      	  	    {
      	  	      if(jpsi_pT < 10.)
			{
			  _h_1D["YAA_pT_mid_020"]->fill(jpsi_pT, pt_weight); // Fig 1 
			  _h_RAA_1D["020_pT_mid_CuCu"]->fill(jpsi_pT); // Fig 3
			}
		    }
      	  	  else if((c > 20.) && (c <= 40.))
      	  	    {
      	  	      if(jpsi_pT < 10.)
      	  		_h_1D["YAA_pT_mid_2040"]->fill(jpsi_pT, pt_weight);
		    }
      	  	  else if((c >= 40.) && (c <= 60.))
      	  	    {
      	  	      if(jpsi_pT < 10.)
      	  		_h_1D["YAA_pT_mid_4060"]->fill(jpsi_pT, pt_weight);
		    }
      	  	  else if((c >= 60.) && (c <= 94.))
      	  	    {
      	  	      if(jpsi_pT < 10.)
      	  		_h_1D["YAA_pT_mid_6094"]->fill(jpsi_pT, pt_weight);
		      if(c > 94.) vetoEvent;
		    }
		} // 443
	    } // for
	 
	  // for pCu RAA vs rapidity, Figure 3
      	  /////////////////////////////////////////////////////////
      	  for(const Particle& p :AllParticles)
      	    {
      	      if(p.pid() == 443)
      		{
		  double jpsi_rap = p.rapidity();
		  
		  if( (c > 0.) && (c <= 20.))
		    _h_RAA_1D_rap["020_rap_all_CuCu"]->fill(jpsi_rap);  // Fig 3

      		}
      	    }

	  // for(const Particle& p : MidParticles)
      	  //   {
	  //     double jpsi_pT = p.pT()/GeV;
	     	    
      	  //     if(p.pid() == 443)
      	  // 	{
	  // 	  if(jpsi_pT < 5.)
	  // 	    {
	  // 	      if( (c > 0.) && (c <= 10.))
	  // 		_h_RAA_1D_cent_pT["010_pT_mid_CuCu"]->fill(jpsi_pT);	 // Fig 4     	  	    	  	    	      	  	    
	  // 	      else if((c > 10.) && (c <= 20.))
	  // 		_h_RAA_1D_cent_pT["1020_pT_mid_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if((c > 20.) && (c <= 30.))
	  // 	      	_h_RAA_1D_cent_pT["2030_pT_mid_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if((c > 30.) && (c <= 40.))
	  // 	      	_h_RAA_1D_cent_pT["3040_pT_mid_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if((c > 40.) && (c <= 50.))
	  // 	      	_h_RAA_1D_cent_pT["4050_pT_mid_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if((c > 50.) && (c <= 60.))
	  // 	      	_h_RAA_1D_cent_pT["5060_pT_mid_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if((c > 60.) && (c <= 94.))
	  // 	      	_h_RAA_1D_cent_pT["6094_pT_mid_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if(c > 94.) vetoEvent;
	  // 	    }
      	  // 	}
	  //   }
	  
	  // for(const Particle& p : FwdParticles)
	  //   {
	  //     double jpsi_pT = p.pT()/GeV;
	    		  
	  //     if(p.pid() == 443)
	  // 	{
	  // 	  if(jpsi_pT < 5.)
	  // 	    {
	  // 	      if( (c > 0.) && (c <= 10.))
	  // 		_h_RAA_1D_cent_pT["010_pT_fwd_CuCu"]->fill(jpsi_pT);	  // Fig 4    	  	    	  	    	      	  	    
	  // 	      else if((c > 10.) && (c <= 20.))
	  // 		_h_RAA_1D_cent_pT["1020_pT_fwd_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if((c > 20.) && (c <= 30.))
	  // 	      	_h_RAA_1D_cent_pT["2030_pT_fwd_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if((c > 30.) && (c <= 40.))
	  // 	      	_h_RAA_1D_cent_pT["3040_pT_fwd_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if((c > 40.) && (c <= 50.))
	  // 	      	_h_RAA_1D_cent_pT["4050_pT_fwd_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if((c > 50.) && (c <= 60.))
	  // 	      	_h_RAA_1D_cent_pT["5060_pT_fwd_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if((c > 60.) && (c <= 70.))
	  // 	      	_h_RAA_1D_cent_pT["6070_pT_fwd_CuCu"]->fill(jpsi_pT);
	  // 	      else if((c > 70.) && (c <= 94.))
	  // 	      	_h_RAA_1D_cent_pT["7094_pT_fwd_CuCu"]->fill(jpsi_pT);	      	  	    	  	    
	  // 	      else if(c > 94.) vetoEvent;
	  // 	    }
	  // 	}
	  //   } // for

	} // if CuCu

      //////////////////////////////////////////////////////////////////////////////
       if(CollSys == pp200)
	 //if(CollSys == CuCu200)
	//////////////////////////////////////////////////////////////////////////////
      	{

	  _c["c_pp"]->fill();
	

      	  // for pp RAA vs rapidity, Figure 3
      	  /////////////////////////////////////////////////////////
      	  for(const Particle& p :AllParticles)
      	    {
	   
      	      if(p.pid() == 443)
      		{
		  double jpsi_rap = p.rapidity();
		  
		  _h_RAA_1D_rap["rap_all_pp"]->fill(jpsi_rap);  //Fig 3

      		}
      	    }
      	  // for pp RAA fwd
      	  /////////////////////////////////////////////////////////
      	  for(const Particle& p : FwdParticles)
      	    {

	      double jpsi_pT = p.pT()/GeV;
	    
      	      if((p.pid() == 443) && (jpsi_pT < 5.))
      		{
		  _h_RAA_1D["pT_fwd_pp"]->fill(jpsi_pT);  // Fig 3
		  _h_RAA_1D_cent_pT["pT_cent_fwd_pp"]->fill(jpsi_pT);	
		}
      	    }
      	  // for pp RAA mid
      	  /////////////////////////////////////////////////////////
      	  for(const Particle& p : MidParticles)
      	    {

	      double jpsi_pT = p.pT()/GeV;
	     
	      if((p.pid() == 443) && (jpsi_pT < 10.))
      		{
		  _h_RAA_1D["pT_mid_pp"]->fill(jpsi_pT);  // Fig 3
		  _h_RAA_1D_cent_pT["pT_cent_mid_pp"]->fill(jpsi_pT);	
		}
      	    }

      	} // if pp


    } // end analyze

    /////////////////////////////////////////////////////////////////////////
    void finalize() {


      // Figure 1  - yields vs pT
      binShift(*_h_1D["YAA_pT_mid_020"]);
      binShift(*_h_1D["YAA_pT_mid_2040"]);
      binShift(*_h_1D["YAA_pT_mid_4060"]);
      binShift(*_h_1D["YAA_pT_mid_6094"]);
      binShift(*_h_1D["YAA_pT_fwd_020"]);
      binShift(*_h_1D["YAA_pT_fwd_2040"]);
      binShift(*_h_1D["YAA_pT_fwd_4060"]);
      binShift(*_h_1D["YAA_pT_fwd_6094"]);
      _h_1D["YAA_pT_mid_020"]->scaleW(1./_c["c_YAA_pT_mid_020"]->sumW());
      _h_1D["YAA_pT_mid_2040"]->scaleW(1./_c["c_YAA_pT_mid_2040"]->sumW());
      _h_1D["YAA_pT_mid_4060"]->scaleW(1./_c["c_YAA_pT_mid_4060"]->sumW());
      _h_1D["YAA_pT_mid_6094"]->scaleW(1./_c["c_YAA_pT_mid_6094"]->sumW());
      _h_1D["YAA_pT_fwd_020"]->scaleW(1./_c["c_YAA_pT_fwd_020"]->sumW());
      _h_1D["YAA_pT_fwd_2040"]->scaleW(1./_c["c_YAA_pT_fwd_2040"]->sumW());
      _h_1D["YAA_pT_fwd_4060"]->scaleW(1./_c["c_YAA_pT_fwd_4060"]->sumW());
      _h_1D["YAA_pT_fwd_6094"]->scaleW(1./_c["c_YAA_pT_fwd_6094"]->sumW());

      // Figure 2 - ptsq vs Npart
      /////////////////////////


      /////////////////////////

      // Figure 3a - RAA mid and fwd vs pT
      binShift(*_h_RAA_1D["020_pT_mid_CuCu"]);
      binShift(*_h_RAA_1D["pT_mid_pp"]);
      _h_RAA_1D["020_pT_mid_CuCu"]->scaleW(1./_c["c_YAA_pT_mid_020"]->sumW());
      _h_RAA_1D["pT_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
      divide(_h_RAA_1D["020_pT_mid_CuCu"], _h_RAA_1D["pT_mid_pp"],_h2D_RAA["RAA_pT_mid_020"]);
      _h2D_RAA["RAA_pT_mid_020"]->scaleY(1./151.8);  // Ncoll from PHENIX AN 638, page 100

      binShift(*_h_RAA_1D["020_pT_fwd_CuCu"]);
      binShift(*_h_RAA_1D["pT_fwd_pp"]);
       _h_RAA_1D["020_pT_fwd_CuCu"]->scaleW(1./_c["c_YAA_pT_fwd_020"]->sumW());
       _h_RAA_1D["pT_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
       divide(_h_RAA_1D["020_pT_fwd_CuCu"], _h_RAA_1D["pT_fwd_pp"],_h2D_RAA["RAA_pT_fwd_020"]);
       _h2D_RAA["RAA_pT_fwd_020"]->scaleY(1./151.8);

      // // Figure 3b - RAA vs rap
      binShift(*_h_RAA_1D_rap["020_rap_all_CuCu"]);
      binShift(*_h_RAA_1D_rap["rap_all_pp"]);
       _h_RAA_1D_rap["020_rap_all_CuCu"]->scaleW(1./_c["c_CuCu_rap_all_020"]->sumW());
       _h_RAA_1D_rap["rap_all_pp"]->scaleW(1./_c["c_pp"]->sumW());
       divide(_h_RAA_1D_rap["020_rap_all_CuCu"], _h_RAA_1D_rap["rap_all_pp"],_h2D_RAA_rap["RAA_rap_all_020"]);
       _h2D_RAA_rap["RAA_rap_all_020"]->scaleY(1./151.8);  // Ncoll from PHENIX AN 638, page 100

      // Figure 4 - RAA vs Npart
      /////////////////////
      // 0 -10 %
      // _h_RAA_1D_cent_pT["010_pT_fwd_CuCu"]->scaleW(1./_c["c_CuCu_010"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["010_pT_fwd_CuCu"],_h_RAA_1D_cent_pT["pT_cent_fwd_pp"],_h2D_RAA_cent_pT["RAA_pT_fwd_010"]);
      // _h2D_RAA_cent_pT["RAA_pT_fwd_010"]->scaleY(1./98.2); 

      // _h_RAA_1D_cent_pT["010_pT_mid_CuCu"]->scaleW(1./_c["c_CuCu_010"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["010_pT_mid_CuCu"],_h_RAA_1D_cent_pT["pT_cent_mid_pp"],_h2D_RAA_cent_pT["RAA_pT_mid_010"]);
      // _h2D_RAA_cent_pT["RAA_pT_mid_010"]->scaleY(1./98.2); 

      // // 10 - 20 %
      // _h_RAA_1D_cent_pT["1020_pT_fwd_CuCu"]->scaleW(1./_c["c_CuCu_1020"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["1020_pT_fwd_CuCu"],_h_RAA_1D_cent_pT["pT_cent_fwd_pp"],_h2D_RAA_cent_pT["RAA_pT_fwd_1020"]);
      // _h2D_RAA_cent_pT["RAA_pT_fwd_1020"]->scaleY(1./73.6);

      // _h_RAA_1D_cent_pT["1020_pT_mid_CuCu"]->scaleW(1./_c["c_CuCu_1020"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["1020_pT_mid_CuCu"],_h_RAA_1D_cent_pT["pT_cent_mid_pp"],_h2D_RAA_cent_pT["RAA_pT_mid_1020"]);
      // _h2D_RAA_cent_pT["RAA_pT_mid_1020"]->scaleY(1./73.6); 

      // // // 20 - 30 %
      // _h_RAA_1D_cent_pT["2030_pT_fwd_CuCu"]->scaleW(1./_c["c_CuCu_2030"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["2030_pT_fwd_CuCu"],_h_RAA_1D_cent_pT["pT_cent_fwd_pp"],_h2D_RAA_cent_pT["RAA_pT_fwd_2030"]);
      // _h2D_RAA_cent_pT["RAA_pT_fwd_2030"]->scaleY(1./53.0); 

      // _h_RAA_1D_cent_pT["2030_pT_mid_CuCu"]->scaleW(1./_c["c_CuCu_2030"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["2030_pT_mid_CuCu"],_h_RAA_1D_cent_pT["pT_cent_mid_pp"],_h2D_RAA_cent_pT["RAA_pT_mid_2030"]);
      // _h2D_RAA_cent_pT["RAA_pT_mid_2030"]->scaleY(1./53.0); 

      // // // 30 - 40 %
      // _h_RAA_1D_cent_pT["3040_pT_fwd_CuCu"]->scaleW(1./_c["c_CuCu_3040"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["3040_pT_fwd_CuCu"],_h_RAA_1D_cent_pT["pT_cent_fwd_pp"],_h2D_RAA_cent_pT["RAA_pT_fwd_3040"]);
      // _h2D_RAA_cent_pT["RAA_pT_fwd_3040"]->scaleY(1./37.3); 

      // _h_RAA_1D_cent_pT["3040_pT_mid_CuCu"]->scaleW(1./_c["c_CuCu_3040"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["3040_pT_mid_CuCu"],_h_RAA_1D_cent_pT["pT_cent_mid_pp"],_h2D_RAA_cent_pT["RAA_pT_mid_3040"]);
      // _h2D_RAA_cent_pT["RAA_pT_mid_3040"]->scaleY(1./37.3); 

      // // // 40 - 50 %
      // _h_RAA_1D_cent_pT["4050_pT_fwd_CuCu"]->scaleW(1./_c["c_CuCu_4050"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["4050_pT_fwd_CuCu"],_h_RAA_1D_cent_pT["pT_cent_fwd_pp"],_h2D_RAA_cent_pT["RAA_pT_fwd_4050"]);
      // _h2D_RAA_cent_pT["RAA_pT_fwd_4050"]->scaleY(1./25.4); 

      // _h_RAA_1D_cent_pT["4050_pT_mid_CuCu"]->scaleW(1./_c["c_CuCu_4050"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["4050_pT_mid_CuCu"],_h_RAA_1D_cent_pT["pT_cent_mid_pp"],_h2D_RAA_cent_pT["RAA_pT_mid_4050"]);
      // _h2D_RAA_cent_pT["RAA_pT_mid_4050"]->scaleY(1./25.4); 

      // // // 50 - 60 %
      // _h_RAA_1D_cent_pT["5060_pT_fwd_CuCu"]->scaleW(1./_c["c_CuCu_5060"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["5060_pT_fwd_CuCu"],_h_RAA_1D_cent_pT["pT_cent_fwd_pp"],_h2D_RAA_cent_pT["RAA_pT_fwd_5060"]);
      // _h2D_RAA_cent_pT["RAA_pT_fwd_5060"]->scaleY(1./16.7); 

      // _h_RAA_1D_cent_pT["5060_pT_mid_CuCu"]->scaleW(1./_c["c_CuCu_5060"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["5060_pT_mid_CuCu"],_h_RAA_1D_cent_pT["pT_cent_mid_pp"],_h2D_RAA_cent_pT["RAA_pT_mid_5060"]);
      // _h2D_RAA_cent_pT["RAA_pT_mid_5060"]->scaleY(1./16.7); 

      // // 60 - 94 %, 60 - 70%
      // _h_RAA_1D_cent_pT["6070_pT_fwd_CuCu"]->scaleW(1./_c["c_CuCu_6070"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["6070_pT_fwd_CuCu"],_h_RAA_1D_cent_pT["pT_cent_fwd_pp"],_h2D_RAA_cent_pT["RAA_pT_fwd_6070"]);
      // _h2D_RAA_cent_pT["RAA_pT_fwd_6070"]->scaleY(1./10.4); 

      // _h_RAA_1D_cent_pT["6094_pT_mid_CuCu"]->scaleW(1./_c["c_CuCu_6094"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["6094_pT_mid_CuCu"],_h_RAA_1D_cent_pT["pT_cent_mid_pp"],_h2D_RAA_cent_pT["RAA_pT_mid_6094"]);
      // _h2D_RAA_cent_pT["RAA_pT_mid_6094"]->scaleY(1./6.4); 

      // // // 70 - 94% at fwd only
      // _h_RAA_1D_cent_pT["7094_pT_fwd_CuCu"]->scaleW(1./_c["c_CuCu_7094"]->sumW());
      // _h_RAA_1D_cent_pT["pT_cent_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
      // divide(_h_RAA_1D_cent_pT["7094_pT_fwd_CuCu"],_h_RAA_1D_cent_pT["pT_cent_fwd_pp"],_h2D_RAA_cent_pT["RAA_pT_fwd_7094"]);
      // _h2D_RAA_cent_pT["RAA_pT_fwd_7094"]->scaleY(1./4.7);


      // ///////////////////////// 
      // // add together the pT bins to get pT integrated for each centrality histogrma
      // // then use set point outside of a loop and manually set each bin in the Npart 15,1,1 and 14,1,1 hitograms with the corresponding RAA value
            
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_mid_010") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  
      // 		{
      // 		  // cout << "Hello World" << endl;
      // 		  cout << point.y() << endl;

      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		      cout << "Sum RAA pT: " << sum_RAA_pT << endl;
      // 		    }
      // 		}
      // 	    }
	  //  _h2D_RAA_mid_Npart["RAA_mid_Npart"]->point(centBins3[name]).setY(sum_RAA_pT);  // set points inside one histogram
	  //	  _h2D_RAA_mid_Npart["RAA_mid_010"]->point(centBins3[name]).setYErr( sqrt( sum_RAA_pTErr) );
       //	}
      // ### 1 mid
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_mid_1020") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.) 
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_mid_Npart["RAA_mid_1020"]->point(centBins3[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_mid_Npart["RAA_mid_1020"]->point(centBins3[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 2 mid
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_mid_2030") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for mid RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_mid_Npart["RAA_mid_2030"]->point(centBins3[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_mid_Npart["RAA_mid_2030"]->point(centBins3[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 3 mid
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_mid_3040") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for mid RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_mid_Npart["RAA_mid_3040"]->point(centBins3[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_mid_Npart["RAA_mid_3040"]->point(centBins3[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 4 mid
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_mid_4050") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for mid RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_mid_Npart["RAA_mid_4050"]->point(centBins3[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_mid_Npart["RAA_mid_4050"]->point(centBins3[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 5 mid
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_mid_5060") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for mid RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_mid_Npart["RAA_mid_5060"]->point(centBins3[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_mid_Npart["RAA_mid_5060"]->point(centBins3[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 6 mid
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_mid_6094") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for mid RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_mid_Npart["RAA_mid_6094"]->point(centBins3[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_mid_Npart["RAA_mid_6094"]->point(centBins3[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 7 mid
      // /////////////////////////////////////////////////////
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_fwd_010") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for fwd RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_010"]->point(centBins4[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_010"]->point(centBins4[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 1 fwd
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_fwd_1020") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for fwd RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_1020"]->point(centBins4[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_1020"]->point(centBins4[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 2 fwd
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_fwd_2030") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for fwd RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_2030"]->point(centBins4[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_2030"]->point(centBins4[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 3 fwd
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_fwd_3040") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for fwd RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_3040"]->point(centBins4[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_3040"]->point(centBins4[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 4 fwd
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_fwd_4050") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for fwd RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_4050"]->point(centBins4[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_4050"]->point(centBins4[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 5 fwd
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_fwd_5060") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for fwd RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_5060"]->point(centBins4[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_5060"]->point(centBins4[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 6 fwd
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_fwd_6070") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for fwd RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_6070"]->point(centBins4[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_6070"]->point(centBins4[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 7 fwd
      // for(auto element : _h2D_RAA_cent_pT) // for each centrality histogram, sum over pT
      // 	{
      // 	  string name = element.first;
      // 	  if(name.find("RAA_pT_fwd_7094") != std::string::npos) continue;

      // 	  YODA::Scatter2D yoda_RAA = *element.second;

      // 	  double sum_RAA_pT = 0.;
      // 	  double sum_RAA_pTErr = 0.;
      
      // 	  for(auto &point : yoda_RAA.points())
      // 	    {
      // 	      if(point.x() > 0. && point.x() < 5.)  // for fwd RAA
      // 		{
      // 		  if(!isnan(point.y()))
      // 		    {
      // 		      sum_RAA_pT += point.y();
      // 		      sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
      // 		    }
      // 		}
      // 	    }
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_7094"]->point(centBins4[name]).setY(sum_RAA_pT);
      // 	  _h2D_RAA_fwd_Npart["RAA_fwd_7094"]->point(centBins4[name]).setYErr( sqrt( sum_RAA_pTErr) );
      // 	}
      // // ### 8 fwd
    }  // end finalize
    
    /////////////////////////////////////////////////////////////////////////
    
    // FIGURE 1
    map<string, Histo1DPtr> _h_1D;  // for tables 1-8 inv yield
    map<string, CounterPtr> _c;  // for tables 1-8  inv yield
    
    // FIGURE 2
    //map<string, Scatter2DPtr> hRaaNpart;
 
    // FIGURE 3
    map<string, int> centBins1;
    map<string, int> centBins2;
    map<string, Histo1DPtr> _h_RAA_1D;
    map<string, Scatter2DPtr> _h2D_RAA;
    map<string, Histo1DPtr> _h_RAA_1D_rap;
    map<string, Scatter2DPtr> _h2D_RAA_rap;   

    // FIGURE 4
    map<string, int> centBins3;
    map<string, int> centBins4;
    map<string, Histo1DPtr> _h_RAA_1D_cent_pT;
    map<string, Scatter2DPtr> _h2D_RAA_cent_pT;   
    map<string, Scatter2DPtr> _h_RAA_1D_cent; 
    map<string, Scatter2DPtr> _h2D_RAA_mid_Npart;
    map<string, Scatter2DPtr> _h2D_RAA_fwd_Npart;
     
    /////////////////////////////////////////////////////////////////////////
    string beamOpt = "NONE";
    enum CollisionSystem {pp200, CuCu200};
    CollisionSystem CollSys;
    /////////////////////////////////////////////////////////////////////////

    }; // public analysis


    DECLARE_RIVET_PLUGIN(PHENIX_2008_I776624);

  } // namespace Rivet



  // 2D scatter connected to the data
      //  book( _h2D_RAA_mid_Npart["RAA_mid_Npart"], 14,1,1);  
      //  book(_h2D_RAA_mid_Npart["RAA_fwd_Npart", 15,1,1);  
 
