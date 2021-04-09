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
    void init() {

      const UnstableParticles up_fwd(Cuts::absrap < 2.2 && Cuts::absrap > 1.2 && Cuts::abspid == 443 && Cuts::abscharge == 0);  // fwd rap
      declare(up_fwd,"up_fwd");
      const UnstableParticles up_mid(Cuts::absrap < 0.35 && Cuts::abspid == 443 && Cuts::abscharge == 0); // mid rap
      declare(up_mid,"up_mid");
      const UnstableParticles up_all(Cuts::abspid == 443 && Cuts::abscharge == 0); // no rapiidty cuts here
      declare(up_mid,"up_all");

      beamOpt = getOption<string>("beam","NONE");
      
      if(beamOpt=="PP200") CollSys = pp200;
      else if(beamOpt=="AUAU200") CollSys = AuAu200;
     
      if(!(CollSys == pp200)) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

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
      book(_c["c_YAA_pT_mid_020"], "c_YAA_pT_mid_020");  
      book(_c["c_YAA_pT_mid_2040"], "c_YAA_pT_mid_2040");  
      book(_c["c_YAA_pT_mid_4060"], "c_YAA_pT_mid_4060");  
      book(_c["c_YAA_pT_mid_6094"], "c_YAA_pT_mid_6094");  
     
      book(_c["c_YAA_pT_fwd_020"], "c_YAA_pT_fwd_020");
      book(_c["c_YAA_pT_fwd_2040"], "c_YAA_pT_fwd_2040");
      book(_c["c_YAA_pT_fwd_4060"], "c_YAA_pT_fwd_4060");
      book(_c["c_YAA_pT_fwd_6094"], "c_YAA_pT_fwd_6094");

      book(_c["c_pp_pT_fwd"], "c_pp_pT_fwd");
      book(_c["c_pp_pT_mid"], "c_pp_pT_mid");

      book(_c["c_pp_rap_all"], "c_pp_rap_all");
      book(_c["c_AuAu_rap_all_020"], "c_AuAu_rap_all_020");

     //  /////////////////////////////////////////////////////////////
   //    // FIGURE 2
   //    // x-axis for AuAu ptsq Table 9 as function of Npart - up to 5 GeV
   //    vector<double> centBins3{0.0, 20.0, 40.0, 60.0};   // no 60-94 provided for mid rapidity
   //    book(_h1D_npart3["ptsq_mid_cent"], "ptsq_mid_cent", centBins3);  
   //    // x-axis for AuAu ptsq Table 10
   //    vector<double> centBins4{0.0, 20.0, 40.0, 60.0, 94.0};
   //    book(_h1D_npart4["ptsq_fwd_cent"], "ptsq_fwd_cent", centBins4);
   
   // // 2D scatter connected to the data
   //    book( _h2D_ptsq_npart_mid["ptsq_npart_mid"], 9,1,1);  
   //    book( _h2D_ptsq_npart_fwd["ptsq_npart_fwd"], 10,1,1);  
   //    /////////////////////////////////////////////////////////////////////////////////////////

      // FIGURE 3 - Only 0-20% Events shown in figures
    
      string refnameRaa1 = mkAxisCode(11,1,1);
      const Scatter2D& refdataRaa1 =refData(refnameRaa1);
      book(_h_RAA_1D["020_pT_mid_AuAu"], refnameRaa1 + "_AuAu200", refdataRaa1);
      book(_h_RAA_1D["pT_mid_pp"], refnameRaa1 + "_pp200", refdataRaa1);
      book(_h2D_RAA["RAA_pT_mid_020"], refnameRaa1);

      string refnameRaa2 = mkAxisCode(12,1,1);
      const Scatter2D& refdataRaa2 =refData(refnameRaa2);
      book(_h_RAA_1D["020_pT_fwd_AuAu"], refnameRaa2 + "_AuAu200", refdataRaa2);
      book(_h_RAA_1D["pT_fwd_pp"], refnameRaa2 + "_pp200", refdataRaa2);
      book(_h2D_RAA["RAA_pT_fwd_020"], refnameRaa2);

      centBins.insert(pair<string, int>("RAA_pT_mid_020",0));
      centBins.insert(pair<string, int>("RAA_pT_fwd_020",1));
 
      string refnameRaa3 = mkAxisCode(13,1,1);
      const Scatter2D& refdataRaa3 =refData(refnameRaa3);
      book(_h_RAA_1D_rap["020_rap_all_AuAu"], refnameRaa3 + "_AuAu200", refdataRaa3);
      book(_h_RAA_1D_rap["rap_all_pp"], refnameRaa3 + "_pp200", refdataRaa3);
      book(_h2D_RAA_rap["RAA_rap_all_020"], refnameRaa3);

      centBins.insert(pair<string, int>("RAA_rap_all_020",0));

 	/////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////

      // FIGURE 4
   
      book(_c["c_AuAu_010"], "c_AuAu_010");
      book(_c["c_AuAu_1020"], "c_AuAu_1020");
      book(_c["c_AuAu_2030"], "c_AuAu_2030");
      book(_c["c_AuAu_3040"], "c_AuAu_3040");
      book(_c["c_AuAu_4050"], "c_AuAu_4050");
      book(_c["c_AuAu_5060"], "c_AuAu_5060");
      book(_c["c_AuAu_6094"], "c_AuAu_6094");
      book(_c["c_AuAu_6070"], "c_AuAu_6070");
      book(_c["c_AuAu_7094"], "c_AuAu_7094");


      vector<double> pTBins{0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0};  
      book(_h_RAA_1D_cent["010_pT_mid_AuAu"], "010_pT_mid_AuAu", pTBins);
      book(_h_RAA_1D_cent["1020_pT_mid_AuAu"], "1020_pT_mid_AuAu", pTBins);
      book(_h_RAA_1D_cent["2030_pT_mid_AuAu"], "2030_pT_mid_AuAu", pTBins);
      book(_h_RAA_1D_cent["3040_pT_mid_AuAu"], "3040_pT_mid_AuAu", pTBins);
      book(_h_RAA_1D_cent["4050_pT_mid_AuAu"], "4050_pT_mid_AuAu", pTBins);
      book(_h_RAA_1D_cent["5060_pT_mid_AuAu"], "5060_pT_mid_AuAu", pTBins);
      book(_h_RAA_1D_cent["6094_pT_mid_AuAu"], "6094_pT_mid_AuAu", pTBins);

      book(_h_RAA_1D_cent["010_pT_fwd_AuAu"], "010_pT_fwd_AuAu", pTBins);
      book(_h_RAA_1D_cent["1020_pT_fwd_AuAu"], "1020_pT_fwd_AuAu", pTBins);
      book(_h_RAA_1D_cent["2030_pT_fwd_AuAu"], "2030_pT_fwd_AuAu", pTBins);
      book(_h_RAA_1D_cent["3040_pT_fwd_AuAu"], "3040_pT_fwd_AuAu", pTBins);
      book(_h_RAA_1D_cent["4050_pT_fwd_AuAu"], "4050_pT_fwd_AuAu", pTBins);
      book(_h_RAA_1D_cent["5060_pT_fwd_AuAu"], "5060_pT_fwd_AuAu", pTBins);
      book(_h_RAA_1D_cent["6070_pT_fwd_AuAu"], "6070_pT_fwd_AuAu", pTBins);
      book(_h_RAA_1D_cent["7094_pT_fwd_AuAu"], "7094_pT_fwd_AuAu", pTBins);

      book(_h2D_RAA_cent["RAA_pT_mid_010"], "RAA_pT_mid_010", pTBins);
      book(_h2D_RAA_cent["RAA_pT_fwd_010"], "RAA_pT_fwd_010", pTBins);



      // string refnameRaa6 = mkAxisCode(11,1,1);
      // const Scatter2D& refdataRaa6 =refData(refnameRaa6);
      // book(_h_RAA_1D_cent["010_pT_mid_AuAu"], refnameRaa6 + "_AuAu200", refdataRaa6);
      // book(_h2D_RAA_cent["RAA_pT_mid_010"], refnameRaa6);

      // string refnameRaa7 = mkAxisCode(11,1,1);
      // const Scatter2D& refdataRaa7 =refData(refnameRaa7);
      // book(_h_RAA_1D_cent["1020_pT_mid_AuAu"], refnameRaa7 + "_AuAu200", refdataRaa7);
      // book(_h2D_RAA_cent["RAA_pT_mid_1020"], refnameRaa7);

      // string refnameRaa8 = mkAxisCode(11,1,1);
      // const Scatter2D& refdataRaa8 =refData(refnameRaa8);
      // book(_h_RAA_1D_cent["2030_pT_mid_AuAu"], refnameRaa8 + "_AuAu200", refdataRaa8);
      // book(_h2D_RAA_cent["RAA_pT_mid_2030"], refnameRaa8);

      // string refnameRaa9 = mkAxisCode(11,1,1);
      // const Scatter2D& refdataRaa9 =refData(refnameRaa9);
      // book(_h_RAA_1D_cent["3040_pT_mid_AuAu"], refnameRaa9 + "_AuAu200", refdataRaa9);
      // book(_h2D_RAA_cent["RAA_pT_mid_3040"], refnameRaa9);

      // string refnameRaa10 = mkAxisCode(11,1,1);
      // const Scatter2D& refdataRaa10 =refData(refnameRaa10);
      // book(_h_RAA_1D_cent["4050_pT_mid_AuAu"], refnameRaa10 + "_AuAu200", refdataRaa10);
      // book(_h2D_RAA_cent["RAA_pT_mid_4050"], refnameRaa10);

      // string refnameRaa11 = mkAxisCode(11,1,1);
      // const Scatter2D& refdataRaa11 =refData(refnameRaa11);
      // book(_h_RAA_1D_cent["5060_pT_mid_AuAu"], refnameRaa11 + "_AuAu200", refdataRaa11);
      // book(_h2D_RAA_cent["RAA_pT_mid_5060"], refnameRaa11);

      // string refnameRaa12 = mkAxisCode(11,1,1);
      // const Scatter2D& refdataRaa12 =refData(refnameRaa12);
      // book(_h_RAA_1D_cent["6094_pT_mid_AuAu"], refnameRaa12 + "_AuAu200", refdataRaa12);
      // book(_h2D_RAA_cent["RAA_pT_mid_6094"], refnameRaa12);

      // centBins1.insert(pair<string, int>("RAA_pT_mid_010",0));
      // centBins1.insert(pair<string, int>("RAA_pT_mid_1020",1));
      // centBins1.insert(pair<string, int>("RAA_pT_mid_2030",2));
      // centBins1.insert(pair<string, int>("RAA_pT_mid_3040",3));
      // centBins1.insert(pair<string, int>("RAA_pT_mid_4050",4));
      // centBins1.insert(pair<string, int>("RAA_pT_mid_5060",5));
      // centBins1.insert(pair<string, int>("RAA_pT_mid_6094",6));

      // string refnameRaa13 = mkAxisCode(12,1,1);
      // const Scatter2D& refdataRaa13 =refData(refnameRaa13);
      // book(_h_RAA_1D_cent["010_pT_fwd_AuAu"], refnameRaa13 + "_AuAu200", refdataRaa13);
      // book(_h2D_RAA_cent["RAA_pT_fwd_010"], refnameRaa13);

      // string refnameRaa14 = mkAxisCode(12,1,1);
      // const Scatter2D& refdataRaa14 =refData(refnameRaa14);
      // book(_h_RAA_1D_cent["1020_pT_fwd_AuAu"], refnameRaa14 + "_AuAu200", refdataRaa14);
      // book(_h2D_RAA_cent["RAA_pT_fwd_1020"], refnameRaa14);

      // string refnameRaa15 = mkAxisCode(12,1,1);
      // const Scatter2D& refdataRaa15 =refData(refnameRaa15);
      // book(_h_RAA_1D_cent["2030_pT_fwd_AuAu"], refnameRaa15 + "_AuAu200", refdataRaa15);
      // book(_h2D_RAA_cent["RAA_pT_fwd_2030"], refnameRaa15);

      // string refnameRaa16 = mkAxisCode(12,1,1);
      // const Scatter2D& refdataRaa16 =refData(refnameRaa16);
      // book(_h_RAA_1D_cent["3040_pT_fwd_AuAu"], refnameRaa16 + "_AuAu200", refdataRaa16);
      // book(_h2D_RAA_cent["RAA_pT_fwd_3040"], refnameRaa16);

      // string refnameRaa17 = mkAxisCode(12,1,1);
      // const Scatter2D& refdataRaa17 =refData(refnameRaa17);
      // book(_h_RAA_1D_cent["4050_pT_fwd_AuAu"], refnameRaa17 + "_AuAu200", refdataRaa17);
      // book(_h2D_RAA_cent["RAA_pT_fwd_4050"], refnameRaa17);

      // string refnameRaa18 = mkAxisCode(12,1,1);
      // const Scatter2D& refdataRaa18 =refData(refnameRaa18);
      // book(_h_RAA_1D_cent["5060_pT_fwd_AuAu"], refnameRaa18 + "_AuAu200", refdataRaa18);
      // book(_h2D_RAA_cent["RAA_pT_fwd_5060"], refnameRaa18);

      // string refnameRaa19 = mkAxisCode(12,1,1);
      // const Scatter2D& refdataRaa19 =refData(refnameRaa19);
      // book(_h_RAA_1D_cent["6070_pT_fwd_AuAu"], refnameRaa19 + "_AuAu200", refdataRaa19);
      // book(_h2D_RAA_cent["RAA_pT_fwd_6070"], refnameRaa19);

      // string refnameRaa20 = mkAxisCode(12,1,1);
      // const Scatter2D& refdataRaa20 =refData(refnameRaa20);
      // book(_h_RAA_1D_cent["7094_pT_fwd_AuAu"], refnameRaa20 + "_AuAu200", refdataRaa20);
      // book(_h2D_RAA_cent["RAA_pT_fwd_7094"], refnameRaa20);

      // centBins2.insert(pair<string, int>("RAA_pT_fwd_010",0));
      // centBins2.insert(pair<string, int>("RAA_pT_fwd_1020",1));
      // centBins2.insert(pair<string, int>("RAA_pT_fwd_2030",2));
      // centBins2.insert(pair<string, int>("RAA_pT_fwd_3040",3));
      // centBins2.insert(pair<string, int>("RAA_pT_fwd_4050",4));
      // centBins2.insert(pair<string, int>("RAA_pT_fwd_5060",5));
      // centBins2.insert(pair<string, int>("RAA_pT_fwd_6070",6));
      // centBins2.insert(pair<string, int>("RAA_pT_fwd_7094",7));
    
    
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////

      // for pT integrated results using SetPoint
      string refnameRaa4 = mkAxisCode(14,1,1);
      const Scatter2D& refdataRaa4 =refData(refnameRaa4);
      book(_h_RAA_1D_cent["Cent_pT_mid_AuAu"], refnameRaa4 + "_AuAu200", refdataRaa4);
      book(_h_RAA_1D_cent["Cent_pT_mid_pp"], refnameRaa4 + "_pp200", refdataRaa4);
      book(_h2D_RAA_Npart["RAA_pT_mid_Npart"], refnameRaa4);

      string refnameRaa5 = mkAxisCode(15,1,1);
      const Scatter2D& refdataRaa5 =refData(refnameRaa5);
      book(_h_RAA_1D_cent["Cent_pT_fwd_AuAu"], refnameRaa5 + "_AuAu200", refdataRaa5);
      book(_h_RAA_1D_cent["Cent_pT_fwd_pp"], refnameRaa5 + "_pp200", refdataRaa5);
      book(_h2D_RAA_Npart["RAA_pT_fwd_Npart"], refnameRaa5);

      centBins.insert(pair<string, int>("RAA_pT_mid_Npart",0));
      centBins.insert(pair<string, int>("RAA_pT_fwd_Npart",1));
 
    }

    /////////////////////////////////////////////////////////////////////////
    void analyze(const Event& event) {

     
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();
    
      Particles FwdParticles = applyProjection<UnstableParticles>(event,"up_fwd").particles();  // muons
      Particles MidParticles = applyProjection<UnstableParticles>(event,"up_mid").particles();  // electrons
      Particles AllParticles = applyProjection<UnstableParticles>(event,"up_all").particles();  // electrons

      if(CollSys==AuAu200)
      	{
      	  // AuAu event counters 
      	  /////////////////////////////////////////////////////////
      	  if((c > 0.) && (c <= 20.))
      	    {
	      _c["c_YAA_pT_fwd_020"]->fill();
	      _c["c_YAA_pT_mid_020"]->fill();
	      _c["c_AuAu_rap_all_020"]->fill();
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
	  
	  if((c > 0.) && (c <= 10.))
	    _c["c_AuAu_010"]->fill();
	  if((c > 10.) && (c <= 20.))
	    _c["c_AuAu_1020"]->fill();
	  if((c > 20.) && (c <= 30.))
	    _c["c_AuAu_2030"]->fill();
	  if((c > 30.) && (c <= 40.))
	    _c["c_AuAu_3040"]->fill();
	  if((c > 40.) && (c <= 50.))
	    _c["c_AuAu_4050"]->fill();
	  if((c > 50.) && (c <= 60.))
	    _c["c_AuAu_5060"]->fill();
	  if((c > 60.) && (c <= 70.))
	    _c["c_AuAu_6070"]->fill();
	  if((c > 60.) && (c <= 94.))
	    _c["c_AuAu_6094"]->fill();
	  if((c > 70.) && (c <= 94.))
	    _c["c_AuAu_7094"]->fill();
	  
  
      	  // // AuAu fwd J/psi
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
			  _h_1D["YAA_pT_fwd_020"]->fill(jpsi_pT, pt_weight);
			  _h_RAA_1D["020_pT_fwd_AuAu"]->fill(jpsi_pT);
			}
		      // if(jpsi_pT < 5.)
      	  	      // 	_h1D_npart3["ptsq_fwd_cent"]->fill();
      	  	    }
      	  	  else if((c > 20.) && (c <= 40.))
      	  	    {
      	  	      if(jpsi_pT < 6.)
      	  		_h_1D["YAA_pT_fwd_2040"]->fill(jpsi_pT, pt_weight);
      	  	      // if(jpsi_pT < 5.)
      	  	      // 	_h1D_npart3["ptsq_fwd_cent"]->fill();
      	  	    }
      	  	  else if((c > 40.) && (c <= 60.))
      	  	    {
      	  	      if(jpsi_pT < 6.)
      	  		_h_1D["YAA_pT_fwd_4060"]->fill(jpsi_pT, pt_weight);
      	  	      // if(jpsi_pT < 5.)
      	  	      // 	_h1D_npart3["ptsq_fwd_cent"]->fill();
      	  	    }
      	  	  else if((c > 60.) && (c <= 94.))
      	  	    {
      	  	      if(jpsi_pT < 6.)
      	  		_h_1D["YAA_pT_fwd_6094"]->fill(jpsi_pT, pt_weight);
      	  	    }
      	  	  else if(c > 94.) vetoEvent;
      	  	}
      	    }

      	  // // AuAu mid J/psi
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
			  _h_1D["YAA_pT_mid_020"]->fill(jpsi_pT, pt_weight);
			_h_RAA_1D["020_pT_mid_AuAu"]->fill(jpsi_pT);
			}
      	  	      // if(jpsi_pT < 5.)
      	  	      // 	_h1D_npart4["ptsq_mid_cent"]->fill();
      	  	    }
      	  	  else if((c > 20.) && (c <= 40.))
      	  	    {
      	  	      if(jpsi_pT < 10.)
      	  		_h_1D["YAA_pT_mid_2040"]->fill(jpsi_pT, pt_weight);
      	  	      // if(jpsi_pT < 5.)
      	  	      // 	_h1D_npart4["ptsq_mid_cent"]->fill();
      	  	    }
      	  	  else if((c >= 40.) && (c <= 60.))
      	  	    {
      	  	      if(jpsi_pT < 10.)
      	  		_h_1D["YAA_pT_mid_4060"]->fill(jpsi_pT, pt_weight);
		      //  if(jpsi_pT < 5.)
		      //	_h1D_npart4["ptsq_mid_cent"]->fill();
      	  	    }
      	  	  else if((c >= 60.) && (c <= 94.))
      	  	    {
      	  	      if(jpsi_pT < 10.)
      	  		_h_1D["YAA_pT_mid_6094"]->fill(jpsi_pT, pt_weight);
      	  	      // if(jpsi_pT < 5.)
      	  	      // 	_h1D_npart4["ptsq_mid_cent"]->fill();
      	  	    }
      	  	  if(c > 94.) vetoEvent;
      	  	}
      	    }
	 
	  // for pAu RAA vs rapidity, Figure 3
      	  /////////////////////////////////////////////////////////
      	  for(const Particle& p :AllParticles)
      	    {
      	      if(p.pid() == 443)
      		{
		  double jpsi_rap = p.rapidity();
		  
		  _h_RAA_1D_rap["020_rap_all_AuAu"]->fill(jpsi_rap); 

      		}
      	    }

	  // for FIgure 4
	  for(const Particle& p : MidParticles)
      	    {
	      double jpsi_pT = p.pT()/GeV;
	     	    
      	      if(p.pid() == 443)
      	  	{
		  if(jpsi_pT < 10.)
		    {
		      if( (c > 0.) && (c <= 10.))
			_h_RAA_1D_cent["010_pT_mid_AuAu"]->fill(jpsi_pT);	      	  	    	  	    	      	  	    
		      else if((c > 10.) && (c <= 20.))
			_h_RAA_1D_cent["1020_pT_mid_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if((c > 20.) && (c <= 30.))
		      	_h_RAA_1D_cent["2030_pT_mid_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if((c > 30.) && (c <= 40.))
		      	_h_RAA_1D_cent["3040_pT_mid_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if((c > 40.) && (c <= 50.))
		      	_h_RAA_1D_cent["4050_pT_mid_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if((c > 50.) && (c <= 60.))
		      	_h_RAA_1D_cent["5060_pT_mid_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if((c > 60.) && (c <= 94.))
		      	_h_RAA_1D_cent["6094_pT_mid_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if(c > 94.) vetoEvent;
		    }
      	  	}
	    }
	  
	  for(const Particle& p : FwdParticles)
	    {
	      double jpsi_pT = p.pT()/GeV;
	    		  
	      if(p.pid() == 443)
		{
		  if(jpsi_pT < 6.)
		    {
		       if( (c > 0.) && (c <= 10.))
			_h_RAA_1D_cent["010_pT_fwd_AuAu"]->fill(jpsi_pT);	      	  	    	  	    	      	  	    
		      else if((c > 10.) && (c <= 20.))
			_h_RAA_1D_cent["1020_pT_fwd_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if((c > 20.) && (c <= 30.))
		      	_h_RAA_1D_cent["2030_pT_fwd_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if((c > 30.) && (c <= 40.))
		      	_h_RAA_1D_cent["3040_pT_fwd_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if((c > 40.) && (c <= 50.))
		      	_h_RAA_1D_cent["4050_pT_fwd_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if((c > 50.) && (c <= 60.))
		      	_h_RAA_1D_cent["5060_pT_fwd_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if((c > 60.) && (c <= 70.))
		      	_h_RAA_1D_cent["6070_pT_fwd_AuAu"]->fill(jpsi_pT);
		      else if((c > 70.) && (c <= 94.))
		      	_h_RAA_1D_cent["7094_pT_fwd_AuAu"]->fill(jpsi_pT);	      	  	    	  	    
		      else if(c > 94.) vetoEvent;
		    }
		}
	    }
	} // if AuAu

      if(CollSys == pp200)
      	{

	    _c["c_pp"]->fill();
	

      	  // for pp RAA vs rapidity, Figure 3
      	  /////////////////////////////////////////////////////////
      	  for(const Particle& p :AllParticles)
      	    {
	   
      	      if(p.pid() == 443)
      		{
		  double jpsi_rap = p.rapidity();
		  
		  _h_RAA_1D_rap["rap_all_pp"]->fill(jpsi_rap);  // mid

      		}
      	    }
      	  // for pp RAA fwd
      	  /////////////////////////////////////////////////////////
      	  for(const Particle& p : FwdParticles)
      	    {

	      double jpsi_pT = p.pT()/GeV;
	    
      	      if((p.pid() == 443) && (jpsi_pT < 6.))
      		{
		  //	  _h1D_npart8_pp["pp_fwd_npart"]->fill();
      		  _h1D_pT_fwd_pp["pT_bins5_fwd_pp"]->fill(jpsi_pT);
		  _h_RAA_1D["pT_fwd_pp"]->fill(jpsi_pT);
		}
      	    }
      	  // for pp RAA mid
      	  /////////////////////////////////////////////////////////
      	  for(const Particle& p : MidParticles)
      	    {

	      double jpsi_pT = p.pT()/GeV;
	     
	      if((p.pid() == 443) && (jpsi_pT < 10.))
      		{
		  //  _h1D_npart7_pp["pp_mid_npart"]->fill();
      		  _h1D_pT_mid_pp["pT_bins5_mid_pp"]->fill(jpsi_pT);
		  _h_RAA_1D["pT_mid_pp"]->fill(jpsi_pT);
		}
      	    }

      	} // if pp


    } // end analyze

    /////////////////////////////////////////////////////////////////////////
    void finalize() {


      // Figure 1  - yields vs pT
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
      _h_RAA_1D["020_pT_mid_AuAu"]->scaleW(1./_c["c_YAA_pT_mid_020"]->sumW());
      _h_RAA_1D["pT_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
      divide(_h_RAA_1D["020_pT_mid_AuAu"], _h_RAA_1D["pT_mid_pp"],_h2D_RAA["RAA_pT_mid_020"]);
      _h2D_RAA["RAA_pT_mid_020"]->scaleY(1./151.8);  // Ncoll from PHENIX AN 638, page 100

      _h_RAA_1D["020_pT_fwd_AuAu"]->scaleW(1./_c["c_YAA_pT_fwd_020"]->sumW());
      _h_RAA_1D["pT_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
      divide(_h_RAA_1D["020_pT_fwd_AuAu"], _h_RAA_1D["pT_fwd_pp"],_h2D_RAA["RAA_pT_fwd_020"]);
      _h2D_RAA["RAA_pT_fwd_020"]->scaleY(1./151.8);

      // Figure 3b - RAA vs rap
      _h_RAA_1D_rap["020_rap_all_AuAu"]->scaleW(1./_c["c_AuAu_rap_all_020"]->sumW());
      _h_RAA_1D_rap["rap_all_pp"]->scaleW(1./_c["c_pp"]->sumW());
      divide(_h_RAA_1D_rap["020_rap_all_AuAu"], _h_RAA_1D_rap["rap_all_pp"],_h2D_RAA_rap["RAA_rap_all_020"]);
      _h2D_RAA_rap["RAA_rap_all_020"]->scaleY(1./151.8);  // Ncoll from PHENIX AN 638, page 100

      // Figure 4 - RAA vs Npart
      /////////////////////
      // 0 -10 %
      _h_RAA_1D_cent["010_pT_fwd_AuAu"]->scaleW(1./_c["c_AuAu_010"]->sumW());
      _h_RAA_1D_cent["pT_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
       divide(_h_RAA_1D_cent["010_pT_fwd_AuAu"],_h_RAA_1D_cent["pT_fwd_pp"],_h2D_RAA_cent["RAA_pT_fwd_010"]);
       _h2D_RAA_cent["RAA_pT_fwd_010"]->scaleY(1./98.2); 

       _h_RAA_1D_cent["010_pT_mid_AuAu"]->scaleW(1./_c["c_AuAu_010"]->sumW());
      _h_RAA_1D_cent["pT_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
       divide(_h_RAA_1D_cent["010_pT_mid_AuAu"],_h_RAA_1D_cent["pT_mid_pp"],_h2D_RAA_cent["RAA_pT_mid_010"]);
       _h2D_RAA_cent["RAA_pT_mid_010"]->scaleY(1./98.2); 

      //  // 10 - 20 %
//        _h_RAA_1D_cent["1020_pT_fwd_AuAu"]->scaleW(1./_c["c_AuAu_1020"]->sumW());
//       _h_RAA_1D_cent["pT_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["1020_pT_fwd_AuAu"],_h_RAA_1D_cent["pT_fwd_pp"],_h2D_RAA_cent["RAA_pT_fwd_1020"]);
//        _h2D_RAA_cent["RAA_pT_fwd_1020"]->scaleY(1./73.6);

//        _h_RAA_1D_cent["1020_pT_mid_AuAu"]->scaleW(1./_c["c_AuAu_1020"]->sumW());
//       _h_RAA_1D_cent["pT_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["1020_pT_mid_AuAu"],_h_RAA_1D_cent["pT_mid_pp"],_h2D_RAA_cent["RAA_pT_mid_1020"]);
//        _h2D_RAA_cent["RAA_pT_mid_1020"]->scaleY(1./73.6); 

//        // 20 - 30 %
//        _h_RAA_1D_cent["2030_pT_fwd_AuAu"]->scaleW(1./_c["c_AuAu_2030"]->sumW());
//       _h_RAA_1D_cent["pT_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["2030_pT_fwd_AuAu"],_h_RAA_1D_cent["pT_fwd_pp"],_h2D_RAA_cent["RAA_pT_fwd_2030"]);
//        _h2D_RAA_cent["RAA_pT_fwd_2030"]->scaleY(1./53.0); 

//        _h_RAA_1D_cent["2030_pT_mid_AuAu"]->scaleW(1./_c["c_AuAu_2030"]->sumW());
//       _h_RAA_1D_cent["pT_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["2030_pT_mid_AuAu"],_h_RAA_1D_cent["pT_mid_pp"],_h2D_RAA_cent["RAA_pT_mid_2030"]);
//        _h2D_RAA_cent["RAA_pT_mid_2030"]->scaleY(1./53.0); 

//        // 30 - 40 %
//        _h_RAA_1D_cent["3040_pT_fwd_AuAu"]->scaleW(1./_c["c_AuAu_3040"]->sumW());
//       _h_RAA_1D_cent["pT_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["3040_pT_fwd_AuAu"],_h_RAA_1D_cent["pT_fwd_pp"],_h2D_RAA_cent["RAA_pT_fwd_3040"]);
//        _h2D_RAA_cent["RAA_pT_fwd_3040"]->scaleY(1./37.3); 

//        _h_RAA_1D_cent["3040_pT_mid_AuAu"]->scaleW(1./_c["c_AuAu_3040"]->sumW());
//       _h_RAA_1D_cent["pT_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["3040_pT_mid_AuAu"],_h_RAA_1D_cent["pT_mid_pp"],_h2D_RAA_cent["RAA_pT_mid_3040"]);
//        _h2D_RAA_cent["RAA_pT_mid_3040"]->scaleY(1./37.3); 

//        // 40 - 50 %
//        _h_RAA_1D_cent["4050_pT_fwd_AuAu"]->scaleW(1./_c["c_AuAu_4050"]->sumW());
//       _h_RAA_1D_cent["pT_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["4050_pT_fwd_AuAu"],_h_RAA_1D_cent["pT_fwd_pp"],_h2D_RAA_cent["RAA_pT_fwd_4050"]);
//        _h2D_RAA_cent["RAA_pT_fwd_4050"]->scaleY(1./25.4); 

//        _h_RAA_1D_cent["4050_pT_mid_AuAu"]->scaleW(1./_c["c_AuAu_4050"]->sumW());
//       _h_RAA_1D_cent["pT_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["4050_pT_mid_AuAu"],_h_RAA_1D_cent["pT_mid_pp"],_h2D_RAA_cent["RAA_pT_mid_4050"]);
//        _h2D_RAA_cent["RAA_pT_mid_4050"]->scaleY(1./25.4); 

//        // 50 - 60 %
//        _h_RAA_1D_cent["5060_pT_fwd_AuAu"]->scaleW(1./_c["c_AuAu_5060"]->sumW());
//       _h_RAA_1D_cent["pT_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["5060_pT_fwd_AuAu"],_h_RAA_1D_cent["pT_fwd_pp"],_h2D_RAA_cent["RAA_pT_fwd_5060"]);
//        _h2D_RAA_cent["RAA_pT_fwd_5060"]->scaleY(1./16.7); 

//        _h_RAA_1D_cent["5060_pT_mid_AuAu"]->scaleW(1./_c["c_AuAu_5060"]->sumW());
//       _h_RAA_1D_cent["pT_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["5060_pT_mid_AuAu"],_h_RAA_1D_cent["pT_mid_pp"],_h2D_RAA_cent["RAA_pT_mid_5060"]);
//        _h2D_RAA_cent["RAA_pT_mid_5060"]->scaleY(1./16.7); 

//        // 60 - 94 %, 60 - 70%
//        _h_RAA_1D_cent["6070_pT_fwd_AuAu"]->scaleW(1./_c["c_AuAu_6070"]->sumW());
//       _h_RAA_1D_cent["pT_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["6070_pT_fwd_AuAu"],_h_RAA_1D_cent["pT_fwd_pp"],_h2D_RAA_cent["RAA_pT_fwd_6070"]);
//        _h2D_RAA_cent["RAA_pT_fwd_6070"]->scaleY(1./10.4); 

//        _h_RAA_1D_cent["6094_pT_mid_AuAu"]->scaleW(1./_c["c_AuAu_6094"]->sumW());
//       _h_RAA_1D_cent["pT_mid_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["6094_pT_mid_AuAu"],_h_RAA_1D_cent["pT_mid_pp"],_h2D_RAA_cent["RAA_pT_mid_6094"]);
//        _h2D_RAA_cent["RAA_pT_mid_6094"]->scaleY(1./6.4); 

//        // 70 - 94% at fwd only
//        _h_RAA_1D_cent["7094_pT_fwd_AuAu"]->scaleW(1./_c["c_AuAu_7094"]->sumW());
//       _h_RAA_1D_cent["pT_fwd_pp"]->scaleW(1./_c["c_pp"]->sumW());
//        divide(_h_RAA_1D_cent["7094_pT_fwd_AuAu"],_h_RAA_1D_cent["pT_fwd_pp"],_h2D_RAA_cent["RAA_pT_fwd_7094"]);
//        _h2D_RAA_cent["RAA_pT_fwd_7094"]->scaleY(1./4.7);


// ///////////////////////// 
//       // add together the pT bins to get pT integrated for each centrality histogrma
    
//       // then use set point outside of a loop and manually set each bin in the Npart 15,1,1 and 14,1,1 hitograms with the corresponding RAA value
//        // then set the point in the 14, 1,1 and 15,1,1 RAA hsitograms using code below
       
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_mid_010") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 5.)  
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_mid_010"]->point(centBins1[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_mid_010"]->point(centBins1[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 1 mid
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_mid_1020") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 5.) 
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_mid_1020"]->point(centBins1[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_mid_1020"]->point(centBins1[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 2 mid
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_mid_2030") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 5.)  // for mid RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_mid_2030"]->point(centBins1[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_mid_2030"]->point(centBins1[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 3 mid
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_mid_3040") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 5.)  // for mid RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_mid_3040"]->point(centBins1[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_mid_3040"]->point(centBins1[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 4 mid
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_mid_4050") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 5.)  // for mid RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_mid_4050"]->point(centBins1[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_mid_4050"]->point(centBins1[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 5 mid
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_mid_5060") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 5.)  // for mid RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_mid_5060"]->point(centBins1[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_mid_5060"]->point(centBins1[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 6 mid
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_mid_6094") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 5.)  // for mid RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_mid_6094"]->point(centBins1[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_mid_6094"]->point(centBins1[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 7 mid
//        /////////////////////////////////////////////////////
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_fwd_010") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 6.)  // for fwd RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_fwd_010"]->point(centBins2[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_fwd_010"]->point(centBins2[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 1 fwd
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_fwd_1020") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 6.)  // for fwd RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_fwd_1020"]->point(centBins2[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_fwd_1020"]->point(centBins2[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 2 fwd
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_fwd_2030") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 6.)  // for fwd RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_fwd_2030"]->point(centBins2[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_fwd_2030"]->point(centBins2[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 3 fwd
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_fwd_3040") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 6.)  // for fwd RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_fwd_3040"]->point(centBins2[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_fwd_3040"]->point(centBins2[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 4 fwd
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_fwd_4050") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 6.)  // for fwd RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_fwd_4050"]->point(centBins2[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_fwd_4050"]->point(centBins2[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 5 fwd
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_fwd_5060") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 6.)  // for fwd RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_fwd_5060"]->point(centBins2[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_fwd_5060"]->point(centBins2[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 6 fwd
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_fwd_6070") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 6.)  // for fwd RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_fwd_6070"]->point(centBins2[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_fwd_6070"]->point(centBins2[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 7 fwd
//        for(auto element : _h2D_RAA_cent) // for each centrality histogram, sum over pT
// 	 {
// 	   string name = element.first;
// 	   if(name.find("RAA_pT_fwd_7094") != std::string::npos) continue;

// 	   YODA::Scatter2D yoda_RAA = *element.second;

// 	   double sum_RAA_pT = 0.;
// 	   double sum_RAA_pTErr = 0.;
      
// 	   for(auto &point : yoda_RAA.points())
// 	     {
// 	       if(point.x() > 0. && point.x() < 6.)  // for fwd RAA
// 		 {
// 		   if(!isnan(point.y()))
// 		     {
// 		       sum_RAA_pT += point.y();
// 		       sum_RAA_pTErr += point.yErrPlus()*point.yErrPlus();
// 		     }
// 		 }
// 	     }
// 	   _h2D_RAA_cent["RAA_pT_fwd_7094"]->point(centBins2[name]).setY(sum_RAA_pT);
// 	   _h2D_RAA_cent["RAA_pT_fwd_7094"]->point(centBins2[name]).setYErr( sqrt( sum_RAA_pTErr) );
// 	 }
//        // ### 7 fwd


    }  // end finalize
    
    /////////////////////////////////////////////////////////////////////////
    
    // FIGURE 1
    map<string, Histo1DPtr> _h_1D;  // for tables 1-8 inv yield
    map<string, CounterPtr> _c;  // for tables 1-8  inv yield
    
    // FIGURE 2
    //map<string, Scatter2DPtr> hRaaNpart;
    map<string, int> centBins;
    map<string, int> centBins1;
    map<string, int> centBins2;
    map<string, Histo1DPtr> _h_RAA_1D;
    map<string, Scatter2DPtr> _h2D_RAA;
 
    // FIGURE 3
    map<string, Histo1DPtr> _h1D_pT_mid_pp;  // denom for table 11  RAA
    map<string, Histo1DPtr> _h1D_pT_fwd_pp;  // denom for table 12  RAA
    map<string, Histo1DPtr> _h1D_rap_pp;  // denom for table 13 RAA
    map<string, Scatter2DPtr> _h2D_RAA_pt_fwd;   
    map<string, Scatter2DPtr> _h2D_RAA_pt_mid;   
    map<string, Histo1DPtr> _h_RAA_1D_rap;
    map<string, Scatter2DPtr> _h2D_RAA_rap;   
    map<string, Histo1DPtr> _h_RAA_1D_cent;
    map<string, Scatter2DPtr> _h2D_RAA_cent;   

    // FIGURE 4
    // map<string, Histo1DPtr> _h1D_npart7;  //  num for table 11-12 RAA
    // map<string, Histo1DPtr> _h1D_npart8;  //num for table 11-12 RAA
    // map<string, Histo1DPtr> _h1D_npart7_pp;  // for denom of RAA
    // map<string, Histo1DPtr> _h1D_npart8_pp;  // for denom RAA
    map<string, Scatter2DPtr> _h2D_RAA_Npart;
     
    /////////////////////////////////////////////////////////////////////////
    string beamOpt;
    enum CollisionSystem {pp200, AuAu200};
   CollisionSystem CollSys;
    /////////////////////////////////////////////////////////////////////////

  }; // public analysis


  DECLARE_RIVET_PLUGIN(PHENIX_2008_I776624);

} // namespace Rivet



