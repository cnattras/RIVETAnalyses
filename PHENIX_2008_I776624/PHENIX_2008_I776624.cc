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

      beamOpt = getOption<string>("beam","NONE");
      
      if(beamOpt=="PP200") collSys = pp200;
      else if(beamOpt=="CUCU200") collSys = CuCu200;
     
      if(!(collSys == pp200)) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      //declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      /////////////////////////////////////////////////////////////
     //  FIGURE 1 (Tables 1-8)
      book(_h_1D["YAA_pT_mid_020"], 1, 1, 1);  
      book(_h_1D["YAA_pT_mid_2040"], 2, 1, 1);  
      book(_h_1D["YAA_pT_mid_4060"], 3, 1, 1); 
      book(_h_1D["YAA_pT_mid_6094"], 4, 1, 1); 
     
      book(_h_1D["YAA_pT_fwd_020"], 5, 1, 1); 
      book(_h_1D["YAA_pT_fwd_2040"], 6, 1, 1); 
      book(_h_1D["YAA_pT_fwd_4060"], 7, 1, 1); 
      book(_h_1D["YAA_pT_fwd_6094"], 8, 1, 1); 

     // counter histograms
      book(_c["c_YAA_mid_020"], "c_YAA_mid_020");  
      book(_c["c_YAA_mid_2040"], "c_YAA_mid_2040");  
      book(_c["c_YAA_mid_4060"], "c_YAA_mid_4060");  
      book(_c["c_YAA_mid_6094"], "c_YAA_mid_6094");  
     
      book(_c["c_YAA_fwd_020"], "c_YAA_fwd_020");
      book(_c["c_YAA_fwd_2040"], "c_YAA_fwd_2040");
      book(_c["c_YAA_fwd_4060"], "c_YAA_fwd_4060");
      book(_c["c_YAA_fwd_6094"], "c_YAA_fwd_6094");

      /////////////////////////////////////////////////////////////
      // FIGURE 2 
      // x-axis for CuCu ptsq Table 9
      std::vector<double> centBins3 = {0., 20., 40., 60.};   // no 60-94 provided for mid rapidity
      book(_h1D_npart3, "ptsq_mid_cent", 3, centBins3);   // Antoniio slack instructions
      // x-axis for CuCu ptsq Table 10
      std::vector<double> centBins4 = {0., 20., 40., 60., 94.};
      book(_h1D_npart4, "ptsq_fwd_cent", 4, centBins4);
   
      /////////////////////////////////////////////////////////////////////////////////////////

      // FIGURE 3
      // x-axis for Tables 11 and 12 - pp vs. pT
      std::vector<double> pT_bins5_pp = {0., 1., 2., 3., 4., 5.};  
      book(_h1D_pT_pp, "pT_bins5_pp", 5, pT_bins5_pp);
  
      //  x-axis for Table 13 - pp vs. rapidity
      std::vector<double> pT_bins5_ = {0., 1., 2., 3., 4., 5.};  
      book(_h1D_rap_pp, "rap_bins9_pp", 9, rap_bins9_rap);
 
      // 2D scatter connected to the data
      book( _h2D_RAA_pt_mid["RAA_pT_mid_020"],11,1,1);  
      book( _h2D_RAA_pt_fwd["RAA_pT_fwd_020"],12,1,1);  
      book( _h2D_RAA_rap["RAA_rap_020"],13,1,1);  
 
 	/////////////////////////////////////////////////////////////////////////////////////////

      // FIGURE 4
      // x-axis for CuCu Table 14 mid
      std::vector<double> centBins7 = {0., 10., 20., 30., 40., 50., 60., 94.};  // dummy 1D histo
      book(_h1D_npart7, "CuCu_mid_npart", 7, centBins7);
      // x-axis for CuCu Table 15 fwd
      std::vector<double> centBins8 = {0., 10., 20., 30., 40., 50., 60., 70., 94.}; // dummy 1D histo
      book(_h1D_npart8, "CuCu_fwd_npart", 8, centBins8);
   
       // x-axis for pp Tble 14
      std::vector<double> centBins7_pp = {0., 10., 20., 30., 40., 50., 60., 94.};   // dummy 1D histo
      book(_h1D_npart7_pp, "pp_mid_npart", 7, centBins7_pp);
      // x-axis for pp Table 15
      std::vector<double> centBins8_pp = {0., 10., 20., 30., 40., 50., 60., 70., 94.};   // dummy 1D histo
      book(_h1D_npart8_pp, "pp_fwd_npart", 8, centBins8_pp);

       // 2D scatter connected to the data
      book( _h2D_RAA_npart_mid["RAA_npart_mid"],14,1,1);  
      book( _h2D_RAA_npart_fwd["RAA_npart_fwd"],15,1,1);  
    }

    /////////////////////////////////////////////////////////////////////////
    void analyze(const Event& event) {

     
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();

      Particles FwdParticles = applyProjection<UnstableParticles>(event,"upFwd").particles();  // muons
      Particles MidParticles = applyProjection<UnstableParticles>(event,"upMid").particles();  // electrons

      if(collSys==CuCu200)
      	{

      	  // Fill counters (no associated particle types?)
      	  if(c < 20.)
      	    {
      	      _c["c_YAA_fwd_020"]->fill();
      	      _c["c_YAA_mid_020"]->fill();
      	    }
      	  else if(c >= 20. && c < 40.)
      	    {
      	      _c["c_YAA_fwd_2040"]->fill();
      	      _c["c_YAA_mid_2040"]->fill();
      	    }
      	  else if(c >= 40. && c < 60.)
      	    {
      	      _c["c_YAA_fwd_4060"]->fill();
      	      _c["c_YAA_mid_4060"]->fill();
      	    }
      	  else if(c >= 60. && c < 94.)
      	    {
      	      _c["c_YAA_fwd_6094"]->fill();
      	      _c["c_YAA_mid_6094"]->fill();
      	    }
      	  if(c > 94.) vetoEvent;


    	  // fill hisograms ( with associated particle types )
      for(const Particle& p : FwdParticles)
    	{
    	  if(c < 20. && p.pid() == 443)
    	    _h_1D["YAA_pT_fwd_020"]->fill(p.pT()/GeV);
    	  else if(c >= 20. && c < 40. && p.pid() == 443)
    	    _h_1D["YAA_pT_fwd_2040"]->fill(p.pT()/GeV);
    	  else if(c >= 40. && c < 60. && p.pid() == 443)
    	    _h_1D["YAA_pT_fwd_4060"]->fill(p.pT()/GeV);
    	  else if(c >= 60. && c < 94. && p.pid() == 443)
    	    _h_1D["YAA_pT_fwd_6094"]->fill(p.pT()/GeV);
    	  if(c > 94.) vetoEvent;
    	}
      return;

      for(const Particle& p : MidParticles)
    	{
    	  if(c < 20. && p.pid() == 443)
    	    _h_1D["YAA_pT_mid_020"]->fill(p.pT()/GeV);
    	  else if(c >= 20. && c < 40. && p.pid() == 443)
    	    _h_1D["YAA_pT_mid_2040"]->fill(p.pT()/GeV);
    	  else if(c >= 40. && c < 60. && p.pid() == 443)
    	    _h_1D["YAA_pT_mid_4060"]->fill(p.pT()/GeV);
    	  else if(c >= 60. && c < 94. && p.pid() == 443)
    	    _h_1D["YAA_pT_mid_6094"]->fill(p.pT()/GeV);
    	  if(c > 94.) vetoEvent;
    	}
      return;
    } // if CuCu


      // Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();

      // for(const Particle& p : fsParticles) 
      // 	{
      // 	  if(c < 10. && p.pid() == 321) _h["AAAA"]->fill(p.pT()/GeV);
      // 	}


    } // end analyze

    /////////////////////////////////////////////////////////////////////////
    void finalize() {

  

    }

    /////////////////////////////////////////////////////////////////////////

    // FIGURE 1
    map<string, Histo1DPtr> _h_1D;  // for tables 1-8 inv yield
    map<string, CounterPtr>  _c;  // for tables 1-8  inv yield

    // FIGURE 2
    map<string, Histo1DPtr> _h1D_npart3;  // for table 9 ptsq 
    map<string, Histo1DPtr> _h1D_npart4;   // for table 10 ptsq 
   
    // FIGURE 3
    map<string, Histo1DPtr>_h1D_pT_pp;  // denom for table 11-12 RAA
    map<string, Histo1DPtr>_h1D_rap_pp;  // denom for table 13 RAA
    map<string, Scatter2DPtr> _h2D_RAA_pt_fwd;   
    map<string, Scatter2DPtr> _h2D_RAA_pt_mid;   
    map<string, Scatter2DPtr> _h2D_RAA_rap;   

    // FIGURE 4
    map<string, Histo1DPtr> _h1D_npart7;  //  num of for table 11-12 RAA
    map<string, Histo1DPtr> _h1D_npart8;  //num of for table 11-12 RAA
    map<string, Histo1DPtr> _h1D_npart7_pp;  // for denom of RAA
    map<string, Histo1DPtr> _h1D_npart8_pp;  // for denom RAA
     map<string, Scatter2DPtr> _h2D_RAA_npart_fwd;   
    map<string, Scatter2DPtr> _h2D_RAA_npart_mid;   
   
    /////////////////////////////////////////////////////////////////////////
    string beamOpt;
    enum CollisionSystem {pp200, CuCu200};
   CollisionSystem collSys;
    /////////////////////////////////////////////////////////////////////////

  }; // public analysis


  DECLARE_RIVET_PLUGIN(PHENIX_2008_I776624);

} // namespace Rivet


