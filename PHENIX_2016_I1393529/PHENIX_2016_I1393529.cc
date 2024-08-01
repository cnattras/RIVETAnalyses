// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2016_I1393529 : public Analysis {
  public:
    
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2016_I1393529);
    
    
    /// @name Analysis methods
    //@{
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
    
    /// Book histograms and initialise projections before the run
    void init() {
     
      beamOpt = getOption<string>("beam", "NONE");
 
      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
      
      const FinalState fs(Cuts::abseta < 0.35 && Cuts::pT > 0.0*GeV && Cuts::pT < 20.0*GeV);
      declare(fs, "fs");

      //counters 
      book(_c["pp"], "_pp");      
      book(_c["AuAu"], "_AuAu");      
      
      //Figure 1
      book(_h["InvYield_charm"], 1, 1, 1);
      book(_h["InvYield_bottom"], 1, 1, 2);
      
      //+++++++++++++++++++//
      //Figure 2
      string refname_1 = mkAxisCode(2, 1, 1);
      const Scatter2D& refdata_1 = refData(2, 1, 1);

      book(_h["bfracn"], refname_1 + "_num", refdata_1);
      book(_h["bfracd"], refname_1 + "_den", refdata_1);
      book(_s["bfrac"], refname_1, true);

      //+++++++++++++++++++//
      //Figure 3
      string refname_2 = mkAxisCode(3, 1, 1);
      const Scatter2D& refdata_2 = refData(3, 1, 1);

      book(_h["RAA_c2en"], refname_2 + "_num", refdata_2);
      book(_h["RAA_c2ed"], refname_2 + "_den", refdata_2);
      book(_s["RAA_c2e"], refname_2, true);

      
      string refname_3 = mkAxisCode(3, 1, 2);
      const Scatter2D& refdata_3 = refData(3, 1, 2);
      
      book(_h["RAA_b2en"], refname_3 + "_num", refdata_3);
      book(_h["RAA_b2ed"], refname_3 + "_den", refdata_3);
      book(_s["RAA_b2e"], refname_3, true);
      
      //+++++++++++++++++++//
      string refname_4 = mkAxisCode(4, 1, 1);
      const Scatter2D& refdata_4 = refData(4, 1, 1);
      
      book(_h["RAA_ratio_c2en"], refname_4 + "_c2e_num", refdata_4);
      book(_h["RAA_ratio_c2ed"], refname_4 + "_c2e_den", refdata_4);
      book(_s["RAA_ratio_c2e"], refname_4 + "_c2e");
      
      book(_h["RAA_ratio_b2en"], refname_4 + "_b2e_num", refdata_4);
      book(_h["RAA_ratio_b2ed"], refname_4 + "_b2e_den", refdata_4);
      book(_s["RAA_ratio_b2e"], refname_4 + "_b2e");
      book(_s["RAA_ratio"], refname_4);
      
      
    }
    
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
      const double c=cent();
//      cout<<"cent: "<<c<<endl;
      
      const ParticlePair& beam = beams();
      double NN=0;

      if (beamOpt == "NONE") {
        if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
	  {
	    NN = 197.;
	    if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = AuAu200;
	    //if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 5)) beamName += "62GeV";
	  }
        if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
	  {
	    NN = 1.;
	    if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = pp_200;
	  }
      }
      else if (beamOpt == "AUAU200") collSys = AuAu200;
      else if (beamOpt == "PP200") collSys = pp_200;      

      Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();
      
      for(const Particle& p : fsParticles) 
	{
	  if(p.fromCharm()){
	    _h["InvYield_charm"]->fill(p.pT()/GeV);
	    if(collSys==AuAu200){
	      _c["AuAu"]->fill();
	      _h["RAA_c2en"]->fill(p.pT()/GeV);
	      _h["RAA_ratio_c2en"]->fill(p.pT()/GeV);
	    }
	    else if(collSys==pp_200){
	      _c["pp"]->fill();
	      _h["RAA_c2ed"]->fill(p.pT()/GeV);
	      _h["RAA_ratio_c2ed"]->fill(p.pT()/GeV);
	    }
	  }
	  if(p.fromBottom()){
	    _h["InvYield_bottom"]->fill(p.pT()/GeV);
	    _h["bfracn"]->fill(p.pT()/GeV);
	    if(collSys==AuAu200){
	      _c["AuAu"]->fill();
	      _h["RAA_b2en"]->fill(p.pT()/GeV);
	      _h["RAA_ratio_b2en"]->fill(p.pT()/GeV);
	    }
	    else if(collSys==pp_200){
	      _c["pp"]->fill();
	      _h["RAA_b2ed"]->fill(p.pT()/GeV);
	      _h["RAA_ratio_b2ed"]->fill(p.pT()/GeV);
	    }
	  }
	  if(p.fromCharm() || p.fromBottom()) _h["bfracd"]->fill(p.pT()/GeV);
	  
	  
	}
      
      
    }
    
    
    /// Normalise histograms etc., after the run
    void finalize() {
      bool has_AuAu {_c["AuAu"]->sumW() > 0 };
      bool has_pp {_c["pp"]->sumW() > 0 };
      
      double scale = 1./(2*M_PI);
      if(has_pp){
	//denominators for the RAA Plots -- Fig 3
  binShift(*_h["RAA_c2ed"]);
  binShift(*_h["RAA_b2ed"]);
	_h["RAA_c2ed"]->scaleW(scale/_c["pp"]->sumW());
	_h["RAA_b2ed"]->scaleW(scale/_c["pp"]->sumW());
      }
      else if(has_AuAu){
	//Figure 1
  binShift(*_h["InvYield_charm"]);
  binShift(*_h["InvYield_bottom"]);
	_h["InvYield_charm"]->scaleW(scale/_c["AuAu"]->sumW());
	_h["InvYield_bottom"]->scaleW(scale/_c["AuAu"]->sumW());
	
	//numerator and denominator for Figure 2
  binShift(*_h["bfracn"]);
  binShift(*_h["bfracd"]);
	_h["bfracn"]->scaleW(scale/_c["AuAu"]->sumW());
	_h["bfracd"]->scaleW(scale/_c["AuAu"]->sumW());
	
	//numerators for the RAA Plots -- Fig 3
  binShift(*_h["RAA_c2en"]);
  binShift(*_h["RAA_b2en"]);
	_h["RAA_c2en"]->scaleW(scale/_c["AuAu"]->sumW());
	_h["RAA_b2en"]->scaleW(scale/_c["AuAu"]->sumW());
	
      }

      
      if(has_pp && has_AuAu){
	//Ratio plots
	//Figure 2
	divide(_h["bfracn"], _h["bfracd"], _s["bfrac"]);
	
	//Figure 3 a and b
	divide(_h["RAA_c2en"], _h["RAA_c2ed"], _s["RAA_c2e"]);
	divide(_h["RAA_b2en"], _h["RAA_b2ed"], _s["RAA_b2e"]);
	
	//Figure 4
  binShift(*_h["RAA_ratio_c2en"]);
  binShift(*_h["RAA_ratio_c2ed"]);
  binShift(*_h["RAA_ratio_b2en"]);
  binShift(*_h["RAA_ratio_b2ed"]);
	divide(_h["RAA_ratio_c2en"], _h["RAA_ratio_c2ed"], _s["RAA_ratio_c2e"]);
	divide(_h["RAA_ratio_b2en"], _h["RAA_ratio_b2ed"], _s["RAA_ratio_b2e"]);
	DivideScatter2D(_s["RAA_b2e"], _s["RAA_c2e"], _s["RAA_ratio"]);
      }
      
      
    }
    void DivideScatter2D(Scatter2DPtr s1, Scatter2DPtr s2, Scatter2DPtr s)
    {
      for(unsigned int i =0; i<s2->numPoints(); i++){
	if(s2->point(i).y() == 0)
	  {
	    s->addPoint(s2->point(i).x(), std::numeric_limits<double>::quiet_NaN());
	    continue; 
	  }
	double yErr = (s1->point(i).y()/s2->point(i).y())*std::sqrt(std::pow(s1->point(i).yErrPlus()/s1->point(i).y(), 2) + std::pow(s2->point(i).yErrPlus()/s2->point(i).y(), 2));
	s->addPoint(s2->point(i).x(), s1->point(i).y()/s2->point(i).y(), s1->point(i).xErrPlus(), yErr);
      }    
    }
    
    /// @name Histograms
    //@}
    map<string, Histo1DPtr> _h; 
    map<string, Profile1DPtr> _p; 
    map<string, CounterPtr> _c; 
    map<string, Scatter2DPtr> _s;
    enum CollisionSystem {pp_200, AuAu200};
    CollisionSystem collSys;
    string beamOpt = "NONE";
    //@}

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2016_I1393529);

}
