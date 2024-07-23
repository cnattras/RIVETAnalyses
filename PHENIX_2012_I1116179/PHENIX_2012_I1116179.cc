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
  class PHENIX_2012_I1116179 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2012_I1116179);

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

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
    
    beamOpt = getOption<string>("beam", "NONE");
    const ParticlePair& beam = beams();
          
          if (beamOpt == "NONE") {
          if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
          {
            float NN = 197.;
            if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = AuAu;
          }
          }
          
          if (beamOpt == "AUAU200"){collSys = AuAu;}

	  declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");
    
	  //const FinalState fs(Cuts::abseta < 0.5 && Cuts::pT > 0.15*GeV);
	  //declare(fs, "fs");

	  const PromptFinalState pfs(Cuts::abseta < 0.35 && Cuts::pid == 22);
		declare(pfs, "pfs");

//	  const UnstableParticles up(Cuts::abseta < 0.5 && Cuts::abspid == 421);
//	  declare(up, "up");

	  book(_h["dir_photon_AuAu0005"], 1, 1, 1);
    book(_h["dir_photon_AuAu0510"], 1, 1, 2);
    book(_h["dir_photon_AuAu1015"], 1, 1, 3);
    book(_h["dir_photon_AuAu1520"], 1, 1, 4);
    book(_h["dir_photon_AuAu2030"], 1, 1, 5);
    book(_h["dir_photon_AuAu3040"], 1, 1, 6);
    book(_h["dir_photon_AuAu4050"], 1, 1, 7);
    book(_h["dir_photon_AuAu5060"], 1, 1, 8);
    book(_h["dir_photon_AuAu6092"], 1, 1, 9);
    book(_h["dir_photon_AuAu0092"], 1, 1, 10);

    //book(_h["dir_photon_pp"], 1, 1, 10);
    
	  book(_c["sow_AuAu0005"], "sow_AuAu0005");
	  book(_c["sow_AuAu0510"], "sow_AuAu0510");
	  book(_c["sow_AuAu1015"], "sow_AuAu1015");
	  book(_c["sow_AuAu1520"], "sow_AuAu1520");
	  book(_c["sow_AuAu2030"], "sow_AuAu2030");
	  book(_c["sow_AuAu3040"], "sow_AuAu3040");
	  book(_c["sow_AuAu4050"], "sow_AuAu4050");
	  book(_c["sow_AuAu5060"], "sow_AuAu5060");
	  book(_c["sow_AuAu6092"], "sow_AuAu6092");
	  book(_c["sow_AuAu0092"], "sow_AuAu0092");
	  //book(_c["sow_pp"], "sow_pp");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      const double c = cent();

      const PromptFinalState pfs = apply<PromptFinalState>(event, "pfs");
      const Particles photons = pfs.particles(Cuts::pT > 4.*GeV && Cuts::pT < 22.*GeV);


	  //cout << "c=" << c << endl;
    //for (auto& w : event.weights()) cout << " " << w << endl;

    if(collSys==AuAu)
    {

      if(c > 0. && c < 5.) 
      {
        _c["sow_AuAu0005"]->fill();
      }
      else if(c < 10.) 
      {
        _c["sow_AuAu0510"]->fill();
      }
      else if(c < 15.) 
      {
        _c["sow_AuAu1015"]->fill();
      }
      else if(c < 20.) 
      {
        _c["sow_AuAu1520"]->fill();
      }
      else if(c < 30.) 
      {
        _c["sow_AuAu2030"]->fill();
      }
      else if(c < 40.) 
      {
        _c["sow_AuAu3040"]->fill();
      }
      else if(c < 50.) 
      {
        _c["sow_AuAu4050"]->fill();
      }
      else if(c < 60.) 
      {
        _c["sow_AuAu5060"]->fill();
      }
      else if(c < 92.) 
      {
        _c["sow_AuAu6092"]->fill();
      }
      
      if(c > 0. && c < 92.)
        _c["sow_AuAu0092"]->fill();

      

      for(const Particle& p : photons)
      {
        if(c >= 0. && c < 5.) _h["dir_photon_AuAu0005"]->fill(p.pT()/GeV);
        else if (c < 10.) _h["dir_photon_AuAu0510"]->fill(p.pT()/GeV);
        else if (c < 15.) _h["dir_photon_AuAu1015"]->fill(p.pT()/GeV);
        else if (c < 20.) _h["dir_photon_AuAu1520"]->fill(p.pT()/GeV);
        else if (c < 30.) _h["dir_photon_AuAu2030"]->fill(p.pT()/GeV);
        else if (c < 40.) _h["dir_photon_AuAu3040"]->fill(p.pT()/GeV);
        else if (c < 50.) _h["dir_photon_AuAu4050"]->fill(p.pT()/GeV);
        else if (c < 60.) _h["dir_photon_AuAu5060"]->fill(p.pT()/GeV);
        else if (c < 92.) _h["dir_photon_AuAu6092"]->fill(p.pT()/GeV);
        
        if(c > 0. && c < 92.) 
          _h["dir_photon_AuAu0092"]->fill(p.pT()/GeV);

      }

    }
/*
    else if(beam="PP")
    {
      _c["sow_pp"]->fill();
      for(const Particle& p : pfsParticles) 
      {
        if(p.pid() == 22) _h["dir_photon_pp"]->fill(p.pT()/GeV);
      }

    }
      
*/
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      
      binShift(*_h["dir_photon_AuAu0005"]);
      binShift(*_h["dir_photon_AuAu0510"]);
      binShift(*_h["dir_photon_AuAu1015"]);
      binShift(*_h["dir_photon_AuAu1520"]);
      binShift(*_h["dir_photon_AuAu2030"]);
      binShift(*_h["dir_photon_AuAu3040"]);
      binShift(*_h["dir_photon_AuAu4050"]);
      binShift(*_h["dir_photon_AuAu5060"]);
      binShift(*_h["dir_photon_AuAu6092"]);
      binShift(*_h["dir_photon_AuAu0092"]);

      _h["dir_photon_AuAu0005"]->scaleW(1./_c["sow_AuAu0005"]->sumW());
      _h["dir_photon_AuAu0510"]->scaleW(1./_c["sow_AuAu0510"]->sumW());
      _h["dir_photon_AuAu1015"]->scaleW(1./_c["sow_AuAu1015"]->sumW());
      _h["dir_photon_AuAu1520"]->scaleW(1./_c["sow_AuAu1520"]->sumW());
      _h["dir_photon_AuAu2030"]->scaleW(1./_c["sow_AuAu2030"]->sumW());
      _h["dir_photon_AuAu3040"]->scaleW(1./_c["sow_AuAu3040"]->sumW());
      _h["dir_photon_AuAu4050"]->scaleW(1./_c["sow_AuAu4050"]->sumW());
      _h["dir_photon_AuAu5060"]->scaleW(1./_c["sow_AuAu5060"]->sumW());
      _h["dir_photon_AuAu6092"]->scaleW(1./_c["sow_AuAu6092"]->sumW());
      _h["dir_photon_AuAu0092"]->scaleW(1./_c["sow_AuAu0092"]->sumW());

      
    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
	  map<string, Scatter2DPtr> _s;
    string beamOpt;
    enum CollisionSystem {AuAu};
    CollisionSystem collSys;
    //@}


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2012_I1116179);

}
