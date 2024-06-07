// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2011_I886590 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2011_I886590);

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

    // Function to get bin centers for scaling
    bool getDeltaPt(YODA::Histo1D hist, double pT, double &deltaPt)
    {
        if(pT > hist.xMin() && pT < hist.xMax())
        {
        	deltaPt = hist.bin(hist.binIndexAt(pT)).xMid();
                return true;
        }
        else return false;
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const ALICE::PrimaryParticles fsPI(Cuts::abseta < 0.35 && Cuts::pT > 0.3*GeV && Cuts::pT < 3*GeV);
      declare(fsPI, "fsPI");

      const ALICE::PrimaryParticles fsK(Cuts::abseta < 0.35 && Cuts::pT > 0.4*GeV && Cuts::pT < 2*GeV);
      declare(fsK, "fsK");

      const ALICE::PrimaryParticles fsP(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV && Cuts::pT < 4.5*GeV);
      declare(fsP, "fsP");

      // Book histograms
      // Histos from HEPdata at 200GeV
      book(_h["xsec_piplus_200"], 1, 1, 1);
      book(_h["xsec_piminus_200"], 1, 1, 2);

      book(_h["xsec_kplus_200"], 2, 1, 1);
      book(_h["xsec_kminus_200"], 2, 1, 2);

      book(_h["xsec_p_noFD_200_1"], 3, 1, 1);
      book(_h["xsec_p_noFD_200_2"], 9, 1, 1);
      book(_h["xsec_pbar_noFD_200_1"], 3, 1, 2);
      book(_h["xsec_pbar_noFD_200_2"], 9, 1, 2);

      book(_h["xsec_p_withFD_200_1"], 4, 1, 1);
      book(_h["xsec_p_withFD_200_2"], 10, 1, 1);
      book(_h["xsec_pbar_withFD_200_1"], 4, 1, 2);
      book(_h["xsec_pbar_withFD_200_2"], 10, 1, 2);

      // Histos from HEPdata at 62.4GeV
      book(_h["xsec_piplus_62"], 5, 1, 1);
      book(_h["xsec_piminus_62"], 5, 1, 2);

      book(_h["xsec_kplus_62"], 6, 1, 1);
      book(_h["xsec_kminus_62"], 6, 1, 2);

      book(_h["xsec_p_noFD_62"], 7, 1, 1);
      book(_h["xsec_pbar_noFD_62"], 7, 1, 2);

      book(_h["xsec_p_withFD_62"], 8, 1, 1);
      book(_h["xsec_pbar_withFD_62"], 8, 1, 2);

      // Counter histos for event weights
      book(_c["sow_pp200"], "_sow_pp200");
      book(_c["sow_pp62"], "_sow_pp62");

      //Cross sections
      book(_c["xsec_pp200"], "xsec_pp200");
      book(_c["xsec_pp62"], "xsec_pp62");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles fsPIParticles = applyProjection<ALICE::PrimaryParticles>(event,"fsPI").particles();
      Particles fsKParticles = applyProjection<ALICE::PrimaryParticles>(event,"fsK").particles();
      Particles fsPParticles = applyProjection<ALICE::PrimaryParticles>(event,"fsP").particles();

      pair<double,double> cs = HepMCUtils::crossSection(*event.genEvent());

       // Beam options
      beamOpt = getOption<string>("beam", "NONE");
      const ParticlePair& beam = beams();
      
      if(beamOpt == "NONE"){

      //check for proton proton
      if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
      {
  			int NN = 1.;
  			if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3))
        {
          collsys = pp200;
          _c["sow_pp200"]->fill();
          _c["xsec_pp200"]->fill(cs.first);
        }
  			if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 1E-3))
        {
          collsys = pp62;
          _c["sow_pp62"]->fill();
          _c["xsec_pp62"]->fill(cs.first);
        }
  		}else{
        std::cerr << "Error: Not PP beam" << std::endl;
      }
    }
     
      
      if (beamOpt == "PP200") collsys = pp200;
      if (beamOpt == "PP62") collsys = pp62;


      // Pions
      for( const Particle& pPI : fsPIParticles)
      {
	      	double PtPI = pPI.pT()/GeV;
                double PtPI_weight = 1./(2.*M_PI*0.7);
                double deltaPtPI = 0;

		// Fill histos 200GeV



		if (collsys == pp200)
		{


			if(getDeltaPt(*_h["xsec_piplus_200"],PtPI,deltaPtPI))
			{
      				if(pPI.pid() == 211)
              {
                PtPI_weight /= deltaPtPI;
                _h["xsec_piplus_200"]->fill(pPI.pT()/GeV, PtPI_weight);
              }
			}

			if(getDeltaPt(*_h["xsec_piminus_200"],PtPI,deltaPtPI))
			{
				if(pPI.pid() == -211)
        {
          PtPI_weight /= deltaPtPI;
          _h["xsec_piminus_200"]->fill(pPI.pT()/GeV, PtPI_weight);
        }
			}
		}

		// Fill histos 62.4GeV
		if (collsys == pp62)
                {

			if(getDeltaPt(*_h["xsec_piplus_62"],PtPI,deltaPtPI))
			{
                		if(pPI.pid() == 211)
                    {
                      PtPI_weight /= deltaPtPI;
                      _h["xsec_piplus_62"]->fill(pPI.pT()/GeV, PtPI_weight);
                    }
			}

			if(getDeltaPt(*_h["xsec_piminus_62"],PtPI,deltaPtPI))
			{

                		if(pPI.pid() == -211)
                    {
                      PtPI_weight /= deltaPtPI;
                      _h["xsec_piminus_62"]->fill(pPI.pT()/GeV, PtPI_weight);
                    }
			}
		}

      }

      // Kaons
      for( const Particle& pK : fsKParticles)
      {
                double PtK = pK.pT()/GeV;
                double PtK_weight = 1./(2.*M_PI*0.7);
                double deltaPtK = 0;

		// Fill histos 200GeV
		if (collsys == pp200)
		{

                        if(getDeltaPt(*_h["xsec_kplus_200"],PtK,deltaPtK))
                        {

                		if(pK.pid() == 321)
                    {
                      PtK_weight /= deltaPtK;
                      _h["xsec_kplus_200"]->fill(pK.pT()/GeV, PtK_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_kminus_200"],PtK,deltaPtK))
                        {

                		if(pK.pid() == -321)
                    {
                      PtK_weight /= deltaPtK;
                      _h["xsec_kminus_200"]->fill(pK.pT()/GeV, PtK_weight);
                    }
			}
		}

		// Fill histos 62.4GeV
		if (collsys == pp62)
		{

                        if(getDeltaPt(*_h["xsec_kplus_62"],PtK,deltaPtK))
                        {

                		if(pK.pid() == 321)
                    {
                      PtK_weight /= deltaPtK;
                      _h["xsec_kplus_62"]->fill(pK.pT()/GeV, PtK_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_kminus_62"],PtK,deltaPtK))
                        {

                		if(pK.pid() == -321)
                    {
                      PtK_weight /= deltaPtK;
                      _h["xsec_kminus_62"]->fill(pK.pT()/GeV, PtK_weight);
                    }
			}
		}

      }

      // Protons and antiprotons
      for( const Particle& pP : fsPParticles)
      {
                double PtP = pP.pT()/GeV;
                double PtP_weight = 1./(2.*M_PI*0.7);
                double deltaPtP = 0;

		// Fill histos 200GeV
		if (collsys == pp200)
		{

                        if(getDeltaPt(*_h["xsec_p_noFD_200_1"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_noFD_200_1"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_p_noFD_200_2"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_noFD_200_2"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_noFD_200_1"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_noFD_200_1"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_noFD_200_2"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_noFD_200_2"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_p_withFD_200_1"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_withFD_200_1"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_p_withFD_200_2"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_withFD_200_2"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_withFD_200_1"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_withFD_200_1"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_withFD_200_2"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_withFD_200_2"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}
		}

		// Fill histos 62.4GeV
		if (collsys == pp62)
		{

                        if(getDeltaPt(*_h["xsec_p_noFD_62"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_noFD_62"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_noFD_62"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_noFD_62"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_p_withFD_62"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_withFD_62"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_withFD_62"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_withFD_62"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}
		}

      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      

      binShift(*_h["xsec_piplus_200"]);
      binShift(*_h["xsec_piminus_200"]);
      binShift(*_h["xsec_kplus_200"]);
      binShift(*_h["xsec_kminus_200"]);
      binShift(*_h["xsec_p_noFD_200_1"]);
      binShift(*_h["xsec_p_noFD_200_2"]);
      binShift(*_h["xsec_pbar_noFD_200_1"]);
      binShift(*_h["xsec_pbar_noFD_200_2"]);
      binShift(*_h["xsec_p_withFD_200_1"]);
      binShift(*_h["xsec_p_withFD_200_2"]);
      binShift(*_h["xsec_pbar_withFD_200_1"]);
      binShift(*_h["xsec_pbar_withFD_200_2"]);
      binShift(*_h["xsec_piplus_62"]);
      binShift(*_h["xsec_piminus_62"]);
      binShift(*_h["xsec_kplus_62"]);
      binShift(*_h["xsec_kminus_62"]);
      binShift(*_h["xsec_p_noFD_62"]);
      binShift(*_h["xsec_pbar_noFD_62"]);
      binShift(*_h["xsec_p_withFD_62"]);
      binShift(*_h["xsec_pbar_withFD_62"]);

      if (_c["sow_pp200"]->sumW() != 0){
        double xsec200 = (_c["xsec_pp200"]->sumW()/millibarn)/_c["sow_pp200"]->sumW();
        _h["xsec_piplus_200"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
        _h["xsec_piminus_200"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
        _h["xsec_kplus_200"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
        _h["xsec_kminus_200"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
        _h["xsec_p_noFD_200_1"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
        _h["xsec_p_noFD_200_2"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
        _h["xsec_pbar_noFD_200_1"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
        _h["xsec_pbar_noFD_200_2"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
        _h["xsec_p_withFD_200_1"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
        _h["xsec_p_withFD_200_2"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
        _h["xsec_pbar_withFD_200_1"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
        _h["xsec_pbar_withFD_200_2"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
      }else{
        std::cerr << "Error: Divide by zero encountered." << std::endl;
      }
		  

      if (_c["sow_pp62"]->sumW() != 0){
        double xsec62 = (_c["xsec_pp62"]->sumW()/millibarn)/_c["sow_pp62"]->sumW();
        if(_c["sow_pp200"]->sumW() !=0){
          _h["xsec_piplus_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
          _h["xsec_piminus_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
          _h["xsec_kplus_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
          _h["xsec_kminus_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
          _h["xsec_p_noFD_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
          _h["xsec_pbar_noFD_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
          _h["xsec_p_withFD_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
          _h["xsec_pbar_withFD_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
        }else{
          std::cerr << "Error: Divide by zero encountered. Cannot Scale Histogram" << std::endl;
        }
      }else{
        std::cerr << "Error: Divide by zero encountered. Cannot define xsec62" << std::endl;
      }
        
		  


    }

    //@}

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;
    string beamOpt;
    enum CollisionSystem {pp200, pp62};
    CollisionSystem collsys;
    //@}

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2011_I886590);

}
