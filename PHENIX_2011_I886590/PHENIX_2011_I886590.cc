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

    // Function to get bin centers for scaling
    bool getDeltaPt(YODA::Histo1D hist, double pT, double &deltaPt)
    {
    	//cout << "pT: " << pT << endl;
        if(pT > hist.xMin() && pT < hist.xMax())
        {
        	deltaPt = hist.bin(hist.binIndexAt(pT)).xMid();
                //cout << "DeltapT: " << deltaPt << endl;
                return true;
        }
        else return false;
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      //cout << "Made it into init\n";

      // Initialise and register projections
      const ALICE::PrimaryParticles fsPI(Cuts::abseta < 0.35 && Cuts::pT > 0.3*GeV && Cuts::pT < 3*GeV);
      declare(fsPI, "fsPI");

      const ALICE::PrimaryParticles fsK(Cuts::abseta < 0.35 && Cuts::pT > 0.4*GeV && Cuts::pT < 2*GeV);
      declare(fsK, "fsK");

      const ALICE::PrimaryParticles fsP(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV && Cuts::pT < 4.5*GeV);
      declare(fsP, "fsP");

      beamOpt = getOption<string>("beam", "NONE");
      if (beamOpt == "pp200") collsys = pp200;
      if (beamOpt == "pp62") collsys = pp62;

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

      book(_h["xsec_p_noFD_62_1"], 7, 1, 1);
      book(_h["xsec_pbar_noFD_62_1"], 7, 1, 2);

      book(_h["xsec_p_withFD_62_1"], 8, 1, 1);
      book(_h["xsec_pbar_withFD_62_1"], 8, 1, 2);

      // Counter histos for event weights
      book(_c["sow_pp200"], "sow_pp200");
      book(_c["sow_pp62"], "sow_pp62");

      //cout << "Made it past init\n";

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      //cout << "Made it into analyze\n";

      Particles fsPIParticles = applyProjection<ALICE::PrimaryParticles>(event,"fsPI").particles();
      Particles fsKParticles = applyProjection<ALICE::PrimaryParticles>(event,"fsK").particles();
      Particles fsPParticles = applyProjection<ALICE::PrimaryParticles>(event,"fsP").particles();
     
      // Pions
      for( const Particle& pPI : fsPIParticles)
      {
		
		//const double pPIWeight = 1.0 / pPI.pt() / 2. / M_PI;
		double partPt = pPI.pT()/GeV;
		double pt_weight = 1./(2.*M_PI*0.7);
		double deltaPt = 0;
		// Fill histos 200GeV
		if (collsys == pp200)
		{
			if(getDeltaPt(*_h["xsec_piplus_200"],partPt,deltaPt))
			{
			pt_weight /= deltaPt;
			_c["sow_pp200"]->fill();
      			if(pPI.pid() == 211) _h["xsec_piplus_200"]->fill(pPI.pT()/GeV, pt_weight);
			if(pPI.pid() == -211) _h["xsec_piminus_200"]->fill(pPI.pT()/GeV, 1.0);
			}
			//if(pPI.pid() == 211) _h["piplus_200_temp"]->fill(pPI.pT()/GeV);
		}

		// Fill histos 62.4GeV
		if (collsys == pp62)
                {
			_c["sow_pp62"]->fill();
                	if(pPI.pid() == 211) _h["xsec_piplus_62"]->fill(pPI.pT()/GeV, pt_weight);
                	if(pPI.pid() == -211) _h["xsec_piminus_62"]->fill(pPI.pT()/GeV, 1.0);
			//if(pPI.pid() == 211) _h["piplus_62_temp"]->fill(pPI.pT()/GeV);
		}

      }

      // Kaons
      for( const Particle& pK : fsKParticles)
      {
	
		// Fill histos 200GeV
		if (collsys == pp200)
		{
                	if(pK.pid() == 321) _h["xsec_kplus_200"]->fill(pK.pT()/GeV, 1.0);
                	if(pK.pid() == -321) _h["xsec_kminus_200"]->fill(pK.pT()/GeV, 1.0);
		}

		// Fill histos 62.4GeV
		if (collsys == pp62)
		{
                	if(pK.pid() == 321) _h["xsec_kplus_62"]->fill(pK.pT()/GeV, 1.0);
                	if(pK.pid() == -321) _h["xsec_kminus_62"]->fill(pK.pT()/GeV, 1.0);
		}

      }

      // Protons and antiprotons
      for( const Particle& pP : fsPParticles)
      {

		// Fill histos 200GeV
		if (collsys == pp200)
		{
                	if(pP.pid() == 2212) _h["xsec_p_noFD_200_1"]->fill(pP.pT()/GeV, 1.0);
                	if(pP.pid() == 2212) _h["xsec_p_noFD_200_2"]->fill(pP.pT()/GeV, 1.0);
                	if(pP.pid() == -2212) _h["xsec_pbar_noFD_200_1"]->fill(pP.pT()/GeV, 1.0);
                	if(pP.pid() == -2212) _h["xsec_pbar_noFD_200_2"]->fill(pP.pT()/GeV, 1.0);

                	if(pP.pid() == 2212) _h["xsec_p_withFD_200_1"]->fill(pP.pT()/GeV, 1.0);
                	if(pP.pid() == 2212) _h["xsec_p_withFD_200_2"]->fill(pP.pT()/GeV, 1.0);
                	if(pP.pid() == -2212) _h["xsec_pbar_withFD_200_1"]->fill(pP.pT()/GeV, 1.0);
                	if(pP.pid() == -2212) _h["xsec_pbar_withFD_200_2"]->fill(pP.pT()/GeV, 1.0);
		}

		// Fill histos 62.4GeV
		if (collsys == pp62)
		{
                	if(pP.pid() == 2212) _h["xsec_p_noFD_62_1"]->fill(pP.pT()/GeV, 1.0);
                	if(pP.pid() == -2212) _h["xsec_pbar_noFD_62_1"]->fill(pP.pT()/GeV, 1.0);

                	if(pP.pid() == 2212) _h["xsec_p_withFD_62_1"]->fill(pP.pT()/GeV, 1.0);
                	if(pP.pid() == -2212) _h["xsec_pbar_withFD_62_1"]->fill(pP.pT()/GeV, 1.0);
		}

      }

      //cout << "Made it past analyze\n";

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //cout << "Made it into finalize\n";

      // Method 1
      //normalize(_h["xsec_piplus_200"]); // normalize to unity
      
      // Method 2
      //scale(_h["xsec_piplus_200"], 1.0/ *_Nevt_after_cuts);

      // Method 3
      //const double s = 1./sow->sumW();
      //scale(_h["xsec_piplus_200"], s);

      // Method 4
      //normalize(_h["xsec_piplus_200"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)

      // Method 5
      //scale(_h["xsec_piplus_200"], crossSection()/picobarn/sow->sumW()); // norm to generated cross-section in pb (after cuts)

      // Method 6
      if (collsys == pp200)
      {
		_h["xsec_piplus_200"]->scaleW(crossSection()*1.E-9/_c["sow_pp200"]->sumW());
		//_h["piplus_200_temp"]->scaleW(1.0/_c["sow_pp200"]->sumW());
      }

      if (collsys == pp62)
      {
		_h["xsec_piplus_62"]->scaleW(crossSection()*1.E-9/_c["sow_pp62"]->sumW());
		//_h["piplus_62_temp"]->scaleW(1.0/_c["sow_pp62"]->sumW());
      }

      //divide(_h["_xec_piplus_200"], _h["_xsec_piplus_62"], _s["piplus_ratio"]);

      //cout << "Made it past finalize\n";

    }

    //@}

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;
    string beamOpt = "";
    enum CollisionSystem {pp200, pp62};
    CollisionSystem collsys;
    //CounterPtr _Nevt_after_cuts;
    //CounterPtr sow;
    //@}

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2011_I886590);

}
