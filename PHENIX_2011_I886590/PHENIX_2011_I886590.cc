// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2011_I886590 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2011_I886590);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fsPI(Cuts::abseta < 0.35 && Cuts::pT > 0.3*GeV && Cuts::pT < 3*GeV);
      declare(fsPI, "fsPI");

      const FinalState fsK(Cuts::abseta < 0.35 && Cuts::pT > 0.4*GeV && Cuts::pT < 2*GeV);
      declare(fsK, "fsK");

      const FinalState fsP(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV && Cuts::pT < 4.5*GeV);
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
      book(_h["xsec_piplus_624"], 5, 1, 1);
      book(_h["xsec_piminus_624"], 5, 1, 2);

      book(_h["xsec_kplus_624"], 6, 1, 1);
      book(_h["xsec_kminus_624"], 6, 1, 2);

      book(_h["xsec_p_noFD_624_1"], 7, 1, 1);
      book(_h["xsec_pbar_noFD_624_1"], 7, 1, 2);

      book(_h["xsec_p_withFD_624_1"], 8, 1, 1);
      book(_h["xsec_pbar_withFD_624_1"], 8, 1, 2);

      //book(_Nevt_after_cuts,"Nevt_after_cuts");
      //book(sow, "sow");
      book(_c["sow_pp200"], "sow_pp200");
      book(_c["sow_pp62"], "sow_pp62");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles fsPIParticles = applyProjection<FinalState>(event,"fsPI").particles();
      Particles fsKParticles = applyProjection<FinalState>(event,"fsK").particles();
      Particles fsPParticles = applyProjection<FinalState>(event,"fsP").particles();

      //_Nevt_after_cuts->fill();
      //sow->fill();
     
      // Pions
      for( const Particle& pPI : fsPIParticles)
      {
		
		const double pPIWeight = 1.0 / pPI.pt() / 2. / M_PI;
		// Fill histos 200GeV
		if (collsys == pp200)
		{
			_c["sow_pp200"]->fill();
      			if(pPI.pid() == 211) _h["xsec_piplus_200"]->fill(pPI.pT()/GeV, pPIWeight);
			if(pPI.pid() == -211) _h["xsec_piminus_200"]->fill(pPI.pT()/GeV, 1.0);
		}

		// Fill histos 62.4GeV
		if (collsys == pp62)
                {
			_c["sow_pp62"]->fill();
                	if(pPI.pid() == 211) _h["xsec_piplus_624"]->fill(pPI.pT()/GeV, 1.0);
                	if(pPI.pid() == -211) _h["xsec_piminus_624"]->fill(pPI.pT()/GeV, 1.0);
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
                	if(pK.pid() == 321) _h["xsec_kplus_624"]->fill(pK.pT()/GeV, 1.0);
                	if(pK.pid() == -321) _h["xsec_kminus_624"]->fill(pK.pT()/GeV, 1.0);
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
                	if(pP.pid() == 2212) _h["xsec_p_noFD_624_1"]->fill(pP.pT()/GeV, 1.0);
                	if(pP.pid() == -2212) _h["xsec_pbar_noFD_624_1"]->fill(pP.pT()/GeV, 1.0);

                	if(pP.pid() == 2212) _h["xsec_p_withFD_624_1"]->fill(pP.pT()/GeV, 1.0);
                	if(pP.pid() == -2212) _h["xsec_pbar_withFD_624_1"]->fill(pP.pT()/GeV, 1.0);
		}

      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h["xsec_piplus_200"]); // normalize to unity
      //scale(_h["xsec_piplus_200"], 1.0/ *_Nevt_after_cuts);
      //const double s = 1./sow->sumW();
      //scale(_h["xsec_piplus_200"], s);
      //normalize(_h["xsec_piplus_200"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
      //scale(_h["xsec_piplus_200"], crossSection()/picobarn/sow->sumW()); // norm to generated cross-section in pb (after cuts)

      if (collsys == pp200)
      {
		_h["xsec_piplus_200"]->scaleW(1.0/_c["sow_pp200"]->sumW());
      }

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
