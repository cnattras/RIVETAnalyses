// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/AliceCommon.hh"

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Projections/AliceCommon.hh"

#include "Centrality/RHICCentrality.hh" //external header for Centrality calculation
#include <math.h>
#define _USE_MATH_DEFINES
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {


  class PHENIX_2008_I778168 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    PHENIX_2008_I778168()
      : Analysis("PHENIX_2008_I778168")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
	
	 std::initializer_list<int> pdgIds = {111};  // Pion 0

      const PrimaryParticles fs(pdgIds, Cuts::abseta < 0.35 && Cuts::abscharge == 0);
      declare(fs, "fs");

      beamOpt = getOption<string>("beam","NONE");


      if(beamOpt=="PP200") collSys = pp200;
      else if(beamOpt=="AUAU") collSys = AuAu;
	
	if(collSys != pp200) declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      /// @todo Initialise and register projections here

      /// @todo Book histograms here, e.g.:
      // _h_XXXX = bookProfile1D(1, 1, 1);
      // _h_YYYY = bookHistogram1D(2, 1, 1);

 //Yields_________________
      book(hPion0Pt["ptyieldsc05"], 1, 1, 3);
      book(hPion0Pt["ptyieldsc10"], 1, 1, 1);
      book(hPion0Pt["ptyieldsc20"], 1, 1, 2);
      book(hPion0Pt["ptyieldsc30"], 2, 1, 1);
      book(hPion0Pt["ptyieldsc40"], 2, 1, 2);
      book(hPion0Pt["ptyieldsc50"], 3, 1, 1);
      book(hPion0Pt["ptyieldsc60"], 3, 1, 2);
      book(hPion0Pt["ptyieldsc70"], 4, 1, 1);
      book(hPion0Pt["ptyieldsc80"], 5, 1, 1);
      book(hPion0Pt["ptyieldsc92"], 5, 1, 2);
      book(hPion0Pt["ptyieldscall"], 1, 1, 3);

      book(sow["sow_c05"],"sow_c05");
      book(sow["sow_c10"],"sow_c10");
      book(sow["sow_c20"],"sow_c20");
      book(sow["sow_c30"],"sow_c30");
      book(sow["sow_c40"],"sow_c40");
      book(sow["sow_c50"],"sow_c50");
      book(sow["sow_c60"],"sow_c60");
      book(sow["sow_c70"],"sow_c70");
      book(sow["sow_c80"],"sow_c80");
      book(sow["sow_c92"],"sow_c92");
      book(sow["sow_call"],"sow_call");

//RAA _______________________________
      string refnameRaa1 = mkAxisCode(6,1,3);
            const Scatter2D& refdataRaa1 =refData(refnameRaa1);
      book(hPion0Pt["c5Pt_AuAu"], refnameRaa1 + "_AuAu", refdataRaa1);
      book(hPion0Pt["c5Pt_pp"], refnameRaa1 + "_pp", refdataRaa1);
      book(hRaa["Raa_c05_AuAu"], refnameRaa1);

      string refnameRaa2 = mkAxisCode(6,1,1);
            const Scatter2D& refdataRaa2 =refData(refnameRaa2);
      book(hPion0Pt["c10Pt_AuAu"], refnameRaa2 + "_AuAu", refdataRaa2);
      book(hPion0Pt["c10Pt_pp"], refnameRaa2 + "_pp", refdataRaa2);
      book(hRaa["Raa_c010_AuAu"], refnameRaa2);

      string refnameRaa3 = mkAxisCode(6,1,1);
            const Scatter2D& refdataRaa3 =refData(refnameRaa3);
      book(hPion0Pt["c40Pt_AuAu"], refnameRaa3 + "_AuAu", refdataRaa3);
      book(hPion0Pt["c40Pt_pp"], refnameRaa3 + "_pp", refdataRaa3);
      book(hRaa["Raa_c2040_AuAu"], refnameRaa3);

      string refnameRaa4 = mkAxisCode(5,1,3);
            const Scatter2D& refdataRaa4 =refData(refnameRaa4);
      book(hPion0Pt["c60Pt_AuAu39"], refnameRaa4 + "_AuAu39", refdataRaa4);
      book(hPion0Pt["c60Pt39_pp"], refnameRaa4 + "_pp39", refdataRaa4);
      book(hRaa["Raa_c4060_AuAu39"], refnameRaa4);

      string refnameRaa5 = mkAxisCode(5,1,4);
            const Scatter2D& refdataRaa5 =refData(refnameRaa5);
      book(hPion0Pt["c86Pt_AuAu39"], refnameRaa5 + "_AuAu39", refdataRaa5);
      book(hPion0Pt["c86Pt39_pp"], refnameRaa5 + "_pp39", refdataRaa5);
      book(hRaa["Raa_c6086_AuAu39"], refnameRaa5);

      string refnameRaa6 = mkAxisCode(6,1,2);
            const Scatter2D& refdataRaa6 =refData(refnameRaa6);
      book(hPion0Pt["callPt_AuAu39"], refnameRaa6 + "_AuAu39", refdataRaa6);
      book(hPion0Pt["callPt39_pp"], refnameRaa6 + "_pp39", refdataRaa6);
      book(hRaa["Raa_minbias_AuAu39"], refnameRaa6);

      string refnameRaa7 = mkAxisCode(7,1,1);
            const Scatter2D& refdataRaa7 =refData(refnameRaa7);
      book(hPion0Pt["c10Pt_AuAu62"], refnameRaa7 + "_AuAu62", refdataRaa7);
      book(hPion0Pt["c10Pt62_pp"], refnameRaa7 + "_pp62", refdataRaa7);
      book(hRaa["Raa_c010_AuAu62"], refnameRaa7);
      
      string refnameRaa8 = mkAxisCode(7,1,2);
            const Scatter2D& refdataRaa8 =refData(refnameRaa8);
      book(hPion0Pt["c20Pt_AuAu62"], refnameRaa8 + "_AuAu62", refdataRaa8);
      book(hPion0Pt["c20Pt62_pp"], refnameRaa8 + "_pp62", refdataRaa8);
      book(hRaa["Raa_c1020_AuAu62"], refnameRaa8);

      string refnameRaa9 = mkAxisCode(7,1,3);
            const Scatter2D& refdataRaa9 =refData(refnameRaa9);
      book(hPion0Pt["c40Pt_AuAu62"], refnameRaa9 + "_AuAu62", refdataRaa9);
      book(hPion0Pt["c40Pt62_pp"], refnameRaa9 + "_pp62", refdataRaa9);
      book(hRaa["Raa_c2040_AuAu62"], refnameRaa9);

      string refnameRaa10 = mkAxisCode(7,1,4);
            const Scatter2D& refdataRaa10 =refData(refnameRaa10);
      book(hPion0Pt["c60Pt_AuAu62"], refnameRaa10 + "_AuAu62", refdataRaa10);
      book(hPion0Pt["c60Pt62_pp"], refnameRaa10 + "_pp62", refdataRaa10);
      book(hRaa["Raa_c4060_AuAu62"], refnameRaa10);

      string refnameRaa11 = mkAxisCode(8,1,1);
            const Scatter2D& refdataRaa11 =refData(refnameRaa11);
      book(hPion0Pt["c86Pt_AuAu62"], refnameRaa11 + "_AuAu62", refdataRaa11);
      book(hPion0Pt["c86Pt62_pp"], refnameRaa11 + "_pp62", refdataRaa11);
      book(hRaa["Raa_c6086_AuAu62"], refnameRaa11);

      string refnameRaa12 = mkAxisCode(7,1,5);
            const Scatter2D& refdataRaa12 =refData(refnameRaa12);
      book(hPion0Pt["callPt_AuAu62"], refnameRaa12 + "_AuAu62", refdataRaa12);
      book(hPion0Pt["callPt62_pp"], refnameRaa12 + "_pp62", refdataRaa12);
      book(hRaa["Raa_minbias_AuAu62"], refnameRaa12);

      book(sow["sow_pp39"],"sow_pp39");
      book(sow["sow_pp62"],"sow_pp62");


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// @todo Do the event by event analysis here

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here

      // scale(_h_YYYY, crossSection()/sumOfWeights()); # norm to cross section
      // normalize(_h_YYYY); # normalize to unity

    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    //@{
    AIDA::IProfile1D *_h_XXXX;
    AIDA::IHistogram1D *_h_YYYY;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PHENIX_2008_I778168);

}
