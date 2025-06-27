// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
//#include "Centralities/RHICCentrality.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2016_I1394433 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2016_I1394433);


    /// @name Analysis methods
    //@{

    void init() {

        beamOpt = getOption<string>("beam", "NONE");

        const ParticlePair& beam = beams();
        int NN = 0;

        if (beamOpt == "NONE") {
        
            if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
          NN = 197.;
          if (fuzzyEquals(sqrtS()/GeV, 7.7*
          NN, 5)) collSys = AuAu7;
          if (fuzzyEquals(sqrtS()/GeV, 14.5*NN, 5)) collSys = AuAu14;
          if (fuzzyEquals(sqrtS()/GeV, 19.6*NN, 5)) collSys = AuAu19;
          if (fuzzyEquals(sqrtS()/GeV, 27*NN, 5)) collSys = AuAu27;
          if (fuzzyEquals(sqrtS()/GeV, 39*NN, 5)) collSys = AuAu39;
          if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 5)) collSys = AuAu62;
          if (fuzzyEquals(sqrtS()/GeV, 130*NN, 5)) collSys = AuAu130;
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = AuAu200;
      }
      if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000290630)
      {
          NN = 63.;
          if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 5)) collSys = CuCu62;
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = CuCu200;
      }
      if ((beam.first.pid() == 1000290630 && beam.second.pid() == 1000791970) || (beam.first.pid() == 1000791970 && beam.second.pid() == 1000290630) )
      {
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = CuAu200;
      }
      if (beam.first.pid() == 1000922380 && beam.second.pid() == 1000922380)
      {
        NN = 238.;
          if (fuzzyEquals(sqrtS()/GeV, 193*NN, 5)) collSys = UU193;
      }
      if ((beam.first.pid() == 1000010020 && beam.second.pid() == 1000791970) || (beam.first.pid() == 1000791970 && beam.second.pid() == 1000010020))
      {
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = dAu200;
      }
      if ((beam.first.pid() == 1000020030 && beam.second.pid() == 1000791970) || (beam.first.pid() == 1000791970 && beam.second.pid() == 1000020030))
      {
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = HeAu200;
      }
          
        }
    
    if (beamOpt == "AUAU200") collSys = AuAu200;
    else if (beamOpt == "AUAU130") collSys = AuAu130;
    else if (beamOpt == "AUAU62") collSys = AuAu62;
    else if (beamOpt == "AUAU39") collSys = AuAu39;
    else if (beamOpt == "AUAU27") collSys = AuAu27;
    else if (beamOpt == "AUAU19") collSys = AuAu19;
    else if (beamOpt == "AUAU14") collSys = AuAu14;
    else if (beamOpt == "AUAU7") collSys = AuAu7;
    else if (beamOpt == "CUCU200") collSys = CuCu200;
    else if (beamOpt == "CUCU62") collSys = CuCu62;
    else if (beamOpt == "UU193") collSys = UU193;
    else if (beamOpt == "DAU200") collSys = dAu200;
    else if (beamOpt == "HEAU200") collSys = HeAu200;
    else if (beamOpt == "CUAU200") collSys = CuAu200;


    declareCentrality(RHICCentrality("PHENIX"),"RHIC_2019_CentralityCalibration:exp=PHENIX","CMULT","CMULT");
        
     //Au+Au collisions

        book(_hist_AuAu_E_200, "d01-x01-y03", refData(1, 1, 3));
        book(_hist_AuAu_Ch_200, "d02-x01-y03", refData(2, 1, 3));
        book(_hist_AuAu_E_130, "d03-x01-y03", refData (3, 1, 3));
        book(_hist_AuAu_Ch_130, "d04-x01-y03", refData (4, 1, 3));
        book(_hist_AuAu_E_62, "d05-x01-y03", refData( 5, 1, 3));
        book(_hist_AuAu_Ch_62, "d06-x01-y03", refData (6, 1, 3));
        book(_hist_AuAu_E_39, "d07-x01-y03", refData (7, 1, 3));
        book(_hist_AuAu_Ch_39, "d08-x01-y03", refData (8, 1, 3));
        book(_hist_AuAu_E_27, "d09-x01-y03", refData(9,1,3));
        book(_hist_AuAu_Ch_27, "d10-x01-y03", refData(10,1,3));
        book(_hist_AuAu_E_19, "d11-x01-y03", refData(11,1,3));
        book(_hist_AuAu_Ch_19, "d12-x01-y03", refData(12,1,3));
        book(_hist_AuAu_E_14, "d13-x01-y03", refData(13,1,3));
        book(_hist_AuAu_Ch_14, "d14-x01-y03", refData(14,1,3));
        book(_hist_AuAu_E_7, "d15-x01-y03", refData(15,1,3));
        book(_hist_AuAu_Ch_7, "d16-x01-y03", refData(16,1,3));
      
     //Cu+Cu collisions   
      
        book(_hist_CuCu_E_200, "d17-x01-y03", refData(17,1,3));
        book(_hist_CuCu_Ch_200, "d18-x01-y03", refData(18,1,3));
        book(_hist_CuCu_E_62, "d19-x01-y03", refData(19,1,3));
        book(_hist_CuCu_Ch_62, "d20-x01-y03", refData(20,1,3));

    //Cu+Au collisions

        book(_hist_CuAu_E_200, "d21-x01-y03", refData(21,1,3));
        book(_hist_CuAu_Ch_200, "d22-x01-y03", refData(22,1,3));

    //U+U collisions

        book(_hist_UU_E_193, "d23-x01-y03", refData(23,1,3));
        book(_hist_UU_Ch_193, "d24-x01-y03", refData(24,1,3));

    //d+Au collisions

        book(_hist_dAu_E_200, "d25-x01-y03", refData(25,1,3));
        book(_hist_dAu_Ch_200, "d26-x01-y03", refData(26,1,3));

    //He+Au collisions
    
        book(_hist_HeAu_E_200, "d27-x01-y03", refData(27,1,3));
        book(_hist_HeAu_Ch_200, "d28-x01-y03", refData(28,1,3));

    // Final State projection
      const FinalState fs(Cuts::abseta < 0.35);
      declare(fs, "fs");
      
        
     declare(ALICE::PrimaryParticles(Cuts::abseta < 0.35 && Cuts::pT > 0.0*MeV && Cuts::abscharge > 0), "APRIM");
    }


    void analyze(const Event& event)
    {
        int nchcounter = 0;  
        double sumET = 0.; 
        
        //get charged particles
        Particles chargedParticles = apply<ALICE::PrimaryParticles>(event,"APRIM").particles();

        //get fs particles
        Particles fsParticles = apply<FinalState>(event,"fs").particles();
        
        //get centrality info
        const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
        const double c = cent();
        if (c > 60) vetoEvent;
        
        //loop over particles to get total charged particles & Et
        nchcounter = chargedParticles.size();

        for(const Particle& part : fsParticles)
        {
            sumET += part.Et()/GeV;
        }
        
        
        //fill hists for AuAU collisions
        if(collSys == AuAu200)
        {
            _hist_AuAu_Ch_200->fill(c,nchcounter/0.7);
            _hist_AuAu_E_200->fill(c,sumET/0.7);
        }
        else if(collSys == AuAu130)
        {
            _hist_AuAu_Ch_130->fill(c,nchcounter/0.7);
            _hist_AuAu_E_130->fill(c,sumET/0.7);
        }
        else if(collSys == AuAu62)
        {
            _hist_AuAu_Ch_62->fill(c,nchcounter/0.7);
            _hist_AuAu_E_62->fill(c,sumET/0.7);
        }
        else if(collSys == AuAu39)
        {
            _hist_AuAu_Ch_39->fill(c,nchcounter/0.7);
            _hist_AuAu_E_39->fill(c,sumET/0.7);
        }
        else if(collSys == AuAu27)
        {
            _hist_AuAu_Ch_27->fill(c,nchcounter/0.7);
            _hist_AuAu_E_27->fill(c,sumET/0.7);
        }
        else if(collSys == AuAu19)
        {
            _hist_AuAu_Ch_19->fill(c,nchcounter/0.7);
            _hist_AuAu_E_19->fill(c,sumET/0.7);
        }
        else if(collSys == AuAu14)
        {
            _hist_AuAu_Ch_14->fill(c,nchcounter/0.7);
            _hist_AuAu_E_14->fill(c,sumET/0.7);
        }
        else if(collSys == AuAu7)
        {
            _hist_AuAu_Ch_7->fill(c,nchcounter/0.7);
            _hist_AuAu_E_7->fill(c,sumET/0.7);
        }

        //fill hists for CuCu collisions

        if(collSys == CuCu200)
        {
            _hist_CuCu_Ch_200->fill(c,nchcounter/0.7);
            _hist_CuCu_E_200->fill(c,sumET/0.7);
        }
        else if(collSys == CuCu62)
        {
            _hist_CuCu_Ch_62->fill(c,nchcounter/0.7);
            _hist_CuCu_E_62->fill(c,sumET/0.7);
        }

         //fill hists for CuAu collisions

         if(collSys == CuAu200)
        {
            _hist_CuAu_Ch_200->fill(c,nchcounter/0.7);
            _hist_CuAu_E_200->fill(c,sumET/0.7);
        }

        //fill hists for UU collisions

        if(collSys == UU193)
        {
            _hist_UU_Ch_193->fill(c,nchcounter/0.7);
            _hist_UU_E_193->fill(c,sumET/0.7);
        }

        //fill hists for dAu collisions

        if(collSys == dAu200)
        {
            _hist_dAu_Ch_200->fill(c,nchcounter/0.7);
            _hist_dAu_E_200->fill(c,sumET/0.7);
        }

        //fill hists for HeAu collisions

        if(collSys == HeAu200)
        {
            _hist_HeAu_Ch_200->fill(c,nchcounter/0.7);
            _hist_HeAu_E_200->fill(c,sumET/0.7);
        }
    }

    void finalize()
    {

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;

    //AuAu collision
    Profile1DPtr _hist_AuAu_Ch_200;
    Profile1DPtr _hist_AuAu_Ch_130;
    Profile1DPtr _hist_AuAu_Ch_62;
    Profile1DPtr _hist_AuAu_Ch_39;
    Profile1DPtr _hist_AuAu_Ch_27;
    Profile1DPtr _hist_AuAu_Ch_19;
    Profile1DPtr _hist_AuAu_Ch_14;
    Profile1DPtr _hist_AuAu_Ch_7;
    Profile1DPtr _hist_AuAu_E_200;
    Profile1DPtr _hist_AuAu_E_130;
    Profile1DPtr _hist_AuAu_E_62;
    Profile1DPtr _hist_AuAu_E_39;
    Profile1DPtr _hist_AuAu_E_27;
    Profile1DPtr _hist_AuAu_E_19;
    Profile1DPtr _hist_AuAu_E_14;
    Profile1DPtr _hist_AuAu_E_7;

    //CuCu collision
    Profile1DPtr _hist_CuCu_Ch_200;
    Profile1DPtr _hist_CuCu_Ch_62;
    Profile1DPtr _hist_CuCu_E_200;
    Profile1DPtr _hist_CuCu_E_62;

    //CuAu collision
    Profile1DPtr _hist_CuAu_Ch_200;
    Profile1DPtr _hist_CuAu_E_200;

    //UU collision
    Profile1DPtr _hist_UU_Ch_193;
    Profile1DPtr _hist_UU_E_193;

    //dAu collision
    Profile1DPtr _hist_dAu_Ch_200;
    Profile1DPtr _hist_dAu_E_200;

    //HeAu collision
    Profile1DPtr _hist_HeAu_Ch_200;
    Profile1DPtr _hist_HeAu_E_200;
    //@}
    string beamOpt;
    enum CollisionSystem {NONE, AuAu7,AuAu14,AuAu19,AuAu27,AuAu39,AuAu62,AuAu130,AuAu200,CuCu62,CuCu200,CuAu200,UU193,dAu200,HeAu200,pp200};
    CollisionSystem collSys;


  };


  RIVET_DECLARE_PLUGIN(PHENIX_2016_I1394433);

}
