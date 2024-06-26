// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
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
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2016_I1394433);


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
          if (fuzzyEquals(sqrtS()/GeV, 7.7*NN, 1E-3)) collSys = AuAu7;
          if (fuzzyEquals(sqrtS()/GeV, 14.5*NN, 1E-3)) collSys = AuAu14;
          if (fuzzyEquals(sqrtS()/GeV, 19.6*NN, 1E-3)) collSys = AuAu19;
          if (fuzzyEquals(sqrtS()/GeV, 27*NN, 1E-3)) collSys = AuAu27;
          if (fuzzyEquals(sqrtS()/GeV, 39*NN, 1E-3)) collSys = AuAu39;
          if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 1E-3)) collSys = AuAu62;
          if (fuzzyEquals(sqrtS()/GeV, 130*NN, 1E-3)) collSys = AuAu130;
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSys = AuAu200;
      }
      if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000290630)
      {
          NN = 63.;
          if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 1E-3)) collSys = CuCu62;
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSys = CuCu200;
      }
      if (beam.first.pid() == 1000290630 && beam.second.pid() == 1000791970)
      {
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSys = CuAu200;
      }
      if (beam.first.pid() == 1000922380 && beam.second.pid() == 1000922380)
      {
          if (fuzzyEquals(sqrtS()/GeV, 193*NN, 1E-3)) collSys = UU193;
      }
      if (beam.first.pid() == 1000010020 && beam.second.pid() == 1000791970)
      {
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSys = DAu200;
      }
      if (beam.first.pid() == 1000020030 && beam.second.pid() == 1000791970)
      {
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSys = HeAu200;
      }
      if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
      {
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) collSys = pp200;
      }    
        }

    declareCentrality(RHICCentrality("PHENIX"),"RHIC_2019_CentralityCalibration:exp=PHENIX","CMULT","CMULT");
        
     //Au+Au collisions

        book(_h["hist_E_200"], 1, 1, 3);
        book(_h["hist_Ch_200"], 2, 1, 3);
        book(_h["hist_E_130"], 3, 1, 3);
        book(_h["hist_Ch_130"], 4, 1, 3);
        book(_h["hist_E_62.4"], 5, 1, 3);
        book(_h["hist_Ch_62.4"], 6, 1, 3);
        book(_h["hist_E_39"], 7, 1, 3);
        book(_h["hist_Ch_39"], 8, 1, 3);
        book(_h["hist_E_27"],9,1,3);
        book(_h["hist_Ch_27"],10,1,3);
        book(_h["hist_E_19.6"],11,1,3);
        book(_h["hist_Ch_19.6"],12,1,3);
        book(_h["hist_E_14.5"],13,1,3);
        book(_h["hist_Ch_14.5"],14,1,3);
        book(_h["hist_E_7.7"],15,1,3);
        book(_h["hist_Ch_7.7"],16,1,3);
      
     //Cu+Cu collisions   
      
        book(_h["hist_E_200_Cu"],17,1,3);
        book(_h["hist_Ch_200_Cu"],18,1,3);
        book(_h["hist_E_62.4_Cu"],19,1,3);
        book(_h["hist_Ch_62.4_Cu"],20,1,3);

    //Cu+Au collisions

        book(_h["hist_E_200_Cu_Au"],21,1,3);
        book(_h["hist_Ch_200_Cu_Au"],22,1,3);

    //U+U collisions

        book(_h["hist_E_193_UU"],23,1,3);
        book(_h["hist_Ch_193_UU"],24,1,3);

    //d+Au collisions

        book(_h["hist_E_200_dAu"],25,1,3);
        book(_h["hist_Ch_200_dAu"],26,1,3);

    //He+Au collisions
    
        book(_h["hist_E_200_He_A"],27,1,3);
        book(_h["hist_Ch_200_He_A"],28,1,3);
      
        
     declare(ALICE::PrimaryParticles(Cuts::abseta < 0.35 && Cuts::pT > 0.0*MeV && Cuts::abscharge > 0), "APRIM");
    }


    void analyze(const Event& event)
    {
        const ParticlePair& beam = beams();
        double NN = 0.;
        string thisbeam = "Empty";
        
        //track beams
        if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
        {
          thisbeam = "AuAu";
          NN = 197.;
          if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3)) thisbeam += "200GeV";
          if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 1E-3)) thisbeam += "62GeV";
          if (fuzzyEquals(sqrtS()/GeV, 39*NN, 1E-3)) thisbeam += "39GeV";
        }
        
        int nch200counter = 0; int nch62counter = 0;
        int nch39counter = 0; double sumET200 = 0.;
        double sumET62 = 0.; double sumET39 = 0.;
        
        //get charged particles
        Particles chargedParticles = applyProjection<ALICE::PrimaryParticles>(event,"APRIM").particles();
        
        //get centrality info
        const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
        const double c = cent();
        if (c > 60) vetoEvent;
        
        //loop over particles to get total charged particles & Et
        for(const Particle& p : chargedParticles)
        {
            const int id = abs(p.pid());
            if(id==11 || id==13) continue; //exclude electron & muon
            
            if(thisbeam == "AuAu200GeV")
            {
                nch200counter ++;
                sumET200 += p.Et()/GeV;
            }
            else if(thisbeam == "AuAu62GeV")
            {
                nch62counter ++;
                sumET62 += p.Et()/GeV;
            }
            else if(thisbeam == "AuAu39GeV")
            {
                nch39counter ++;
                sumET39 += p.Et()/GeV;
            }
        }
        
        
        //fill hists
        if(thisbeam == "AuAu200GeV")
        {
            _hist_Ch_200->fill(c,nch200counter/0.7);
            _hist_E_200->fill(c,sumET200/0.7);
        }
        else if(thisbeam == "AuAu62GeV")
        {
            _hist_Ch_62->fill(c,nch62counter/0.7);
            _hist_E_62->fill(c,sumET62/0.7);
        }
        else if(thisbeam == "AuAu39GeV")
        {
            _hist_Ch_39->fill(c,nch39counter/0.7);
            _hist_E_39->fill(c,sumET39/0.7);
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
    Profile1DPtr _hist_Ch_200;
    Profile1DPtr _hist_Ch_62;
    Profile1DPtr _hist_Ch_39;
    Profile1DPtr _hist_E_200;
    Profile1DPtr _hist_E_62;
    Profile1DPtr _hist_E_39;

    //@}
    string beamOpt;
    enum CollisionSystem {AuAu7,AuAu14,AuAu19,AuAu27,AuAu39,AuAu62,AuAu130,AuAu200,CuCu62,CuCu200,CuAu200,UU193,DAu200,HeAu200,pp200};
    CollisionSystem collSys;


  };


  DECLARE_RIVET_PLUGIN(PHENIX_2016_I1394433);

}
