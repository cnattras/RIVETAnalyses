// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "Rivet/Math/Math.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2025_I2820229 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2025_I2820229);

// Helper to compute deltaPhi in [-π, π]
static double deltaPhi(double phi1, double phi2) {
    double dphi = phi1 - phi2;
    while (dphi > M_PI) dphi -= 2*M_PI;
    while (dphi <= -M_PI) dphi += 2*M_PI;
    return dphi;
}


    /// Book histograms and initialise projections before the run
    void init() {

// Initialise and register projections
      //John, this is sometihng for you to investigate.  We're supposed ot look at jets with |eta|<0.15.  I'm not sure what the right selection is here.  I've opted to be a bit generous.
      //I also want to flag the "Final State" distinction.  This takes all final state particles.  No feeddown correction, for instance.  It shouldn't matter too much for jets, BUT it maybe slightly does since these jets are small.
      const FinalState fs(Cuts::abseta < 0.5);
        //In case "NONE" is given as option
      const ParticlePair& beam = beams();

      beamOpt = getOption<string>("beam","NONE");
        

      if (beamOpt == "NONE") {
      if (beam.first.pid() == 2212 && beam.second.pid() == 2212){
        float NN = 1.;
        if (fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) collSys = pp;
      }
      }

      if (beamOpt =="PP200") collSys = pp;
      

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.3
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, JetAlg::ANTIKT, 0.3, JetMuons::NONE, JetInvisibles::NONE);
      declare(jetfs, "jets");


      // Book histograms
      // specify custom binning
     
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      // book(h["jetpt"], 1, 1, 1);//this is not real data!  It's the mean pT
      book(h["jetcross"], 1, 1, 1);

      book(h["zg910"], 2, 1, 1);//bin 0
      book(h["zg1012"], 2, 1, 2);//bin 1
      book(h["zg1214"], 2, 1, 3);//bin 2
      book(h["zg1417"], 2, 1, 4);//bin3
      book(h["zg1720"], 2, 1, 5);//bin 4
      book(h["zg2024"], 2, 1, 6);//bin 5
      book(h["zg2429"], 2, 1, 7);//bin 6
           
      book(h["Xi910"], 3, 1, 1);//bin 0

      book(h["Xi1012"], 4, 1, 1);//bin 1

      book(h["Xi1214"], 5, 1, 1);//bin 2
      book(h["Xi1417"], 5, 1, 2);//bin 3
      book(h["Xi1720"], 5, 1, 3);//bin 4
      book(h["Xi2024"], 5, 1, 4);//bin 5
      book(h["Xi2429"], 5, 1, 5);//bin 6

      book(h["jT910"], 6, 1, 1);//bin 0
      book(h["jT1012"], 6, 1, 2);//bin 1
      book(h["jT1214"], 6, 1, 3);//bin 2
      book(h["jT1417"], 6, 1, 4);//bin 3
      book(h["jT1720"], 6, 1, 5);//bin 4
      book(h["jT2024"], 7, 1, 1);//bin 5
      book(h["jT2429"], 7, 1, 2);//bin 6

      book(h["R910"], 8, 1, 1);//bin 0
      book(h["R1012"], 8, 1, 2);//bin 1
      book(h["R1214"], 8, 1, 3);//bin 2
      book(h["R1417"], 8, 1, 4);//bin 3
      book(h["R1720"], 8, 1, 5);//bin 4
      book(h["R2024"], 8, 1, 6);//bin 5
      book(h["R2429"], 8, 1, 7);//bin 6
      // Book counters
      for (int i = 0; i < NPTBINS; ++i) {
        book(c["ptbin" + std::to_string(i)], "ptbin" + std::to_string(i));
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 5*GeV);
      if (collSys == pp){
      // Loop over all jets
      for (const Jet& jet : jets) {
          if(jet.abseta()<0.15){
            h["jetcross"]->fill(jet.pT());

            if(jet.pT()>PTBINS[0]){
              
              //fill counters
              if(jet.pT()>PTBINS[0]&& jet.pT()<PTBINS[1]){//bin 0
                c["ptbin0"]->fill();
              }
              if(jet.pT()>PTBINS[1]&& jet.pT()<PTBINS[2]){//bin 1
                c["ptbin1"]->fill();
              }
              if(jet.pT()>PTBINS[2]&& jet.pT()<PTBINS[3]){//bin 2
                c["ptbin2"]->fill();
              }
              if(jet.pT()>PTBINS[3]&& jet.pT()<PTBINS[4]){//bin 3
                c["ptbin3"]->fill();
              }
              if(jet.pT()>PTBINS[4]&& jet.pT()<PTBINS[5]){//bin 4
                c["ptbin4"]->fill();
              }
              if(jet.pT()>PTBINS[5]&& jet.pT()<PTBINS[6]){//bin 5
                c["ptbin5"]->fill();
              }
              if(jet.pT()>PTBINS[6]&& jet.pT()<PTBINS[7]){//bin 6
                c["ptbin6"]->fill();
              }


            for (const Particle& p : jet.particles()) {
              double xi = - log( p.p() / jet.p());
              //double dr = deltaR(p, jet);
              double deta = p.eta() - jet.eta();
              double dphi = deltaPhi(jet.phi(), p.phi());
              double dr =  sqrt(deta*deta + dphi*dphi);
              //There are edge cases where jT is infinite so we protect against that when filling the histos
              double jTOverpTjet = jet.p3().cross(p.p3()).mod()/jet.p()/jet.pT();
              //cout<<"p "<<p.p()<<" pjet "<<jet.p()<<" jTOverpTjet "<< jet.p3().cross(p.p3()).mod()/jet.p()<<" jTOverpTjet/pTjet "<<jT<<endl;
              // cout<<"xi "<<xi<<" dr "<<dr<<endl;
              if(jet.pT()>PTBINS[0]&& jet.pT()<PTBINS[1]){//bin 0
                h["Xi910"]->fill(xi);
                h["R910"]->fill(dr);
                h["jT910"]->fill(jTOverpTjet,1.0/jTOverpTjet);
              }
              if(jet.pT()>PTBINS[1]&& jet.pT()<PTBINS[2]){//bin 1
                h["Xi1012"]->fill(xi);
                h["R1012"]->fill(dr);
                h["jT1012"]->fill(jTOverpTjet,1.0/jTOverpTjet);
              }
              if(jet.pT()>PTBINS[2]&& jet.pT()<PTBINS[3]){//bin 2
                h["Xi1214"]->fill(xi);
                h["R1214"]->fill(dr);
                h["jT1214"]->fill(jTOverpTjet,1.0/jTOverpTjet);
              }
              if(jet.pT()>PTBINS[3]&& jet.pT()<PTBINS[4]){//bin 3
                h["Xi1417"]->fill(xi);
                h["R1417"]->fill(dr);
                h["jT1417"]->fill(jTOverpTjet,1.0/jTOverpTjet);
              }
              if(jet.pT()>PTBINS[4]&& jet.pT()<PTBINS[5]){//bin 4
                h["Xi1720"]->fill(xi);
                h["R1720"]->fill(dr);
                h["jT1720"]->fill(jTOverpTjet,1.0/jTOverpTjet);
              }
              if(jet.pT()>PTBINS[5]&& jet.pT()<PTBINS[6]){//bin 5
                h["Xi2024"]->fill(xi);
                h["R2024"]->fill(dr);
                h["jT2024"]->fill(jTOverpTjet,1.0/jTOverpTjet);
              }
              if(jet.pT()>PTBINS[6]&& jet.pT()<PTBINS[7]){//bin 6
                h["Xi2429"]->fill(xi);
                h["R2429"]->fill(dr);
                h["jT2429"]->fill(jTOverpTjet,1.0/jTOverpTjet);
              }
            }
          

              //John, I will TOTALLY admit to writing this with ChatGPT.  Please double check!
              fastjet::contrib::SoftDrop sd(beta, zcut);

              // Apply SoftDrop grooming
              fastjet::PseudoJet fjJet = jet.pseudojet();
              fastjet::PseudoJet sdJet = sd(fjJet);

              // Check if the jet has substructure (was groomed)
              if (sdJet.pieces().size() < 2) continue;

              // Calculate zg
              double pt1 = sdJet.pieces()[0].pt();
              double pt2 = sdJet.pieces()[1].pt();
              double zg = std::min(pt1, pt2) / (pt1 + pt2);
              if(jet.pT()>PTBINS[0]&& jet.pT()<PTBINS[1]){//bin 0
                h["zg910"]->fill(zg);
              }
              if(jet.pT()>PTBINS[1]&& jet.pT()<PTBINS[2]){//bin 1
                h["zg1012"]->fill(zg);
              }
              if(jet.pT()>PTBINS[2]&& jet.pT()<PTBINS[3]){//bin 2
                h["zg1214"]->fill(zg);
              }
              if(jet.pT()>PTBINS[3]&& jet.pT()<PTBINS[4]){//bin 3
                h["zg1417"]->fill(zg);
              }
              if(jet.pT()>PTBINS[4]&& jet.pT()<PTBINS[5]){//bin 4
                h["zg1720"]->fill(zg);
              }
              if(jet.pT()>PTBINS[5]&& jet.pT()<PTBINS[6]){//bin 5
                h["zg2024"]->fill(zg);
              }
              if(jet.pT()>PTBINS[6]&& jet.pT()<PTBINS[7]){//bin 6
                h["zg2429"]->fill(zg);
              }
            
          }
}
      }

    }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

       // factor 0.3 needed because it is differential in deta
       // const double sf =0.3* crossSection()/millibarn;

      double normalization = crossSection()/millibarn/sumOfWeights()/(0.3); 
       // std::cout<< "sf: " << sf << std::endl;
      scale(h["jetcross"], normalization);


      for (int i = 0; i < NPTBINS; ++i) {
        scale(h["zg" + std::to_string(static_cast<int>(PTBINS[i])) + std::to_string(static_cast<int>(PTBINS[i + 1]))], 1.0 / (c["ptbin" + std::to_string(i)]->sumW()));
        scale(h["Xi" + std::to_string(static_cast<int>(PTBINS[i])) + std::to_string(static_cast<int>(PTBINS[i + 1]))], 1.0 / (c["ptbin" + std::to_string(i)]->sumW()));
        scale(h["R" + std::to_string(static_cast<int>(PTBINS[i])) + std::to_string(static_cast<int>(PTBINS[i + 1]))], 1.0 / (c["ptbin" + std::to_string(i)]->sumW()));
        scale(h["jT" + std::to_string(static_cast<int>(PTBINS[i])) + std::to_string(static_cast<int>(PTBINS[i + 1]))], 1.0 / (c["ptbin" + std::to_string(i)]->sumW()));
      }

 
      //normalize(h["XXXX"]); // normalize to unity
      //normalize(h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in pb (no cuts)
      //scale(h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> h;
    map<string, Profile1DPtr> p;
    map<string, CounterPtr> c;
    string beamOpt;
    enum CollisionSystem {NONE, pp};
    CollisionSystem collSys;
    /// @}
  private:

      int NPTBINS = 7;
      double PTBINS[8] =  {  9.0,    10.0,    12.0,    14.5,    17.5,    20.5,    24.5,    29.0};

        // Define SoftDrop parameters
        const double beta = 0.0;
        const double zcut = 0.1;

  };


  RIVET_DECLARE_PLUGIN(PHENIX_2025_I2820229);

}
