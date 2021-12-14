// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2020_I1755387 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2020_I1755387);


    /// @name Analysis methods
    ///@{
    double GetRho(const Jets& jetsSet, int nLeadJetExclud)
    {
        if(!jetsSet.size()) return 0.;

        Jets jets = sortBy(jetsSet,cmpMomByPt);
        vector<double> jetPtDensityVector;

        for(auto jet : jets)
        {
            if(nLeadJetExclud == 0)
            {
                    double jetPtDensity = (jet.pT()/GeV)/(jet.pseudojet().area());
                    jetPtDensityVector.push_back(jetPtDensity);
            }
            else nLeadJetExclud--;

        }

        std::sort(jetPtDensityVector.begin(), jetPtDensityVector.end());

        double nMediam = jetPtDensityVector.size();
        double rho = 0.;
        int index = ceil(nMediam/2.)-1;

        if(index < 0)
        {
            printf("WARNING! Only leading jets in the event! Cannot calculate rho. \n");
            return 0.;
        }

        if(int(nMediam)%2 == 1)
        {
            rho = jetPtDensityVector[index]; // median
        }
        else
        {
            rho = 0.5*(jetPtDensityVector[index] + jetPtDensityVector[index+1]); // median
        }


        return rho;
    }

    /// Book histograms and initialise projections before the run
    void init() {

      // Centralities
	declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Initialise and register projections
	//const ALICE::PrimaryParticles aprim(Cuts::abseta < 0.9 && Cuts::abscharge > 0);
        const ALICE::PrimaryParticles aprim((Cuts::abscharge > 0 && Cuts::pT > 0.15*GeV && Cuts::abseta < 0.9) || (Cuts::abscharge == 0 && Cuts::pT > 0.3*GeV && Cuts::abseta < 0.7));
	declare(aprim, "aprim");

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 0.9);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      fastjet::AreaType fjAreaType = fastjet::active_area_explicit_ghosts;
      fastjet::GhostedAreaSpec fjGhostAreaSpec = fastjet::GhostedAreaSpec(1., 1, 0.005, 1., 0.1, 1e-100);

      FastJets jet01(fs, FastJets::ANTIKT, 0.1, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet01, "jets01");
      //FastJets jet02(fs, FastJets::ANTIKT, 0.2, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      //declare(jet02, "jets02");
      fjAreaDef02 = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
      FastJets jet02(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.2, fjAreaDef02, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet02, "jets02");
      FastJets jet03(fs, FastJets::ANTIKT, 0.3, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet03, "jets03");
      //FastJets jet04(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      //declare(jet04, "jets04");
      fjAreaDef04 = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
      FastJets jet04(fs, fastjet::JetAlgorithm::antikt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.4, fjAreaDef04, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet04, "jets04");
      FastJets jet05(fs, FastJets::ANTIKT, 0.5, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet05, "jets05");
      FastJets jet06(fs, FastJets::ANTIKT, 0.6, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jet06, "jets06");

      //Background in PbPb
      fjAreaDef02KT = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
      FastJets jetsKTR02FJ(fs, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.2, fjAreaDef02KT, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetsKTR02FJ, "jetsKTR02FJ");

      fjAreaDef04KT = new fastjet::AreaDefinition(fjGhostAreaSpec, fjAreaType);
      FastJets jetsKTR04FJ(fs, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::pt_scheme, 0.4, fjAreaDef04KT, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetsKTR04FJ, "jetsKTR04FJ");

      // Book histograms
      // specify custom binning
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)

      book(_h["ppspectraR0.1"], 1, 1, 1);
      book(_h["ppspectraR0.2"], 2, 1, 1);
      book(_h["ppspectraR0.3"], 3, 1, 1);
      book(_h["ppspectraR0.4"], 4, 1, 1);
      book(_h["ppspectraR0.5"], 5, 1, 1);
      book(_h["ppspectraR0.6"], 6, 1, 1);
//      book(_s["ppratioR0.1divR0.2"], 13, 1, 1);
//      book(_s["ppratioR0.1divR0.3"], 14, 1, 1);
//      book(_s["ppratioR0.1divR0.4"], 15, 1, 1);
//      book(_s["ppratioR0.1divR0.5"], 16, 1, 1);
//      book(_s["ppratioR0.1divR0.6"], 17, 1, 1);
//      book(_s["ppratioR0.2divR0.3"], 18, 1, 1);
//      book(_s["ppratioR0.2divR0.4"], 19, 1, 1);
//      book(_s["ppratioR0.2divR0.5"], 20,1, 1);
//      book(_s["ppratioR0.2divR0.6"], 21, 1, 1);
      book(_h["ppspectra5GeVleadtrackR0.2"], 22, 1, 1);
      book(_h["pbspectra5GeVleadtrackR0.2"], 23, 1, 1);
      book(_h["ppspectra7GeVleadtrackR0.4"], 24, 1, 1);
      book(_h["pbspectra7GeVleadtrackR0.4"], 25, 1, 1);
//      book(_s["ppcrossleadtrackbias5div0"], 26, 1, 1);
//      book(_s["ppcrossleadtrackbias7div0"], 27, 1, 1);
      book(_s["ppcrossleadtrackbias7div5"], 28, 1, 1);
//      book(_s["pbcross7to5leadcharge"], 29, 1, 1);
      // book(_s["jetRaa5GeVleadtrackR0.2"], 30, 1, 1);
      // book(_s["jetRaa7GeVleadtrackR0.4"], 31, 1, 1);
      // book(_s["jetRaa5GeVpbleadtrackR0.2"], 32, 1, 1);
      // book(_s["jetRaa7GeVpbleadtrackR0.4"], 33, 1, 1);

string refname13 = mkAxisCode(13, 1, 1);
const Scatter2D& refdata13 = refData(refname13);
book(_h["ppratioR0.1divR0.2_R0.1"], refname13 + "_R0.1", refdata13);
book(_h["ppratioR0.1divR0.2_R0.2"], refname13 + "_R0.2", refdata13);
book(_s["ppratioR0.1divR0.2"], refname13);

string refname14 = mkAxisCode(14, 1, 1);
const Scatter2D& refdata14 = refData(refname14);
book(_h["ppratioR0.1divR0.3_R0.1"], refname14 + "_R0.1", refdata14);
book(_h["ppratioR0.1divR0.3_R0.3"], refname14 + "_R0.3", refdata14);
book(_s["ppratioR0.1divR0.3"], refname14);

string refname15 = mkAxisCode(15, 1, 1);
const Scatter2D& refdata15 = refData(refname15);
book(_h["ppratioR0.1divR0.4_R0.1"], refname15 + "_R0.1", refdata15);
book(_h["ppratioR0.1divR0.4_R0.4"], refname15 + "_R0.4", refdata15);
book(_s["ppratioR0.1divR0.4"], refname15);

string refname16 = mkAxisCode(16, 1, 1);
const Scatter2D& refdata16 = refData(refname16);
book(_h["ppratioR0.1divR0.5_R0.1"], refname16 + "_R0.1", refdata16);
book(_h["ppratioR0.1divR0.5_R0.5"], refname16 + "_R0.5", refdata16);
book(_s["ppratioR0.1divR0.5"], refname16);

string refname17 = mkAxisCode(17, 1, 1);
const Scatter2D& refdata17 = refData(refname17);
book(_h["ppratioR0.1divR0.6_R0.1"], refname17 + "_R0.1", refdata17);
book(_h["ppratioR0.1divR0.6_R0.6"], refname17 + "_R0.6", refdata17);
book(_s["ppratioR0.1divR0.6"], refname17);

string refname18 = mkAxisCode(18, 1, 1);
const Scatter2D& refdata18 = refData(refname18);
book(_h["ppratioR0.2divR0.3_R0.2"], refname18 + "_R0.2", refdata18);
book(_h["ppratioR0.2divR0.3_R0.3"], refname18 + "_R0.3", refdata18);
book(_s["ppratioR0.2divR0.3"], refname18);

string refname19 = mkAxisCode(19, 1, 1);
const Scatter2D& refdata19 = refData(refname19);
book(_h["ppratioR0.2divR0.4_R0.2"], refname19 + "_R0.2", refdata19);
book(_h["ppratioR0.2divR0.4_R0.4"], refname19 + "_R0.4", refdata19);
book(_s["ppratioR0.2divR0.4"], refname19);

string refname20 = mkAxisCode(20, 1, 1);
const Scatter2D& refdata20 = refData(refname20);
book(_h["ppratioR0.2divR0.5_R0.2"], refname20 + "_R0.2", refdata20);
book(_h["ppratioR0.2divR0.5_R0.5"], refname20 + "_R0.5", refdata20);
book(_s["ppratioR0.2divR0.5"], refname20);

string refname21 = mkAxisCode(21, 1, 1);
const Scatter2D& refdata21 = refData(refname21);
book(_h["ppratioR0.2divR0.6_R0.2"], refname21 + "_R0.2", refdata21);
book(_h["ppratioR0.2divR0.6_R0.6"], refname21 + "_R0.6", refdata21);
book(_s["ppratioR0.2divR0.6"], refname21);

string refname26 = mkAxisCode(26, 1, 1);
const Scatter2D& refdata26 = refData(refname26);
book(_h["ppcrossleadtrackbias5"], refname26 + "_5GeV_R0.2", refdata26);
book(_s["ppcrossleadtrackbias5div0"], refname26);

string refname27 = mkAxisCode(27, 1, 1);
const Scatter2D& refdata27 = refData(refname27);
book(_h["ppcrossleadtrackbias7"], refname27 + "_7GeV_R0.4", refdata27);
book(_s["ppcrossleadtrackbias7div0"], refname27);

string refname29 = mkAxisCode(29, 1, 1);
const Scatter2D& refdata29 = refData(refname29);
book(_h["pbcrossleadtrackbias7R0.2"], refname29 + "_7GeV_R0.2", refdata29);
book(_h["pbcrossleadtrackbias5R0.2"], refname29 + "_5GeV_R0.2", refdata29);
book(_s["pbcrossleadtrackbias7div5"], refname29);


string refname30 = mkAxisCode(30, 1, 1);
const Scatter2D& refdata30 = refData(refname30);
book(_h["pbscaledspectrumtrackbias5R0.2Table30Figure6"], refname30 + "_5GeV_R0.2pb", refdata30);
book(_h["ppcrossleadtrackbias5R0.2Table30Figure6"], refname30 + "_5GeV_R0.2pp", refdata30);
book(_s["RAAbias5R0.2"], refname30);

string refname31 = mkAxisCode(31, 1, 1);
const Scatter2D& refdata31 = refData(refname31);
book(_h["pbscaledspectrumtrackbias7R0.4Table31Figure6"], refname31 + "_7GeV_R0.4pb", refdata31);
book(_h["ppcrossleadtrackbias7R0.4Table31Figure6"], refname31 + "_7GeV_R0.4pp", refdata31);
book(_s["RAAbias7R0.4"], refname31);



string refname32 = mkAxisCode(32, 1, 1);
const Scatter2D& refdata32 = refData(refname32);
book(_h["pbscaledspectrumtrackbias5R0.2Table32Figure6"], refname32 + "_5GeV_R0.2pb", refdata32);
book(_h["ppcrossleadtrackbias0R0.2Table32Figure6"], refname32 + "_0GeV_R0.2pp", refdata32);
book(_s["RAAbias0and5R0.2"], refname32);

string refname33 = mkAxisCode(33, 1, 1);
const Scatter2D& refdata33 = refData(refname33);
book(_h["pbscaledspectrumtrackbias7R0.4Table33Figure6"], refname33 + "_7GeV_R0.4pb", refdata33);
book(_h["ppcrossleadtrackbias0R0.4Table33Figure6"], refname33 + "_0GeV_R0.4pp", refdata33);
book(_s["RAAbias0and7R0.4"], refname33);

      book(_c["sow"], "sow");
      book(_c["sowPBPB"], "sowPBPB");
      book(_c["ppXSec"], "ppXSec");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // beam generation

	const ParticlePair& beam = beams();
            string CollSystem = "Empty";

        if (beam.first.pid() == 1000822080 && beam.second.pid() == 1000822080)
	{
    	CollSystem = "PBPB";
	}

	if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
	{
    	CollSystem = "PP";
	}


      // Retrieve clustered jets, sorted by pT, with a minimum pT cut

        const ALICE::PrimaryParticles aprim = apply<ALICE::PrimaryParticles>(event, "aprim");
	const Particles ALICEparticles = aprim.particles();

        FastJets FJjets02 = apply<FastJets>(event, "jets02");
        FastJets FJjets04 = apply<FastJets>(event, "jets04");

        FJjets02.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
        FJjets04.calc(ALICEparticles); //give ALICE primary particles to FastJet projection

        Jets jets02 = FJjets02.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.5); //get jets (ordered by pT)
        Jets jets04 = FJjets04.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.3); //get jets (ordered by pT)

	if (CollSystem == "PBPB") {
		const CentralityProjection& centProj = apply<CentralityProjection>(event,"V0M");
		const double cent = centProj();
		if (cent >= 10) vetoEvent;
                _c["sowPBPB"]->fill();

                FastJets FJjets02KT = apply<FastJets>(event, "jetsKTR02FJ");
                FastJets FJjets04KT = apply<FastJets>(event, "jetsKTR04FJ");

                FJjets02KT.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
                FJjets04KT.calc(ALICEparticles); //give ALICE primary particles to FastJet projection

                Jets jets02KT = FJjets02KT.jetsByPt(Cuts::abseta < 0.5); //get jets (ordered by pT)
                Jets jets04KT = FJjets04KT.jetsByPt(Cuts::abseta < 0.3); //get jets (ordered by pT)

                double rho02 = GetRho(jets02KT, 2);
                double rho04 = GetRho(jets04KT, 2);

	        for(auto jet : jets02)
       		 {
       		         if(jet.particles(Cuts::pT > 5.*GeV).size() > 0){
                        _h["pbspectra5GeVleadtrackR0.2"]->fill(jet.pT()/GeV - rho02*jet.pseudojet().area());//histo 23
                        _h["pbcrossleadtrackbias5R0.2"]->fill(jet.pT()/GeV - rho02*jet.pseudojet().area());//histo 29
                        _h["pbscaledspectrumtrackbias5R0.2Table30Figure6"]->fill(jet.pT()/GeV - rho02*jet.pseudojet().area());//histo 30
                        _h["pbscaledspectrumtrackbias5R0.2Table32Figure6"]->fill(jet.pT()/GeV - rho02*jet.pseudojet().area());//histo 32
                }
                     if(jet.particles(Cuts::pT > 7.*GeV).size() > 0){
                        _h["pbcrossleadtrackbias7R0.2"]->fill(jet.pT()/GeV - rho02*jet.pseudojet().area());//histo 29
                    }
		}
                for(auto jet : jets04)
                 {
                         if(jet.particles(Cuts::pT > 7.*GeV).size() > 0){
                        _h["pbspectra7GeVleadtrackR0.4"]->fill(jet.pT()/GeV - rho04*jet.pseudojet().area());//histo 25
                        _h["pbscaledspectrumtrackbias7R0.4Table31Figure6"]->fill(jet.pT()/GeV - rho04*jet.pseudojet().area());//histo 31
                        _h["pbscaledspectrumtrackbias7R0.4Table33Figure6"]->fill(jet.pT()/GeV - rho04*jet.pseudojet().area());//histo 33
                }
                }
		return;
		}

        //pp cross-sections
        pair<double,double> cs = HepMCUtils::crossSection(*event.genEvent());
        _c["ppXSec"]->fill(cs.first);

	FastJets FJjets01 = apply<FastJets>(event, "jets01");
	FastJets FJjets03 = apply<FastJets>(event, "jets03");
	FastJets FJjets05 = apply<FastJets>(event, "jets05");
	FastJets FJjets06 = apply<FastJets>(event, "jets06");

	FJjets01.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
        FJjets03.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
        FJjets05.calc(ALICEparticles); //give ALICE primary particles to FastJet projection
        FJjets06.calc(ALICEparticles); //give ALICE primary particles to FastJet projection

	Jets jets01 = FJjets01.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.6); //get jets (ordered by pT)
        Jets jets03 = FJjets03.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.4); //get jets (ordered by pT)
        Jets jets05 = FJjets05.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.2); //get jets (ordered by pT)
        Jets jets06 = FJjets06.jetsByPt(Cuts::pT >= 20.*GeV && Cuts::abseta < 0.1); //get jets (ordered by pT)

	_c["sow"]->fill();

	for(auto jet : jets01)
	{

		_h["ppspectraR0.1"]->fill(jet.pT()/GeV);
		_h["ppratioR0.1divR0.2_R0.1"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.3_R0.1"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.4_R0.1"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.5_R0.1"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.6_R0.1"]->fill(jet.pT()/GeV);

	}

        for(auto jet : jets02)
        {
		if(jet.particles(Cuts::pT > 5.*GeV).size() > 0){
			_h["ppcrossleadtrackbias5"]->fill(jet.pT()/GeV);
 			_h["ppspectra5GeVleadtrackR0.2"]->fill(jet.pT()/GeV);
            _h["ppcrossleadtrackbias5R0.2Table30Figure6"]->fill(jet.pT()/GeV);//Table 30
		}
                if(jet.particles(Cuts::pT > 7.*GeV).size() > 0){
                        _h["ppcrossleadtrackbias7"]->fill(jet.pT()/GeV);
               	}
	        _h["ppspectraR0.2"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.2_R0.2"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.3_R0.2"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.4_R0.2"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.5_R0.2"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.6_R0.2"]->fill(jet.pT()/GeV);
                _h["ppcrossleadtrackbias0R0.2Table32Figure6"]->fill(jet.pT()/GeV);//Table 32

    }

        for(auto jet : jets03)
        {
                _h["ppspectraR0.3"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.3_R0.3"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.3_R0.3"]->fill(jet.pT()/GeV);

	 }

        for(auto jet : jets04)
        {
                if(jet.particles(Cuts::pT > 7.*GeV).size() > 0){
                        _h["ppspectra7GeVleadtrackR0.4"]->fill(jet.pT()/GeV);
                        _h["ppcrossleadtrackbias7R0.4Table31Figure6"]->fill(jet.pT()/GeV);//Table 31
                    }

		 _h["ppspectraR0.4"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.4_R0.4"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.4_R0.4"]->fill(jet.pT()/GeV);
                _h["ppcrossleadtrackbias0R0.4Table33Figure6"]->fill(jet.pT()/GeV);

        }

        for(auto jet : jets05)
        {
                _h["ppspectraR0.5"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.5_R0.5"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.5_R0.5"]->fill(jet.pT()/GeV);

	 }

        for(auto jet : jets06)
        {
                _h["ppspectraR0.6"]->fill(jet.pT()/GeV);
                _h["ppratioR0.1divR0.6_R0.6"]->fill(jet.pT()/GeV);
                _h["ppratioR0.2divR0.6_R0.6"]->fill(jet.pT()/GeV);
        }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

        //if (_c["sow"]->sumW()>0){//pp only histograms

        double ppXSec = _c["ppXSec"]->sumW()/_c["sow"]->sumW();

	_h["ppspectraR0.1"]->scaleW((ppXSec/millibarn)/(1.2*_c["sow"]->sumW()));
        _h["ppspectraR0.2"]->scaleW((ppXSec/millibarn)/(1.*_c["sow"]->sumW()));
        _h["ppspectraR0.3"]->scaleW((ppXSec/millibarn)/(0.8*_c["sow"]->sumW()));
        _h["ppspectraR0.4"]->scaleW((ppXSec/millibarn)/(0.6*_c["sow"]->sumW()));
        _h["ppspectraR0.5"]->scaleW((ppXSec/millibarn)/(0.4*_c["sow"]->sumW()));
        _h["ppspectraR0.6"]->scaleW((ppXSec/millibarn)/(0.2*_c["sow"]->sumW()));

        _h["ppspectra5GeVleadtrackR0.2"]->scaleW((ppXSec/millibarn)/(1.*_c["sow"]->sumW()));
        _h["ppspectra7GeVleadtrackR0.4"]->scaleW((ppXSec/millibarn)/(0.6*_c["sow"]->sumW()));

        _h["ppratioR0.1divR0.2_R0.1"]->scaleW(1./1.2);
        _h["ppratioR0.1divR0.3_R0.1"]->scaleW(1./1.2);
        _h["ppratioR0.1divR0.4_R0.1"]->scaleW(1./1.2);
        _h["ppratioR0.1divR0.5_R0.1"]->scaleW(1./1.2);
        _h["ppratioR0.1divR0.6_R0.1"]->scaleW(1./1.2);
        _h["ppratioR0.1divR0.3_R0.3"]->scaleW(1./0.8);
        _h["ppratioR0.2divR0.3_R0.3"]->scaleW(1./0.8);
        _h["ppratioR0.1divR0.4_R0.4"]->scaleW(1./0.6);
        _h["ppratioR0.2divR0.4_R0.4"]->scaleW(1./0.6);
        _h["ppratioR0.1divR0.5_R0.5"]->scaleW(1./0.4);
        _h["ppratioR0.2divR0.5_R0.5"]->scaleW(1./0.4);
        _h["ppratioR0.1divR0.6_R0.6"]->scaleW(1./0.2);
        _h["ppratioR0.2divR0.6_R0.6"]->scaleW(1./0.2);


	divide(_h["ppratioR0.1divR0.2_R0.1"], _h["ppratioR0.1divR0.2_R0.2"], _s["ppratioR0.1divR0.2"]);
        divide(_h["ppratioR0.1divR0.3_R0.1"], _h["ppratioR0.1divR0.3_R0.3"], _s["ppratioR0.1divR0.3"]);
        divide(_h["ppratioR0.1divR0.4_R0.1"], _h["ppratioR0.1divR0.4_R0.4"], _s["ppratioR0.1divR0.4"]);
        divide(_h["ppratioR0.1divR0.5_R0.1"], _h["ppratioR0.1divR0.5_R0.5"], _s["ppratioR0.1divR0.5"]);
        divide(_h["ppratioR0.1divR0.6_R0.1"], _h["ppratioR0.1divR0.6_R0.6"], _s["ppratioR0.1divR0.6"]);
        divide(_h["ppratioR0.2divR0.3_R0.2"], _h["ppratioR0.2divR0.3_R0.3"], _s["ppratioR0.2divR0.3"]);
        divide(_h["ppratioR0.2divR0.4_R0.2"], _h["ppratioR0.2divR0.4_R0.4"], _s["ppratioR0.2divR0.4"]);
        divide(_h["ppratioR0.2divR0.5_R0.2"], _h["ppratioR0.2divR0.5_R0.5"], _s["ppratioR0.2divR0.5"]);
        divide(_h["ppratioR0.2divR0.6_R0.2"], _h["ppratioR0.2divR0.6_R0.6"], _s["ppratioR0.2divR0.6"]);

        //Figure 5 left
        _h["ppcrossleadtrackbias5"]->scaleW((ppXSec/millibarn)/_c["sow"]->sumW());
        divide(_h["ppcrossleadtrackbias5"], _h["ppspectraR0.2"], _s["ppcrossleadtrackbias5div0"]);
        _h["ppcrossleadtrackbias7"]->scaleW((ppXSec/millibarn)/_c["sow"]->sumW());
        divide(_h["ppcrossleadtrackbias7"], _h["ppspectraR0.2"], _s["ppcrossleadtrackbias7div0"]);
        divide(_h["ppcrossleadtrackbias7"], _h["ppcrossleadtrackbias5"], _s["ppcrossleadtrackbias7div5"]);

        //Figure 6 histos
        _h["ppcrossleadtrackbias5R0.2Table30Figure6"]->scaleW((ppXSec/millibarn)/_c["sow"]->sumW());
        _h["ppcrossleadtrackbias7R0.4Table31Figure6"]->scaleW((ppXSec/millibarn)/(0.6*_c["sow"]->sumW()));
        _h["ppcrossleadtrackbias0R0.2Table32Figure6"]->scaleW((ppXSec/millibarn)/_c["sow"]->sumW());
        _h["ppcrossleadtrackbias0R0.4Table33Figure6"]->scaleW((ppXSec/millibarn)/(0.6*_c["sow"]->sumW()));
        //}

        //if (_c["sowPBPB"]->sumW()>0){//pb only histograms
        divide(_h["pbcrossleadtrackbias7R0.2"], _h["pbcrossleadtrackbias5R0.2"], _s["pbcrossleadtrackbias7div5"]);//Fig. 5 right

        float TAA = 23.07;
    _h["pbspectra5GeVleadtrackR0.2"]->scaleW(1.0/(TAA*_c["sowPBPB"]->sumW()));//Table 23, figure 4
    _h["pbspectra7GeVleadtrackR0.4"]->scaleW(1.0/(0.6*TAA*_c["sowPBPB"]->sumW()));//Table 25, figure 5
    //Figure 6 histos
    _h["pbscaledspectrumtrackbias5R0.2Table30Figure6"]->scaleW(1.0/(TAA*_c["sowPBPB"]->sumW()));//Table 30, figure 6
    _h["pbscaledspectrumtrackbias7R0.4Table31Figure6"]->scaleW(1.0/(0.6*TAA*_c["sowPBPB"]->sumW()));//Table 31, figure 6
    _h["pbscaledspectrumtrackbias5R0.2Table32Figure6"]->scaleW(1.0/(TAA*_c["sowPBPB"]->sumW()));//Table 32, figure 6
    _h["pbscaledspectrumtrackbias7R0.4Table33Figure6"]->scaleW(1.0/(0.6*TAA*_c["sowPBPB"]->sumW()));//Table 33, figure 6



	//}
//pb and pp histograms
        //if (_c["sow"]->sumW()>0 && _c["sowPBPB"]->sumW()>0){//pp only histograms
//Figure 6
        divide(_h["pbscaledspectrumtrackbias5R0.2Table30Figure6"], _h["ppcrossleadtrackbias5R0.2Table30Figure6"], _s["RAAbias5R0.2"]);//Table 30
        divide(_h["pbscaledspectrumtrackbias7R0.4Table31Figure6"], _h["ppcrossleadtrackbias7R0.4Table31Figure6"], _s["RAAbias7R0.4"]);//Table 31
        divide(_h["pbscaledspectrumtrackbias5R0.2Table32Figure6"], _h["ppcrossleadtrackbias0R0.2Table32Figure6"], _s["RAAbias0and5R0.2"]);//Table 32
        divide(_h["pbscaledspectrumtrackbias7R0.4Table33Figure6"], _h["ppcrossleadtrackbias0R0.4Table33Figure6"], _s["RAAbias0and7R0.4"]);//Table 33

        //}

	}
    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;

    fastjet::AreaDefinition *fjAreaDef02;
    fastjet::AreaDefinition *fjAreaDef04;

    fastjet::AreaDefinition *fjAreaDef02KT;
    fastjet::AreaDefinition *fjAreaDef04KT;

    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2020_I1755387);

}
