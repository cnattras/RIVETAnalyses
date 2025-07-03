// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "../Centralities/RHICCentrality.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2013_I1127262 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2013_I1127262);

    //adding a function to get the acc of the event plant dectector
    //so that reaction plant dependency can be calculated
    double GetEventPlaneDetectorAcc(int n, const FinalState pf, vector<Cut> etaRxP, int nPhiSections)
    {
            double QIn = 0.;
            double QRn = 0.;

            for(auto eta : etaRxP)
            {
                    for(int iphi = 0; iphi < nPhiSections; iphi++)
                    {
                            Cut aCut = eta && Cuts::phi > iphi*(2.*M_PI/nPhiSections) && Cuts::phi < (iphi+1)*(2.*M_PI/nPhiSections);
                            Particles particles = pf.particles(aCut);
                            int weight = particles.size();
                            for(const Particle& p : particles)
                            {
                                    QIn += weight*sin(n*p.phi());
                                    QRn += weight*cos(n*p.phi());
                            }
                    }
            }

            double eventPlane = mapAngle0To2Pi((1./n)*atan2(QIn,QRn));

            return eventPlane;
    }

    //adding this funcion to help with RxP calculations
    string Form(double number, int precision)
    {
            std::stringstream stream;
            stream << std::fixed << std::setprecision(precision) << number;

            return stream.str();
    }

    double Resolution(double chi)
    {
        double A = sqrt(M_PI/2.)/2.;

        double funcRes = (A*chi*exp(-chi*chi/4.))*( BesselI0(chi *chi /4.) + BesselI1(chi*chi /4.));

        return funcRes;
    }

    void FillVn(Profile1DPtr vnHisto, const Particles& particles, double eventPlane, int n)
    {
        for(const Particle &p : particles)
        {
            vnHisto->fill(p.pT()/GeV, cos(n*(p.phi() - eventPlane)));
        }
    }

    void FillRAA(Histo1DPtr RAAHisto, const Particles& particles, double eventPlane, int n)
    {
    	for(const Particle &p : particles){
		RAAHisto->fill(p.pT()/GeV, cos(n*(p.phi() - eventPlane)));
	}
    }

    /*void FillRAA(Estimate1DPtr RAAHisto, const Particles& particles, double eventPlane, int n)
    {
        for(const Particle &p : particles){
                RAAHisto->fill(p.pT()/GeV, cos(n*(p.phi() - eventPlane)));
        }
    }*/

    double CalculateChi(double res)
    {
	double chi = 2.;
        double delta = 1.;
        double con = sqrt(M_PI/2.)/2.;
        for ( int i = 0; i < 15; i++)
        {
            chi = (( con*chi*exp(-chi*chi/4.))*( BesselIn(0, chi *chi /4.) + BesselIn(1, chi*chi /4.) )  < res) ? chi + delta : chi - delta ;
            delta = delta / 2.;
        }

        return chi;
    }

    double BesselIn(int n, double x)
    {
	int acc = 40; //accuracy
	double bigN = 1.e10;
            double smallN = 1.e-10;

            if(n < 0){
                    throw UserError("bessel n has to be larger than zero: return 0");
                    return 0;
            }

            if(n == 0) return BesselI0(x);
            if(n == 1) return BesselI1(x);

            if(x == 0) return 0;
            if(abs(x) > bigN) return 0;

            double tox = 2/abs(x);
            double bip = 0.;
            double bim = 0.;
            double bi  = 1.;
            double besselIn = 0.;
            int init = 2*((n + int(sqrt(float(acc*n)))));

            for(int j = init; j>=1; j--)
            {
                    bim = bip + (j*tox*bi);
                    bip = bi;
                    bi  = bim;
		    if(abs(bi) > bigN)
                    {
                            besselIn *= smallN;
                            bi *= smallN;
                            bip *= smallN;
                    }

                    if(j==n) besselIn=bip;
            }

            besselIn *= BesselI0(x)/bi;
	    if((x < 0) && (n%2 == 1)) besselIn = -besselIn;

            return besselIn;
    }

    double BesselI0(double x)
    {
	    double p[7] = {1.0, 3.5156229, 3.0899424, 1.2067492, 0.2659732, 0.360768e-1, 0.45813e-2};
            double q[9] = {0.39894228, 0.1328592e-1, 0.225319e-2, -0.157565e-2, 0.916281e-2, -0.2057706e-1, 0.2635537e-1, -0.1647633e-1, 0.392377e-2};

            double absx = abs(x);
            double y = 0.;
            double besselI0 = 0.;

	    if (absx < 3.75){
                    y = x/3.75;
                    y *= y;
                    besselI0 = p[0]+y*(p[1]+y*(p[2]+y*(p[3]+y*(p[4]+y*(p[5]+y*p[6])))));
            }
            else{
                    y = 3.75/absx;
                    besselI0 = (exp(absx)/sqrt(absx))*(q[0]+y*(q[1]+y*(q[2]+y*(q[3]+y*(q[4]+y*(q[5]+y*(q[6]+y*(q[7]+y*q[8]))))))));
            }

            return besselI0;
    }

    double BesselI1(double x)
    {
	    double p[7] = {0.5, 0.87890594, 0.51498869, 0.15084934, 0.2658733e-1, 0.301532e-2, 0.32411e-3};
            double q[9] = {0.39894228, -0.3988024e-1, -0.362018e-2, 0.163801e-2, -0.1031555e-1, 0.2282967e-1, -0.2895312e-1, 0.1787654e-1, -0.420059e-2};
            double absx = abs(x);
            double y = 0.;
            double besselI1 = 0.;

	    if (absx < 3.75){
                    y = x/3.75;
                    y *= y;
                    besselI1 = x*(p[0]+y*(p[1]+y*(p[2]+y*(p[3]+y*(p[4]+y*(p[5]+y*p[6]))))));
            }
            else{
                    y = 3.75/absx;
                    besselI1 = (exp(absx)/sqrt(absx))*(q[0]+y*(q[1]+y*(q[2]+y*(q[3]+y*(q[4]+y*(q[5]+y*(q[6]+y*(q[7]+y*q[8]))))))));
            }

            if (x < 0) besselI1 = -besselI1;

            return besselI1;
    }

    int GetNDeltaPhi(double evPlane, Particle p)
    {
            int n = 0;

            evPlane = mapAngle0To2Pi(evPlane);

            double deltaPhi = mapAngle0ToPi(evPlane - p.phi());

            if(deltaPhi < M_PI/2.)
            {
                    n = floor(deltaPhi/(M_PI/12.));
            }
            else
            {
                    deltaPhi = mapAngle0ToPi(evPlane + M_PI - p.phi());
                    n = floor(deltaPhi/(M_PI/12.));
            }

            return n;
    }





    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      const UnstableParticles ufs(Cuts::abseta < 0.35 && Cuts::pT > 1*GeV && Cuts::pid == 111);
      declare(ufs, "ufs");

      //booking and delcaring RxP profiles for calculations
      const FinalState RxP(Cuts::abseta > 1. && Cuts::abseta < 2.8);
      declare(RxP, "RxP");

      const FinalState RxPPos(Cuts::eta > 1. && Cuts::eta < 2.8);
      declare(RxPPos, "RxPPos");

      const FinalState RxPNeg(Cuts::eta < -1. && Cuts::eta > -2.8);
      declare(RxPNeg, "RxPNeg");

      book(_p["RxPcosPos"], "RxPcosPos", 6, 0., 6.);

      for(int icent = 0; icent < 6; icent++)
      {
              string refnameCent = mkAxisCode(4, 1, icent+1);
              const Estimate1D& refdataCent = refData(refnameCent);
              book(_h["RAA_pt_AuAu" + to_string(icent)], refnameCent + "_AuAu" + to_string(icent), refdataCent);
              book(_c["sow_AuAu" + to_string(icent)], "sow_AuAu" + to_string(icent));



              for(int idelta = 0; idelta < 6; idelta++)
              {
                      string refname = mkAxisCode(idelta+4, 1, icent+1);
                      //const Estimate1D& refdata = refData(refname);
                      book(_s["RAA_pt_" + to_string(idelta+4) + "1" + to_string(icent+1)], refname);
              }

      }

      string refnamePP = mkAxisCode(4, 1, 1);
      const Estimate1D& refdataPP = refData(refnamePP);
      book(_h["RAA_pt_pp"], refnamePP + "_pp", refdataPP);
      book(_c["sow_pp"], "sow_pp");

      for(unsigned int icent = 0; icent < v2centBins.size()-1; icent++)
      {
            string v2string = "v2_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);
            book(_p[v2string], v2string, refdataPP);
      }

      for(unsigned int ipt = 0; ipt < refdataPP.numBins(); ipt++)
      {
              for(int idelta = 0; idelta < 6; idelta++)
              {
                      for(int icent = 0; icent < 6; icent++)
                      {
                              book(_c["DeltaPhi" + to_string(idelta) + "_pt" + to_string(ipt) + "_cent" + to_string(icent)], "DeltaPhi" + to_string(idelta) + "_pt" + to_string(ipt) + "_cent" + to_string(icent));
                      }

              }
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      double nNucleons = 0.;
      string CollSystem = "Empty";
      const ParticlePair& beam = beams();

      if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
              CollSystem = "AuAu";
              nNucleons = 197.;
              if (fuzzyEquals(sqrtS()/GeV, 200*nNucleons, 5)) CollSystem += "200GeV";
      }
      if (beam.first.pid() == 2212 && beam.second.pid() == 2212)
      {
              CollSystem = "pp";
              nNucleons = 1.;
              if (fuzzyEquals(sqrtS()/GeV, 200*nNucleons, 5)) CollSystem += "200GeV";
      }

      const UnstableParticles& ufs = apply<UnstableParticles>(event, "ufs");
      const Particles particles = ufs.particles();

      //Calcualte Reaction Plane Dependency:

      //filling the RAA histograms for AuAu and pp RAA dependecy calculations
      if(CollSystem == "AuAu200GeV")
      {
              const CentralityProjection& centProj = apply<CentralityProjection>(event,"CMULT");
              double c = centProj();

              if(c >= 60.) vetoEvent;

              //get reaction plane positive and negative final state values
              const FinalState& RxP = apply<FinalState>(event, "RxP");
              const FinalState& RxPPos = apply<FinalState>(event, "RxPPos");
              const FinalState& RxPNeg = apply<FinalState>(event, "RxPNeg");

              //Inner and Outer rings of North and South sections of the RxP dectector
              vector<Cut> etaRxP = {Cuts::eta > 1. && Cuts::eta < 1.5, Cuts::eta > 1.5 && Cuts::eta < 2.8, Cuts::eta < -1. && Cuts::eta > -1.5, Cuts::eta < -1.5 && Cuts::eta > -2.8};
              int nPhiSections = 12;

              double evPPosNeg = GetEventPlaneDetectorAcc(2, RxP, etaRxP, nPhiSections);

              vector<Cut> etaRxPPos = {Cuts::eta > 1. && Cuts::eta < 1.5, Cuts::eta > 1.5 && Cuts::eta < 2.8};
              vector<Cut> etaRxPNeg = {Cuts::eta < -1. && Cuts::eta > -1.5, Cuts::eta < -1.5 && Cuts::eta > -2.8};

              double evPPos = GetEventPlaneDetectorAcc(2, RxPPos, etaRxPPos, nPhiSections);
              double evPNeg = GetEventPlaneDetectorAcc(2, RxPNeg, etaRxPNeg, nPhiSections);

              _p["RxPcosPos"]->fill(int(floor(c/10))+0.5, cos(2*(evPPos-evPNeg)));

              string v2string = "v2_cent" + Form(floor(c/10)*10., 0) + Form((floor(c/10)*10.)+10., 0);
              FillVn(_p[v2string], particles, evPPosNeg, 2);

              //cout << "sow_AuAu" + to_string(int(floor(c/10.))) << endl;

              _c["sow_AuAu" + to_string(int(floor(c/10.)))]->fill();

              for(const auto& p : particles)
              {
                      int n = GetNDeltaPhi(evPPosNeg, p);
                      int ptbin = binRef.indexAt(p.pt()/GeV);
                      if(ptbin >= 1)
                      {
                              _c["DeltaPhi" + to_string(n) + "_pt" + to_string(ptbin-1)  + "_cent" + to_string(int(floor(c/10.)))]->fill();
                      }
                      _h["RAA_pt_AuAu" + to_string(int(floor(c/10.)))]->fill(p.pT()/GeV);
              }
      }
      else if(CollSystem == "pp200GeV")
      {
              _c["sow_pp"]->fill();
              for(const auto& p : particles)
              {
                      _h["RAA_pt_pp"]->fill(p.pT()/GeV);
              }
      }
   }


    /// Normalise histograms etc., after the run
    void finalize() {


      //Reaction Plane Dependency Claculations:

      int centBin = 0;

      std::vector<double> EPres(7, 0.);

     for(auto bin : _p["RxPcosPos"]->bins())
      {
	if(bin.numEntries() > 0)
        {
                if(bin.xMean() <=0) continue;
		double RxPPosRes = sqrt(bin.xMean());
		double chiRxPPos = CalculateChi(RxPPosRes);
		double res = Resolution(sqrt(2)*chiRxPPos);
		EPres[centBin] = res;
                // _s["ResCent"]->addPoint((centBin*10.)+5., res, 5., 0.);
        }
        centBin++;
      }

      for(unsigned int icent = 0; icent < v2centBins.size()-1; icent++)
      {
               string v2string = "v2_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);
	       if(EPres[icent] > 0) _p[v2string]->scale(1, 1./EPres[icent]);
      }

      _h["RAA_pt_pp"]->scaleW(1./_c["sow_pp"]->sumW());

      const vector<double> Ncoll = {955.4, 602.6, 373.8, 219.8, 120.3, 61.};

      for(int icent = 0; icent < 6; icent++)
      {
              _h["RAA_pt_AuAu" + to_string(icent)]->scaleW(1./(Ncoll[icent]*_c["sow_AuAu" + to_string(icent)]->sumW()));
      }





      for(int idelta = 0; idelta < 6; idelta++)
      {
              for(int icent = 0; icent < 6; icent++)
              {
                      divide(_h["RAA_pt_AuAu" + to_string(icent)], _h["RAA_pt_pp"], _s["RAA_pt_" + to_string(idelta+4) + "1" + to_string(icent+1)]);
                      for(unsigned int ipoint = 0; ipoint < _s["RAA_pt_" + to_string(idelta+4) + "1" + to_string(icent+1)]->numBins(); ipoint++)
                      {
                              double sumDeltaPhi = 0.;
                              for(int jdelta = 0; jdelta < 6; jdelta++)
                              {
                                      sumDeltaPhi += (1/6.)*_c["DeltaPhi" + to_string(jdelta) + "_pt" + to_string(ipoint) + "_cent" + to_string(icent)]->effNumEntries();
                              }
                              if(sumDeltaPhi > 0.) _s["RAA_pt_" + to_string(idelta+4) + "1" + to_string(icent+1)]->bin(ipoint+1).scale(_c["DeltaPhi" + to_string(idelta) + "_pt" + to_string(ipoint) + "_cent" + to_string(icent)]->effNumEntries()/sumDeltaPhi);
                      }
              }

      }

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, Estimate1DPtr> _s;
    map<string, CounterPtr> _c;
    std::vector<double> v2centBins = {0., 10., 20., 30., 40., 50., 60.};
    YODA::Histo1D binRef;
    //@}


  };


  RIVET_DECLARE_PLUGIN(PHENIX_2013_I1127262);

}
