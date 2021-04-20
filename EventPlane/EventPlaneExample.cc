// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "../Centralities/RHICCentrality.hh"
//#include "Rivet/Projections/EventPlane.hh"
#include <cmath>
#include <iostream>
//#include "/usr/include/boost/math/special_functions/modf.hpp"

namespace Rivet {


  /// @brief Add a short analysis description here
  class EventPlaneExample : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(EventPlaneExample);

    double GetEventPlanePtWeight(int n, Particles particles)
    {
            double QIn = 0.;
            double QRn = 0.;

            for(const Particle& p : particles)
            {
                QIn += (p.pT()/GeV)*sin(n*p.phi());
                QRn += (p.pT()/GeV)*cos(n*p.phi());
            }

            double eventPlane = (1./n)*atan2(QIn,QRn);

            if(eventPlane < 0) eventPlane += 2.*M_PI;

            return eventPlane;
    }

    double GetEventPlane(int n, Particles particles)
    {
            double QIn = 0.;
            double QRn = 0.;

            for(const Particle& p : particles)
            {
                QIn += sin(n*p.phi());
                QRn += cos(n*p.phi());
            }

            double eventPlane = (1./n)*atan2(QIn,QRn);

            if(eventPlane < 0) eventPlane += 2.*M_PI;

            return eventPlane;
    }

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

            double eventPlane = (1./n)*atan2(QIn,QRn);

            if(eventPlane < 0) eventPlane += 2.*M_PI;

            return eventPlane;
    }

    double CalculateVn(const Particles& particles, double eventPlane, int n, double minPt, double maxPt, int nBins=-1)
    {
        if (nBins == -1) nBins = sqrt(particles.size());

        std::vector<double> phiBins(nBins,0.0);

        double binWidth = 2*M_PI/nBins;

        for(const Particle &p : particles)
        {
            if((p.pT()/GeV > minPt) && (p.pT()/GeV < maxPt))
            {
                phiBins[static_cast<int>(mapAngle0To2Pi(p.phi() - eventPlane)/binWidth)] += p.pT()/GeV;
            }
        }

        double integral = 0.;
        double binCenter = 0.5*binWidth;
        double Vn = 0.;

        for (int i = 0; i < nBins; i++)
        {
            integral += phiBins[i];
            Vn += phiBins[i]*cos(n*((i*binWidth) + binCenter));
        }

        Vn /= integral;
        return Vn;
    }

    void FillVn(Histo1DPtr vnHisto, const Particles& particles, double eventPlane, double minPt, double maxPt)
    {
        for(const Particle &p : particles)
        {
            if((p.pT()/GeV > minPt) && (p.pT()/GeV < maxPt))
            {
                vnHisto->fill(mapAngle0To2Pi(p.phi() - eventPlane), p.pT()/GeV);
            }
        }
    }

    void FillVn(Histo1DPtr vnHisto, const Particles& particles, double eventPlane)
    {
        for(const Particle &p : particles)
        {
            vnHisto->fill(mapAngle0To2Pi(p.phi() - eventPlane), p.pT()/GeV);
        }
    }

    double GetVn(Histo1DPtr vnHisto, int n)
    {
        double integral = 0.;
        double Vn = 0.;

        for(auto bin : vnHisto->bins())
        {
            integral += vnHisto->sumW();
            Vn += vnHisto->sumW()*cos(n*bin.xMid());
        }

        if(integral > 0.) Vn /= integral;
        else return 0.;

        return Vn;
    }

    double CalculateVnAlt(const Particles& particles, double eventPlane, int n, double minPt, double maxPt)
    {

        double integral = 0.;
        double Vn = 0.;

        for(const Particle &p : particles)
        {
            if ((p.pT()/GeV > minPt) && (p.pT()/GeV < maxPt))
            {
                integral += p.pT()/GeV;
                Vn += (p.pT()/GeV)*cos(n*mapAngle0To2Pi(p.phi() - eventPlane));
            }
        }

        Vn /= integral;
        return Vn;
    }

    double CalculateChi(double res)
    {
        //Implementation from A. M. Poskanzer and S. A. Voloshin, Phys. Rev. C 58, 1671 – Published 1 September 1998
        double chi = 2.;
        double delta = 1.;
        double con = sqrt(M_PI)/2.;
        for ( int i = 0; i < 15; i++)
        {
            chi = (( con*chi*exp(-chi*chi/2.))*( BesselI(0, chi *chi /2.) + BesselI(1, chi*chi /2.) )  < res) ? chi + delta : chi - delta ;
            delta = delta / 2.;
        }

        return chi;
    }

    double BesselI0(double x)
    {
            //From NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)

            //Polynomial parameters
            double p[7] = {1.0, 3.5156229, 3.0899424, 1.2067492, 0.2659732, 0.360768e-1, 0.45813e-2};
            double q[9] = {0.39894228, 0.1328592e-1, 0.225319e-2, -0.157565e-2, 0.916281e-2, -0.2057706e-1, 0.2635537e-1, -0.1647633e-1, 0.392377e-2};

            double absx = abs(x);
            double y = 0.;
            double besselI0 = 0.;

            //Polynomial fit
            if (absx < 3.75)
            {
                    y = x/3.75;
                    y *= y;
                    besselI0 = p[0]+y*(p[1]+y*(p[2]+y*(p[3]+y*(p[4]+y*(p[5]+y*p[6])))));
            }
            else
            {
                    y = 3.75/absx;
                    besselI0 = (exp(absx)/sqrt(absx))*(q[0]+y*(q[1]+y*(q[2]+y*(q[3]+y*(q[4]+y*(q[5]+y*(q[6]+y*(q[7]+y*q[8]))))))));
            }

            return besselI0;
    }

    double BesselI1(double x)
    {
            //From NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)

            //Polynomial parameters
            double p[7] = {0.5, 0.87890594, 0.51498869, 0.15084934, 0.2658733e-1, 0.301532e-2, 0.32411e-3};

            double q[9] = {0.39894228, -0.3988024e-1, -0.362018e-2, 0.163801e-2, -0.1031555e-1, 0.2282967e-1, -0.2895312e-1, 0.1787654e-1, -0.420059e-2};

            double absx = abs(x);

            double y = 0.;
            double besselI1 = 0.;

            //Polynomial fit
            if (absx < 3.75)
            {
                    y = x/3.75;
                    y *= y;
                    besselI1 = x*(p[0]+y*(p[1]+y*(p[2]+y*(p[3]+y*(p[4]+y*(p[5]+y*p[6]))))));
            }
            else
            {
                    y = 3.75/absx;
                    besselI1 = (exp(absx)/sqrt(absx))*(q[0]+y*(q[1]+y*(q[2]+y*(q[3]+y*(q[4]+y*(q[5]+y*(q[6]+y*(q[7]+y*q[8]))))))));
            }

            if (x < 0) besselI1 = -besselI1;

            return besselI1;
    }

    double BesselI(int n,double x)
    {
            //From NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)

            int acc = 40; //accuracy
            double bigN = 1.e10;
            double smallN = 1.e-10;

            if(n < 0)
            {
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

                    // Renormalise to prevent overflows
                    if(abs(bi) > bigN)
                    {
                            besselIn *= smallN;
                            bi *= smallN;
                            bip *= smallN;
                    }

                    if(j==n) besselIn=bip;
            }

            besselIn *= BesselI0(x)/bi; // Normalise with BesselI0(x)

            if((x < 0) && (n%2 == 1)) besselIn = -besselIn;

            return besselIn;
    }

    double Resolution(double chi)
    {
        double A = sqrt(M_PI/2.)/2.;

        double funcRes = (A*chi*exp(-chi*chi/4.))*( BesselI0(chi *chi /4.) + BesselI1(chi*chi /4.));

        return funcRes;
    }

    string Form(double number, int precision)
    {
            std::stringstream stream;
            stream << std::fixed << std::setprecision(precision) << number;

            return stream.str();
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      const FinalState fs(Cuts::abseta < 0.5 && Cuts::pT > 0.150*GeV);
      declare(fs, "fs");

      const FinalState RxP(Cuts::abseta > 1. && Cuts::abseta < 2.8);
      declare(RxP, "RxP");

      const FinalState RxPPos(Cuts::eta > 1. && Cuts::eta < 2.8);
      declare(RxPPos, "RxPPos");

      const FinalState RxPNeg(Cuts::eta < -1. && Cuts::eta > -2.8);
      declare(RxPNeg, "RxPNeg");

      book(_p["RxPcosPos"], "RxPcosPos", 10, 0., 10.);
      book(_s["ResCent"], "ResCent");






      for(int icent = 0; icent < v2centBins.size()-1; icent++)
      {
              string v2string = "v2_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);

              book(_s[v2string], v2string);

              for(int ipt = 0; ipt < v2ptBins.size()-1; ipt++)
              {
                      string phiEP = "phiEP_pt_" + Form(v2ptBins[ipt], 2) + "_" + Form(v2ptBins[ipt+1], 2) + v2string;
                      book(_h[phiEP], phiEP,  36, 0, 2*M_PI);
              }
      }



    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      //add calculation of reaction plane angle
      const double c = cent();

      if(c >= 50.) vetoEvent;

      const FinalState& fs = apply<FinalState>(event, "fs");

      const FinalState& RxP = apply<FinalState>(event, "RxP");
      const FinalState& RxPPos = apply<FinalState>(event, "RxPPos");
      const FinalState& RxPNeg = apply<FinalState>(event, "RxPNeg");

      //Inner and Outer rings of the North and South sections of the RxP detector
      vector<Cut> etaRxP = {Cuts::eta > 1. && Cuts::eta < 1.5, Cuts::eta > 1.5 && Cuts::eta < 2.8, Cuts::eta < -1. && Cuts::eta > -1.5, Cuts::eta < -1.5 && Cuts::eta > -2.8};
      int nPhiSections = 12;

      double evPPosNeg = GetEventPlaneDetectorAcc(2, RxP, etaRxP, nPhiSections);

      vector<Cut> etaRxPPos = {Cuts::eta > 1. && Cuts::eta < 1.5, Cuts::eta > 1.5 && Cuts::eta < 2.8};
      vector<Cut> etaRxPNeg = {Cuts::eta < -1. && Cuts::eta > -1.5, Cuts::eta < -1.5 && Cuts::eta > -2.8};

      double evPPos = GetEventPlaneDetectorAcc(2, RxPPos, etaRxPPos, nPhiSections);
      double evPNeg = GetEventPlaneDetectorAcc(2, RxPNeg, etaRxPNeg, nPhiSections);

      double deltaEP = abs(evPPos-evPNeg);

      if(deltaEP > M_PI) deltaEP = (2.*M_PI) - deltaEP;

      //std::floor(double a) returns the largest integer value smaller than a
      _p["RxPcosPos"]->fill(int(floor(c/10))+0.5, cos(2*(deltaEP)));

      //Particles particles = fs.particles();

      for(int ipt = 0; ipt < v2ptBins.size()-1; ipt++)
      {
              string v2string = "v2_cent" + Form(floor(c/10)*10., 0) + Form((floor(c/10)*10.)+10., 0);
              string phiEP = "phiEP_pt_" + Form(v2ptBins[ipt], 2) + "_" + Form(v2ptBins[ipt+1], 2) + v2string;
              //FillVn(_h[phiEP], particles, evPPosNeg, v2ptBins[i], v2ptBins[i+1]);
              Particles particles = fs.particles(Cuts::pT > v2ptBins[ipt]*GeV && Cuts::pT < v2ptBins[ipt+1]*GeV);

              FillVn(_h[phiEP], particles, evPPosNeg);
      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {

            int centBin = 0;

            std::vector<double> EPres(5, 0.);

            for(auto bin : _p["RxPcosPos"]->bins())
            {
                    //if(bin.numEntries() > 0 && bin.mean() > 0.)
                    if(bin.numEntries() > 0)
                    {
                            double RxPPosRes = sqrt(bin.mean()); //resolution of one of the RxP sides (positive eta)
                            double chiRxPPos = CalculateChi(RxPPosRes); //Chi of one of the RxP sides (positive eta)
                            //Using the approx chiRxP(pos+neg) = sqrt(2)*chiRxPPos
                            double res = Resolution(sqrt(2)*chiRxPPos); //resolution of RxP(pos+neg)
                            EPres[centBin] = res;
                            _s["ResCent"]->addPoint((centBin*10.)+5., res, 5., 0.);
                    }
                    centBin++;
            }

            for(int icent = 0; icent < v2centBins.size()-1; icent++)
            {
                    string v2string = "v2_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);

                    for(int ipt = 0; ipt < v2ptBins.size()-1; ipt++)
                    {
                            string phiEP = "phiEP_pt_" + Form(v2ptBins[ipt], 2) + "_" + Form(v2ptBins[ipt+1], 2) + v2string;
                            double v2raw = GetVn(_h[phiEP], 2);
                            double v2 = 0.;
                            if(EPres[icent] > 0) v2 = v2raw / EPres[icent]; //v2 = v2_raw / EP_resolution
                            double binCenter = (v2ptBins[ipt]+v2ptBins[ipt+1])/2.;
                            _s[v2string]->addPoint(binCenter, v2, v2ptBins[ipt+1]-binCenter, 0.);
                    }
            }


    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, Scatter2DPtr> _s;
    std::vector<double> v2ptBins = {0.25, 0.5, 7.5, 1., 1.25, 1.5, 1.75, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 8.};
    std::vector<double> v2centBins = {0., 10., 20., 30., 40., 50.};
    //@}


  };


  DECLARE_RIVET_PLUGIN(EventPlaneExample);

}
