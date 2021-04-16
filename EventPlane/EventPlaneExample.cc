// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
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

            return eventPlane;
    }

    double CalculateVn(const Particles& particles, double eventPlane, int nth, double minPt, int nBins=-1)
    {
        if (nBins == -1) nBins = sqrt(particles.size());

        std::vector<double> phiBins(nBins,0.0);

        double binWidth = 2*M_PI/nBins;

        for(const Particle &p : particles)
        {
            if (p.pT()/GeV < minPt)
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
            Vn += phiBins[i]*cos(nth*((i*binWidth) + binCenter));
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
            chi = (( con*chi*exp(-chi*chi/2.))*( BesselI0(chi *chi /2.) + BesselI1(chi*chi /2.) )  < res) ? chi + delta : chi - delta ;
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
        double con = sqrt(M_PI)/2.;

        double funcRes = (con*chi*exp(-chi*chi/2.))*( BesselI0(chi *chi /2.) + BesselI1(chi*chi /2.));

        return funcRes;
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      const FinalState fs(Cuts::abseta < 0.5 && Cuts::pT > 0.150*GeV);
      declare(fs, "fs");



      /*
      Cut cutsBBCP = Cuts::eta < 0.5 && Cuts::eta > 0. && Cuts::pT > 0.150*GeV;

      EventPlane bbcP(2,cutsBBCP);
      declare(bbcP, "bbcP");

      Cut cutsBBCN = Cuts::eta > -0.5 && Cuts::eta < 0. && Cuts::pT > 0.150*GeV;

      EventPlane bbcN(2,cutsBBCN);
      declare(bbcN, "bbcN");
      */
      /*
      cout << "CalculateChi: " << CalculateChi(0.2) << endl;

      double chi = CalculateChi(0.2);

      cout << "Resolution: " << Resolution(chi) << endl;

      cout << "CalculateChi: " << CalculateChi(0.2/sqrt(2.)) << endl;

      chi = CalculateChi(0.2/sqrt(2.));

      cout << "Resolution: " << Resolution(chi)*sqrt(2.) << endl;
      */
      book(hPhi, "hPhi", 36, 0, 2*M_PI);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const FinalState& fs = apply<FinalState>(event, "fs");

      Particles particles = fs.particles();

      /*
      EventPlane bbcP = apply<EventPlane>(event,"bbcP");
      double epBBCP = bbcP.getNthOrderAngle(2);

      EventPlane bbcN = apply<EventPlane>(event,"bbcN");
      double epBBCN = bbcN.getNthOrderAngle(2);

      cout << "epBBCP: " << epBBCP << " epBBCN: " << epBBCN << endl;
      */

      //cout << "nparticles: " << particles.size() << endl;

      for(Particle& p : particles)
      {
              hPhi->fill(p.phi());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {



    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr hPhi;
    //@}


  };


  DECLARE_RIVET_PLUGIN(EventPlaneExample);

}
