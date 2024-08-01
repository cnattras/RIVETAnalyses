// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "../Centralities/RHICCentrality.hh"
#include <stdio.h>

namespace Rivet {


  class Correlator {

    private:
      std::vector<int> _indices;
      string _collSystemAndEnergy;
      pair<double,double> _centrality;
      pair<double,double> _triggerRange;
      pair<double,double> _associatedRange;
	pair<int,int> _RxnPlaneAngleRange;
      vector<int> _pid;
	int _RxnPlaneAngle = 0;
	int _eventPlaneMethod = 0;
      bool _noCentrality = false;
      bool _noAssoc = false;
	Histo1DPtr _deltaPhi;
	CounterPtr _counter;
        CounterPtr _cTriggers;


    public:

	vector<int> findIndicies(){return _indices; }

	Correlator(int index0, int index1, int index2) { _indices = {index0, index1, index2};}
	Correlator(int index0, int index1) {_indices = {index0, index1};}
	Correlator(int index0) {_indices = {index0};}
	Correlator(std::vector<int> vindex) {_indices = vindex;}

      void SetCollSystemAndEnergy(string s){ _collSystemAndEnergy = s; }
      void SetCentrality(double cmin, double cmax){ _centrality = make_pair(cmin, cmax); }
      void SetNoCentrality(){ _noCentrality = true; }
      void SetNoAssoc(){ _noAssoc = true; }
      void SetTriggerRange(double tmin, double tmax){ _triggerRange = make_pair(tmin, tmax); }
      void SetAssociatedRange(double amin, double amax){ _associatedRange = make_pair(amin, amax); }
      void SetPID(std::initializer_list<int> pid){ _pid = pid; }
	void SetCorrelationFunction(Histo1DPtr cf){ _deltaPhi = cf; }
	void SetCounter(CounterPtr c){ _counter = c; }
	void SetRxnPlaneAngle(int rxnMin, int rxnMax){ _RxnPlaneAngleRange = make_pair(rxnMin,rxnMax); }
	void SetEventPlaneMethod (int eventPlane) { _eventPlaneMethod = eventPlane; }
        void SetTriggerCounter(CounterPtr c){ _cTriggers = c; }

      string GetCollSystemAndEnergy(){ return _collSystemAndEnergy; }
      pair<double,double> GetCentrality(){ return _centrality; }
      double GetCentralityMin(){ return _centrality.first; }
      double GetCentralityMax(){ return _centrality.second; }
      pair<double,double> GetTriggerRange(){ return _triggerRange; }
      double GetTriggerRangeMin(){ return _triggerRange.first; }
      double GetTriggerRangeMax(){ return _triggerRange.second; }
      pair<double,double> GetAssociatedRange(){ return _associatedRange; }
      double GetAssociatedRangeMin(){ return _associatedRange.first; }
      double GetAssociatedRangeMax(){ return _associatedRange.second; }
	pair<int,int> GetRxnPlaneAngle() { return _RxnPlaneAngleRange; }
	int GetEventPlaneMethod() { return _eventPlaneMethod; }
      vector<int> GetPID(){ return _pid; }
	double GetWeight(){ return _counter->sumW(); }
	Histo1DPtr GetCorrelationFunction(){ return _deltaPhi; }
	CounterPtr GetCounter(){ return _counter; }

	double GetDeltaPhi(Particle pAssoc, Particle pTrig)
	{
		double dPhi = deltaPhi(pTrig, pAssoc, true);
		//this does NOT rotate the delta phi to be in a given range

	        if(dPhi < -M_PI/2.)
	        {
	            dPhi += 2.*M_PI;
	        }
	        else if(dPhi > 3.*M_PI/2.)
	        {
	            dPhi -= 2*M_PI;
	        }

	        return dPhi;
	}


	void AddCorrelation(Particle pTrig, Particle pAssoc)
	{
		double dPhi = GetDeltaPhi(pTrig, pAssoc);
		_deltaPhi->fill(dPhi);
	}

	void AddWeight()
	{
		_counter->fill();
	}

	int GetIndex(int i){ return _indices[i]; }
      string GetFullIndex()
      {
          string fullIndex = "";
  	  for(int index : _indices)
 	  {
		fullIndex += to_string(index);
	  }
          return fullIndex;
      }

      void AddTrigger()
	{
		_cTriggers->fill();
	}

      void Normalize(double weight = 1.)
	{
		if(_cTriggers->effNumEntries()*_counter->sumW() > 0) _deltaPhi->scaleW((weight*_counter->effNumEntries())/(_cTriggers->effNumEntries()*_counter->sumW()));
	}

      bool CheckCollSystemAndEnergy(string s){ return _collSystemAndEnergy.compare(s) == 0 ? true : false; }
      bool CheckCentrality(double cent){ return ((cent>_centrality.first && cent<_centrality.second) || _noCentrality == true) ? true : false; }
      bool CheckTriggerRange(double tpt){ return (tpt>_triggerRange.first && tpt<_triggerRange.second) ? true : false; }
      bool CheckAssociatedRange(double apt){ return ((apt>_associatedRange.first && apt<_associatedRange.second) || _noAssoc == true) ? true : false; }
      bool CheckAssociatedRangeMaxTrigger(double apt, double tpt){ return (apt>_associatedRange.first && apt<tpt) ? true : false; }
      bool CheckPID(std::initializer_list<int> pid)
      {
          bool inList = false;

          for(int id : pid)
          {
              auto it = std::find(_pid.begin(), _pid.end(), id);

              if(it != _pid.end())
              {
                  inList = true;
                  break;
              }
          }
          return inList;
      }

      bool CheckConditions(string s, double cent, double tpt, double apt)
      {
        if(!CheckConditions(s, cent, tpt)) return false;
        if(!CheckAssociatedRange(apt)) return false;

        return true;

      }

      bool CheckConditions(string s, double cent, double tpt)
      {
        if(!CheckConditions(s, cent)) return false;
        if(!CheckTriggerRange(tpt)) return false;

        return true;

      }

      bool CheckConditions(string s, double cent)
      {
        if(!CheckCollSystemAndEnergy(s)) return false;
        if(!CheckCentrality(cent)) return false;

        return true;

      }

      bool CheckConditionsMaxTrigger(string s, double cent, double tpt, double apt)
      {
        if(!CheckCollSystemAndEnergy(s)) return false;
        if(!CheckCentrality(cent)) return false;
        if(!CheckTriggerRange(tpt)) return false;
        if(!CheckAssociatedRangeMaxTrigger(apt,tpt)) return false;

        return true;

      }
     };


  /// @brief Add a short analysis description here
  class PHENIX_2019_I1658594 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2019_I1658594);


    /// @name Analysis methods
    //@{

	Histo1DPtr SubtractBackgroundZYAM(Histo1DPtr histo)
	{
		YODA::Histo1D hist = *histo;

		double minValue = sqrt(-2);
	        double binWidth = 0.;
	        int minValueEntries = 0.;

        	for(auto &bin : hist.bins())
	        {
	            if(std::isnan(minValue))
	            {
	                minValue = bin.sumW();
	                binWidth = bin.width();
	                minValueEntries = bin.numEntries();
	            }
	            if(bin.sumW()/bin.width() < minValue/binWidth)
	            {
	                minValue = bin.sumW();
	                binWidth = bin.width();
	                minValueEntries = bin.numEntries();
	            }
        	}

		//if the minValue is already zero, do not subtract
		if(minValue == 0 || minValueEntries==0) return histo;

	        hist.reset();

	        for(auto &bin : hist.bins())
	        {
	            bin.fillBin((minValue*bin.width())/(minValueEntries*binWidth), minValueEntries);
	        }

	        *histo = YODA::subtract(*histo, hist);

	        return histo;
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

            double eventPlane = mapAngle0To2Pi((1./n)*atan2(QIn,QRn));

            return eventPlane;
    }

    void FillVn(Profile1DPtr vnHisto, const Particles& particles, double eventPlane, int n)
    {
        for(const Particle &p : particles)
        {
            vnHisto->fill(p.pT()/GeV, cos(n*(p.phi() - eventPlane)));
        }
    }

    double CalculateChi(double res)
    {
        //Implementation from A. M. Poskanzer and S. A. Voloshin, Phys. Rev. C 58, 1671 – Published 1 September 1998
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

    double BesselIn(int n, double x)
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

    /// Book histograms and initialise projections before the run
    void init() {

      const ChargedFinalState cfs(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV);
      declare(cfs, "CFS");

      declareCentrality(RHICCentrality("PHENIX"), "RHIC_2019_CentralityCalibration:exp=PHENIX", "CMULT", "CMULT");

      const FinalState RxP(Cuts::abseta > 1. && Cuts::abseta < 2.8);
      declare(RxP, "RxP");

      const FinalState RxPPos(Cuts::eta > 1. && Cuts::eta < 2.8);
      declare(RxPPos, "RxPPos");

      const FinalState RxPNeg(Cuts::eta < -1. && Cuts::eta > -2.8);
      declare(RxPNeg, "RxPNeg");

      book(_p["RxPcosPosv2"], "RxPcosPosv2", 10, 0., 10.);
      book(_p["RxPcosPosv3"], "RxPcosPosv3", 10, 0., 10.);
      book(_p["RxPcosPosv4"], "RxPcosPosv4", 10, 0., 10.);

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      //const FinalState fs(Cuts::abseta < 4.9);


	//initialize iterators and varibles
	int minCent=0, maxCent=0, minPlane=0, maxPlane=0, a=0, e=0, iterator=1;
	float min_pT=0, max_pT=0, min_pA=0, max_pA=0, minV=0, maxV=0, b=0, c=0, d=0;
	char corrName[200], bookName[200], corrNameTrigger[200];
/*
	//fig 6
	minCent=0, maxCent=40, minPlane=2, maxPlane=5;
	for(e=minPlane; e<=maxPlane; e++){
		for(a=minCent; a<=maxCent; a+=10){
			if(e!=5){
				snprintf(corrName,200,"CounterFig6EventPlane%iCent%iTo%i",e,a,a+10);
                                snprintf(corrNameTrigger,200,"CounterFig6EventPlane%iCent%iTo%i%s",e,a,a+10,"Triggers");
				snprintf(bookName,200,"Fig6EventPlane%iCent%iTo%i",e,a,a+10);
			}
			else{
				snprintf(corrName,200,"CounterFig6EventPlane%i{Psi2}Cent%iTo%i",e,a,a+10);
                                snprintf(corrNameTrigger,200,"CounterFig6EventPlane%i{Psi2}Cent%iTo%i%s",e,a,a+10,"Triggers");
				snprintf(bookName,200,"Fig6EventPlane%i{Psi2}Cent%iTo%i",e,a,a+10);
			}
			book(_h[bookName],1,1,iterator);
			book(_c[corrName], corrName);
                        book(_c[corrNameTrigger], corrNameTrigger);

			Correlator corrFig6(a,e);
                        corrFig6.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig6.SetCentrality(a,a+10);
                        corrFig6.SetEventPlaneMethod(e);
                        corrFig6.SetCorrelationFunction(_h[bookName]);
                        corrFig6.SetCounter(_c[corrName]);
                        corrFig6.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig6);
                        iterator++;
		}
	}
	iterator=1;
*/
	//fig 12
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=0.5, max_pA=4, minV=1, maxV=1;
	b=2;
	for(a=minCent; a<=maxCent; a+=10){
		for(c=min_pA; c<=max_pA; c*=2){
			snprintf(corrName,200,"CounterFig12Cent%iTo%iPtA%2.1fTo%2.1f",a,a+10,c,c*2);
                        snprintf(corrNameTrigger,200,"CounterFig12Cent%iTo%iPtA%2.1fTo%2.1f%s",a,a+10,c,c*2, "Triggers");
			snprintf(bookName,200,"Fig12Cent%iTo%iPtA%2.1fTo%2.1f",a,a+10,c,c*2);
			book(_h[bookName],2,1,iterator);
			book(_c[corrName], corrName);
                        book(_c[corrNameTrigger], corrNameTrigger);

			Correlator corrFig12(a,(int) c*10);
	 		corrFig12.SetCollSystemAndEnergy("AuAu200GeV");
	 		corrFig12.SetCentrality(a,a+10);
	 		corrFig12.SetTriggerRange(b,(b==4)?10:b*2);
	 		corrFig12.SetAssociatedRange(c,(c==4)?10:c*2);
	 		corrFig12.SetCorrelationFunction(_h[bookName]);
	 		corrFig12.SetCounter(_c[corrName]);
                        corrFig12.SetTriggerCounter(_c[corrNameTrigger]);
			Correlators.push_back(corrFig12);
			iterator++;
		}
	}
	iterator=1;


	//fig 15
	minCent=0, maxCent=40, min_pT=1, max_pT=2, min_pA=0.5, max_pA=2, minV=1, maxV=1;
	for(a=minCent; a<=maxCent; a+=10){
		for(b=min_pT; b<=max_pT; b*=2){
			for(c=min_pA; c<=max_pA && c!=b ; c*=2){
				snprintf(corrName,200,"CounterFig15Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
					a,a+10,b,b*2,c,c*2);
				snprintf(bookName,200,"Fig15Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
					a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
				 snprintf(corrNameTrigger,200,"CounterFig15Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f%s",
                                        a,a+10,b,b*2,c,c*2, "Trigger");
				book(_h[bookName],3,1,iterator);
				book(_c[corrName], corrName);
				book(_c[corrNameTrigger], corrNameTrigger);

				Correlator corrFig15(a,(int)b*10,(int)c*10);
				corrFig15.SetCollSystemAndEnergy("AuAu200GeV");
				corrFig15.SetCentrality(a,a+10);
				corrFig15.SetTriggerRange(b,(b==4)?10:b*2);
				corrFig15.SetAssociatedRange(c,(c==4)?10:c*2);
				corrFig15.SetCorrelationFunction(_h[bookName]);
				corrFig15.SetCounter(_c[corrName]);
				corrFig15.SetTriggerCounter(_c[corrNameTrigger]);
				Correlators.push_back(corrFig15);
				iterator++;
			}
		}
	}
	iterator=1;

	//Fig 18
	minCent=0, maxCent=40, min_pT=2, max_pT=4, min_pA=1, max_pA=2, minV=1, maxV=1;
	for(int sides=0; sides<2; sides++){
		for(a=minCent; a<=maxCent; a+=10){
                	for(b=min_pT; b<=max_pT; b*=2){
                        	for((b==4) ? c=min_pA*2 : c=min_pA; c<=max_pA; c*=2){
                                	if(sides==0){
						snprintf(bookName,200,
							"NearSideFig18Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
							a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrName,200,
                                                	"CounterNearSideFig18Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
	                                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrNameTrigger,200,
                                                        "CounterNearSideFig18Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f%s",
                                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2, "Triggers");
					}
					else{
						snprintf(bookName,200,
							"FarSideFig18Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
							a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrName,200,
	                                                "CounterFarSideFig18Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
        	                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrNameTrigger,200,
                                                        "CounterFarSideFig18Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f%s",
                                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2, "Triggers");
					}
                                	book(_h[bookName],4,1,iterator);
					book(_c[corrName], corrName);
					book(_c[corrNameTrigger], corrNameTrigger);

	                                Correlator corrFig18(a,(int)b*10,(int)c*10);
        	                        corrFig18.SetCollSystemAndEnergy("AuAu200GeV");
                	                corrFig18.SetCentrality(a,a+10);
	                                if(b==0.5 || b==1 || b==2) corrFig18.SetTriggerRange(b,b*2);
        	                        else corrFig18.SetTriggerRange(b,10);
                	                if(c==0.5 || c==1 || c==2) corrFig18.SetAssociatedRange(c,c*2);
                        	        else corrFig18.SetAssociatedRange(c,10);
					corrFig18.SetCorrelationFunction(_h[bookName]);
	                                corrFig18.SetCounter(_c[corrName]);
					corrFig18.SetTriggerCounter(_c[corrNameTrigger]);
        	                        Correlators.push_back(corrFig18);
                	                iterator++;
				}
			}
		}
	}
	iterator=1;

	//Fig 20
	minCent=0, maxCent=40, min_pT=2, max_pT=4, min_pA=1, max_pA=2, minV=1, maxV=1;
        for(int sides=0; sides<2; sides++){
                for(a=minCent; a<=maxCent; a+=10){
                        for(b=min_pT; b<=max_pT; b*=2){
                                for((b==4) ? c=min_pA*2 : c=min_pA; c<=max_pA; c*=2){
                                	if(sides==0){
						snprintf(bookName,200,
							"NearSideFig20Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
							a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrName,200,
	                                                "CounterNearSideFig20Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
        	                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						 snprintf(corrNameTrigger,200,
                                                        "CounterNearSideFig20Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f%s",
                                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2, "Trigger");
					}
                                	else{
						snprintf(bookName,200,
							"FarSideFig20Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
							a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						snprintf(corrName,200,
	                                                "CounterFarSideFig20Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
        	                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
						 snprintf(corrNameTrigger,200,
                                                        "CounterFarSideFig20Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f%s",
                                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2, "Trigger");

					}
                                        book(_h[bookName],5,1,iterator);
					book(_c[corrName], corrName);
					book(_c[corrNameTrigger], corrNameTrigger);

                                        Correlator corrFig20(a,(int)b*10,(int)c*10);
                                        corrFig20.SetCollSystemAndEnergy("AuAu200GeV");
                                        corrFig20.SetCentrality(a,a+10);
                                        if(b==0.5 || b==1 || b==2) corrFig20.SetTriggerRange(b,b*2);
                                        else corrFig20.SetTriggerRange(b,10);
                                        if(c==0.5 || c==1 || c==2) corrFig20.SetAssociatedRange(c,c*2);
                                        else corrFig20.SetAssociatedRange(c,10);
					corrFig20.SetCorrelationFunction(_h[bookName]);
                                        corrFig20.SetCounter(_c[corrName]);
					corrFig20.SetTriggerCounter(_c[corrNameTrigger]);
                                        Correlators.push_back(corrFig20);
                                        iterator++;
				}
			}
		}
	}
	iterator=1;

	//Fig 21
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=0.5, max_pA=4, minV=1, maxV=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(b=min_pT; b<=max_pT; b*=2){
                        for(c=min_pA; c<=max_pA; c*=2){
                                snprintf(corrName,200,"CounterFig21Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
				snprintf(bookName,200,
					"Fig21Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
					a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2);
				snprintf(corrNameTrigger,200,"CounterFig21Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f%s",
                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2, "Triggers");
                                book(_h[bookName],6,1,iterator);
				book(_c[corrName], corrName);
				book(_c[corrNameTrigger], corrNameTrigger);

                                Correlator corrFig21(a,(int)c*10);
                                corrFig21.SetCollSystemAndEnergy("AuAu200GeV");
                                corrFig21.SetCentrality(a,a+10);
				corrFig21.SetTriggerRange(b,(b==4)?10:b*2);
				corrFig21.SetAssociatedRange(c,(c==4)?10:c*2);
				corrFig21.SetCorrelationFunction(_h[bookName]);
                                corrFig21.SetCounter(_c[corrName]);
				corrFig21.SetTriggerCounter(_c[corrNameTrigger]);
                                Correlators.push_back(corrFig21);
				iterator++;
			}
		}
	}
	iterator=1;

	//Fig22
	minCent=0, maxCent=40, min_pT=1, max_pT=2, min_pA=0.5, max_pA=2, minV=1, maxV=1;
	for(a=minCent; a<=maxCent; a+=10){
        	for(b=min_pT; b<=max_pT; b*=2){
        		for(c=min_pA; c<=max_pA && !(c>b); c*=2){
                                snprintf(corrName,200,
					"CounterFig22Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2);
                                snprintf(bookName,200,
                                        "Fig22Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f",
                                        a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2);
				snprintf(corrNameTrigger,200,
                                        "CounterFig22Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1f%s",
                                        a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2, "Triggers");
                                book(_h[bookName],7,1,iterator);
				book(_c[corrName], corrName);
				book(_c[corrNameTrigger], corrNameTrigger);

                                Correlator corrFig22(a,(int)c*10);
                                corrFig22.SetCollSystemAndEnergy("AuAu200GeV");
                                corrFig22.SetCentrality(a,a+10);
                                corrFig22.SetTriggerRange(b,(b==4)?10:b*2);
                                corrFig22.SetAssociatedRange(c,(c==4)?10:c*2);
				corrFig22.SetCorrelationFunction(_h[bookName]);
                                corrFig22.SetCounter(_c[corrName]);
				corrFig22.SetTriggerCounter(_c[corrNameTrigger]);
                                Correlators.push_back(corrFig22);
                                iterator++;
                        }
                }
        }
        iterator=1;

	//Fig23
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
	for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
			snprintf(corrName,200,
				"CounterFig23Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig23Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig23Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "Trigger");
                        book(_h[bookName],8,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig23(a,(int)d);
                        corrFig23.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig23.SetCentrality(a,a+10);
                        corrFig23.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig23.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig23.SetRxnPlaneAngle(d,d+1);
			corrFig23.SetCorrelationFunction(_h[bookName]);
			corrFig23.SetCounter(_c[corrName]);
			corrFig23.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig23);
                        iterator++;
                }
        }
        iterator=1;



	//Fig24
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig24Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig24Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(corrNameTrigger,200,
                                "CounterFig24Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "Trigger");
                        book(_h[bookName],9,1,iterator);
                        book(_c[corrName], corrName);
                        book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig24(a,(int)d);
                        corrFig24.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig24.SetCentrality(a,a+10);
                        corrFig24.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig24.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig24.SetRxnPlaneAngle(d,d+1);
                        corrFig24.SetCorrelationFunction(_h[bookName]);
                        corrFig24.SetCounter(_c[corrName]);
                        corrFig24.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig24);
                        iterator++;
                }
        }
        iterator=1;

	//Fig25
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig25Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig25Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(corrNameTrigger,200,
                                "CounterFig25Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "Trigger");
                        book(_h[bookName],10,1,iterator);
                        book(_c[corrName], corrName);
                        book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig25(a,(int)d);
                        corrFig25.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig25.SetCentrality(a,a+10);
                        corrFig25.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig25.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig25.SetRxnPlaneAngle(d,d+1);
                        corrFig25.SetCorrelationFunction(_h[bookName]);
                        corrFig25.SetCounter(_c[corrName]);
                        corrFig25.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig25);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 26
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig26Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig26Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig26Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],11,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig26(a,(int)d);
                        corrFig26.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig26.SetCentrality(a,a+10);
                        corrFig26.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig26.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig26.SetRxnPlaneAngle(d,d+1);
			corrFig26.SetCorrelationFunction(_h[bookName]);
                        corrFig26.SetCounter(_c[corrName]);
			corrFig26.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig26);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 27
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig27Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig27Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig27Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],12,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig27(a,(int)d);
                        corrFig27.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig27.SetCentrality(a,a+10);
                        corrFig27.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig27.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig27.SetRxnPlaneAngle(d,d+1);
			corrFig27.SetCorrelationFunction(_h[bookName]);
                        corrFig27.SetCounter(_c[corrName]);
			corrFig27.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig27);
                        iterator++;
                }
        }
        iterator=1;


	//Fig 28
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig28Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig28Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig28Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],13,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig28(a,(int)d);
                        corrFig28.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig28.SetCentrality(a,a+10);
                        corrFig28.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig28.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig28.SetRxnPlaneAngle(d,d+1);
			corrFig28.SetCorrelationFunction(_h[bookName]);
                        corrFig28.SetCounter(_c[corrName]);
			corrFig28.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig28);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 29
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig29Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig29Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			 snprintf(corrNameTrigger,200,
                                "CounterFig29Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1,"_Triggers");
                        book(_h[bookName],14,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig29(a,(int)d);
                        corrFig29.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig29.SetCentrality(a,a+10);
                        corrFig29.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig29.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig29.SetRxnPlaneAngle(d,d+1);
			corrFig29.SetCorrelationFunction(_h[bookName]);
                        corrFig29.SetCounter(_c[corrName]);
			corrFig29.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig29);
                        iterator++;
                }
        }
        iterator=1;


	//Fig 30
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig30Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig30Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig30Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],15,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig30(a,(int)d);
                        corrFig30.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig30.SetCentrality(a,a+10);
                        corrFig30.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig30.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig30.SetRxnPlaneAngle(d,d+1);
			corrFig30.SetCorrelationFunction(_h[bookName]);
                        corrFig30.SetCounter(_c[corrName]);
			corrFig30.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig30);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 31
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig31Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig31Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			 snprintf(corrNameTrigger,200,
                                "CounterFig31Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Trigger");
                        book(_h[bookName],16,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig31(a,(int)d);
                        corrFig31.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig31.SetCentrality(a,a+10);
                        corrFig31.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig31.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig31.SetRxnPlaneAngle(d,d+1);
			corrFig31.SetCorrelationFunction(_h[bookName]);
                        corrFig31.SetCounter(_c[corrName]);
			corrFig31.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig31);
                        iterator++;
                }
        }
        iterator=1;


	//Fig 32
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig32Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig32Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			 snprintf(corrNameTrigger,200,
                                "CounterFig32Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Trigger");
                        book(_h[bookName],17,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig32(a,(int)d);
                        corrFig32.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig32.SetCentrality(a,a+10);
                        corrFig32.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig32.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig32.SetRxnPlaneAngle(d,d+1);
			corrFig32.SetCorrelationFunction(_h[bookName]);
                        corrFig32.SetCounter(_c[corrName]);
			corrFig32.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig32);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 33
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig33Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig33Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,d,d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig33Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Trigger");
                        book(_h[bookName],18,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig33(a,(int)d);
                        corrFig33.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig33.SetCentrality(a,a+10);
                        corrFig33.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig33.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig33.SetRxnPlaneAngle(d,d+1);
			corrFig33.SetCorrelationFunction(_h[bookName]);
                        corrFig33.SetCounter(_c[corrName]);
			corrFig33.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig33);
                        iterator++;
                }
        }
        iterator=1;

	//Fig 34
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig34Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig34Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			 snprintf(corrNameTrigger,200,
                                "CounterFig34Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Trigger");
                        book(_h[bookName],19,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig34(a,(int)d);
                        corrFig34.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig34.SetCentrality(a,a+10);
                        corrFig34.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig34.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig34.SetRxnPlaneAngle(d,d+1);
			corrFig34.SetCorrelationFunction(_h[bookName]);
                        corrFig34.SetCounter(_c[corrName]);
			corrFig34.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig34);
                        iterator++;
                }
        }
        iterator=1;

	//Psi2_Cor2-4x1-2GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(35)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig(35)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig(35)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],20,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig35(a,(int)d);
                        corrFig35.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig35.SetCentrality(a,a+10);
                        corrFig35.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig35.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig35.SetRxnPlaneAngle(d,d+1);
			corrFig35.SetCorrelationFunction(_h[bookName]);
                        corrFig35.SetCounter(_c[corrName]);
			corrFig35.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig35);
                        iterator++;
                }
        }
        iterator=1;


	//Psi2_Cor2-4x2-4GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(36)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig36Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig(36)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],21,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig36(a,(int)d);
                        corrFig36.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig36.SetCentrality(a,a+10);
                        corrFig36.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig36.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig36.SetRxnPlaneAngle(d,d+1);
			corrFig36.SetCorrelationFunction(_h[bookName]);
                        corrFig36.SetCounter(_c[corrName]);
			corrFig36.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig36);
                        iterator++;
                }
        }
        iterator=1;

	//Psi2_Cor4-10x2-4GeV
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(37)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig(37)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			 snprintf(corrNameTrigger,200,
                                "CounterFig(37)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],22,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig37(a,(int)d);
                        corrFig37.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig37.SetCentrality(a,a+10);
                        corrFig37.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig37.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig37.SetRxnPlaneAngle(d,d+1);
			corrFig37.SetCorrelationFunction(_h[bookName]);
                        corrFig37.SetCounter(_c[corrName]);
			corrFig37.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig37);
                        iterator++;
                }
        }
        iterator=1;


	//Psi2_PTY2-4x1-2GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(38)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig(38)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig(38)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],23,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig38(a,(int)d);
                        corrFig38.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig38.SetCentrality(a,a+10);
                        corrFig38.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig38.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig38.SetRxnPlaneAngle(d,d+1);
			corrFig38.SetCorrelationFunction(_h[bookName]);
                        corrFig38.SetCounter(_c[corrName]);
			corrFig38.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig38);
                        iterator++;
                }
        }
        iterator=1;

	//Psi2_PTY2-4x2-4GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(39)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig(39)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			 snprintf(corrNameTrigger,200,
                                "CounterFig(39)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],24,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig39(a,(int)d);
                        corrFig39.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig39.SetCentrality(a,a+10);
                        corrFig39.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig39.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig39.SetRxnPlaneAngle(d,d+1);
			corrFig39.SetCorrelationFunction(_h[bookName]);
                        corrFig39.SetCounter(_c[corrName]);
			corrFig39.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig39);
                        iterator++;
                }
        }
        iterator=1;


	//Psi2_PTY4-10x2-4GeV
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(40)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%iTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig(40)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig(40)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%iTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],25,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);


                        Correlator corrFig40(a,(int)d);
                        corrFig40.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig40.SetCentrality(a,a+10);
                        corrFig40.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig40.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig40.SetRxnPlaneAngle(d,d+1);
			corrFig40.SetCorrelationFunction(_h[bookName]);
                        corrFig40.SetCounter(_c[corrName]);
			corrFig40.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig40);
                        iterator++;
                }
        }
        iterator=1;

	//Psi3_Cor2-4x1-2GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(41)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig(41)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			 snprintf(corrNameTrigger,200,
                                "CounterFig(41)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],26,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig41(a,(int)d);
                        corrFig41.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig41.SetCentrality(a,a+10);
                        corrFig41.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig41.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig41.SetRxnPlaneAngle(d,d+1);
			corrFig41.SetCorrelationFunction(_h[bookName]);
                        corrFig41.SetCounter(_c[corrName]);
			corrFig41.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig41);
                        iterator++;
                }
        }
        iterator=1;

	//Psi3_Cor2-4x2-4GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(42)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig(42)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig(42)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],27,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig42(a,(int)d);
                        corrFig42.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig42.SetCentrality(a,a+10);
                        corrFig42.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig42.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig42.SetRxnPlaneAngle(d,d+1);
			corrFig42.SetCorrelationFunction(_h[bookName]);
                        corrFig42.SetCounter(_c[corrName]);
			corrFig42.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig42);
                        iterator++;
                }
        }
        iterator=1;

	//Psi3_Cor4-10x2-4GeV
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(43)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig(43)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig(43)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],28,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig43(a,(int)d);
                        corrFig43.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig43.SetCentrality(a,a+10);
                        corrFig43.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig43.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig43.SetRxnPlaneAngle(d,d+1);
			corrFig43.SetCorrelationFunction(_h[bookName]);
                        corrFig43.SetCounter(_c[corrName]);
			corrFig43.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig43);
                        iterator++;
                }
        }
        iterator=1;


	//Psi3_PTY2-4x1-2GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=1, max_pA=1, minV=-4, maxV=-1;
        b=2,c=1;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(44)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%iTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig(44)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig(44)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%iTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],29,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig44(a,(int)d);
                        corrFig44.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig44.SetCentrality(a,a+10);
                        corrFig44.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig44.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig44.SetRxnPlaneAngle(d,d+1);
			corrFig44.SetCorrelationFunction(_h[bookName]);
                        corrFig44.SetCounter(_c[corrName]);
			corrFig44.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig44);
                        iterator++;
                }
        }
        iterator=1;

	//Psi3_PTY2-4x2-4GeV
	minCent=0, maxCent=40, min_pT=2, max_pT=2, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=2,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(45)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig(45)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig(45)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],30,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig45(a,(int)d);
                        corrFig45.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig45.SetCentrality(a,a+10);
                        corrFig45.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig45.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig45.SetRxnPlaneAngle(d,d+1);
			corrFig45.SetCorrelationFunction(_h[bookName]);
                        corrFig45.SetCounter(_c[corrName]);
			corrFig45.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig45);
                        iterator++;
                }
        }
        iterator=1;


	//Psi3_PTY4-10x2-4GeV
	minCent=0, maxCent=40, min_pT=4, max_pT=4, min_pA=2, max_pA=2, minV=-4, maxV=-1;
        b=4,c=2;
        for(a=minCent; a<=maxCent; a+=10){
                for(d=minV; d<=maxV; d++){
                        snprintf(corrName,200,
                                "CounterFig(46)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1);
                        snprintf(bookName,200,
                                "Fig(46)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi",
                                a,a+10,b,(b==4)?10:b*2,c,(b==4)?10:c*2,(int)d,(int)d+1);
			snprintf(corrNameTrigger,200,
                                "CounterFig(46)Cent%iTo%iPtT%2.1fTo%2.1fPtA%2.1fTo%2.1fAngle%i/piTo%i/pi%s",
                                a,a+10,b,(b==4)?10:b*2,c,(c==4)?10:c*2,(int)d,(int)d+1, "_Triggers");
                        book(_h[bookName],31,1,iterator);
			book(_c[corrName], corrName);
			book(_c[corrNameTrigger], corrNameTrigger);

                        Correlator corrFig46(a,(int)d);
                        corrFig46.SetCollSystemAndEnergy("AuAu200GeV");
                        corrFig46.SetCentrality(a,a+10);
                        corrFig46.SetTriggerRange(b,(b==4)?10:b*2);
                        corrFig46.SetAssociatedRange(c,(c==4)?10:c*2);
                        corrFig46.SetRxnPlaneAngle(d,d+1);
			corrFig46.SetCorrelationFunction(_h[bookName]);
                        corrFig46.SetCounter(_c[corrName]);
			corrFig46.SetTriggerCounter(_c[corrNameTrigger]);
                        Correlators.push_back(corrFig46);
                        iterator++;
                }
        }
        iterator=1;

        for(unsigned int icent = 0; icent < v2centBins.size()-1; icent++)
        {
                string v2string = "v2_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);
                book(_p[v2string], 1, 1, 1+icent);

                string v3string = "v3_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);
                book(_p[v3string], 1, 1, 6+icent);

                string v4string = "v4_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);
                book(_p[v4string], 1, 1, 11+icent);

                string v4ep2string = "v4ep2_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);
                book(_p[v4ep2string], 1, 1, 16+icent);
        }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const CentralityProjection& cent = apply<CentralityProjection>(event, "CMULT");
      //add calculation of reaction plane angle
      const double c = cent();
      const ParticlePair& beam = beams();
      string CollSystem = "Empty";

      //add pp collision eventually
	if (beam.first.pid() == 1000791970 && beam.second.pid() == 1000791970)
      {
          CollSystem = "AuAu200GeV";
          //if(fuzzyEquals(sqrtS()/GeV, 200*NN, 5)) CollSystem += "200GeV";
      }

      if(CollSystem == "AuAu200GeV" && c > 50)
      {
        vetoEvent;
      }

      //cout << c << endl;
      for(Correlator& corr : Correlators)
      {
        if(!corr.CheckCollSystemAndEnergy(CollSystem)) continue;
        if(!corr.CheckCentrality(c)) continue;
        corr.AddWeight();
      }

      //Correlator corr = Correlators[0];

      for(auto pTrig : cfs.particles())
      {
        for (Correlator &corr : Correlators)
        {
          if (!corr.CheckCollSystemAndEnergy(CollSystem)) continue;
          if (!corr.CheckCentrality(c)) continue;
          //cout << "hi" << '\n';
          if (!corr.CheckTriggerRange(pTrig.pT() / GeV)) continue;
            ///Add a function to check reaction plane angle and make sure that the
            //if (!corr.CheckXiRange(log(pTrig.pT()/ pAssoc.pT()))) continue;
          //somthing here is stoping the pp
          //cout << c << '\n';
          //cout << "hi" << '\n';
          corr.AddTrigger();
        }
        for (auto pAssoc : cfs.particles())
        {
          for (Correlator &corr : Correlators)
          {
            if (!corr.CheckCollSystemAndEnergy(CollSystem)) continue;
            if (!corr.CheckCentrality(c)) continue;
            if (!corr.CheckTriggerRange(pTrig.pT() / GeV)) continue;
            if (!corr.CheckAssociatedRange(pAssoc.pT() / GeV)) continue;
            ///Add a function to check reaction plane angle and make sure that the
            //if (!corr.CheckXiRange(log(pTrig.pT()/ pAssoc.pT()))) continue;
            corr.AddCorrelation(pTrig, pAssoc);
          }
        }
      }

	const FinalState& RxP = apply<FinalState>(event, "RxP");
      const FinalState& RxPPos = apply<FinalState>(event, "RxPPos");
      const FinalState& RxPNeg = apply<FinalState>(event, "RxPNeg");

	//Inner and Outer rings of the North and South sections of the RxP detector
	vector<Cut> etaRxP = {Cuts::eta > 1. && Cuts::eta < 1.5, Cuts::eta > 1.5 && Cuts::eta < 2.8, Cuts::eta < -1. && Cuts::eta > -1.5, Cuts::eta < -1.5 && Cuts::eta > -2.8};
      int nPhiSections = 12;

      double evPPosNeg = GetEventPlaneDetectorAcc(2, RxP, etaRxP, nPhiSections);

      vector<Cut> etaRxPPos = {Cuts::eta > 1. && Cuts::eta < 1.5, Cuts::eta > 1.5 && Cuts::eta < 2.8};
      vector<Cut> etaRxPNeg = {Cuts::eta < -1. && Cuts::eta > -1.5, Cuts::eta < -1.5 && Cuts::eta > -2.8};

      //event plane calcaulted with the dectector in the positive/negative absolute rapidity
      double evPPos = GetEventPlaneDetectorAcc(2, RxPPos, etaRxPPos, nPhiSections);
      double evPNeg = GetEventPlaneDetectorAcc(2, RxPNeg, etaRxPNeg, nPhiSections);

	//std::floor(double a) returns the largest integer value smaller than a
	 _p["RxPcosPosv2"]->fill(int(floor(c/10))+0.5, cos(2*(evPPos-evPNeg)));

         //EP3 and Resolution

         double evPPosNeg3 = GetEventPlaneDetectorAcc(3, RxP, etaRxP, nPhiSections);

         double evPPos3 = GetEventPlaneDetectorAcc(3, RxPPos, etaRxPPos, nPhiSections);
         double evPNeg3 = GetEventPlaneDetectorAcc(3, RxPNeg, etaRxPNeg, nPhiSections);

         _p["RxPcosPosv3"]->fill(int(floor(c/10))+0.5, cos(3*(evPPos3-evPNeg3)));

         //EP4 and Resolution

         double evPPosNeg4 = GetEventPlaneDetectorAcc(4, RxP, etaRxP, nPhiSections);

         double evPPos4 = GetEventPlaneDetectorAcc(4, RxPPos, etaRxPPos, nPhiSections);
         double evPNeg4 = GetEventPlaneDetectorAcc(4, RxPNeg, etaRxPNeg, nPhiSections);

         _p["RxPcosPosv4"]->fill(int(floor(c/10))+0.5, cos(4*(evPPos4-evPNeg4)));


         Particles particles = cfs.particles();

         string v2string = "v2_cent" + Form(floor(c/10)*10., 0) + Form((floor(c/10)*10.)+10., 0);
         FillVn(_p[v2string], particles, evPPosNeg, 2);

         string v3string = "v3_cent" + Form(floor(c/10)*10., 0) + Form((floor(c/10)*10.)+10., 0);
         FillVn(_p[v3string], particles, evPPosNeg3, 3);

         string v4string = "v4_cent" + Form(floor(c/10)*10., 0) + Form((floor(c/10)*10.)+10., 0);
         FillVn(_p[v4string], particles, evPPosNeg4, 4);

         string v4ep2string = "v4ep2_cent" + Form(floor(c/10)*10., 0) + Form((floor(c/10)*10.)+10., 0);
         FillVn(_p[v4ep2string], particles, evPPosNeg, 4);



    }


    /// Normalise histograms etc., after the run
    void finalize() {
	int i=1;
	for(Correlator& corr : Correlators)
	{
		corr.Normalize();
		Histo1DPtr h = corr.GetCorrelationFunction();
		if((i>=11&&i<=60)||(i>=161&&i<=280)||(i>=401&&i<=460)||(i>=521&&i<=580)){
	  	   //cout << i << " " << corr.GetTriggerRange() << "x" << corr.GetAssociatedRange() << " " << corr.findIndicies() << '\n';
			h = SubtractBackgroundZYAM(h);
		}
		i++;
	}

	int centBin = 0;

            std::vector<double> EPres(5, 0.);

            for(auto bin : _p["RxPcosPosv2"]->bins())
            {
		if(bin.numEntries() > 0)
                    {
                         double RxPPosRes;
                         if(bin.mean()>0)
                         {
                                 RxPPosRes = sqrt(bin.mean());
                                 double chiRxPPos = CalculateChi(RxPPosRes);
        			 double res = Resolution(sqrt(2)*chiRxPPos);
                                 EPres[centBin] = res;
                         }
                         else EPres[centBin] = 0;

                    }
                    centBin++;
            }

            centBin = 0;
            std::vector<double> EPres3(5, 0.);

            for(auto bin : _p["RxPcosPosv3"]->bins())
            {
		if(bin.numEntries() > 0)
                    {
                         double RxPPosRes;
			 if(bin.mean()>0)
                         {

                                 RxPPosRes = sqrt(bin.mean());
                                 double chiRxPPos = CalculateChi(RxPPosRes);
        			 double res = Resolution(sqrt(2)*chiRxPPos);
                                 EPres3[centBin] = res;
                         }
                         else EPres3[centBin] = 0;

                    }
                    centBin++;
            }

            centBin = 0;
            std::vector<double> EPres4(5, 0.);

            for(auto bin : _p["RxPcosPosv4"]->bins())
            {
		if(bin.numEntries() > 0)
                    {
                        double RxPPosRes;
			if (bin.mean()>0)
                        {
                                RxPPosRes = sqrt(bin.mean());
                                double chiRxPPos = CalculateChi(RxPPosRes);
                                double res = Resolution(sqrt(2)*chiRxPPos);
                                EPres4[centBin] = res;
                        }
                        else EPres4[centBin] = 0;

                    }
                    centBin++;
            }

            for(unsigned int icent = 0; icent < v2centBins.size()-1; icent++)
            {
                    string v2string = "v2_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);
                    if(EPres[icent]>0) _p[v2string]->scaleY(1./EPres[icent]);
                    else throw UserError("EPres[icent] is less than/equal to 0, scaling not occuring");

                    string v3string = "v3_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);
                    if(EPres3[icent]>0) _p[v3string]->scaleY(1./EPres3[icent]);
                    else throw UserError("EPres3[icent] is less than/equal to 0, scaling not occuring");

                    string v4string = "v4_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);
                    if(EPres4[icent]>0) _p[v4string]->scaleY(1./EPres4[icent]);
                    else throw UserError("EPres4[icent] is less than/equal to 0, scaling not occuring");

                    string v4ep2string = "v4ep2_cent" + Form(v2centBins[icent], 0) + Form(v2centBins[icent+1], 0);
                    if(EPres[icent]>0) _p[v4ep2string]->scaleY(1./EPres[icent]);
                    else throw UserError("EPres[icent] is less than/equal to 0, scaling not occuring");

            }

/*
<<<<<<< HEAD
=======


	         for(Correlator& corr : Correlators)
      {
              corr.Normalize();
      }
>>>>>>> 968959d13c8f8d35f16d92923403a26295e69023*/
    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, Scatter2DPtr> _s;
    map<string, CounterPtr> _c;
    std::vector<double> v2centBins = {0., 10., 20., 30., 40., 50.};
    //@}

    vector<Correlator> Correlators;

    enum CollisionSystem {pp, AuAu};
    CollisionSystem collSys;

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2019_I1658594);

}
