// -*- C++ -*-
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Projections/FastJets.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include "Rivet/Projections/BackgroundRho.hh"

namespace Rivet { 

    //Constructors
    BackgroundRho::BackgroundRho(const FastJets& fastjets, const int nLeadJetExclud, const double jetAreaCut, const Cut jetCuts) :
    _nLeadJetExclud(nLeadJetExclud), _jetAreaCut(jetAreaCut), _jetCuts(jetCuts)
    {
        setName("BackgroundRho");
        _jetProjName = "fastjets";
        declare(fastjets,_jetProjName);
    }
    
    BackgroundRho::BackgroundRho(const int nLeadJetExclud, const double jetAreaCut) :
    _nLeadJetExclud(nLeadJetExclud), _jetAreaCut(jetAreaCut), _jetCuts(Cuts::open())
    {
        setName("BackgroundRho");
        _jetProjName = "";
    }
    
    void BackgroundRho::project(const Event& e)
    {
        Jets jets;
        
        if(_jetProjName != "") jets = apply<FastJets>(e, _jetProjName).jetsByPt(_jetCuts);
        
        _rho = getRho(jets);
    }
    
    CmpState BackgroundRho::compare(const Projection& p) const
    {
        const BackgroundRho& other = dynamic_cast<const BackgroundRho&>(p);
        if(_nLeadJetExclud != other._nLeadJetExclud) return CmpState::NEQ;
        if(_jetAreaCut != other._jetAreaCut) return CmpState::NEQ;
        
        return CmpState::EQ;
    }
    
    double BackgroundRho::getRho()
    {
        return _rho;
    }
    
    double BackgroundRho::getRho(Jets jetsSet)
    {
        Jets jets = sortBy(jetsSet,cmpMomByPt);
        vector<double> jetPtDensityVector;
        
        for(auto jet : jets)
        {
            if(jet.pseudojet().area() < _jetAreaCut) continue;
                                    
            double jetPtDensity = (jet.pT()/GeV)/(jet.pseudojet().area());
            jetPtDensityVector.push_back(jetPtDensity);
            
        }
        
        double nMediam = jetPtDensityVector.size() - _nLeadJetExclud; //
        double rho = 0.;
        int index = ceil(nMediam/2.)-1;
        
        if(index < 0)
        {
            //MSG_INFO("Only leading jets in the event! Cannot calculate rho.");
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
    
    double BackgroundRho::getLocalRho(double phi, ParticleVn pvn, EventPlane ep, double jetR)
    {
        std::vector<double> epAngle = ep.getAngleVector();
        std::vector<double> vn = pvn.getVnVector();
        
        std::vector<int> nthOrder = ep.getOrderVector();
        
        double deltaPhi;
        
        _localRho = 1.;
        
        for(unsigned int n = 0; n < nthOrder.size(); n++)
        {
            deltaPhi = mapAngle0To2Pi(phi-epAngle[n]);
            _localRho += 2.*cos(nthOrder[n]*deltaPhi)*getNormalization(nthOrder[n], jetR);
        }
        
        _localRho *= _rho;
        
        return _localRho;
    }
    
    double BackgroundRho::getNormalization(double nth, double jetR)
    {
        double norm = sin(nth*jetR)/(nth*jetR);
        return norm;
    }
    
}
