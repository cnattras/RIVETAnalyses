// -*- C++ -*-
#ifndef RIVET_BackgroundRho_HH
#define RIVET_BackgroundRho_HH

#include "Rivet/Projection.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParticleVn.hh"
#include "Rivet/Projections/EventPlane.hh"


namespace Rivet
{

    /// @brief Project Rho and Rho(phi)
    class BackgroundRho : public Projection
    {
        public:
            
            BackgroundRho(const FastJets& fastjets, const int nLeadJetExclud, const double jetAreaCut, const Cut jetCuts);
            BackgroundRho(const int nLeadJetExclud, const double jetAreaCut);
        
            /// Clone on the heap.
            DEFAULT_RIVET_PROJ_CLONE(BackgroundRho);
        
            double getRho();
            double getRho(Jets jets);
            
            double getLocalRho(double phi, ParticleVn pvn, EventPlane ep, double jetR);
            double getNormalization(double nth, double jetR);
        
        
        protected:
            
            /// Apply the projection to the event.
            void project(const Event& e);
            
            CmpState compare(const Projection& p) const;
            
            int _nLeadJetExclud;
            double _jetAreaCut;
            Cut _jetCuts;
            double _rho;
            double _localRho;
            string _jetProjName;
        
    };
}

#endif
