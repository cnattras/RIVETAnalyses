// -*- C++ -*-
#ifndef RIVET_RHICCENTRALITY_HH
#define RIVET_RHICCENTRALITY_HH

#include "Rivet/Projections/PercentileProjection.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include "Rivet/Tools/Cuts.hh"
#include <map>

namespace Rivet {

class RHICCentrality: public SingleValueProjection {

public:

  /// Default constructor.
  RHICCentrality(const std::string& name) {
    Cut c = Cuts::open();
    if(name == "STAR")
    {
        c = Cuts::abseta < 0.5;
        expName = "RHIC";
    }
    else if(name == "PHENIX")
    {
        c = Cuts::abseta > 3.1 && Cuts::abseta < 3.9;
        expName = "RHIC";
    }
    else if(name == "sPHENIX")
    {
        c = Cuts::abseta > 3.51 && Cuts::abseta < 4.61;
        expName = "RHIC";
    }
    else if(name == "CMS")
    {
        c = Cuts::abseta > 3 && Cuts::abseta < 5;
        expName = "CMS";
    }
    
    setName(name);
    declare(ChargedFinalState(c), "CFS");

}

  //Rivet 3.1
  //DEFAULT_RIVET_PROJ_CLONE(RHICCentrality);
  //Rivet 4.1 
  std::unique_ptr<Projection> clone() const override {
  return std::make_unique<RHICCentrality>(*this);
  }

  /// @BRIEF Add a new centality estimate.
  ///
  /// The SingelValueProjection, @a p, should return a value between 0
  /// and 100, and the @a pname should be one of "REF", "GEN", "IMP",
  /// "USR", or "RAW", as described above.
  void add(const SingleValueProjection & p, string pname) {
    _projNames.push_back(pname);
    declare(p, pname);
  }

  /// Perform all internal projections.
  /// Perform the projection
    void project(const Event& e) {
      clear();
      double estimate = 0;
      
      
      if(expName == "RHIC")
      {
          estimate = apply<FinalState>(e, "CFS").particles().size();
          
      }
      else if(expName == "CMS")
      {
          //Rivet 3.1
          //const Particles& particles = applyProjection<ChargedFinalState>(e, "CFS").particles();
          //Rivet 4.1
          const Particles& particles = apply<ChargedFinalState>(e, "CFS").particles();
          for(const Particle& p : particles)
          {
              estimate += p.Et();
          }
              
          
      }
      
      //cout << "Estimate: " << estimate << endl;

      setValue(estimate);
    }

  /// Cheek if no internal projections have been added.
  bool empty() const {
    return _projNames.empty();
  }

  /// Return the percentile of the @a i'th projection.
  ///
  /// Note that operator() will return the zero'th projection.
  double operator[](int i) const {
    return _values[i];
  }

  // Standard comparison function.
  /// Compare projections.
    CmpState compare(const Projection& p) const {
      return mkNamedPCmp(p, "CFS");
    }

  /// The list of names of the internal projections.
  vector<string> projections() const {
    return _projNames;
  }

private:

  /// The list of names of the internal projections.
  vector<string> _projNames;

  /// The list of percentiles resulting from the last projection.
  vector<double> _values;
  
  string expName = "";

};

}

#endif
