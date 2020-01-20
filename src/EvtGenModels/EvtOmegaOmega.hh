#ifndef EVTOMEGAOMEGA_HH
#define EVTOMEGAOMEGA_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector3R.hh"


class EvtParticle;

class EvtOmegaOmega:public  EvtDecayAmp  {

public:

    EvtOmegaOmega() {}
    virtual ~EvtOmegaOmega();

    std::string getName();
    EvtDecayBase* clone();

    void initProbMax();
    void init();
    const double ff(const double f0, const double alpha, const double beta, const double gamma, EvtVector4R qqq);
    void HadronicAmp( EvtParticle* parent, 
                                 EvtParticle* child, 
                                 EvtVector4C* T,
                                 const int i, 
                                 const int j );
    EvtVector4C EvtSigmaCurrent(const EvtDiracSpinor &d, const EvtDiracSpinor &dp) ;

    void decay(EvtParticle *p); 
    double f1a, f1b, f1c, f1d, f2a, f2b, f2c, f2d, g1a, g1b, g1c, g1d, g2a, g2b, g2c, g2d;
};

#endif
