#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <iostream>
#include <string>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenModels/EvtXiccXicNpi.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include "EvtGenModels/EvtWHad.hh"

EvtXiccXicNpi::~EvtXiccXicNpi() {
    f1a = 1;
    f1b = 1;
    f1c = 1;
    f1d = 1;
    f2a = 1;
    f2b = 1;
    f2c = 1;
    f2d = 1;
    g1a = 1;
    g1b = 1;
    g1c = 1;
    g1d = 1;
    g2a = 1;
    g2b = 1;
    g2c = 1;
    g2d = 1;
}

std::string EvtXiccXicNpi::getName() {

    return "XICCXICNPI";

}

EvtDecayBase* EvtXiccXicNpi::clone() {

    return new EvtXiccXicNpi;

}

void EvtXiccXicNpi::init() {

    // check that there are 0 arguments
    checkNArg(0);

    checkSpinParent(EvtSpinType::DIRAC);

    checkSpinDaughter(0, EvtSpinType::DIRAC);

    EvtId parnum = getParentId();
    EvtId daughtnum = getDaug(0);
    
    // get FF parameters
    if (parnum == EvtPDL::getId("Xi_cc++") && daughtnum == EvtPDL::getId("Xi_c+")) {
        // \Xi_{cc}^{++}\Xi_{c}^{+}
        f1a = 0.914249;
        f1b = 0.347935;
        f1c = 0.0817908;
        f1d = 0.0186894;
        f2a = 0.0116363;
        f2b = 1.307;
        f2c = 0.512884;
        f2d = 0.149521;
        g1a = 0.258429;
        g1b = 0.208055;
        g1c = 0.0262441;
        g1d = -0.000279583;
        g2a = -0.0608048;
        g2b = 0.364131;
        g2c = 0.288789;
        g2d = -0.160034;
    } else {
        std::cout << "Wrong baryons\n";
        ::abort();
    };
}

void EvtXiccXicNpi::initProbMax() {
    EvtId parnum = getParentId();
    EvtId daughtnum1 = getDaug(0);
    EvtId daughtnum2 = getDaug(1);
    EvtId daughtnum3 = getDaug(2);
    size_t n = getNDaug();
    

    
    if(n==4) {// 3pi case
        if (parnum == EvtPDL::getId("Xi_cc++") && daughtnum1 == EvtPDL::getId("Xi_c+")) {
            setProbMax(6500);
            return;
        }
    }
    else     if(n==6) {// 5pi case
        if (parnum == EvtPDL::getId("Xi_cc++") && daughtnum1 == EvtPDL::getId("Xi_c+")) {
            setProbMax(7000);
            return;
        }
    };

    std::cout<<" Wrong mode"<<std::endl;
    return;
}
    

void EvtXiccXicNpi::HadronicAmp(EvtParticle* parent,
        EvtParticle* child,
        EvtVector4C* T,
        const int i,
        const int j) {
    const EvtDiracSpinor Sfinal = child->spParent(i);
    const EvtDiracSpinor Sinit = parent->sp(j);
    const EvtVector4R pp = parent->getP4Lab();
    const EvtVector4R cp = child->getP4Lab();
    const double pm = parent->mass();
    const double cm = child->mass();
    // \bar{u} \gamma^{\mu} u
    T[0] = EvtLeptonVCurrent(Sfinal, Sinit);

    // \bar{u} \gamma^{\mu}\gamma^{5} u
    T[1] = EvtLeptonACurrent(Sfinal, Sinit);

    // \bar{u} u
    T[2] = EvtLeptonSCurrent(Sfinal, Sinit) * ((pp + cp) / pm);

    // \bar{u} \gamma^{5} u
    T[3] = EvtLeptonPCurrent(Sfinal, Sinit) * ((pp + cp) / pm);

    return;
}

const double EvtXiccXicNpi::ff(const double f0, const double alpha, const double beta, const double gamma, EvtVector4R qqq) {
    return f0 * (1 + alpha * qqq * qqq + beta * qqq * qqq * qqq * qqq + gamma * qqq * qqq * qqq * qqq * qqq * qqq);
}

void EvtXiccXicNpi::decay(EvtParticle *b1) {
    static EvtId TAUM = EvtPDL::getId("tau-");

    b1->initializePhaseSpace(getNDaug(), getDaugs());

    EvtParticle *b2;
    b2 = b1->getDaug(0);



    const EvtVector4R pb1 = b1->getP4Lab();
    const EvtVector4R pb2 = b2->getP4Lab();
    EvtVector4R q = pb1 - pb2;
    const double m1 = b1->mass();
    const double m2 = b2->mass();


    const double f1 = ff(f1a, f1b, f1c, f1d, q);
    const double f2 = ff(f2a, f2b, f2c, f2d, q);
    const double g1 = ff(g1a, g1b, g1c, g1d, q);
    const double g2 = ff(g2a, g2b, g2c, g2d, q);

    EvtVector4C H[2][2]; // vector current
    EvtVector4C T[6];
    // Hadronic current
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            HadronicAmp(b1, b2, T, i, j);

            H[i][j] = (f1 * T[0] - g1 * T[1] + f2 * T[2] - f2 * (m1 + m2) / m1 * T[0] - g2 * T[3] - g2 * (m1 - m2) / m1 * T[1]);

        }
    }
    EvtId check1 = getDaug(1);


    //l nu_l case
    if (check1 == EvtPDL::getId("e+") || check1 == EvtPDL::getId("e-")) {
        EvtParticle *l, *nul;
        l = b1->getDaug(1);
        nul = b1->getDaug(2);
        //Leptonic current
        EvtVector4C lep[2];
        lep[0] = EvtLeptonVACurrent(l->spParent(0), nul->spParentNeutrino());
        lep[1] = EvtLeptonVACurrent(l->spParent(1), nul->spParentNeutrino());

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    vertex(i, j, k, H[i][j] * lep[k]);
                }
            }
        }
    }
    // 2pi case
    else if(getNDaug()==3) {
        EvtParticle *pi1 = b1->getDaug(1),
                *pi2 = b1->getDaug(2);
        EvtVector4R ppi1 = pi1->getP4Lab(),
                ppi2 = pi2->getP4Lab();
        EvtWHad qqqqqq = EvtWHad();
        EvtVector4C lep = qqqqqq.WCurrent(ppi1, ppi2);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                vertex(i, j, H[i][j] * lep);
            }
        }
    }
    //3pi case
    else if (getNDaug() == 4) {
        EvtParticle *pi1, *pi2, *pi3;

        pi1 = b1->getDaug(1);
        pi2 = b1->getDaug(2);
        pi3 = b1->getDaug(3);

        EvtVector4R ppi1, ppi2, ppi3; //Импульсы пионов
        ppi1 = pi1 -> getP4Lab();
        ppi2 = pi2 -> getP4Lab();
        ppi3 = pi3 -> getP4Lab();
        EvtWHad qqqqqq = EvtWHad();
        EvtVector4C lep = qqqqqq.WCurrent(ppi1, ppi2, ppi3);

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                vertex(i, j, H[i][j] * lep);
            }
        }
    }
    //5pi case
    else if (getNDaug() == 6) {
        EvtParticle *pi1, *pi2, *pi3, *pi4, *pi5;
        pi1 = b1->getDaug(1);
        pi2 = b1->getDaug(2);
        pi3 = b1->getDaug(3);
        pi4 = b1->getDaug(4);
        pi5 = b1->getDaug(5);
        EvtVector4R ppi1, ppi2, ppi3, ppi4, ppi5; //Импульсы пионов
        ppi1 = pi1 -> getP4Lab();
        ppi2 = pi2 -> getP4Lab();
        ppi3 = pi3 -> getP4Lab();
        ppi4 = pi4 -> getP4Lab();
        ppi5 = pi5 -> getP4Lab();
        EvtWHad qqqqqqq = EvtWHad();
        EvtVector4C lep = qqqqqqq.WCurrent(ppi1, ppi2, ppi3, ppi4, ppi5);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                vertex(i, j, H[i][j] * lep);
            }
        }
    }
    else {
        std::cout<<"Wrong final state!"<<std::endl;
        ::abort();
    }

    return;
}