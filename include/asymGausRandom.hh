//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * Copyright@2019 Vi Ho Phong, email: phong@ribf.riken.jp           *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications.                    *
// ********************************************************************
//

//! \file asymGausRandom.hh
//! \brief A class to generate asymetric uniform random number

#ifndef asymGausRandom_h
#define asymGausRandom_h 1

#include "TTree.h"
#include "TMath.h"
//#include "TUnuran.h"
//#include "TUnuranContDist.h"
#include "TF1.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooBifurGauss.h"
#include "RooRandom.h"
using namespace RooFit;

class asymGausRandom
{
public:
    asymGausRandom();
    virtual ~asymGausRandom();
    Double_t generate(Double_t meanVal, Double_t sigmaVal);
    Double_t generate(Double_t meanVal, Double_t sigmaLVal,Double_t sigmaRVal);
//    void SetSeed(Int_t seedno){unr->SetSeed(seedno);}
    void Init(Int_t seedno);
private:
//    TUnuran* unr;
//    TUnuranContDist* dist;
    TF1* fAsymGaus;
    Double_t asymGausFcn(Double_t *x, Double_t *par);

    RooRealVar* xx;
    RooRealVar* mean;
    RooRealVar* sigmaL;
    RooRealVar* sigmaR;
    RooBifurGauss *bfGaus;
};


#endif
