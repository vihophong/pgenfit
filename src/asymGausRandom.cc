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
/// \file asymGausRandom.cc
/// \brief Implementation of the asymGausRandom class

#include "asymGausRandom.hh"

asymGausRandom::asymGausRandom()
{
//    unr = new TUnuran;
    xx=new RooRealVar("xx","xx",5,40);
    mean=new RooRealVar("mean","mean",10,5,20);
    sigmaL=new RooRealVar("sigmaL","sigmaL",1,0.5,2);
    sigmaR=new RooRealVar("sigmaR","sigmaR",5,1,4);
}

asymGausRandom::~asymGausRandom()
{

}

Double_t asymGausRandom::asymGausFcn(Double_t *x, Double_t *par)
{
    //! par[0] -mean; par[1]- sigmalow; par[2] - sigmahigh
    Double_t ret = 0;
    if (x[0]<par[0]){
        ret = TMath::Gaus(x[0],par[0],par[1]);
    }else{
        ret = TMath::Gaus(x[0],par[0],par[2]);
    }
    return ret;
}


void asymGausRandom::Init(Int_t seedno){
    RooRandom::randomGenerator()->SetSeed(seedno);
    bfGaus=new RooBifurGauss("bigaus","Bi Gaussian PDF",*xx,*mean,*sigmaL,*sigmaR);
//    unr->SetSeed(seedno);
    fAsymGaus = new TF1("fAsymGaus",this,&asymGausRandom::asymGausFcn,-10,10,3,"asymGausRandom","asymGausFcn");
//    dist = new TUnuranContDist(fAsymGaus);
}

Double_t asymGausRandom::generate(Double_t meanVal, Double_t sigmaVal){
//    fAsymGaus->SetParameter(0,0);
//    fAsymGaus->SetParameter(1,1);
//    fAsymGaus->SetParameter(2,1);
//    unr->Init(*dist,"method=arou");
//    return unr->Sample()*sigmaVal+meanVal;
    //std::cout<<"EEEE"<<meanVal<<"\t"<<sigmaVal<<std::endl;
    xx->setRange(meanVal-sigmaVal*10,meanVal+sigmaVal*10);
    mean->setRange(meanVal/10,meanVal*10);mean->setVal(meanVal);
    sigmaL->setRange(sigmaVal/10,sigmaVal*10);sigmaL->setVal(sigmaVal);
    sigmaR->setRange(sigmaVal/10,sigmaVal*10);sigmaR->setVal(sigmaVal);
    RooDataSet *data = bfGaus->generate(*xx,1);
    Double_t returnVal= data->get(0)->getRealValue("xx");
    delete data;
    return returnVal;
}

Double_t asymGausRandom::generate(Double_t meanVal, Double_t sigmaLVal,Double_t sigmaRVal){
//    fAsymGaus->SetParameter(0,0);
//    fAsymGaus->SetParameter(1,1);
//    fAsymGaus->SetParameter(2,sigmaRVal/sigmaLVal);
//    unr->Init(*dist,"method=arou");
//    return unr->Sample()*sigmaLVal+meanVal;
    xx->setRange(meanVal-sigmaLVal*10,meanVal+sigmaRVal*10);
    mean->setRange(meanVal/10,meanVal*10);mean->setVal(meanVal);
    sigmaL->setRange(sigmaLVal/10,sigmaLVal*10);sigmaL->setVal(sigmaLVal);
    sigmaR->setRange(sigmaRVal/10,sigmaRVal*10);sigmaR->setVal(sigmaRVal);
    RooDataSet *data = bfGaus->generate(*xx,1);
    Double_t returnVal= data->get(0)->getRealValue("xx");
    delete data;
    return returnVal;
}
