
/// \file
/// \ingroup tutorial_fit
/// \notebook
/// Combined (simultaneous) fit of two histogram with separate functions
/// and some common parameters
///
/// See http://root.cern.ch/phpBB3//viewtopic.php?f=3&t=11740#p50908
/// for a modified version working with Fumili or GSLMultiFit
///
/// N.B. this macro must be compiled with ACliC
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Lorenzo Moneta

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TFile.h"
#include <fstream>


const Int_t knrri=9;
const Int_t kmaxxpar=5;
const Int_t kmaxnddecay=10;
const Int_t kmaxxpaths=100;

Double_t neuueff=0.605;//changed to 62 %
Double_t neuueff_mean=0.605;
Double_t neuueff_err=0.053;


Bool_t reject=false;
Double_t rejecttrange=0.075;//first 0.075 ms


//! variable
Int_t npathhs=100;
Int_t ndeccay[kmaxxpaths];
Int_t decaymmap[kmaxxpaths][kmaxnddecay];
Int_t nneeu[kmaxxpaths][kmaxnddecay];


//! Global Bateaman function
Double_t corefcn(Int_t ndeccay,Int_t* decaymmap,Int_t* nneeu, Double_t* b1n,Double_t* b2n,Double_t* lamda,Double_t N0,Double_t t){
    Double_t fcnret=0;

    Double_t factor1=1.;
    //only parrent decay p2n
    for (int i=0;i<ndeccay-1;i++){
        if (nneeu[i]==0){
            factor1=factor1 * (1-b1n[decaymmap[i]-1]-b2n[decaymmap[i]-1])*lamda[decaymmap[i]-1];
        }else if (nneeu[i]==1){
            factor1=factor1 * b1n[decaymmap[i]-1]*lamda[decaymmap[i]-1];
        }else{
            factor1=factor1 * b2n[decaymmap[i]-1]*lamda[decaymmap[i]-1];
        }
    }

    Double_t factor2=0;
    for (int i=0;i<ndeccay;i++){
        Double_t factor2i=exp(-lamda[decaymmap[i]-1]*t);
        Double_t factor2ij=1;
        for (int j=0;j<ndeccay;j++)
            if (j!=i) factor2ij=factor2ij*(lamda[decaymmap[j]-1]-lamda[decaymmap[i]-1]);
        factor2=factor2+factor2i/factor2ij;
    }

    fcnret=factor1*N0*factor2;
    return fcnret;
}

//! Global function
Double_t fcn_gen(Double_t *x, Double_t *par) {
    Double_t bkg=par[knrri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;

    Double_t* pn=&par[knrri];
    Double_t* lamda=par;
    Double_t p2n[knrri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knrri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knrri*2+1]/par[0];

    //! Parent nuclei
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0]);

    for (Int_t i=0;i<npathhs;i++){
        returnval+=lamda[decaymmap[i][ndeccay[i]-1]-1]*corefcn(ndeccay[i],decaymmap[i],nneeu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejecttrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}


//! Global function
Double_t fcn_gen_bkg(Double_t *x, Double_t *par) {
    Double_t bkg=par[knrri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;
    return returnval;
}

Double_t fcn_gen_des(Double_t *x, Double_t *par) {
    Double_t bkg=par[knrri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;

    Double_t* pn=&par[knrri];
    Double_t* lamda=par;
    Double_t p2n[knrri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knrri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knrri*2+1]/par[0];

    //! Parent nuclei
    //returnval+=lamda[0]*N0*exp(-lamda[0]*x[0]);

    for (Int_t i=0;i<npathhs;i++){
        returnval+=lamda[decaymmap[i][ndeccay[i]-1]-1]*corefcn(ndeccay[i],decaymmap[i],nneeu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejecttrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}


Double_t fcn_gennuc1(Double_t *x, Double_t *par) {
    Double_t bkg=par[knrri*2+2];
    if (x[0]<0) return bkg;
    Double_t returnval=bkg;

    Double_t* pn=&par[knrri];
    Double_t* lamda=par;
    Double_t p2n[knrri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knrri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knrri*2+1]/par[0];

    //! Parent nuclei
    returnval+=lamda[0]*N0*exp(-lamda[0]*x[0]);

    //for (Int_t i=0;i<npathhs;i++){
        //returnval+=lamda[decaymmap[i][ndeccay[i]-1]-1]*corefcn(ndeccay[i],decaymmap[i],nneeu[i],pn,p2n,lamda,N0,x[0]);
    //}

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejecttrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gennuc2(Double_t *x, Double_t *par) {
    Double_t bkg=par[knrri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knrri];
    Double_t* lamda=par;
    Double_t p2n[knrri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knrri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knrri*2+1]/par[0];


    for (Int_t i=0;i<npathhs;i++){
        if (i==0)
        returnval+=lamda[decaymmap[i][ndeccay[i]-1]-1]*corefcn(ndeccay[i],decaymmap[i],nneeu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejecttrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gennuc3(Double_t *x, Double_t *par) {
    Double_t bkg=par[knrri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knrri];
    Double_t* lamda=par;
    Double_t p2n[knrri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knrri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knrri*2+1]/par[0];


    for (Int_t i=0;i<npathhs;i++){
        if (i==3)
        returnval+=lamda[decaymmap[i][ndeccay[i]-1]-1]*corefcn(ndeccay[i],decaymmap[i],nneeu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejecttrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gennuc6(Double_t *x, Double_t *par) {
    Double_t bkg=par[knrri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knrri];
    Double_t* lamda=par;
    Double_t p2n[knrri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knrri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knrri*2+1]/par[0];


    for (Int_t i=0;i<npathhs;i++){
        if (i==4)
        returnval+=lamda[decaymmap[i][ndeccay[i]-1]-1]*corefcn(ndeccay[i],decaymmap[i],nneeu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejecttrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}


Double_t fcn_gennuc4(Double_t *x, Double_t *par) {
    Double_t bkg=par[knrri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knrri];
    Double_t* lamda=par;
    Double_t p2n[knrri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knrri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knrri*2+1]/par[0];


    for (Int_t i=0;i<npathhs;i++){
        if (i==1)
        returnval+=lamda[decaymmap[i][ndeccay[i]-1]-1]*corefcn(ndeccay[i],decaymmap[i],nneeu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejecttrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gennuc5(Double_t *x, Double_t *par) {
    Double_t bkg=par[knrri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knrri];
    Double_t* lamda=par;
    Double_t p2n[knrri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knrri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knrri*2+1]/par[0];


    for (Int_t i=0;i<npathhs;i++){
        if (i==5||i==6)
        returnval+=lamda[decaymmap[i][ndeccay[i]-1]-1]*corefcn(ndeccay[i],decaymmap[i],nneeu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejecttrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

Double_t fcn_gennuc9(Double_t *x, Double_t *par) {
    Double_t bkg=par[knrri*2+2];
    if (x[0]<0) return 0;
    Double_t returnval=0;

    Double_t* pn=&par[knrri];
    Double_t* lamda=par;
    Double_t p2n[knrri]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    p2n[0]=par[knrri*2];//special for p2n of parrent, to be develope later
    Double_t N0=par[knrri*2+1]/par[0];

    for (Int_t i=0;i<npathhs;i++){
        if (i==7||i==8)
        returnval+=lamda[decaymmap[i][ndeccay[i]-1]-1]*corefcn(ndeccay[i],decaymmap[i],nneeu[i],pn,p2n,lamda,N0,x[0]);
    }

    //! reject 50 ms bump
    if (reject && x[0] > 0 && x[0] < rejecttrange) {
       TF1::RejectPoint();
       return 0;
    }
    return returnval;
}

TF1* getFunction(Double_t*parms, Double_t irejecttrange=0.05, Double_t upperlimit=20, Double_t neutronefficiency=0.62){
    rejecttrange=irejecttrange;
    neuueff=neutronefficiency;
    neuueff_mean=neutronefficiency;
    //! define decay map
    npathhs=12;
    //! path 1(go for ri2)
    ndeccay[0]=2;
    decaymmap[0][0]=1;decaymmap[0][1]=2;
    nneeu[0][0]=0;
    //! path 2(go for ri4)
    ndeccay[1]=3;
    decaymmap[1][0]=1;decaymmap[1][1]=2;decaymmap[1][2]=4;
    nneeu[1][0]=0;nneeu[1][1]=0;
    //! path 3(go for ri7)
    ndeccay[2]=4;
    decaymmap[2][0]=1;decaymmap[2][1]=2;decaymmap[2][2]=4;decaymmap[2][3]=7;
    nneeu[2][0]=0;nneeu[2][1]=0;nneeu[2][2]=0;
    //! path4(go for ri3)
    ndeccay[3]=2;
    decaymmap[3][0]=1;decaymmap[3][1]=3;
    nneeu[3][0]=1;
    //!path5(go for ri6)
    ndeccay[4]=2;
    decaymmap[4][0]=1;decaymmap[4][1]=6;
    nneeu[4][0]=2;

    //! path6(go for ri5-route 1)
    ndeccay[5]=3;
    decaymmap[5][0]=1;decaymmap[5][1]=2;decaymmap[5][2]=5;
    nneeu[5][0]=0;nneeu[5][1]=1;
    //! path7(go for ri5-route 2)
    ndeccay[6]=3;
    decaymmap[6][0]=1;decaymmap[6][1]=3;decaymmap[6][2]=5;
    nneeu[6][0]=1;nneeu[6][1]=0;

    //! path8 (go for ri9-route 1)
    ndeccay[7]=3;
    decaymmap[7][0]=1;decaymmap[7][1]=3;decaymmap[7][2]=9;
    nneeu[7][0]=1;nneeu[7][1]=1;
    //! path9 (go for ri9-route 2)
    ndeccay[8]=3;
    decaymmap[8][0]=1;decaymmap[8][1]=6;decaymmap[8][2]=9;
    nneeu[8][0]=2;nneeu[8][1]=0;

    //! path10 (go for ri8-route 1)
    ndeccay[9]=4;
    decaymmap[9][0]=1;decaymmap[9][1]=2;decaymmap[9][2]=4;decaymmap[9][3]=8;
    nneeu[9][0]=0;nneeu[9][1]=0;nneeu[9][2]=1;

    //! path11 (go for ri8-route 2)
    ndeccay[10]=4;
    decaymmap[10][0]=1;decaymmap[10][1]=2;decaymmap[10][2]=5;decaymmap[10][3]=8;
    nneeu[10][0]=0;nneeu[10][1]=1;nneeu[10][2]=0;

    //! path12 (go for ri8-route 3)
    ndeccay[11]=4;
    decaymmap[11][0]=1;decaymmap[11][1]=3;decaymmap[11][2]=5;decaymmap[11][3]=8;
    nneeu[11][0]=1;nneeu[11][1]=0;nneeu[11][2]=0;

    //!******************************************Define All BETA decay function
    TF1* fB=new TF1("fB",fcn_gen,rejecttrange,upperlimit,knrri*2+3);
    fB->SetNpx(2000);
    fB->SetLineWidth(2);
    fB->SetLineWidth(0);

    //!****************************************** All BETA decay function
    //! read input and get parameters
    for (int i=0;i<knrri*2+1;i++){
        fB->FixParameter(i,parms[i]);
    }
    fB->SetParameter(knrri*2+1,parms[knrri*2+1]);
    //fB->SetParameter(knrri*2+2,parms[knrri*2+2]);
    fB->FixParameter(knrri*2+2,parms[knrri*2+2]);

    fB->SetParLimits(knrri*2+1,parms[knrri*2+1]/5,parms[knrri*2+1]*5);
    //fB->SetParLimits(knrri*2+2,parms[knrri*2+2]/5,parms[knrri*2+2]*5);
    return fB;
}

TF1* getDaughterFunction(Double_t*parms, Double_t irejecttrange=0.05, Double_t upperlimit=20, Double_t neutronefficiency=0.62){
    rejecttrange=irejecttrange;
    neuueff=neutronefficiency;
    neuueff_mean=neutronefficiency;
    //! define decay map
    npathhs=12;
    //! path 1(go for ri2)
    ndeccay[0]=2;
    decaymmap[0][0]=1;decaymmap[0][1]=2;
    nneeu[0][0]=0;
    //! path 2(go for ri4)
    ndeccay[1]=3;
    decaymmap[1][0]=1;decaymmap[1][1]=2;decaymmap[1][2]=4;
    nneeu[1][0]=0;nneeu[1][1]=0;
    //! path 3(go for ri7)
    ndeccay[2]=4;
    decaymmap[2][0]=1;decaymmap[2][1]=2;decaymmap[2][2]=4;decaymmap[2][3]=7;
    nneeu[2][0]=0;nneeu[2][1]=0;nneeu[2][2]=0;
    //! path4(go for ri3)
    ndeccay[3]=2;
    decaymmap[3][0]=1;decaymmap[3][1]=3;
    nneeu[3][0]=1;
    //!path5(go for ri6)
    ndeccay[4]=2;
    decaymmap[4][0]=1;decaymmap[4][1]=6;
    nneeu[4][0]=2;

    //! path6(go for ri5-route 1)
    ndeccay[5]=3;
    decaymmap[5][0]=1;decaymmap[5][1]=2;decaymmap[5][2]=5;
    nneeu[5][0]=0;nneeu[5][1]=1;
    //! path7(go for ri5-route 2)
    ndeccay[6]=3;
    decaymmap[6][0]=1;decaymmap[6][1]=3;decaymmap[6][2]=5;
    nneeu[6][0]=1;nneeu[6][1]=0;

    //! path8 (go for ri9-route 1)
    ndeccay[7]=3;
    decaymmap[7][0]=1;decaymmap[7][1]=3;decaymmap[7][2]=9;
    nneeu[7][0]=1;nneeu[7][1]=1;
    //! path9 (go for ri9-route 2)
    ndeccay[8]=3;
    decaymmap[8][0]=1;decaymmap[8][1]=6;decaymmap[8][2]=9;
    nneeu[8][0]=2;nneeu[8][1]=0;

    //! path10 (go for ri8-route 1)
    ndeccay[9]=4;
    decaymmap[9][0]=1;decaymmap[9][1]=2;decaymmap[9][2]=4;decaymmap[9][3]=8;
    nneeu[9][0]=0;nneeu[9][1]=0;nneeu[9][2]=1;

    //! path11 (go for ri8-route 2)
    ndeccay[10]=4;
    decaymmap[10][0]=1;decaymmap[10][1]=2;decaymmap[10][2]=5;decaymmap[10][3]=8;
    nneeu[10][0]=0;nneeu[10][1]=1;nneeu[10][2]=0;

    //! path12 (go for ri8-route 3)
    ndeccay[11]=4;
    decaymmap[11][0]=1;decaymmap[11][1]=3;decaymmap[11][2]=5;decaymmap[11][3]=8;
    nneeu[11][0]=1;nneeu[11][1]=0;nneeu[11][2]=0;

    //!******************************************Define All BETA decay function
    TF1* fBdaugter=new TF1("fBdaugter",fcn_gen_des,rejecttrange,upperlimit,knrri*2+3);
    fBdaugter->SetNpx(2000);
    fBdaugter->SetLineWidth(2);
    fBdaugter->SetLineColor(8);

    //!****************************************** All BETA decay function
    //! read input and get parameters
    for (int i=0;i<knrri*2+1;i++){
        fBdaugter->FixParameter(i,parms[i]);
    }
    fBdaugter->SetParameter(knrri*2+1,parms[knrri*2+1]);
    fBdaugter->SetParameter(knrri*2+2,parms[knrri*2+2]);
    fBdaugter->SetParLimits(knrri*2+1,parms[knrri*2+1]/5,parms[knrri*2+1]*5);
    fBdaugter->SetParLimits(knrri*2+2,parms[knrri*2+2]/5,parms[knrri*2+2]*5);
    return fBdaugter;
}

TF1* getParentFunction(Double_t*parms, Double_t irejecttrange=0.05, Double_t upperlimit=20, Double_t neutronefficiency=0.62){
    rejecttrange=irejecttrange;
    neuueff=neutronefficiency;
    neuueff_mean=neutronefficiency;
    //! define decay map
    npathhs=12;
    //! path 1(go for ri2)
    ndeccay[0]=2;
    decaymmap[0][0]=1;decaymmap[0][1]=2;
    nneeu[0][0]=0;
    //! path 2(go for ri4)
    ndeccay[1]=3;
    decaymmap[1][0]=1;decaymmap[1][1]=2;decaymmap[1][2]=4;
    nneeu[1][0]=0;nneeu[1][1]=0;
    //! path 3(go for ri7)
    ndeccay[2]=4;
    decaymmap[2][0]=1;decaymmap[2][1]=2;decaymmap[2][2]=4;decaymmap[2][3]=7;
    nneeu[2][0]=0;nneeu[2][1]=0;nneeu[2][2]=0;
    //! path4(go for ri3)
    ndeccay[3]=2;
    decaymmap[3][0]=1;decaymmap[3][1]=3;
    nneeu[3][0]=1;
    //!path5(go for ri6)
    ndeccay[4]=2;
    decaymmap[4][0]=1;decaymmap[4][1]=6;
    nneeu[4][0]=2;

    //! path6(go for ri5-route 1)
    ndeccay[5]=3;
    decaymmap[5][0]=1;decaymmap[5][1]=2;decaymmap[5][2]=5;
    nneeu[5][0]=0;nneeu[5][1]=1;
    //! path7(go for ri5-route 2)
    ndeccay[6]=3;
    decaymmap[6][0]=1;decaymmap[6][1]=3;decaymmap[6][2]=5;
    nneeu[6][0]=1;nneeu[6][1]=0;

    //! path8 (go for ri9-route 1)
    ndeccay[7]=3;
    decaymmap[7][0]=1;decaymmap[7][1]=3;decaymmap[7][2]=9;
    nneeu[7][0]=1;nneeu[7][1]=1;
    //! path9 (go for ri9-route 2)
    ndeccay[8]=3;
    decaymmap[8][0]=1;decaymmap[8][1]=6;decaymmap[8][2]=9;
    nneeu[8][0]=2;nneeu[8][1]=0;

    //! path10 (go for ri8-route 1)
    ndeccay[9]=4;
    decaymmap[9][0]=1;decaymmap[9][1]=2;decaymmap[9][2]=4;decaymmap[9][3]=8;
    nneeu[9][0]=0;nneeu[9][1]=0;nneeu[9][2]=1;

    //! path11 (go for ri8-route 2)
    ndeccay[10]=4;
    decaymmap[10][0]=1;decaymmap[10][1]=2;decaymmap[10][2]=5;decaymmap[10][3]=8;
    nneeu[10][0]=0;nneeu[10][1]=1;nneeu[10][2]=0;

    //! path12 (go for ri8-route 3)
    ndeccay[11]=4;
    decaymmap[11][0]=1;decaymmap[11][1]=3;decaymmap[11][2]=5;decaymmap[11][3]=8;
    nneeu[11][0]=1;nneeu[11][1]=0;nneeu[11][2]=0;

    //!******************************************Define All BETA decay function
    TF1* fBparent=new TF1("fBparent",fcn_gennuc1,rejecttrange,upperlimit,knrri*2+3);
    fBparent->SetNpx(2000);
    fBparent->SetLineWidth(2);
    fBparent->SetLineColor(8);

    //!****************************************** All BETA decay function
    //! read input and get parameters
    for (int i=0;i<knrri*2+1;i++){
        fBparent->FixParameter(i,parms[i]);
    }
    fBparent->SetParameter(knrri*2+1,parms[knrri*2+1]);
    fBparent->SetParameter(knrri*2+2,parms[knrri*2+2]);
    fBparent->SetParLimits(knrri*2+1,parms[knrri*2+1]/5,parms[knrri*2+1]*5);
    fBparent->SetParLimits(knrri*2+2,parms[knrri*2+2]/5,parms[knrri*2+2]*5);
    return fBparent;
}


TF1* getBkgFunction(Double_t*parms, Double_t irejecttrange=0.05, Double_t upperlimit=20, Double_t neutronefficiency=0.62){
    rejecttrange=irejecttrange;
    neuueff=neutronefficiency;
    neuueff_mean=neutronefficiency;
    //! define decay map
    npathhs=12;
    //! path 1(go for ri2)
    ndeccay[0]=2;
    decaymmap[0][0]=1;decaymmap[0][1]=2;
    nneeu[0][0]=0;
    //! path 2(go for ri4)
    ndeccay[1]=3;
    decaymmap[1][0]=1;decaymmap[1][1]=2;decaymmap[1][2]=4;
    nneeu[1][0]=0;nneeu[1][1]=0;
    //! path 3(go for ri7)
    ndeccay[2]=4;
    decaymmap[2][0]=1;decaymmap[2][1]=2;decaymmap[2][2]=4;decaymmap[2][3]=7;
    nneeu[2][0]=0;nneeu[2][1]=0;nneeu[2][2]=0;
    //! path4(go for ri3)
    ndeccay[3]=2;
    decaymmap[3][0]=1;decaymmap[3][1]=3;
    nneeu[3][0]=1;
    //!path5(go for ri6)
    ndeccay[4]=2;
    decaymmap[4][0]=1;decaymmap[4][1]=6;
    nneeu[4][0]=2;

    //! path6(go for ri5-route 1)
    ndeccay[5]=3;
    decaymmap[5][0]=1;decaymmap[5][1]=2;decaymmap[5][2]=5;
    nneeu[5][0]=0;nneeu[5][1]=1;
    //! path7(go for ri5-route 2)
    ndeccay[6]=3;
    decaymmap[6][0]=1;decaymmap[6][1]=3;decaymmap[6][2]=5;
    nneeu[6][0]=1;nneeu[6][1]=0;

    //! path8 (go for ri9-route 1)
    ndeccay[7]=3;
    decaymmap[7][0]=1;decaymmap[7][1]=3;decaymmap[7][2]=9;
    nneeu[7][0]=1;nneeu[7][1]=1;
    //! path9 (go for ri9-route 2)
    ndeccay[8]=3;
    decaymmap[8][0]=1;decaymmap[8][1]=6;decaymmap[8][2]=9;
    nneeu[8][0]=2;nneeu[8][1]=0;

    //! path10 (go for ri8-route 1)
    ndeccay[9]=4;
    decaymmap[9][0]=1;decaymmap[9][1]=2;decaymmap[9][2]=4;decaymmap[9][3]=8;
    nneeu[9][0]=0;nneeu[9][1]=0;nneeu[9][2]=1;

    //! path11 (go for ri8-route 2)
    ndeccay[10]=4;
    decaymmap[10][0]=1;decaymmap[10][1]=2;decaymmap[10][2]=5;decaymmap[10][3]=8;
    nneeu[10][0]=0;nneeu[10][1]=1;nneeu[10][2]=0;

    //! path12 (go for ri8-route 3)
    ndeccay[11]=4;
    decaymmap[11][0]=1;decaymmap[11][1]=3;decaymmap[11][2]=5;decaymmap[11][3]=8;
    nneeu[11][0]=1;nneeu[11][1]=0;nneeu[11][2]=0;

    //!******************************************Define All BETA decay function
    TF1* fBbkg=new TF1("fBbkg",fcn_gen_bkg,rejecttrange,upperlimit,knrri*2+3);
    fBbkg->SetNpx(2000);
    fBbkg->SetLineWidth(2);
    fBbkg->SetLineColor(8);

    //!****************************************** All BETA decay function
    //! read input and get parameters
    for (int i=0;i<knrri*2+1;i++){
        fBbkg->FixParameter(i,parms[i]);
    }
    fBbkg->SetParameter(knrri*2+1,parms[knrri*2+1]);
    fBbkg->SetParameter(knrri*2+2,parms[knrri*2+2]);
    fBbkg->SetParLimits(knrri*2+1,parms[knrri*2+1]/5,parms[knrri*2+1]*5);
    fBbkg->SetParLimits(knrri*2+2,parms[knrri*2+2]/5,parms[knrri*2+2]*5);
    return fBbkg;
}



