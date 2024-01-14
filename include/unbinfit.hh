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

//! \file unbinfit.hh
//! \brief Definition of the unbinfit class

#ifndef unbinfit_h
#define unbinfit_h 1

//! Unbinned fitting class.

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"

#include "RooDataSet.h"

#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooPlot.h"

#include "RooAddPdf.h"
#include "TString.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFitResult.h"

#include "RooCurve.h"
#include "RooHist.h"

#include "TStyle.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "TStopwatch.h"
#include "TRandom3.h"

#include "common.hh"
#include "fitF.hh"

#include "decaypath.hh"


#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"

#include "asymGausRandom.hh"

class asymGausRandom;

class unbinfit
{
  public:
    unbinfit();
    virtual ~unbinfit();
    void Init(char* inputParms, char* inputData);
    void setInputEffParms(char* inputEffParms){sprintf(finputEffParms,"%s",inputEffParms);}

    void SetParameters();

    void setMinEffMC(Double_t effin){fmineffMC = effin;}
    void setMaxEffMC(Double_t effin){fmaxeffMC = effin;}

    void setOutputFile(char* outputData){sprintf(foutputData,"%s",outputData);}
    void setEntriesLimit(Long64_t entries){fnentrieslimit = entries;}
    void setStartTime(double deadtime){p_deadtime=deadtime;}
    void setTimeRange(double timerange){p_timerange = timerange;}
    void Run();
    void RunBinFit();

    void generateRoofitEvaluate();

    void fitBackground(Int_t opt=0);

    void initFitParameters();
    void setNormalFit();//decide parameter is fix or not
    void setExernalContrainFit();

    void prepareData();
    void SetRandomSeed(Int_t seednoin){seedno = seednoin;}
    void prepareMonteCarloData(int nevents);
    void doFit();

    void calculateUpperLimit();//calculate upper limit
    void plotResults();
    void plotResultsMore(Int_t opt=0);
    void writeResults();
    void writeOutputTree(){foutputtree->Write();}
    void closeOutputFile(){fout->Close();}

    //! stuff for MC set
    void setNumberOfMC(Int_t ntimes){fnMC=ntimes;}

    void setCentralParameters();
    void bookOutputTree();
    void getParameters();
    void printCurrentParameters();
    void setValParameters();
    void generateMC();
    void writeResultsMC();

    void setNBinHists(Int_t nbins){
        nbinsHB=nbins;
        nbinsHSB=nbins;
        nbinsHSB2=nbins;
    }

    void setdataHistzero(RooHist* datahist){
        for (Int_t i=0;i<datahist->GetN();i++){if (datahist->GetY()[i]==0) {
                datahist->SetPointEYhigh(i,0);
                datahist->SetPointEYlow(i,0);
            }
        }
    }

 private:
    void setModel();
    void calculateChiSquare(Int_t opt=0);
    void writeFitComponents();
    char* finputParms;
    char* finputData;
    char finputEffParms[500];
    char foutputData[500];
    TFile* fout;
    Long64_t fnentrieslimit;
    Double_t p_deadtime;
    Double_t p_timerange;

    Int_t ncpu;

    Int_t ffitopt;

    decaypath* fdecaypath;

    //! 2 dimensions parameters
    RooRealVar *x;
    RooCategory *y;

    //! background parameters
    RooRealVar* xbkg;
    RooRealVar* bkg1nratio;
    RooRealVar* bkg2nratio;

    RooRealVar* slope1pos;
    RooRealVar* slope2pos;
    RooRealVar* slope3pos;

    fitFbkg* bkgmodelneg;
    fitFbkg* bkgmodelpos;

    //! fit paramters
    // Declare all parameters
    RooAbsReal* p[kmaxparms];//fdecaypath->getNMember()*5+4];
    RooRealVar* pvar[kmaxparms];//[fdecaypath->getNMember()*5+4];
    RooRealVar* nbkg;
    RooRealVar* nsig;

    //! constrains set
    RooArgSet *externalconstrains;


    //! constrains set
    fitF* totdecaymodel;
    RooAddPdf* final_pdf;


    //! data sets
    TTree* tree;
    TTree* treeb;
    RooDataSet* databkg; //background data
    RooDataSet* data;

    //! fit results
    RooFitResult* fitres;

    //! stuffs for MC generation
    TRandom3* rseed;
    asymGausRandom* rseedA;

    Int_t seedno;
    Int_t fnMC;
    TTree* foutputtree;

    Double_t pVal[kmaxparms];
    Double_t pCentralVal[kmaxparms];
    Double_t pValError[kmaxparms];
    Double_t pValErrorHi[kmaxparms];
    Int_t ispVary[kmaxparms];
    Int_t ipVal[kmaxparms];

    Double_t nbkgVal;
    Double_t nbkgCentralVal;
    Double_t nbkgError;
    Double_t bkg1nratioVal;
    Double_t bkg1nratioCentralVal;
    Double_t bkg1nratioError;
    Double_t bkg2nratioVal;
    Double_t bkg2nratioCentralVal;
    Double_t bkg2nratioError;
    Double_t slope1posVal;
    Double_t slope2posVal;
    Double_t slope3posVal;
    Double_t slope1posCentralVal;
    Double_t slope2posCentralVal;
    Double_t slope3posCentralVal;
    Double_t slope1posError;
    Double_t slope2posError;
    Double_t slope3posError;

    Int_t fitStatus;
    Int_t fitCovQual;
    Int_t fitNumInvalidNLL;
    Double_t fitEdm;
    Double_t fitMinNll;
    Double_t chiSquareNDF;
    Double_t chiSquareNDF1n;
    Double_t chiSquareNDF2n;

    Double_t nsigVal;
    Double_t nsigCentralVal;
    Double_t nsigError;

    TStopwatch *fStopWatch;//for debuging
    Double_t fMCGenTime;
    Double_t fFitTime;

    Int_t nbinsHB;
    Int_t nbinsHSB;
    Int_t nbinsHSB2;

    TH1F* hB;
    TH1F* hSB;
    TH1F* hSB2;

    TF1* fB;
    TF1* fSB;
    TF1* fSB2;
    TF1* fB_bkgneg;
    TF1* fSB_bkgneg;
    TF1* fSB2_bkgneg;
    TF1* fB_bkgpos;
    TF1* fSB_bkgpos;
    TF1* fSB2_bkgpos;


    fitF* totdecaymodelforplot;
    TF1* fB_parent;
    TF1* fB_daugter;
    TF1* fSB_parent;
    TF1* fSB_daugter;
    TF1* fSB2_parent;
    TF1* fSB2_daugter;

    TF1* fSB_c1;
    TF1* fSB_c2;
    TF1* fSB_c3;
    TF1* fSB_c23;
    TF1* fSB2_c1;
    TF1* fSB2_c2;
    TF1* fSB2_c3;
    TF1* fSB2_c4;
    TF1* fSB2_c134;

    RooCurve* model0nCurve;
    RooHist* model0nHist;
    RooCurve* model1nCurve;
    RooHist* model1nHist;
    RooCurve* model2nCurve;
    RooHist* model2nHist;

    RooCurve* modelbkg0nCurve;
    RooHist* modelbkg0nHist;
    RooCurve* modelbkg1nCurve;
    RooHist* modelbkg1nHist;
    RooCurve* modelbkg2nCurve;
    RooHist* modelbkg2nHist;

    Double_t plotrangelow;
    Double_t plotrangehi;

    Double_t fmineffMC;
    Double_t fmaxeffMC;

    Double_t nsig_hB_firstbin;
    Double_t binfitparms[kmaxparms];
    Double_t binfitbkgparms[6];

    //! stuffs for upper limits calculations
    Double_t N0b;
    Double_t N0b1n;
    Double_t N0b2n;
    Double_t N0b3n;

    Double_t N0b1n_bkg;
    Double_t N0b2n_bkg;
    Double_t N0b3n_bkg;


    Double_t binw;

    struct GlobalChi2 {
       GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                    ROOT::Math::IMultiGenFunction & f2,
                    ROOT::Math::IMultiGenFunction & f3) :
          fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3)
       {
           std::ifstream pathfile("path.txt");
           pathfile>>nri;
           pathfile.close();
       }

       double operator() (const double *par) const {
           double p1[nri*5+10];
           for (int i = 0; i < nri*5+8; ++i) p1[i] = par[i];
           p1[nri*5+8]=par[nri*5+8];
           p1[nri*5+9]=par[nri*5+9];
           double p2[nri*5+10];
           for (int i = 0; i < nri*5+8; ++i) p2[i] = par[i];
           p2[nri*5+8]=par[nri*5+10];
           p2[nri*5+9]=par[nri*5+11];
           double p3[nri*5+10];
           for (int i = 0; i < nri*5+8; ++i) p3[i] = par[i];
           p3[nri*5+8]=par[nri*5+12];
           p3[nri*5+9]=par[nri*5+13];
          return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3);
       }
       const  ROOT::Math::IMultiGenFunction * fChi2_1;
       const  ROOT::Math::IMultiGenFunction * fChi2_2;
       const  ROOT::Math::IMultiGenFunction * fChi2_3;

       Int_t nri;

    };

};

#endif

