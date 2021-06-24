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




class unbinfit
{
  public:
    unbinfit();
    virtual ~unbinfit();
    void Init(char* inputParms, char* inputData);
    void setInputEffParms(char* inputEffParms){sprintf(finputEffParms,"%s",inputEffParms);}

    void setOutputFile(char* outputData){sprintf(foutputData,"%s",outputData);}

    void setStartTime(double deadtime){p_deadtime=deadtime;}

    void fitBackground(Int_t opt=0);
 private:
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
    fitFbkgflat* bkgmodelneg;
    fitFbkgflat* bkgmodelpos;

    //! data sets
    TTree* tree;
    TTree* treeb;
    RooDataSet* databkg; //background data
    RooDataSet* data;
};

#endif

