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
/// \file unbinfit.cc
/// \brief Implementation of the unbinfit class

#include "unbinfit.hh"
#include "TSystem.h"
#include <iostream>



using namespace RooFit;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unbinfit::unbinfit()
{
    ffitopt=0;// only 0 - fix parameter; 1 - constrain parameters

    fnentrieslimit=ENTRYLIMIT;

    p_deadtime=STARTFIT;
    p_timerange=10;

    ncpu=NCPUS_UNBINFIT;

    finputData=new char[1000];
    finputParms=new char[1000];

    fdecaypath=new decaypath;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unbinfit::~unbinfit()
{
    delete finputData;
    delete finputParms;
    delete fdecaypath;
    delete tree;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::Init(char* inputParms, char* inputData)
{

    sprintf(finputParms,"%s",inputParms);
    sprintf(finputData,"%s",inputData);
    std::clog<< __PRETTY_FUNCTION__ <<"read input files:"<<
               std::endl<<finputParms<<
               std::endl<<finputData<<std::endl;

    fdecaypath->Init(finputParms);
    fdecaypath->makePath();
    fdecaypath->printMember();
    fdecaypath->printPath();
    fdecaypath->writePath();
    fdecaypath->drawPath((char*)TString("outdecayroutes.root").Data());

    //check if dummy input file
    if (gSystem->AccessPathName(finputData))
        return;
    // Import tree to total decay data
    TFile *f=TFile::Open(finputData);
    f->GetObject("tree",tree);
    f->GetObject("treebw",treeb);
    if (fnentrieslimit>0) tree->SetEntries(fnentrieslimit);

    //!*****************************************
    //! Prepare X,y
    //! *****************************************
    x=new RooRealVar("x","x",p_deadtime,p_timerange) ;
    xbkg=new RooRealVar("x","x",-p_timerange,-p_deadtime) ;
    // define discrete variable y
    y=new RooCategory("y","y");
    y->defineType("0neu",0);
    y->defineType("1neu",1);
    y->defineType("2neu",2);

    // prepare data set for bkg fit
    databkg=new RooDataSet("databkg","databkg",RooArgSet(*xbkg,*y),Import(*tree));
    databkg->Print() ;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::fitBackground(Int_t opt)
{
    // Fit and get positive backgrounds
    //!*****************************************
    //! Prepare and fit negative background function
    //! *****************************************
    // bkg parameters
    bkg1nratio=new RooRealVar("bkg1nratio","bkg1nratio",0.5,0,1) ;
    bkg2nratio=new RooRealVar("bkg2nratio","bkg2nratio",0.5,0,1) ;

    bkgmodelneg= new fitFbkgflat("bkgmodel","bkgmodel",*xbkg,*y,*bkg1nratio,*bkg2nratio);
    bkgmodelneg->fitTo(*databkg,NumCPU(ncpu),Save()) ;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
