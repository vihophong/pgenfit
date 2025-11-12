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

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/FrequentistCalculator.h"
#include "TGraphAsymmErrors.h"
using namespace RooFit;
using namespace RooStats;


#define LONG_FIT_RANGE 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unbinfit::unbinfit()
{
    ffitopt=0;// only 0 - fix parameter; 1 - constrain parameters

    chiSquareNDF=0;
    chiSquareNDF1n=0;
    chiSquareNDF2n=0;
    fnentrieslimit=ENTRYLIMIT;

    p_deadtime=STARTFIT;
#ifdef LONG_FIT_RANGE
    p_timerange=400;
#else
    p_timerange=10;
#endif
    nbinsHB=200;
    nbinsHSB=200;
    nbinsHSB2=200;

    ncpu=NCPUS_UNBINFIT;

    finputData=new char[1000];
    finputParms=new char[1000];

    fdecaypath=new decaypath;

    fStopWatch=new TStopwatch;
    fMCGenTime=0;
    fFitTime=0;

    seedno=4357;
    fnMC=0;

    for (int i=0;i<kmaxparms;i++)ipVal[i]=i;
    fitStatus=-9999;
    fitCovQual=-9999;
    fitNumInvalidNLL=-9999;
    fitEdm=-9999;
    fitMinNll-9999;

    plotrangelow=-1.5;
    plotrangehi=5;
#ifdef LONG_FIT_RANGE
    plotrangelow=-100;
    plotrangehi=400;
#else
    plotrangehi=10;
#endif
    fmineffMC = 0.475349;
    fmaxeffMC = 0.664527;

    significance = -9999;
    upperLimit = -9999;
    fitres = 0;
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
    //rseed=new TRandom3();
    rseed=new TRandom3(seedno);

    rseedA = new asymGausRandom();
    rseedA->Init(seedno);

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

    const char* evn_nentrieslimit = std::getenv("ENTRYLIMIT");
    if (evn_nentrieslimit){
        fnentrieslimit = atoi(evn_nentrieslimit);
        cout<<"ENTRYLIMIT = "<<fnentrieslimit<<endl;
    }
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
    t12databkg=new RooDataSet("t12databkg","t12databkg",RooArgSet(*xbkg,*y),Import(*tree));
    t12databkg->Print() ;

    //!*****************************************
    //! Setup histograms for unbinned fit
    //! *****************************************

    tree->Draw(Form("x>>hB(%d,%f,%f)",nbinsHB,-p_timerange,p_timerange));
    tree->Draw(Form("x>>hSB(%d,%f,%f)",nbinsHSB,-p_timerange,p_timerange),"y==1");
    tree->Draw(Form("x>>hSB2(%d,%f,%f)",nbinsHSB2,-p_timerange,p_timerange),"y==2");
    char tempchar1[500];
    sprintf(tempchar1,"hB");
    hB=(TH1F*) gDirectory->Get(tempchar1);
    nsig_hB_firstbin=hB->GetBinContent(hB->FindBin(p_deadtime));
    sprintf(tempchar1,"hSB");
    hSB=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hSB2");
    hSB2=(TH1F*) gDirectory->Get(tempchar1);

    tree->Draw(Form("x>>hA(%d,%f,%f)",nbinsHB,-p_timerange,p_timerange),"z>0");

    sprintf(tempchar1,"hA");
    hA=(TH1F*) gDirectory->Get(tempchar1);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::fitBackground(Int_t opt)
{
    // Fit and get positive backgrounds
    //!*****************************************
    //! Prepare and fit negative background function
    //! *****************************************

    // bkg parameters
    Double_t ini_nnbkg=tree->Draw("",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    Double_t ini_nnbkg1n = tree->Draw("",Form("x<%f&&x>%f&&y==1",-p_deadtime,-p_timerange),"goff");
    Double_t ini_nnbkg2n = tree->Draw("",Form("x<%f&&x>%f&&y==2",-p_deadtime,-p_timerange),"goff");

    cout<<"Background ratios = "<<ini_nnbkg1n/ini_nnbkg<<"\t"<<ini_nnbkg2n/ini_nnbkg1n<<endl;

    bkg1nratio=new RooRealVar("bkg1nratio","bkg1nratio",ini_nnbkg1n/ini_nnbkg,ini_nnbkg1n/ini_nnbkg/5,ini_nnbkg1n/ini_nnbkg*5) ;
    bkg2nratio=new RooRealVar("bkg2nratio","bkg2nratio",ini_nnbkg2n/ini_nnbkg1n,ini_nnbkg2n/ini_nnbkg1n/5,ini_nnbkg2n/ini_nnbkg1n*5) ;

    RooRealVar slope1("slope1","slope1",0.,-0.1,0.1) ;
    RooRealVar slope2("slope2","slope2",0.,-0.1,0.1) ;
    RooRealVar slope3("slope3","slope3",0.,-0.1,0.1) ;
#ifdef FLAT_BACKGROUNDS
    slope1.setConstant();
    slope2.setConstant();
    slope3.setConstant();
#endif

    // bkg pdf
    bkgmodelneg= new fitFbkg("bkgmodel","bkgmodel",*xbkg,*y,*bkg1nratio,*bkg2nratio,slope1,slope2,slope3);
    bkgmodelnegT12= new RooPolynomial("bkgmodelnegT12","bkgmodelnegT12",*xbkg,RooArgList(slope1));
    // fit background
#ifdef GPUMODE
    if (ffitopt!=2)
        if (opt==0) bkgmodelneg->fitTo(*databkg,BatchMode("cuda"),Save()) ;

    bool flag_fit_data_empty = false;
    if (databkg->numEntries()==0){
        flag_fit_data_empty=true;
    }else{
        if (opt==0) {
            //            bkgmodelnegT12->fitTo(*databkg,BatchMode("cuda"),Save()) ;
            bkgmodelnegT12->fitTo(*databkg,NumCPU(ncpu),Save(kTRUE)) ;
            xframe4 = xbkg->frame(Title("all fit bkg2")) ;
            databkg->plotOn(xframe4,Binning(nbinsHB/2,-p_timerange,0),RooFit::Name("bkg0n2")) ;
            bkgmodelnegT12->plotOn(xframe4,RooFit::Name("bkg0nmodel"));
        }
    }


#else
    if (opt==0) bkgmodelneg->fitTo(*databkg,NumCPU(ncpu),Save()) ;
    if (opt==0) bkgmodelnegT12->fitTo(*databkg,NumCPU(ncpu),Save()) ;
#endif

    //!*****************************************
    //! Prepare positive background function
    //! *****************************************
    // set positive background parms
    slope1pos=new RooRealVar("slope1pos","slope1pos",-slope1.getVal(),-0.1,0.1) ;
    slope2pos=new RooRealVar("slope2pos","slope2pos",-slope2.getVal(),-0.1,0.1) ;
    slope3pos=new RooRealVar("slope3pos","slope3pos",-slope3.getVal(),-0.1,0.1) ;

    slope1pos->setError(slope1.getError());
    slope2pos->setError(slope2.getError());
    slope3pos->setError(slope3.getError());

    // bkg pdf positive
    bkgmodelpos=new fitFbkg("bkgmodelpos","bkgmodelpos",*x,*y,*bkg1nratio,*bkg2nratio,*slope1pos,*slope2pos,*slope3pos);
    bkgmodelposT12= new RooPolynomial("bkgmodelposT12","bkgmodelposT12",*x,RooArgList(*slope1pos));


    //    //! bin fit background negative
    fB_bkgneg=new TF1("fB_bkgneg","pol1",-p_timerange,0);
    fSB_bkgneg=new TF1("fSB_bkgneg","pol1",-p_timerange,0);
    fSB2_bkgneg=new TF1("fSB2_bkgneg","pol1",-p_timerange,0);
#ifdef FLAT_BACKGROUNDS
    fB_bkgneg->FixParameter(1,0.);
    fSB_bkgneg->FixParameter(1,0.);
    fSB2_bkgneg->FixParameter(1,0.);
#endif
    fB_bkgneg->FixParameter(1,slope3.getVal());
    fSB_bkgneg->FixParameter(1,slope1.getVal());
    fSB2_bkgneg->FixParameter(1,slope2.getVal());


    hB->Fit(fB_bkgneg,"LEQR0+","goff");
    hSB->Fit(fSB_bkgneg,"LEQR0+","goff");
    hSB2->Fit(fSB2_bkgneg,"LEQR0+","goff");


    fB_bkgpos=new TF1("fB_bkgpos","pol1",p_deadtime,p_timerange);
    fSB_bkgpos=new TF1("fSB_bkgpos","pol1",p_deadtime,p_timerange);
    fSB2_bkgpos=new TF1("fSB2_bkgpos","pol1",p_deadtime,p_timerange);
    fB_bkgpos->FixParameter(0,fB_bkgneg->GetParameter(0));
    fSB_bkgpos->FixParameter(0,fSB_bkgneg->GetParameter(0));
    fSB2_bkgpos->FixParameter(0,fSB2_bkgneg->GetParameter(0));
    fB_bkgpos->FixParameter(1,-fB_bkgneg->GetParameter(1));
    fSB_bkgpos->FixParameter(1,-fSB_bkgneg->GetParameter(1));
    fSB2_bkgpos->FixParameter(1,-fSB2_bkgneg->GetParameter(1));


    binfitbkgparms[0]=fB_bkgpos->GetParameter(0);
    binfitbkgparms[1]=fB_bkgpos->GetParameter(1);
    binfitbkgparms[2]=fSB_bkgpos->GetParameter(0);
    binfitbkgparms[3]=fSB_bkgpos->GetParameter(1);
    binfitbkgparms[4]=fSB2_bkgpos->GetParameter(0);
    binfitbkgparms[5]=fSB2_bkgpos->GetParameter(1);

    fA_bkgneg=new TF1("fA_bkgneg","pol0",-p_timerange,0);
    hA->Fit(fA_bkgneg,"LEQR0+","goff");
    fA_bkgpos=new TF1("fA_bkgpos","pol0",0,p_timerange);
    fA_bkgpos->FixParameter(0,fA_bkgneg->GetParameter(0));
}

//....oooOO0OOooo........ofooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::initFitParameters()
{
    //!*****************************************
    //! Initialize all parameters
    //! *****************************************
    // Initialize decay parameters
    for (int i=0;i<fdecaypath->getNMember();i++){
        p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),fdecaypath->getMember(i)->decay_lamda,fdecaypath->getMember(i)->decay_lamdalow,fdecaypath->getMember(i)->decay_lamdaup);
        pvar[i]=(RooRealVar*) p[i];
        pvar[i]->setError(fdecaypath->getMember(i)->decay_lamdaerr);
        pValErrorHi[i] = fdecaypath->getMember(i)->decay_lamdaerrhi;

        p[fdecaypath->getNMember()+i]=new RooRealVar(Form("p%d",fdecaypath->getNMember()+i),Form("p%d",fdecaypath->getNMember()+i),fdecaypath->getMember(i)->decay_p1n,fdecaypath->getMember(i)->decay_p1nlow,fdecaypath->getMember(i)->decay_p1nup);
        pvar[fdecaypath->getNMember()+i]=(RooRealVar*) p[fdecaypath->getNMember()+i];
        pvar[fdecaypath->getNMember()+i]->setError(fdecaypath->getMember(i)->decay_p1nerr);
        pValErrorHi[fdecaypath->getNMember()+i] = fdecaypath->getMember(i)->decay_p1nerrhi;

        p[fdecaypath->getNMember()*2+i]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*2+i),Form("p%d",fdecaypath->getNMember()*2+i),fdecaypath->getMember(i)->decay_p2n,fdecaypath->getMember(i)->decay_p2nlow,fdecaypath->getMember(i)->decay_p2nup);
        pvar[fdecaypath->getNMember()*2+i]=(RooRealVar*) p[fdecaypath->getNMember()*2+i];
        pvar[fdecaypath->getNMember()*2+i]->setError(fdecaypath->getMember(i)->decay_p2nerr);
        pValErrorHi[fdecaypath->getNMember()*2+i] = fdecaypath->getMember(i)->decay_p2nerrhi;

        //! population ratio, replaced with the alpha branching
        p[fdecaypath->getNMember()*3+i]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*3+i),Form("p%d",fdecaypath->getNMember()*3+i),fdecaypath->getMember(i)->decay_abr/100.,fdecaypath->getMember(i)->decay_abrlow/100.,fdecaypath->getMember(i)->decay_abrup/100.);
        pvar[fdecaypath->getNMember()*3+i]=(RooRealVar*) p[fdecaypath->getNMember()*3+i];
        pvar[fdecaypath->getNMember()*3+i]->setError(fdecaypath->getMember(i)->decay_abrerr/100.);
        if (ffitopt==2){
            p[fdecaypath->getNMember()*4+i]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*4+i),Form("p%d",fdecaypath->getNMember()*4+i),0.5,0.,10.);
            pvar[fdecaypath->getNMember()*4+i]=(RooRealVar*) p[fdecaypath->getNMember()*4+i];
            if (i>1)
                pvar[fdecaypath->getNMember()*4+i]->setConstant();
        }else{
            p[fdecaypath->getNMember()*4+i]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*4+i),Form("p%d",fdecaypath->getNMember()*4+i),fdecaypath->getMember(i)->neueff,fdecaypath->getMember(i)->neuefflow,fdecaypath->getMember(i)->neueffup);
            pvar[fdecaypath->getNMember()*4+i]=(RooRealVar*) p[fdecaypath->getNMember()*4+i];
            if (fdecaypath->getMember(i)->neuefferr<0){
                pvar[fdecaypath->getNMember()*4+i]->SetTitle((char*)"-");
                pvar[fdecaypath->getNMember()*4+i]->setError(-fdecaypath->getMember(i)->neuefferr);
            }else{
                pvar[fdecaypath->getNMember()*4+i]->setError(fdecaypath->getMember(i)->neuefferr);
            }
        }

    }
    // initialize initial activity and set fixed to 1
    p[fdecaypath->getNMember()*5]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5),Form("p%d",fdecaypath->getNMember()*5),1,0,2);
    pvar[fdecaypath->getNMember()*5]=(RooRealVar*) p[fdecaypath->getNMember()*5];
    pvar[fdecaypath->getNMember()*5]->setError(0);
    pvar[fdecaypath->getNMember()*5]->setConstant();

    // Get histograms to calculate random coincidence parameters
    char tempchar1[1000];
    sprintf(tempchar1,"hdecay");
    TH1F* hdecay=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecay1nbwd");
    TH1F* hdecay1nbwd=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecaygt0nbwd");
    TH1F* hdecaygt0nbwd=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hdecay2nbwd");
    TH1F* hdecay2nbwd=(TH1F*) gDirectory->Get(tempchar1);

    // Calculate random coincidence paramters
    Double_t n1nbwd=(Double_t) hdecay1nbwd->GetEntries();
    Double_t gt0nbwd=(Double_t) hdecaygt0nbwd->GetEntries();
    Double_t n2nbwd=(Double_t) hdecay2nbwd->GetEntries();
    Double_t nball=(Double_t) hdecay->GetEntries();

    Double_t randcoinf1n=n1nbwd/nball;
    Double_t randcoinfgt0n=gt0nbwd/nball;
    Double_t randcoinf2n=n2nbwd/nball;

    Double_t randcoinf1nerr=n1nbwd/nball*TMath::Sqrt(1/n1nbwd+1/nball);
    Double_t randcoinfgt0nerr=gt0nbwd/nball*TMath::Sqrt(1/gt0nbwd+1/nball);
    Double_t randcoinf2nerr=n2nbwd/nball*TMath::Sqrt(1/n2nbwd+1/nball);

    // Initialize random coincicence parameters
    p[fdecaypath->getNMember()*5+1]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5+1),Form("p%d",fdecaypath->getNMember()*5+1),randcoinf1n,0,1);
    p[fdecaypath->getNMember()*5+2]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5+2),Form("p%d",fdecaypath->getNMember()*5+2),randcoinfgt0n,0,1);
    p[fdecaypath->getNMember()*5+3]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5+3),Form("p%d",fdecaypath->getNMember()*5+3),randcoinf2n,0,1);
    pvar[fdecaypath->getNMember()*5+1]=(RooRealVar*) p[fdecaypath->getNMember()*5+1];
    pvar[fdecaypath->getNMember()*5+2]=(RooRealVar*) p[fdecaypath->getNMember()*5+2];
    pvar[fdecaypath->getNMember()*5+3]=(RooRealVar*) p[fdecaypath->getNMember()*5+3];

    pvar[fdecaypath->getNMember()*5+1]->setError(randcoinf1nerr);
    pvar[fdecaypath->getNMember()*5+2]->setError(randcoinfgt0nerr);
    pvar[fdecaypath->getNMember()*5+3]->setError(randcoinf2nerr);

    // Initialize correction factors for beta and neutron efficiency of parent (1n,2n) decay
    std::string line;
    std::ifstream infile(finputEffParms);
    Int_t nlinesread = 0;
    Double_t be,b1ne,b2ne,n1n2ne;
    Double_t err_be,err_b1ne,err_b2ne,err_n1n2ne, err_n1n2ne_hi;
    Int_t isvary_be,isvary_b1ne,isvary_b2ne,isvary_n1n2ne;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (line[0]=='#') continue;
        if (!(iss >> be >> err_be >> b1ne >> err_b1ne >> b2ne >> err_b2ne >> n1n2ne >> err_n1n2ne>>err_n1n2ne_hi)) break;
        nlinesread++;
    }

    if (be<0) {isvary_be=1;be=-be;}else{isvary_be=0;}
    if (b1ne<0){isvary_b1ne=1;b1ne=-b1ne;}else{isvary_b1ne=0;}
    if (b2ne<0){isvary_b2ne=1;b2ne=-b2ne;}else{isvary_b2ne=0;}
    if (n1n2ne<0){isvary_n1n2ne=1;n1n2ne=-n1n2ne;}else{isvary_n1n2ne=0;}

    std::cout<<"read-in efficiency factors:"<<std::endl;
    std::cout<<be<<"\t"<<err_be<<"\t"<<b1ne<<"\t"<<err_b1ne<<"\t"<<b2ne<<"\t"<<err_b2ne<<"\t"<<n1n2ne<<"\t"<<err_n1n2ne<<"\t"<<err_n1n2ne_hi<<"\t"<<std::endl;

    if (ffitopt==2){
        pvar[fdecaypath->getNMember()*4]->setVal(1.);
        pvar[fdecaypath->getNMember()*4]->setError(0.);
        pvar[fdecaypath->getNMember()*4]->setConstant();
        pvar[fdecaypath->getNMember()*4+1]->setVal(be);
        pvar[fdecaypath->getNMember()*4+1]->setError(err_be);
        pvar[fdecaypath->getNMember()*4+1]->setConstant();
    }

    p[fdecaypath->getNMember()*5+4]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5+4),Form("p%d",fdecaypath->getNMember()*5+4),be,0,1);
    p[fdecaypath->getNMember()*5+5]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5+5),Form("p%d",fdecaypath->getNMember()*5+5),b1ne,0,1);
    p[fdecaypath->getNMember()*5+6]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5+6),Form("p%d",fdecaypath->getNMember()*5+6),b2ne,0,1);
    p[fdecaypath->getNMember()*5+7]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5+7),Form("p%d",fdecaypath->getNMember()*5+7),n1n2ne,0,1);
    pvar[fdecaypath->getNMember()*5+4]=(RooRealVar*) p[fdecaypath->getNMember()*5+4];
    pvar[fdecaypath->getNMember()*5+5]=(RooRealVar*) p[fdecaypath->getNMember()*5+5];
    pvar[fdecaypath->getNMember()*5+6]=(RooRealVar*) p[fdecaypath->getNMember()*5+6];
    pvar[fdecaypath->getNMember()*5+7]=(RooRealVar*) p[fdecaypath->getNMember()*5+7];
    pvar[fdecaypath->getNMember()*5+4]->setError(err_be);
    pvar[fdecaypath->getNMember()*5+5]->setError(err_b1ne);
    pvar[fdecaypath->getNMember()*5+6]->setError(err_b2ne);
    if (err_n1n2ne<0){
        pvar[fdecaypath->getNMember()*5+7]->setError(-err_n1n2ne);
        pvar[fdecaypath->getNMember()*5+7]->setAsymError(-err_n1n2ne,err_n1n2ne_hi);
        pvar[fdecaypath->getNMember()*5+7]->SetTitle((char*)"-");
    }else{
        pvar[fdecaypath->getNMember()*5+7]->setError(err_n1n2ne);
        pvar[fdecaypath->getNMember()*5+7]->setAsymError(err_n1n2ne,err_n1n2ne_hi);
    }


    if (isvary_be==0) pvar[fdecaypath->getNMember()*5+4]->setConstant(kTRUE);
    if (isvary_b1ne==0) pvar[fdecaypath->getNMember()*5+5]->setConstant(kTRUE);
    if (isvary_b2ne==0) pvar[fdecaypath->getNMember()*5+6]->setConstant(kTRUE);
    if (isvary_n1n2ne==0) pvar[fdecaypath->getNMember()*5+7]->setConstant(kTRUE);


    // set background/signal counts
    Double_t nnsig=tree->Draw("",Form("x>%f&&x<%f",p_deadtime,p_timerange),"goff");
    Double_t nnbkg=tree->Draw("",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    nnsig=nnsig-nnbkg;
    if (nnsig<0)
        nnsig = 10.;
    nbkg=new RooRealVar("nbkg","nbkg",nnbkg,nnbkg/3,nnbkg*3);
    nsig=new RooRealVar("nsig","nsig",nnsig,0,nnsig*5);
    //    nsig=new RooRealVar("nsig","nsig",5.00e+04,0,nnsig*5);

    cout<<"NSIG = "<<nnsig<<endl;
    nbkg->setError(TMath::Sqrt(nnbkg));

    //    if (ffitopt==2){
    //        for (int i=0;i<fdecaypath->getNMember()*4+2;i++){
    //            pT12[i] = p[i];
    //        }
    //    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::setNormalFit()
{
    //!******************************************
    //! Decide whether paramters are fixed
    //! *****************************************
    // decay parameters
    for (int i=0;i<fdecaypath->getNMember();i++){
        if (fdecaypath->getMember(i)->is_decay_hl_fix!=0)
            pvar[i]->setConstant(kTRUE);
        if (fdecaypath->getMember(i)->is_decay_p1n_fix!=0)
            pvar[fdecaypath->getNMember()+i]->setConstant(kTRUE);
        if (fdecaypath->getMember(i)->is_decay_p2n_fix!=0)
            pvar[fdecaypath->getNMember()*2+i]->setConstant(kTRUE);
        //population ratio, replaced with the alpha branching
        if (fdecaypath->getMember(i)->is_decay_abr_fix!=0)
            pvar[fdecaypath->getNMember()*3+i]->setConstant(kTRUE);
        //        if (fdecaypath->getMember(i)->is_neueff_fix!=0)
        //            pvar[fdecaypath->getNMember()*4+i]->setConstant(kTRUE);
    }

    // others parameters
    slope1pos->setConstant();
    slope2pos->setConstant();
    slope3pos->setConstant();
    // backgrounds
    nbkg->setConstant(kTRUE);
    // background ratio and slope
    bkg1nratio->setConstant();
    bkg2nratio->setConstant();

    // random coincicence parameters
    pvar[fdecaypath->getNMember()*5+1]->setConstant(kTRUE);//randcoinf1n
    pvar[fdecaypath->getNMember()*5+2]->setConstant(kTRUE);//randcoinfgt0n
    pvar[fdecaypath->getNMember()*5+3]->setConstant(kTRUE);//randcoinf2n

    //! call set model
    setModel();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::setExernalContrainFit()
{
    //! if parameter is zero (like P2n), skip the constrains
    for (int i=0;i<fdecaypath->getNMember();i++){
        if (fdecaypath->getMember(i)->is_decay_hl_fix==2)
            pvar[i]->setConstant(kTRUE);
        if (fdecaypath->getMember(i)->is_decay_p1n_fix==2)
            pvar[fdecaypath->getNMember()+i]->setConstant(kTRUE);
        if (fdecaypath->getMember(i)->is_decay_p2n_fix==2)
            pvar[fdecaypath->getNMember()*2+i]->setConstant(kTRUE);
        if (fdecaypath->getMember(i)->is_population_ratio_fix==2)
            pvar[fdecaypath->getNMember()*3+i]->setConstant(kTRUE);
        //        if (fdecaypath->getMember(i)->is_neueff_fix==2){
        //            pvar[fdecaypath->getNMember()*4+i]->setConstant(kTRUE);
        //        }
    }
    // set slope parameters to constant for the moment
    //!******************************************
    //! Define roogaussian for error propagation
    //! *****************************************
    externalconstrains=new RooArgSet;

    RooGaussian* nbkgconstr=new RooGaussian("nbkgconstr","nbkgconstr",*nbkg,RooConst(nbkg->getVal()),RooConst(nbkg->getError()));
    RooGaussian* bkg1nratiocnstr=new RooGaussian("bkg1nratiocnstr","bkg1nratiocnstr",*bkg1nratio,RooConst(bkg1nratio->getVal()),RooConst(bkg1nratio->getError()));
    RooGaussian* bkg2nratiocnstr=new RooGaussian("bkg2nratiocnstr","bkg2nratiocnstr",*bkg2nratio,RooConst(bkg2nratio->getVal()),RooConst(bkg2nratio->getError()));

    RooGaussian* slope1poscnstr=new RooGaussian("slope1poscnstr","slope1poscnstr",*slope1pos,RooConst(slope1pos->getVal()),RooConst(slope1pos->getError()));
    RooGaussian* slope2poscnstr=new RooGaussian("slope2poscnstr","slope2poscnstr",*slope2pos,RooConst(slope2pos->getVal()),RooConst(slope2pos->getError()));
    RooGaussian* slope3poscnstr=new RooGaussian("slope3poscnstr","slope3poscnstr",*slope3pos,RooConst(slope3pos->getVal()),RooConst(slope3pos->getError()));

    // decay paramters
    RooGaussian* pconstr[fdecaypath->getNMember()*5+4];
    for (int i=0;i<fdecaypath->getNMember();i++){
        pconstr[i]=new RooGaussian(Form("p%dconstr",i),Form("p%dconstr",i),*p[i],RooConst(pvar[i]->getVal()),RooConst(fdecaypath->getMember(i)->decay_lamdaerr));
        pconstr[fdecaypath->getNMember()+i]=new RooGaussian(Form("p%dconstr",fdecaypath->getNMember()+i),Form("p%dconstr",fdecaypath->getNMember()+i),*p[fdecaypath->getNMember()+i],RooConst(pvar[fdecaypath->getNMember()+i]->getVal()),RooConst(fdecaypath->getMember(i)->decay_p1nerr));
        pconstr[fdecaypath->getNMember()*2+i]=new RooGaussian(Form("p%dconstr",fdecaypath->getNMember()*2+i),Form("p%dconstr",fdecaypath->getNMember()*2+i),*p[fdecaypath->getNMember()*2+i],RooConst(pvar[fdecaypath->getNMember()*2+i]->getVal()),RooConst(fdecaypath->getMember(i)->decay_p2nerr));
        pconstr[fdecaypath->getNMember()*3+i]=new RooGaussian(Form("p%dconstr",fdecaypath->getNMember()*3+i),Form("p%dconstr",fdecaypath->getNMember()*3+i),*p[fdecaypath->getNMember()*3+i],RooConst(pvar[fdecaypath->getNMember()*3+i]->getVal()),RooConst(fdecaypath->getMember(i)->population_ratioerr));
        pconstr[fdecaypath->getNMember()*4+i]=new RooGaussian(Form("p%dconstr",fdecaypath->getNMember()*4+i),Form("p%dconstr",fdecaypath->getNMember()*4+i),*p[fdecaypath->getNMember()*4+i],RooConst(pvar[fdecaypath->getNMember()*4+i]->getVal()),RooConst(0.1));
    }
    // random coincidence paramters
    pconstr[fdecaypath->getNMember()*5+1]=new RooGaussian(Form("p%dconstr",fdecaypath->getNMember()*5+1),Form("p%dconstr",fdecaypath->getNMember()*5+1),*p[fdecaypath->getNMember()*5+1],RooConst(pvar[fdecaypath->getNMember()*5+1]->getVal()),RooConst(pvar[fdecaypath->getNMember()*5+1]->getError()));
    pconstr[fdecaypath->getNMember()*5+2]=new RooGaussian(Form("p%dconstr",fdecaypath->getNMember()*5+2),Form("p%dconstr",fdecaypath->getNMember()*5+2),*p[fdecaypath->getNMember()*5+2],RooConst(pvar[fdecaypath->getNMember()*5+2]->getVal()),RooConst(pvar[fdecaypath->getNMember()*5+2]->getError()));
    pconstr[fdecaypath->getNMember()*5+3]=new RooGaussian(Form("p%dconstr",fdecaypath->getNMember()*5+3),Form("p%dconstr",fdecaypath->getNMember()*5+3),*p[fdecaypath->getNMember()*5+3],RooConst(pvar[fdecaypath->getNMember()*5+3]->getVal()),RooConst(pvar[fdecaypath->getNMember()*5+3]->getError()));

    // add to contrains set
    for (int i=0;i<fdecaypath->getNMember()*5;i++) {
        if (!pvar[i]->isConstant()&&i!=0&&i!=fdecaypath->getNMember()&&i!=fdecaypath->getNMember()*2) // let decay parameters of parent nuclei free
            externalconstrains->add(*pconstr[i]);
    }

    // backgrounds parameters
    externalconstrains->add(*nbkgconstr);
    externalconstrains->add(*bkg1nratiocnstr);
    externalconstrains->add(*bkg2nratiocnstr);

    //! fix slope parameters for the moment
    //externalconstrains->add(*slope1poscnstr);
    //externalconstrains->add(*slope2poscnstr);
    //externalconstrains->add(*slope3poscnstr);
    slope1pos->setConstant();
    slope2pos->setConstant();
    slope3pos->setConstant();

    //random coincicence parameters
    externalconstrains->add(*pconstr[fdecaypath->getNMember()*5+1]);
    externalconstrains->add(*pconstr[fdecaypath->getNMember()*5+2]);
    externalconstrains->add(*pconstr[fdecaypath->getNMember()*5+3]);

    //! call set model
    setModel();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::setModel()
{
    //!******************************************
    //! Construct final fit model
    //! *****************************************
    //!
    totdecaymodel=new fitF("totdecaymodel","totdecaymodel",*x,*y,p);
    final_pdf=new RooAddPdf("final_pdf","final pdf",RooArgList(*totdecaymodel,*bkgmodelpos),RooArgList(*nsig,*nbkg));
    totdecaymodelT12=new fitF_T12("totdecaymodelT12","totdecaymodelT12",*x,p);
    final_pdfT12=new RooAddPdf("final_pdfT12","final_pdfT12 ",RooArgList(*totdecaymodelT12,*bkgmodelposT12),RooArgList(*nsig,*nbkg));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::prepareData()
{
    //!******************************************
    //! Prepare data set for fitting forward correlated data
    //! *****************************************
    data=new RooDataSet("data","data",RooArgSet(*x,*y),Import(*tree)) ;
    data->Print() ;
    t12data=new RooDataSet("t12data","t12data",RooArgSet(*x),Import(*tree)) ;
    t12data->Print() ;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::bookOutputTree()
{
    fout=new TFile(foutputData,"recreate");
    foutputtree=new TTree("treeout","treeout");
    foutputtree->Branch("pVal",pVal,Form("pVal[%d]/D",kmaxparms));
    foutputtree->Branch("pCentralVal",pCentralVal,Form("pCentralVal[%d]/D",kmaxparms));
    foutputtree->Branch("pValError",pValError,Form("pValError[%d]/D",kmaxparms));
    foutputtree->Branch("ispVary",ispVary,Form("ispVary[%d]/I",kmaxparms));
    foutputtree->Branch("ipVal",ipVal,Form("ipVal[%d]/I",kmaxparms));

    foutputtree->Branch("nsigCentralVal",&nsigCentralVal,"nsigCentralVal/D");
    foutputtree->Branch("nsigVal",&nsigVal,"nsigVal/D");
    foutputtree->Branch("nbkgVal",&nbkgVal,"nbkgVal/D");
    foutputtree->Branch("nbkgError",&nbkgError,"nbkgError/D");
    foutputtree->Branch("bkg1nratioVal",&bkg1nratioVal,"bkg1nratioVal/D");
    foutputtree->Branch("bkg1nratioError",&bkg1nratioError,"bkg1nratioError/D");
    foutputtree->Branch("bkg2nratioVal",&bkg2nratioVal,"bkg2nratioVal/D");
    foutputtree->Branch("bkg2nratioError",&bkg2nratioError,"bkg2nratioError/D");
    foutputtree->Branch("slope1posVal",&slope1posVal,"slope1posVal/D");
    foutputtree->Branch("slope1posVal",&slope1posVal,"slope1posVal/D");
    foutputtree->Branch("slope1posVal",&slope1posVal,"slope1posVal/D");
    foutputtree->Branch("slope1posError",&slope1posError,"slope1posError/D");
    foutputtree->Branch("slope2posError",&slope2posError,"slope2posError/D");
    foutputtree->Branch("slope3posError",&slope3posError,"slope3posError/D");

    foutputtree->Branch("fFitTime",&fFitTime,"fFitTime/D");

    foutputtree->Branch("fitStatus",&fitStatus,"fitStatus/I");
    foutputtree->Branch("fitCovQual",&fitCovQual,"fitCovQual/I");
    foutputtree->Branch("fitNumInvalidNLL",&fitNumInvalidNLL,"fitNumInvalidNLL/I");
    foutputtree->Branch("fitEdm",&fitEdm,"fitEdm/D");
    foutputtree->Branch("fitMinNll",&fitMinNll,"fitMinNll/D");
    foutputtree->Branch("chiSquareNDF",&chiSquareNDF,"chiSquareNDF/D");
    foutputtree->Branch("chiSquareNDF1n",&chiSquareNDF1n,"chiSquareNDF1n/D");
    foutputtree->Branch("chiSquareNDF2n",&chiSquareNDF2n,"chiSquareNDF2n/D");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::prepareMonteCarloData(int nevents)
{
    //!******************************************
    //! Toy MC dataset (just for testing)
    //! *****************************************
    fStopWatch->Start();
    data = final_pdf->generate(RooArgSet(*x,*y),nevents) ;
    data->Print() ;
    fStopWatch->Stop();
    fMCGenTime=fStopWatch->RealTime();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::getParameters()
{
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        pVal[i]=pvar[i]->getVal();
        pValError[i]=pvar[i]->getError();
    }
    nsigVal=nsig->getVal();
    nsigError=nsig->getError();

    nbkgVal=nbkg->getVal();
    nbkgError=nbkg->getError();
    bkg1nratioVal=bkg1nratio->getVal();
    bkg1nratioError=bkg1nratio->getError();
    bkg2nratioVal=bkg2nratio->getVal();
    bkg2nratioError=bkg2nratio->getError();

    slope1posVal=slope1pos->getVal();
    slope2posVal=slope2pos->getVal();
    slope3posVal=slope3pos->getVal();

    //slope1posError=slope1pos->getError();
    //slope2posError=slope2pos->getError();
    //slope3posError=slope3pos->getError();
    slope1posError=0;
    slope2posError=0;
    slope3posError=0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::setCentralParameters()
{
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        pCentralVal[i]=pvar[i]->getVal();
        if (pvar[i]->isConstant()) ispVary[i]=0;
        else ispVary[i]=1;
    }
    nsigCentralVal=nsig->getVal();
    nbkgCentralVal=nbkg->getVal();
    bkg1nratioCentralVal=bkg1nratio->getVal();
    bkg2nratioCentralVal=bkg2nratio->getVal();

    slope1posCentralVal=slope1pos->getVal();
    slope2posCentralVal=slope2pos->getVal();
    slope3posCentralVal=slope3pos->getVal();

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::printCurrentParameters()
{
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        if (i<fdecaypath->getNMember())
            std::cout<<"lambda"<<std::endl;
        else if (i<fdecaypath->getNMember()*2)
            std::cout<<"p1n"<<std::endl;
        else if (i<fdecaypath->getNMember()*3)
            std::cout<<"p2n"<<std::endl;
        else if (i<fdecaypath->getNMember()*4)
            std::cout<<"alphabranch"<<std::endl;
        else if (i<fdecaypath->getNMember()*5)
            std::cout<<"neueff"<<std::endl;
        std::cout<<"parms"<<i<<"=\t"<<pvar[i]->getVal()<<" +/- "<<pvar[i]->getError()<<"\tisfix="<<pvar[i]->isConstant()<<std::endl;
    }
    std::cout<<"nbkg=\t"<<nbkg->getVal()<<" +/- "<<nbkg->getError()<<std::endl;
    std::cout<<"bkg1nratio=\t"<<bkg1nratio->getVal()<<" +/- "<<bkg1nratio->getError()<<std::endl;
    std::cout<<"bkg2nratio=\t"<<bkg2nratio->getVal()<<" +/- "<<bkg2nratio->getError()<<std::endl;
    std::cout<<"slope1pos=\t"<<slope1pos->getVal()<<" +/- "<<slope1pos->getError()<<std::endl;
    std::cout<<"slope2pos=\t"<<slope2pos->getVal()<<" +/- "<<slope2pos->getError()<<std::endl;
    std::cout<<"slope3pos=\t"<<slope3pos->getVal()<<" +/- "<<slope3pos->getError()<<std::endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::setValParameters()
{
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        if (ffitopt==2){//fit with half-life only
            if (i!=fdecaypath->getNMember()*4 && i!=fdecaypath->getNMember()*4+1){
                pvar[i]->setVal(pVal[i]);
            }else if (i==fdecaypath->getNMember()*4+1){
                pvar[i]->setVal(pVal[fdecaypath->getNMember()*4+1]);
            }
        }else{
            pvar[i]->setVal(pVal[i]);
        }
    }
    nbkg->setVal(nbkgVal);
    bkg1nratio->setVal(bkg1nratioVal);
    bkg2nratio->setVal(bkg2nratioVal);
    slope1pos->setVal(slope1posVal);
    slope2pos->setVal(slope2posVal);
    slope3pos->setVal(slope3posVal);

    //! for binfit
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){//decay parameters
        binfitparms[i]=pVal[i];
    }
    binfitparms[fdecaypath->getNMember()*5]=nsig_hB_firstbin;//initial activity
    binfitparms[fdecaypath->getNMember()*5+8]=binfitbkgparms[0];//bkg hdecay offset
    binfitparms[fdecaypath->getNMember()*5+9]=binfitbkgparms[1];//bkg hdecay slope
    binfitparms[fdecaypath->getNMember()*5+10]=binfitbkgparms[2];//bkg hdecay1n offset
    binfitparms[fdecaypath->getNMember()*5+11]=binfitbkgparms[3];//bkg hdecay1n slope
    binfitparms[fdecaypath->getNMember()*5+12]=binfitbkgparms[4];//bkg hdecay2n offset
    binfitparms[fdecaypath->getNMember()*5+13]=binfitbkgparms[5];//bkg hdecay2n slope
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::generateMC()
{
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        if ((pvar[i]->isConstant())&&pValError[i]!=0) {
            if(i>=fdecaypath->getNMember()*4&&i<fdecaypath->getNMember()*5){//neutron efficiency parameter
                if (strcmp(pvar[i]->GetTitle(),(char*)"-")==0){//generate uniform random distribution from mean-error_low to mean+error_high
                    pVal[i]=pCentralVal[i]-pValError[i]+rseed->Rndm()*((pCentralVal[i]+fdecaypath->getMember(i-fdecaypath->getNMember()*4)->neuefferrhi)-(pCentralVal[i]-pValError[i]));
                }else{// Special function for Asymetric error
                    pVal[i]=rseedA->generate(pCentralVal[i],pValError[i],fdecaypath->getMember(i-fdecaypath->getNMember()*4)->neuefferrhi);
                }
            }else{//other parameters
                if (i<fdecaypath->getNMember()*3){//HL, P1n and P2n errors
                    //                    pVal[i]=rseedA->generate(pCentralVal[i],pValError[i],pValErrorHi[i]);
                    if (pValErrorHi[i]==0){
                        pVal[i]=rseed->Rndm();
                        //                        cout<<"AAAA"<<pCentralVal[i]<<"\t"<<pValError[i]<<"\t"<<pValErrorHi[i]<<endl;
                    }else{
                        pVal[i]=rseedA->generate(pCentralVal[i],pValError[i],pValErrorHi[i]);
                    }
                }else{
                    if (strcmp(pvar[i]->GetTitle(),(char*)"-")==0){//generate uniform random distribution
                        if (i==fdecaypath->getNMember()*5+7)//2n efficiency of parent
                            pVal[i]=pCentralVal[i]-pValError[i]+rseed->Rndm()*(pValError[i]+pvar[i]->getAsymErrorHi());
                        else
                            pVal[i]=pCentralVal[i]-pValError[i]+rseed->Rndm()*pValError[i]*2;
                    }else{
                        if (i==fdecaypath->getNMember()*5+7){//2n efficiency of parent
                            pVal[i]=rseedA->generate(pCentralVal[i],pValError[i],pvar[i]->getAsymErrorHi());
                        }else{
                            pVal[i]=rseedA->generate(pCentralVal[i],pValError[i]);
                        }
                    }
                }
            }
        }else{
            pVal[i]=pCentralVal[i];
        }

        //#ifdef PARENT_NEUEFF_UNIFORM
        //        //! neueff of parent randomly distributed from 40 to 68%
        //        if (i==fdecaypath->getNMember()*4) pVal[i]=fmineffMC+rseed->Rndm()*(fmaxeffMC-fmineffMC);
        //#endif
    }
    nsigVal=nsigCentralVal;//reset central parameters after 1 fit


    //    nbkgVal = nbkgCentralVal;
    //    bkg1nratioVal = bkg1nratioCentralVal;
    //    bkg2nratioVal = bkg2nratioCentralVal;
    //nbkgVal=rseed->Gaus(nbkgCentralVal,nbkgError);
    //bkg1nratioVal=rseed->Gaus(bkg1nratioCentralVal,bkg1nratioError);
    //bkg2nratioVal=rseed->Gaus(bkg2nratioCentralVal,bkg2nratioError);
    if (nbkgError>0&&bkg1nratioError>0&&bkg2nratioError>0){
        nbkgVal=rseedA->generate(nbkgCentralVal,nbkgError);
        bkg1nratioVal=rseedA->generate(bkg1nratioCentralVal,bkg1nratioError);
        bkg2nratioVal=rseedA->generate(bkg2nratioCentralVal,bkg2nratioError);
    }

    slope1posVal=slope1posCentralVal;
    slope2posVal=slope2posCentralVal;
    slope3posVal=slope3posCentralVal;
    //slope1posVal=rseed->Gaus(slope1posCentralVal,slope1posError);
    //slope2posVal=rseed->Gaus(slope2posCentralVal,slope2posError);
    //slope3posVal=rseed->Gaus(slope3posCentralVal,slope3posError);

    //slope1posVal=rseedA->generate(slope1posCentralVal,slope1posError);
    //slope2posVal=rseedA->generate(slope2posCentralVal,slope2posError);
    //slope3posVal=rseedA->generate(slope3posCentralVal,slope3posError);

    //! for binfit
    //binfitbkgparms[0]=rseed->Gaus(fB_bkgpos->GetParameter(0),fB_bkgneg->GetParError(0));
    //binfitbkgparms[2]=rseed->Gaus(fSB_bkgpos->GetParameter(0),fSB_bkgneg->GetParError(0));
    //binfitbkgparms[4]=rseed->Gaus(fSB2_bkgpos->GetParameter(0),fSB2_bkgneg->GetParError(0));
    //binfitbkgparms[1]=rseed->Gaus(fB_bkgpos->GetParameter(1),fB_bkgneg->GetParError(1));
    //binfitbkgparms[3]=rseed->Gaus(fSB_bkgpos->GetParameter(1),fSB_bkgneg->GetParError(1));
    //binfitbkgparms[5]=rseed->Gaus(fSB2_bkgpos->GetParameter(1),fSB2_bkgneg->GetParError(1));

    binfitbkgparms[0]=rseedA->generate(fB_bkgpos->GetParameter(0),fB_bkgneg->GetParError(0));
    binfitbkgparms[2]=rseedA->generate(fSB_bkgpos->GetParameter(0),fSB_bkgneg->GetParError(0));
    binfitbkgparms[4]=rseedA->generate(fSB2_bkgpos->GetParameter(0),fSB2_bkgneg->GetParError(0));
    //    binfitbkgparms[1]=rseedA->generate(fB_bkgpos->GetParameter(1),fB_bkgneg->GetParError(1));
    //    binfitbkgparms[3]=rseedA->generate(fSB_bkgpos->GetParameter(1),fSB_bkgneg->GetParError(1));
    //    binfitbkgparms[5]=rseedA->generate(fSB2_bkgpos->GetParameter(1),fSB2_bkgneg->GetParError(1));
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::doFit()
{
    //!******************************************
    //! Perform the fit
    //! *****************************************
    fStopWatch->Clear();
    fStopWatch->Start();

    if (ffitopt==0){
#ifdef GPUMODE
        fitres=final_pdf->fitTo(*data,BatchMode("cuda"),Save(kTRUE),PrintLevel(3));
#else
        fitres=final_pdf->fitTo(*data,NumCPU(ncpu),Save(kTRUE),PrintLevel(3));
#endif
    }else if (ffitopt==1){//with external constrain
#ifdef GPUMODE
        fitres=final_pdf->fitTo(*data,ExternalConstraints(*externalconstrains),BatchMode("cuda"),Save(kTRUE),PrintLevel(3));
#else
        fitres=final_pdf->fitTo(*data,ExternalConstraints(*externalconstrains),NumCPU(ncpu),Save(kTRUE),PrintLevel(3));
#endif
    }else{//fit T1/2 only, gpu mode not support for now

        fitres=final_pdfT12->fitTo(*t12data,NumCPU(ncpu),Save(kTRUE),PrintLevel(3));
        //        fitres=final_pdfT12->fitTo(*t12data,BatchMode("cuda"));
    }
    fStopWatch->Stop();
    fFitTime=fStopWatch->RealTime();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::writeResultsMC()
{
    getParameters();
    if (fitres){
        fitStatus=fitres->status();
        fitCovQual=fitres->covQual();
        fitNumInvalidNLL=fitres->numInvalidNLL();
        fitEdm=fitres->edm();
        fitMinNll=fitres->minNll();
    }
    foutputtree->Fill();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::plotResults()
{
    //!******************************************
    //! Outputs
    //! *****************************************
    TCanvas *c1 = new TCanvas("c1","c1",1200, 800);
    c1->Divide(2,3);
    c1->cd(1);
    RooPlot* xframe0 = x->frame(Title("all fit")) ;
    if (ffitopt==2){
        t12data->plotOn(xframe0,Binning(nbinsHB/2,0.,p_timerange),RooFit::Name("data0n")) ;
        binw=p_timerange*2/nbinsHB;
        final_pdfT12->plotOn(xframe0,RooFit::Name("data0nmodel")) ;
        final_pdfT12->plotOn(xframe0,Components("bkgmodelposT12"),RooFit::Name("bkg0nPos")) ;
    }else{
        data->plotOn(xframe0,Binning(nbinsHB/2,0.,p_timerange),RooFit::Name("data0n")) ;
        binw=p_timerange*2/nbinsHB;
        final_pdf->plotOn(xframe0,RooFit::Name("data0nmodel")) ;
    }
    model0nCurve=(RooCurve*)xframe0->getCurve("data0nmodel");
    model0nHist=(RooHist*)xframe0->getHist("data0n");
    modelbkg0nCurvePositive=(RooCurve*)xframe0->getCurve("bkg0nPos");
    modelbkg0nCurvePositive->SetLineWidth(0);
    xframe0->Draw() ;
    //    c1->cd(3);
    //    RooPlot* xframe1 = x->frame(Title("1 neutron fit")) ;
    //    data->plotOn(xframe1,Cut("y==y::1neu"),Binning(nbinsHSB/2),RooFit::Name("data1n")) ;
    //    final_pdf->plotOn(xframe1,Slice(*y,"1neu"),RooFit::Name("data1nmodel")) ;
    //    //bkgmodelpos->plotOn(xframe1,Slice(*y,"1neu"),RooFit::Name("bkg1nposmodel")) ;
    //    xframe1->Draw() ;
    //    c1->cd(5);
    //    RooPlot* xframe2 = x->frame(Title("2 neutron fit")) ;
    //    data->plotOn(xframe2,Cut("y==y::2neu"),Binning(nbinsHSB2/2),RooFit::Name("data2n")) ;
    //    final_pdf->plotOn(xframe2,Slice(*y,"2neu"),RooFit::Name("data2nmodel")) ;
    //    //bkgmodelpos->plotOn(xframe2,Slice(*y,"2neu"),RooFit::Name("bkg2nposmodel")) ;
    //    xframe2->Draw() ;

    c1->cd(2);
    //    RooPlot* xframe3 = xbkg->frame(Title("all fit bkg")) ;
    //    databkg->plotOn(xframe3,Binning(nbinsHB/2),RooFit::Name("bkg0n")) ;

    //    bkgmodelnegT12->plotOn(xframe3,RooFit::Name("bkg0nmodel"));
    //    xframe3->Draw() ;
    modelbkg0nHist=(RooHist*)xframe4->getHist("bkg0n2");
    modelbkg0nCurve=(RooCurve*)xframe4->getCurve("bkg0nmodel");
    for (int i=0;i<modelbkg0nCurve->GetN();i++){
        modelbkg0nCurve->SetPoint(i,modelbkg0nCurve->GetPointX(i),modelbkg0nCurvePositive->GetPointY(i));
    }
    //    modelbkg0nCurve->SetLineWidth(0);
    xframe4->Draw();

    //    c1->cd(4);
    //    RooPlot* xframe4 = xbkg->frame(Title("1 neutron fit bkg")) ;
    //    databkg->plotOn(xframe4,Cut("y==y::1neu"),Binning(nbinsHSB/2),RooFit::Name("bkg1n")) ;
    //    xframe4->Draw() ;

    //    c1->cd(6);
    //    RooPlot* xframe5 = xbkg->frame(Title("2 neutron fit bkg")) ;
    //    databkg->plotOn(xframe5,Cut("y==y::2neu"),Binning(nbinsHSB2/2),RooFit::Name("bkg2n")) ;
    //    xframe5->Draw() ;
    c1->Write();

    //    model1nCurve=(RooCurve*)xframe1->getCurve("data1nmodel");
    //    model1nHist=(RooHist*)xframe1->getHist("data1n");
    //    model2nCurve=(RooCurve*)xframe2->getCurve("data2nmodel");
    //    model2nHist=(RooHist*)xframe2->getHist("data2n");

    //    modelbkg1nHist=(RooHist*)xframe3->getHist("bkg1n");
    //    modelbkg2nHist=(RooHist*)xframe5->getHist("bkg2n");

    model0nHist->Write();
    //    model1nHist->Write();
    //    model2nHist->Write();
    model0nCurve->Write();
    //    model1nCurve->Write();
    //    model2nCurve->Write();

    modelbkg0nHist->Write();\
    //    modelbkg1nHist->Write();
    //    modelbkg2nHist->Write();


    //    //!******************************************
    //    //! Plotting components using bin fit
    //    //! *****************************************
    //    //! fix bkg parms for binned model from result of backward fit
    Double_t bkg1nratioval=bkg1nratio->getVal();
    Double_t bkg2nratioval=bkg2nratio->getVal();
    Double_t slope1posval=slope1pos->getVal();
    Double_t slope2posval=slope2pos->getVal();
    Double_t slope3posval=slope3pos->getVal();
    Double_t a0=slope3posval;
    Double_t b0=(2*nbkg->getVal()-a0*p_timerange*p_timerange)/2/p_timerange*(p_timerange/nbinsHB*2);
    b0 = modelbkg0nCurve->GetPointY(1);
    Double_t a1=slope1posval*bkg1nratioval;
    Double_t b1=(2*nbkg->getVal()*bkg1nratioval-a1*p_timerange*p_timerange)/2/p_timerange*(p_timerange/nbinsHB*2);
    Double_t a2=slope2posval*bkg1nratioval*bkg2nratioval;
    Double_t b2=(2*nbkg->getVal()*bkg1nratioval*bkg2nratioval-a2*p_timerange*p_timerange)/2/p_timerange*(p_timerange/nbinsHB*2);

    fB_bkgpos->FixParameter(0,b0);
    //    fSB_bkgpos->FixParameter(0,b1);
    //    fSB2_bkgpos->FixParameter(0,b2);
    fB_bkgpos->FixParameter(1,a0);
    //    fSB_bkgpos->FixParameter(1,a1);
    //    fSB2_bkgpos->FixParameter(1,a2);

    fB_bkgneg->FixParameter(0,b0);
    //    fSB_bkgneg->FixParameter(0,b1);
    //    fSB2_bkgneg->FixParameter(0,b2);
    //    fB_bkgneg->FixParameter(1,-a0);
    //    fB_bkgneg->FixParameter(1,-a1);
    //    fB_bkgneg->FixParameter(1,-a2);

    totdecaymodelforplot=new fitF_T12("totdecaymodelforplot","totdecaymodelforplot",*x,p);
    totdecaymodelforplot->initPath();
    fB=new TF1("fB",totdecaymodelforplot,&fitF_T12::fcndecay,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay");
    //    //! set fix parameters
    nsig_hB_firstbin=model0nCurve->Eval(p_deadtime)-fB_bkgpos->Eval(p_deadtime);
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        if (i!=fdecaypath->getNMember()*4) {
            fB->FixParameter(i,pvar[i]->getVal());
        }else{
            fB->SetParameter(fdecaypath->getNMember()*4,nsig_hB_firstbin);
            fB->SetParLimits(fdecaypath->getNMember()*4,nsig_hB_firstbin/10,nsig_hB_firstbin*10);
        }
    }

    //    b0 = modelbkg0nCurve->GetPointY(0);
    fB->FixParameter(fdecaypath->getNMember()*5+8,b0);
    fB->FixParameter(fdecaypath->getNMember()*5+9,a0);
    //    fB->Write();
    //    fSB=new TF1("fSB",totdecaymodelforplot,&fitF_T12::fcndecay1n,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay1n");
    //    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
    //        fSB->FixParameter(i,pvar[i]->getVal());
    //    }
    //    fSB->FixParameter(fdecaypath->getNMember()*5+8,b1);
    //    fSB->FixParameter(fdecaypath->getNMember()*5+9,a1);

    //    fSB2=new TF1("fSB2",totdecaymodelforplot,&fitF_T12::fcndecay2n,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay2n");
    //    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
    //        fSB2->FixParameter(i,pvar[i]->getVal());
    //    }
    //    fSB2->FixParameter(fdecaypath->getNMember()*5+8,b2);
    //    fSB2->FixParameter(fdecaypath->getNMember()*5+9,a2);

    model0nCurve->Fit(fB,"LEQR","goff");
    model0nCurve->GetFunction("fB")->SetLineWidth(0);

    //    //fB->FixParameter(fdecaypath->getNMember()*5,nsig_hB_firstbin);
    //    fSB->FixParameter(fdecaypath->getNMember()*5,fB->GetParameter(fdecaypath->getNMember()*5));
    //    fSB2->FixParameter(fdecaypath->getNMember()*5,fB->GetParameter(fdecaypath->getNMember()*5));

    //    //! construct parent/daugters decay components for plotting demonstration
    fB_parent=new TF1("fB_parent",totdecaymodelforplot,&fitF_T12::fcndecay_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay_parent");
    fB_parentnobkg=new TF1("fB_parentnobkg",totdecaymodelforplot,&fitF_T12::fcndecay_parentnobkg,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay_parentnobkg");
    fB_daugter=new TF1("fB_daugter",totdecaymodelforplot,&fitF_T12::fcndecay_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay_daugter");
    fB_daugternobkg=new TF1("fB_daugternobkg",totdecaymodelforplot,&fitF_T12::fcndecay_daugternobkg,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay_daugternobkg");
    fB_alpha=new TF1("fB_alpha",totdecaymodelforplot,&fitF_T12::fcndecayAlpha,0,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecayAlpha");
    fB_alphanobkg=new TF1("fB_alphanobkg",totdecaymodelforplot,&fitF_T12::fcndecayAlphanobkg,0,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecayAlphanobkg");
    //    fB_parent->Write();
    //    fB->Write();
    //    fSB_parent=new TF1("fSB_parent",totdecaymodelforplot,&fitF_T12::fcndecay1n_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay1n_parent");
    //    fSB_daugter=new TF1("fSB_daugter",totdecaymodelforplot,&fitF_T12::fcndecay1n_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay1n_daugter");

    //    fSB2_parent=new TF1("fSB2_parent",totdecaymodelforplot,&fitF_T12::fcndecay2n_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay2n_parent");
    //    fSB2_daugter=new TF1("fSB2_daugter",totdecaymodelforplot,&fitF_T12::fcndecay2n_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay2n_daugter");

    //    fSB_c1=new TF1("fSB_c1",totdecaymodelforplot,&fitF_T12::fcndecay1n_c1,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay1n_c1");
    //    fSB_c2=new TF1("fSB_c2",totdecaymodelforplot,&fitF_T12::fcndecay1n_c2,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay1n_c2");
    //    fSB_c3=new TF1("fSB_c3",totdecaymodelforplot,&fitF_T12::fcndecay1n_c3,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay1n_c3");
    //    fSB_c23=new TF1("fSB_c23",totdecaymodelforplot,&fitF_T12::fcndecay1n_c23,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay1n_c23");
    //    fSB2_c1=new TF1("fSB2_c1",totdecaymodelforplot,&fitF_T12::fcndecay2n_c1,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay2n_c1");
    //    fSB2_c2=new TF1("fSB2_c2",totdecaymodelforplot,&fitF_T12::fcndecay2n_c2,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay2n_c2");
    //    fSB2_c3=new TF1("fSB2_c3",totdecaymodelforplot,&fitF_T12::fcndecay2n_c3,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay2n_c3");
    //    fSB2_c4=new TF1("fSB2_c4",totdecaymodelforplot,&fitF_T12::fcndecay2n_c4,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay2n_c4");
    //    fSB2_c134=new TF1("fSB2_c134",totdecaymodelforplot,&fitF_T12::fcndecay2n_c134,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF_T12","fcndecay2n_c134");

    writeFitComponents();
    plotResultsMore();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::writeResults()
{
    char tempstr[500];
    sprintf(tempstr,"%s.txt",foutputData);
    std::ofstream ofs(tempstr,std::ios::app);
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        if (!pvar[i]->isConstant())
            ofs<<i<<"\t"<<pvar[i]->getVal()<<"\t"<<pvar[i]->getError()<<std::endl;
    }
    ofs<<"nsig = "<<nsig->getVal()<<"\tnbkg = "<<nbkg->getVal()<<std::endl;
    ofs<<"chisquare/NDF = "<<chiSquareNDF<<"\t"<<chiSquareNDF1n<<"\t"<<chiSquareNDF2n<<std::endl;
    ofs<<"FitTime = "<<fFitTime<<endl;
    if (fitres)
        fitres->Print();
    std::cout<<pvar[0]->getVal()<<std::endl;
    std::cout<<"Time for MC generation = "<<fMCGenTime<<std::endl;
    std::cout<<"Time for Fitting = "<<fFitTime<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::RunBinFit()
{
    fitBackground(1);
    initFitParameters();
    if (ffitopt==0)
        setNormalFit();
    else
        setExernalContrainFit();
    printCurrentParameters();

    prepareData();
    bookOutputTree();
    //! setup stuffs for MC fits
    setCentralParameters();
    getParameters();


    totdecaymodel->initPath();
    fB=new TF1("fB",totdecaymodel,&fitF::fcndecay,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay");
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        fB->FixParameter(i,pCentralVal[i]);
    }
    fB->FixParameter(fdecaypath->getNMember()*5+8,100);
    fB->FixParameter(fdecaypath->getNMember()*5+9,0);


    fSB=new TF1("fSB",totdecaymodel,&fitF::fcndecay1n,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay1n");
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        fSB->FixParameter(i,pCentralVal[i]);
    }
    fSB->FixParameter(fdecaypath->getNMember()*5+8,10);
    fSB->FixParameter(fdecaypath->getNMember()*5+9,0);

    fSB2=new TF1("fSB2",totdecaymodel,&fitF::fcndecay2n,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay2n");
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        fSB2->FixParameter(i,pCentralVal[i]);
    }
    fSB2->FixParameter(fdecaypath->getNMember()*5+8,1);
    fSB2->FixParameter(fdecaypath->getNMember()*5+9,0);

    //! set fix parameters
    fB->FixParameter(fdecaypath->getNMember()*5,nsig_hB_firstbin);
    fSB->FixParameter(fdecaypath->getNMember()*5,nsig_hB_firstbin);
    fSB2->FixParameter(fdecaypath->getNMember()*5,nsig_hB_firstbin);
    cout<<"Eval fB = "<<fB->Eval(p_timerange/2+p_deadtime/2)<<endl;
    cout<<"Eval fSB = "<<fSB->Eval(p_timerange/2+p_deadtime/2)<<endl;
    cout<<"Eval fSB2 = "<<fSB2->Eval(p_timerange/2+p_deadtime/2)<<endl;

    ROOT::Math::WrappedMultiTF1 wfB(*fB,1);
    ROOT::Math::WrappedMultiTF1 wfSB(*fSB,1);
    ROOT::Math::WrappedMultiTF1 wfSB2(*fSB2,1);

    ROOT::Fit::DataOptions opt;
    //! limit within the fitting range
    opt.fUseRange  =true;

    ROOT::Fit::DataRange rangeB;
    rangeB.SetRange(p_deadtime,p_timerange);
    ROOT::Fit::BinData dataB(opt,rangeB);
    ROOT::Fit::FillData(dataB, hB);

    ROOT::Fit::DataRange rangeSB;
    rangeSB.SetRange(p_deadtime,p_timerange);
    ROOT::Fit::BinData dataSB(opt,rangeSB);
    ROOT::Fit::FillData(dataSB, hSB);

    ROOT::Fit::DataRange rangeSB2;
    rangeSB2.SetRange(p_deadtime,p_timerange);
    ROOT::Fit::BinData dataSB2(opt,rangeSB2);
    ROOT::Fit::FillData(dataSB2, hSB2);


    ROOT::Fit::PoissonLLFunction chi2_B(dataB, wfB);
    ROOT::Fit::PoissonLLFunction chi2_SB(dataSB, wfSB);
    ROOT::Fit::PoissonLLFunction chi2_SB2(dataSB2, wfSB2);

    GlobalChi2 globalChi2(chi2_B, chi2_SB, chi2_SB2);

    ROOT::Fit::Fitter fitter;

    setValParameters();//initiate parameters
    fitter.Config().SetParamsSettings(fdecaypath->getNMember()*5+14,binfitparms);
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        if (pvar[i]->isConstant()&&i!=fdecaypath->getNMember()*5){
            fitter.Config().ParSettings(i).Fix();
        }else{
            if (i==fdecaypath->getNMember()*5)
                fitter.Config().ParSettings(i).SetLimits(0,nsig_hB_firstbin*20);
            else
                fitter.Config().ParSettings(i).SetLimits(pvar[i]->getMin(),pvar[i]->getMax());
        }
    }
    for (int i=fdecaypath->getNMember()*5+8;i<fdecaypath->getNMember()*5+14;i++) fitter.Config().ParSettings(i).Fix();//fix background parameters

    fitter.Config().SetMinimizer("Minuit2","Migrad");
    //fitter.Config().SetMinosErrors();

    fitter.FitFCN(fdecaypath->getNMember()*5+14,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);
    fitter.Result().Print(std::cout);
    const Double_t* resultparcenter=fitter.Result().GetParams();
    const Double_t* resulterrcenter=fitter.Result().GetErrors();

    Int_t nfreeparms = 0;
    for (Int_t i=0;i<fdecaypath->getNMember()*5+8;i++){
        pVal[i]=resultparcenter[i];
        if (!pvar[i]->isConstant()){
            pValError[i]=resulterrcenter[i];
            nfreeparms++;
        }
    }

    fitStatus=fitter.Result().Status();
    fitCovQual=nfreeparms+1;//fitter.Result().Ndf();
    fitNumInvalidNLL=fitter.Result().NCalls();
    fitEdm=fitter.Result().Edm();
    fitMinNll=fitter.Result().MinFcnValue();
    //foutputtree->Fill();

    //! construct parent/daugters decay components for plotting demonstration
    fB_parent=new TF1("fB_parent",totdecaymodel,&fitF::fcndecay_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay_parent");
    fB_daugter=new TF1("fB_daugter",totdecaymodel,&fitF::fcndecay_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay_daugter");

    fSB_parent=new TF1("fSB_parent",totdecaymodel,&fitF::fcndecay1n_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay1n_parent");
    fSB_daugter=new TF1("fSB_daugter",totdecaymodel,&fitF::fcndecay1n_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay1n_daugter");

    fSB2_parent=new TF1("fSB2_parent",totdecaymodel,&fitF::fcndecay2n_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay2n_parent");
    fSB2_daugter=new TF1("fSB2_daugter",totdecaymodel,&fitF::fcndecay2n_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay2n_daugter");

    fSB_c1=new TF1("fSB_c1",totdecaymodel,&fitF::fcndecay1n_c1,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay1n_c1");
    fSB_c2=new TF1("fSB_c2",totdecaymodel,&fitF::fcndecay1n_c2,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay1n_c2");
    fSB_c3=new TF1("fSB_c3",totdecaymodel,&fitF::fcndecay1n_c3,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay1n_c3");
    fSB_c23=new TF1("fSB_c23",totdecaymodel,&fitF::fcndecay1n_c23,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay1n_c23");
    fSB2_c1=new TF1("fSB2_c1",totdecaymodel,&fitF::fcndecay2n_c1,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay2n_c1");
    fSB2_c2=new TF1("fSB2_c2",totdecaymodel,&fitF::fcndecay2n_c2,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay2n_c2");
    fSB2_c3=new TF1("fSB2_c3",totdecaymodel,&fitF::fcndecay2n_c3,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay2n_c3");
    fSB2_c4=new TF1("fSB2_c4",totdecaymodel,&fitF::fcndecay2n_c4,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay2n_c4");
    fSB2_c134=new TF1("fSB2_c134",totdecaymodel,&fitF::fcndecay2n_c134,p_deadtime,p_timerange,fdecaypath->getNMember()*5+10,"fitF","fcndecay2n_c134");
    writeFitComponents();
    binw=hB->GetXaxis()->GetBinWidth(1);
    hB->Write();
    hSB->Write();
    hSB2->Write();
    plotResultsMore(1);

    char tempstr[500];
    sprintf(tempstr,"%s.txt",foutputData);
    std::ofstream ofs(tempstr,std::ios::app);
    for (int i=0;i<fdecaypath->getNMember()*5+8;i++){
        if (!pvar[i]->isConstant())
            ofs<<i<<"\t"<<pVal[i]<<"\t"<<pValError[i]<<std::endl;
    }
    ofs<<"nsig = "<<pVal[fdecaypath->getNMember()*5]<<"\tnbkg = "<<fdecaypath->getNMember()*5+8<<std::endl;
    ofs<<"chisquare/NDF = "<<chiSquareNDF<<"\t"<<chiSquareNDF1n<<"\t"<<chiSquareNDF2n<<std::endl;
    for (int i=0;i<fnMC;i++){
        generateMC();
        setValParameters();
        fitter.Config().SetParamsSettings(fdecaypath->getNMember()*5+14,binfitparms);
        fitter.FitFCN(fdecaypath->getNMember()*5+14,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);
        fitter.Result().Print(std::cout);
        const Double_t* resultpar=fitter.Result().GetParams();
        const Double_t* resulterr=fitter.Result().GetErrors();

        //const Double_t* resulterr=fitter.Result().GetErrors();
        Int_t nfreeparms = 0;
        for (Int_t i=0;i<fdecaypath->getNMember()*5+8;i++){
            pVal[i]=resultpar[i];
            if (!pvar[i]->isConstant()){
                pValError[i]=resulterr[i];
                nfreeparms++;
            }
        }
        fitStatus=fitter.Result().Status();
        fitCovQual=nfreeparms+1;//fitter.Result().Ndf();
        fitNumInvalidNLL=fitter.Result().NCalls();
        fitEdm=fitter.Result().Edm();
        fitMinNll=fitter.Result().MinFcnValue();
        calculateChiSquare(1);
        foutputtree->Fill();
    }
    writeOutputTree();
    closeOutputFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::writeFitComponents()
{
    for (int i=0;i<fdecaypath->getNMember()*5+10;i++){
        fB_parent->FixParameter(i,fB->GetParameter(i));
        fB_parentnobkg->FixParameter(i,fB->GetParameter(i));
        fB_daugter->FixParameter(i,fB->GetParameter(i));
        fB_daugternobkg->FixParameter(i,fB->GetParameter(i));
        fB_alpha->FixParameter(i,fB->GetParameter(i));
        fB_alphanobkg->FixParameter(i,fB->GetParameter(i));
        //        fSB_parent->FixParameter(i,fSB->GetParameter(i));
        //        fSB_daugter->FixParameter(i,fSB->GetParameter(i));
        //        fSB2_parent->FixParameter(i,fSB2->GetParameter(i));
        //        fSB2_daugter->FixParameter(i,fSB2->GetParameter(i));

        //        fSB_c1->FixParameter(i,fSB->GetParameter(i));
        //        fSB_c2->FixParameter(i,fSB->GetParameter(i));
        //        fSB_c3->FixParameter(i,fSB->GetParameter(i));
        //        fSB_c23->FixParameter(i,fSB->GetParameter(i));
        //        fSB2_c1->FixParameter(i,fSB2->GetParameter(i));
        //        fSB2_c2->FixParameter(i,fSB2->GetParameter(i));
        //        fSB2_c3->FixParameter(i,fSB2->GetParameter(i));
        //        fSB2_c4->FixParameter(i,fSB2->GetParameter(i));
        //        fSB2_c134->FixParameter(i,fSB2->GetParameter(i));
    }

    fB->SetNpx(nbinsHB*10);
    fB_bkgneg->SetNpx(nbinsHB*10);
    fB_bkgpos->SetNpx(nbinsHB*10);
    fB_parent->SetNpx(nbinsHB*10);
    fB_parentnobkg->SetNpx(nbinsHB*10);
    fB_daugter->SetNpx(nbinsHB*10);
    fB_daugternobkg->SetNpx(nbinsHB*10);
    fB_alpha->SetNpx(nbinsHB*10);
    fB_alphanobkg->SetNpx(nbinsHB*10);

    fA_all = new TF1("fA_all", [&](double *x, double *p) {
        return fB_alphanobkg->Eval(x[0]) + fA_bkgpos->Eval(x[0]);
    }, 0, p_timerange, 0);
    fA_all->SetNpx(nbinsHB*10);
    fA_all->SetLineWidth(3);
    fA_all->SetLineColor(7);

    fA_bkgneg->SetNpx(nbinsHB*10);
    fA_bkgneg->SetLineWidth(3);
    fA_bkgneg->SetLineColor(7);

    //    fSB->SetNpx(nbinsHSB*10);
    //    fSB_bkgneg->SetNpx(nbinsHSB*10);
    //    fSB_bkgpos->SetNpx(nbinsHSB*10);
    //    fSB_parent->SetNpx(nbinsHSB*10);
    //    fSB_daugter->SetNpx(nbinsHSB*10);
    //    fSB2->SetNpx(nbinsHSB2*10);
    //    fSB2_bkgneg->SetNpx(nbinsHSB2*10);
    //    fSB2_bkgpos->SetNpx(nbinsHSB2*10);
    //    fSB2_parent->SetNpx(nbinsHSB2*10);
    //    fSB2_daugter->SetNpx(nbinsHSB2*10);
    //    fSB_c1->SetNpx(nbinsHB*10);
    //    fSB_c2->SetNpx(nbinsHB*10);
    //    fSB_c3->SetNpx(nbinsHB*10);
    //    fSB_c23->SetNpx(nbinsHB*10);
    //    fSB2_c1->SetNpx(nbinsHSB2*10);
    //    fSB2_c2->SetNpx(nbinsHSB2*10);
    //    fSB2_c3->SetNpx(nbinsHSB2*10);
    //    fSB2_c4->SetNpx(nbinsHSB2*10);
    //    fSB2_c134->SetNpx(nbinsHSB2*10);

    fB->Write();
    //    fSB->Write();
    //    fSB2->Write();
    fB_bkgneg->Write();
    //    fSB_bkgneg->Write();
    //    fSB2_bkgneg->Write();
    fB_bkgpos->Write();
    //    fSB_bkgpos->Write();
    //    fSB2_bkgpos->Write();

    fB_parent->Write();
    fB_parentnobkg->Write();
    fB_daugter->Write();
    fB_daugternobkg->Write();
    fB_alpha->Write();
    fB_alphanobkg->Write();
    fA_all->Write();
    hA->Write();
    //    fSB_parent->Write();
    //    fSB_daugter->Write();
    //    fSB2_parent->Write();
    //    fSB2_daugter->Write();

    //    fSB_c1->Write();
    //    fSB_c2->Write();
    //    fSB_c3->Write();
    //    fSB_c23->Write();
    //    fSB2_c1->Write();
    //    fSB2_c2->Write();
    //    fSB2_c3->Write();
    //    fSB2_c4->Write();
    //    fSB2_c134->Write();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::calculateChiSquare(Int_t opt){
    Double_t chisquare = 0;
    Double_t chisquare1n = 0;
    Double_t chisquare2n = 0;
    //! 0n
    if (opt==0){
        for (Int_t i=0;i<model0nHist->GetN();i++){
            Double_t xi=model0nHist->GetX()[i];
            Double_t yi=model0nHist->GetY()[i];
            Double_t yeval=model0nCurve->Eval(xi);
            Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
            chisquare+=chisquarei;
        }
        chisquare=2*chisquare;
        if (fitres)
            chiSquareNDF=chisquare/(model0nHist->GetN()-fitres->floatParsFinal().getSize());
    }else{
        Int_t k=0;
        for (Int_t i=0;i<hB->GetNbinsX();i++){
            Double_t xi=hB->GetBinCenter(i+1);
            Double_t yi=hB->GetBinContent(i+1);
            if (xi>p_deadtime){
                Double_t yeval=fB->Eval(xi);
                Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
                chisquare+=chisquarei;
                k++;
            }
        }
        chisquare=2*chisquare;
        chiSquareNDF=chisquare/(k-fitCovQual);
    }

    //!1n
    if (opt==0){
        for (Int_t i=0;i<model1nHist->GetN();i++){
            Double_t xi=model1nHist->GetX()[i];
            Double_t yi=model1nHist->GetY()[i];
            Double_t yeval=model1nCurve->Eval(xi);
            Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
            chisquare1n+=chisquarei;
        }
        chisquare1n=2*chisquare1n;
        if (fitres)
            chiSquareNDF1n=chisquare1n/(model1nHist->GetN()-fitres->floatParsFinal().getSize());
    }else{
        Int_t k=0;
        for (Int_t i=0;i<hSB->GetNbinsX();i++){
            Double_t xi=hSB->GetBinCenter(i+1);
            Double_t yi=hSB->GetBinContent(i+1);
            if (xi>p_deadtime){
                Double_t yeval=fSB->Eval(xi);
                Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
                chisquare1n+=chisquarei;
                k++;
            }
        }
        chisquare1n=2*chisquare1n;
        chiSquareNDF1n=chisquare1n/(k-fitCovQual);
    }

    //!2n
    if (opt==0){
        for (Int_t i=0;i<model2nHist->GetN();i++){
            Double_t xi=model2nHist->GetX()[i];
            Double_t yi=model2nHist->GetY()[i];
            Double_t yeval=model2nCurve->Eval(xi);
            Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
            chisquare2n+=chisquarei;
        }
        chisquare2n=2*chisquare2n;
        if (fitres)
            chiSquareNDF2n=chisquare2n/(model2nHist->GetN()-fitres->floatParsFinal().getSize());
    }else{
        Int_t k=0;
        for (Int_t i=0;i<hSB2->GetNbinsX();i++){
            Double_t xi=hSB2->GetBinCenter(i+1);
            Double_t yi=hSB2->GetBinContent(i+1);
            if (xi>p_deadtime){
                Double_t yeval=fSB2->Eval(xi);
                Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
                chisquare2n+=chisquarei;
                k++;
            }
        }
        chisquare2n=2*chisquare2n;
        chiSquareNDF2n=chisquare2n/(k-fitCovQual);
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::calculateUpperLimit()
{


    //! Add 2023, Mar 1: integraged range reduced to only up to 20 times of half-life
    Double_t integratedTimeRange = TMath::Log(2)/fB_parent->GetParameter(0)*20;
    if (integratedTimeRange>p_timerange)
        integratedTimeRange = p_timerange;
    //! old
    //    Double_t integratedTimeRange = p_timerange;
    //! write P1n,P2n and P3n upper limit
    char tempstr[500];
    sprintf(tempstr,"%s.txt",foutputData);
    std::ofstream ofs(tempstr,std::ios::app);
    N0b=fB_parent->Eval(0.)/pvar[0]->getVal()/binw-fB_parent->Eval(integratedTimeRange)/pvar[0]->getVal()/binw;
    Double_t N0b2 = fB_parent->Eval(0.)/pvar[0]->getVal()/binw;
    Double_t Nparent= fB_parent->Integral(0.,50.);
    Double_t Ndaugter = fB_daugter->Integral(0.,50.);

    Double_t staterr_cnt =  fB_parent->Eval(0.)*nsig->getError()/nsig->getVal();
    Double_t dN0b = fB_parent->Eval(0.) / (fB->GetParameter(0) * binw)*sqrt(pow(staterr_cnt / fB_parent->Eval(0.), 2) + pow(pvar[0]->getError() / pvar[0]->getVal(), 2));
    Double_t dN0b_not12err = fB_parent->Eval(0.) / (fB->GetParameter(0) * binw)*staterr_cnt / fB_parent->Eval(0.);
    ofs<<N0b<<"\t"<<N0b2<<"\t"<<dN0b<<"\t"<<N0b*nsig->getError()/nsig->getVal()<<"\t"<<dN0b_not12err<<"\t"<<nsig->getVal()*Nparent/(Nparent+Ndaugter)<<endl;

    //    N0b1n=tree->Draw("",Form("x>%f&&x<%f&&y==1",p_deadtime,integratedTimeRange),"goff")-
    //            tree->Draw("",Form("x>%f&&x<%f&&y==1",-integratedTimeRange,-p_deadtime),"goff")-
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==1",p_deadtime,integratedTimeRange),"goff")+
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==1",-integratedTimeRange,-p_deadtime),"goff");
    //    N0b1n_bkg=tree->Draw("",Form("x>%f&&x<%f&&y==1",-integratedTimeRange,-p_deadtime),"goff")+
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==1",p_deadtime,integratedTimeRange),"goff")-
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==1",-integratedTimeRange,-p_deadtime),"goff");
    //    N0b2n=tree->Draw("",Form("x>%f&&x<%f&&y==2",p_deadtime,integratedTimeRange),"goff")-
    //            tree->Draw("",Form("x>%f&&x<%f&&y==2",-integratedTimeRange,-p_deadtime),"goff")-
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==2",p_deadtime,integratedTimeRange),"goff")+
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==2",-integratedTimeRange,-p_deadtime),"goff");
    //    N0b2n_bkg=tree->Draw("",Form("x>%f&&x<%f&&y==2",-integratedTimeRange,-p_deadtime),"goff")+
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==2",p_deadtime,integratedTimeRange),"goff")-
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==2",-integratedTimeRange,-p_deadtime),"goff");
    //    N0b3n=tree->Draw("",Form("x>%f&&x<%f&&y==3",p_deadtime,integratedTimeRange),"goff")-
    //            tree->Draw("",Form("x>%f&&x<%f&&y==3",-integratedTimeRange,-p_deadtime),"goff")-
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==3",p_deadtime,integratedTimeRange),"goff")+
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==3",-integratedTimeRange,-p_deadtime),"goff");
    //    N0b3n_bkg=tree->Draw("",Form("x>%f&&x<%f&&y==3",-integratedTimeRange,-p_deadtime),"goff")+
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==3",p_deadtime,integratedTimeRange),"goff")-
    //            treeb->Draw("",Form("x>%f&&x<%f&&y==3",-integratedTimeRange,-p_deadtime),"goff");
    //    ofs<<N0b<<"\t"<<N0b1n<<"\t"<<N0b2n<<"\t"<<N0b3n<<"\t"<<N0b1n_bkg<<"\t"<<N0b2n_bkg<<"\t"<<
    //         N0b3n_bkg<<"\t"<<pVal[fdecaypath->getNMember()*4]<<"\t"<<pValError[fdecaypath->getNMember()*4]<<"\t"<<
    //      pVal[fdecaypath->getNMember()*5+7]<<"\t"<<pValError[fdecaypath->getNMember()*5+7]<<"\t"<<integratedTimeRange<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::plotResultsMore(Int_t opt)
{

    calculateUpperLimit();

    gStyle->SetOptStat(0);
    TCanvas* c0n=new TCanvas("c0n","c0n",1200,800);
    TPad *pad1_c0n = new TPad("pad1_c0n","pad1_c0n",0,0.3,1,1);
    TPad *pad2_c0n = new TPad("pad2_c0n","pad2_c0n",0,0,1,0.3);
    pad1_c0n->SetTopMargin(0.09);
    pad1_c0n->SetBottomMargin(0.);
    pad1_c0n->SetBorderMode(0);
    //pad1_c0n->SetLogy();
    pad2_c0n->SetTopMargin(0.);
    pad2_c0n->SetBottomMargin(0.4);
    pad2_c0n->SetBorderMode(0);
    pad1_c0n->Draw();
    pad2_c0n->Draw();
    pad1_c0n->cd();
    Int_t npremove=0;
    if (opt==0){
        TH1F* hdummyc0n=new TH1F("hdummy0n","",20,plotrangelow,plotrangehi);
        hdummyc0n->Draw();
        hdummyc0n->GetYaxis()->SetRangeUser(0.,model0nHist->GetYaxis()->GetXmax());
        hdummyc0n->GetYaxis()->SetTitleSize(0.06);
        hdummyc0n->GetYaxis()->SetTitleOffset(0.58);
        hdummyc0n->GetYaxis()->SetTitle("Counts");
        hdummyc0n->GetXaxis()->SetTitle("t_{#beta} - t_{ion} (s)");
        hdummyc0n->GetYaxis()->SetLabelSize(0.05);
        hdummyc0n->GetXaxis()->SetTitleSize(0.06);
        hdummyc0n->GetXaxis()->SetLabelSize(0.05);

        //! for unbin fit
        for (Int_t i=0;i<model0nCurve->GetN();i++) if (model0nCurve->GetX()[i]<p_deadtime) npremove++; else break;
        for (Int_t i=0;i<npremove;i++) model0nCurve->RemovePoint(0);
        model0nHist->SetMarkerSize(1.2);
        modelbkg0nHist->SetMarkerSize(1.2);
        model0nHist->Draw("sameP");
        modelbkg0nHist->Draw("sameP");
        model0nCurve->SetLineWidth(3);
        model0nCurve->SetLineColor(4);
        model0nCurve->Draw("same");
    }else{
        hB->SetMarkerStyle(20);
        hB->SetMarkerSize(1.);
        hB->Draw("P0E");
        hB->GetXaxis()->SetRangeUser(plotrangelow,plotrangehi);
        hB->GetYaxis()->SetTitleSize(0.06);
        hB->GetYaxis()->SetTitleOffset(0.58);
        hB->GetYaxis()->SetTitle("Counts");
        hB->GetYaxis()->SetLabelSize(0.05);
        fB->SetLineWidth(3);
        fB->SetLineColor(4);
        fB->Draw("same");
    }

    fB_parent->SetLineWidth(3);
    fB_parent->SetLineColor(2);
    fB_parent->Draw("same");
    fB_daugter->SetLineWidth(3);
    fB_daugter->SetLineColor(6);
    fB_daugter->Draw("same");

    fB_alpha->SetLineWidth(3);
    fB_alpha->SetLineColor(7);
    fB_alpha->Draw("same");

    fB_alphanobkg->SetLineWidth(3);
    fB_alphanobkg->SetLineColor(7);
    //    fB_alphanobkg->Draw("same");

    fB_bkgneg->SetLineWidth(3);
    fB_bkgneg->SetLineColor(4);
    fB_bkgneg->Draw("same");
    pad1_c0n->Draw();

    pad2_c0n->cd();
    //!Calculate and residual plot - X2/ndf
    Double_t chisquare=0;
    Double_t xres[10000];
    Double_t yres[10000];
    Double_t yreserr[10000];
    Double_t yreserrlow[10000];

    TGraphAsymmErrors * resplot_0n;
    TGraphAsymmErrors * resplotbkg_0n;

    for (Int_t i=0;i<model0nHist->GetN();i++){
        Double_t xi=model0nHist->GetX()[i];
        if (xi<p_deadtime)
            continue;
        Double_t yi=model0nHist->GetY()[i];
        Double_t yeval=model0nCurve->Eval(xi);
        Double_t reldev=yeval-yi;
        xres[i]=xi;
        yres[i]=reldev;
        yreserr[i]=model0nHist->GetEYlow()[i];//sqrt(yi+yeval);
        yreserrlow[i]=model0nHist->GetEYhigh()[i];//sqrt(yi+yeval);
        Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
        chisquare+=chisquarei;
    }
    chisquare=2*chisquare;
    if (fitres){
        cout<<"ndf="<<fitres->floatParsFinal().getSize()<<endl;
        chiSquareNDF=chisquare/(model0nHist->GetN()-fitres->floatParsFinal().getSize());
    }
    cout<<"chisquare/ndf="<<chiSquareNDF<<endl;
    //        resplot_0n=new TGraphErrors(model0nHist->GetN(),xres,yres,0,yreserr);

    resplot_0n=new TGraphAsymmErrors(model0nHist->GetN(),xres,yres,0,0,yreserrlow,yreserr);
    for (Int_t i=0;i<modelbkg0nHist->GetN();i++){
        Double_t xi=modelbkg0nHist->GetX()[i];
        if (xi>-p_deadtime)
            continue;
        Double_t yi=modelbkg0nHist->GetY()[i];
        Double_t yeval=fB_bkgneg->Eval(xi);
        Double_t reldev=yeval-yi;
        xres[i]=xi;
        yres[i]=reldev;
        yreserr[i]=modelbkg0nHist->GetEYlow()[i];//sqrt(yi+yeval);
        yreserrlow[i]=modelbkg0nHist->GetEYhigh()[i];//sqrt(yi+yeval);
    }
    resplotbkg_0n=new TGraphAsymmErrors(modelbkg0nHist->GetN(),xres,yres,0,0,yreserrlow,yreserr);


    TH1F* hdummyc0nres=new TH1F("hdummyc0nres","",20,plotrangelow,plotrangehi);
    hdummyc0nres->SetLineColor(1);
    hdummyc0nres->SetLineWidth(3);
    hdummyc0nres->Draw();
    hdummyc0nres->GetYaxis()->SetRangeUser(resplot_0n->GetYaxis()->GetXmin(),resplot_0n->GetYaxis()->GetXmax());

    hdummyc0nres->GetYaxis()->SetLabelSize(0.12);
    hdummyc0nres->GetYaxis()->SetTitle("fit - data (counts)");
    hdummyc0nres->GetYaxis()->SetTitleSize(0.12);
    hdummyc0nres->GetYaxis()->SetTitleOffset(0.29);

    hdummyc0nres->GetXaxis()->SetLabelSize(0.14);
    hdummyc0nres->GetXaxis()->SetTitle("t_{#beta} - t_{ion} (s)");
    hdummyc0nres->GetXaxis()->SetTitleSize(0.17);
    hdummyc0nres->GetXaxis()->SetTitleOffset(0.95);

    resplot_0n->SetMarkerStyle(20);
    resplotbkg_0n->SetMarkerStyle(20);

    resplot_0n->SetLineColor(2);
    resplotbkg_0n->SetLineColor(2);
    resplot_0n->SetMarkerSize(1.2);
    resplotbkg_0n->SetMarkerSize(1.2);

    resplot_0n->Draw("sameP");
    resplotbkg_0n->Draw("sameP");
    pad2_c0n->Draw();
    c0n->Write();


    TCanvas* c1n=new TCanvas("c1n","c1n",1400,700);
    TPad *pad1_c1n = new TPad("pad1_c1n","pad1_c1n",0,0.3,1,1);
    TPad *pad2_c1n = new TPad("pad2_c1n","pad2_c1n",0,0,1,0.3);
    pad1_c1n->SetTopMargin(0.09);
    pad1_c1n->SetBottomMargin(0.00001);
    pad1_c1n->SetBorderMode(0);
    //pad1_c1n->SetLogy();
    pad2_c1n->SetTopMargin(0.00001);
    pad2_c1n->SetBottomMargin(0.4);
    pad2_c1n->SetBorderMode(0);
    pad1_c1n->Draw();
    pad2_c1n->Draw();
    pad1_c1n->cd();



    hA->SetMarkerStyle(20);
    hA->SetMarkerSize(1.);
    hA->SetMarkerColor(2);
    hA->Draw("P0E");

    hA->GetXaxis()->SetRangeUser(plotrangelow,plotrangehi);
    hA->GetYaxis()->SetTitleSize(0.06);
    hA->GetYaxis()->SetTitleOffset(0.58);
    hA->GetYaxis()->SetTitle("Counts");
    hA->GetYaxis()->SetLabelSize(0.05);
    fA_all->Draw("same");
    fA_bkgneg->Draw("same");


    pad1_c1n->Draw();
    pad2_c1n->cd();
    Double_t chisquare1n=0;

    TGraphErrors * resplot_1n;
    TGraphErrors * resplotbkg_1n;
    Int_t k=0;
    for (Int_t i=0;i<hA->GetNbinsX();i++){
        Double_t xi=hA->GetBinCenter(i+1);
        Double_t yi=hA->GetBinContent(i+1);
        if (xi>p_deadtime){
            Double_t yeval=fA_all->Eval(xi);
            Double_t reldev=yeval-yi;
            xres[k]=xi;
            yres[k]=reldev;
            yreserr[k]=sqrt(yi+yeval);
            Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
            chisquare1n+=chisquarei;
            k++;
        }
    }
    chisquare1n=2*chisquare1n;
    chiSquareNDF1n=chisquare1n/(k+fitCovQual);
    cout<<"chisquare/NDF 1n="<<chiSquareNDF1n<<endl;
    resplot_1n=new TGraphErrors(k,xres,yres,0,yreserr);
    k=0;
    for (Int_t i=0;i<hA->GetNbinsX();i++){
        Double_t xi=hA->GetBinCenter(i+1);
        Double_t yi=hA->GetBinContent(i+1);
        if (xi<0){
            Double_t yeval=fA_bkgneg->Eval(xi);
            Double_t reldev=yeval-yi;
            xres[k]=xi;
            yres[k]=reldev;
            yreserr[k]=sqrt(yi+yeval);
            k++;
        }
    }
    resplotbkg_1n=new TGraphErrors(k,xres,yres,0,yreserr);



    TH1F* hdummyc1nres=new TH1F("hdummyc1nres","",20,plotrangelow,plotrangehi);
    hdummyc1nres->SetLineColor(1);
    hdummyc1nres->SetLineWidth(2);
    hdummyc1nres->Draw();
    hdummyc1nres->GetYaxis()->SetRangeUser(resplot_1n->GetYaxis()->GetXmin(),resplot_1n->GetYaxis()->GetXmax());

    hdummyc1nres->GetYaxis()->SetLabelSize(0.12);
    hdummyc1nres->GetYaxis()->SetTitle("fit - data (counts)");
    hdummyc1nres->GetYaxis()->SetTitleSize(0.12);
    hdummyc1nres->GetYaxis()->SetTitleOffset(0.29);

    hdummyc1nres->GetXaxis()->SetLabelSize(0.14);
    hdummyc1nres->GetXaxis()->SetTitle("t_{#beta} - t_{ion} (s)");
    hdummyc1nres->GetXaxis()->SetTitleSize(0.17);
    hdummyc1nres->GetXaxis()->SetTitleOffset(0.95);

    resplot_1n->SetMarkerStyle(20);
    resplotbkg_1n->SetMarkerStyle(20);

    resplot_1n->SetLineColor(2);
    resplotbkg_1n->SetLineColor(2);

    resplot_1n->Draw("sameP");
    resplotbkg_1n->Draw("sameP");
    pad2_c1n->Draw();
    c1n->Write();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::Run()
{

    fitBackground();
    initFitParameters();
    //    printCurrentParameters();
    if (ffitopt==0 || ffitopt==2){
        setNormalFit();
    }else{
        setExernalContrainFit();
    }
    printCurrentParameters();
    prepareData();
    //    //prepareMonteCarloData(1);

    bookOutputTree();
    //! setup stuffs for MC fits
    setCentralParameters();
    getParameters();

    //! perform first fit
    doFit();
    plotResults();
    writeResults();
    writeResultsMC();

    for (int i=0;i<fnMC;i++){
        generateMC();
        setValParameters();
        printCurrentParameters();
        doFit();
        writeResultsMC();
        //        calculateChiSquare();
    }

    writeOutputTree();
    closeOutputFile();

    //    //Found that this ROOSTAT calculation change final value.
    //    double nsigval = nsig->getVal();
    //    const char* plcflag = std::getenv("PLC");
    //    if (plcflag){
    //        //Manual asymptotic Z-value estimation.
    //        double nll_sb = fitres->minNll(); // log(L_s+b)
    //        nsig->setVal(0);//
    //        nsig->setConstant(kTRUE);
    //        if (pvar[0]->isConstant()){
    //            pvar[0]->setMin(pvar[0]->getVal()-0.0001);
    //            pvar[0]->setMax(pvar[0]->getVal()+0.0001);
    //            pvar[0]->setConstant(kFALSE);
    //        }
    //        fitres=final_pdfT12->fitTo(*t12data,NumCPU(ncpu),Save(kTRUE),PrintLevel(3));
    //        double nll_b = fitres->minNll(); // log(L_b)
    //        double q0 = 2.0 * (nll_b - nll_sb);  // likelihood ratio
    //        if (q0 < 0) q0 = 0;
    //        double Z = sqrt(q0);
    //        double significance3 = Z;
    ////        nbkg->setConstant(false);
    ////        nbkg->setMin(nbkg->getVal()-nbkg->getError()*3);
    ////        nbkg->setMax(nbkg->getVal()+nbkg->getError()*3);
    //        nsig->setConstant(kFALSE);
    //        nsig->setVal(nsigval);
    //        ProfileLikelihoodCalculator plc(*t12data, *final_pdfT12, RooArgSet(*nsig,*nbkg));
    //        plc.SetConfidenceLevel(0.683);  // 1-sigma region
    //        LikelihoodInterval* interval = plc.GetInterval();
    //        upperLimit = interval->UpperLimit(*nsig);
    //        double lower = interval->LowerLimit(*nsig);
    //        // Rough-Rough estimation of Significance Relies indirectly on Wilks' theorem
    //        double significance1 = nsigval / ((upperLimit - lower) / 2.0);  // Approximate significance
    //        double significance2 = upperLimit / nsig->getError();  // Approximate significance 2

    //        //Buit in asymptotic Z-value estimation of ProfileLikelihoodCalculator
    //        RooArgSet nullparams("nullparams");
    //        nullparams.addClone(*nsig);
    //        nullparams.setRealValue(nsig->GetName(), 0);
    //        plc.SetNullParameters(nullparams);
    //        std::cout << "Perform Test of Hypothesis : null Hypothesis is " << nsig->GetName() << 0
    //                  << std::endl;
    //        auto result = plc.GetHypoTest();
    //        std::cout << "\n>>>> Hypotheis Test Result \n";
    //        result->Print();
    //        significance = result->Significance();
    //        cout<<significance<<"\t"<<significance3<<"\t"<<significance1<<"\t"<<significance2<<endl;
    //        char tempstr[500];
    //        sprintf(tempstr,"%s.txt",foutputData);
    //        std::ofstream ofs(tempstr,std::ios::app);
    //        ofs<<"Significane = "<<significance<<"\t"<<significance1<<"\t"<<significance2<<endl;

    //    }

    //    const char* MCflag = std::getenv("MC");
    //    if (MCflag){
    //        // Step 8: Create a workspace and ModelConfig
    //        RooWorkspace w("w", true);
    //        w.import(*final_pdfT12);
    //        w.import(*t12data);
    //        w.defineSet("obs", "x");
    //        w.defineSet("poi", "nsig");


    //        ModelConfig* mc = new ModelConfig("MyModel",&w);
    //        mc->SetPdf(*final_pdfT12);
    //        // Set observables
    ////        RooArgSet* obs = new RooArgSet(*x);
    //        mc->SetObservables(*x);
    //        // Set parameter of interest
    ////        RooArgSet* poi = new RooArgSet(*nsig);
    //        mc->SetParametersOfInterest(*nsig);
    //        // Set nuisance parameters (optional)
    ////        RooArgSet* nuis = new RooArgSet(*nbkg);
    //        mc->SetNuisanceParameters(*nbkg);
    //        mc->SetSnapshot(*nsig);
    //        ModelConfig* bModel = new ModelConfig(*mc);
    //        bModel->SetName("BackgroundOnly");
    //        nsig->setVal(0);
    //        bModel->SetSnapshot(*nsig);  // Force signal=0 for null hypothesis

    ////        ModelConfig* bModel = new ModelConfig("MyBModel",&w);
    ////        bModel->SetPdf(*final_pdfT12);
    ////        bModel->SetObservables(*x);
    ////        bModel->SetParametersOfInterest(*nsig);
    ////        mc->SetNuisanceParameters(*nbkg);
    ////        bModel->SetSnapshot(*nsig);
    ////        nsig->setVal(0);

    //        FrequentistCalculator fc(*t12data, *mc, *bModel);  // Model config for both signal+background and background-only

    //        fc.SetToys(1000,500);
    //        RooStats::HypoTestResult* result = fc.GetHypoTest();
    //        result->Print();
    //    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::generateRoofitEvaluate()
{
    path* fpath=new path;
    std::ifstream pathfile("path.txt");
    pathfile>>fpath->nri;
    pathfile>>fpath->npaths;
    for (int i=0;i<fpath->npaths;i++){
        pathfile>>fpath->ndecay[i];
        pathfile>>fpath->ispathhasflow[i];
        for (int j=0;j<fpath->ndecay[i];j++){
            pathfile>>fpath->decaymap[i][j]>>fpath->nneu[i][j];
        }
    }
    pathfile>>fpath->nisomers;
    for (int i=0;i<fpath->nisomers;i++){
        pathfile>>fpath->isomer_gs_index[i]>>fpath->isomer_ex_index[i];
    }
    pathfile.close();


    //! Generate functions for printing...
    std::ofstream ofnc("ofnc.txt");

    ofnc<<"#include <fitF.hh>"<<std::endl;
    ofnc<<"#include <RooAbsReal.h>"<<std::endl;
    ofnc<<"#include <RooAbsCategory.h>"<<std::endl;
    ofnc<<"#include <math.h>"<<std::endl;
    ofnc<<"#include <TMath.h>"<<std::endl;
    ofnc<<"#include <TStopwatch.h>"<<std::endl;
    ofnc<<"#ifdef EVAL_FAST"<<std::endl;
    ofnc<<"Double_t fitF::evaluate() const"<<std::endl;
    ofnc<<"{"<<std::endl;
    ofnc<<"double t = x;"<<std::endl;

    for (Int_t k=0;k<fpath->nri;k++){
        ofnc<<"double l"<<k<<"=(*p["<<k<<"]);"<<std::endl;
    }
    for (Int_t k=0;k<fpath->nri;k++){
        ofnc<<"double e"<<k<<"=exp(-l"<<k<<"*t);"<<std::endl;
    }
    for (Int_t k=0;k<fpath->nri;k++){
        ofnc<<"double p1n"<<k<<"=(*p["<<k+fpath->nri<<"]);"<<std::endl;
    }
    for (Int_t k=0;k<fpath->nri;k++){
        ofnc<<"double p2n"<<k<<"=(*p["<<k+fpath->nri*2<<"]);"<<std::endl;
    }
    for (Int_t k=0;k<fpath->nri;k++){
        ofnc<<"double py"<<k<<"=(*p["<<k+fpath->nri*3<<"]);"<<std::endl;
    }
    for (Int_t k=0;k<fpath->nri;k++){
        ofnc<<"double ne"<<k<<"=(*p["<<k+fpath->nri*4<<"]);"<<std::endl;
    }
#ifdef ISOMER_SUM_UNITY
    for (Int_t i=0;i<fpath->nisomers;i++){
        ofnc<<"py"<<fpath->isomer_ex_index[i]<<"=1-py"<<fpath->isomer_gs_index[i]<<";"<<std::endl;
    }
#endif
    ofnc<<"double be=*p["<<fpath->nri*5+4<<"];"<<std::endl;
    ofnc<<"double b1ne=*p["<<fpath->nri*5+5<<"];"<<std::endl;
    ofnc<<"double b2ne=*p["<<fpath->nri*5+6<<"];"<<std::endl;
    ofnc<<"double n1n2ne=*p["<<fpath->nri*5+7<<"];"<<std::endl;

    ofnc<<"double N0=*p["<<fpath->nri*5<<"]/l0;"<<std::endl;
    ofnc<<"double fparentdecay=l0*N0*e0;"<<std::endl;

    //! all decay function
    ofnc<<"double fdecay=fparentdecay*be;"<<std::endl;
    for (Int_t k=0;k<fpath->npaths;k++){
#ifdef PATHFLOW
        if (fpath->ispathhasflow[k]){
#endif
            // corefunction
            ofnc<<"double f"<<k<<"=l"<<fpath->decaymap[k][fpath->ndecay[k]-1]<<"*";
            for (int i=0;i<fpath->ndecay[k]-1;i++){
                if (fpath->nneu[k][i]==0){
                    ofnc<<"py"<<fpath->decaymap[k][i+1]<<"*(1-p1n"<<fpath->decaymap[k][i]<<"-p2n"<<fpath->decaymap[k][i]<<")*l"<<fpath->decaymap[k][i]<<"*";
                }else if (fpath->nneu[k][i]==1){
                    ofnc<<"py"<<fpath->decaymap[k][i+1]<<"*p1n"<<fpath->decaymap[k][i]<<"*l"<<fpath->decaymap[k][i]<<"*";
                }else{
                    ofnc<<"py"<<fpath->decaymap[k][i+1]<<"*p2n"<<fpath->decaymap[k][i]<<"*l"<<fpath->decaymap[k][i]<<"*";
                }
            }

            ofnc<<"(";
            for (int i=0;i<fpath->ndecay[k];i++){
                ofnc<<"e"<<fpath->decaymap[k][i]<<"/(";
                for (int j=0;j<fpath->ndecay[k];j++){
                    if (j!=i) {
                        ofnc<<"(l"<<fpath->decaymap[k][j]<<"-l"<<fpath->decaymap[k][i]<<")*";
                    }
                }
                ofnc<<"1)";
                ofnc<<"+";
            }
            ofnc<<"0)";
            ofnc<<"*N0;"<<std::endl;
            //end corefunction
#ifdef PATHFLOW
        }
#endif
    }
    ofnc<<"fdecay+=";
    for (Int_t k=0;k<fpath->npaths;k++){
#ifdef PATHFLOW
        if (fpath->ispathhasflow[k]){
#endif
            ofnc<<"f"<<k<<"+";
#ifdef PATHFLOW
        }
#endif
    }
    ofnc<<"0;"<<std::endl;
    ofnc<<std::endl;

    ofnc<<"double randcoinf2n=*p["<<fpath->nri*5+3<<"];"<<std::endl;
    ofnc<<"double randcoinfgt0n=*p["<<fpath->nri*5+2<<"];"<<std::endl;
    ofnc<<"double randcoinf1n=*p["<<fpath->nri*5+1<<"];"<<std::endl;


    //! calculation for 1neu
    ofnc<<"double fdecay1n=fparentdecay*(be*randcoinf1n+b1ne*ne0*p1n0*(1-randcoinf1n-randcoinfgt0n)+b2ne*2*(n1n2ne*(1-n1n2ne))*p2n0*(1-randcoinf1n-randcoinfgt0n)-b2ne*n1n2ne*n1n2ne*p2n0*randcoinf1n);"<<std::endl;

    for (Int_t i=0;i<fpath->npaths;i++){
#ifdef PATHFLOW
        if (fpath->ispathhasflow[i]){
#endif
            ofnc<<"fdecay1n+=f"<<i<<"*(randcoinf1n+ne"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*p1n"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*(1-randcoinf1n-randcoinfgt0n)+2*(ne"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*(1-ne"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"))*p2n"
               <<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*(1-randcoinf1n-randcoinfgt0n)-ne"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*ne"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*p2n"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*randcoinf1n);"<<std::endl;
#ifdef PATHFLOW
        }
#endif
    }

    //! calculation for 2neu
    ofnc<<"double fdecay2n=fparentdecay*(b2ne*n1n2ne*n1n2ne*p2n0*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n*be+b1ne*ne0*p1n0*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+b2ne*2*(n1n2ne*(1-n1n2ne))*p2n0*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));"<<std::endl;
    for (Int_t i=0;i<fpath->npaths;i++){
#ifdef PATHFLOW
        if (fpath->ispathhasflow[i]){
#endif
            ofnc<<"fdecay2n+=f"<<i<<"*(ne"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*ne"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*p2n"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*p1n"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*(1-ne"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"))*p2n"<<fpath->decaymap[i][fpath->ndecay[i]-1]<<"*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));"<<std::endl;
#ifdef PATHFLOW
        }
#endif
    }
    ofnc<<std::endl;
    ofnc<<"double ret; if (y==0) ret= fdecay-fdecay1n-fdecay2n; else if (y==1) ret= fdecay1n; else ret= fdecay2n;"<<std::endl;
    ofnc<<"return ret;"<<std::endl;
    ofnc<<"}"<<std::endl;
    ofnc<<"#endif"<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
