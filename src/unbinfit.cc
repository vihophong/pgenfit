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

    p_timerange_plus=20;
    nbinsHB=300;
    nbinsHSB=300;
    nbinsHSB2=300;

    ncpu=16;

    finputData=new char[1000];
    finputParms=new char[1000];

    fdecaypath=new decaypath;

    fStopWatch=new TStopwatch;
    fMCGenTime=0;
    fFitTime=0;

    //rseed=new TRandom3;
    rseed=new TRandom3(0);
    fnMC=0;

    fitStatus=-9999;
    fitCovQual=-9999;
    fitNumInvalidNLL=-9999;
    fitEdm=-9999;
    fitMinNll-9999;
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
    if (fnentrieslimit>0) tree->SetEntries(fnentrieslimit);

    //!*****************************************
    //! Setup histograms for unbinned fit
    //! *****************************************

    tree->Draw(Form("x>>hB(%d,%f,%f)",nbinsHB,-p_timerange,p_timerange_plus));
    tree->Draw(Form("x>>hSB(%d,%f,%f)",nbinsHSB,-p_timerange,p_timerange_plus),"y==1");
    tree->Draw(Form("x>>hSB2(%d,%f,%f)",nbinsHSB2,-p_timerange,p_timerange_plus),"y==2");
    char tempchar1[500];
    sprintf(tempchar1,"hB");
    hB=(TH1F*) gDirectory->Get(tempchar1);
    nsig_hB_firstbin=hB->GetBinContent(hB->FindBin(p_deadtime));
    sprintf(tempchar1,"hSB");
    hSB=(TH1F*) gDirectory->Get(tempchar1);
    sprintf(tempchar1,"hSB2");
    hSB2=(TH1F*) gDirectory->Get(tempchar1);

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
void unbinfit::fitBackground()
{
    // Fit and get positive backgrounds
    //!*****************************************
    //! Prepare and fit negative background function
    //! *****************************************
    // bkg parameters
    bkg1nratio=new RooRealVar("bkg1nratio","bkg1nratio",0.5,0,1) ;
    bkg2nratio=new RooRealVar("bkg2nratio","bkg2nratio",0.5,0,1) ;

    RooRealVar slope1("slope1","slope1",0.,-0.1,0.1) ;
    RooRealVar slope2("slope2","slope2",0.,-0.1,0.1) ;
    RooRealVar slope3("slope3","slope3",0.,-0.1,0.1) ;

    // bkg pdf
    fitFbkg bkgmodel("bkgmodel","bkgmodel",*xbkg,*y,*bkg1nratio,*bkg2nratio,slope1,slope2,slope3);
    // fit background
    bkgmodel.fitTo(*databkg,NumCPU(ncpu),Save()) ;

    //!*****************************************
    //! Prepare positive background function
    //! *****************************************
    // set positive background parms
#ifdef FLAT_BACKGROUNDS
    slope1pos=new RooRealVar("slope1pos","slope1pos",0,-0.1,0.1) ;
    slope2pos=new RooRealVar("slope2pos","slope2pos",0,-0.1,0.1) ;
    slope3pos=new RooRealVar("slope3pos","slope3pos",0,-0.1,0.1) ;
#else
    slope1pos=new RooRealVar("slope1pos","slope1pos",-slope1.getVal(),-0.1,0.1) ;
    slope2pos=new RooRealVar("slope2pos","slope2pos",-slope2.getVal(),-0.1,0.1) ;
    slope3pos=new RooRealVar("slope3pos","slope3pos",-slope3.getVal(),-0.1,0.1) ;
#endif


    slope1pos->setError(slope1.getError());
    slope2pos->setError(slope2.getError());
    slope3pos->setError(slope3.getError());

    // bkg pdf positive
    bkgmodelpos=new fitFbkg("bkgmodelpos","bkgmodelpos",*x,*y,*bkg1nratio,*bkg2nratio,*slope1pos,*slope2pos,*slope3pos);


    //! bin fit background negative
    fB_bkgneg=new TF1("fB_bkgneg","pol1",-p_timerange,0);
    fSB_bkgneg=new TF1("fSB_bkgneg","pol1",-p_timerange,0);
    fSB2_bkgneg=new TF1("fSB2_bkgneg","pol1",-p_timerange,0);
    //fB_bkgneg->FixParameter(1,slope3.getVal());
    //fSB_bkgneg->FixParameter(1,slope1.getVal());
    //fSB2_bkgneg->FixParameter(1,slope2.getVal());
    hB->Fit(fB_bkgneg,"LEQR0+","goff");
    hSB->Fit(fSB_bkgneg,"LEQR0+","goff");
    hSB2->Fit(fSB2_bkgneg,"LEQR0+","goff");

    fB_bkgpos=new TF1("fB_bkgpos","pol1",p_deadtime,p_timerange_plus);
    fSB_bkgpos=new TF1("fSB_bkgpos","pol1",p_deadtime,p_timerange_plus);
    fSB2_bkgpos=new TF1("fSB2_bkgpos","pol1",p_deadtime,p_timerange_plus);
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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

        p[fdecaypath->getNMember()+i]=new RooRealVar(Form("p%d",fdecaypath->getNMember()+i),Form("p%d",fdecaypath->getNMember()+i),fdecaypath->getMember(i)->decay_p1n,fdecaypath->getMember(i)->decay_p1nlow,fdecaypath->getMember(i)->decay_p1nup);
        pvar[fdecaypath->getNMember()+i]=(RooRealVar*) p[fdecaypath->getNMember()+i];
        pvar[fdecaypath->getNMember()+i]->setError(fdecaypath->getMember(i)->decay_p1nerr);

        p[fdecaypath->getNMember()*2+i]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*2+i),Form("p%d",fdecaypath->getNMember()*2+i),fdecaypath->getMember(i)->decay_p2n,fdecaypath->getMember(i)->decay_p2nlow,fdecaypath->getMember(i)->decay_p2nup);
        pvar[fdecaypath->getNMember()*2+i]=(RooRealVar*) p[fdecaypath->getNMember()*2+i];
        pvar[fdecaypath->getNMember()*2+i]->setError(fdecaypath->getMember(i)->decay_p1nerr);

        p[fdecaypath->getNMember()*3+i]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*3+i),Form("p%d",fdecaypath->getNMember()*3+i),fdecaypath->getMember(i)->population_ratio,fdecaypath->getMember(i)->population_ratiolow,fdecaypath->getMember(i)->population_ratioup);
        pvar[fdecaypath->getNMember()*3+i]=(RooRealVar*) p[fdecaypath->getNMember()*3+i];
        pvar[fdecaypath->getNMember()*3+i]->setError(fdecaypath->getMember(i)->population_ratioerr);

        p[fdecaypath->getNMember()*4+i]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*4+i),Form("p%d",fdecaypath->getNMember()*4+i),fdecaypath->getMember(i)->neueff,fdecaypath->getMember(i)->neuefflow,fdecaypath->getMember(i)->neueffup);
        pvar[fdecaypath->getNMember()*4+i]=(RooRealVar*) p[fdecaypath->getNMember()*4+i];
        pvar[fdecaypath->getNMember()*4+i]->setError(fdecaypath->getMember(i)->neuefferr);
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

    // set background/signal counts
    Double_t nnsig=tree->Draw("",Form("x>%f&&x<%f",p_deadtime,p_timerange),"goff");
    Double_t nnbkg=tree->Draw("",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    nnsig=nnsig-nnbkg;

    nbkg=new RooRealVar("nbkg","nbkg",nnbkg,nnbkg/3,nnbkg*3);
    nsig=new RooRealVar("nsig","nsig",nnsig,0,nnsig*10);
    nbkg->setError(TMath::Sqrt(nnbkg));
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
        if (fdecaypath->getMember(i)->is_population_ratio_fix!=0)
            pvar[fdecaypath->getNMember()*3+i]->setConstant(kTRUE);
        if (fdecaypath->getMember(i)->is_neueff_fix!=0)
            pvar[fdecaypath->getNMember()*4+i]->setConstant(kTRUE);
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
        if (fdecaypath->getMember(i)->is_neueff_fix==2){
            //cout<<"EEEEE "<<fdecaypath->getMember(i)->neueff<<endl;
            pvar[fdecaypath->getNMember()*4+i]->setConstant(kTRUE);
        }
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
    totdecaymodel=new fitF("totdecaymodel","totdecaymodel",*x,*y,p);
    final_pdf=new RooAddPdf("final_pdf","final pdf",RooArgList(*totdecaymodel,*bkgmodelpos),RooArgList(*nsig,*nbkg));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::prepareData()
{
    //!******************************************
    //! Prepare data set for fitting forward correlated data
    //! *****************************************
    data=new RooDataSet("data","data",RooArgSet(*x,*y),Import(*tree)) ;
    data->Print() ;
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
    for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
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
    for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
        pCentralVal[i]=pvar[i]->getVal();
        ipVal[i]=i;
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
    for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
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
    for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
        pvar[i]->setVal(pVal[i]);
    }
    nbkg->setVal(nbkgVal);
    bkg1nratio->setVal(bkg1nratioVal);
    bkg2nratio->setVal(bkg2nratioVal);
    slope1pos->setVal(slope1posVal);
    slope2pos->setVal(slope2posVal);
    slope3pos->setVal(slope3posVal);

    //! for binfit
    for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
        binfitparms[i]=pVal[i];
    }
    binfitparms[fdecaypath->getNMember()*5]=nsig_hB_firstbin;//initial activity
    binfitparms[fdecaypath->getNMember()*5+4]=binfitbkgparms[0];//bkg hdecay offset
    binfitparms[fdecaypath->getNMember()*5+5]=binfitbkgparms[1];//bkg hdecay slope
    binfitparms[fdecaypath->getNMember()*5+6]=binfitbkgparms[2];//bkg hdecay1n offset
    binfitparms[fdecaypath->getNMember()*5+7]=binfitbkgparms[3];//bkg hdecay1n slope
    binfitparms[fdecaypath->getNMember()*5+8]=binfitbkgparms[4];//bkg hdecay2n offset
    binfitparms[fdecaypath->getNMember()*5+9]=binfitbkgparms[5];//bkg hdecay2n slope
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::generateMC()
{
    for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
        if ((pvar[i]->isConstant())&&pValError[i]!=0) pVal[i]=rseed->Gaus(pCentralVal[i],pValError[i]);
        else pVal[i]=pCentralVal[i];

#ifdef PARENT_NEUEFF_UNIFORM
        //! neueff of parent randomly distributed from 40 to 68%
        if (i==fdecaypath->getNMember()*4) pVal[i]=0.4+rseed->Rndm()*(0.68-0.4);
#endif
    }
    nsigVal=nsigCentralVal;//reset central parameters after 1 fit

    nbkgVal=rseed->Gaus(nbkgCentralVal,nbkgError);
    bkg1nratioVal=rseed->Gaus(bkg1nratioCentralVal,bkg1nratioError);
    bkg2nratioVal=rseed->Gaus(bkg2nratioCentralVal,bkg2nratioError);

    slope1posVal=slope1posCentralVal;
    slope2posVal=slope2posCentralVal;
    slope3posVal=slope3posCentralVal;
    //slope1posVal=rseed->Gaus(slope1posCentralVal,slope1posError);
    //slope2posVal=rseed->Gaus(slope2posCentralVal,slope2posError);
    //slope3posVal=rseed->Gaus(slope3posCentralVal,slope3posError);


    //! for binfit
    //no variation background for now

    //binfitbkgparms[0]=rseed->Gaus(fB_bkgpos->GetParameter(0),fB_bkgneg->GetParError(0));
    //binfitbkgparms[1]=rseed->Gaus(fB_bkgpos->GetParameter(1),fB_bkgneg->GetParError(1));
    //binfitbkgparms[2]=rseed->Gaus(fSB_bkgpos->GetParameter(0),fB_bkgneg->GetParError(0));
    //binfitbkgparms[3]=rseed->Gaus(fSB_bkgpos->GetParameter(1),fB_bkgneg->GetParError(1));
    //binfitbkgparms[4]=rseed->Gaus(fSB2_bkgpos->GetParameter(0),fB_bkgneg->GetParError(0));
    //binfitbkgparms[5]=rseed->Gaus(fSB2_bkgpos->GetParameter(1),fB_bkgneg->GetParError(1));

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::doFit()
{
    //!******************************************
    //! Perform the fit
    //! *****************************************
    fStopWatch->Clear();
    fStopWatch->Start();

    if (ffitopt==0)
        fitres=final_pdf->fitTo(*data,NumCPU(ncpu),Save(kTRUE),PrintLevel(3));
    else
        fitres=final_pdf->fitTo(*data,ExternalConstraints(*externalconstrains),NumCPU(ncpu),Save(kTRUE),PrintLevel(3));
    fStopWatch->Stop();
    fFitTime=fStopWatch->RealTime();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::writeResultsMC()
{
    getParameters();
    fitStatus=fitres->status();
    fitCovQual=fitres->covQual();
    fitNumInvalidNLL=fitres->numInvalidNLL();
    fitEdm=fitres->edm();
    fitMinNll=fitres->minNll();
    foutputtree->Fill();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::plotResults()
{
    //!******************************************
    //! Outputs
    //! *****************************************
    TCanvas *c1 = new TCanvas("c1","c1",1200, 800);
    c1->Divide(3,1);
    c1->cd(1);
    RooPlot* xframe0 = x->frame(Title("0 neutron fit")) ;
    data->plotOn(xframe0,Cut("y==y::0neu"),Binning(500)) ;
    final_pdf->plotOn(xframe0,Slice(*y,"0neu")) ;
    //RooPlot* xframe0 = x->frame(Title("all decay fit")) ;
    //data->plotOn(xframe0,Binning(500)) ;
    //final_pdf->plotOn(xframe0) ;
    xframe0->Draw() ;
    c1->cd(2);
    RooPlot* xframe1 = x->frame(Title("1 neutron fit")) ;
    data->plotOn(xframe1,Cut("y==y::1neu"),Binning(500)) ;
    final_pdf->plotOn(xframe1,Slice(*y,"1neu")) ;
    xframe1->Draw() ;
    c1->cd(3);
    RooPlot* xframe2 = x->frame(Title("2 neutron fit")) ;
    data->plotOn(xframe2,Cut("y==y::2neu"),Binning(500)) ;
    final_pdf->plotOn(xframe2,Slice(*y,"2neu")) ;
    xframe2->Draw() ;
    c1->Write();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void unbinfit::writeResults()
{
    char tempstr[500];
    sprintf(tempstr,"%s.txt",foutputData);
    std::ofstream ofs(tempstr,std::ios::app);
    for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
        if (!pvar[i]->isConstant())
            ofs<<i<<"\t"<<pvar[i]->getVal()<<"\t"<<pvar[i]->getError()<<std::endl;
    }
    fitres->Print();
    std::cout<<"Time for MC generation = "<<fMCGenTime<<std::endl;
    std::cout<<"Time for Fitting = "<<fFitTime<<std::endl;
}


void unbinfit::RunBinFit()
{
    fitBackground();

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
   fB=new TF1("fB",totdecaymodel,&fitF::fcndecay,p_deadtime,p_timerange_plus,fdecaypath->getNMember()*5+6,"fitF","fcndecay");
   for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
       fB->FixParameter(i,pCentralVal[i]);
   }
   fB->FixParameter(fdecaypath->getNMember()*5+4,100);
   fB->FixParameter(fdecaypath->getNMember()*5+5,0);


   fSB=new TF1("fSB",totdecaymodel,&fitF::fcndecay1n,p_deadtime,p_timerange_plus,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n");
   for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
       fSB->FixParameter(i,pCentralVal[i]);
   }
   fSB->FixParameter(fdecaypath->getNMember()*5+4,10);
   fSB->FixParameter(fdecaypath->getNMember()*5+5,0);

   fSB2=new TF1("fSB2",totdecaymodel,&fitF::fcndecay2n,p_deadtime,p_timerange_plus,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n");
   for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
       fSB2->FixParameter(i,pCentralVal[i]);
   }
   fSB2->FixParameter(fdecaypath->getNMember()*5+4,1);
   fSB2->FixParameter(fdecaypath->getNMember()*5+5,0);

   //! set fix parameters
   fB->FixParameter(fdecaypath->getNMember()*5,nsig_hB_firstbin);
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
   fitter.Config().SetParamsSettings(fdecaypath->getNMember()*5+10,binfitparms);
   for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
       if (pvar[i]->isConstant()&&i!=fdecaypath->getNMember()*5){
           fitter.Config().ParSettings(i).Fix();
       }else{
           if (i==fdecaypath->getNMember()*5)
               fitter.Config().ParSettings(i).SetLimits(0,nsig_hB_firstbin*20);
           else
               fitter.Config().ParSettings(i).SetLimits(pvar[i]->getMin(),pvar[i]->getMax());
       }
   }
   for (int i=fdecaypath->getNMember()*5+4;i<fdecaypath->getNMember()*5+10;i++) fitter.Config().ParSettings(i).Fix();//fix background parameters


   fitter.Config().SetMinimizer("Minuit2","Migrad");
   //fitter.Config().SetMinosErrors();

   fitter.FitFCN(fdecaypath->getNMember()*5+10,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);
   fitter.Result().Print(std::cout);
   const Double_t* resultparcenter=fitter.Result().GetParams();
   const Double_t* resulterrcenter=fitter.Result().GetErrors();
   for (Int_t i=0;i<fdecaypath->getNMember()*5+4;i++){
       pVal[i]=resultparcenter[i];
       if (!pvar[i]->isConstant())
        pValError[i]=resulterrcenter[i];
   }
   fitStatus=fitter.Result().Status();
   fitCovQual=fitter.Result().Ndf();
   fitNumInvalidNLL=fitter.Result().NCalls();
   fitEdm=fitter.Result().Edm();
   fitMinNll=fitter.Result().MinFcnValue();
   foutputtree->Fill();

   char tempstr[500];
   sprintf(tempstr,"%s.txt",foutputData);
   std::ofstream ofs(tempstr,std::ios::app);
   for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
       if (!pvar[i]->isConstant())
           ofs<<i<<"\t"<<pVal[i]<<"\t"<<pValError[i]<<std::endl;
   }


   for (int i=0;i<fnMC;i++){
       generateMC();
       setValParameters();
       fitter.Config().SetParamsSettings(fdecaypath->getNMember()*5+10,binfitparms);
       fitter.FitFCN(fdecaypath->getNMember()*5+10,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);
       fitter.Result().Print(std::cout);
       const Double_t* resultpar=fitter.Result().GetParams();
       const Double_t* resulterr=fitter.Result().GetErrors();
       for (Int_t i=0;i<fdecaypath->getNMember()*5+4;i++){
           pVal[i]=resultpar[i];
           if (!pvar[i]->isConstant())
            pValError[i]=resulterrcenter[i];
       }
       fitStatus=fitter.Result().Status();
       fitCovQual=fitter.Result().Ndf();
       fitNumInvalidNLL=fitter.Result().NCalls();
       fitEdm=fitter.Result().Edm();
       fitMinNll=fitter.Result().MinFcnValue();
       foutputtree->Fill();
   }

   hB->Write();
   hSB->Write();
   hSB2->Write();
   fB->Write();
   fSB->Write();
   fSB2->Write();
   fB_bkgneg->Write();
   fSB_bkgneg->Write();
   fSB2_bkgneg->Write();
   fB_bkgpos->Write();
   fSB_bkgpos->Write();
   fSB2_bkgpos->Write();

   writeOutputTree();
   closeOutputFile();
}


void unbinfit::Run()
{
    fitBackground();

    initFitParameters();
    if (ffitopt==0)
        setNormalFit();
    else
        setExernalContrainFit();
    printCurrentParameters();
    prepareData();
    //prepareMonteCarloData(1);


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
    }
    writeOutputTree();
    closeOutputFile();

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

    ofnc<<"double N0=*p["<<fpath->nri*5<<"]/l0;"<<std::endl;
    ofnc<<"double fparentdecay=l0*N0*e0;"<<std::endl;

    //! all decay function
    ofnc<<"double fdecay=fparentdecay;"<<std::endl;
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
    //! random coinc of beta decay of parent
    //ofnc<<"double fdecay1n=fparentdecay*randcoinf1n;"<<std::endl;
    //! decay with 1 neutron of parent
    //ofnc<<"fdecay1n+=ne0*p1n0*fparentdecay*(1-randcoinf1n-randcoinfgt0n);"<<std::endl;
    //! decay with 1 neutron of parent from p2n
    //ofnc<<"fdecay1n+=2*(ne0*(1-ne0))*p2n0*fparentdecay*(1-randcoinf1n-randcoinfgt0n);"<<std::endl;
    //! decay with 2 neutron of parent (not random 1 neutron)
    //ofnc<<"fdecay1n-=ne0*ne0*p2n0*fparentdecay*randcoinf1n;"<<std::endl;

    ofnc<<"double fdecay1n=fparentdecay*(randcoinf1n+ne0*p1n0*(1-randcoinf1n-randcoinfgt0n)+2*(ne0*(1-ne0))*p2n0*(1-randcoinf1n-randcoinfgt0n)-ne0*ne0*p2n0*randcoinf1n);"<<std::endl;
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
    //! decay with 2 neutron from P2n of parent
    //ofnc<<"double fdecay2n=ne0*ne0*p2n0*fparentdecay*(1-randcoinf2n-randcoinfgt0n);"<<std::endl;
    //! random coinc of beta decay of parent
    //ofnc<<"fdecay2n+=fparentdecay*randcoinf2n;"<<std::endl;
    //! random 1n decay of parent
    //ofnc<<"fdecay2n+=ne0*p1n0*fparentdecay*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);"<<std::endl;
    //! decay with 1 neutron from P2n of parent - randomly correlated
    //ofnc<<"fdecay2n+=2*(ne0*(1-ne0))*p2n0*fparentdecay*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n);"<<std::endl;
    ofnc<<"double fdecay2n=fparentdecay*(ne0*ne0*p2n0*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne0*p1n0*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne0*(1-ne0))*p2n0*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));"<<std::endl;
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
