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
    nbinsHB=200;
    nbinsHSB=200;
    nbinsHSB2=200;

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

    plotrangelow=-0.5;
    plotrangehi=3;

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
    bkgmodelneg= new fitFbkg("bkgmodel","bkgmodel",*xbkg,*y,*bkg1nratio,*bkg2nratio,slope1,slope2,slope3);
    // fit background
    bkgmodelneg->fitTo(*databkg,NumCPU(ncpu),Save()) ;

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
    for (int i=0;i<fdecaypath->getNMember()*5+4;i++){//decay parameters
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
    c1->Divide(2,3);
    c1->cd(1);
    RooPlot* xframe0 = x->frame(Title("all fit")) ;
    data->plotOn(xframe0,Binning(nbinsHB/2),RooFit::Name("data0n")) ;
    final_pdf->plotOn(xframe0,RooFit::Name("data0nmodel")) ;
    //bkgmodelpos->plotOn(xframe0,RooFit::Name("bkg0nposmodel")) ;
    xframe0->Draw() ;
    c1->cd(3);
    RooPlot* xframe1 = x->frame(Title("1 neutron fit")) ;
    data->plotOn(xframe1,Cut("y==y::1neu"),Binning(nbinsHSB/2),RooFit::Name("data1n")) ;
    final_pdf->plotOn(xframe1,Slice(*y,"1neu"),RooFit::Name("data1nmodel")) ;
    //bkgmodelpos->plotOn(xframe1,Slice(*y,"1neu"),RooFit::Name("bkg1nposmodel")) ;
    xframe1->Draw() ;
    c1->cd(5);
    RooPlot* xframe2 = x->frame(Title("2 neutron fit")) ;
    data->plotOn(xframe2,Cut("y==y::2neu"),Binning(nbinsHSB2/2),RooFit::Name("data2n")) ;
    final_pdf->plotOn(xframe2,Slice(*y,"2neu"),RooFit::Name("data2nmodel")) ;
    //bkgmodelpos->plotOn(xframe2,Slice(*y,"2neu"),RooFit::Name("bkg2nposmodel")) ;
    xframe2->Draw() ;

    c1->cd(2);
    RooPlot* xframe3 = xbkg->frame(Title("all fit bkg")) ;
    databkg->plotOn(xframe3,Binning(nbinsHB/2),RooFit::Name("bkg0n")) ;
    xframe3->Draw() ;

    c1->cd(4);
    RooPlot* xframe4 = xbkg->frame(Title("1 neutron fit bkg")) ;
    databkg->plotOn(xframe4,Cut("y==y::1neu"),Binning(nbinsHSB/2),RooFit::Name("bkg1n")) ;
    xframe4->Draw() ;

    c1->cd(6);
    RooPlot* xframe5 = xbkg->frame(Title("2 neutron fit bkg")) ;
    databkg->plotOn(xframe5,Cut("y==y::2neu"),Binning(nbinsHSB2/2),RooFit::Name("bkg2n")) ;
    xframe5->Draw() ;
    c1->Write();

    model0nCurve=(RooCurve*)xframe0->getCurve("data0nmodel");
    model0nHist=(RooHist*)xframe0->getHist("data0n");
    model1nCurve=(RooCurve*)xframe1->getCurve("data1nmodel");
    model1nHist=(RooHist*)xframe1->getHist("data1n");
    model2nCurve=(RooCurve*)xframe2->getCurve("data2nmodel");
    model2nHist=(RooHist*)xframe2->getHist("data2n");

    modelbkg0nHist=(RooHist*)xframe3->getHist("bkg0n");
    modelbkg1nHist=(RooHist*)xframe4->getHist("bkg1n");
    modelbkg2nHist=(RooHist*)xframe5->getHist("bkg2n");

    model0nHist->Write();
    model1nHist->Write();
    model2nHist->Write();
    model0nCurve->Write();
    model1nCurve->Write();
    model2nCurve->Write();

    modelbkg0nHist->Write();
    modelbkg1nHist->Write();
    modelbkg2nHist->Write();

    //!******************************************
    //! Plotting components using bin fit
    //! *****************************************
    //! fix bkg parms for binned model from result of backward fit
    Double_t bkg1nratioval=bkg1nratio->getVal();
    Double_t bkg2nratioval=bkg2nratio->getVal();
    Double_t slope1posval=slope1pos->getVal();
    Double_t slope2posval=slope2pos->getVal();
    Double_t slope3posval=slope3pos->getVal();
    Double_t a0=slope3posval;
    Double_t b0=(2*nbkg->getVal()-a0*p_timerange*p_timerange)/2/p_timerange*(p_timerange/nbinsHB*2);
    Double_t a1=slope1posval*bkg1nratioval;
    Double_t b1=(2*nbkg->getVal()*bkg1nratioval-a1*p_timerange*p_timerange)/2/p_timerange*(p_timerange/nbinsHB*2);
    Double_t a2=slope2posval*bkg1nratioval*bkg2nratioval;
    Double_t b2=(2*nbkg->getVal()*bkg1nratioval*bkg2nratioval-a2*p_timerange*p_timerange)/2/p_timerange*(p_timerange/nbinsHB*2);

    fB_bkgpos->FixParameter(0,b0);
    fSB_bkgpos->FixParameter(0,b1);
    fSB2_bkgpos->FixParameter(0,b2);
    fB_bkgpos->FixParameter(1,a0);
    fSB_bkgpos->FixParameter(1,a1);
    fSB2_bkgpos->FixParameter(1,a2);

    fB_bkgneg->FixParameter(0,b0);
    fSB_bkgneg->FixParameter(0,b1);
    fSB2_bkgneg->FixParameter(0,b2);
    fB_bkgneg->FixParameter(1,-a0);
    fB_bkgneg->FixParameter(1,-a1);
    fB_bkgneg->FixParameter(1,-a2);

    totdecaymodelforplot=new fitF("totdecaymodelforplot","totdecaymodelforplot",*x,*y,p);
    totdecaymodelforplot->initPath();
    fB=new TF1("fB",totdecaymodelforplot,&fitF::fcndecay,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay");
    for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
        fB->FixParameter(i,pvar[i]->getVal());
    }
    fB->FixParameter(fdecaypath->getNMember()*5+4,b0);
    fB->FixParameter(fdecaypath->getNMember()*5+5,a0);


    fSB=new TF1("fSB",totdecaymodelforplot,&fitF::fcndecay1n,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n");
    for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
        fSB->FixParameter(i,pvar[i]->getVal());
    }
    fSB->FixParameter(fdecaypath->getNMember()*5+4,b1);
    fSB->FixParameter(fdecaypath->getNMember()*5+5,a1);

    fSB2=new TF1("fSB2",totdecaymodelforplot,&fitF::fcndecay2n,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n");
    for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
        fSB2->FixParameter(i,pvar[i]->getVal());
    }
    fSB2->FixParameter(fdecaypath->getNMember()*5+4,b2);
    fSB2->FixParameter(fdecaypath->getNMember()*5+5,a2);

    //! set fix parameters
    nsig_hB_firstbin=model0nCurve->Eval(p_deadtime)-fB_bkgpos->Eval(p_deadtime);
    fB->FixParameter(fdecaypath->getNMember()*5,nsig_hB_firstbin);
    fSB->FixParameter(fdecaypath->getNMember()*5,nsig_hB_firstbin);
    fSB2->FixParameter(fdecaypath->getNMember()*5,nsig_hB_firstbin);

    //! construct parent/daugters decay components for plotting demonstration
    fB_parent=new TF1("fB_parent",totdecaymodelforplot,&fitF::fcndecay_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay_parent");
    fB_daugter=new TF1("fB_daugter",totdecaymodelforplot,&fitF::fcndecay_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay_daugter");

    fSB_parent=new TF1("fSB_parent",totdecaymodelforplot,&fitF::fcndecay1n_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_parent");
    fSB_daugter=new TF1("fSB_daugter",totdecaymodelforplot,&fitF::fcndecay1n_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_daugter");

    fSB2_parent=new TF1("fSB2_parent",totdecaymodelforplot,&fitF::fcndecay2n_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_parent");
    fSB2_daugter=new TF1("fSB2_daugter",totdecaymodelforplot,&fitF::fcndecay2n_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_daugter");

    fSB_c1=new TF1("fSB_c1",totdecaymodelforplot,&fitF::fcndecay1n_c1,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_c1");
    fSB_c2=new TF1("fSB_c2",totdecaymodelforplot,&fitF::fcndecay1n_c2,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_c2");
    fSB_c3=new TF1("fSB_c3",totdecaymodelforplot,&fitF::fcndecay1n_c3,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_c3");
    fSB_c23=new TF1("fSB_c23",totdecaymodelforplot,&fitF::fcndecay1n_c23,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_c23");
    fSB2_c1=new TF1("fSB2_c1",totdecaymodelforplot,&fitF::fcndecay2n_c1,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_c1");
    fSB2_c2=new TF1("fSB2_c2",totdecaymodelforplot,&fitF::fcndecay2n_c2,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_c2");
    fSB2_c3=new TF1("fSB2_c3",totdecaymodelforplot,&fitF::fcndecay2n_c3,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_c3");
    fSB2_c4=new TF1("fSB2_c4",totdecaymodelforplot,&fitF::fcndecay2n_c4,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_c4");
    fSB2_c134=new TF1("fSB2_c134",totdecaymodelforplot,&fitF::fcndecay2n_c134,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_c134");

    writeFitComponents();
    plotResultsMore();
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
    ofs<<"nsig="<<nsig->getVal()<<"\tnbkg="<<nbkg->getVal()<<std::endl;
    fitres->Print();
    std::cout<<"Time for MC generation = "<<fMCGenTime<<std::endl;
    std::cout<<"Time for Fitting = "<<fFitTime<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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
   fB=new TF1("fB",totdecaymodel,&fitF::fcndecay,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay");
   for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
       fB->FixParameter(i,pCentralVal[i]);
   }
   fB->FixParameter(fdecaypath->getNMember()*5+4,100);
   fB->FixParameter(fdecaypath->getNMember()*5+5,0);


   fSB=new TF1("fSB",totdecaymodel,&fitF::fcndecay1n,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n");
   for (int i=0;i<fdecaypath->getNMember()*5+4;i++){
       fSB->FixParameter(i,pCentralVal[i]);
   }
   fSB->FixParameter(fdecaypath->getNMember()*5+4,10);
   fSB->FixParameter(fdecaypath->getNMember()*5+5,0);

   fSB2=new TF1("fSB2",totdecaymodel,&fitF::fcndecay2n,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n");
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

   //! construct parent/daugters decay components for plotting demonstration
   fB_parent=new TF1("fB_parent",totdecaymodel,&fitF::fcndecay_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay_parent");
   fB_daugter=new TF1("fB_daugter",totdecaymodel,&fitF::fcndecay_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay_daugter");

   fSB_parent=new TF1("fSB_parent",totdecaymodel,&fitF::fcndecay1n_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_parent");
   fSB_daugter=new TF1("fSB_daugter",totdecaymodel,&fitF::fcndecay1n_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_daugter");

   fSB2_parent=new TF1("fSB2_parent",totdecaymodel,&fitF::fcndecay2n_parent,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_parent");
   fSB2_daugter=new TF1("fSB2_daugter",totdecaymodel,&fitF::fcndecay2n_daugter,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_daugter");

   fSB_c1=new TF1("fSB_c1",totdecaymodel,&fitF::fcndecay1n_c1,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_c1");
   fSB_c2=new TF1("fSB_c2",totdecaymodel,&fitF::fcndecay1n_c2,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_c2");
   fSB_c3=new TF1("fSB_c3",totdecaymodel,&fitF::fcndecay1n_c3,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_c3");
   fSB_c23=new TF1("fSB_c23",totdecaymodel,&fitF::fcndecay1n_c23,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay1n_c23");
   fSB2_c1=new TF1("fSB2_c1",totdecaymodel,&fitF::fcndecay2n_c1,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_c1");
   fSB2_c2=new TF1("fSB2_c2",totdecaymodel,&fitF::fcndecay2n_c2,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_c2");
   fSB2_c3=new TF1("fSB2_c3",totdecaymodel,&fitF::fcndecay2n_c3,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_c3");
   fSB2_c4=new TF1("fSB2_c4",totdecaymodel,&fitF::fcndecay2n_c4,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_c4");
   fSB2_c134=new TF1("fSB2_c134",totdecaymodel,&fitF::fcndecay2n_c134,p_deadtime,p_timerange,fdecaypath->getNMember()*5+6,"fitF","fcndecay2n_c134");
   writeFitComponents();
   hB->Write();
   hSB->Write();
   hSB2->Write();
   plotResultsMore(1);

   for (int i=0;i<fnMC;i++){
       generateMC();
       setValParameters();
       fitter.Config().SetParamsSettings(fdecaypath->getNMember()*5+10,binfitparms);
       fitter.FitFCN(fdecaypath->getNMember()*5+10,globalChi2,0,dataB.Size()+dataSB.Size()+dataSB2.Size(),false);
       fitter.Result().Print(std::cout);
       const Double_t* resultpar=fitter.Result().GetParams();
       //const Double_t* resulterr=fitter.Result().GetErrors();
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
   writeOutputTree();


   closeOutputFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::writeFitComponents()
{
    for (int i=0;i<fdecaypath->getNMember()*5+6;i++){
        fB_parent->FixParameter(i,fB->GetParameter(i));
        fB_daugter->FixParameter(i,fB->GetParameter(i));
        fSB_parent->FixParameter(i,fSB->GetParameter(i));
        fSB_daugter->FixParameter(i,fSB->GetParameter(i));
        fSB2_parent->FixParameter(i,fSB2->GetParameter(i));
        fSB2_daugter->FixParameter(i,fSB2->GetParameter(i));

        fSB_c1->FixParameter(i,fSB->GetParameter(i));
        fSB_c2->FixParameter(i,fSB->GetParameter(i));
        fSB_c3->FixParameter(i,fSB->GetParameter(i));
        fSB_c23->FixParameter(i,fSB->GetParameter(i));
        fSB2_c1->FixParameter(i,fSB2->GetParameter(i));
        fSB2_c2->FixParameter(i,fSB2->GetParameter(i));
        fSB2_c3->FixParameter(i,fSB2->GetParameter(i));
        fSB2_c4->FixParameter(i,fSB2->GetParameter(i));
        fSB2_c134->FixParameter(i,fSB2->GetParameter(i));
    }

    fB->SetNpx(nbinsHB*10);
    fB_bkgneg->SetNpx(nbinsHB*10);
    fB_bkgpos->SetNpx(nbinsHB*10);
    fB_parent->SetNpx(nbinsHB*10);
    fB_daugter->SetNpx(nbinsHB*10);
    fSB->SetNpx(nbinsHSB*10);
    fSB_bkgneg->SetNpx(nbinsHSB*10);
    fSB_bkgpos->SetNpx(nbinsHSB*10);
    fSB_parent->SetNpx(nbinsHSB*10);
    fSB_daugter->SetNpx(nbinsHSB*10);
    fSB2->SetNpx(nbinsHSB2*10);
    fSB2_bkgneg->SetNpx(nbinsHSB2*10);
    fSB2_bkgpos->SetNpx(nbinsHSB2*10);
    fSB2_parent->SetNpx(nbinsHSB2*10);
    fSB2_daugter->SetNpx(nbinsHSB2*10);
    fSB_c1->SetNpx(nbinsHB*10);
    fSB_c2->SetNpx(nbinsHB*10);
    fSB_c3->SetNpx(nbinsHB*10);
    fSB_c23->SetNpx(nbinsHB*10);
    fSB2_c1->SetNpx(nbinsHSB2*10);
    fSB2_c2->SetNpx(nbinsHSB2*10);
    fSB2_c3->SetNpx(nbinsHSB2*10);
    fSB2_c4->SetNpx(nbinsHSB2*10);
    fSB2_c134->SetNpx(nbinsHSB2*10);

    fB->Write();
    fSB->Write();
    fSB2->Write();
    fB_bkgneg->Write();
    fSB_bkgneg->Write();
    fSB2_bkgneg->Write();
    fB_bkgpos->Write();
    fSB_bkgpos->Write();
    fSB2_bkgpos->Write();

    fB_parent->Write();
    fB_daugter->Write();
    fSB_parent->Write();
    fSB_daugter->Write();
    fSB2_parent->Write();
    fSB2_daugter->Write();

    fSB_c1->Write();
    fSB_c2->Write();
    fSB_c3->Write();
    fSB_c23->Write();
    fSB2_c1->Write();
    fSB2_c2->Write();
    fSB2_c3->Write();
    fSB2_c4->Write();
    fSB2_c134->Write();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void unbinfit::plotResultsMore(Int_t opt)
{
    gStyle->SetOptStat(0);

    /**
        canvas 0n
    */

    TCanvas* c0n=new TCanvas("c0n","c0n",1200,800);
    TPad *pad1_c0n = new TPad("pad1_c0n","pad1_c0n",0,0.3,1,1);
    TPad *pad2_c0n = new TPad("pad2_c0n","pad2_c0n",0,0,1,0.3);
    pad1_c0n->SetTopMargin(0.09);
    pad1_c0n->SetBottomMargin(0.00001);
    pad1_c0n->SetBorderMode(0);
    //pad1_c0n->SetLogy();
    pad2_c0n->SetTopMargin(0.00001);
    pad2_c0n->SetBottomMargin(0.4);
    pad2_c0n->SetBorderMode(0);
    pad1_c0n->Draw();
    pad2_c0n->Draw();
    pad1_c0n->cd();
    Int_t npremove=0;
    if (opt==0){
        TH1F* hdummyc0n=new TH1F("hdummy0n","",20,plotrangelow,plotrangehi);
        hdummyc0n->Draw();
        hdummyc0n->GetYaxis()->SetRangeUser(model0nHist->GetYaxis()->GetXmin(),model0nHist->GetYaxis()->GetXmax());
        hdummyc0n->GetYaxis()->SetTitleSize(0.06);
        hdummyc0n->GetYaxis()->SetTitleOffset(0.58);
        hdummyc0n->GetYaxis()->SetTitle("Counts");
        hdummyc0n->GetYaxis()->SetLabelSize(0.05);
        //! for unbin fit
        for (Int_t i=0;i<model0nCurve->GetN();i++) if (model0nCurve->GetX()[i]<p_deadtime) npremove++; else break;
        for (Int_t i=0;i<npremove;i++) model0nCurve->RemovePoint(0);
        model0nHist->SetMarkerSize(0.8);
        modelbkg0nHist->SetMarkerSize(0.8);
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
    fB_parent->SetLineStyle(9);
    fB_parent->Draw("same");
    fB_daugter->SetLineWidth(3);
    fB_daugter->SetLineColor(6);
    fB_daugter->SetLineStyle(10);
    fB_daugter->Draw("same");
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

    TGraphErrors * resplot_0n;
    TGraphErrors * resplotbkg_0n;

    if (opt==0){
        for (Int_t i=0;i<model0nHist->GetN();i++){
            Double_t xi=model0nHist->GetX()[i];
            Double_t yi=model0nHist->GetY()[i];
            Double_t yeval=model0nCurve->Eval(xi);
            Double_t reldev=yeval-yi;
            xres[i]=xi;
            yres[i]=reldev;
            yreserr[i]=sqrt(yi+yeval);
            Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
            chisquare+=chisquarei;
        }
        chisquare=2*chisquare;
        cout<<"ndf="<<fitres->floatParsFinal().getSize()<<endl;
        cout<<"chisquare/ndf="<<chisquare/(nbinsHB*2+fitres->floatParsFinal().getSize())<<endl;
        resplot_0n=new TGraphErrors(model0nHist->GetN(),xres,yres,0,yreserr);
        for (Int_t i=0;i<modelbkg0nHist->GetN();i++){
            Double_t xi=modelbkg0nHist->GetX()[i];
            Double_t yi=modelbkg0nHist->GetY()[i];
            Double_t yeval=fB_bkgneg->Eval(xi);
            Double_t reldev=yeval-yi;
            xres[i]=xi;
            yres[i]=reldev;
            yreserr[i]=sqrt(yi+yeval);
        }
        resplotbkg_0n=new TGraphErrors(modelbkg0nHist->GetN(),xres,yres,0,yreserr);
    }else{
        Int_t k=0;
        for (Int_t i=0;i<hB->GetNbinsX();i++){
            Double_t xi=hB->GetBinCenter(i+1);
            Double_t yi=hB->GetBinContent(i+1);
            if (xi>p_deadtime){
                Double_t yeval=fB->Eval(xi);
                Double_t reldev=yeval-yi;
                xres[k]=xi;
                yres[k]=reldev;
                yreserr[k]=sqrt(yi+yeval);
                Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
                chisquare+=chisquarei;
                k++;
            }
        }
        chisquare=2*chisquare;
        cout<<"k="<<k<<endl;
        cout<<"ndf="<<fitCovQual<<endl;
        cout<<"chisquare/ndf="<<chisquare/(k+fitCovQual)<<endl;
        resplot_0n=new TGraphErrors(k,xres,yres,0,yreserr);
        k=0;
        for (Int_t i=0;i<hB->GetNbinsX();i++){
           Double_t xi=hB->GetBinCenter(i+1);
           Double_t yi=hB->GetBinContent(i+1);
           if (xi<0){
               Double_t yeval=fB_bkgneg->Eval(xi);
               Double_t reldev=yeval-yi;
               xres[k]=xi;
               yres[k]=reldev;
               yreserr[k]=sqrt(yi+yeval);
               k++;
           }
       }
        resplotbkg_0n=new TGraphErrors(k,xres,yres,0,yreserr);
    }

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

    resplot_0n->Draw("sameP");
    resplotbkg_0n->Draw("sameP");
    pad2_c0n->Draw();
    c0n->Write();


    //!canvas 1n

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

    if (opt==0){
        TH1F* hdummyc1n=new TH1F("hdummy1n","",20,plotrangelow,plotrangehi);
        hdummyc1n->Draw();
        hdummyc1n->GetYaxis()->SetRangeUser(model1nHist->GetYaxis()->GetXmin(),model1nHist->GetYaxis()->GetXmax());

        hdummyc1n->GetYaxis()->SetTitleSize(0.06);
        hdummyc1n->GetYaxis()->SetTitleOffset(0.58);
        hdummyc1n->GetYaxis()->SetTitle("Counts");

        hdummyc1n->GetYaxis()->SetLabelSize(0.05);
        //! remove several points at the begining for the model curves
        npremove=0;
        for (Int_t i=0;i<model1nCurve->GetN();i++) if (model1nCurve->GetX()[i]<p_deadtime) npremove++; else break;
        for (Int_t i=0;i<npremove;i++) model1nCurve->RemovePoint(0);
        model1nHist->SetMarkerSize(0.8);
        modelbkg1nHist->SetMarkerSize(0.8);
        setdataHistzero(model1nHist);
        setdataHistzero(modelbkg1nHist);
        model1nHist->Draw("sameP");
        modelbkg1nHist->Draw("sameP");
        model1nCurve->SetLineWidth(3);
        model1nCurve->SetLineColor(4);
        model1nCurve->Draw("same");
    }else{
        hSB->SetMarkerStyle(20);
        hSB->SetMarkerSize(1.);
        hSB->Draw("P0E");
        hSB->GetXaxis()->SetRangeUser(plotrangelow,plotrangehi);
        hSB->GetYaxis()->SetTitleSize(0.06);
        hSB->GetYaxis()->SetTitleOffset(0.58);
        hSB->GetYaxis()->SetTitle("Counts");
        hSB->GetYaxis()->SetLabelSize(0.05);
        fSB->SetLineWidth(3);
        fSB->SetLineColor(4);
        fSB->Draw("same");
    }


    fSB_c23->SetLineWidth(3);
    fSB_c23->SetLineColor(2);
    fSB_c23->SetLineStyle(9);
    fSB_c23->Draw("same");
    fSB_c1->SetLineWidth(3);
    fSB_c1->SetLineColor(6);
    fSB_c1->SetLineStyle(10);
    fSB_c1->Draw("same");
    fSB_bkgneg->SetLineWidth(3);
    fSB_bkgneg->SetLineColor(4);
    fSB_bkgneg->Draw("same");
    pad1_c1n->Draw();

    pad2_c1n->cd();
    Double_t chisquare1n=0;

    TGraphErrors * resplot_1n;
    TGraphErrors * resplotbkg_1n;
    if (opt==0){
        for (Int_t i=0;i<model1nHist->GetN();i++){
            Double_t xi=model1nHist->GetX()[i];
            Double_t yi=model1nHist->GetY()[i];
            Double_t yeval=model1nCurve->Eval(xi);
            Double_t reldev=yeval-yi;
            xres[i]=xi;
            yres[i]=reldev;
            yreserr[i]=sqrt(yi+yeval);
            Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
            chisquare1n+=chisquarei;
        }
        cout<<"chisquare1n="<<chisquare1n<<endl;
        resplot_1n=new TGraphErrors(model1nHist->GetN(),xres,yres,0,yreserr);
        for (Int_t i=0;i<modelbkg1nHist->GetN();i++){
            Double_t xi=modelbkg1nHist->GetX()[i];
            Double_t yi=modelbkg1nHist->GetY()[i];
            Double_t yeval=fSB_bkgneg->Eval(xi);
            Double_t reldev=yeval-yi;
            xres[i]=xi;
            yres[i]=reldev;
            yreserr[i]=sqrt(yi+yeval);
        }
        resplotbkg_1n=new TGraphErrors(modelbkg1nHist->GetN(),xres,yres,0,yreserr);
    }else{
        Int_t k=0;
        for (Int_t i=0;i<hSB->GetNbinsX();i++){
            Double_t xi=hSB->GetBinCenter(i+1);
            Double_t yi=hSB->GetBinContent(i+1);
            if (xi>p_deadtime){
                Double_t yeval=fSB->Eval(xi);
                Double_t reldev=yeval-yi;
                xres[k]=xi;
                yres[k]=reldev;
                yreserr[k]=sqrt(yi+yeval);
                Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
                chisquare1n+=chisquarei;
                k++;
            }
        }
        cout<<"chisquare1n="<<chisquare1n<<endl;
        resplot_1n=new TGraphErrors(k,xres,yres,0,yreserr);

        k=0;
        for (Int_t i=0;i<hSB->GetNbinsX();i++){
           Double_t xi=hSB->GetBinCenter(i+1);
           Double_t yi=hSB->GetBinContent(i+1);
           if (xi<0){
               Double_t yeval=fSB_bkgneg->Eval(xi);
               Double_t reldev=yeval-yi;
               xres[k]=xi;
               yres[k]=reldev;
               yreserr[k]=sqrt(yi+yeval);
               k++;
           }
       }
        resplotbkg_1n=new TGraphErrors(k,xres,yres,0,yreserr);
    }


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

    //!canvas 2n

    TCanvas* c2n=new TCanvas("c2n","c2n",1400,700);
    TPad *pad1_c2n = new TPad("pad1_c2n","pad1_c2n",0,0.3,1,1);
    TPad *pad2_c2n = new TPad("pad2_c2n","pad2_c2n",0,0,1,0.3);
    pad1_c2n->SetTopMargin(0.09);
    pad1_c2n->SetBottomMargin(0.00001);
    pad1_c2n->SetBorderMode(0);
    //pad1_c2n->SetLogy();
    pad2_c2n->SetTopMargin(0.00001);
    pad2_c2n->SetBottomMargin(0.4);
    pad2_c2n->SetBorderMode(0);
    pad1_c2n->Draw();
    pad2_c2n->Draw();
    pad1_c2n->cd();
    if (opt==0){
        TH1F* hdummyc2n=new TH1F("hdummy2n","",20,plotrangelow,plotrangehi);
        hdummyc2n->Draw();
        hdummyc2n->GetYaxis()->SetRangeUser(model2nHist->GetYaxis()->GetXmin(),model2nHist->GetYaxis()->GetXmax());

        hdummyc2n->GetYaxis()->SetTitleSize(0.06);
        hdummyc2n->GetYaxis()->SetTitleOffset(0.58);
        hdummyc2n->GetYaxis()->SetTitle("Counts");

        hdummyc2n->GetYaxis()->SetLabelSize(0.05);
        //! remove several points at the begining for the model curves
        npremove=0;
        for (Int_t i=0;i<model2nCurve->GetN();i++) if (model2nCurve->GetX()[i]<p_deadtime) npremove++; else break;
        for (Int_t i=0;i<npremove;i++) model2nCurve->RemovePoint(0);
        model2nHist->SetMarkerSize(0.8);
        modelbkg2nHist->SetMarkerSize(0.8);
        setdataHistzero(model2nHist);
        setdataHistzero(modelbkg2nHist);
        model2nHist->Draw("sameP");
        modelbkg2nHist->Draw("sameP");
        model2nCurve->SetLineWidth(3);
        model2nCurve->SetLineColor(4);
        model2nCurve->Draw("same");
    }else{
        hSB2->SetMarkerStyle(20);
        hSB2->SetMarkerSize(1.);
        hSB2->Draw("P0E");
        hSB2->GetXaxis()->SetRangeUser(plotrangelow,plotrangehi);
        hSB2->GetYaxis()->SetTitleSize(0.06);
        hSB2->GetYaxis()->SetTitleOffset(0.58);
        hSB2->GetYaxis()->SetTitle("Counts");
        hSB2->GetYaxis()->SetLabelSize(0.05);
        fSB2->SetLineWidth(3);
        fSB2->SetLineColor(4);
        fSB2->Draw("same");
    }

    fSB2_c2->SetLineWidth(3);
    fSB2_c2->SetLineColor(2);
    fSB2_c2->SetLineStyle(9);
    fSB2_c2->Draw("same");
    fSB2_c134->SetLineWidth(3);
    fSB2_c134->SetLineColor(6);
    fSB2_c134->SetLineStyle(10);
    fSB2_c134->Draw("same");
    fSB2_bkgneg->SetLineWidth(3);
    fSB2_bkgneg->SetLineColor(4);
    fSB2_bkgneg->Draw("same");
    pad1_c2n->cd();

    pad2_c2n->cd();
    Double_t chisquare2n=0;
    TGraphErrors * resplot_2n;
    TGraphErrors * resplotbkg_2n;
    if (opt==0){
        for (Int_t i=0;i<model2nHist->GetN();i++){
            Double_t xi=model2nHist->GetX()[i];
            Double_t yi=model2nHist->GetY()[i];
            Double_t yeval=model2nCurve->Eval(xi);
            Double_t reldev=yeval-yi;
            xres[i]=xi;
            yres[i]=reldev;
            yreserr[i]=sqrt(yi+yeval);
            Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
            chisquare2n+=chisquarei;
        }
        resplot_2n=new TGraphErrors(model2nHist->GetN(),xres,yres,0,yreserr);
        for (Int_t i=0;i<modelbkg2nHist->GetN();i++){
            Double_t xi=modelbkg2nHist->GetX()[i];
            Double_t yi=modelbkg2nHist->GetY()[i];
            Double_t yeval=fSB2_bkgneg->Eval(xi);
            Double_t reldev=yeval-yi;
            xres[i]=xi;
            yres[i]=reldev;
            yreserr[i]=sqrt(yi+yeval);
        }
        resplotbkg_2n=new TGraphErrors(modelbkg2nHist->GetN(),xres,yres,0,yreserr);
    }else{
        Int_t k=0;
        for (Int_t i=0;i<hSB2->GetNbinsX();i++){
            Double_t xi=hSB2->GetBinCenter(i+1);
            Double_t yi=hSB2->GetBinContent(i+1);
            if (xi>p_deadtime){
                Double_t yeval=fSB2->Eval(xi);
                Double_t reldev=yeval-yi;
                xres[k]=xi;
                yres[k]=reldev;
                yreserr[k]=sqrt(yi+yeval);
                Double_t chisquarei=yeval-yi+yi*TMath::Log(yi/yeval);
                chisquare2n+=chisquarei;
                k++;
            }
        }
        resplot_2n=new TGraphErrors(k,xres,yres,0,yreserr);
        k=0;
        for (Int_t i=0;i<hSB2->GetNbinsX();i++){
           Double_t xi=hSB2->GetBinCenter(i+1);
           Double_t yi=hSB2->GetBinContent(i+1);
           if (xi<0){
               Double_t yeval=fSB2_bkgneg->Eval(xi);
               Double_t reldev=yeval-yi;
               xres[k]=xi;
               yres[k]=reldev;
               yreserr[k]=sqrt(yi+yeval);
               k++;
           }
       }
        resplotbkg_2n=new TGraphErrors(k,xres,yres,0,yreserr);
    }

    TH1F* hdummyc2nres=new TH1F("hdummyc2nres","",20,plotrangelow,plotrangehi);
    hdummyc2nres->SetLineColor(1);
    hdummyc2nres->SetLineWidth(2);
    hdummyc2nres->Draw();
    hdummyc2nres->GetYaxis()->SetRangeUser(resplot_2n->GetYaxis()->GetXmin(),resplot_2n->GetYaxis()->GetXmax());

    hdummyc2nres->GetYaxis()->SetLabelSize(0.12);
    hdummyc2nres->GetYaxis()->SetTitle("fit - data (counts)");
    hdummyc2nres->GetYaxis()->SetTitleSize(0.12);
    hdummyc2nres->GetYaxis()->SetTitleOffset(0.29);

    hdummyc2nres->GetXaxis()->SetLabelSize(0.14);
    hdummyc2nres->GetXaxis()->SetTitle("t_{#beta} - t_{ion} (s)");
    hdummyc2nres->GetXaxis()->SetTitleSize(0.17);
    hdummyc2nres->GetXaxis()->SetTitleOffset(0.95);

    resplot_2n->SetMarkerStyle(20);
    resplotbkg_2n->SetMarkerStyle(20);

    resplot_2n->SetLineColor(2);
    resplotbkg_2n->SetLineColor(2);

    resplot_2n->Draw("sameP");
    resplotbkg_2n->Draw("sameP");
    pad2_c2n->Draw();
    c2n->Write();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
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
