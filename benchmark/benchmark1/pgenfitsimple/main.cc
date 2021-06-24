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

#include <iostream>
#include <fstream>
#include "fitF.hh"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "TRandom3.h"
#include "TFile.h"
#include "RooGaussian.h"
#include "RooConstVar.h"


#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAddition.h"
#include "RooDataHist.h"
#include "RooPoisson.h"
#include "RooPlot.h"

#include "RooNLLVar.h"

#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooStats/ProofConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooProfileLL.h"

#include "RooHist.h"
#include "RooCurve.h"

#include <fstream>

#include "decaypath.hh"

using namespace RooFit ;
using namespace RooStats;

int main(int argc, char *argv[])
{
    Double_t p_deadtime=0;
    Double_t p_timerange=10;
    // Import tree to total decay data
    TFile *f=TFile::Open("outhist.root");
    TTree* tree;
    f->GetObject("tree",tree);
    tree->SetEntries(5000);

    //! define variables
    RooRealVar x("x","x",p_deadtime,p_timerange) ;
    RooRealVar xbkg("x","x",-p_timerange,-p_deadtime) ;

    //! define discrete variable y
    RooCategory y("y","y");
    y.defineType("0neu",0);
    y.defineType("1neu",1);
    y.defineType("2neu",2);

    RooRealVar bkg1nratio("bkg1nratio","bkg1nratio",0.5,0,1) ;
    RooRealVar bkg2nratio("bkg2nratio","bkg2nratio",0.5,0,1) ;

    RooRealVar slope1("slope1","slope1",0.,-0.1,0.1) ;
    RooRealVar slope2("slope2","slope2",0.,-0.1,0.1) ;
    RooRealVar slope3("slope3","slope3",0.,-0.1,0.1) ;


    //! bkg pdf
    fitFbkg bkgmodel("bkgmodel","bkgmodel",xbkg,y,bkg1nratio,bkg2nratio,slope1,slope2,slope3);

    //! prepare data set for bkg fit
    RooDataSet* data2=new RooDataSet("data","data",RooArgSet(xbkg,y),Import(*tree)) ;
    data2->Print() ;

    //! fit background
    bkgmodel.fitTo(*data2,NumCPU(8),Save()) ;

    RooRealVar slope1pos("slope1pos","slope1pos",-slope1.getVal(),-0.1,0.1) ;
    RooRealVar slope2pos("slope2pos","slope2pos",-slope2.getVal(),-0.1,0.1) ;
    RooRealVar slope3pos("slope3pos","slope3pos",-slope3.getVal(),-0.1,0.1) ;

    //! set background ratio and slope as constant
    bkg1nratio.setConstant();
    bkg2nratio.setConstant();

    slope1pos.setConstant();
    slope2pos.setConstant();
    slope3pos.setConstant();

    //! bkg pdf positive
    fitFbkg bkgmodelpos("bkgmodelpos","bkgmodelpos",x,y,bkg1nratio,bkg2nratio,slope1pos,slope2pos,slope3pos);

    //! get random coincidence
    TH1F* hdecayconly1neub;
    TH1F* hdecaycmorethan0neub;
    TH1F* hdecayconly2neub;

    TH1F* hdecayall;

    f->GetObject("hdecayconly1neub",hdecayconly1neub);
    f->GetObject("hdecaycmorethan0neub",hdecaycmorethan0neub);
    f->GetObject("hdecayconly2neub",hdecayconly2neub);

    f->GetObject("hdecay",hdecayall);

    Double_t  nball=hdecayall->GetEntries();
    Double_t n1nbwd=hdecayconly1neub->GetEntries();
    Double_t gt0nbwd=hdecaycmorethan0neub->GetEntries();
    Double_t n2nbwd=hdecayconly2neub->GetEntries();

    //! set background/signal counts
    Double_t nnsig=tree->Draw("",Form("x>%f",p_deadtime),"goff");
    Double_t nnbkg=tree->Draw("",Form("x<%f",-p_deadtime),"goff");
    nnsig=nnsig-nnbkg;

    RooRealVar nbkg("nbkg","nbkg",nnbkg,0,nnbkg*10) ;
    RooRealVar nsig("nsig","nsig",nnsig,0,nnsig*10) ;
    nbkg.setConstant();

    //! prepare data set for fitting forward correlated data
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(x,y),Import(*tree)) ;
    data->Print() ;

    //! prepare fit paramters
    // Declare all parameters
    RooAbsReal* p[kmaxparms];
    RooRealVar* pvar[kmaxparms];

    decaypath* fdecaypath=new decaypath;

    fdecaypath->Init((char*)"benchmarkIn134.txt");
    fdecaypath->makePath();
    fdecaypath->printMember();
    fdecaypath->printPath();
    fdecaypath->writePath();
    fdecaypath->drawPath((char*)TString("outdecayroutes.root").Data());

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
    std::ifstream infile((char*)"effparmsbenchmark.txt");
    Int_t nlinesread = 0;
    Double_t be,b1ne,b2ne,n1n2ne;
    Double_t err_be,err_b1ne,err_b2ne,err_n1n2ne;
    Int_t isvary_be,isvary_b1ne,isvary_b2ne,isvary_n1n2ne;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (line[0]=='#') continue;
        if (!(iss >> be >> err_be >> b1ne >> err_b1ne >> b2ne >> err_b2ne >> n1n2ne >> err_n1n2ne)) break;
        nlinesread++;
    }

    if (be<0) {isvary_be=1;be=-be;}else{isvary_be=0;}
    if (b1ne<0){isvary_b1ne=1;b1ne=-b1ne;}else{isvary_b1ne=0;}
    if (b2ne<0){isvary_b2ne=1;b2ne=-b2ne;}else{isvary_b2ne=0;}
    if (n1n2ne<0){isvary_n1n2ne=1;n1n2ne=-n1n2ne;}else{isvary_n1n2ne=0;}

    std::cout<<"read-in efficiency factors:"<<std::endl;
    std::cout<<be<<"\t"<<err_be<<"\t"<<b1ne<<"\t"<<err_b1ne<<"\t"<<b2ne<<"\t"<<err_b2ne<<"\t"<<n1n2ne<<"\t"<<err_n1n2ne<<"\t"<<std::endl;

    p[fdecaypath->getNMember()*5+4]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5+4),Form("p%d",fdecaypath->getNMember()*5+4),be,be-10*err_be,be+10*err_be);
    p[fdecaypath->getNMember()*5+5]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5+5),Form("p%d",fdecaypath->getNMember()*5+5),b1ne,b1ne-10*err_b1ne,b1ne+10*err_b1ne);
    p[fdecaypath->getNMember()*5+6]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5+6),Form("p%d",fdecaypath->getNMember()*5+6),b2ne,b2ne-10*err_b2ne,b2ne+10*err_b2ne);
    p[fdecaypath->getNMember()*5+7]=new RooRealVar(Form("p%d",fdecaypath->getNMember()*5+7),Form("p%d",fdecaypath->getNMember()*5+7),n1n2ne,n1n2ne-10*err_n1n2ne,n1n2ne+10*err_n1n2ne);
    pvar[fdecaypath->getNMember()*5+4]=(RooRealVar*) p[fdecaypath->getNMember()*5+4];
    pvar[fdecaypath->getNMember()*5+5]=(RooRealVar*) p[fdecaypath->getNMember()*5+5];
    pvar[fdecaypath->getNMember()*5+6]=(RooRealVar*) p[fdecaypath->getNMember()*5+6];
    pvar[fdecaypath->getNMember()*5+7]=(RooRealVar*) p[fdecaypath->getNMember()*5+7];
    pvar[fdecaypath->getNMember()*5+4]->setError(err_be);
    pvar[fdecaypath->getNMember()*5+5]->setError(err_b1ne);
    pvar[fdecaypath->getNMember()*5+6]->setError(err_b2ne);
    pvar[fdecaypath->getNMember()*5+7]->setError(err_n1n2ne);

    if (isvary_be==0) pvar[fdecaypath->getNMember()*5+4]->setConstant(kTRUE);
    if (isvary_b1ne==0) pvar[fdecaypath->getNMember()*5+5]->setConstant(kTRUE);
    if (isvary_b2ne==0) pvar[fdecaypath->getNMember()*5+6]->setConstant(kTRUE);
    if (isvary_n1n2ne==0) pvar[fdecaypath->getNMember()*5+7]->setConstant(kTRUE);

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
    slope1pos.setConstant();
    slope2pos.setConstant();
    slope3pos.setConstant();
    // backgrounds
    nbkg.setConstant(kTRUE);
    // background ratio and slope
    bkg1nratio.setConstant();
    bkg2nratio.setConstant();

    // random coincicence parameters
    pvar[fdecaypath->getNMember()*5+1]->setConstant(kTRUE);//randcoinf1n
    pvar[fdecaypath->getNMember()*5+2]->setConstant(kTRUE);//randcoinfgt0n
    pvar[fdecaypath->getNMember()*5+3]->setConstant(kTRUE);//randcoinf2n

    fitF totdecaymodel("totdecaymodel","totdecaymodel",x,y,p);
    RooAddPdf final_pdf("final_pdf","final pdf",RooArgList(totdecaymodel,bkgmodelpos),RooArgList(nsig,nbkg));

    std::ofstream ofsparms((char*)"oparms.txt",ios::out);
    pvar[0]->setMin(6.93147e-12);
    pvar[0]->setMax(23.3407);
    for (Int_t i=0;i<fdecaypath->getNMember()*5+7;i++){
        ofsparms<<pvar[i]->getVal()<<"\t"<<pvar[i]->getError()<<"\t"<<pvar[i]->getMin()<<"\t"<<pvar[i]->getMax()<<"\t"<<pvar[i]->isConstant()<<endl;
    }
    // Fit
    RooFitResult* fitres;
    fitres=final_pdf.fitTo(*data,NumCPU(8),Save()) ;

    return 0;
}
