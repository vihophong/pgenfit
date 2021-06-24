/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
///  'MULTIDIMENSIONAL MODELS' RooFit tutorial macro #312
///
///  Performing fits in multiple (disjoint) ranges in one or more dimensions
///
/// \macro_output
/// \macro_code
/// \author 07/2008 - Wouter Verkerke


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
#include "fitF.cxx"
#include "fitFbkg.cxx"
#include "fitFdaughters.cxx"
#include "fitFparent.cxx"
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
#include "histofunction.C"

#include <fstream>

using namespace RooFit ;
using namespace RooStats;

Double_t p_deadtime=0;
Double_t p_timerange=10;

Int_t ncpus=24;

using namespace std;


void setdataHistzero(RooHist* datahist){
    for (Int_t i=0;i<datahist->GetN();i++){if (datahist->GetY()[i]==0) {
            datahist->SetPointEYhigh(i,0);
            datahist->SetPointEYlow(i,0);
        }
    }
}


void getparms(Double_t* parms,Double_t* parmserr, Double_t* parmsmax, Double_t* parmsmin, Int_t* isparmsfix,char* infile)
{
    std::ifstream ifs(infile);
    Double_t decayparms[knri][3];
    Double_t decayparms_err[knri][3];
    Double_t decayparms_max[knri][3];
    Double_t decayparms_min[knri][3];
    Int_t decayparms_isfix[knri][3];

    std::string tempstr;
    Int_t temp;
    Int_t rino=0;

    std::getline(ifs,tempstr);
    for (int i=0;i<knri;i++){
        ifs>>tempstr>>temp>>decayparms[rino][0]>>decayparms_err[rino][0]>>decayparms[rino][1]
                >>decayparms_err[rino][1]>>decayparms[rino][2]>>decayparms_err[rino][2]
                >>decayparms_isfix[rino][0]>>decayparms_min[rino][0]>>decayparms_max[rino][0]
                >>decayparms_isfix[rino][1]>>decayparms_min[rino][1]>>decayparms_max[rino][1]
                >>decayparms_isfix[rino][2]>>decayparms_min[rino][2]>>decayparms_max[rino][2];
        //cout<<tempstr<<"\t"<<temp<<"\t"<<decayparms[rino][0]<<"\t"<<decayparms_min[rino][0]<<"\t"<<decayparms_max[rino][0]<<endl;
        rino++;
    }

    for (int i=0;i<knri;i++){
        for (int j=0;j<3;j++){
            if (j==0){//for half-life
                decayparms_err[i][j]=TMath::Log(2)/decayparms[i][j]/decayparms[i][j]*decayparms_err[i][j];
                decayparms[i][j]=TMath::Log(2)/decayparms[i][j];
                Double_t maxdecayrate=TMath::Log(2)/decayparms_min[i][j];
                Double_t mindecayrate=TMath::Log(2)/decayparms_max[i][j];
                decayparms_min[i][j]=mindecayrate;
                decayparms_max[i][j]=maxdecayrate;
                //cout<<i<<"\t"<<decayparms[i][j]<<"\t"<<decayparms_min[i][j]<<endl;
            }
        }
    }


    //! calculate output
    for (int i=0;i<knri;i++){
        for (int j=0;j<2;j++){
            isparmsfix[j*knri+i]=decayparms_isfix[i][j];
            parms[j*knri+i]=decayparms[i][j];
            parmserr[j*knri+i]=decayparms_err[i][j];
            parmsmin[j*knri+i]=decayparms_min[i][j];
            parmsmax[j*knri+i]=decayparms_max[i][j];
        }
    }

    //! read p2n parrent
    parms[knri*2]=decayparms[0][2];
    parmserr[knri*2]=decayparms_err[0][2];
    parms[knri*2]=decayparms[0][2];
    parmsmin[knri*2]=decayparms_min[0][2];
    parmsmax[knri*2]=decayparms_max[0][2];
    isparmsfix[knri*2]=decayparms_isfix[0][2];
}

void mlhfitexpv1(char* infile,char* parmsfile,char*outfile,Double_t inputneueff,Double_t inputneueff_err,Double_t beginfit=0,Int_t fitopt=0,Int_t nbinsplot=200,Double_t plotrangelow=-1,Double_t plotrangehi=3)// fit option 0: no constrains(central value) fit option 1: constrains for sys error estimation; fit option 2: constrains for background only
{    
    p_deadtime=beginfit;
    //! Get paramters
    Double_t parms[knri*2+6];
    Double_t parmserr[knri*2+6];
    Double_t parmsmax[knri*2+6];
    Double_t parmsmin[knri*2+6];
    Int_t isparmsfix[knri*2+6];
    getparms(parms,parmserr,parmsmax,parmsmin,isparmsfix,parmsfile);
    for (int i=0;i<knri*2+1;i++){
        cout<<"parm "<<i<<"\t"<<parms[i]<<"\t"<<parmserr[i]<<"\t"<<parmsmin[i]<<"\t"<<parmsmax[i]<<"\t"<<isparmsfix[i]<<endl;
    }

    // Import tree to total decay data
    TFile *f=TFile::Open(infile);
    TTree* tree;
    f->GetObject("tree",tree);
    tree->SetEntries(5000);
    //tree->SetEntries(20000);

    //! define variables
    RooRealVar x("x","x",p_deadtime,p_timerange) ;
    RooRealVar xbkg("x","x",-p_timerange,-p_deadtime) ;

    cout<<"TIMERANGE: "<<p_deadtime<<"\t"<<p_timerange<<endl;

    //! define discrete variable y
    RooCategory y("y","y");
    y.defineType("0neu",0);
    y.defineType("1neu",1);
    y.defineType("2neu",2);


    //! estimate background
    Double_t nnbkg=tree->Draw("",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    Double_t nnbkg1n = tree->Draw("",Form("x<%f&&x>%f&&y==1",-p_deadtime,-p_timerange),"goff");
    Double_t nnbkg2n = tree->Draw("",Form("x<%f&&x>%f&&y==2",-p_deadtime,-p_timerange),"goff");

    cout<<"Background ratios = "<<nnbkg1n/nnbkg<<"\t"<<nnbkg2n/nnbkg1n<<endl;

    RooRealVar bkg1nratio("bkg1nratio","bkg1nratio",nnbkg1n/nnbkg,nnbkg1n/nnbkg/5,nnbkg1n/nnbkg*5) ;
    RooRealVar bkg2nratio("bkg2nratio","bkg2nratio",nnbkg2n/nnbkg1n,nnbkg2n/nnbkg1n/5,nnbkg2n/nnbkg1n*5) ;


    //! bkg pdf
    fitFbkg bkgmodel("bkgmodel","bkgmodel",xbkg,y,bkg1nratio,bkg2nratio);

    //! prepare data set for bkg fit
    RooDataSet* data2=new RooDataSet("data","data",RooArgSet(xbkg,y),Import(*tree)) ;
    data2->Print() ;

    //! fit background
    bkgmodel.fitTo(*data2,NumCPU(ncpus),Save()) ;

    fitFbkg bkgmodelfw("bkgmodelfw","bkgmodelfw",x,y,bkg1nratio,bkg2nratio);

    RooRealVar* p[knri*2+4];
    // decay parameters
    for (int i=0;i<knri*2+1;i++){
        p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),parms[i],parmsmin[i],parmsmax[i]);
    }

    //! get random coincidence

    TH1F* hdecayconly1neub;
    TH1F* hdecaycmorethan0neub;
    TH1F* hdecayconly2neub;
    TH1F* hdecayall;

    f->GetObject("hdecay1nbwd",hdecayconly1neub);
    f->GetObject("hdecaygt0nbwd",hdecaycmorethan0neub);
    f->GetObject("hdecay2nbwd",hdecayconly2neub);
    f->GetObject("hdecay",hdecayall);

    Double_t  nball=(Double_t) hdecayall->GetEntries();
    Double_t n1nbwd=(Double_t) hdecayconly1neub->GetEntries();
    Double_t gt0nbwd=(Double_t) hdecaycmorethan0neub->GetEntries();
    Double_t n2nbwd=(Double_t) hdecayconly2neub->GetEntries();

    // random coincicence parameters
    parms[knri*2+1]=n1nbwd/nball;
    parmserr[knri*2+1]=n1nbwd/nball*TMath::Sqrt(1/n1nbwd+1/nball);
    parmsmax[knri*2+1]=n1nbwd/nball*10;
    parmsmin[knri*2+1]=0;//n1nbwd/nball/10;
    isparmsfix[knri*2+1]=1;

    parms[knri*2+2]=gt0nbwd/nball;
    parmserr[knri*2+2]=gt0nbwd/nball*TMath::Sqrt(1/gt0nbwd+1/nball);
    parmsmax[knri*2+2]=gt0nbwd/nball*10;
    parmsmin[knri*2+2]=0;//gt0nbwd/nball/10;
    isparmsfix[knri*2+2]=1;

    parms[knri*2+3]=n2nbwd/nball;
    parmserr[knri*2+3]=n2nbwd/nball*TMath::Sqrt(1/n2nbwd+1/nball);
    parmsmax[knri*2+3]=n2nbwd/nball*10;
    parmsmin[knri*2+3]=0;//n2nbwd/nball/10;
    isparmsfix[knri*2+3]=1;

    p[knri*2+1]=new RooRealVar("p19","p19",parms[knri*2+1],parmsmin[knri*2+1],parmsmax[knri*2+1]);
    p[knri*2+2]=new RooRealVar("p20","p20",parms[knri*2+2],parmsmin[knri*2+2],parmsmax[knri*2+2]);
    p[knri*2+3]=new RooRealVar("p21","p21",parms[knri*2+3],parmsmin[knri*2+3],parmsmax[knri*2+3]);


    cout<<"p"<<knri*2+1<<"\t val="<<p[knri*2+1]->getVal()<<endl;
    cout<<"p"<<knri*2+2<<"\t val="<<p[knri*2+2]->getVal()<<endl;
    cout<<"p"<<knri*2+3<<"\t val="<<p[knri*2+3]->getVal()<<endl;

    RooArgSet externalconstrains;

    //! define roogaussian for error propagation
    RooGaussian* pconstr[knri*2+4];

    for (int i=0;i<knri*2+4;i++){
        pconstr[i]=new RooGaussian(Form("p%dconstr",i),Form("p%dconstr",i),*p[i],RooConst(parms[i]),RooConst(parmserr[i]));
        if (isparmsfix[i]!=0){
            if (fitopt==0) {
                //if (i<20)
                p[i]->setConstant();
                cout<<"set Constant for parameter p"<<i<<"\tmean="<<parms[i]<<"\tstd="<<parmserr[i]<<endl;
            }else if (fitopt==2){
                if (i<knri*2+1) {
                    p[i]->setConstant();
                    cout<<"set Constant for parameter p"<<i<<"\tmean="<<parms[i]<<"\tstd="<<parmserr[i]<<endl;
                }else{
                    externalconstrains.add(*pconstr[i]);
                    cout<<"set constrain model for p"<<i<<"\tmean="<<parms[i]<<"\tstd="<<parmserr[i]<<endl;
                }
            }else{
                if (parmserr[i]==0&&parmserr[i]==0){
                    cout<<"set Constant for parameter p"<<i<<"\tmean="<<parms[i]<<"\tstd="<<parmserr[i]<<endl;
                    p[i]->setConstant();
                }else{
                    if (isparmsfix[i]==1) {
                        externalconstrains.add(*pconstr[i]);
                        cout<<"set constrain model for p"<<i<<"\tmean="<<parms[i]<<"\tstd="<<parmserr[i]<<endl;
                    }else{
                        p[i]->setConstant();
                        cout<<"set Constant for parameter p"<<i<<"\tmean="<<parms[i]<<"\tstd="<<parmserr[i]<<endl;
                    }
                }
            }
        }
    }

    p[0]->setMin(6.93147e-12);
    p[0]->setMax(23.3407);
    std::ofstream ofsparms((char*)"oparms.txt",ios::out);
    for (Int_t i=0;i<knri*2+4;i++){
        ofsparms<<p[i]->getVal()<<"\t"<<p[i]->getError()<<"\t"<<p[i]->getMin()<<"\t"<<p[i]->getMax()<<"\t"<<"\t"<<p[i]->isConstant()<<endl;
    }
    //! neutron detection efficinecy
    //RooRealVar neueff("neueff","neueff",0.613384,0.,1.) ;
    RooRealVar neueff("neueff","neueff",inputneueff,0.,1.) ;
    RooGaussian* neueffconstr=new RooGaussian("neueffconstr","neueffconstr",neueff,RooConst(neueff.getVal()),RooConst(inputneueff_err));
    cout<<"Neutron efficiency: "<<neueff.getVal()<<"\t"<<inputneueff_err<<endl;

    if (fitopt==1) externalconstrains.add(*neueffconstr);
    else neueff.setConstant();

    //! sig pdf
    fitF totdecaymodel("totdecaymodel","totdecaymodel",x,y,neueff,*p[0],*p[1],*p[2],*p[3],*p[4],*p[5],*p[6],*p[7],*p[8],*p[9],*p[10],*p[11],*p[12],*p[13],*p[14],*p[15],*p[16],*p[17],*p[18],*p[19],*p[20],*p[21]);

    //! set background/signal counts
    Double_t nnsig=tree->Draw("",Form("x<%f&&x>%f",p_timerange,p_deadtime),"goff");
    nnsig=nnsig-nnbkg;
    cout<<"NSIG = "<<nnsig<<"-NBKG = "<<nnbkg<<endl;

    RooRealVar nbkg("nbkg","nbkg",nnbkg,0,nnbkg*3) ;
    RooRealVar nsig("nsig","nsig",nnsig,nnsig*0.5,nnsig*1.5) ;
    // nbackgrounds and bkg ratio from previous background fit
    RooGaussian* nbkgconstr=new RooGaussian("nbkgconstr","nbkgconstr",nbkg,RooConst(nbkg.getVal()),RooConst(TMath::Sqrt(nbkg.getVal())));
    RooGaussian* bkg1nratiocnstr=new RooGaussian("bkg1nratiocnstr","bkg1nratiocnstr",bkg1nratio,RooConst(bkg1nratio.getVal()),RooConst(bkg1nratio.getError()));
    RooGaussian* bkg2nratiocnstr=new RooGaussian("bkg2nratiocnstr","bkg2nratiocnstr",bkg2nratio,RooConst(bkg2nratio.getVal()),RooConst(bkg2nratio.getError()));

    externalconstrains.add(*nbkgconstr);
    externalconstrains.add(*bkg1nratiocnstr);
    externalconstrains.add(*bkg2nratiocnstr);

    RooAddPdf final_pdf("final_pdf","final pdf",RooArgList(totdecaymodel,bkgmodelfw),RooArgList(nsig,nbkg));
    //! prepare data set for fitting forward correlated data
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(x,y),Import(*tree)) ;
    data->Print() ;
    // Fit
    RooFitResult* fitres;

    if (fitopt==0){
        neueff.setConstant();
        nbkg.setConstant();
        bkg1nratio.setConstant();
        bkg2nratio.setConstant();
        fitres=final_pdf.fitTo(*data,Extended(kTRUE),NumCPU(ncpus),Save(kTRUE)) ;
    }else{
        fitres=final_pdf.fitTo(*data,Minos(kTRUE), Extended(kTRUE), ExternalConstraints(externalconstrains),NumCPU(ncpus),Save(kTRUE)) ;
    }


//    fitFparent parentmodel("parentmodel","parentmodel",x,y,neueff,*p[0],*p[1],*p[2],*p[3],*p[4],*p[5],*p[6],*p[7],*p[8],*p[9],*p[10],*p[11],*p[12],*p[13],*p[14],*p[15],*p[16],*p[17],*p[18],*p[19],*p[20],*p[21]);
//    fitFdaughters daughtmodel("daughtmodel","daughtmodel",x,y,neueff,*p[0],*p[1],*p[2],*p[3],*p[4],*p[5],*p[6],*p[7],*p[8],*p[9],*p[10],*p[11],*p[12],*p[13],*p[14],*p[15],*p[16],*p[17],*p[18],*p[19],*p[20],*p[21]);

//    Double_t igparent=parentmodel.createIntegral(x)->getVal();
//    Double_t igdaughter=daughtmodel.createIntegral(x)->getVal();

//    RooRealVar nsigparent("nsigparent","nsigparent",nsig.getVal()*igparent/(igparent+igdaughter),nsig.getVal()*igparent/(igparent+igdaughter)*0.5,nsig.getVal()*igparent/(igparent+igdaughter)*1.5) ;
//    RooAddPdf final_pdfparent("final_pdfparent","final_pdfparent pdf",RooArgList(parentmodel,bkgmodelfw),RooArgList(nsigparent,nbkg));

    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1","c1",1200, 800);
    c1->Divide(2,3);

    c1->cd(1);
    RooPlot* xframe0 = x.frame(Title("0 neutron fit")) ;
    data->plotOn(xframe0,Cut("y==y::0neu"),Binning(nbinsplot),RooFit::Name("data0n")) ;
    final_pdf.plotOn(xframe0,Slice(y,"0neu"),LineColor(kRed),RooFit::Name("data0nmodel")) ;
    //data->plotOn(xframe0,Binning(nbinsplot),RooFit::Name("data0n")) ;
    //final_pdf.plotOn(xframe0,LineColor(kRed),RooFit::Name("data0nmodel")) ;
    //final_pdf.plotOn(xframe0,Components(bkgmodelfw),LineColor(kGreen),RooFit::Name("data0nmodelbgfw")) ;
    xframe0->Draw();

    c1->cd(3);
    RooPlot* xframe1 = x.frame(Title("1 neutron fit")) ;
    data->plotOn(xframe1,Cut("y==y::1neu"),Binning(nbinsplot),RooFit::Name("data1n")) ;
    final_pdf.plotOn(xframe1,Slice(y,"1neu"),LineColor(kRed),RooFit::Name("data1nmodel")) ;
    xframe1->Draw() ;

    c1->cd(5);
    RooPlot* xframe2 = x.frame(Title("2 neutron fit")) ;
    data->plotOn(xframe2,Cut("y==y::2neu"),Binning(nbinsplot),RooFit::Name("data2n")) ;
    final_pdf.plotOn(xframe2,Slice(y,"2neu"),LineColor(kRed),RooFit::Name("data2nmodel")) ;    
    xframe2->Draw() ;

    RooCurve* model0nCurve=(RooCurve*)xframe0->getCurve("data0nmodel");
    //RooCurve* model0nbgfwCurve=(RooCurve*)xframe0->getCurve("data0nmodelbgfw");

    RooHist* model0nHist=(RooHist*)xframe0->getHist("data0n");
    RooCurve* model1nCurve=(RooCurve*)xframe1->getCurve("data1nmodel");
    RooHist* model1nHist=(RooHist*)xframe1->getHist("data1n");
    RooCurve* model2nCurve=(RooCurve*)xframe2->getCurve("data2nmodel");
    RooHist* model2nHist=(RooHist*)xframe2->getHist("data2n");

    c1->cd(2);
    RooPlot* xframe3 = xbkg.frame(Title("0 neutron fit")) ;
    data2->plotOn(xframe3,Binning(nbinsplot),RooFit::Name("bkg0n")) ;
    bkgmodel.plotOn(xframe3,RooFit::Name("bkg0nmodel")) ;
    xframe3->Draw() ;
    c1->cd(4);
    RooPlot* xframe4 = xbkg.frame(Title("1 neutron fit")) ;
    data2->plotOn(xframe4,Cut("y==y::1neu"),Binning(nbinsplot),RooFit::Name("bkg1n")) ;
    bkgmodel.plotOn(xframe4,Slice(y,"1neu"),RooFit::Name("bkg1nmodel")) ;
    xframe4->Draw() ;
    c1->cd(6);
    RooPlot* xframe5 = xbkg.frame(Title("2 neutron fit")) ;
    data2->plotOn(xframe5,Cut("y==y::2neu"),Binning(nbinsplot),RooFit::Name("bkg2n")) ;
    bkgmodel.plotOn(xframe5,Slice(y,"2neu"),RooFit::Name("bkg2nmodel")) ;
    xframe5->Draw() ;

    RooCurve* modelbkg0nCurve=(RooCurve*)xframe3->getCurve("bkg0nmodel");
    RooHist* modelbkg0nHist=(RooHist*)xframe3->getHist("bkg0n");
    RooCurve* modelbkg1nCurve=(RooCurve*)xframe4->getCurve("bkg1nmodel");
    RooHist* modelbkg1nHist=(RooHist*)xframe4->getHist("bkg1n");
    RooCurve* modelbkg2nCurve=(RooCurve*)xframe5->getCurve("bkg2nmodel");
    RooHist* modelbkg2nHist=(RooHist*)xframe5->getHist("bkg2n");


    //!Calculate and residual plot - X2/ndf
    Double_t chisquare=0;
    Double_t xres[10000];
    Double_t yres[10000];
    Double_t yreserr[10000];
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
        //cout<<"reldev "<<xi<<"\t"<<reldev<<endl;
    }
    chisquare=2*chisquare;
    cout<<"ndf="<<fitres->floatParsFinal().getSize()<<endl;
    cout<<"chisquare/ndf="<<chisquare/(nbinsplot+fitres->floatParsFinal().getSize())<<endl;
    TGraphErrors * resplot_0n=new TGraphErrors(model0nHist->GetN(),xres,yres,0,yreserr);

    Double_t chisquare1n=0;
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
    TGraphErrors * resplot_1n=new TGraphErrors(model1nHist->GetN(),xres,yres,0,yreserr);

    Double_t chisquare2n=0;
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
    TGraphErrors * resplot_2n=new TGraphErrors(model2nHist->GetN(),xres,yres,0,yreserr);


    for (Int_t i=0;i<modelbkg0nHist->GetN();i++){
        Double_t xi=modelbkg0nHist->GetX()[i];
        Double_t yi=modelbkg0nHist->GetY()[i];
        Double_t yeval=modelbkg0nCurve->Eval(xi);
        Double_t reldev=yeval-yi;
        xres[i]=xi;
        yres[i]=reldev;
        yreserr[i]=sqrt(yi+yeval);
    }
    TGraphErrors * resplotbkg_0n=new TGraphErrors(modelbkg0nHist->GetN(),xres,yres,0,yreserr);

    for (Int_t i=0;i<modelbkg1nHist->GetN();i++){
        Double_t xi=modelbkg1nHist->GetX()[i];
        Double_t yi=modelbkg1nHist->GetY()[i];
        Double_t yeval=modelbkg1nCurve->Eval(xi);
        Double_t reldev=yeval-yi;
        xres[i]=xi;
        yres[i]=reldev;
        yreserr[i]=sqrt(yi+yeval);
    }
    TGraphErrors * resplotbkg_1n=new TGraphErrors(modelbkg1nHist->GetN(),xres,yres,0,yreserr);

    for (Int_t i=0;i<modelbkg2nHist->GetN();i++){
        Double_t xi=modelbkg2nHist->GetX()[i];
        Double_t yi=modelbkg2nHist->GetY()[i];
        Double_t yeval=model2nCurve->Eval(xi);
        Double_t reldev=yeval-yi;
        xres[i]=xi;
        yres[i]=reldev;
        yreserr[i]=sqrt(yi+yeval);
    }
    TGraphErrors * resplotbkg_2n=new TGraphErrors(modelbkg2nHist->GetN(),xres,yres,0,yreserr);


    ofstream str("out.txt",ios::app);
    str<<outfile<<"\t"<<TMath::Log(2)/p[0]->getVal()*1000<<"\t"<<p[9]->getVal()*100<<"\t"<<p[18]->getVal()*100<<"\t"<<log(2)/p[0]->getVal(0)/p[0]->getVal(0)*p[0]->getError()*1000<<"\t"<<p[9]->getError()*100<<"\t"<<p[18]->getError()*100<<"\t"<<chisquare/(nbinsplot+fitres->floatParsFinal().getSize())<<"\t"<<fitres->status()<<endl;

    /**
        Make another canvases for thesis

        canvas 0n
    */
    //!
    TCanvas* c0n=new TCanvas("c0n","c0n",1400,700);

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
    TH1F* hdummyc0n=new TH1F("hdummy0n","",20,plotrangelow,plotrangehi);
    hdummyc0n->Draw();
    hdummyc0n->GetYaxis()->SetRangeUser(model0nHist->GetYaxis()->GetXmin(),model0nHist->GetYaxis()->GetXmax());

    hdummyc0n->GetYaxis()->SetTitleSize(0.06);
    hdummyc0n->GetYaxis()->SetTitleOffset(0.58);
    hdummyc0n->GetYaxis()->SetTitle("Counts");

    hdummyc0n->GetYaxis()->SetLabelSize(0.05);

    //! remove several points at the begining for the model curves
    Int_t npremove=0;
    for (Int_t i=0;i<model0nCurve->GetN();i++) if (model0nCurve->GetX()[i]<p_deadtime) npremove++; else break;
    for (Int_t i=0;i<npremove;i++) model0nCurve->RemovePoint(0);
    //! remove several points at the end of the background model curves
    npremove=0;
    for (Int_t i=modelbkg0nCurve->GetN()-1;i>0;i--) if (modelbkg0nCurve->GetX()[i]>-p_deadtime) npremove++; else break;
    for (Int_t i=0;i<npremove;i++) modelbkg0nCurve->RemovePoint(modelbkg0nCurve->GetN()-1);

    model0nHist->SetMarkerSize(0.8);
    modelbkg0nHist->SetMarkerSize(0.8);

    model0nHist->Draw("sameP");
    model0nCurve->Draw("same");
    modelbkg0nHist->Draw("sameP");
    modelbkg0nCurve->SetLineColor(6);
    modelbkg0nCurve->SetLineStyle(10);
    modelbkg0nCurve->Draw("same");



    //! fake drawing
    Double_t parmhist[knri*2+3];
    for (Int_t i=0;i<knri*2+1;i++) parmhist[i]=p[i]->getVal();
    parmhist[knri*2+2]=modelbkg0nCurve->Eval((-p_timerange-p_deadtime)/2);
    parmhist[knri*2+1]=model0nCurve->Eval(p_deadtime)-parmhist[knri*2+2];
    TF1* tfdecay=getFunction(parmhist,p_deadtime,p_timerange,neueff.getVal());
    model0nCurve->Fit(tfdecay,"EQR+","goff");

    for (Int_t i=0;i<knri*2+3;i++) parmhist[i]=tfdecay->GetParameter(i);
    TF1* tfdaughterdecay=getDaughterFunction(parmhist,p_deadtime,p_timerange,neueff.getVal());
    TF1* tfparentdecay=getParentFunction(parmhist,p_deadtime,p_timerange,neueff.getVal());
    TF1* tfbkg=getBkgFunction(parmhist,p_deadtime,p_timerange,neueff.getVal());
    pad1_c0n->cd();
    //tfdecay->SetLineColor(3);
    //tfdecay->Draw("same");
    tfdaughterdecay->SetLineColor(4);
    tfdaughterdecay->Draw("same");
    tfparentdecay->SetLineColor(5);
    tfparentdecay->Draw("same");
    tfbkg->SetLineColor(6);
    tfbkg->SetLineStyle(10);
    tfbkg->Draw("same");
    pad1_c0n->Draw();


    pad2_c0n->cd();
    TH1F* hdummyc0nres=new TH1F("hdummyc0nres","",20,plotrangelow,plotrangehi);
    hdummyc0nres->SetLineColor(1);
    hdummyc0nres->SetLineWidth(2);
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


    /**
        Make another canvases for thesis

        canvas 1n
    */
    //!
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
    //! remove several points at the end of the background model curves
    npremove=0;
    for (Int_t i=modelbkg1nCurve->GetN()-1;i>0;i--) if (modelbkg1nCurve->GetX()[i]>-p_deadtime) npremove++; else break;
    for (Int_t i=0;i<npremove;i++) modelbkg1nCurve->RemovePoint(modelbkg1nCurve->GetN()-1);

    model1nHist->SetMarkerSize(0.8);
    modelbkg1nHist->SetMarkerSize(0.8);

    setdataHistzero(model1nHist);
    setdataHistzero(modelbkg1nHist);
    model1nHist->Draw("sameP");
    model1nCurve->Draw("same");
    modelbkg1nHist->Draw("sameP");
    modelbkg1nCurve->Draw("same");
    pad2_c1n->cd();
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

    /**
        Make another canvases for thesis

        canvas 2n
    */
    //!
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
    //! remove several points at the end of the background model curves
    npremove=0;
    for (Int_t i=modelbkg2nCurve->GetN()-1;i>0;i--) if (modelbkg2nCurve->GetX()[i]>-p_deadtime) npremove++; else break;
    for (Int_t i=0;i<npremove;i++) modelbkg2nCurve->RemovePoint(modelbkg2nCurve->GetN()-1);

    model2nHist->SetMarkerSize(0.8);
    modelbkg2nHist->SetMarkerSize(0.8);

    setdataHistzero(model2nHist);
    setdataHistzero(modelbkg2nHist);
    model2nHist->Draw("sameP");
    model2nCurve->Draw("same");
    modelbkg2nHist->Draw("sameP");
    modelbkg2nCurve->Draw("same");
    pad2_c2n->cd();
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

    TFile* outfileroot=new TFile(outfile,"recreate");

    model0nCurve->Write();
    model1nCurve->Write();
    model2nCurve->Write();

    model0nHist->Write();
    model1nHist->Write();
    model2nHist->Write();

    modelbkg0nCurve->Write();
    modelbkg1nCurve->Write();
    modelbkg2nCurve->Write();

    modelbkg0nHist->Write();
    modelbkg1nHist->Write();
    modelbkg2nHist->Write();

    resplot_0n->SetName("plot0nres");
    resplot_0n->Write();
    resplot_1n->SetName("plot1nres");
    resplot_1n->Write();
    resplot_2n->SetName("plot2nres");
    resplot_2n->Write();

    resplotbkg_0n->SetName("plotbkg0nres");
    resplotbkg_1n->SetName("plotbkg1nres");
    resplotbkg_2n->SetName("plotbkg2nres");

    resplotbkg_0n->Write();
    resplotbkg_1n->Write();
    resplotbkg_2n->Write();

    c1->Write();
    c0n->Write();
    c1n->Write();
    c2n->Write();

    tfdecay->Write();
    outfileroot->Close();
}


