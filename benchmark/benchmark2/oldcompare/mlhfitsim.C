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
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooCategory.h"
#include "TRandom3.h"
#include "TFile.h"
#include "fitF.cxx"
#include "fitFbkg.cxx"
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

#include <fstream>

using namespace RooFit ;
using namespace RooStats;

Double_t p_deadtime=0;
Double_t p_timerange=10;

using namespace std;


void getparms(Double_t* parms,Double_t* parmserr, Double_t* parmsmax, Double_t* parmsmin, Bool_t* isparmsfix,char* infile,Double_t nsig=10.)
{
    std::ifstream ifs(infile);
    Int_t rino;
    Double_t temp;
    Bool_t flagfix[knri][3];
    Double_t decayparms[knri][3];
    Double_t decayparms_err[knri][3];

    for (int i=0;i<knri;i++){
        ifs>>rino;
        for (int j=0;j<3;j++){
            ifs>>temp;
            if (temp>=0){
                flagfix[i][j]=true;
                decayparms[i][j]=temp;
            }else{
                flagfix[i][j]=false;
                decayparms[i][j]=(-temp);
            }
            ifs>>temp;decayparms_err[i][j]=temp;
        }
    }

    for (int i=0;i<knri;i++){
        for (int j=0;j<3;j++){
            if (j==0){//for half-life
                decayparms_err[i][j]=TMath::Log(2)/decayparms[i][j]/decayparms[i][j]*decayparms_err[i][j];
                decayparms[i][j]=TMath::Log(2)/decayparms[i][j];
            }
        }
    }


    //! calculate output
    for (int i=0;i<knri;i++){
        for (int j=0;j<2;j++){
            isparmsfix[j*knri+i]=flagfix[i][j];
            parms[j*knri+i]=decayparms[i][j];
            parmserr[j*knri+i]=decayparms_err[i][j];
            parmsmax[j*knri+i]=decayparms[i][j]+decayparms_err[i][j]*nsig;
            if ((decayparms[i][j]-decayparms_err[i][j]*nsig)>0)
                parmsmin[j*knri+i]=decayparms[i][j]-decayparms_err[i][j]*nsig;
            else
                parmsmin[j*knri+i]=TMath::Log(2)/100000000000;
            //for pn
            if (j!=0) {parmsmin[j*knri+i]=0;parmsmax[j*knri+i]=1;}
        }
    }

    //! read p2n parrent
    parms[knri*2]=decayparms[0][2];
    parmserr[knri*2]=decayparms_err[0][2];
    parms[knri*2]=decayparms[0][2];
    parmsmin[knri*2]=0;parmsmax[knri*2]=1;
    isparmsfix[knri*2]=flagfix[0][2];

    //! read initial activity, background and random coincidence factor
    for (int i=0;i<5;i++){
        ifs>>parms[knri*2+i+1]>>parmsmin[knri*2+i+1]>>parmsmax[knri*2+i+1];
        parmserr[knri*2+i+1]=parms[knri*2+i+1]-parmsmin[knri*2+i+1];
        cout<<knri*2+i+1<<"\t"<<parms[knri*2+i+1]<<"\t"<<parmsmin[knri*2+i+1]<<"\t"<<parmsmax[knri*2+i+1]<<endl;
        isparmsfix[knri*2+i+1]=false;
        if (i==4) {
            parmsmin[knri*2+i+1]=0;parmsmax[knri*2+i+1]=1;
            if (parms[knri*2+i+1]>=0){
                isparmsfix[knri*2+i+1]=true;
            }else{
                parms[knri*2+i+1]=-parms[knri*2+i+1];
                isparmsfix[knri*2+i+1]=false;
            }
        }
    }
}

void mlhfitsim()
{
    //! Get paramters

    Double_t parms[knri*2+6];
    Double_t parmserr[knri*2+6];
    Double_t parmsmax[knri*2+6];
    Double_t parmsmin[knri*2+6];
    Bool_t isparmsfix[knri*2+6];

    getparms(parms,parmserr,parmsmax,parmsmin,isparmsfix,"parmsfit.txt",5);
    for (int i=0;i<knri*2+6;i++){
        cout<<"parm "<<i<<"\t"<<parms[i]<<"\t"<<parmserr[i]<<"\t"<<parmsmin[i]<<"\t"<<parmsmax[i]<<"\t"<<isparmsfix[i]<<endl;
    }

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


    //! estimate background
    Double_t nnbkg=tree->Draw("",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    Double_t nnbkg1n = tree->Draw("",Form("x<%f&&x>%f&&y==1",-p_deadtime,-p_timerange),"goff");
    Double_t nnbkg2n = tree->Draw("",Form("x<%f&&x>%f&&y==2",-p_deadtime,-p_timerange),"goff");

    cout<<"Background ratios = "<<nnbkg1n/nnbkg<<"\t"<<nnbkg2n/nnbkg1n<<endl;

    RooRealVar bkg1nratio("bkg1nratio","bkg1nratio",nnbkg1n/nnbkg,nnbkg1n/nnbkg/5,nnbkg1n/nnbkg*5) ;
    RooRealVar bkg2nratio("bkg2nratio","bkg2nratio",nnbkg2n/nnbkg1n,nnbkg2n/nnbkg1n/5,nnbkg2n/nnbkg1n*5) ;

    RooRealVar slope1("slope1","slope1",0.,-0.1,0.1) ;
    RooRealVar slope2("slope2","slope2",0.,-0.1,0.1) ;
    RooRealVar slope3("slope3","slope3",0.,-0.1,0.1) ;
    slope1.setConstant(kTRUE);
    slope2.setConstant(kTRUE);
    slope3.setConstant(kTRUE);


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

    RooRealVar* p[knri*2+4];
    // decay parameters
    for (int i=0;i<knri*2+1;i++){
        p[i]=new RooRealVar(Form("p%d",i),Form("p%d",i),parms[i],parmsmin[i],parmsmax[i]);
        if (isparmsfix[i]==1){
            cout<<"fixed value p"<<i<<"\tval="<<parms[i]<<endl;
            p[i]->setConstant(kTRUE);
        }else{
            cout<<"variable value p"<<i<<"\tval "<<parms[i]<<"\tmin="<<parmsmin[i]<<"\tmax="<<parmsmax[i]<<endl;
        }
    }

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

    // random coincicence parameters
    p[knri*2+1]=new RooRealVar("p19","p19",n1nbwd/nball,n1nbwd/nball/3,n1nbwd/nball*3);
    p[knri*2+2]=new RooRealVar("p20","p20",gt0nbwd/nball,gt0nbwd/nball/3,gt0nbwd/nball*3);
    p[knri*2+3]=new RooRealVar("p21","p21",n2nbwd/nball,n2nbwd/nball/3,n2nbwd/nball*3);
    p[knri*2+1]->setConstant(kTRUE);
    p[knri*2+2]->setConstant(kTRUE);
    p[knri*2+3]->setConstant(kTRUE);

    //! neutron detection efficinecy
    //RooRealVar neueff("neueff","neueff",0.613384,0.,1.) ;
    RooRealVar neueff("neueff","neueff",0.62,0.,1.) ;
    neueff.setConstant(kTRUE);

    //! sig pdf
    fitF totdecaymodel("totdecaymodel","totdecaymodel",x,y,neueff,*p[0],*p[1],*p[2],*p[3],*p[4],*p[5],*p[6],*p[7],*p[8],*p[9],*p[10],*p[11],*p[12],*p[13],*p[14],*p[15],*p[16],*p[17],*p[18],*p[19],*p[20],*p[21]);

    //! set background/signal counts
    Double_t nnsig=tree->Draw("",Form("x>%f&&x<%f",p_deadtime,p_timerange),"goff");
    //Double_t nnbkg=tree->Draw("",Form("x<%f&&x>%f",-p_deadtime,-p_timerange),"goff");
    nnsig=nnsig-nnbkg;

    RooRealVar nbkg("nbkg","nbkg",nnbkg,0,nnbkg*3) ;
    RooRealVar nsig("nsig","nsig",nnsig*0.5,nnsig*1.5) ;
    nbkg.setConstant();

    RooAddPdf final_pdf("final_pdf","final pdf",RooArgList(totdecaymodel,bkgmodelpos),RooArgList(nsig,nbkg));

    std::ofstream ofsparms((char*)"oparms.txt",ios::out);
    for (Int_t i=0;i<knri*2+4;i++){
        ofsparms<<p[i]->getVal()<<"\t"<<p[i]->getError()<<"\t"<<p[i]->getMin()<<"\t"<<p[i]->getMax()<<"\t"<<"\t"<<p[i]->isConstant()<<endl;
    }
    //! prepare data set for fitting forward correlated data
    RooDataSet* data=new RooDataSet("data","data",RooArgSet(x,y),Import(*tree)) ;
    data->Print() ;

    //return;
    // Fit
    RooFitResult* fitres;
    fitres=final_pdf.fitTo(*data,NumCPU(8),Save()) ;


    TCanvas *c1 = new TCanvas("c1","c1",1200, 800);
    c1->Divide(2,3);
    c1->cd(1);
    RooPlot* xframe0 = x.frame(Title("0 neutron fit")) ;
    data->plotOn(xframe0,Cut("y==y::0neu"),Binning(500)) ;
    final_pdf.plotOn(xframe0,Slice(y,"0neu")) ;
    final_pdf.plotOn(xframe0,Components(bkgmodelpos),LineColor(kGreen),RooFit::Name("data0nmodelbgfw")) ;
    xframe0->Draw() ;
    c1->cd(3);
    RooPlot* xframe1 = x.frame(Title("1 neutron fit")) ;
    data->plotOn(xframe1,Cut("y==y::1neu"),Binning(500)) ;
    final_pdf.plotOn(xframe1,Slice(y,"1neu")) ;
    xframe1->Draw() ;
    c1->cd(5);
    RooPlot* xframe2 = x.frame(Title("2 neutron fit")) ;
    data->plotOn(xframe2,Cut("y==y::2neu"),Binning(500)) ;
    final_pdf.plotOn(xframe2,Slice(y,"2neu")) ;
    xframe2->Draw() ;



    c1->cd(2);
    RooPlot* xframe3 = xbkg.frame(Title("0 neutron fit")) ;
    data2->plotOn(xframe3,Cut("y==y::0neu"),Binning(500)) ;
    bkgmodel.plotOn(xframe3,Slice(y,"0neu")) ;
    xframe3->Draw() ;
    c1->cd(4);
    RooPlot* xframe4 = xbkg.frame(Title("1 neutron fit")) ;
    data2->plotOn(xframe4,Cut("y==y::1neu"),Binning(500)) ;
    bkgmodel.plotOn(xframe4,Slice(y,"1neu")) ;
    xframe4->Draw() ;
    c1->cd(6);
    RooPlot* xframe5 = xbkg.frame(Title("2 neutron fit")) ;
    data2->plotOn(xframe5,Cut("y==y::2neu"),Binning(500)) ;
    bkgmodel.plotOn(xframe5,Slice(y,"2neu")) ;
    xframe5->Draw() ;




    ofstream str("out.txt",ios::app);
    str<<p[0]->getVal()<<"\t"<<p[9]->getVal()<<"\t"<<p[18]->getVal()<<"\t"<<p[0]->getError()<<"\t"<<p[9]->getError()<<"\t"<<p[18]->getError()<<endl;


}
