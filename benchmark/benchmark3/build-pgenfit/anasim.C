#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TAxis.h"

void anasim(char* inputfile)
{
    TFile* f=TFile::Open(inputfile);
    TTree* treecorr=(TTree*)f->Get("treecorr");
    //! all
    treecorr->Draw("beta.T-ion.T>>hall(2000,-10,10)");
    TH1F* hall=(TH1F*) gDirectory->Get("hall");
    treecorr->Draw("beta.T-ion.T>>hreal(2000,-10,10)","beta.evt==ion.evt");
    TH1F* hreal=(TH1F*) gDirectory->Get("hreal");
    treecorr->Draw("beta.T-ion.T>>hparent(2000,-10,10)","beta.evt==ion.evt&&beta.id==0");
    TH1F* hparent=(TH1F*) gDirectory->Get("hparent");
    treecorr->Draw("beta.T-ion.T>>hdaughter(2000,-10,10)","beta.evt==ion.evt&&beta.id>0");
    TH1F* hdaughter=(TH1F*) gDirectory->Get("hdaughter");
    treecorr->Draw("beta.T-ion.T>>hbkg(2000,-10,10)","beta.evt!=ion.evt");
    TH1F* hbkg=(TH1F*) gDirectory->Get("hbkg");


    //! 1n
    treecorr->Draw("beta.T-ion.T>>h1n(2000,-10,10)","multf==1");
    TH1F* h1n=(TH1F*) gDirectory->Get("h1n");
    treecorr->Draw("beta.T-ion.T>>h1nreal(2000,-10,10)","multf==1&&beta.evt==ion.evt");
    TH1F* h1nreal=(TH1F*) gDirectory->Get("h1nreal");
    treecorr->Draw("beta.T-ion.T>>h1n1(2000,-10,10)","multf==1&&beta.evt==ion.evt&&beta.evt!=evtf");
    TH1F* h1n1=(TH1F*) gDirectory->Get("h1n1");
    treecorr->Draw("beta.T-ion.T>>h1n2(2000,-10,10)","multf==1&&beta.evt==ion.evt&&beta.evt==evtf&&beta.mode==1");
    TH1F* h1n2=(TH1F*) gDirectory->Get("h1n2");
    treecorr->Draw("beta.T-ion.T>>h1n3(2000,-10,10)","multf==1&&beta.evt==ion.evt&&beta.evt==evtf&&beta.mode==2");
    TH1F* h1n3=(TH1F*) gDirectory->Get("h1n3");

    //! 2n
    treecorr->Draw("beta.T-ion.T>>h2n(2000,-10,10)","multf==2");
    TH1F* h2n=(TH1F*) gDirectory->Get("h2n");
    treecorr->Draw("beta.T-ion.T>>h2nreal(2000,-10,10)","multf==2&&beta.evt==ion.evt");
    TH1F* h2nreal=(TH1F*) gDirectory->Get("h2nreal");
    treecorr->Draw("beta.T-ion.T>>h2n1(2000,-10,10)","multf==2&&beta.evt==ion.evt&&beta.evt!=evtf[0]&&beta.evt!=evtf[1]");
    TH1F* h2n1=(TH1F*) gDirectory->Get("h2n1");
    treecorr->Draw("beta.T-ion.T>>h2n2(2000,-10,10)","multf==2&&beta.evt==ion.evt&&beta.evt==evtf[0]&&beta.evt==evtf[1]&&modef[0]==2&&modef[1]==2");
    TH1F* h2n2=(TH1F*) gDirectory->Get("h2n2");
    treecorr->Draw("beta.T-ion.T>>h2n3(2000,-10,10)","multf==2&&beta.evt==ion.evt&&((beta.evt==evtf[0]&&beta.evt!=evtf[1]&&modef[0]==1)||(beta.evt!=evtf[0]&&beta.evt==evtf[1]&&modef[1]==1))");
    TH1F* h2n3=(TH1F*) gDirectory->Get("h2n3");
    treecorr->Draw("beta.T-ion.T>>h2n4(2000,-10,10)","multf==2&&beta.evt==ion.evt&&((beta.evt==evtf[0]&&beta.evt!=evtf[1]&&modef[0]==2)||(beta.evt!=evtf[0]&&beta.evt==evtf[1]&&modef[1]==2))");
    TH1F* h2n4=(TH1F*) gDirectory->Get("h2n4");

    //! Canvas
    TCanvas *c1=new TCanvas("c1","c1",900,700);
    c1->cd();
    hall->SetLineColor(1);
    hreal->SetLineColor(2);
    hparent->SetLineColor(3);
    hdaughter->SetLineColor(4);
    hbkg->SetLineColor(5);
    hall->Draw();hall->SetMinimum(0);
    hreal->Draw("same");
    hparent->Draw("same");
    hdaughter->Draw("same");
    hbkg->Draw("same");

    TCanvas *c2=new TCanvas("c2","c2",900,700);
    c2->cd();
    h1n->SetLineColor(1);
    h1nreal->SetLineColor(2);
    h1n1->SetLineColor(3);
    h1n2->SetLineColor(4);
    h1n3->SetLineColor(5);
    h1n->Draw();h1n->SetMinimum(0);
    h1nreal->Draw("same");
    h1n1->Draw("same");
    h1n2->Draw("same");
    h1n3->Draw("same");

    TCanvas *c3=new TCanvas("c3","c3",900,700);
    c3->cd();
    h2n->SetLineColor(1);
    h2nreal->SetLineColor(2);
    h2n1->SetLineColor(3);
    h2n2->SetLineColor(4);
    h2n3->SetLineColor(5);
    h2n4->SetLineColor(6);
    h2n->Draw();h2n->SetMinimum(0);
    h2nreal->Draw("same");
    h2n1->Draw("same");
    h2n2->Draw("same");
    h2n3->Draw("same");
    h2n4->Draw("same");
}
