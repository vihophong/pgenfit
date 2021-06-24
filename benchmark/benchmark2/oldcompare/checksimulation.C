#include <math.h>
#include <stdlib.h>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <map>

typedef struct {
    double T; 	 // Calibrated time
    double Tcorr; //correlated time
    double x,y,z;// number of pixel for AIDA, or number of tube for BELEN
    int type;
    int type2;
    int evt;
} datatype;

void checksimulation()
{
    TRandom3 r;
    //! read back data
    double windowlow=10;
    double windowup=10;
    Double_t deltaxy=4.;
    cout<<"deltaxy = "<<deltaxy<<endl;


    //! veto gate neutron - beam
    double vetoneutronbeamlow=0;
    double vetoneutronbeamup=400000./1e9;

    //! window gate beta - neutron
    double windowbetaneutronlow=400000./1e9;
    double windowbetaneutronup=400000./1e9;//0

    //! artificial dead time after implantation

    double tdead=0.002;//2ms
    //Int_t tdead=0;//2ms


    std::multimap < double,std::pair< unsigned int, datatype> > fimplantMap;
    std::multimap < double,std::pair< unsigned int, datatype> >::iterator fimplantMap_it;

    std::multimap < double, unsigned int> fbetaMap;
    std::multimap < double, unsigned int>::iterator fbetaMap_it;

    std::multimap < double,std::pair< unsigned int, datatype> > fneuMap;
    std::multimap < double,std::pair< unsigned int, datatype> >::iterator fneuMap_it;

    datatype readbackion;
    datatype readbackbeta;
    datatype readbackneu;

    TTree* readbacktreebeta=0;
    TTree* readbacktreeion=0;
    TTree* readbacktreeneu=0;

    TFile* f = new TFile("outtree.root");
    f->GetObject("ion",readbacktreeion);
    f->GetObject("beta",readbacktreebeta);
    f->GetObject("neutron",readbacktreeneu);

    readbacktreeion->SetBranchAddress("ion",&readbackion);
    readbacktreebeta->SetBranchAddress("beta",&readbackbeta);
    readbacktreeneu->SetBranchAddress("neutron",&readbackneu);


    TFile* fout=new TFile("outhist.root","RECREATE");
    Double_t tp1n=0;
    Double_t tp2n=0;
    Double_t tall=0;
    Int_t mult=0;
    Int_t nflag=0;
    Int_t nrealflag=0;
    Int_t nfwd=0;
    Int_t nbwd=0;

    Int_t btype=-1;
    Int_t breal=-1;
    Double_t neu_T=-999999;
    Double_t neuf_T=-999999;
    Double_t neub_T=-999999;

    TTree* treeb=new TTree("treeb","treeb");
    treeb->Branch("x",&tall,"x/D");
    treeb->Branch("nflag",&nflag,"nflag/I");
    treeb->Branch("nrealflag",&nrealflag,"nrealflag/I");
    treeb->Branch("nfwd",&nfwd,"nfwd/I");
    treeb->Branch("nbwd",&nbwd,"nbwd/I");
    treeb->Branch("btype",&btype,"btype/I");
    treeb->Branch("breal",&breal,"breal/I");
    treeb->Branch("neu_T",&neu_T,"neu_T/D");
    treeb->Branch("neuf_T",&neuf_T,"neuf_T/D");
    treeb->Branch("neub_T",&neub_T,"neub_T/D");

    TTree* treep1n=new TTree("treep1n","treep1n");
    treep1n->Branch("x",&tp1n,"x/D");
    TTree* treep2n=new TTree("treep2n","treep2n");
    treep2n->Branch("x",&tp2n,"x/D");

    TTree* treemlh=new TTree("tree","tree");
    treemlh->Branch("x",&tall,"x/D");
    treemlh->Branch("y",&mult,"y/I");




    TTree* treeptmp=new TTree("treeptmp","treeptmp");
    datatype datatmp;
    datatype datatmparr[20];
    treeptmp->Branch("neu",&datatmp,"T/D:Tcorr/D:x/D:y/D:z/D:type/I:type2/I:evt/I");


    TH1F* hdecayc=new TH1F("hdecay","hdecay",2000,-windowlow/1,windowup/1);
    TH1F* hdecaycomponent[17];
    for (Int_t i=0;i<17;i++){
        hdecaycomponent[i]=new TH1F(Form("hdecay%d",i+1),Form("hdecay%d",i+1),2000,-windowlow/1,windowup/1);
    }
    TH1F* hdecaycreal=new TH1F("hdecaycreal","hdecaycreal",2000,-windowlow/1,windowup/1);

    TH1F* hbkggaus=new TH1F("hbkggaus","hbkggaus",2000,-windowlow/1,windowup/1);
    TH1F* hbkguniform=new TH1F("hbkguniform","hbkguniform",2000,-windowlow/1,windowup/1);

    TH2F* hdeltaxy=new TH2F("hdeltaxy","hdeltaxy",200,-4,4,200,-4,4);

    TH1F* hneubeam=new TH1F("hneubeam","hneubeam",2000,-vetoneutronbeamlow,vetoneutronbeamup);
    TH1F* hneubeamreal=new TH1F("hneubeamreal","hneubeamreal",2000,-vetoneutronbeamlow,vetoneutronbeamup);

    TH1F* hneubeta=new TH1F("hneubeta","hneubeta",2000,-windowbetaneutronup,windowbetaneutronlow);
    TH1F* hneubetareal=new TH1F("hneubetareal","hneubetareal",2000,-windowbetaneutronup,windowbetaneutronlow);

    TH1F* hneuneu=new TH1F("hneuneu","hneuneu",2000,-windowbetaneutronup,windowbetaneutronlow);


    TH1F* hdecayc1neu=new TH1F("hdecayc1neu","hdecayc1neu",2000,-windowlow/1,windowup/1);
    TH1F* hdecayc1neub=new TH1F("hdecayc1neub","hdecayc1neub",2000,-windowlow/1,windowup/1);
    TH1F* hdecayconly1neu=new TH1F("hdecayconly1neu","hdecayconly1neu",2000,-windowlow/1,windowup/1);
    TH1F* hdecayconly1neub=new TH1F("hdecayconly1neub","hdecayconly1neub",2000,-windowlow/1,windowup/1);


    TH1F* hdecay1nbwd=new TH1F("hdecay1nbwd","hdecay1nbwd",2000,-windowlow/1,windowup/1);
    TH1F* hdecaygt0nbwd=new TH1F("hdecaygt0nbwd","hdecaygt0nbwd",2000,-windowlow/1,windowup/1);
    TH1F* hdecay2nbwd=new TH1F("hdecay2nbwd","hdecay2nbwd",2000,-windowlow/1,windowup/1);


    TH1F* hdecayconly2neu=new TH1F("hdecayconly2neu","hdecayconly2neu",2000,-windowlow/1,windowup/1);
    TH1F* hdecayconly2neub=new TH1F("hdecayconly2neub","hdecayconly2neub",2000,-windowlow/1,windowup/1);


    TH1F* hdecaycmorethan0neu=new TH1F("hdecaycmorethan0neu","hdecaycmorethan0neu",2000,-windowlow/1,windowup/1);
    TH1F* hdecaycmorethan0neub=new TH1F("hdecaycmorethan0neub","hdecaycmorethan0neub",2000,-windowlow/1,windowup/1);

    TH1F* hdecaycmorethan1neu=new TH1F("hdecaycmorethan1neu","hdecaycmorethan1neu",2000,-windowlow/1,windowup/1);
    TH1F* hdecaycmorethan1neub=new TH1F("hdecaycmorethan1neub","hdecaycmorethan1neub",2000,-windowlow/1,windowup/1);

    TH1F* hdecayc1neureal=new TH1F("hdecayc1neureal","hdecayc1neureal",2000,-windowlow/1,windowup/1);

    TH1F* hdecayc1neurealin2n=new TH1F("hdecayc1neurealin2n","hdecayc1neurealin2n",2000,-windowlow/1,windowup/1);
    TH1F* hdecayc2neuin1neuin2n=new TH1F("hdecayc2neuin1neuin2n","hdecayc1neurealin2n",2000,-windowlow/1,windowup/1);


    TH1F* hdecayc2neuin1neu=new TH1F("hdecayc2neuin1neu","hdecayc2neuin1neu",2000,-windowlow/1,windowup/1);

    TH1F* hdecayc1neuin2neu=new TH1F("hdecayc1neuin2neu","hdecayc1neuin2neu",2000,-windowlow/1,windowup/1);


    TH1F* hdecayc2neu=new TH1F("hdecayc2neu","hdecayc2neu",2000,-windowlow/1,windowup/1);
    TH1F* hdecayc2neureal=new TH1F("hdecayc2neureal","hdecayc2neureal",2000,-windowlow/1,windowup/1);
    TH1F* hdecayc2neub=new TH1F("hdecayc2neub","hdecayc2neub",2000,-windowlow/1,windowup/1);


    TH1F* hdecayc1neureal3=new TH1F("hdecayc1neureal3","hdecayc1neureal3",2000,-windowlow,windowup);
    TH1F* hdecayc1neureal4=new TH1F("hdecayc1neureal4","hdecayc1neureal4",2000,-windowlow,windowup);
    TH1F* hdecayc1neureal8=new TH1F("hdecayc1neureal8","hdecayc1neureal8",2000,-windowlow,windowup);
    TH1F* hdecayc1neureal10=new TH1F("hdecayc1neureal10","hdecayc1neureal10",2000,-windowlow,windowup);
    TH1F* hdecayc1neureal14=new TH1F("hdecayc1neureal14","hdecayc1neureal14",2000,-windowlow,windowup);


    TH1F* betatdiff=new TH1F("betatdiff","betatdiff",20000,0,1000);

    Long64_t nentriesion=readbacktreeion->GetEntries();
    Long64_t nentriesbeta=readbacktreebeta->GetEntries();
    Long64_t nentriesneu=readbacktreeneu->GetEntries();


    for (unsigned int jentry=0;jentry<nentriesion;jentry++){
        readbacktreeion->GetEntry(jentry);
        datatype aidadata;
        aidadata.x=readbackion.x;
        aidadata.y=readbackion.y;
        aidadata.z=readbackion.z;
        aidadata.T=readbackion.T;
        aidadata.type=readbackion.type;
        aidadata.type2=readbackion.type2;
        aidadata.evt=readbackion.evt;
        fimplantMap.insert(make_pair((double)(aidadata.T*1),make_pair(jentry,aidadata)));
    }

    Double_t tsprev=0;
    for (unsigned int jentry=0;jentry<nentriesbeta;jentry++){
        readbacktreebeta->GetEntry(jentry);
        bool flagfill=true;
        if (readbackbeta.T-tsprev>0.4e-3) flagfill=true; else flagfill=false;
        tsprev=readbackbeta.T;

        flagfill=true;
        if (flagfill) fbetaMap.insert(make_pair((double)(readbackbeta.T*1),jentry));
    }

    tsprev=0;
    for (unsigned int jentry=0;jentry<nentriesneu;jentry++){
        readbacktreeneu->GetEntry(jentry);
        bool flagfill=true;
        if (readbackneu.T-tsprev>0.0006) flagfill=true; else flagfill=false;
        tsprev=readbackneu.T;

        datatype neudata;
        neudata.x=readbackneu.x;
        neudata.y=readbackneu.y;
        neudata.z=readbackneu.z;
        neudata.T=readbackneu.T;
        neudata.type=readbackneu.type;
        neudata.type2=0;//default no corr with beam
        neudata.evt=readbackneu.evt;

        flagfill=true;
        if (flagfill) fneuMap.insert(make_pair((double)(neudata.T*1),make_pair(jentry,neudata)));
    }
    cout<<"\ntotal number of implant = "<<fimplantMap.size()<<endl;
    cout<<"total number of beta = "<<fbetaMap.size()<<endl;
    cout<<"total number of neutron = "<<fneuMap.size()<<endl;

    /*
    //! scan neutron correlation
    for (fneuMap_it=fneuMap.begin();fneuMap_it!=fneuMap.end();fneuMap_it++){
        double ts=fneuMap_it->first;
        unsigned int entry = fneuMap_it->second.first;
        datatype neudata=fneuMap_it->second.second;
        double ts1 = ts - vetoneutronbeamup;
        double ts2 = ts + vetoneutronbeamlow;
        double corrts = 0;
        unsigned int correntry = 0;
        double check_time = 0;
        fimplantMap_it = fimplantMap.lower_bound(ts1);
        while(fimplantMap_it!=fimplantMap.end()&&fimplantMap_it->first<ts2){
            corrts =  fimplantMap_it->first;
            correntry = fimplantMap_it->second.first;
            if (corrts!=check_time){
                neudata.type2=1;//has correlation with beam
                hneubeam->Fill(ts-corrts);
                if (neudata.type==0) hneubeamreal->Fill(ts-corrts);
                break;
            }
            fimplantMap_it++;
        }
    }
    */

    int k=0;
    int ktotal=fbetaMap.size();
    Long64_t ncorr = 0;


    double tprev=0;

    for (fbetaMap_it=fbetaMap.begin();fbetaMap_it!=fbetaMap.end();fbetaMap_it++){
        //if (k%1000==0) cout<<k<<"/"<<ktotal<<"\tncorr="<<ncorr<<endl;
        //if (k>100000) break;
        double ts=(double)fbetaMap_it->first;
        //if ((ts-tprev)/1e3<600) continue;
        betatdiff->Fill((ts-tprev)/1e6);
        tprev=ts;
        unsigned int entry = fbetaMap_it->second;
        readbacktreebeta->GetEntry(entry);
        double betax = readbackbeta.x;
        double betay = readbackbeta.y;
        unsigned short betaz = readbackbeta.z;
        int betatype=readbackbeta.type;
        int betaevt=readbackbeta.evt;

        //int lastneutype=-9999;
        //! with neutron forward
        int ndelayedneutronf=0;
        double ts1 = ts - 0;
        double ts2 = ts + windowbetaneutronlow;
        double corrts = 0;
        unsigned int correntry = 0;
        double check_time = 0;
        fneuMap_it = fneuMap.lower_bound(ts1);

        int flagrealneutron=0;

        neu_T=-99999;
        neuf_T=-99999;
        neub_T=-99999;

        while(fneuMap_it!=fneuMap.end()&&fneuMap_it->first<ts2){
            corrts = fneuMap_it->first;
            correntry = fneuMap_it->second.first;
            datatype neudata=fneuMap_it->second.second;
            //if (corrts!=check_time){
            //if (corrts!=check_time&&neudata.type2==0){//Gate on if neutron is not associate with beam
            if (corrts!=check_time){
                hneubeta->Fill(corrts-ts);
                if (betaevt==neudata.evt&&readbackbeta.type==neudata.type){
                    hneubetareal->Fill(corrts-ts);
                    flagrealneutron++;
                }
                neu_T=corrts-ts;
                neuf_T=corrts-ts;
                ndelayedneutronf++;
            }
            fneuMap_it++;
        }

        //! with neutron backward
        int ndelayedneutronb=0;
        ts1 = ts - windowbetaneutronup;
        ts2 = ts + 0;
        corrts = 0;
        correntry = 0;
        check_time = 0;
        fneuMap_it = fneuMap.lower_bound(ts1);
        while(fneuMap_it!=fneuMap.end()&&fneuMap_it->first<ts2){
            corrts = fneuMap_it->first;
            correntry = fneuMap_it->second.first;
            datatype neudata=fneuMap_it->second.second;
            //if (corrts!=check_time){
            //if (corrts!=check_time&&neudata.type2==0){//Gate on if neutron is not associate with beam
            if (corrts!=check_time){
                hneubeta->Fill(corrts-ts);
                neu_T=corrts-ts;
                neub_T=corrts-ts;
                ndelayedneutronb++;
            }
            fneuMap_it++;
        }


        //! with implantation
        int fimplanti=0;
        ts1 = ts - windowup;
        ts2 = ts + windowlow;
        corrts = 0;
        correntry = 0;
        check_time = 0;
        fimplantMap_it = fimplantMap.lower_bound(ts1);

        Int_t flagbeta1=0;
        if (betatype==1||betatype==4||betatype==6) flagbeta1=1;

        while(fimplantMap_it!=fimplantMap.end()&&fimplantMap_it->first<ts2){
            corrts =  (double) fimplantMap_it->first;
            correntry = fimplantMap_it->second.first;
            datatype iont = fimplantMap_it->second.second;

            if (corrts!=check_time&&iont.z==betaz&&iont.x>=0&&iont.y>=0){
                //Double_t fdeltaxy = sqrt((iont.x-betax)*(iont.x-betax)+(iont.y-betay)*(iont.y-betay));
                //if (fdeltaxy>deltaxy){
                if (!((betax-iont.x>=-deltaxy)&&(betax-iont.x<=deltaxy)&&(betay-iont.y>=-deltaxy)&&(betay-iont.y<=deltaxy))){
                    fimplantMap_it++;
                    continue;
                }
                fimplanti ++;
                //! fill data here!
                hdeltaxy->Fill(betax-iont.x,betay-iont.y);
                double tsdiff=ts-corrts;
                tall=tsdiff/1;
                hdecayc->Fill(tall);
                //tall=r.Exp(0.128/log(2));
                //if (tsdiff>=tdead||tsdiff<0) treeb->Fill();
                //treeb->Fill();

                nflag=0;
                nrealflag=0;
                breal=0;
                if (betaevt==iont.evt){
                    if (betatype>=1&&betatype<=17) hdecaycreal->Fill(tall);
                    if (betatype==1||betatype==4||betatype==6) flagbeta1++;
                    if (betatype>=1&&betatype<=17) hdecaycomponent[betatype-1]->Fill(tall);
                    if (flagrealneutron>0){
                        hdecayc1neureal->Fill(tall);
                    }
                    nrealflag=flagrealneutron;
                    if (flagrealneutron&&ndelayedneutronf>1){
                        hdecayc1neurealin2n->Fill(tall);
                    }
                    breal=1;
                }
                if (betatype==20) hbkggaus->Fill(tall);
                if (betatype==21) hbkguniform->Fill(tall);


                //* calculation for getting random coincidence
                if (ndelayedneutronf>0){
                    hdecaycmorethan0neu->Fill(tall);
                    hdecayc1neu->Fill(tall);
                }
                if (ndelayedneutronb>0){
                    hdecaycmorethan0neub->Fill(tall);
                    hdecayc1neub->Fill(tall);
                    hdecaygt0nbwd->Fill(tall);
                }

                if (ndelayedneutronf==1){
                    nflag=1;
                    hdecayconly1neu->Fill(tall);
                }
                if (ndelayedneutronb==1){
                    hdecayconly1neub->Fill(tall);
                    hdecay1nbwd->Fill(tall);
                }

                if (ndelayedneutronf==2){
                    hdecayconly2neu->Fill(tall);                    
                }
                if (ndelayedneutronb==2){
                    hdecayconly2neub->Fill(tall);
                    hdecay2nbwd->Fill(tall);
                }
                nfwd=ndelayedneutronf;
                nbwd=ndelayedneutronb;
                btype=betatype;
                treeb->Fill();

                mult=ndelayedneutronf;
                treemlh->Fill();
            }
            fimplantMap_it++;
        }
        if (flagbeta1==1) cout<<betaevt<<endl;

        if (fimplanti>0) ncorr++;
        k++;
    }
    hdeltaxy->Write();
    hdecayc->Write();
    hdecaycreal->Write();

    for (Int_t i=0;i<17;i++){
        hdecaycomponent[i]->Write();
    }


    hbkggaus->Write();
    hbkguniform->Write();
    TH1F* hbkgother=(TH1F*) hdecayc->Clone("hbkgother");
    hbkgother->Add(hdecaycreal,-1);
    hbkgother->Write();

    hneubeam->Write();
    hneubeta->Write();
    hneubeamreal->Write();
    hneubetareal->Write();

    hdecayc1neu->Write();
    hdecayc2neu->Write();
    hdecayc1neureal->Write();
    hdecayc1neub->Write();
    hdecayc1neuin2neu->Write();

    hdecayconly1neu->Write();
    hdecayconly1neub->Write();

    hdecayconly2neu->Write();
    hdecayconly2neub->Write();


    hdecayc2neureal->Write();
    hdecayc2neub->Write();

    hdecayc1neureal3->Write();
    hdecayc1neureal4->Write();
    hdecayc1neureal8->Write();
    hdecayc1neureal10->Write();
    hdecayc1neureal14->Write();
    hdecayc2neuin1neu->Write();
    hdecayc1neurealin2n->Write();
    hdecayc2neuin1neuin2n->Write();


    hdecaycmorethan0neu->Write();
    hdecaycmorethan0neub->Write();
    hdecaycmorethan1neu->Write();
    hdecaycmorethan1neub->Write();

    hdecay1nbwd->Write();
    hdecaygt0nbwd->Write();
    hdecay2nbwd->Write();


    hneuneu->Write();
    treeptmp->Write();

    betatdiff->Write();



    treeb->Write();
    treep1n->Write();
    treep2n->Write();
    treeptmp->Write();
    treemlh->Write();

    fout->Close();
}
