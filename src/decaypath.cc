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
/// \file decaypath.cc
/// \brief Implementation of the decaypath class

#include "decaypath.hh"
#include <iostream>
#include <fstream>


#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TBox.h"

#include "TLine.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TPad.h"
#include "TFrame.h"
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

decaypath::decaypath():
    flistofdecaymember()
{
    fdecaypath=new path();
    finputParms=new char[1000];
    flistofdecaymember.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

decaypath::~decaypath()
{
    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        delete *flistofdecaymember_it;
    }
    flistofdecaymember.clear();
    delete finputParms;
    delete fdecaypath;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::Init(char* inputParms)
{
    sprintf(finputParms,"%s",inputParms);
    std::clog<< __PRETTY_FUNCTION__<<" read input files:"<<
               std::endl<<finputParms<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::ProcessMember(MemberDef *obj)
{
    char temp[100];
    sprintf(temp,"^{%i}",obj->n+obj->z);
    obj->name=obj->name.Prepend(temp);

    //! further processing
    if (obj->decay_hl<0){
        obj->decay_hl=-obj->decay_hl;
        obj->is_decay_hl_fix=0;
    }else{//exclude decay half-life =0;
        obj->is_decay_hl_fix=1;
    }
    if (obj->decay_hlerr==0){//if error ==0-> fix amd no contrain on parameter
        obj->is_decay_hl_fix=2;
    }

    // convert half-life into activity
    obj->decay_lamda=log(2)/obj->decay_hl;
    obj->decay_lamdaerr=log(2)/obj->decay_hl/obj->decay_hl*obj->decay_hlerr;
    obj->decay_lamdalow=log(2)/obj->decay_hlup;
    obj->decay_lamdaup=log(2)/obj->decay_hllow;
    obj->is_decay_lamda_fix=obj->is_decay_hl_fix;

    if (obj->decay_p1n<0){
        obj->decay_p1n=-obj->decay_p1n;
        obj->is_decay_p1n_fix=0;
    }else if (obj->decay_p1n==0){//exclude decay with pn = 0;
        obj->is_decay_p1n_fix=2;
    }else{
        obj->is_decay_p1n_fix=1;
    }

    if (obj->decay_p2n<0){
        obj->decay_p2n=-obj->decay_p2n;
        obj->is_decay_p2n_fix=0;
    }else if (obj->decay_p2n==0){//exclude decay p2n =0;
        obj->is_decay_p2n_fix=2;
    }else{
        obj->is_decay_p2n_fix=1;
    }

    if (obj->neueff<0){
        obj->neueff=-obj->neueff;
        obj->is_neueff_fix=0;
    }else{
        obj->is_neueff_fix=1;
    }
    if (obj->neuefferr==0){//if error ==0-> fix amd no contrain on parameter
        obj->is_neueff_fix=2;
    }

    // convert pn in % to pn in 1
    obj->decay_p1n=obj->decay_p1n/100;
    obj->decay_p1nerr=obj->decay_p1nerr/100;
    obj->decay_p1nlow=obj->decay_p1nlow/100;
    obj->decay_p1nup=obj->decay_p1nup/100;

    obj->decay_p2n=obj->decay_p2n/100;
    obj->decay_p2nerr=obj->decay_p2nerr/100;
    obj->decay_p2nlow=obj->decay_p2nlow/100;
    obj->decay_p2nup=obj->decay_p2nup/100;

    obj->decay_p0n=1-obj->decay_p1n-obj->decay_p2n;
    obj->decay_p0nerr=sqrt(obj->decay_p1nerr*obj->decay_p1nerr+obj->decay_p2nerr*obj->decay_p2nerr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::CopyMember(MemberDef *source, MemberDef *destination)
{
    destination-> id = source->  id;
    destination-> z = source->  z;
    destination-> n = source->  n;
    destination-> name = source->  name;

    destination-> decay_hl = source->  decay_hl;
    destination-> decay_lamda = source->  decay_lamda;

    destination-> decay_p0n = source->  decay_p0n;
    destination-> decay_p1n = source->  decay_p1n;
    destination-> decay_p2n = source->  decay_p2n;

    destination-> decay_hlerr = source->  decay_hlerr;
    destination-> decay_lamdaerr = source->  decay_lamdaerr;
    destination-> decay_p0nerr = source->  decay_p0nerr;
    destination-> decay_p1nerr = source->  decay_p1nerr;
    destination-> decay_p2nerr = source->  decay_p2nerr;

    destination-> decay_hlup = source->  decay_hlup;
    destination-> decay_lamdaup = source->  decay_lamdaup;
    destination-> decay_p0nup = source->  decay_p0nup;
    destination-> decay_p1nup = source->  decay_p1nup;
    destination-> decay_p2nup = source->  decay_p2nup;

    destination-> decay_hllow = source->  decay_hllow;
    destination-> decay_lamdalow = source->  decay_lamdalow;
    destination-> decay_p0nlow = source->  decay_p0nlow;
    destination-> decay_p1nlow = source->  decay_p1nlow;
    destination-> decay_p2nlow = source->  decay_p2nlow;

    destination-> population_ratio = source->  population_ratio;
    destination-> population_ratioerr = source->  population_ratioerr;
    destination-> population_ratioup = source->  population_ratioup;
    destination-> population_ratiolow = source->  population_ratiolow;

    destination-> neueff = source->  neueff;
    destination-> neuefferr = source->  neuefferr;
    destination-> neueffup = source->  neueffup;
    destination-> neuefflow = source->  neuefflow;

    destination-> is_decay_hl_fix = source->  is_decay_hl_fix;
    destination-> is_decay_lamda_fix = source->  is_decay_lamda_fix;
    destination-> is_decay_p0n_fix = source->  is_decay_p0n_fix;
    destination-> is_decay_p1n_fix = source->  is_decay_p1n_fix;
    destination-> is_decay_p2n_fix = source->  is_decay_p2n_fix;
    destination-> is_neueff_fix = source->  is_neueff_fix;

    destination-> gspatner = source->  gspatner;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::appendvectors(vector< vector<Int_t> >& pathoriginal,vector< vector<Int_t> >& pathnew,Int_t newmember,vector< vector<Int_t> >& nneupathoriginal,vector< vector<Int_t> >& nneupathnew,Int_t nneunewmember)
{

    //! Vectors containing information about the node in the decay network
    // if original vector is empty (parent nuclei)
    if (pathoriginal.size()==0){
        vector<Int_t> row;
        row.push_back(0);
        row.push_back(newmember);
        pathnew.push_back(row);
    }

    // if original vector is not empty
    for (Size_t i=0;i<pathoriginal.size();i++)
    {
        // copy original vector
        vector<Int_t> row;
        for (Size_t j=0;j<pathoriginal[i].size();j++)
        {
            row.push_back(pathoriginal[i][j]);
        }
        // add one more element to the row
        row.push_back(newmember);
        pathnew.push_back(row);
    }

    //! Vectors containing information about how many neutron emitted in a decay
    // if original vector is empty (parent nuclei)
    if (nneupathoriginal.size()==0){
        vector<Int_t> nneurow;
        nneurow.push_back(-1); // no meaning for first row
        nneurow.push_back(nneunewmember);
        nneupathnew.push_back(nneurow);
    }

    // if original vector is not empty
    for (Size_t i=0;i<nneupathoriginal.size();i++)
    {
        // copy original vector
        vector<Int_t> nneurow;
        for (Size_t j=0;j<nneupathoriginal[i].size();j++)
        {
            nneurow.push_back(nneupathoriginal[i][j]);
        }
        // add one more element to the row
        nneurow.push_back(nneunewmember);
        nneupathnew.push_back(nneurow);
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::makePath()
{
    std::clog<< __PRETTY_FUNCTION__<<std::endl;
    //! reading input file and put informations to list of members
    Int_t id=0;
    std::string line;
    std::ifstream infile(finputParms);
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (line[0]=='#') continue;
        Int_t a;
        // decay properies
        MemberDef* obj=new MemberDef();
        obj->id=id;
        //! read info without isomer
//        if (!(iss >> obj->name >> obj->z >> a >> obj->decay_hl >> obj->decay_hlerr >> obj->decay_hllow >> obj->decay_hlup >>
//              obj->decay_p1n >> obj->decay_p1nerr >> obj->decay_p1nlow >> obj->decay_p1nup >>
//              obj->decay_p2n >> obj->decay_p2nerr >> obj->decay_p2nlow >> obj->decay_p2nup)) break;
        obj->population_ratio = 1;
        obj->population_ratioerr = 0;
        obj->population_ratiolow = 0;
        obj->population_ratioup = 2;
        obj->is_population_ratio_fix = 2;

        //! read info with isomer
        if (!(iss >> obj->name >> obj->z >> a >> obj->decay_hl >> obj->decay_hlerr >> obj->decay_hllow >> obj->decay_hlup >>
              obj->decay_p1n >> obj->decay_p1nerr >> obj->decay_p1nlow >> obj->decay_p1nup >>
              obj->decay_p2n >> obj->decay_p2nerr >> obj->decay_p2nlow >> obj->decay_p2nup >>
              obj->neueff >> obj->neuefferr >> obj->neuefflow >> obj->neueffup)) break;
        obj->n=a-obj->z;

        obj->gspatner = -1;

        //! isomeric state
        Bool_t flagisomer = false;
        std::string tempstring(obj->name.Data());
        MemberDef* objisomer;
        if (tempstring.back()=='*') {
            objisomer=new MemberDef();
            CopyMember(obj,objisomer);
            objisomer->id = obj->id + 1;
            cout<<"ISOMER of "<<obj->name<<endl;
            if (!(iss >> objisomer->population_ratio >> objisomer->population_ratioerr >> objisomer->population_ratiolow >> objisomer->population_ratioup >> objisomer->decay_hl >> objisomer->decay_hlerr >> objisomer->decay_hllow >> objisomer->decay_hlup >>
                  objisomer->decay_p1n >> objisomer->decay_p1nerr >> objisomer->decay_p1nlow >> objisomer->decay_p1nup >>
                  objisomer->decay_p2n >> objisomer->decay_p2nerr >> objisomer->decay_p2nlow >> objisomer->decay_p2nup >>
                  objisomer->neueff >> objisomer->neuefferr >> objisomer->neuefflow >> objisomer->neueffup)) break;
            obj->name = TString(tempstring.substr(0,tempstring.length()-1).data());

            if (objisomer->population_ratio>0){
                obj->is_population_ratio_fix = 1;
                objisomer->is_population_ratio_fix = 1;
            }else{
                objisomer->population_ratio = -objisomer->population_ratio;
                obj->is_population_ratio_fix = 0;
                objisomer->is_population_ratio_fix = 0;
            }

            obj->population_ratio = 1 - objisomer->population_ratio;
            obj->population_ratioerr = objisomer->population_ratioerr;
            obj->population_ratioup = 1 - objisomer->population_ratiolow;
            obj->population_ratiolow = 1 - objisomer->population_ratioup;

            obj->gspatner = obj->id;
            flagisomer = true;
        }

        //! ground state
        ProcessMember(obj);
        flistofdecaymember.emplace(flistofdecaymember.end(),obj);
        id++;

        if (flagisomer){
            //! isomeric state
            ProcessMember(objisomer);
            flistofdecaymember.emplace(flistofdecaymember.end(),objisomer);
            id++;
        }
    }

    std::list<MemberDef*>::iterator listofdecaymember_it2;
    //! main loop to make paths
    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        for (listofdecaymember_it2 = flistofdecaymember.begin(); listofdecaymember_it2 != flistofdecaymember.end(); listofdecaymember_it2++)
        {
            if (((*flistofdecaymember_it)->z-(*listofdecaymember_it2)->z)==1){
                if (((*flistofdecaymember_it)->n-(*listofdecaymember_it2)->n)==-1){//p0n
                    appendvectors((*listofdecaymember_it2)->path,(*flistofdecaymember_it)->path,(*flistofdecaymember_it)->id,
                                  (*listofdecaymember_it2)->nneupath,(*flistofdecaymember_it)->nneupath,0);
                }else if (((*flistofdecaymember_it)->n-(*listofdecaymember_it2)->n)==-2){//p1n
                    appendvectors((*listofdecaymember_it2)->path,(*flistofdecaymember_it)->path,(*flistofdecaymember_it)->id,
                                  (*listofdecaymember_it2)->nneupath,(*flistofdecaymember_it)->nneupath,1);
                }else if (((*flistofdecaymember_it)->n-(*listofdecaymember_it2)->n)==-3){//p2n
                    appendvectors((*listofdecaymember_it2)->path,(*flistofdecaymember_it)->path,(*flistofdecaymember_it)->id,
                                  (*listofdecaymember_it2)->nneupath,(*flistofdecaymember_it)->nneupath,2);
                }
            }
        }
    }

    fdecaypath->npaths=0;
    //! passing paths
    // number of ri
    fdecaypath->nri=flistofdecaymember.size();
    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        for (Size_t i=0;i<(*flistofdecaymember_it)->path.size();i++){
            Double_t prevz=0;
            Double_t prevn=0;
            Double_t previd=0;
            Double_t prevp1n=0;
            Double_t prevp2n=0;
            Bool_t isflow=true;
            fdecaypath->ndecay[fdecaypath->npaths]=(*flistofdecaymember_it)->path[i].size();
            for (Size_t j=0;j<(*flistofdecaymember_it)->path[i].size();j++){
                fdecaypath->decaymap[fdecaypath->npaths][(Int_t)j]=(*flistofdecaymember_it)->path[i][j];

                //check if the current member matched
                for (listofdecaymember_it2 = flistofdecaymember.begin(); listofdecaymember_it2 != flistofdecaymember.end(); listofdecaymember_it2++)
                {
                    if ((*flistofdecaymember_it)->path[i][j]==(*listofdecaymember_it2)->id){
                        Double_t presz=(*listofdecaymember_it2)->z;
                        Double_t presn=(*listofdecaymember_it2)->n;
                        Double_t presid=(*listofdecaymember_it2)->id;
                        if (prevz+prevn==presz+presn){
                            if ((1-prevp1n+prevp2n)==0) {
                                isflow=false;
                            }
                        }else if(presz+presn==prevz+prevn-1){
                            if (prevp1n==0) {
                                isflow=false;
                            }
                        }else if((presz+presn==prevz+prevn-2)){
                            if (prevp2n==0) {
                                isflow=false;
                            }
                        }

                        prevz=presz;
                        prevn=presn;
                        prevp1n=(*listofdecaymember_it2)->decay_p1n;
                        prevp2n=(*listofdecaymember_it2)->decay_p2n;
                        previd=presid;
                    }
                }//end of check if the current member matched

            }
            fdecaypath->ispathhasflow[fdecaypath->npaths]=isflow;
            //cout<<"Flow"<<fdecaypath->npaths<<"\t"<<isflow<<endl;

            for (Size_t j=1;j<(*flistofdecaymember_it)->path[i].size();j++){
                fdecaypath->nneu[fdecaypath->npaths][(Int_t)j-1]=(*flistofdecaymember_it)->nneupath[i][j];
            }

            fdecaypath->npaths++;
        }
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::writePath()
{
    std::ofstream pathfile("path.txt");
    pathfile<<fdecaypath->nri<<std::endl;
    pathfile<<fdecaypath->npaths<<std::endl;

    for (int i=0;i<fdecaypath->npaths;i++){
        pathfile<<fdecaypath->ndecay[i]<<std::endl;
        pathfile<<fdecaypath->ispathhasflow[i]<<std::endl;
        for (int j=0;j<fdecaypath->ndecay[i];j++){
            pathfile<<fdecaypath->decaymap[i][j]<<"\t"<<fdecaypath->nneu[i][j]<<std::endl;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::printPath()
{
    std::clog<< __PRETTY_FUNCTION__<<std::endl;
    Int_t npaths=0;

    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        cout<<"********* Go for Isotope "<<(*flistofdecaymember_it)->id<<" ("<<(*flistofdecaymember_it)->name<<")"<<endl;
        for (Size_t i=0;i<(*flistofdecaymember_it)->path.size();i++){
            cout<<"row "<<i<<" = ";
            std::list<MemberDef*>::iterator listofdecaymember_it2;
            for (Size_t j=0;j<(*flistofdecaymember_it)->path[i].size();j++){
                cout<<(*flistofdecaymember_it)->path[i][j]<<" ";
                for (listofdecaymember_it2 = flistofdecaymember.begin(); listofdecaymember_it2 != flistofdecaymember.end(); listofdecaymember_it2++)
                {
                    if ((*flistofdecaymember_it)->path[i][j]==(*listofdecaymember_it2)->id){
                        cout<<"("<<(*listofdecaymember_it2)->name<<")\t to \t";
                    }
                }
            }
            cout<<" end"<<endl;
            npaths++;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::printMember()
{
    std::clog<< __PRETTY_FUNCTION__<<std::endl;
    //! display information of the paths to decay
    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        cout<<(*flistofdecaymember_it)->id<<"\t"<<(*flistofdecaymember_it)->name<<"\t"<<(*flistofdecaymember_it)->z<<"\t"<<
              (*flistofdecaymember_it)->n<<"\t"<<(*flistofdecaymember_it)->n+(*flistofdecaymember_it)->z<<"\t"<<(*flistofdecaymember_it)->decay_hl<<"\t"<<(*flistofdecaymember_it)->decay_lamda<<"\t"<<
              (*flistofdecaymember_it)->decay_p1n<<"\t"<<(*flistofdecaymember_it)->decay_p2n<<"\t"<<(*flistofdecaymember_it)->neueff<<"\t"<<
              (*flistofdecaymember_it)->is_decay_hl_fix<<"\t"<<(*flistofdecaymember_it)->is_decay_p1n_fix<<"\t"<<(*flistofdecaymember_it)->is_decay_p2n_fix<<"\t"<<(*flistofdecaymember_it)->is_neueff_fix<<endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::drawPath(char* outputFileName)
{
    std::clog<< __PRETTY_FUNCTION__<<std::endl;
    Double_t minz=100000;
    Double_t minn=100000;
    Double_t maxz=0;
    Double_t maxn=0;

    Double_t expandZ=3;
    Double_t expandN=3;

    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        if ((*flistofdecaymember_it)->z<minz) minz=(*flistofdecaymember_it)->z;
        if ((*flistofdecaymember_it)->n<minn) minn=(*flistofdecaymember_it)->n;

        if ((*flistofdecaymember_it)->z>maxz) maxz=(*flistofdecaymember_it)->z;
        if ((*flistofdecaymember_it)->n>maxn) maxn=(*flistofdecaymember_it)->n;
    }

    TCanvas* c1=new TCanvas("cc","",900,700) ;
    Double_t xrange[2]={minn-expandN,maxn+expandN};
    Double_t yrange[2]={minz-expandZ,maxz+expandZ};

    Double_t minhalflife=0.0001;//100 ns

    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    gStyle->SetOptStat(0);

    Int_t NumRI;

    Int_t nprot;
    Int_t nneut;
    Int_t nmass;
    Double_t hlval;

    //TCanvas* c1=new TCanvas("c1","",900,700) ;

    NumRI = 5346;

    ifstream fdat;
    fdat.open("FRDM-QRPA12-halflife.txt");

    TH2F *hchart = new TH2F("hist","",185,-0.5,184.5,127,-0.5,126.5);
    for (Int_t i=0; i<NumRI; i++) {
      fdat >> nprot >> nneut >> hlval;
      nmass = nneut+nprot;
      hchart->Fill(nneut,nprot,hlval);
    }

    c1->SetLogz(0);
    hchart->SetTitleSize(0.04);
    hchart->GetXaxis()->SetTitleOffset(1.0);
    hchart->GetYaxis()->SetTitleOffset(1.2);
    hchart->GetYaxis()->CenterTitle();
    hchart->GetXaxis()->SetLabelSize(0.03);
    hchart->GetYaxis()->SetLabelSize(0.03);
    //  hchart->GetXaxis()->SetTitleSize(1.1);
    hchart->GetYaxis()->SetTitle("N_{Proton}");
    hchart->GetXaxis()->SetTitle("N_{Neutron}");

    hchart->GetXaxis()->SetRangeUser(xrange[0],xrange[1]);
    hchart->GetYaxis()->SetRangeUser(yrange[0],yrange[1]);

    hchart->SetMinimum(minhalflife);


    c1->SetLogz();
    hchart->SetLineWidth(10);
    hchart->SetLineColor(1);
    hchart->Draw("COLZ");

    //! draw border of isotopes
    TBox b2;
    b2.SetFillStyle(0);
    b2.SetLineColor(2);
    b2.SetLineWidth(1);
    fdat.seekg(0, ios::beg);
    for (Int_t i=0; i<NumRI; i++) {
      fdat >> nprot >> nneut >> hlval;
      if(nprot>=yrange[0]&&nprot<=yrange[1]&&nneut>=xrange[0]&&nneut<=xrange[1]) b2.DrawBox(nneut-0.5,nprot-0.5,nneut+0.5,nprot+0.5);
    }
    fdat.close();

    //! Drawing magic number
    Double_t dd = 0.5;
    TLine a1;
    //  a1.SetLineWidth(1.5);
    a1.SetLineWidth(3.0);
    a1.SetLineColor(7);

    Int_t magicn[]={8,20,28,50,82,126};

    for (Int_t i=0;i<6;i++){
        a1.DrawLine(magicn[i]-dd,yrange[0]-dd,magicn[i]-dd,yrange[1]+dd); a1.DrawLine(magicn[i]+1-dd,yrange[0]-dd,magicn[i]+1-dd,yrange[1]+dd);
        a1.DrawLine(xrange[0]-dd,magicn[i]-dd,xrange[1]+dd,magicn[i]-dd); a1.DrawLine(xrange[0]-dd,magicn[i]+1-dd,xrange[1]+dd,magicn[i]+1-dd);
    }

    //! Draw stable isotopes
    ifstream fdat7;
    fdat7.open("stable.csv");
    Int_t n7=287;
    Double_t xx7[300];
    Double_t yy7[300];

    for (Int_t i=0; i<n7; i++) {
      fdat7 >> yy7[i] >> xx7[i];
    }
    TGraph *gr7 = new TGraph(n7,xx7,yy7);
    gr7->SetMarkerStyle(21);
    gr7->SetMarkerColor(1);
    gr7->SetMarkerSize(1.8);
    gr7->Draw("PS");

    //! draw paths
    TLatex latex;
    latex.SetTextAlign(12);
    latex.SetTextSize(0.025);

    TArrow arr;
    arr.SetLineColor(6);
    arr.SetFillColor(2);
    //arr.SetLineStyle(2);
    // Plot the flow
    Int_t npaths=0;

    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        Int_t isplotiso=0;
        for (Size_t i=0;i<(*flistofdecaymember_it)->path.size();i++){
            Double_t prevz=0;
            Double_t prevn=0;
            Double_t previd=0;
            Double_t prevp1n=0;
            Double_t prevp2n=0;
            Bool_t isplot=true;

            std::list<MemberDef*>::iterator listofdecaymember_it2;
            for (Size_t j=0;j<(*flistofdecaymember_it)->path[i].size();j++){
                //! draw arrow
                for (listofdecaymember_it2 = flistofdecaymember.begin(); listofdecaymember_it2 != flistofdecaymember.end(); listofdecaymember_it2++)
                {
                    if ((*flistofdecaymember_it)->path[i][j]==(*listofdecaymember_it2)->id){

                        Double_t presz=(*listofdecaymember_it2)->z;
                        Double_t presn=(*listofdecaymember_it2)->n;
                        Double_t presid=(*listofdecaymember_it2)->id;


                        if (prevz+prevn==presz+presn){
                            if ((1-prevp1n+prevp2n)==0) {
                                isplot=false;
                            }
                        }else if(presz+presn==prevz+prevn-1){
                            if (prevp1n==0) {
                                isplot=false;
                            }
                        }else if((presz+presn==prevz+prevn-2)){
                            if (prevp2n==0) {
                                isplot=false;
                            }
                        }

                        if (prevz!=0&&isplot){
                            arr.DrawArrow(prevn,prevz,presn,presz,0.01,">");
                            //! A trick for plotting, draw at the end of each track another arrow
                            arr.DrawArrow(presn,presz,presn-1,presz+1,0.01,">");
                            isplotiso++;
                            //latex.DrawLatex((*listofdecaymember_it2)->n-0.5,(*listofdecaymember_it2)->z,Form("%s",(*listofdecaymember_it2)->name.Data()));
                        }

                        prevz=presz;
                        prevn=presn;
                        prevp1n=(*listofdecaymember_it2)->decay_p1n;
                        prevp2n=(*listofdecaymember_it2)->decay_p2n;
                        previd=presid;
                    }
                }

            }
            npaths++;
        }
    }

    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        if ((*flistofdecaymember_it)->id==0) latex.DrawLatex((*flistofdecaymember_it)->n-0.5,(*flistofdecaymember_it)->z,Form("%s",(*flistofdecaymember_it)->name.Data()));
        else latex.DrawLatex((*flistofdecaymember_it)->n-0.5,(*flistofdecaymember_it)->z,Form("%s",(*flistofdecaymember_it)->name.Data()));
    }

    //! save to file
    c1->SaveAs(outputFileName);
}

