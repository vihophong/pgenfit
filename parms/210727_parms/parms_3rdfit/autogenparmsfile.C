#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <unistd.h>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

#include "TLine.h"
#include "TArrow.h"


#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"

#include "TTree.h"
#include "TRandom3.h"
#include "TFile.h"
#include <fstream>


#include "TCutG.h"
#include "TLatex.h"
#include "TMath.h"
#include "TCanvas.h"

#include "TGraph.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/GSLSimAnMinimizer.h"
#include "Math/Functor.h"

#include <list>

#include <sstream>
#include <string>


using namespace std; 


double deadtime_corr = 1-0.0320731;
double neueff=0.668*(1-0.0320731);//nominal efficiency
double neuefferr=0.02*(1-0.0320731);

string extractIntegerWords(string str) 
{ 
  stringstream ss; 

  /* Storing the whole string into string stream */
  ss << str; 

  /* Running loop till the end of the stream */
  string temp; 
  int found; 
  string temp2="";
  while (!ss.eof()) { 

    /* extracting word by word from stream */
    ss >> temp; 
    /* Checking the given word is integer or not */
    
    if (stringstream(temp) >> found) {
      char tempchars[100];
      sprintf(tempchars,"%i",found);
      temp2+=string(tempchars);
    }      
    /* To save from space at the end of string */
    temp = "";
  }
  return temp2;
} 
string extractElement(string str)
{
  string tempint=extractIntegerWords(str);
  return str.substr(tempint.length(),str.length()-tempint.length());
}
string converttolatex(string str){
  string tempint=extractIntegerWords(str);
  string tempele=extractElement(str);
  string ret="^{"+tempint+"}"+tempele;
    return ret;
}
void checkconverttolatex() 
{ 
  string str = "134In"; 
  string ee=converttolatex(str);
  cout<<ee<<endl;
} 
typedef struct {
    // idendification
    Int_t id;
    Int_t z;
    Int_t n;
    Int_t a;
    TString name;

    // decay properies
    Double_t decay_hl;
    Double_t decay_p0n;//decay to several isomerics states or ground state
    Double_t decay_p1n;
    Double_t decay_p2n;
    Double_t decay_p3n;

    Double_t decay_hlerr;
    Double_t decay_p0nerr;
    Double_t decay_p1nerr;
    Double_t decay_p2nerr;
    Double_t decay_p3nerr;

    Double_t decay_hlerrhi;
    Double_t decay_p1nerrhi;
    Double_t decay_p2nerrhi;
    Double_t decay_p3nerrhi;

    Double_t decay_neueff;
    Double_t decay_neuefferr;
    Double_t decay_neuefferrHi;

    Double_t decay_2neueff;
    Double_t decay_2neuefferr;
    Double_t decay_2neuefferrHi;

    Int_t flag;
} MemberDef;
void CopyMember(MemberDef* source, MemberDef* destination)
{
    destination-> id = source->  id;
    destination-> a = source->  a;
    destination-> z = source->  z;
    destination-> n = source->  n;
    destination-> name = source->  name;

    destination-> decay_hl = source->  decay_hl;

    destination-> decay_p0n = source->  decay_p0n;
    destination-> decay_p1n = source->  decay_p1n;
    destination-> decay_p2n = source->  decay_p2n;
    destination-> decay_p3n = source->  decay_p3n;

    destination-> decay_hlerr = source->  decay_hlerr;
    destination-> decay_p0nerr = source->  decay_p0nerr;
    destination-> decay_p1nerr = source->  decay_p1nerr;
    destination-> decay_p2nerr = source->  decay_p2nerr;
    destination-> decay_p3nerr = source->  decay_p3nerr;

    destination-> decay_hlerrhi = source->  decay_hlerrhi;
    destination-> decay_p1nerrhi = source->  decay_p1nerrhi;
    destination-> decay_p2nerrhi = source->  decay_p2nerrhi;
    destination-> decay_p3nerrhi = source->  decay_p3nerrhi;

    destination-> decay_neueff = source->  decay_neueff;
    destination-> decay_neuefferr = source->  decay_neuefferr;
    destination-> decay_neuefferrHi = source->  decay_neuefferrHi;

    destination-> decay_2neueff = source->  decay_2neueff;
    destination-> decay_2neuefferr = source->  decay_2neuefferr;
    destination-> decay_2neuefferrHi = source->  decay_2neuefferrHi;

    destination-> flag = source->  flag;
}
void plainPlot(TCanvas* c1,Double_t xrange[],Double_t yrange[])
{

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
  //  NumRI = 24;

  //  ofstream fout2;
  //  fout2.open("zzz.dat");

  ifstream fdat;
  fdat.open("FRDM-QRPA12-halflife.txt");

  cout << "Get data" << endl;

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

  //! Drawing r-process path
  TLine a0;
  a0.SetLineWidth(4);
  a0.SetLineStyle(1);
  a0.SetLineColor(3);
  Double_t nn;
  Double_t pp1;
  Double_t pp2;
  ifstream rpathfile("r-process_path.txt");
   while (rpathfile.good()){
       rpathfile>>nn>>pp1>>pp2;
       if (nn>=xrange[0]&&nn<=xrange[1]){
           Bool_t isplot=true;
           if (pp1<yrange[0]) pp1=yrange[0];
           if (pp1>yrange[1]) isplot=false;
           if (pp2<yrange[0]) isplot=false;
           if (pp2>yrange[1]) pp2=yrange[1];
           if (isplot){
            a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5);a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5);
            a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5);a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);
           }
       }
   }

  ifstream fdat7;
  fdat7.open("stable.csv");
  cout << "Stable data " << endl;
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
}

void autogenparmsfile(char* outputfile, Int_t Ainput=138, Int_t Zinput=50)
{
    Int_t Ninput = Ainput-Zinput;

    list<MemberDef*> listofdecaymember;

    //! readthings into arrays
    std::string line;
    std::ifstream infile("listofeval.txt");
    Int_t id=0;
    Int_t ndecaymember=0;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (line[0]=='#') continue;
        // decay properies
        MemberDef* obj=new MemberDef();
        obj->flag=0;
        obj->id=id;
        Int_t isiso;
        Double_t tempnum;
        //! read info with isomer
        if (!(iss >> obj->name >> obj->z >> obj->a >> isiso)) break;
        for (Int_t i=0;i<12;i++){
            if (!(iss >> tempnum)) break;
        }
        if (!(iss >> obj->decay_hl >> obj->decay_hlerr >> obj->decay_p1n >> obj->decay_p1nerr >>
              obj->decay_p2n >> obj->decay_p2nerr >>
              obj->decay_p3n >> obj->decay_p3nerr >> obj->decay_neueff >> obj->decay_neuefferr >> obj->decay_neuefferrHi>>
              obj->decay_2neueff >> obj->decay_2neuefferr >> obj->decay_2neuefferrHi>>
              obj->decay_hlerrhi >> obj->decay_p1nerrhi >>obj->decay_p2nerrhi)) break;

        obj->decay_p0n=100- obj->decay_p1n - obj->decay_p2n - obj->decay_p3n;
        obj->decay_p0nerr = sqrt(obj->decay_p1nerr*obj->decay_p1nerr+obj->decay_p2nerr*obj->decay_p2nerr+obj->decay_p3nerr*obj->decay_p3nerr);

        obj->n=obj->a-obj->z;
        if (isiso==0) {
            listofdecaymember.emplace(listofdecaymember.end(),obj);
            cout<<"ground state in "<<obj->name<<" T1/2= "<<obj->decay_hl<<endl;
            ndecaymember++;
        }
        //else{
            //cout<<"isomeric in "<<obj->name<<" T1/2= "<<obj->decay_hl<<endl;
        //}

        id++;
    }
    TCanvas* cc=new TCanvas("cc","",900,700) ;
    Double_t xrange[2]={70,95};
    Double_t yrange[2]={45,60};
    plainPlot(cc,xrange,yrange);



    //! construct Z scheme

    list<MemberDef*> listofavailablemember;
    list<MemberDef*>::iterator listofdecaymember_it;
    list<MemberDef*>::iterator listofdecaymember_it2;

    TLatex latex;
    latex.SetTextAlign(12);
    latex.SetTextSize(0.025);

    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        string riname=converttolatex(string((*listofdecaymember_it)->name.Data()));
        (*listofdecaymember_it)->name=TString(riname.data());
        latex.DrawLatex((*listofdecaymember_it)->n-0.5,(*listofdecaymember_it)->z,Form("%s",(*listofdecaymember_it)->name.Data()));
        ndecaymember++;
    }
    for (listofdecaymember_it = listofdecaymember.begin(); listofdecaymember_it != listofdecaymember.end(); listofdecaymember_it++)
    {
        if ((*listofdecaymember_it)->z==Zinput&&(*listofdecaymember_it)->n==Ninput) {
            MemberDef* obj=new MemberDef();
            CopyMember(*listofdecaymember_it,obj);
            if (obj->decay_p1n==0) obj->decay_p1n=50;
            if (obj->decay_p2n==0) obj->decay_p2n=25;
            listofavailablemember.emplace(listofavailablemember.end(),obj);
            break;
        }
    }
    //cout<<"***************"<<endl;

    // main loop
    while(1){
        Int_t ndaughter=0;
        for (listofdecaymember_it = listofavailablemember.begin(); listofdecaymember_it != listofavailablemember.end(); listofdecaymember_it++)
        {
            if((*listofdecaymember_it)->flag==0){
                for (listofdecaymember_it2 = listofdecaymember.begin(); listofdecaymember_it2 != listofdecaymember.end(); listofdecaymember_it2++)
                {
                    if (((*listofdecaymember_it2)->z-(*listofdecaymember_it)->z)==1){
                        if (((*listofdecaymember_it2)->n-(*listofdecaymember_it)->n)==-1){//p0n
                            if ((*listofdecaymember_it)->decay_p0n>0){
                                cout<<(*listofdecaymember_it2)->name<<endl;
                                (*listofdecaymember_it)->flag=1;
                                MemberDef* obj=new MemberDef();
                                CopyMember(*listofdecaymember_it2,obj);
                                listofavailablemember.emplace(listofavailablemember.end(),obj);
                                ndaughter++;
                            }
                        }else if (((*listofdecaymember_it2)->n-(*listofdecaymember_it)->n)==-2){//p1n
                            if ((*listofdecaymember_it)->decay_p1n>0){
                                cout<<(*listofdecaymember_it2)->name<<endl;
                                (*listofdecaymember_it)->flag=1;
                                MemberDef* obj=new MemberDef();
                                CopyMember(*listofdecaymember_it2,obj);
                                listofavailablemember.emplace(listofavailablemember.end(),obj);
                                ndaughter++;
                            }
                        }else if (((*listofdecaymember_it2)->n-(*listofdecaymember_it)->n)==-3){//p2n
                            if ((*listofdecaymember_it)->decay_p2n>0){
                                cout<<(*listofdecaymember_it2)->name<<endl;
                                (*listofdecaymember_it)->flag=1;
                                MemberDef* obj=new MemberDef();
                                CopyMember(*listofdecaymember_it2,obj);
                                listofavailablemember.emplace(listofavailablemember.end(),obj);
                                ndaughter++;
                            }
                        }else if (((*listofdecaymember_it2)->n-(*listofdecaymember_it)->n)==-4){//p3n
                            if ((*listofdecaymember_it)->decay_p3n>0){
                                cout<<(*listofdecaymember_it2)->name<<endl;
                                (*listofdecaymember_it)->flag=1;
                                MemberDef* obj=new MemberDef();
                                CopyMember(*listofdecaymember_it2,obj);
                                listofavailablemember.emplace(listofavailablemember.end(),obj);
                                ndaughter++;
                            }
                        }
                    }
                }
            }
        }
        if (ndaughter==0) break;
    }

    cout<<"output"<<endl;
    for (listofdecaymember_it = listofavailablemember.begin(); listofdecaymember_it != listofavailablemember.end(); listofdecaymember_it++)
    {
        cout<<(*listofdecaymember_it)->name<<endl;
    }
    map<int, MemberDef*> mapA;
    map<int, MemberDef*>::iterator mapA_it;
    cout<<"sorting 1st"<<endl;
    for (listofdecaymember_it = listofavailablemember.begin(); listofdecaymember_it != listofavailablemember.end(); listofdecaymember_it++)
    {
        mapA.insert(make_pair((*listofdecaymember_it)->z*200+(*listofdecaymember_it)->n,*listofdecaymember_it));
    }
    for (mapA_it = mapA.begin(); mapA_it != mapA.end(); mapA_it++)
    {
        cout<<(mapA_it->second)->name<<endl;
    }


    map<int, MemberDef*, greater<int> > mapAsort;
    map<int, MemberDef*>::iterator mapA_it2;

    list<MemberDef*> listofavailablemembersorted;
    Int_t prevZ=0;
    Int_t idd=0;
    cout<<"sorting 2nd"<<endl;
    for (mapA_it = mapA.begin(); mapA_it != mapA.end(); mapA_it++)
    {
        /*
        if ((mapA_it->second)->z==prevZ||prevZ==0){
            mapAsort.insert(make_pair((mapA_it->second)->a,mapA_it->second));
            cout<<(mapA_it->second)->name<<"-"<<(mapA_it->second)->z<<"-"<<prevZ<<endl;
        }else{
            //if (prevZ!=0){
                cout<<"*******sep******"<<endl;
                listofavailablemembersorted.emplace(listofavailablemembersorted.end(),mapA_it->second);
                mapAsort.clear();
                mapAsort.insert(make_pair((mapA_it->second)->a,mapA_it->second));
                cout<<(mapA_it->second)->name<<"-"<<(mapA_it->second)->z<<"-"<<prevZ<<endl;
            //}

        }

        if (idd==mapA.size()-1){
            cout<<"*******sep******"<<endl;
        }
        idd++;
        prevZ=(mapA_it->second)->z;
        */

        if ((mapA_it->second)->z==prevZ||prevZ==0){
            mapAsort.insert(make_pair((mapA_it->second)->a,mapA_it->second));
            cout<<(mapA_it->second)->name<<"-"<<(mapA_it->second)->z<<"-"<<prevZ<<endl;
        }
        if ((mapA_it->second)->z!=prevZ&&prevZ!=0){
            cout<<"*******sep******"<<endl;
            for (mapA_it2 = mapAsort.begin(); mapA_it2 != mapAsort.end(); mapA_it2++){
                listofavailablemembersorted.emplace(listofavailablemembersorted.end(),mapA_it2->second);
            }
            mapAsort.clear();
            mapAsort.insert(make_pair((mapA_it->second)->a,mapA_it->second));
            cout<<(mapA_it->second)->name<<"-"<<(mapA_it->second)->z<<"-"<<prevZ<<endl;
        }

        if (idd==mapA.size()-1){
            cout<<"*******sep******"<<endl;
            for (mapA_it2 = mapAsort.begin(); mapA_it2 != mapAsort.end(); mapA_it2++){
                listofavailablemembersorted.emplace(listofavailablemembersorted.end(),mapA_it2->second);
            }
            mapAsort.clear();
        }
        idd++;
        prevZ=(mapA_it->second)->z;
    }
    cout<<"Final"<<endl;

    std::ofstream str(outputfile);
    str<<"# Note should start from #, non-zero upper value indicate varying variable, half-life is in second"<<endl;
    str<<"# RI should be in the sequences of increasing Z and decreasing A for each isotropic row"<<endl;
    str<<"#Name	Z	A	Half-life	Abs_Error_Half-life Abs_Error_Half-life_Hi	lowerHL	upperHL	P1n	Abs_Error_P1n   Abs_Error_P1n_Hi	lowerP1n	upperP1n	P2n	Abs_Error_P2n   Abs_Error_P2n_Hi	lowerP2n	upperP2n	Neu.Eff	Neu.Eff Err	Neu.Eff ErrHi	lowerNeu.Eff	upperNeu.Eff    isomer_ratio	isomer_ratio_err	lower_isomer_ratio	upper_isomer_ratio  ..."<<endl;
    idd=0;
    for (listofdecaymember_it = listofavailablemembersorted.begin(); listofdecaymember_it != listofavailablemembersorted.end(); listofdecaymember_it++)
    {
        MemberDef* obj=*listofdecaymember_it;
        cout<<obj->name<<endl;
        if (idd==0){
            str<<obj->name<<"\t"<<obj->z<<"\t"<<obj->a<<"\t"<<-obj->decay_hl<<"\t"<<obj->decay_hlerr<<"\t"<<obj->decay_hlerrhi<<"\t"<<obj->decay_hl/5<<"\t"<<obj->decay_hl*5<<"\t"<<-obj->decay_p1n<<"\t"<<obj->decay_p1nerr<<"\t"<<obj->decay_p1nerrhi<<"\t"<<0.<<"\t"<<200.<<"\t"<<-obj->decay_p2n<<"\t"<<obj->decay_p2nerr<<"\t"<<obj->decay_p2nerrhi<<"\t"<<0.<<"\t"<<200.<<"\t"<<obj->decay_neueff*deadtime_corr<<"\t"<<obj->decay_neuefferr*deadtime_corr<<"\t"<<obj->decay_neuefferrHi*deadtime_corr<<"\t"<<0.<<"\t"<<1.<<endl;
        }else{
            str<<obj->name<<"\t"<<obj->z<<"\t"<<obj->a<<"\t"<<obj->decay_hl<<"\t"<<obj->decay_hlerr<<"\t"<<obj->decay_hlerrhi<<"\t"<<obj->decay_hl/5<<"\t"<<obj->decay_hl*5<<"\t"<<obj->decay_p1n<<"\t"<<obj->decay_p1nerr<<"\t"<<obj->decay_p1nerrhi<<"\t"<<0.<<"\t"<<200.<<"\t"<<obj->decay_p2n<<"\t"<<obj->decay_p2nerr<<"\t"<<obj->decay_p2nerrhi<<"\t"<<0.<<"\t"<<200.<<"\t"<<obj->decay_neueff*deadtime_corr<<"\t"<<obj->decay_neuefferr*deadtime_corr<<"\t"<<obj->decay_neuefferrHi*deadtime_corr<<"\t"<<0.<<"\t"<<1.<<endl;

        }

        idd++;
    }
    char tmp[500];
    sprintf(tmp,"%s_effparms",outputfile);
    std::ofstream str2(tmp);
    str2<<"#set value to negative: vary from val-err to val+err; set err to 0 mean constant (no MC variation)"<<endl;
    str2<<"#betaEffFactor	err	betaEffFactor_1n	err	betaEffFactor_2n	err	neutronEffFactor_1nvs2n	errLo errHi"<<endl;
    for (listofdecaymember_it = listofavailablemembersorted.begin(); listofdecaymember_it != listofavailablemembersorted.end(); listofdecaymember_it++)
    {
        MemberDef* obj=*listofdecaymember_it;
        str2<<"1	0	1	0	1	0	"<<obj->decay_2neueff*deadtime_corr<<"\t"<<obj->decay_2neuefferr*deadtime_corr<<"\t"<<obj->decay_2neuefferrHi*deadtime_corr<<endl;
        break;
    }
}

