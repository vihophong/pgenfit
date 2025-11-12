// autogen_with_alpha_full.cpp
// Generates parameter files including ?/? branches (with uncertainties) and explores decay graph.
//
// Build example:
//   g++ -O2 -std=c++17 autogen_with_alpha_full.cpp `root-config --cflags --libs` -o autogen_with_alpha_full
//
// Run example:
//   ./autogen_with_alpha_full params.txt 138 50
#define P2NVARY 0
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TLine.h>
#include <TArrow.h>
#include <TList.h>
#include <Math/WrappedMultiTF1.h>
#include <HFitInterface.h>
#include <TRandom3.h>
#include <TCutG.h>
#include <TLatex.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <Math/GSLMinimizer.h>
#include <Math/Functor.h>
#include <Math/GSLSimAnMinimizer.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <unistd.h>
#include <map>
#include <list>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <functional>

using namespace std;

// =================== Config switches ===================
// If true: write P1n/P2n as ABSOLUTE fractions of all decays (?% × conditional Pxn / 100).
// If false: write P1n/P2n as CONDITIONAL on ? (legacy behavior).
static const bool WRITE_ABSOLUTE_PXN = true;

// =================== Global factors (from your original) ===================
double deadtime_corr = 1 - 0.0320731;
double neueff      = 0.668 * (1 - 0.0320731); // nominal efficiency (kept for compatibility)
double neuefferr   = 0.02  * (1 - 0.0320731);

// =================== Helpers ===================
string extractIntegerWords(string str)
{
  stringstream ss; ss << str;
  string temp; int found; string temp2="";
  while (!ss.eof()) {
    ss >> temp;
    if (stringstream(temp) >> found) {
      char tempchars[100];
      sprintf(tempchars,"%i",found);
      temp2 += string(tempchars);
    }
    temp = "";
  }
  return temp2;
}
string extractElement(string str)
{
  string tempint = extractIntegerWords(str);
  return str.substr(tempint.length(), str.length() - tempint.length());
}
string converttolatex(string str){
  string tempint = extractIntegerWords(str);
  string tempele = extractElement(str);
  return "^{" + tempint + "}" + tempele;
}

// =================== Data structure (?/? with uncertainties) ===================
typedef struct {
    // identification
    Int_t id;
    Int_t z;
    Int_t n;
    Int_t a;
    TString name;

    // decay properties (?-delayed neutrons; conditional on ?)
    Double_t decay_hl;
    Double_t decay_p0n;   // computed as 100 - p1n - p2n - p3n (conditional on ?)
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

    // efficiencies
    Double_t decay_neueff;
    Double_t decay_neuefferr;
    Double_t decay_neuefferrHi;

    Double_t decay_2neueff;
    Double_t decay_2neuefferr;
    Double_t decay_2neuefferrHi;

    // absolute ? and ? branching ratios (percent of all decays) with uncertainties
    Double_t decay_bbr;        // beta branch (%)
    Double_t decay_bbrerr;     // D_b_br (|?| or 1?)
    Double_t decay_bbrerrhi;   // D_b_br_Hi (asymmetric +)
    Double_t decay_abr;        // alpha branch (%)
    Double_t decay_abrerr;     // D_a_br
    Double_t decay_abrerrhi;   // D_a_br_Hi

    Int_t flag; // 0 = not processed, 1 = processed
} MemberDef;

void CopyMember(MemberDef* s, MemberDef* d)
{
    d->id = s->id; d->a = s->a; d->z = s->z; d->n = s->n; d->name = s->name;

    d->decay_hl = s->decay_hl;

    d->decay_p0n = s->decay_p0n;
    d->decay_p1n = s->decay_p1n;
    d->decay_p2n = s->decay_p2n;
    d->decay_p3n = s->decay_p3n;

    d->decay_hlerr   = s->decay_hlerr;
    d->decay_p0nerr  = s->decay_p0nerr;
    d->decay_p1nerr  = s->decay_p1nerr;
    d->decay_p2nerr  = s->decay_p2nerr;
    d->decay_p3nerr  = s->decay_p3nerr;

    d->decay_hlerrhi   = s->decay_hlerrhi;
    d->decay_p1nerrhi  = s->decay_p1nerrhi;
    d->decay_p2nerrhi  = s->decay_p2nerrhi;
    d->decay_p3nerrhi  = s->decay_p3nerrhi;

    d->decay_neueff     = s->decay_neueff;
    d->decay_neuefferr  = s->decay_neuefferr;
    d->decay_neuefferrHi= s->decay_neuefferrHi;

    d->decay_2neueff     = s->decay_2neueff;
    d->decay_2neuefferr  = s->decay_2neuefferr;
    d->decay_2neuefferrHi= s->decay_2neuefferrHi;

    d->decay_bbr      = s->decay_bbr;
    d->decay_bbrerr   = s->decay_bbrerr;
    d->decay_bbrerrhi = s->decay_bbrerrhi;
    d->decay_abr      = s->decay_abr;
    d->decay_abrerr   = s->decay_abrerr;
    d->decay_abrerrhi = s->decay_abrerrhi;

    d->flag = s->flag;
}

// =================== Plot helper (unchanged from your base) ===================
void plainPlot(TCanvas* c1, Double_t xrange[], Double_t yrange[])
{
  Double_t minhalflife=0.0001;//100 ns

  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(12);
  gStyle->SetOptStat(0);

  Int_t NumRI = 5346;

  ifstream fdat("FRDM-QRPA12-halflife.txt");
  cout << "Get data" << endl;

  TH2F *hchart = new TH2F("hist","",185,-0.5,184.5,127,-0.5,126.5);

  Int_t nprot, nneut, nmass;
  Double_t hlval;
  for (Int_t i=0; i<NumRI; i++) {
    fdat >> nprot >> nneut >> hlval;
    nmass = nneut + nprot;
    hchart->Fill(nneut, nprot, hlval);
  }

  c1->SetLogz(0);
  hchart->SetTitleSize(0.04);
  hchart->GetXaxis()->SetTitleOffset(1.0);
  hchart->GetYaxis()->SetTitleOffset(1.2);
  hchart->GetYaxis()->CenterTitle();
  hchart->GetXaxis()->SetLabelSize(0.03);
  hchart->GetYaxis()->SetLabelSize(0.03);
  hchart->GetYaxis()->SetTitle("N_{Proton}");
  hchart->GetXaxis()->SetTitle("N_{Neutron}");

  hchart->GetXaxis()->SetRangeUser(xrange[0],xrange[1]);
  hchart->GetYaxis()->SetRangeUser(yrange[0],yrange[1]);

  hchart->SetMinimum(minhalflife);

  c1->SetLogz();
  hchart->SetLineWidth(10);
  hchart->SetLineColor(1);
  hchart->Draw("COLZ");

  // draw isotope borders
  TBox b2; b2.SetFillStyle(0); b2.SetLineColor(2); b2.SetLineWidth(1);
  fdat.clear(); fdat.seekg(0, ios::beg);
  for (Int_t i=0; i<NumRI; i++) {
    fdat >> nprot >> nneut >> hlval;
    if(nprot>=yrange[0]&&nprot<=yrange[1]&&nneut>=xrange[0]&&nneut<=xrange[1])
      b2.DrawBox(nneut-0.5,nprot-0.5,nneut+0.5,nprot+0.5);
  }
  fdat.close();

  // magic numbers
  Double_t dd = 0.5;
  TLine a1; a1.SetLineWidth(3.0); a1.SetLineColor(7);
  Int_t magicn[]={8,20,28,50,82,126};
  for (Int_t i=0;i<6;i++){
      a1.DrawLine(magicn[i]-dd,yrange[0]-dd,magicn[i]-dd,yrange[1]+dd);
      a1.DrawLine(magicn[i]+1-dd,yrange[0]-dd,magicn[i]+1-dd,yrange[1]+dd);
      a1.DrawLine(xrange[0]-dd,magicn[i]-dd,xrange[1]+dd,magicn[i]-dd);
      a1.DrawLine(xrange[0]-dd,magicn[i]+1-dd,xrange[1]+dd,magicn[i]+1-dd);
  }

  // r-process path
  TLine a0; a0.SetLineWidth(4); a0.SetLineStyle(1); a0.SetLineColor(3);
  Double_t nn, pp1, pp2;
  ifstream rpathfile("r-process_path.txt");
  while (rpathfile.good()){
    rpathfile >> nn >> pp1 >> pp2;
    if (!rpathfile.good()) break;
    if (nn>=xrange[0]&&nn<=xrange[1]){
      Bool_t isplot=true;
      if (pp1<yrange[0]) pp1=yrange[0];
      if (pp1>yrange[1]) isplot=false;
      if (pp2<yrange[0]) isplot=false;
      if (pp2>yrange[1]) pp2=yrange[1];
      if (isplot){
        a0.DrawLine(nn-0.5,pp1-0.5,nn+0.5,pp1-0.5); a0.DrawLine(nn-0.5,pp1-0.5,nn-0.5,pp2+0.5);
        a0.DrawLine(nn-0.5,pp2+0.5,nn+0.5,pp2+0.5); a0.DrawLine(nn+0.5,pp1-0.5,nn+0.5,pp2+0.5);
      }
    }
  }

  // stable nuclei
  ifstream fdat7("stable.csv");
  cout << "Stable data " << endl;
  Int_t n7=287;
  Double_t xx7[300]; Double_t yy7[300];
  for (Int_t i=0; i<n7; i++) {
    fdat7 >> yy7[i] >> xx7[i];
  }
  TGraph *gr7 = new TGraph(n7,xx7,yy7);
  gr7->SetMarkerStyle(21);
  gr7->SetMarkerColor(1);
  gr7->SetMarkerSize(1.8);
  gr7->Draw("PS");
}

// =================== Main generator ===================
void autogenparmsfile(const char* outputfile, Int_t Ainput=138, Int_t Zinput=50)
{
    Int_t Ninput = Ainput - Zinput;

    list<MemberDef*> listofdecaymember;

    // read listofeval.txt (with 6 extra columns for ?/? & uncertainties)
    std::ifstream infile("listofeval.txt");
    if (!infile.is_open()){
        std::cerr << "ERROR: cannot open listofeval.txt\n";
        return;
    }

    std::string line;
    Int_t id=0, ndecaymember=0;
    while (std::getline(infile, line))
    {
        if (line.empty()) continue;
        if (line[0]=='#') continue;

        std::istringstream iss(line);
        MemberDef* obj = new MemberDef();
        obj->flag = 0;
        obj->id   = id;

        Int_t isiso;
        Double_t tempnum;

        // name z a isiso
        if (!(iss >> obj->name >> obj->z >> obj->a >> isiso)) {
            delete obj; break;
        }

        // skip 12 numbers (as in your original)
        for (Int_t i=0; i<12; ++i) {
            if (!(iss >> tempnum)) { tempnum = 0; }
        }

        // read main fields (same order as original)
        if (!(iss >> obj->decay_hl >> obj->decay_hlerr
                  >> obj->decay_p1n >> obj->decay_p1nerr
                  >> obj->decay_p2n >> obj->decay_p2nerr
                  >> obj->decay_p3n >> obj->decay_p3nerr
                  >> obj->decay_neueff >> obj->decay_neuefferr >> obj->decay_neuefferrHi
                  >> obj->decay_2neueff >> obj->decay_2neuefferr >> obj->decay_2neuefferrHi
                  >> obj->decay_hlerrhi >> obj->decay_p1nerrhi >> obj->decay_p2nerrhi))
        {
            delete obj; break;
        }

        // read 6 extra columns: b_br, D_b_br, D_b_br_Hi, a_br, D_a_br, D_a_br_Hi
        if (!(iss >> obj->decay_bbr >> obj->decay_bbrerr >> obj->decay_bbrerrhi
                  >> obj->decay_abr >> obj->decay_abrerr >> obj->decay_abrerrhi))
        {
            // fallback defaults if missing
            obj->decay_bbr = 100.0; obj->decay_bbrerr = 0.0; obj->decay_bbrerrhi = 0.0;
            obj->decay_abr = 0.0;   obj->decay_abrerr = 0.0; obj->decay_abrerrhi = 0.0;
        }

        obj->decay_p0n = 100.0 - obj->decay_p1n - obj->decay_p2n - obj->decay_p3n;
        obj->decay_p0nerr = std::sqrt(
            obj->decay_p1nerr*obj->decay_p1nerr +
            obj->decay_p2nerr*obj->decay_p2nerr +
            obj->decay_p3nerr*obj->decay_p3nerr
        );

        obj->n = obj->a - obj->z;

        if (isiso==0) {
            listofdecaymember.emplace(listofdecaymember.end(), obj);
            cout << "ground state in " << obj->name << " T1/2= " << obj->decay_hl << endl;
            ndecaymember++;
        } else {
            // drop isomers for this graph (as in your original)
            delete obj;
        }
        id++;
    }
    infile.close();

    // // visual (optional)
	// // === automatic axis ranges ===
	// Double_t minZ = 9999, maxZ = -9999;
	// Double_t minN = 9999, maxN = -9999;
	// for (auto it = listofdecaymember.begin(); it != listofdecaymember.end(); ++it) {
	//     if ((*it)->z < minZ) minZ = (*it)->z;
	//     if ((*it)->z > maxZ) maxZ = (*it)->z;
	//     if ((*it)->n < minN) minN = (*it)->n;
	//     if ((*it)->n > maxN) maxN = (*it)->n;
	// }
	// // Add a small padding for better visualization
	// Double_t xrange[2] = {minN - 1, maxN + 1};
	// Double_t yrange[2] = {minZ - 1, maxZ + 1};

	// // === plotting ===
	// TCanvas* cc = new TCanvas("cc","",900,700);
	// plainPlot(cc, xrange, yrange);


    // // latexify labels
    // TLatex latex; latex.SetTextAlign(12); latex.SetTextSize(0.025);
    // for (auto it = listofdecaymember.begin(); it != listofdecaymember.end(); ++it) {
    //     string riname = converttolatex(string((*it)->name.Data()));
    //     (*it)->name = TString(riname.data());
    //     latex.DrawLatex((*it)->n-0.5, (*it)->z, Form("%s", (*it)->name.Data()));
    //     ndecaymember++;
    // }

    // seed available members with the requested (Z,N)
    list<MemberDef*> listofavailablemember;
    for (auto it = listofdecaymember.begin(); it != listofdecaymember.end(); ++it) {
        if ((*it)->z==Zinput && (*it)->n==Ninput) {
            MemberDef* obj = new MemberDef();
            CopyMember(*it, obj);
            if (obj->decay_p1n==0) obj->decay_p1n = 50;
            if (obj->decay_p2n==0) obj->decay_p2n = P2NVARY;
            listofavailablemember.emplace(listofavailablemember.end(), obj);
            break;
        }
    }

    // ====== BUILD DAUGHTERS ======
    // Fix: do NOT set flag until ALL possible daughters (? and ?) are tried for that parent.
    while (true) {
        Int_t ndaughter_total = 0;

        for (auto it = listofavailablemember.begin(); it != listofavailablemember.end(); ++it) {
            if ((*it)->flag != 0) continue;

            Int_t ndaughter_this_parent = 0;

            for (auto it2 = listofdecaymember.begin(); it2 != listofdecaymember.end(); ++it2) {
                // ---- ? branches: Z+1, N-1,-2,-3,-4 if beta BR > 0
                if ( ((*it2)->z - (*it)->z) == 1 && (*it)->decay_bbr > 0.0 ) {
                    int dn = (*it2)->n - (*it)->n;

                    if (dn == -1 && (*it)->decay_p0n > 0.0) {        // p0n
                        MemberDef* obj=new MemberDef(); CopyMember(*it2,obj);
                        listofavailablemember.emplace(listofavailablemember.end(),obj);
                        ndaughter_this_parent++; ndaughter_total++;
                    } else if (dn == -2 && (*it)->decay_p1n > 0.0) { // p1n
                        MemberDef* obj=new MemberDef(); CopyMember(*it2,obj);
                        listofavailablemember.emplace(listofavailablemember.end(),obj);
                        ndaughter_this_parent++; ndaughter_total++;
                    } else if (dn == -3 && (*it)->decay_p2n > 0.0) { // p2n
                        MemberDef* obj=new MemberDef(); CopyMember(*it2,obj);
                        listofavailablemember.emplace(listofavailablemember.end(),obj);
                        ndaughter_this_parent++; ndaughter_total++;
                    } else if (dn == -4 && (*it)->decay_p3n > 0.0) { // p3n
                        MemberDef* obj=new MemberDef(); CopyMember(*it2,obj);
                        listofavailablemember.emplace(listofavailablemember.end(),obj);
                        ndaughter_this_parent++; ndaughter_total++;
                    }
                }

                // ---- ? branch: Z-2, N-2 if alpha BR > 0
                if ( ((*it2)->z - (*it)->z) == -2
                  && ((*it2)->n - (*it)->n) == -2
                  && (*it)->decay_abr > 0.0 )
                {
                    MemberDef* obj=new MemberDef(); CopyMember(*it2,obj);
                    listofavailablemember.emplace(listofavailablemember.end(),obj);
                    ndaughter_this_parent++; ndaughter_total++;
                }
            }

            // mark parent processed only AFTER exploring all possible daughters
            if (ndaughter_this_parent > 0) {
                (*it)->flag = 1;
            } else {
                // Optional: warn if nothing matched (could indicate missing daughter in list)
                // std::cerr << "No daughters found for " << (*it)->name
                //           << " (Z=" << (*it)->z << ", N=" << (*it)->n << ")\n";
            }
        }

        if (ndaughter_total == 0) break;
    }

    // ====== SORT FOUND SET ======
    // sort by increasing Z and decreasing A within each Z (your two-step scheme)
    map<int, MemberDef*> mapA; // key: z*200 + n
    for (auto it = listofavailablemember.begin(); it != listofavailablemember.end(); ++it) {
        mapA.insert(make_pair((*it)->z*200 + (*it)->n, *it));
    }

    map<int, MemberDef*, greater<int> > mapAsort; // A descending per Z
    list<MemberDef*> listofavailablemembersorted;
    Int_t prevZ = 0;
    Int_t idd   = 0;

    for (auto it = mapA.begin(); it != mapA.end(); ++it)
    {
        if ((it->second)->z == prevZ || prevZ == 0) {
            mapAsort.insert(make_pair((it->second)->a, it->second));
        }
        if ((it->second)->z != prevZ && prevZ != 0) {
            for (auto it2 = mapAsort.begin(); it2 != mapAsort.end(); ++it2){
                listofavailablemembersorted.emplace(listofavailablemembersorted.end(), it2->second);
            }
            mapAsort.clear();
            mapAsort.insert(make_pair((it->second)->a, it->second));
        }
        if (idd == (int)mapA.size()-1) {
            for (auto it2 = mapAsort.begin(); it2 != mapAsort.end(); ++it2){
                listofavailablemembersorted.emplace(listofavailablemembersorted.end(), it2->second);
            }
            mapAsort.clear();
        }
        idd++;
        prevZ = (it->second)->z;
    }

    // ====== WRITE OUTPUTS ======
    std::ofstream str(outputfile);
    if (!str.is_open()){
        std::cerr << "ERROR: cannot open " << outputfile << " for writing\n";
        return;
    }

    str << "# Note: start comments with '#'. Negative values indicate variation seeds; half-life in seconds.\n";
    str << "# RI rows sorted by increasing Z and decreasing A within each Z.\n";
    str << "#Name\tZ\tA\tHalf-life\tAbs_Error_HL\tAbs_Error_HL_Hi\tlowerHL\tupperHL\t"
        "P1n\tAbs_Error_P1n\tAbs_Error_P1n_Hi\tlowerP1n\tupperP1n\t"
        "P2n\tAbs_Error_P2n\tAbs_Error_P2n_Hi\tlowerP2n\tupperP2n\t"
        "Neu.Eff\tNeu.Eff_Err\tNeu.Eff_ErrHi\tlowerNeu.Eff\tupperNeu.Eff\t"
        "AlphaBR\tD_AlphaBR\tD_AlphaBR_Hi\tlowerAlphaBR\tupperAlphaBR\n";

    auto absP = [&](const MemberDef* o, double p)->double{
        return WRITE_ABSOLUTE_PXN ? (p * o->decay_bbr / 100.0) : p;
    };
// identify the seed we actually requested
const int Zseed = Zinput;
const int Nseed = Ainput - Zinput;

for (auto it = listofavailablemembersorted.begin(); it != listofavailablemembersorted.end(); ++it)
{
    MemberDef* obj = *it;

    double P1 = absP(obj, obj->decay_p1n);
    double P2 = absP(obj, obj->decay_p2n);

    const bool isSeedLine = (obj->z == Zseed && obj->n == Nseed);

    if (isSeedLine) {
        // Seed line: NEGATIVE values + optional trailing token
        str << obj->name << "\t" << obj->z << "\t" << obj->a << "\t"
            << -obj->decay_hl << "\t" << obj->decay_hlerr << "\t" << obj->decay_hlerrhi << "\t"
            << obj->decay_hl/5 << "\t" << obj->decay_hl*5 << "\t"
            << -P1 << "\t" << obj->decay_p1nerr << "\t" << obj->decay_p1nerrhi << "\t"
            << 0. << "\t" << 200. << "\t"
            << ((P2>0) ? -P2 : 0.) << "\t" << obj->decay_p2nerr << "\t" << obj->decay_p2nerrhi << "\t"
            << 0. << "\t" << 200. << "\t"
            << obj->decay_neueff*deadtime_corr << "\t"
            << obj->decay_neuefferr*deadtime_corr << "\t"
            << obj->decay_neuefferrHi*deadtime_corr << "\t"
            << 0. << "\t" << 1. << "\t"
            << obj->decay_abr << "\t" << obj->decay_abrerr << "\t" << obj->decay_abrerrhi << "\t"
            << 0. << "\t" << 200.
            << "\tPARENT"   // optional: helps downstream makePath() pick the true parent
            << "\n";
    } else {
        // Normal lines: POSITIVE values
        str << obj->name << "\t" << obj->z << "\t" << obj->a << "\t"
            << obj->decay_hl << "\t" << obj->decay_hlerr << "\t" << obj->decay_hlerrhi << "\t"
            << obj->decay_hl/5 << "\t" << obj->decay_hl*5 << "\t"
            << P1 << "\t" << obj->decay_p1nerr << "\t" << obj->decay_p1nerrhi << "\t"
            << 0. << "\t" << 200. << "\t"
            << P2 << "\t" << obj->decay_p2nerr << "\t" << obj->decay_p2nerrhi << "\t"
            << 0. << "\t" << 200. << "\t"
            << obj->decay_neueff*deadtime_corr << "\t"
            << obj->decay_neuefferr*deadtime_corr << "\t"
            << obj->decay_neuefferrHi*deadtime_corr << "\t"
            << 0. << "\t" << 1. << "\t"
            << obj->decay_abr << "\t" << obj->decay_abrerr << "\t" << obj->decay_abrerrhi << "\t"
            << 0. << "\t" << 200.
            << "\n";
    }
}
    str.close();

    // efficiency parameter file
    char tmp[512];
    sprintf(tmp, "%s_effparms", outputfile);
    std::ofstream str2(tmp);
    if (!str2.is_open()){
        std::cerr << "ERROR: cannot open " << tmp << " for writing\n";
        return;
    }
    str2 << "# set value negative: vary from val-err to val+err; err=0 => constant\n";
    str2 << "# betaEffFactor\terr\tbetaEffFactor_1n\terr\tbetaEffFactor_2n\terr\tneutronEffFactor_1nvs2n\terrLo\terrHi\n";
    if (!listofavailablemembersorted.empty()){
        MemberDef* obj = listofavailablemembersorted.front();
        str2 << "1\t0\t1\t0\t1\t0\t"
             << obj->decay_2neueff*deadtime_corr << "\t"
             << obj->decay_2neuefferr*deadtime_corr << "\t"
             << obj->decay_2neuefferrHi*deadtime_corr << "\n";
    }
    str2.close();

    cout << "Written: " << outputfile << " and " << tmp << endl;
    // cc->SaveAs("test.root");
}

// =================== CLI main ===================
int main(int argc, char** argv)
{
    if (argc < 2 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <output_file> [Ainput] [Zinput]\n";
        std::cerr << "  Example: " << argv[0] << " params.txt 138 50\n";
        return 1;
    }
    const char* out = argv[1];
    int A = 138, Z = 50;
    if (argc >= 3) A = atoi(argv[2]);
    if (argc >= 4) Z = atoi(argv[3]);

    autogenparmsfile(out, A, Z);
    return 0;
}

