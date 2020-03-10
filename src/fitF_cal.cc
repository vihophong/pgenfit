#include <fitF.hh>
#include <RooAbsReal.h>
#include <RooAbsCategory.h>
#include <math.h>
#include <TMath.h>
#include <TStopwatch.h>

Double_t fitF::fcndecay(Double_t *x, Double_t *par) const
{
    Double_t t=x[0];
    Double_t* l=&par[0];
    Double_t e[fpath->nri];
    for (Int_t i=0;i<fpath->nri;i++) e[i]=TMath::Exp(-l[i]*t);

    Double_t* p1n=&par[fpath->nri];
    Double_t* p2n=&par[fpath->nri2];
    Double_t* py=&par[fpath->nri3];
    Double_t* ne=&par[fpath->nri4];
    Double_t N0=par[fpath->nri5]/l[0];
    Double_t fparent=0;
    Double_t fdecayall=0;
    Double_t fdaugters[fpath->npaths];
    calculateDecay(fdecayall,fparent, fdaugters, l,e,p1n,p2n,py,N0);
    return fdecayall+par[fpath->nri5+4]+par[fpath->nri5+5]*t;
}
Double_t fitF::fcndecay1n(Double_t *x, Double_t *par) const
{
    Double_t t=x[0];
    Double_t* l=&par[0];
    Double_t e[fpath->nri];
    for (Int_t i=0;i<fpath->nri;i++) e[i]=TMath::Exp(-l[i]*t);
    Double_t* p1n=&par[fpath->nri];
    Double_t* p2n=&par[fpath->nri2];
    Double_t* py=&par[fpath->nri3];
    Double_t* ne=&par[fpath->nri4];
    Double_t N0=par[fpath->nri5]/l[0];
    Double_t randcoinfgt0n=*p[fpath->nri5+2];
    Double_t randcoinf1n=*p[fpath->nri5+1];
    Double_t fparent=0;
    Double_t fdecayall=0;
    Double_t fdaugters[fpath->npaths];
    calculateDecay(fdecayall,fparent, fdaugters, l,e,p1n,p2n,py,N0);
    //! fdecay1n
    return calculateDecay1n(fparent,fdaugters,p1n,p2n,ne,randcoinf1n,randcoinfgt0n)+par[fpath->nri5+4]+par[fpath->nri5+5]*t;
}
Double_t fitF::fcndecay2n(Double_t *x, Double_t *par) const
{
    Double_t t=x[0];
    Double_t* l=&par[0];
    Double_t e[fpath->nri];
    for (Int_t i=0;i<fpath->nri;i++) e[i]=TMath::Exp(-l[i]*t);
    Double_t* p1n=&par[fpath->nri];
    Double_t* p2n=&par[fpath->nri2];
    Double_t* py=&par[fpath->nri3];
    Double_t* ne=&par[fpath->nri4];
    Double_t N0=par[fpath->nri5]/l[0];
    Double_t randcoinf2n=*p[fpath->nri5+3];
    Double_t randcoinfgt0n=*p[fpath->nri5+2];
    Double_t randcoinf1n=*p[fpath->nri5+1];
    Double_t fparent=0;
    Double_t fdecayall=0;
    Double_t fdaugters[fpath->npaths];
    calculateDecay(fdecayall,fparent, fdaugters, l,e,p1n,p2n,py,N0);
    //! fdecay2n
    return calculateDecay2n(fparent,fdaugters,p1n,p2n,ne,randcoinf1n,randcoinfgt0n,randcoinf2n)+par[fpath->nri5+4]+par[fpath->nri5+5]*t;
}

void fitF::calculateDecay(Double_t &fdecayall,Double_t &fparent, Double_t* fdaugters, Double_t *l,Double_t *e,Double_t *p1n,Double_t *p2n,Double_t *py,Double_t N0) const
{
    fparent=l[0]*N0*e[0];
    fdecayall=fparent;
    for (Int_t i=0;i<fpath->npaths;i++){
#ifdef PATHFLOW
        if (fpath->ispathhasflow[i]){
#endif
            fdaugters[i]=1.;
            for (int j=0;j<fpath->ndecay[i]-1;j++){
                if (fpath->nneu[i][j]==0){
                    fdaugters[i]=fdaugters[i] * py[fpath->decaymap[i][j+1]]*(1-p1n[fpath->decaymap[i][j]]-p2n[fpath->decaymap[i][j]])*l[fpath->decaymap[i][j]];//branching here!
                }else if (fpath->nneu[i][j]==1){
                    fdaugters[i]=fdaugters[i] * py[fpath->decaymap[i][j+1]]*p1n[fpath->decaymap[i][j]]*l[fpath->decaymap[i][j]];
                }else{
                    fdaugters[i]=fdaugters[i] * py[fpath->decaymap[i][j+1]]*p2n[fpath->decaymap[i][j]]*l[fpath->decaymap[i][j]];
                }
            }

            Double_t factor2=0;
            for (int j=0;j<fpath->ndecay[i];j++){
                Double_t factor2jk=1;
                for (int k=0;k<fpath->ndecay[i];k++){
                    if (k!=j) {
                        factor2jk=factor2jk*(l[fpath->decaymap[i][k]]-l[fpath->decaymap[i][j]]);
                    }
                }
                factor2=factor2+e[fpath->decaymap[i][j]]/factor2jk;
            }
            fdaugters[i]=l[fpath->decaymap[i][fpath->ndecay[i]-1]]*fdaugters[i]*N0*factor2;

            fdecayall+=fdaugters[i];
#ifdef PATHFLOW
        }
#endif
    }
}

Double_t fitF::calculateDecay1n(Double_t fparent, Double_t* fdaugters, Double_t *p1n,Double_t *p2n,Double_t *ne,Double_t randcoinf1n,Double_t randcoinfgt0n) const
{
    //! fdecay1n
    Double_t fdecay1n=fparent*(randcoinf1n+ne[0]*p1n[0]*(1-randcoinf1n-randcoinfgt0n)+2*(ne[0]*(1-ne[0]))*p2n[0]*(1-randcoinf1n-randcoinfgt0n)-ne[0]*ne[0]*p2n[0]*randcoinf1n);
    for (Int_t i=0;i<fpath->npaths;i++){
#ifdef PATHFLOW
        if (fpath->ispathhasflow[i]){
#endif
            fdecay1n+=fdaugters[i]*(randcoinf1n+ne[fpath->decaymap[i][fpath->ndecay[i]-1]]*p1n[fpath->decaymap[i][fpath->ndecay[i]-1]]*(1-randcoinf1n-randcoinfgt0n)+2*(ne[fpath->decaymap[i][fpath->ndecay[i]-1]]*(1-ne[fpath->decaymap[i][fpath->ndecay[i]-1]]))*p2n[fpath->decaymap[i][fpath->ndecay[i]-1]]*(1-randcoinf1n-randcoinfgt0n)-ne[fpath->decaymap[i][fpath->ndecay[i]-1]]*ne[fpath->decaymap[i][fpath->ndecay[i]-1]]*p2n[fpath->decaymap[i][fpath->ndecay[i]-1]]*randcoinf1n);
#ifdef PATHFLOW
        }
#endif
    }
    return fdecay1n;
}
Double_t fitF::calculateDecay2n(Double_t fparent, Double_t* fdaugters, Double_t *p1n,Double_t *p2n,Double_t *ne,Double_t randcoinf1n,Double_t randcoinfgt0n,Double_t randcoinf2n) const
{
    //! fdecay2n
    Double_t fdecay2n=fparent*(ne[0]*ne[0]*p2n[0]*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne[0]*p1n[0]*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne[0]*(1-ne[0]))*p2n[0]*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
    for (Int_t i=0;i<fpath->npaths;i++){
#ifdef PATHFLOW
        if (fpath->ispathhasflow[i]){
#endif
            fdecay2n+=fdaugters[i]*(ne[fpath->decaymap[i][fpath->ndecay[i]-1]]*ne[fpath->decaymap[i][fpath->ndecay[i]-1]]*p2n[fpath->decaymap[i][fpath->ndecay[i]-1]]*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne[fpath->decaymap[i][fpath->ndecay[i]-1]]*p1n[fpath->decaymap[i][fpath->ndecay[i]-1]]*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne[fpath->decaymap[i][fpath->ndecay[i]-1]]*(1-ne[fpath->decaymap[i][fpath->ndecay[i]-1]]))*p2n[fpath->decaymap[i][fpath->ndecay[i]-1]]*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
#ifdef PATHFLOW
        }
#endif
    }
    return fdecay2n;
}
