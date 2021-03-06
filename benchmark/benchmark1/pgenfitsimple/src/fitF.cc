/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...


#include "fitF.hh"


ClassImp(fitF)

 fitF::fitF(const char *name, const char *title,
                        RooAbsReal& _x,
                        RooAbsCategory& _y,
                        RooAbsReal *_pp[]) :
   RooAbsPdf(name,title),
   x("x","x",this,_x),
   y("y","y",this,_y)
 {
    std::ifstream pathfile("path.txt");
    Int_t nri;
    pathfile>>nri;
    pathfile.close();
    //std::cout<<nri<<std::endl;
    for (Int_t i=0;i<nri*5+8;i++){
        p[i]=new RooRealProxy(Form("p%i",i),Form("p%i",i),this,*_pp[i]);
    }
 }


 fitF::fitF(const fitF& other, const char* name) :
   RooAbsPdf(other,name),
   x("x",this,other.x),
   y("y",this,other.y)
 {
     initPath();
     for (Int_t i=0;i<fpath->nri5+8;i++)
         p[i]=new RooRealProxy(Form("p%i",i),this,*other.p[i]);
 }

 void fitF::initPath()
 {
     fpath=new path;
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
     pathfile>>fpath->nisomers;
     for (int i=0;i<fpath->nisomers;i++){
         pathfile>>fpath->isomer_gs_index[i]>>fpath->isomer_ex_index[i];
     }
     pathfile.close();
     fpath->nri2=fpath->nri*2;
     fpath->nri3=fpath->nri*3;
     fpath->nri4=fpath->nri*4;
     fpath->nri5=fpath->nri*5;
//     std::cout<<fpath->nri<<std::endl;
//     std::cout<<fpath->npaths<<std::endl;
//     for (int i=0;i<fpath->npaths;i++){
//         std::cout<<fpath->ndecay[i]<<std::endl;
//         for (int j=0;j<fpath->ndecay[i];j++){
//             std::cout<<fpath->decaymap[i][j]<<"\t"<<fpath->nneu[i][j]<<std::endl;
//         }
//     }
 }

#ifndef EVAL_FAST
 Double_t fitF::evaluate() const
 {
      //TStopwatch stopwatch;
      //stopwatch.Start();
      Double_t l[fpath->nri];
      Double_t e[fpath->nri];
      Double_t p1n[fpath->nri];
      Double_t p2n[fpath->nri];
      Double_t py[fpath->nri];
      Double_t ne[fpath->nri];

      Double_t t=x;
      for (Int_t i=0;i<fpath->nri;i++) {
          l[i]=(*p[i]);
          e[i]=TMath::Exp(-l[i]*t);
          p1n[i]=(*p[fpath->nri+i]);
          p2n[i]=(*p[fpath->nri2+i]);
          py[i]=(*p[fpath->nri3+i]);
          ne[i]=(*p[fpath->nri4+i]);
      }
      Double_t N0=*p[fpath->nri5]/ l[0];


      Double_t be=*p[fpath->nri5+4];//factor of parent's beta efficiency, relative to daugters
      Double_t b1ne=*p[fpath->nri5+5];//factor of parent's beta efficiency to 1 neutron emission, relative to daugters
      Double_t b2ne=*p[fpath->nri5+6];//factor of parent's beta efficiency to 2 neutron emission, relative to daugters
      Double_t n1n2ne=*p[fpath->nri5+7];//factor of parent's neutrons efficiency to 2 neutron emission, relative to 1 neutron emission.

      //! replace for more general expression
      Double_t fparent=0;
      Double_t fdecayall=0;
      Double_t fdaugters[fpath->npaths];
      calculateDecay(fdecayall,fparent, fdaugters, l,e,p1n,p2n,py,N0,be);
      Double_t randcoinf2n=*p[fpath->nri5+3];
      Double_t randcoinfgt0n=*p[fpath->nri5+2];
      Double_t randcoinf1n=*p[fpath->nri5+1];

      //! fdecay1n
      Double_t fdecay1n=calculateDecay1n(fparent,fdaugters,p1n,p2n,ne,randcoinf1n,randcoinfgt0n,be,b1ne,b2ne,n1n2ne);

      //! fdecay2n
      Double_t fdecay2n=calculateDecay2n(fparent,fdaugters,p1n,p2n,ne,randcoinf1n,randcoinfgt0n,randcoinf2n,be,b1ne,b2ne,n1n2ne);

      Double_t ret;
      if (y==0){
          ret = fdecayall-fdecay1n-fdecay2n;
      }else if (y==1){
          ret = fdecay1n;
      }else{
          ret = fdecay2n;
      }
      //std::cout<<"EEEEE="<<ret<<std::endl;
      //stopwatch.Stop();
      //std::cout<<"Time (us) = "<<stopwatch.RealTime()*1e6<<std::endl;
      return ret;
 }
#endif

 ClassImp(fitFbkg)

  fitFbkg::fitFbkg(const char *name, const char *title,
                         RooAbsReal& _x,
                         RooAbsCategory& _y,
                         RooAbsReal& _bkg1,
                         RooAbsReal& _bkg2,
                         RooAbsReal& _slope1,
                         RooAbsReal& _slope2,
                         RooAbsReal& _slope3) :
    RooAbsPdf(name,title),
    x("x","x",this,_x),
    y("y","y",this,_y),
    bkg1("bkg1","bkg1",this,_bkg1),
    bkg2("bkg2","bkg2",this,_bkg2),
    slope1("slope1","slope1",this,_slope1),
    slope2("slope2","slope2",this,_slope2),
    slope3("slope3","slope3",this,_slope3)
  {
  }


  fitFbkg::fitFbkg(const fitFbkg& other, const char* name) :
    RooAbsPdf(other,name),
    x("x",this,other.x),
    y("y",this,other.y),
    bkg1("bkg1",this,other.bkg1),
    bkg2("bkg2",this,other.bkg2),
    slope1("slope1",this,other.slope1),
    slope2("slope2",this,other.slope2),
    slope3("slope3",this,other.slope3)
  {
  }
  Double_t fitFbkg::evaluate() const
  {
      Double_t ret=0;
      if (y==0){
          ret = slope3*x+1-(1+slope1*x)*bkg1-(1+slope2*x)*bkg2*bkg1;
      }else if (y==1){
          ret = bkg1*(1+slope1*x);
      }else{
          ret = bkg2*bkg1*(1+slope2*x);
      }
      return ret ;
  }

  ClassImp(fitFbkgflat)

   fitFbkgflat::fitFbkgflat(const char *name, const char *title,
                          RooAbsReal& _x,
                          RooAbsCategory& _y,
                          RooAbsReal& _bkg1,
                          RooAbsReal& _bkg2) :
     RooAbsPdf(name,title),
     x("x","x",this,_x),
     y("y","y",this,_y),
     bkg1("bkg1","bkg1",this,_bkg1),
     bkg2("bkg2","bkg2",this,_bkg2)
   {
   }


   fitFbkgflat::fitFbkgflat(const fitFbkgflat& other, const char* name) :
     RooAbsPdf(other,name),
     x("x",this,other.x),
     y("y",this,other.y),
     bkg1("bkg1",this,other.bkg1),
     bkg2("bkg2",this,other.bkg2)
   {
   }
   Double_t fitFbkgflat::evaluate() const
   {
       Double_t ret=0;
       if (y==0){
           ret = 1-bkg1-bkg2*bkg1;
       }else if (y==1){
           ret = bkg1;
       }else{
           ret = bkg2*bkg1;
       }
       return ret ;
   }




