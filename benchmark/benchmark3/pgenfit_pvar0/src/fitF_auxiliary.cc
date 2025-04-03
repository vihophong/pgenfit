#include <fitF.hh>
#include <RooAbsReal.h>
#include <RooAbsCategory.h>
#include <math.h>
#include <TMath.h>
#include <TStopwatch.h>
#ifdef EVAL_FAST
Double_t fitF::evaluate() const
{
double t = x;
double l0=(*p[0]);
double l1=(*p[1]);
double l2=(*p[2]);
double l3=(*p[3]);
double l4=(*p[4]);
double l5=(*p[5]);
double l6=(*p[6]);
double l7=(*p[7]);
double l8=(*p[8]);
double e0=exp(-l0*t);
double e1=exp(-l1*t);
double e2=exp(-l2*t);
double e3=exp(-l3*t);
double e4=exp(-l4*t);
double e5=exp(-l5*t);
double e6=exp(-l6*t);
double e7=exp(-l7*t);
double e8=exp(-l8*t);
double p1n0=(*p[9]);
double p1n1=(*p[10]);
double p1n2=(*p[11]);
double p1n3=(*p[12]);
double p1n4=(*p[13]);
double p1n5=(*p[14]);
double p1n6=(*p[15]);
double p1n7=(*p[16]);
double p1n8=(*p[17]);
double p2n0=(*p[18]);
double p2n1=(*p[19]);
double p2n2=(*p[20]);
double p2n3=(*p[21]);
double p2n4=(*p[22]);
double p2n5=(*p[23]);
double p2n6=(*p[24]);
double p2n7=(*p[25]);
double p2n8=(*p[26]);
double py0=(*p[27]);
double py1=(*p[28]);
double py2=(*p[29]);
double py3=(*p[30]);
double py4=(*p[31]);
double py5=(*p[32]);
double py6=(*p[33]);
double py7=(*p[34]);
double py8=(*p[35]);
double ne0=(*p[36]);
double ne1=(*p[37]);
double ne2=(*p[38]);
double ne3=(*p[39]);
double ne4=(*p[40]);
double ne5=(*p[41]);
double ne6=(*p[42]);
double ne7=(*p[43]);
double ne8=(*p[44]);
double be=*p[49];
double b1ne=*p[50];
double b2ne=*p[51];
double n1n2ne=*p[52];
double N0=*p[45]/l0;
double fparentdecay=l0*N0*e0;
double fdecay=fparentdecay*be;
double f0=l1*py1*(1-p1n0-p2n0)*l0*(e0/((l1-l0)*1)+e1/((l0-l1)*1)+0)*N0;
double f1=l2*py2*p1n0*l0*(e0/((l2-l0)*1)+e2/((l0-l2)*1)+0)*N0;
double f2=l3*py3*p2n0*l0*(e0/((l3-l0)*1)+e3/((l0-l3)*1)+0)*N0;
double f3=l4*py1*(1-p1n0-p2n0)*l0*py4*(1-p1n1-p2n1)*l1*(e0/((l1-l0)*(l4-l0)*1)+e1/((l0-l1)*(l4-l1)*1)+e4/((l0-l4)*(l1-l4)*1)+0)*N0;
double f4=l5*py1*(1-p1n0-p2n0)*l0*py5*p1n1*l1*(e0/((l1-l0)*(l5-l0)*1)+e1/((l0-l1)*(l5-l1)*1)+e5/((l0-l5)*(l1-l5)*1)+0)*N0;
double f5=l5*py2*p1n0*l0*py5*(1-p1n2-p2n2)*l2*(e0/((l2-l0)*(l5-l0)*1)+e2/((l0-l2)*(l5-l2)*1)+e5/((l0-l5)*(l2-l5)*1)+0)*N0;
double f7=l6*py2*p1n0*l0*py6*p1n2*l2*(e0/((l2-l0)*(l6-l0)*1)+e2/((l0-l2)*(l6-l2)*1)+e6/((l0-l6)*(l2-l6)*1)+0)*N0;
double f8=l6*py3*p2n0*l0*py6*(1-p1n3-p2n3)*l3*(e0/((l3-l0)*(l6-l0)*1)+e3/((l0-l3)*(l6-l3)*1)+e6/((l0-l6)*(l3-l6)*1)+0)*N0;
double f9=l7*py1*(1-p1n0-p2n0)*l0*py4*(1-p1n1-p2n1)*l1*py7*(1-p1n4-p2n4)*l4*(e0/((l1-l0)*(l4-l0)*(l7-l0)*1)+e1/((l0-l1)*(l4-l1)*(l7-l1)*1)+e4/((l0-l4)*(l1-l4)*(l7-l4)*1)+e7/((l0-l7)*(l1-l7)*(l4-l7)*1)+0)*N0;
double f11=l8*py1*(1-p1n0-p2n0)*l0*py5*p1n1*l1*py8*(1-p1n5-p2n5)*l5*(e0/((l1-l0)*(l5-l0)*(l8-l0)*1)+e1/((l0-l1)*(l5-l1)*(l8-l1)*1)+e5/((l0-l5)*(l1-l5)*(l8-l5)*1)+e8/((l0-l8)*(l1-l8)*(l5-l8)*1)+0)*N0;
double f12=l8*py2*p1n0*l0*py5*(1-p1n2-p2n2)*l2*py8*(1-p1n5-p2n5)*l5*(e0/((l2-l0)*(l5-l0)*(l8-l0)*1)+e2/((l0-l2)*(l5-l2)*(l8-l2)*1)+e5/((l0-l5)*(l2-l5)*(l8-l5)*1)+e8/((l0-l8)*(l2-l8)*(l5-l8)*1)+0)*N0;
fdecay+=f0+f1+f2+f3+f4+f5+f7+f8+f9+f11+f12+0;

double randcoinf2n=*p[48];
double randcoinfgt0n=*p[47];
double randcoinf1n=*p[46];
double fdecay1n=fparentdecay*(be*randcoinf1n+b1ne*ne0*p1n0*(1-randcoinf1n-randcoinfgt0n)+b2ne*2*(n1n2ne*(1-n1n2ne))*p2n0*(1-randcoinf1n-randcoinfgt0n)-b2ne*n1n2ne*n1n2ne*p2n0*randcoinf1n);
fdecay1n+=f0*(randcoinf1n+ne1*p1n1*(1-randcoinf1n-randcoinfgt0n)+2*(ne1*(1-ne1))*p2n1*(1-randcoinf1n-randcoinfgt0n)-ne1*ne1*p2n1*randcoinf1n);
fdecay1n+=f1*(randcoinf1n+ne2*p1n2*(1-randcoinf1n-randcoinfgt0n)+2*(ne2*(1-ne2))*p2n2*(1-randcoinf1n-randcoinfgt0n)-ne2*ne2*p2n2*randcoinf1n);
fdecay1n+=f2*(randcoinf1n+ne3*p1n3*(1-randcoinf1n-randcoinfgt0n)+2*(ne3*(1-ne3))*p2n3*(1-randcoinf1n-randcoinfgt0n)-ne3*ne3*p2n3*randcoinf1n);
fdecay1n+=f3*(randcoinf1n+ne4*p1n4*(1-randcoinf1n-randcoinfgt0n)+2*(ne4*(1-ne4))*p2n4*(1-randcoinf1n-randcoinfgt0n)-ne4*ne4*p2n4*randcoinf1n);
fdecay1n+=f4*(randcoinf1n+ne5*p1n5*(1-randcoinf1n-randcoinfgt0n)+2*(ne5*(1-ne5))*p2n5*(1-randcoinf1n-randcoinfgt0n)-ne5*ne5*p2n5*randcoinf1n);
fdecay1n+=f5*(randcoinf1n+ne5*p1n5*(1-randcoinf1n-randcoinfgt0n)+2*(ne5*(1-ne5))*p2n5*(1-randcoinf1n-randcoinfgt0n)-ne5*ne5*p2n5*randcoinf1n);
fdecay1n+=f7*(randcoinf1n+ne6*p1n6*(1-randcoinf1n-randcoinfgt0n)+2*(ne6*(1-ne6))*p2n6*(1-randcoinf1n-randcoinfgt0n)-ne6*ne6*p2n6*randcoinf1n);
fdecay1n+=f8*(randcoinf1n+ne6*p1n6*(1-randcoinf1n-randcoinfgt0n)+2*(ne6*(1-ne6))*p2n6*(1-randcoinf1n-randcoinfgt0n)-ne6*ne6*p2n6*randcoinf1n);
fdecay1n+=f9*(randcoinf1n+ne7*p1n7*(1-randcoinf1n-randcoinfgt0n)+2*(ne7*(1-ne7))*p2n7*(1-randcoinf1n-randcoinfgt0n)-ne7*ne7*p2n7*randcoinf1n);
fdecay1n+=f11*(randcoinf1n+ne8*p1n8*(1-randcoinf1n-randcoinfgt0n)+2*(ne8*(1-ne8))*p2n8*(1-randcoinf1n-randcoinfgt0n)-ne8*ne8*p2n8*randcoinf1n);
fdecay1n+=f12*(randcoinf1n+ne8*p1n8*(1-randcoinf1n-randcoinfgt0n)+2*(ne8*(1-ne8))*p2n8*(1-randcoinf1n-randcoinfgt0n)-ne8*ne8*p2n8*randcoinf1n);
double fdecay2n=fparentdecay*(b2ne*n1n2ne*n1n2ne*p2n0*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n*be+b1ne*ne0*p1n0*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+b2ne*2*(n1n2ne*(1-n1n2ne))*p2n0*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
fdecay2n+=f0*(ne1*ne1*p2n1*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne1*p1n1*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne1*(1-ne1))*p2n1*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
fdecay2n+=f1*(ne2*ne2*p2n2*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne2*p1n2*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne2*(1-ne2))*p2n2*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
fdecay2n+=f2*(ne3*ne3*p2n3*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne3*p1n3*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne3*(1-ne3))*p2n3*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
fdecay2n+=f3*(ne4*ne4*p2n4*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne4*p1n4*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne4*(1-ne4))*p2n4*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
fdecay2n+=f4*(ne5*ne5*p2n5*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne5*p1n5*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne5*(1-ne5))*p2n5*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
fdecay2n+=f5*(ne5*ne5*p2n5*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne5*p1n5*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne5*(1-ne5))*p2n5*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
fdecay2n+=f7*(ne6*ne6*p2n6*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne6*p1n6*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne6*(1-ne6))*p2n6*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
fdecay2n+=f8*(ne6*ne6*p2n6*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne6*p1n6*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne6*(1-ne6))*p2n6*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
fdecay2n+=f9*(ne7*ne7*p2n7*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne7*p1n7*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne7*(1-ne7))*p2n7*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
fdecay2n+=f11*(ne8*ne8*p2n8*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne8*p1n8*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne8*(1-ne8))*p2n8*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));
fdecay2n+=f12*(ne8*ne8*p2n8*(1-randcoinf2n-randcoinfgt0n)+randcoinf2n+ne8*p1n8*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n)+2*(ne8*(1-ne8))*p2n8*(randcoinf1n*(1-randcoinfgt0n)-randcoinf2n));

double ret; if (y==0) ret= fdecay-fdecay1n-fdecay2n; else if (y==1) ret= fdecay1n; else ret= fdecay2n;
return ret;
}
#endif
