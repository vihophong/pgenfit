import ROOT
import matplotlib.pyplot as plt
import numpy as np
import array


ftable  = np.load("parms/f_table_wilkinson74.npy",allow_pickle='TRUE')
gr = ROOT.TGraph("parms/neueff.txt","%lg %lg")
def f_func(x,p):
    Z = int(p[0])
    E = p[1]-x[0]
    W = E/0.510998918+1 #unit of mec2
    if (W**2-1)<0:
        p0 = 0
    else:
        p0 = np.sqrt(W**2-1)
    if (W+p0)<0:
        fZ0 = 1./60.*(2*W**4-9*W**2-8)*p0 +1/4*W
    else:
        fZ0 = 1./60.*(2*W**4-9*W**2-8)*p0 +1/4*W*np.log(W+p0)
    An = ftable[Z-2]["An"]
    An_i=[]
    if (E<0.079):
        An_i.append(An[0])
        An_i.append(An[1])
        An_i.append(An[2])
        An_i.append(An[3])
    elif (E>=0.079 and E<0.501):
        An_i.append(An[4])
        An_i.append(An[5])
        An_i.append(An[6])
        An_i.append(An[7])
    elif (E>=0.501 and E<3.162):
        An_i.append(An[8])
        An_i.append(An[9])
        An_i.append(An[10])
        An_i.append(An[11])
    elif (E>=3.162 and E<12.589):
        An_i.append(An[12])
        An_i.append(An[13])
        An_i.append(An[14])
        An_i.append(An[15])
    elif (E>=12.589):
        An_i.append(An[16])
        An_i.append(An[17])
        An_i.append(An[18])
        An_i.append(An[19])
    if (E<0):
        E = 0.001
    S = np.exp(An_i[0]+An_i[1]*np.log(E)+An_i[2]*(np.log(E)**2)+An_i[3]*(np.log(E)**3))
    return fZ0*S          

def f_func_num(x,p):
    return f_func(x,p)*(x[0]-p[2])

def evalEmean(Qb,Sn,ZZ):  
    f_d = ROOT.TF1("f_func",f_func,0.,Qb,2);
    f_d.SetParameter(0,ZZ)
    f_d.SetParameter(1,Qb)
    f_n = ROOT.TF1("f_func_num",f_func_num,0.,Qb,3);
    f_n.SetParameter(0,ZZ)
    f_n.SetParameter(1,Qb)
    f_n.SetParameter(2,Sn)
    return f_n.Integral(Sn,Qb,0.01)/f_d.Integral(Sn,Qb,0.01)

def evalEmeanVar(Qb,Sn,dQb,dSn,ZZ,n=1000):
    retVal = []
    r = ROOT.TRandom3()
    for i in range(n):
        Sn_i = r.Gaus(Sn,dSn)
        Qb_i = r.Gaus(Qb,dQb)
        retVal.append(evalEmean(Qb_i,Sn_i,ZZ))
    retVal = array.array('f',retVal)
    meanfit = ROOT.TMath.Mean(n,retVal)
    mean = evalEmean(Qb,Sn,ZZ)
    sigma = ROOT.TMath.RMS(n,retVal)
    return mean,sigma,meanfit

def cal_eff(E):
    return gr.Eval(E)

def cal_all(Qb,Sn,dQb,dSn,e_factor,Z):
    mean,sigma,meanfit = evalEmeanVar(Qb,Sn,dQb,dSn,Z)
    e_mean = mean * e_factor
    e_mean_fit = meanfit * e_factor
    e_sigma = sigma * e_factor
    e_low  = e_mean - np.sqrt((e_mean-e_mean/2)**2+e_sigma**2)
    if (e_low<0):
        e_low  = e_mean/2
    e_hi = e_mean + np.sqrt((e_mean*2-e_mean)**2+e_sigma**2)
    return e_mean,cal_eff(e_mean),e_low,cal_eff(e_low),e_hi,cal_eff(e_hi),e_mean_fit,e_sigma


#note all values are in keV
def cal1n2neff(nuclide,Z, Qb,D_Qb,S1n,D_S1n,S2n,D_S2n,Qb1n,Qb2n):
    ee_mean=0.;eff_mean=0.668;e_low=0.;eff_low=0.668-0.02;e_hi=0.;eff_hi=0.668+0.02;e_mean_fit=0.;e_sigma=0.
    if (Qb>0):
        if (Qb1n>0):
            e_mean,eff_mean,e_low,eff_low,e_hi,eff_hi,e_mean_fit,e_sigma = cal_all(Qb/1000.,S1n/1000.,D_Qb/1000.,D_S1n/1000.,1.,Z)
    ee_mean2=0.;eff_mean2=0.668;e_low2=0.;eff_low2=0.668-0.02;e_hi2=0.;eff_hi2=0.668+0.02;e_mean_fit2=0.;e_sigma2=0.
    if (Qb>0):
        if (Qb2n>0):
            e_mean2,eff_mean2,e_low2,eff_low2,e_hi2,eff_hi2,e_mean_fit2,e_sigma2 = cal_all(Qb/1000.,S2n/1000.,D_Qb/1000.,D_S2n/1000.,0.5,Z)

    effmean_norm = 66.8 / 100.
    deffmean_norm  = 2. / 100.
    if (e_mean<=0):
        eff_mean = effmean_norm
        eff_low = effmean_norm + deffmean_norm
        eff_hi = effmean_norm - deffmean_norm
    eff_m = eff_mean - eff_hi
    eff_p = eff_low - eff_mean
    if (e_mean2<=0):
        eff_mean2 = effmean_norm
        eff_low2 = effmean_norm + deffmean_norm
        eff_hi2 = effmean_norm - deffmean_norm
    eff_m2 = eff_mean2 - eff_hi2
    eff_p2 = eff_low2 - eff_mean2

    if (eff_m<0):
        print(nuclide,"eff_m <0, ",e_mean)
        eff_mean = effmean_norm
        eff_m = deffmean_norm
        eff_p = deffmean_norm
        
    if (eff_p<0):
        print(nuclide,"eff_p <0, ",e_mean)   
        eff_mean = effmean_norm
        eff_m = deffmean_norm
        eff_p = deffmean_norm      

    if (eff_m2<0):
        print(nuclide,"eff_m2 <0, ",e_mean2)
        eff_mean2 = effmean_norm
        eff_m2 = deffmean_norm
        eff_p2 = deffmean_norm
        
    if (eff_p2<0):
        print(nuclide,"eff_p2 <0, ",e_mean2)
        eff_mean2 = effmean_norm
        eff_m2 = deffmean_norm
        eff_p2 = deffmean_norm
    return eff_mean,eff_m,eff_p,eff_mean2,eff_p2,eff_m2
