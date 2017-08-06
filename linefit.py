import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column
from glob import glob
import time
import multiprocessing as mp
import argparse
from astropy.cosmology import FlatLambdaCDM
#you can also use pre-defined parameters, e.g.: from astropy.cosmology import WMAP7
import astropy.units as u
import pdb
from scipy import optimize

#define the cosmology (if you import WMAP7, you don't need this line)
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3) 

import matplotlib.pyplot as plt

def read_spec(file):
    sp=fits.open(file)
    header = sp[1].header
    wcs = WCS(header)
    
    #read flux and flux err
    spec=sp[3].data
    specerr=sp[6].data/np.sqrt(sp[8].data)
    #convert pixel values to wavelength
    #make index array
    index = np.arange(spec.size) # or np.arange(header['NAXIS1'])
    specwl = wcs.wcs_pix2world(index[:,np.newaxis], 0)
    specwl=np.reshape(specwl, specwl.size) #reshape to have a 1d array
    
    return header,spec,specerr,specwl

def find_ha_arrs(specwl,spec,specerr):
    #find the wavelength, flux, and flux err range for Ha
    sind1 = np.where(specwl == 6555.)[0]
    find1 = np.where(specwl == 6575.)[0]
    start = np.where(specwl == 6400.)[0]
    sind0 = 0.
    if spec[start] != 0: sind0=start
    else:
        for i in range(200):
            if spec[start+i] !=0: break
            sind0=start+i
    find0 = np.where(specwl == 6530.)[0]
    sind2 = np.where(specwl == 6600.)[0]
    final = np.where(specwl == 6680.)[0]
    find2 = 0
    if spec[final] != 0: find2=final
    else:
        for i in reversed(range(200)):
            if spec[final+i] !=0: break
            find2=final+i
    sind0=sind0[0]; sind1=sind1[0]; sind2=sind2[0]
    find0=find0[0]; find1=find1[0]; find2=find2[0]
    w0=specwl[sind0:find0]; w1=specwl[sind1:find1]; w2=specwl[sind2:find2]
    wind_ha=np.append(np.append(w0,w1),w2)
    s0=spec[sind0:find0]; s1=spec[sind1:find1]; s2=spec[sind2:find2]
    spec_ha0=np.append(np.append(s0,s1),s2)
    e0=specerr[sind0:find0]; e1=specerr[sind1:find1]; e2=specerr[sind2:find2]
    e=np.append(np.append(e0,e1),e2)

    return wind_ha,spec_ha0,e

def find_hb_arrs(specwl,spec,specerr):
    #find the wavelength, flux, and flux err range for Hb
    sind1 = np.where(specwl == 4750.)[0]
    find1 = np.where(specwl == 4950.)[0]
    sind1=sind1[0]; find1=find1[0]
    wind_hb=specwl[sind1:find1]
    spec_hb0=spec[sind1:find1]
    e=specerr[sind1:find1]
    return wind_hb,spec_hb0,e

def find_oiii_arrs(specwl,spec,specerr):
    #find the wavelength, flux, and flux err range for OIII5008
    sind1 = np.where(specwl == 4979.)[0]
    find1 = np.where(specwl == 5060.)[0]
    sind1=sind1[0]; find1=find1[0]
    wind_oiii=specwl[sind1:find1]
    spec_oiii0=spec[sind1:find1]
    e=specerr[sind1:find1]
    return wind_oiii,spec_oiii0,e

def find_oii_arrs(specwl,spec,specerr):
    #find the wavelength, flux, and flux err range for OII3727
    sind1 = np.where(specwl == 3675.)[0]
    find1 = np.where(specwl == 3775.)[0]
    sind1=sind1[0]; find1=find1[0]
    wind_oii=specwl[sind1:find1]
    spec_oii0=spec[sind1:find1]
    e=specerr[sind1:find1]
    return wind_oii,spec_oii0,e

def fit_gauss(p0,wave,spec,e):
    '''
    Fit a gaussian function with a continuum to data.
    
    Input
    =====
    p0: initial guess for the guassian parameters in this form: [amplitude,center,width,continuum]
    wave: wavelength array over which the fit should be done
    spec: flux array corresponding to the wavelength array
    e: flux error array
    
    Output
    ======
    p1: the best fit parameters in this form: [amplitude,center,width,continuum]
    '''
    
    #define a gaussian function and an error function for the xi^2 minimization
    gauss = lambda p, t : p[3] + p[0]*np.exp(-(t-p[1])**2/(2*p[2]**2))
    errfunc = lambda p, t, y, err: (y-gauss(p,t))/err
    
    #minimize the sum of squares of errfunc
    p1, success = optimize.leastsq(errfunc, p0[:], args=(wave,spec,e))
    if success > 4.: print('# Not a good fit, success = ',success)
        
    return p1

def fit_double_gauss(p0,wave,spec,e):
    '''
    Fit a double-gaussian function with a continuum to data.
    
    Input
    =====
    p0: initial guess for the guassian parameters in this form:
    [amplitude1,center1,width1,amplitude2,center2,width2,continuum]
    wave: wavelength array over which the fit should be done
    spec: flux array corresponding to the wavelength array
    e: flux error array
    
    Output
    ======
    p1: the best fit parameters in this form: [amplitude,center,width,continuum]
    '''
    
    #define a gaussian function and an error function for the xi^2 minimization
    gauss = lambda p, t : p[6] + p[0]*np.exp(-(t-p[1])**2/(2*p[2]**2)) + p[3]*np.exp(-(t-p[4])**2/(2*p[5]**2))
    errfunc = lambda p, t, y, err: (y-gauss(p,t))/err
    
    #minimize the sum of squares of errfunc
    p1, success = optimize.leastsq(errfunc, p0[:], args=(wave,spec,e))
    if success > 4.: print('# Not a good fit, success = ',success)
        
    return p1

def ha(header, spec, specerr, specwl):
    wind_ha,spec_ha0,e=find_ha_arrs(specwl,spec,specerr)
    if ((np.log10(spec_ha0.max()) > 2.) | (np.log10(spec_ha0.max()) < -1.)): normfac=1./spec_ha0.max()
    else: normfac=1.
        
    p0 = [.4, 6564., 1., .1]  # Initial guess for the parameters [amplitude,center,width,continuum]
    
    pert=1000
    lha_arr=np.zeros(pert)
    gauss = lambda p, t : p[3] + p[0]*np.exp(-(t-p[1])**2/(2*p[2]**2))
    errfunc = lambda p, t, y, err: (y-gauss(p,t))/err
    for p in range(pert):
        spec_ha = e*np.random.randn(spec_ha0.size)+spec_ha0
        p1=fit_gauss(p0,wind_ha,spec_ha*normfac,e*normfac)
        #area under the line:
        area=np.trapz(gauss(p1,np.arange(6555,6575,.1))-p1[3],np.arange(6555,6575,.1))
        lha_arr[p] = area

    #p1, success = optimize.leastsq(errfunc, p0[:], args=(wind_ha,spec_ha0*normfac,e*normfac))
    #plt.clf
    #plt.figure()
    #plt.errorbar(wind_ha,spec_ha0,e,fmt='.',capsize=2,markersize=3,ecolor='grey')
    #plt.plot(wind_ha,gauss(p1,wind_ha), '-',color='orange')
    
    lha=np.mean(lha_arr[lha_arr > 0.])/normfac
    #print('L(Ha) = {0:3.3}'.format(lha))
    #print('Detection at {0:4.4}'.format(lha_arr.mean()/np.std(lha_arr)),' sigma')
       
       
    #Balmer correction
    #print('L(Ha_abs)/L(Ha) = {0:4.4}'.format(header['haabs']/lha))
    lha = lha + header['haabs']
    #print('L(Ha)+L(Ha_abs) = {0:3.4}'.format(lha),' +/- {0:0.3}'.format(lha_arr.std()))
     
    sfrha = 0.079 * lha / 1.8       #in Chabrier

    print('#')
    print('# L(Ha), L(Ha_abs), L(Ha)_err')
    print('{0:e}'.format(lha),' {0:e}'.format(header['haabs']),' {0:e}'.format(lha_arr.std()))

    return lha,header['haabs'],lha_arr.std()

def hb(header, spec, specerr, specwl):
    
    wind_hb,spec_hb0,e=find_hb_arrs(specwl,spec,specerr)
    if ((np.log10(spec_hb0.max()) > 2.) | (np.log10(spec_hb0.max()) < -1.)): normfac=1./spec_hb0.max()
    else: normfac=1.
    
    p0 = [1., 4863., 1., .1]  # Initial guess for the parameters [amplitude,center,width,continuum]
    
    pert=1000
    lhb_arr=np.zeros(pert)
    gauss = lambda p, t : p[3] + p[0]*np.exp(-(t-p[1])**2/(2*p[2]**2))
    errfunc = lambda p, t, y, err: (y-gauss(p,t))/err
    for p in range(pert):
        spec_hb = e*np.random.randn(spec_hb0.size)+spec_hb0
        p1=fit_gauss(p0,wind_hb,spec_hb*normfac,e*normfac)
        #area under the line:
        area=np.trapz(gauss(p1,np.arange(4760,4940,.1))-p1[3],np.arange(4760,4940,.1))
        lhb_arr[p] = area

    lhb=np.mean(lhb_arr[lhb_arr > 0.])/normfac
    #print('L(Hb) = {0:3.3}'.format(lhb))
    #print('Detection at {0:4.4}'.format(lhb_arr.mean()/np.std(lhb_arr)),' sigma')


    #Balmer correction
    #print('L(Hb_abs)/L(Hb) = {0:4.4}'.format(header['hbabs']/lhb))
    lhb = lhb + header['hbabs']
    #print('L(Hb)+L(Hb_abs) = {0:3.4}'.format(lhb),' +/- {0:0.3}'.format(lhb_arr.std()))
    
    print('#')
    print('# L(Hb), L(Hb_abs), L(Hb)_err')
    print('{0:e}'.format(lhb),' {0:e}'.format(header['hbabs']),' {0:e}'.format(lhb_arr.std()))
    return lhb,header['hbabs'],lhb_arr.std()

def oiii(header, spec, specerr, specwl):
    
    wind_oiii,spec_oiii0,e=find_oiii_arrs(specwl,spec,specerr)
    if ((np.log10(spec_oiii0.max()) > 2.) | (np.log10(spec_oiii0.max()) < -1.)): normfac=1./spec_oiii0.max()
    else: normfac=1.
    
    p0 = [1., 5008., 1., .1]  # Initial guess for the parameters [amplitude,center,width,continuum]
    
    pert=1000
    loiii_arr=np.zeros(pert)
    gauss = lambda p, t : p[3] + p[0]*np.exp(-(t-p[1])**2/(2*p[2]**2))
    errfunc = lambda p, t, y, err: (y-gauss(p,t))/err
    for p in range(pert):
        spec_oiii = e*np.random.randn(spec_oiii0.size)+spec_oiii0
        p1=fit_gauss(p0,wind_oiii,spec_oiii*normfac,e*normfac)
        #area under the line:
        area=np.trapz(gauss(p1,np.arange(4979,5050,.1))-p1[3],np.arange(4979,5050,.1))
        loiii_arr[p] = area

    loiii=np.mean(loiii_arr[loiii_arr > 0.])/normfac
  
    print('#')
    print('# L(OIII), L(OIII)_err')
    print('{0:e}'.format(loiii),' {0:e}'.format(loiii_arr.std()))
    return loiii,loiii_arr.std()

def oii(header, spec, specerr, specwl):
    
    wind_oii,spec_oii0,e=find_oii_arrs(specwl,spec,specerr)
    if ((np.log10(spec_oii0.max()) > 2.) | (np.log10(spec_oii0.max()) < -1.)): normfac=1./spec_oii0.max()
    else: normfac=1.
    
    p0 = [1., 3727., 1., 1., 3729., 1., .1]  # Initial guess for the parameters [amplitude1,center1,width1,amplitude2,center2,width2,continuum]
    
    pert=1000
    loii_arr=np.zeros(pert)
    gauss = lambda p, t : p[6] + p[0]*np.exp(-(t-p[1])**2/(2*p[2]**2)) + p[3]*np.exp(-(t-p[4])**2/(2*p[5]**2))
    errfunc = lambda p, t, y, err: (y-gauss(p,t))/err
    for p in range(pert):
        spec_oii = e*np.random.randn(spec_oii0.size)+spec_oii0
        p1=fit_gauss(p0,wind_oii,spec_oii*normfac,e*normfac)
        #area under the line:
        area=np.trapz(gauss(p1,np.arange(3690,3745,.1))-p1[6],np.arange(3690,3745,.1))
        loii_arr[p] = area

    loii=np.mean(loii_arr[loii_arr > 0.])/normfac
  
    print('#')
    print('# L(OII), L(OII)_err')
    print('{0:e}'.format(loii),' {0:e}'.format(loii_arr.std()))
    return loii,loii_arr.std()

def main(file,*arg):
    import sys
    file=sys.argv[1]
    header, spec, specerr, specwl = read_spec(file)

    if sys.argv[2] == 'Ha': ha(header, spec, specerr, specwl)
    if sys.argv[2] == 'Hb': hb(header, spec, specerr, specwl)
    if sys.argv[2] == 'O3': oiii(header, spec, specerr, specwl)
    if sys.argv[2] == 'O2': oii(header, spec, specerr, specwl)
    if len(sys.argv) == 4:
        if sys.argv[3] == 'Ha': ha(header, spec, specerr, specwl)
        if sys.argv[3] == 'Hb': hb(header, spec, specerr, specwl)

import sys
file=sys.argv[1]
main(file)
