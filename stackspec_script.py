import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column
from glob import glob
import time
import multiprocessing as mp
import argparse
from astropy.cosmology import FlatLambdaCDM
#you can also use pre-defined parameters, e.g.:
#from astropy.cosmology import WMAP7
import astropy.units as u
import pdb
#define the cosmology (if you import WMAP7, you don't need this line)
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3) 

def clip(vec,sigclip,cthresh):
  '''
  
  '''
  
  old=vec.std()
  ii=(np.abs(vec-vec.mean()) <= sigclip*old)
  clipvec=vec[ii]
  #sigma=[]#dblarr(1)
  sigma=clipvec.std()
  veci=clipvec
  while np.abs(sigma-old)/old > cthresh: 
     old = sigma
     ii=(np.abs(vec-veci.mean()) <= sigclip*veci.std())
     veci=vec[ii]
     sigma=veci.std()
  avg=veci.mean()
  dev=veci.std()
  return avg


def removebad(spec):
    '''
    This function can be used to determine the range of good pixel values, 
    such that the zero-value pixels in the beginning and the end of the 
    spectra is removed.
    
    Input
    =====
    spectra: an array of the spec flux values
    
    Output
    ======
    A list of indices that include non-zero spec values (i.e., exludes
    zero specs from the beginning and the end of the spectra)
    '''
    
    if spec[0] != 0: bindx=0
    else:
        for i,s in enumerate(spec):
            if s != 0: break
            bindx=i+10
            
    arr=range(len(spec))
    if spec[-1] != 0: eindx=len(spec)-1
    else:
        for i in reversed(arr):
            if spec[i] != 0: break
            eindx=i-10
    
    
    return np.arange(bindx,eindx+1)

def read_spec(file):
    sp=fits.open(file)
    header = sp[1].header
    wcs = WCS(header)
    
    #read flux and flux err
    spec=sp[1].data
    specerr=sp[2].data
    #convert pixel values to wavelength
    #make index array
    index = np.arange(len(spec)) # or np.arange(header['NAXIS1'])
    specwl = wcs.wcs_pix2world(index[:,np.newaxis], 0)
    specwl=np.reshape(specwl, len(specwl)) #reshape to have a 1d array
    
    return header,spec,specerr,specwl

###################### MAIN CODE ########################
def main(specpath,tblpath,obj_ind,outfile,*normto):
    tbl=fits.getdata(tblpath, ext=1)

#input arrays
    masklist=tbl['maskname'][obj_ind]
    idlist=tbl['id'][obj_ind]
    zarray=tbl['z_mosfire'][obj_ind]
    ha6565_lum=tbl['ha6565_lum'][obj_ind]
    hb4863_lum=tbl['hb4863_lum'][obj_ind]
    luvcorr=tbl['luv_ebmv_smc'][obj_ind]
    ha6565_abs_flux=tbl['ha6565_abs_flux'][obj_ind]
    hb4863_abs_flux=tbl['hb4863_abs_flux'][obj_ind]
    ha6565_lum_err=tbl['ha6565_lum_err'][obj_ind]
    hb4863_lum_err=tbl['hb4863_lum_err'][obj_ind]
    ebmv=tbl['ebmv_smc'][obj_ind]
    fesc=np.exp(-5.214* ebmv**0.312) #f_escape based on Reddy+2016b
    fesc[np.isfinite(fesc) == False] = 0.
    fesc[fesc == 1.] = .4 #to avoid unrealistically high xi_ion
    idliststr=[str(e) for e in idlist]

    specfile = specpath + masklist.strip()+'.*.'+idliststr+'.ell.1d.fits'

#counting the number of files
    speccount=0
    for i,f in enumerate(specfile):
        filelisttemp=glob(f)
        for n in filelisttemp:
            speccount += 1

    #making a list of file paths
    filelist=['' for i in range(speccount)]
    zlist=np.zeros(speccount)
    halumlist=np.zeros(speccount)
    hblumlist=np.zeros(speccount)
    uvlumlist=np.zeros(speccount)
    ha6565_abs_fluxlist=np.zeros(speccount)
    hb4863_abs_fluxlist=np.zeros(speccount)
    fesclist=np.zeros(speccount)
    specind=0
    for i,f in enumerate(specfile):
        filelisttemp=glob(f)
        for n in filelisttemp:
            filelist[specind]=n
            zlist[specind]=zarray[i]
            halumlist[specind] = ha6565_lum[i]
            hblumlist[specind] = hb4863_lum[i]
            uvlumlist[specind] = luvcorr[i]
            ha6565_abs_fluxlist[specind] = ha6565_abs_flux[i]
            hb4863_abs_fluxlist[specind] = hb4863_abs_flux[i]
            fesclist[specind] = fesc[i]
            specind += 1
    import sys
    if sys.argv[5] == 'Ha': norm = 1./halumlist ; normbalm = 1./ha6565_lum
    if sys.argv[5] == 'Hb': norm = 1./hblumlist ; normbalm = 1./hb4863_lum
    if sys.argv[5] == 'UV': norm = 1./uvlumlist ; normbalm = 1./luvcorr
    if sys.argv[5] == 'UV_fesc': norm = (1.+fesclist)/uvlumlist ; normbalm = (1.+fesc)/luvcorr
    if sys.argv[5] == 'none': norm = np.ones(len(halumlist)) ; normbalm = np.ones(len(ha6565_lum))
    if (sys.argv[5] != 'Ha') & (sys.argv[5] != 'Hb') & (sys.argv[5] != 'UV') & (sys.argv[5] != 'none') & (sys.argv[5] != 'UV_fesc'):
        print('normto keyword should be set to one of these: "Ha","Hb","UV", "UV_fesc", or "none" ')

#make a grid of wavelength with the desired resolution
    wavemin = 3250
    wavemax = 10000
    delwave = 0.5                 #in AA
    gridwl = np.arange(wavemin,wavemax,delwave)
    nwave = len(gridwl)
    specall = np.zeros((nwave,speccount))
    specerrall = np.zeros((nwave,speccount))

    t0=time.time()
    print('# Stacking ',speccount,' spectra of ',len(idlist),' objects')
    for i,specfile in enumerate(filelist):
 
        header,spec,specerr,specwl = read_spec(specfile)
    
    #cut the bad beginning and end of the spectra
        goodpix=removebad(spec)
        spec=spec[goodpix]
        specerr=specerr[goodpix]
        specwl=specwl[goodpix]
    
        ldist=cosmo.luminosity_distance(zlist[i]).to(u.cm)
        lspec =1e-40 *ldist*ldist*4*np.pi*(1+zlist[i]) * spec
        lspecerr =1e-40 *ldist*ldist*4*np.pi*(1+zlist[i]) * specerr
    
        lspec = lspec * norm[i]
        lspecerr =  lspecerr * norm[i]
    
    #calculate rest-frame wavelength
        specwl /= 1.+zlist[i]
    
    #interpolate to the new wavelength grid
        UNDEF = -999.
        gridspec = np.interp(gridwl,specwl,lspec,left=UNDEF,right=UNDEF)
        gridspecerr = np.interp(gridwl,specwl,lspecerr,left=UNDEF,right=UNDEF)
    
    #take into account the resampling for the error spectrum
        errfac = np.sqrt(header['CDELT1']/(1.+zlist[i])/delwave)
        gridspecerr = gridspecerr*errfac
                   
    #remove sky lines (remove those with error > 3. * median(err))
        mederr = np.median(gridspecerr[gridspecerr > 0.])
        keep = (np.abs(gridspecerr) < mederr*3.)
        gridspec[keep == False] = UNDEF
        #pdb.set_trace()
    #assign the spetra and error spectra to arrays:
        specall[:,i]=gridspec
        specerrall[:,i]=gridspecerr
        
    #declare arrays
    wt_avg=np.zeros(nwave)
    nwt_avg=np.zeros(nwave)
    med=np.zeros(nwave)
    clip_avg=np.zeros(nwave)
    wt_err=np.zeros(nwave)
    disp=np.zeros(nwave)
    z=np.zeros(nwave)
    num=np.zeros(nwave)
    #stacking
    for j in range(nwave):
        speccol = specall[j,:]
        specerrcol = specerrall[j,:]
        good=np.where(speccol != UNDEF)
        ngood=len(speccol[good])
        if ngood == 0: continue
        if ngood == 1:
            speccol=speccol[good]
            specerrcol = specerrcol[good]
            zcol=zlist[good]
            weight=1/(specerrcol*specerrcol)
            wt_avg[j] = np.nansum(speccol*weight)/np.nansum(weight)
            nwt_avg[j] = np.nanmean(speccol)
            med[j] = speccol
            wt_err[j] = np.sqrt(1./np.nansum(weight))
            disp[j] = 0.
            z[j] = zcol.mean()
            num[j] = ngood
        if ngood > 1:
            speccol=speccol[good]
            specerrcol = specerrcol[good]
            zcol=zlist[good]
        
            weight=1/(specerrcol*specerrcol)
            wt_avg[j] = np.nansum(speccol*weight)/np.nansum(weight)
            nwt_avg[j] = np.nanmean(speccol)
            med[j] = np.nanmedian(speccol)
            clip_avg[j] = clip(speccol,3.,.01)
            wt_err[j] = np.sqrt(1./np.nansum(weight))
            disp[j] = np.nanstd(speccol)
            z[j] = zcol.mean()
            num[j] = ngood


    #Measure weighted average Balmer absorption
    ldistarr=cosmo.luminosity_distance(zarray).to(u.cm)
    goodha=np.where((ha6565_abs_flux != -999.) & (ha6565_lum_err != 0.) & (np.isfinite(ha6565_abs_flux))); n1=len(ha6565_lum_err[goodha])
    goodhb=np.where((hb4863_abs_flux != -999.) & (np.isfinite(hb4863_abs_flux))); n2=len(hb4863_abs_flux[goodhb])
    print('# Number of objs with Balmer absorption of two lines: ',n1,' and ',n2)
    lhaabs = 1e-40 *ldistarr*ldistarr*4*np.pi * ha6565_abs_flux * normbalm
    lhbabs = 1e-40 *ldistarr*ldistarr*4*np.pi * hb4863_abs_flux * normbalm
    habalm = (np.ma.masked_invalid(lhaabs[goodha]/(ha6565_lum_err[goodha]*ha6565_lum_err[goodha])).sum()/
              np.ma.masked_invalid(1/(ha6565_lum_err[goodha]*ha6565_lum_err[goodha])).sum() )
    hbbalm = (np.ma.masked_invalid(lhbabs[goodhb]/(hb4863_lum_err[goodhb]*hb4863_lum_err[goodhb])).sum()/
              np.ma.masked_invalid(1/(hb4863_lum_err[goodhb]*hb4863_lum_err[goodhb])).sum() )

    #defining the columns for the output
    #wt_avg_col=fits.Column(name='wt_avg',format='D',array=wt_avg)
    #nwt_avg_col=fits.Column(name='nwt_avg',format='D',array=nwt_avg)
    #med_col=fits.Column(name='med',format='D',array=med)
    #wt_err_col=fits.Column(name='wt_err',format='D',array=wt_err)
    #disp_col=fits.Column(name='disp',format='D',array=disp)
    #z_col=fits.Column(name='z',format='D',array=z)
    #num_col=fits.Column(name='num',format='D',array=num)

    #cols = fits.ColDefs([wt_avg_col,nwt_avg_col,med_col,wt_err_col,disp_col,z_col,num_col])
    #define the output file
    #out = fits.BinTableHDU.from_columns(cols)
    
    out=fits.PrimaryHDU()
    hdr=out.header

    #making the header of the output file
    out.header['UNITS'] = '1.d40 erg/s/A'
    out.header['CTYPE1'] = 'LINEAR'
    out.header['CRPIX1'] = 1.0
    out.header['CRVAL1'] = wavemin
    out.header['CDELT1'] = delwave
    out.header['CD1_1'] = delwave
    out.header['haabs'] = habalm.value
    out.header['hbabs'] = hbbalm.value
    out.header['uvcor'] = np.median(uvlumlist)
    out.header['uvcorerr'] = uvlumlist.std()/len(uvlumlist)
    out.header['COMMENT'] = 'Ext 1: weighted average'
    out.header['COMMENT'] = 'Ext 2: error'
    out.header['COMMENT'] = 'Ext 3: unweighted average'
    out.header['COMMENT'] = 'Ext 4: median'
    out.header['COMMENT'] = 'Ext 5: 3sigma-clipped average'
    out.header['COMMENT'] = 'Ext 6: standard deviation in each wavelength bin'
    out.header['COMMENT'] = 'Ext 7: average redshift of objs contributing to each wavelength bin'
    out.header['COMMENT'] = 'Ext 8: number of objs in each wavelength bin'
    for i in range(len(idlist)):
        hdrspecname = masklist[i]+'_'+str(idlist[i])
        out.header['COMMENT'] = hdrspecname

    #Writing the output
    #out.writeto('test_py.fits', clobber=True)
    outname=outfile
    out.writeto(outname, clobber=True)
    hdr['OBJ']='wt_avg'
    fits.append(outname,wt_avg,hdr)
    hdr['OBJ']='wt_err'
    fits.append(outname,wt_err,hdr)
    hdr['OBJ']='nwt_avg'
    fits.append(outname,nwt_avg,hdr)
    hdr['OBJ']='med'
    fits.append(outname,med,hdr)
    hdr['OBJ']='clip_avg'
    fits.append(outname,clip_avg,hdr)
    hdr['OBJ']='disp'
    fits.append(outname,disp,hdr)
    hdr['OBJ']='z'
    fits.append(outname,z,hdr)
    hdr['OBJ']='num'
    fits.append(outname,num,hdr)

    print('# Stacking took: {0:3.2f}'.format((time.time()-t0)/60.),' min.')

#import cProfile
#cProfile.run('main()')
import sys
specpath=sys.argv[1]
tblpath=sys.argv[2]
objfile=sys.argv[3]
outfile=sys.argv[4]

#calculate number of lines in the file
with open(objfile) as f:
    n=sum(1 for _ in f)

#put the input in arr array
arr=np.zeros(n)
for i,line in enumerate(open(objfile)):
    arr[i]=line.strip()
#covert 0,1 arr to a boolean array
obj_ind = arr == 1

print('# Spectra is normalized to ',sys.argv[5])
print('# Stacked spectra file : ',sys.argv[4])

main(specpath,tblpath,obj_ind,outfile)
