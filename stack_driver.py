from __future__ import print_function

def driver(*args):
    '''
    the first argument is for the spectra stack normalization, and it can be one of these three options: 'Ha', 'Hb', 'UV', 'none'

    the second argument shows which line you want to fit a gaussian to: 'Ha', 'Hb', 'oiii', 'oii'
    if you want a second line, you add another argument.
    Example:
    > python3 stack_driver.py 'Ha' 'Ha' 'Hb' #normilize the spectra to Ha luminosity, stack it, and then measure Ha and Hb lines.
    
    modify the paths and the criteria for objects to use in the stacks below.
    '''
    import subprocess
    import sys
    import os
    import numpy as np
    from astropy.io import fits
    import time
    t00=time.time()

    #paths and source files
    specpath='/Path/Data/'
    tblpath='/Path/mosdef0d.fits'
    outfile='/Path/hanhb_normtouv.fits'
    lineoutfile='/Path/linefit_hanhb_normtouv.txt'
    
    #objects to stack
    tbl=fits.getdata(tblpath, ext=1)
    obj_ind=((tbl['z_mosfire'] < 2.) & (tbl['ha6565_lum'] != -999.))
    print('Stacking ',tbl['o32'][obj_ind].size,' objs')

    obj_ind=obj_ind*1. #convert boolean to 0,1 array
    
    def writetxt(filename, x):
        with open(filename,'w') as f:
            for a in x: print(a, file=f)
            
    writetxt("obj_ind.txt",obj_ind)
    objfile="obj_ind.txt"

    l="python3 /PathToCode/stackspec_script.py "+specpath+' '+tblpath+' '+objfile+' '+outfile
    
    if sys.argv[1] == 'Ha': os.system(l+' "Ha" > '+lineoutfile)
    if sys.argv[1] == 'Hb': os.system(l+' "Hb" > '+lineoutfile)
    if sys.argv[1] == 'UV': os.system(l+' "UV" > '+lineoutfile)
    if sys.argv[1] == 'none': os.system(l+' "none" > '+lineoutfile)
    
    l="python3 /PathToCode/linefit.py "+outfile
    if len(sys.argv) == 3:
        if sys.argv[2] == 'Ha': os.system(l+' "Ha" >> '+lineoutfile)
        if sys.argv[2] == 'Hb': os.system(l+' "Hb" >> '+lineoutfile)
        if sys.argv[2] == 'oiii': os.system(l+' "oiii" >> '+lineoutfile)
        if sys.argv[2] == 'oii': os.system(l+' "oii" >> '+lineoutfile)
     
    if len(sys.argv) == 4:
        if (sys.argv[2] == 'Ha') & (sys.argv[3] == 'Hb'): os.system(l+' "Ha" "Hb" >> '+lineoutfile)
        if (sys.argv[2] == 'Ha') & (sys.argv[3] == 'oiii'): os.system(l+' "Ha" "oiii" >> '+lineoutfile)
        if (sys.argv[2] == 'Ha') & (sys.argv[3] == 'oii'): os.system(l+' "Ha" "oii" >> '+lineoutfile)
        if (sys.argv[2] == 'Hb') & (sys.argv[3] == 'oiii'): os.system(l+' "Hb" "oiii" >> '+lineoutfile)
        if (sys.argv[2] == 'Hb') & (sys.argv[3] == 'oii'): os.system(l+' "Hb" "oii" >> '+lineoutfile)
        if (sys.argv[2] == 'oiii') & (sys.argv[3] == 'oii'): os.system(l+' "oiii" "oii" >> '+lineoutfile)

    l="rm "+objfile
    os.system(l)
    print('The program took: {0:3.2f}'.format((time.time()-t00)/60.),' min.')
    print('~~~~ Bye!! ~~~~')
      
driver()
