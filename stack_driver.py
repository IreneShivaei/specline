def driver(*args):
    '''
    the first argument is for the spectra stack normalization, and it can be one of these three options: 'Ha', 'Hb', 'UV', 'none'

    the second argument shows which line you want to fit a gaussian to: 'Ha', 'Hb'
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
    import pdb
    t00=time.time()

    #paths and source files
    specpath='/Volumes/irene.shivaei/MOSDEF/Data/v0.2+0.4/'
    tblpath='/Users/ireneshivaei/irene.shivaei/MOSDEF/xi_ion/xi.fits'
    outfile='/Users/ireneshivaei/irene.shivaei/MOSDEF/xi_ion/o32/py/hanhb_normtouv_o32_loglt0.fits'
    lineoutfile='/Users/ireneshivaei/irene.shivaei/MOSDEF/xi_ion/o32/py/linefit_hanhb_normtouv_o32_loglt0.txt'
    
    #objects to stack
    tbl=fits.getdata('/Users/ireneshivaei/irene.shivaei/MOSDEF/xi_ion/xi.fits', ext=1)
    obj_ind=((tbl['ha6565_detflag'] == 0) & (tbl['hb4863_detflag'] != -999.) & ((tbl['agnflag'] == 0) | (tbl['agnflag'] == 5) | (tbl['agnflag'] == 6)) &
             (tbl['o32'] > 0.) & (np.log10(tbl['o32']) < 0.))
    print(tbl['o32'][obj_ind].size)
    #pdb.set_trace()
    obj_ind=obj_ind*1. #convert boolean to 0,1 array
    
    def writetxt(filename, x):
        with open(filename,'w') as f:
            for a in x: print(a, file=f)
            
    writetxt("obj_ind.txt",obj_ind)
    objfile="obj_ind.txt"

    l="python3 /Users/ireneshivaei/irene.shivaei/MOSDEF/stackspec/stackspec_script.py "+specpath+' '+tblpath+' '+objfile+' '+outfile
    
    if sys.argv[1] == 'Ha': os.system(l+' "Ha" > '+lineoutfile)
    if sys.argv[1] == 'Hb': os.system(l+' "Hb" > '+lineoutfile)
    if sys.argv[1] == 'UV': os.system(l+' "UV" > '+lineoutfile)
    if sys.argv[1] == 'none': os.system(l+' "none" > '+lineoutfile)
    
    l="python3 /Users/ireneshivaei/irene.shivaei/MOSDEF/linefit/linefit_ha.py "+outfile
    if len(sys.argv) == 3:
        if sys.argv[2] == 'Ha': os.system(l+' "Ha" >> '+lineoutfile)
        if sys.argv[2] == 'Hb': os.system(l+' "Hb" >> '+lineoutfile)
     
    if len(sys.argv) == 4:
        if (sys.argv[2] == 'Ha') & (sys.argv[3] == 'Hb'): os.system(l+' "Ha" "Hb" >> '+lineoutfile)

    l="rm "+objfile
    os.system(l)
    print('The program took: {0:3.2f}'.format((time.time()-t00)/60.),' min.')
    print('~~~~ Bye!! ~~~~')
      
driver()
