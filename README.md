# Making Composite Spectrum and Fitting Emission Lines

The python scripts presented here are to stack multiple spectra, make a composite spectrum, and fit optical emission lines with Gaussian profiles to calculate the line luminosities. 

The details of the stacking method is discussed at *Shivaei et al. (2017, ApJ, in press)*.

Contributions are very welcomed!


## Code Details

### stack_driver.py

This is the driver that calls the two main codes [stackspec_script.py](#stackspec_scriptpy) and [linefit.py](#linefitpy).

### [stackspec_script.py](#stackspec_scriptpy)

This code stacks the spectra. The output is a fits file with eight extensions:

* Ext 1: weighted average
* Ext 2: error
* Ext 3: unweighted average
* Ext 4: median
* Ext 5: 3 sigma-clipped average
* Ext 6: standard deviation in each wavelength bin
* Ext 7: average redshift of objs contributing to each wavelength bin
* Ext 8: number of objs in each wavelength bin

### [linefit.py](#linefitpy)

This code fits the H-alpha, H-beta, and [OIII]5007 lines with single Gaussian profiles, and the [OII]3727 line with a double Gaussian profile to calculate line luminosities and their errors.

## How to run the code

* The code uses `python 3.6.1`

* Required libraries: `numpy 1.12.1` and `astropy 1.3.2`. Other standard libraries include glob, time, argparse, multiprocessing, subprocess, sys, os.

The driver code is [stack_driver.py](#stack_driverpy). Before running this code you need to modify a few lines:

1. lines 21-24: First you need to specify paths for output and source files
	* specpath: path to the directory that contains the spectra
	* tblpath: path of the catalog of your objects (should be a fits file). For this specific code, the catalog needs to have these entries: 
		* 'maskname' (name of the observed mask), 
		* 'id' (object ID), 'z_mosfire' (redshift), 
		* 'ha6565_lum' (H-alpha luminosity), 
		* 'hb4863_lum' (H-beta luminosity), 
		* 'ha6565_lum_err' (H-alpha luminosity error), 
		* 'hb4863_lum_err' (H-beta luminosity error), 
		* 'luv' (dust-corrected UV luminosity), 
		* 'ha6565_abs_flux' (H-alpha absorption flux), 
		* 'hb4863_abs_flux' (H-beta absorption flux),
		
		In case you want to modify/remove any of these entries, edit lines 90-99 in [stackspec_script.py](#stackspec_scriptpy)
        
   	* outfile: path for the output composite spectrum (should be a fits file: ~/Example.fits)
   	* lineoutfile: path for the output text file that contains the line luminosity measurements (a text file, e.g. ~/Example.txt)

2. line 28:
	* Here you set the criteria for the objects you want to draw from the main catalog (tblpath above) and stack. For example:
	
    ```
    obj_ind=((tbl['z_mosfire'] < 2.) & (tbl['ha6565_lum'] != -999.))
    ```
    
3. lines 40 and 47:
	* Edit the paths to point to where you saved the stackspec_script.py (line 40) and linefit.py (line 47) scripts.


Once you modify the driver code, you are ready to run it through your terminal. Type in the terminal:

```
python stack_driver.py [arg1] [arg2] [arg3-optional]
```

The first argument `[arg1]` indicates the stack spectra normalization. You have four choices: 

* 'Ha' : each individual spectra will be first normalized to H-alpha luminosity, taken from the objects catalog (tblpath, see above).
* 'Hb' : same as above, but normalized to H-beta luminosity.
* 'UV' : same as above, but normalized to UV luminosity.
*  'none' : no normalization will be used.

You can define other normalization terms in line 134 of [stackspec_script.py](#stackspec_scriptpy).

The second argument `[arg2]` indicates which line you want to measure its luminosity at the end. You have these options: 
* 'Ha' : H-alpha 6565 ang 
* 'Hb' : 'H-beta' 4863 ang
* 'oiii' : [OIII] 5008 ang
* 'oii' : [OII] 3727 ang doublet

If you want to fit two lines, you can add the third optional argument `[arg3-optional]`.

*You can only call these combinations (the order is important!):*

* `'Ha' 'Hb'`
* `'Ha' 'oiii'`
* `'Ha' 'oii'`
* `'Hb' 'oiii'`
* `'Hb' 'oii'`
* `'oiii' 'oii'`

**Example:**

```
python stack_driver.py 'Ha' 'Ha' 'Hb' 
```
This will normilize the individual spectra to H-alpha luminosity, stack them, and measure the Ha and Hb luminosities from the composite spectrum.



## Acknowledgment

The details of the stacking method is discussed at *Shivaei et al. (2017, ApJ, in press)*.