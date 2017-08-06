# SPECLINE

The main goal is to stack multiple spectra, make a composite spectrum, and fit optical emission lines with Gaussian profiles to calculate their luminosities. 

## Code Details

### stack_driver.py

This is the driver that calls the two main codes [stackspec_script.py](#stackspec_scriptpy) and [linefit.py](#linefitpy).

### stackspec_script.py

This code stacks the spectra. The output has eight extensions:

* Ext 1: weighted average
* Ext 2: error
* Ext 3: unweighted average
* Ext 4: median
* Ext 5: 3 sigma-clipped average
* Ext 6: standard deviation in each wavelength bin
* Ext 7: average redshift of objs contributing to each wavelength bin
* Ext 8: number of objs in each wavelength bin

### linefit.py

This code fits H-alpha, H-beta, and [OIII]5007 lines with a single Gaussian profile, and [OII]3727 line with a double Gaussian profile to calculate line luminosities and their errors