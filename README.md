# IDL_Absorption_line_fit


FIT_ABS_LINES:
IDL function to fit absorption (/emission) line with Voigt Profile.

FUNCTIONS CALLED:
	1.	MPFITFUN
	2.	CHEBEVAL

CALLING SEQUENCE:
result = fit_abs_lines(w,f,e,line=line[,yfit=yfit,parname=parname,dof=dof,chisq=chisq,quiet=quiet,c_order=c_order])


INPUTS:
	w 	- wavevector
	f 	- flux
	e 	- 1 sigma error in flux
	line	- STRUCTURE
		.x0 	- [Array] - line centers
		.gamma 	- [Array] - damping constant
		.nx 	- [Array] - column density
		.fosc 	- [Array] - oscillator strength
		.b 	- [Array] - doppler broadening parameter
		.typ 	- [Array] - -1 for absorption / +1 for emission
	OPTIONAL:
	quiet 	- 1 to suppress MPFIT print statements DEFAULT = 1
	c_order - Order of CHEBYSHEV polynomial for continuum DEFAULT=3

OUTPUT:
	result	- BEST FIT values
	OPTIONAL:
	yfit 	- [Array] of best fit to Voigt profile
	parname	- [Array] Voigt Profile parameter names passed to MPFIT
	dof 	- Degrees of Freedom returned by MPFIT
	chisq 	- CHI-SQUARE calculated

FIT_ABS_BATCH:
IDL code to run fit_abs_lines for batch of files.

CALLING SEQUENCE:
.run fit_abs_batch

DESCRIPTION: User can place their SED data (.fits format) in a 'data/SED' directory. Output is written to a CSV file saved 'abs_fit.csv' in 'data/' directory.
Plots are saved in 'plots/abs_line/' in PNG format. FITS file structure to be .wave for wavevector, .flux for flux and .sigma for 1 sigma error in flux. Users may modify these for their needs.
