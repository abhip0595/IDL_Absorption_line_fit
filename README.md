# IDL_Absorption_line_Voigt_fit


FIT_ABS_LINES:
IDL function to fit absorption (/emission) line with Voigt Profile.

FUNCTIONS CALLED:
	1.	MPFITFUN
	2.	CHEBEVAL

CALLING SEQUENCE
result = fit_abs_lines(w,f,e,line=line[,yfit=yfit,parname=parname,dof=dof,chisq=chisq,quiet=quiet,c_order=c_order])

w - wavevector - i/p
f - flux - i/p
e - 1 sigma error in flux - i/p
line - STRUCTURE: - i/p
	.x0 - [Array] - line centers
	.gamma - [Array] - damping constant
	.nx - [Array] - column density
	.fosc - [Array] - oscillator strength
	.b - [Array] - doppler broadening parameter
	.typ - [Array] - -1 for absorption / +1 for emission

OPTIONAL:
	yfit - [Array] of best fit to Voigt profile - o/p
	parname	- [Array] Voigt Profile parameter names passed to MPFIT	- o/p
	dof - Degrees of Freedom returned by MPFIT - o/p
	chisq - CHI-SQUARE calculated - o/p
	quiet - 1 to suppress MPFIT print statements DEFAULT = 1 - i/p
	c_order - Order of CHEBYSHEV polynomial for continuum DEFAULT=3	- i/p
