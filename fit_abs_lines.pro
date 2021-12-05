function voigt_fit,x,p,tau=tau
	ec = 1.602176634e-19	; electron charge in C
	me = 9.10938356e-31	; electron mass in kg
	cv = 2.99792458e18	; light speed in vacuum (As-1)

	c_order = p[0]	; Order of chebyshev polynomial
	nl = p[c_order+1]
	const = p[1:c_order]
	tau = fltarr(n_elements(x))

	for i = 0,nl-1 do begin
		t = c_order+ 2 + i*6
		x0 = p[t+0]
		fosc = p[t+1]
		gamma = p[t+2]
		nx = 10^p[t+3]
		b = p[t+4]/1e13
		ty = p[t+5]
		ki = ec^2*sqrt(!pi)*fosc*x0/me/cv/b
		ai = x0*gamma/4./!pi/b
		ld = b*x0/cv
		xi = (x-x0)/ld
		ti = ki*nx*ai*voigt(ai,xi)
		tau = tau + ty*ti
	endfor

	prof = exp(tau) * chebeval(x,const)
	return,prof
end

function fit_abs_lines,w,f,e,line=line,yfit=yfit,parname=parname,$
    dof=dof,chisq=chisq,quiet=quiet,c_order=c_order

	nl = n_elements(line.x0)    ; Number of line centers
	npar = (c_order+1) + (nl*6) + 1
	par = replicate({value:0.d,fixed:0,limited:[0,0],$
   	limits:[0.d,0.d],parname:''},npar)

	par[0].value = c_order
	par[0].parname = 'ORDER'
	par[0].fixed = 1
	for i=1,c_order do begin
		par[i].value = 1.
		par[i].parname = 'C'+strtrim(string(i-1),1)
	endfor
	par[c_order+1].value = nl
	par[c_order+1].parname = 'N_LINES'
	par[c_order+1].fixed = 1

	for i = 0,nl-1 do begin
		x = c_order+ 2 + 6*i
		par[x + 0].value = line[i].x0
		par[x + 0].parname = 'LINE_CENTER'
		par[x + 0].limited[0:1] = [1,1]
		par[x + 0].limits[0:1] = [line[i].x0-10.,line[i].x0+10.]
		par[x + 1].value = line[i].fosc
		par[x + 1].parname = 'FOSC'
		par[x + 1].fixed = 1
		par[x + 2].value = line[i].gamma
		par[x + 2].parname = 'GAMMA'
		par[x + 2].fixed = 1
		par[x + 3].value = line[i].nx
		par[x + 3].parname = 'NH'
		par[x + 3].limited[0:1] = 1
		par[x + 3].limits[0:1] = [0.,50.]
		par[x + 4].value = line[i].b
		par[x + 4].parname = 'B_KMPS'
		par[x + 5].value = line[i].typ
		par[x + 5].fixed = 1
		par[x + 5].parname = 'LINE_TYPE'
	endfor
	parname = par.parname

	result = mpfitfun('voigt_fit', w, f, e,  perror = dparms,$
    yfit = yfit,PARINFO = par,/nan,quiet=quiet,dof=dof)

	chisq = total((f-yfit)^2/e^2,/nan)
	return,result
end
