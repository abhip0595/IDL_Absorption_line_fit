
  openw,lun,'data/abs_fit.csv',/get_lun   ;output file
  ;Put your data here: 'data/SED/'
  files = file_search('data/SED/*.fits',count=nf)

  if nf gt 0 then begin
    stars = strsplit(files,'.',/extract)
    stars = stars.ToArray()
    stars = stars[*,0]
    stars = strsplit(stars,path_sep(),/extract)
    stars = stars.ToArray()
    stars = stars[*,3]

  for si = 0,nf-1 do begin
    star_name	= stars[si]
    c_order = 3   ; Order of polynomial
    xlim = [1150,1300]    ; Wavelength range
    lines = [1215.7]    ; line centers
    gamma = [4.6986e8]    ; gamma
    fosc = [4.1641e-1]    ; fosc
    nx = [20.]    ; column density
    b = [10.]   ; broading parameter
    typ	=	[-1]    ; +1: emission -1: absorption
    rb = 1    ; 1: to rebin data to increase SNR
    quiet = 1   ; Set 1 to supress MPFIT
    pr = 1    ; Set 1 to print result
    sp = 1    ; Set 1 to save plot
    line = replicate({x0:0.,gamma:0.,nx:0.,fosc:0.,$
          b:0.,typ:+1},n_elements(lines))
    for i = 0,n_elements(lines)-1 do begin
      line[i].x0 = lines[i]
      line[i].gamma = gamma[i]
      line[i].nx = nx[i]
      line[i].fosc = fosc[i]
      line[i].b = b[i]
      line[i].typ = typ[i]
    endfor
    im	= mrdfits('data/SED/'+star_name+'.fits',1,hdr)
    w = im.wave
    f	= im.flux
    e	= im.sigma

    if rb eq 1 then begin
      q = where(w ge xlim[0]-100. and w le xlim[1]+100. and (f*e) gt 0.)
      w = w[q]
      f = f[q]
      e = e[q]
      wx =	40.
      nx = (max(w)-min(w))/wx
      for i = 0,nx-1 do begin
        xi = min(w)+i*wx-wx/2.
        q	= where(w ge xi and w le xi+wx,nq)
        if nq gt 0 then begin
          mf	= max(f[q],/nan)
          f[q]=f[q]/mf
          e[q]=e[q]/mf
        endif
      endfor
      q = where(w ge xlim[0] and w le xlim[1] and (f*e) gt 0.)
      w = w[q]
      f = f[q]
      e = e[q]
      endif else begin
        q = where(w ge xlim[0] and w le xlim[1] and (f*e) gt 0.)
        w = w[q]
        f = f[q]
        e = e[q]
        mf = max(f,/nan)
        f = f/mf
        e = e/mf
  endelse

  yrange = [min(f-e,/nan),max(f+e,/nan)]
  angstrom = '!6!sA!r!u!9 %!6!n'
  xtitle = 'Wavelength' + ' ( '+angstrom+' )'
  cgplot,w,f,xtitle=xtitle,ytitle='Normalised Flux',$
          Title=star_name,xrange=xlim,col=cgcolor('black'),yrange=yrange
  oploterror,w,f,e,col=cgcolor('black')

  abs = fit_abs_lines(w,f,e,line=line,yfit=yfit,parname=parname,$
  											dof=dof,chisq=chisq,quiet=quiet,c_order=c_order)

  oplot,w,yfit,col=255,thick=2
  write_png,'plots/abs_line/'+star_name+'.png',tvrd(/true)

  if pr eq 1 then begin
    for i=0,n_elements(parname)-1 do begin
      str = parname[i]+' = '+strtrim(string(abs[i]),1)
      print,str
    endfor
    print,'CHISQR = '+strtrim(string(chisq/dof),1)
  endif

  if si eq 0 then begin
    str = 'PROGRAM STAR'
    for i=0,n_elements(parname)-1 do str = str + ','+parname[i]
    str = str+',REDUCED_CHISQ,DOF'
  printf,lun,str
  endif

  str = star_name
  for i=0,n_elements(parname)-1 do str = str + ','+strtrim(string(abs[i]),1)

  str = str+','+strtrim(string(chisq/dof),1)+','+strtrim(string(dof),1)
  printf,lun,str

  endfor
  free_lun,lun
  endif

end
