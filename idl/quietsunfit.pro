@~/idl/ftsread.pro
@~/idl/ver.pro
@~/idl/line_min.pro
@~/idl/min_line2.pro
@~/idl/whereq

pro quietsunfit,spectrumy,spectrumx=spectrumx,prefilter=prefilter,halfwidth=halfwidth
  if not keyword_set(spectrumx) then begin ; Wavelength calibration
     print,'Please, enter a wavelength range (Angstroms)'
     print,'First wavelength (A):'
     read,x0
     print,'Last wavelength (A):'
     read,x1
     x=[x0,x1]
  endif else begin
     x=spectrumx
  endelse
  
  ftsread,atlasy,fix(min(x))-1,fix(max(x))+1,xlam=atlasx
  atlasy=atlasy/1e4
  if not keyword_set(spectrumx) then begin ; Wavelength calibration
     print,'Click on two spectral lines to use for wavelength calibration'
     plot,findgen(n_elements(atlasy)),atlasy
     cursor,x,y,/up
     done=0
     i=2
     previmin=x
     while done eq 0 and x-i gt 0 and x+i lt n_elements(atlasy) do begin
        imin=(whereq(atlasy,min(atlasy[x-i:x+i])))[0]
        if imin eq previmin then done=1
        previmin=imin
        i=i+1
     endwhile
     ver,imin,lin=2
     xa1=line_min(atlasy[imin-5:imin+5])+imin-5
     xa1=atlasx[fix(xa1)]+(atlasx[fix(xa1)+1]-atlasx[fix(xa1)])*(xa1-fix(xa1))
     cursor,x,y,/up
     done=0
     i=2
     previmin=x
     while done eq 0 and x-i gt 0 and x+i lt n_elements(atlasy) do begin
        imin=(whereq(atlasy,min(atlasy[x-i:x+i])))[0]
        if imin eq previmin then done=1
        previmin=imin
        i=i+1
     endwhile
     ver,imin,lin=2
     xa2=line_min(atlasy[imin-5:imin+5])+imin-5
     xa2=atlasx[fix(xa2)]+(atlasx[fix(xa2)+1]-atlasx[fix(xa2)])*(xa2-fix(xa2))
     plot,atlasx,atlasy
     ver,xa1,lin=2
     ver,xa2,lin=2
     print,'Click on the same two spectral lines here'
     plot,spectrumy
     cursor,x,y,/up
     done=0
     i=2
     previmin=x
     while done eq 0 and x-i gt 0 and x+i lt n_elements(atlasy) do begin
        imin=(whereq(spectrumy,min(spectrumy[x-i:x+i])))[0]
        if imin eq previmin then done=1
        previmin=imin
        i=i+1
     endwhile
     xb1=line_min(spectrumy[imin-5:imin+5])+imin-5
     ver,xb1,lin=2
     cursor,x,y,/up
     done=0
     i=2
     previmin=x
     while done eq 0 and x-i gt 0 and x+i lt n_elements(atlasy) do begin
        imin=(whereq(spectrumy,min(spectrumy[x-i:x+i])))[0]
        if imin eq previmin then done=1
        previmin=imin
        i=i+1
     endwhile
     ver,imin,lin=2
     xb2=line_min(spectrumy[imin-5:imin+5])+imin-5
     ver,xb2,lin=2
     spectrumx=(findgen(n_elements(spectrumy))-xb1)*(xa2-xa1)/(xb2-xb1) + xa1
     plot,spectrumx,spectrumy
  endif


  if min(atlasx) gt min(spectrumx) or max(atlasx) lt max(spectrumx) then begin
     ftsread,atlasy,fix(min(spectrumx))-1,fix(max(spectrumx))+1,xlam=atlasx
     atlasy=atlasy/1e4
  endif
  ind=where(atlasx ge min(spectrumx) and atlasx le max(spectrumx))
  ampatlas=max(atlasy[ind])-min(atlasy[ind])
  ampspec=max(spectrumy)-min(spectrumy)
  spectrumy2=(spectrumy-min(spectrumy))*ampatlas/ampspec+min(atlasy)

  plot,spectrumx,spectrumy2
  oplot,atlasx,atlasy,lin=2,th=2

  ind=where(atlasx ge min(spectrumx) and atlasx le max(spectrumx))
  xx=atlasx[ind]
  yy=interpol(spectrumy2,spectrumx,xx)

  if not keyword_set(prefilter) then begin
     yvec=yy/atlasy[ind]
     yvec=smooth(yvec,n_elements(yvec)/20)
     result=poly_fit(findgen(n_elements(xx)),yvec,3,yfit=prefilter,/double)

     plot,xx,yy
     oplot,xx,prefilter,lin=2,th=3
     print,'Divide by this prefilter shape? (y/n)'
     ans=''
     while ans ne 'y' and ans ne 'n' do begin
        read,ans
     endwhile

     if ans eq 'y' then $
        yy=yy/prefilter

     plot,xx,yy
     oplot,atlasx,atlasy,lin=2

  endif else prefilter=yy*0.+1.
     
  if not keyword_set(halfwidth) then begin
     print,'Enter convolution kernel half width at half maximum (0 for no convolution)'
     print,'Units are Angstroms'
     read,halfwidth
  endif
  sigma=halfwidth/1.17740
  if halfwidth gt 0 then begin
      dx=median(deriv(atlasx))
     nkernel=fix(sigma/dx*5)*2
     nkernel2=fix(sigma/dx*5)
     kernelx=(findgen(nkernel)-nkernel2)*dx
     kernely=exp(-((kernelx)^2)/(2.*sigma)^2)
     kernely=kernely/total(kernely)

     convlvatlas=convol(atlasy,kernely)
     convlvatlas[0:nkernel2]=convlvatlas[nkernel2]
     convlvatlas[n_elements(convlvatlas)-nkernel2:n_elements(convlvaltas)-1]=convlvatlas[n_elements(convlvatlas)-nkernel2]
  endif else convlvatlas=atlasy

  plot,xx,yy
  oplot,atlasx,convlvatlas,lin=2

  tmp=interpol(prefilter,xx,spectrumx)
  prefilter=tmp
  plot,spectrumx,spectrumy/prefilter
  ampspec=max(spectrumy/prefilter)-min(spectrumy/prefilter)

  ampatlas=max(convlvatlas)-min(convlvatlas)
  norm=ampspec/ampatlas
;  oplot,atlasx,(convlvatlas-min(convlvatlas))*ampspec/ampatlas+min(spectrumy),lin=2
  oplot,atlasx,( convlvatlas +min(spectrumy/prefilter)/norm-min(convlvatlas) )*norm,lin=2
  print,'Width (sigma, Angstroms)=',sigma
  print,'Multiplicative constant=',norm
  print,'Additive constant=',min(spectrumy/prefilter)/norm-min(convlvatlas)
  
  save,file='fit.idl',/var,/compress
  
end
