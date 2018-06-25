; Take a model in NICOLE binary format and interpolate it to a new
; grid. NOTE: The model is smoothed with a 2-pixel window to 
; remove discontinuities if the keyword smooth is present

pro interpolate_model,filein,fileout,newscale,smooth=smoot

if filein eq fileout then begin
   print,'Error. Input and output files must be different'
   stop
endif

nzout=n_elements(newscale)
openr,inunit,filein,/get_lun,/swap_if_big_endian     
openw,outunit,fileout,/get_lun,/swap_if_big

tmp=lon64arr(13)
readu,inunit,tmp                ; Ignore first record (signature)
if tmp(0) ne 3328834590979877230 or tmp(1) ne 2314885530823516726 then begin
   print,'Error. Input file is not in v2.6 format'
   stop
endif   
nzin=tmp(3)
tmp=lonarr(26)
free_lun,inunit
openr,inunit,filein,/get_lun,/swap_if_big_endian   
readu,inunit,tmp                ; Ignore first record (signature)
nx=tmp(4)
ny=tmp(5)
free_lun,inunit
openr,inunit,filein,/get_lun,/swap_if_big_endian     

print,'nx,ny=',nx,ny

tmp=lon64arr(22*nzout+11+92)
tmp(0)=3328834590979877230
tmp(1)=2314885530823516726
tmp(2)=ny*long64(2)^32 + nx
tmp(3)=nzout
writeu,outunit,tmp

print,'nzin:',nzin,' nzout:',nzout

tmp=dblarr(22*nzin+11+92)
tmpout=dblarr(22*nzout+11+92)
readu,inunit,tmp                ; Ignore first record (signature)
for ix=0,nx-1 do for iy=0,ny-1 do begin
   readu,inunit,tmp
   zin=tmp(0:nzin-1)
   tauin=tmp(nzin:2*nzin-1)
   tmpout(0:nzout-1)=0.
   tmpout(nzout:2*nzout-1)=newscale
   for ivar=2,21 do begin
      yin=tmp(ivar*nzin:(ivar+1)*nzin-1)
      if (ivar ne 3 and ivar ne 4 and ivar ne 5) then begin ; normal variables
         yout=interpol(yin,tauin,newscale)
      endif else begin          ; log interpolation
         yout(0:nzout-1)=exp(interpol(alog(yin),tauin,newscale))
      endelse
      for ipoint=0,nzout-2 do begin
         ind=where( tauin ge newscale(ipoint) and tauin le newscale(ipoint+1))
         if (n_elements(ind) ge 2) then $
            yout[ipoint]=mean( yin[ind] )
      endfor
      if keyword_set(smoot) then yout=smooth(yout,2)
      tmpout(ivar*nzout:(ivar+1)*nzout-1)=yout
   endfor
   tmpout(22*nzout:*)=tmp(22*nzin:*)
   writeu,outunit,tmpout
endfor

free_lun,inunit
free_lun,outunit

end

