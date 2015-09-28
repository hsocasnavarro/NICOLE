; Function to read a spectral profile from NICOLE binary file
;
; The profile is returned in the variables StkI, StkQ, StkU and StkV,
;
; Note that if the array size is 1 in either dimension, IDL will
; sometimes drop that dimension when operating with the
; array. However, when when feeding NICOLE with IDL variables it's
; very important that both dimensions are present. Otherwise it will
; not be possible to determine which dimension (spectral or spatial)
; has only one element and the Python preprocessor will crash.
;
;
; Inputs:
;   filename: Profile file in NICOLE native format
;   outfile (optional): if present, create an IDL savefile with the
;     individual variables, in a form that NICOLE can accept as input
;
; Outputs:
;   StkI, StkQ, StkU, StkV: Stokes profiles. Note that I is the
;      function output whereas Q, U and V are output
;      arguments. Dimensions are nx,ny,nlambda
;   nlambda (optional): Number of wavelengths in the profile
;   npix (optional): Number of profiles contained in the file
;   nx, ny (optional): Horizontal dimension of models contained in the file
;   formatversion (optional): String with the version of the file
;      format (currently either '1.6' or '2.3')

function read_profile,filename,stkq,stku,stkv,outfile=outfile, $
                      nlambda=nl,npix=npix,nx=nx,ny=ny,$
                      formatversion=formatversion

  ; Check file signature
  openr,inunit,filename,/get_lun,/swap_if_big_endian
  tmp=lon64arr(4)
  readu,inunit,tmp
  free_lun,inunit
  nl=tmp(3)
  formatversion='xxx'
  if tmp(0) eq 3328553116003166574 and tmp(1) eq 2314885530823713334 then begin
     formatversion='1.6'
     npix=tmp(2)
     nx=1l
     ny=npix
  endif
  if tmp(0) eq 3328834590979877230 and tmp(1) eq 2314885530823713331 then begin
     formatversion='2.3'
     openr,inunit,filename,/get_lun,/swap_if_big_endian
     tmp=lonarr(26)
     readu,inunit,tmp
     free_lun,inunit
     nx=tmp(4)
     ny=tmp(5)
     npix=nx*ny
  endif
  if formatversion eq 'xxx' then begin
     print,'File is not a valid profile file:',filename
     stop
  endif

  stki=fltarr(nx,ny,nl)
  stki=reform(stki,nx,ny,nl)
  stkq=reform(stki,nx,ny,nl)
  stku=reform(stki,nx,ny,nl)
  stkv=reform(stki,nx,ny,nl)

  openr,inunit,filename,/get_lun,/swap_if_big_endian
  tmp=dblarr(4l*nl)
  readu,inunit,tmp ; Ignore first record (signature)
  ON_IOERROR,ers
  for ix=0l,nx-1 do for iy=0l,ny-1 do begin
     readu,inunit,tmp
ers:     stki(ix,iy,*)=tmp(lindgen(nl)*4l)
     stkq(ix,iy,*)=tmp(lindgen(nl)*4l+1l)
     stku(ix,iy,*)=tmp(lindgen(nl)*4l+2l)
     stkv(ix,iy,*)=tmp(lindgen(nl)*4l+3l)
  endfor
  free_lun,inunit

  if keyword_set(outfile) then begin
     save,file=outfile,stki,stkq,stku,stkv
  endif

  return,stki
end
