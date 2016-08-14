; Function to read the wavelength grid used in the last calculation
; It is read from __input.dat_1
;
; Inputs:
;    region (optional): If present, set to an integer with the desired
;      wavelength region to read. Other regions will be ignored
;    index (optional): If present, an array is returned with the
;      region of each wavelength
;
; Outputs:
;    wavelengths in A
;
;
function read_wavelengths,region=region,index=index,verbose=verbose

  str=''
  openr,iunit,'__input.dat_1',/get_lun
  while (strpos(str,'! nregions') eq -1) do readf,iunit,str
  reads,str,nregions
  if (keyword_set(verbose)) then print,'Number of regions:',nregions
  lambda=[0.]
  index=[0]
  for iregion=0,nregions-1 do begin
     readf,iunit,str
     if (keyword_set(verbose)) then print,str
     reads,str,lambda0,step,nlambda
     readit=0
     if not keyword_set(region) then readit=1
     if keyword_set(region) then if region eq iregion then readit=1
     if readit eq 1 then begin
        lambda=[lambda,findgen(nlambda)*step+lambda0]
        index=[index,intarr(nlambda)+iregion]
     endif
  endfor
  if (n_elements(lambda) gt 1) then begin
     index=index[1:*]
     lambda=lambda[1:*]
  endif

  return,lambda
  free_lun,iunit

end

