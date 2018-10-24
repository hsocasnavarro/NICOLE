; Function to read a model atmosphere from NICOLE in binary format
;
; Inputs:
;   filename: Model file in NICOLE native format
;   outfile (optional): if present, create an IDL savefile with the
;     individual variables (not in the form of a structure), in a
;     form that NICOLE can accept as input
;
; Outputs:
;   structure with all the model variables. Each variable has
;      dimensions (nx,ny,nz)
;   nz (optional): Number of grid points in the vertical direction
;   npix (optional): Number of models contained in the file
;   nx, ny (optional): Horizontal dimension of models contained in the file
;   chisq (optional): Array of dimensions nx,ny with the values of
;      Chi-square. Note that this value must be set to an existing
;      variable in order to read the Chi-square
;   formatversion (optional): String with the version of the file
;      format (currently either '1.6' or '2.3')
;   mask (optional): If the file maskinvert.dat exists, it is read and
;      the mask returned in this array
;   chrom (optional): If this keyword is present, add chromospheric
;      temperature increase as defined by model parameters chrom_x
;      and chrom_y to temperature stratification


function read_model,filename,outfile=outfile,nz=nz,npix=npix,chisq=chisq, $
                    nx=nx,ny=ny,formatversion=formatversion,mask=mask, $
                    chrom=chrom

  ; Check file signature
  openr,inunit,filename,/get_lun,/swap_if_big_endian
  tmp=lon64arr(13)
  readu,inunit,tmp
  free_lun,inunit
  nz=tmp(3)
  formatversion='xxx'
  if tmp(0) eq 3328553116003166574 and tmp(1) eq 2314885530823516726 then begin
     formatversion='1.6'
     npix=tmp(2)
     nx=1l
     ny=npix
  endif
  if tmp(0) eq 3328834590979877230 and tmp(1) eq 2314885530823516723 then begin
     formatversion='2.3'
     openr,inunit,filename,/get_lun,/swap_if_big_endian
     tmp=lonarr(26)
     readu,inunit,tmp
     free_lun,inunit
     nx=tmp(4)
     ny=tmp(5)
     npix=long(nx*ny)
  endif
  if tmp(0) eq 3328834590979877230 and tmp(1) eq 2314885530823516726 then begin
     formatversion='2.6'
     openr,inunit,filename,/get_lun,/swap_if_big_endian
     tmp=lonarr(26)
     readu,inunit,tmp
     free_lun,inunit
     nx=tmp(4)
     ny=tmp(5)
     npix=long(nx*ny)
  endif
  if tmp(0) eq 4049129056382445934 and tmp(1) eq 2314885530823504944 then begin
     formatversion='18.04'
     openr,inunit,filename,/get_lun,/swap_if_big_endian
     tmp=lonarr(26)
     readu,inunit,tmp
     free_lun,inunit
     nx=tmp(4)
     ny=tmp(5)
     npix=long(nx*ny)
  endif
  if formatversion eq 'xxx' then begin
     print,'File is not a valid model file:',filename
     stop
  endif

  m={z:fltarr(nx,ny,nz), tau:fltarr(nx,ny,nz),  t:fltarr(nx,ny,nz), $
     gas_p:fltarr(nx,ny,nz), rho:fltarr(nx,ny,nz), el_p:fltarr(nx,ny,nz), $
     v_los:fltarr(nx,ny,nz), v_mic:fltarr(nx,ny,nz), $
     b_los_z:fltarr(nx,ny,nz), b_los_x:fltarr(nx,ny,nz), $
     b_los_y:fltarr(nx,ny,nz), b_x:fltarr(nx,ny,nz), $
     b_y:fltarr(nx,ny,nz), b_z:fltarr(nx,ny,nz), $
     v_x:fltarr(nx,ny,nz), v_y:fltarr(nx,ny,nz), v_z:fltarr(nx,ny,nz), $
     nH: fltarr(nx,ny,nz), nHminus: fltarr(nx,ny,nz), nHplus: fltarr(nx,ny,nz),$
     nH2: fltarr(nx,ny,nz), nh2plus: fltarr(nx,ny,nz), $
     v_mac: fltarr(nx,ny), stray_frac: fltarr(nx,ny),$
     chrom_x: fltarr(nx,ny), chrom_y: fltarr(nx,ny),$
     keep_el_p: fltarr(nx,ny), keep_gas_p: fltarr(nx,ny), $
     keep_rho: fltarr(nx,ny), keep_nH: fltarr(nx,ny), $
     keep_nHminus: fltarr(nx,ny), keep_nHplus: fltarr(nx,ny), $
     keep_nH2: fltarr(nx,ny), keep_nh2plus: fltarr(nx,ny), $
     ffactor: fltarr(nx,ny), abundance: fltarr(nx,ny,92)}

  openr,inunit,filename,/get_lun,/swap_if_big_endian

  if formatversion eq '2.3' then begin
     tmp=dblarr(17*nz+3)
     readu,inunit,tmp           ; Ignore first record (signature)
     ON_IOERROR,ers
     for ix=0l,nx-1 do for iy=0l,ny-1 do begin
        if not eof(inunit) then readu,inunit,tmp
        ers: m.z(ix,iy,*)=tmp(0*nz:(0+1)*nz-1)
        m.tau(ix,iy,*)=tmp(1*nz:(1+1)*nz-1)
        m.t(ix,iy,*)=tmp(2*nz:(2+1)*nz-1)
        m.gas_p(ix,iy,*)=tmp(3*nz:(3+1)*nz-1)
        m.rho(ix,iy,*)=tmp(4*nz:(4+1)*nz-1)
        m.el_p(ix,iy,*)=tmp(5*nz:(5+1)*nz-1)
        m.v_los(ix,iy,*)=tmp(6*nz:(6+1)*nz-1)
        m.v_mic(ix,iy,*)=tmp(7*nz:(7+1)*nz-1)
        m.b_los_z(ix,iy,*)=tmp(8*nz:(8+1)*nz-1)
        m.b_los_x(ix,iy,*)=tmp(9*nz:(9+1)*nz-1)
        m.b_los_y(ix,iy,*)=tmp(10*nz:(10+1)*nz-1)
        m.b_x(ix,iy,*)=tmp(11*nz:(11+1)*nz-1)
        m.b_y(ix,iy,*)=tmp(12*nz:(12+1)*nz-1)
        m.b_z(ix,iy,*)=tmp(13*nz:(13+1)*nz-1)
        m.v_x(ix,iy,*)=tmp(14*nz:(14+1)*nz-1)
        m.v_y(ix,iy,*)=tmp(15*nz:(15+1)*nz-1)
        m.v_z(ix,iy,*)=tmp(16*nz:(16+1)*nz-1)
        m.v_mac(ix,iy)=tmp(17*nz)
        m.stray_frac(ix,iy)=tmp(17*nz+1)
        m.ffactor(ix,iy)=tmp(17*nz+2)
        m.chrom_x(ix,iy)=0.
        m.chrom_y(ix,iy)=0.
     endfor
     ON_IOERROR,NULL
 endif

  if formatversion eq '2.6' then begin
     tmp=dblarr(22*nz+11+92)
     readu,inunit,tmp           ; Ignore first record (signature)
     ON_IOERROR,ers3
     for ix=0l,nx-1 do for iy=0l,ny-1 do begin
        if not eof(inunit) then readu,inunit,tmp
        ers3: m.z(ix,iy,*)=tmp(0*nz:(0+1)*nz-1)
        m.tau(ix,iy,*)=tmp(1*nz:(1+1)*nz-1)
        m.t(ix,iy,*)=tmp(2*nz:(2+1)*nz-1)
        m.gas_p(ix,iy,*)=tmp(3*nz:(3+1)*nz-1)
        m.rho(ix,iy,*)=tmp(4*nz:(4+1)*nz-1)
        m.el_p(ix,iy,*)=tmp(5*nz:(5+1)*nz-1)
        m.v_los(ix,iy,*)=tmp(6*nz:(6+1)*nz-1)
        m.v_mic(ix,iy,*)=tmp(7*nz:(7+1)*nz-1)
        m.b_los_z(ix,iy,*)=tmp(8*nz:(8+1)*nz-1)
        m.b_los_x(ix,iy,*)=tmp(9*nz:(9+1)*nz-1)
        m.b_los_y(ix,iy,*)=tmp(10*nz:(10+1)*nz-1)
        m.b_x(ix,iy,*)=tmp(11*nz:(11+1)*nz-1)
        m.b_y(ix,iy,*)=tmp(12*nz:(12+1)*nz-1)
        m.b_z(ix,iy,*)=tmp(13*nz:(13+1)*nz-1)
        m.v_x(ix,iy,*)=tmp(14*nz:(14+1)*nz-1)
        m.v_y(ix,iy,*)=tmp(15*nz:(15+1)*nz-1)
        m.v_z(ix,iy,*)=tmp(16*nz:(16+1)*nz-1)
        m.nH(ix,iy,*)=tmp(17*nz:(17+1)*nz-1)
        m.nHminus(ix,iy,*)=tmp(18*nz:(18+1)*nz-1)
        m.nHplus(ix,iy,*)=tmp(19*nz:(19+1)*nz-1)
        m.nH2(ix,iy,*)=tmp(20*nz:(20+1)*nz-1)
        m.nh2plus(ix,iy,*)=tmp(21*nz:(21+1)*nz-1)
        m.v_mac(ix,iy)=tmp(22*nz)
        m.stray_frac(ix,iy)=tmp(22*nz+1)
        m.ffactor(ix,iy)=tmp(22*nz+2)
        m.keep_el_p(ix,iy)=tmp(22*nz+3)
        m.keep_gas_p(ix,iy)=tmp(22*nz+4)
        m.keep_rho(ix,iy)=tmp(22*nz+5)
        m.keep_nH(ix,iy)=tmp(22*nz+6)
        m.keep_nHminus(ix,iy)=tmp(22*nz+7)
        m.keep_nHplus(ix,iy)=tmp(22*nz+8)
        m.keep_nH2(ix,iy)=tmp(22*nz+9)
        m.keep_nh2plus(ix,iy)=tmp(22*nz+10)
        m.abundance(ix,iy,0:91)=tmp(22*nz+10+1:22*nz+10+92)
        m.chrom_x(ix,iy)=0.
        m.chrom_y(ix,iy)=0.
     endfor
     ON_IOERROR,NULL
 endif

  if formatversion eq '18.04' then begin
     tmp=dblarr(22*nz+13+92)
     readu,inunit,tmp           ; Ignore first record (signature)
     ON_IOERROR,ers4
     for ix=0l,nx-1 do for iy=0l,ny-1 do begin
        if not eof(inunit) then readu,inunit,tmp
        ers4: m.z(ix,iy,*)=tmp(0*nz:(0+1)*nz-1)
        m.tau(ix,iy,*)=tmp(1*nz:(1+1)*nz-1)
        m.t(ix,iy,*)=tmp(2*nz:(2+1)*nz-1)
        m.gas_p(ix,iy,*)=tmp(3*nz:(3+1)*nz-1)
        m.rho(ix,iy,*)=tmp(4*nz:(4+1)*nz-1)
        m.el_p(ix,iy,*)=tmp(5*nz:(5+1)*nz-1)
        m.v_los(ix,iy,*)=tmp(6*nz:(6+1)*nz-1)
        m.v_mic(ix,iy,*)=tmp(7*nz:(7+1)*nz-1)
        m.b_los_z(ix,iy,*)=tmp(8*nz:(8+1)*nz-1)
        m.b_los_x(ix,iy,*)=tmp(9*nz:(9+1)*nz-1)
        m.b_los_y(ix,iy,*)=tmp(10*nz:(10+1)*nz-1)
        m.b_x(ix,iy,*)=tmp(11*nz:(11+1)*nz-1)
        m.b_y(ix,iy,*)=tmp(12*nz:(12+1)*nz-1)
        m.b_z(ix,iy,*)=tmp(13*nz:(13+1)*nz-1)
        m.v_x(ix,iy,*)=tmp(14*nz:(14+1)*nz-1)
        m.v_y(ix,iy,*)=tmp(15*nz:(15+1)*nz-1)
        m.v_z(ix,iy,*)=tmp(16*nz:(16+1)*nz-1)
        m.nH(ix,iy,*)=tmp(17*nz:(17+1)*nz-1)
        m.nHminus(ix,iy,*)=tmp(18*nz:(18+1)*nz-1)
        m.nHplus(ix,iy,*)=tmp(19*nz:(19+1)*nz-1)
        m.nH2(ix,iy,*)=tmp(20*nz:(20+1)*nz-1)
        m.nh2plus(ix,iy,*)=tmp(21*nz:(21+1)*nz-1)
        m.v_mac(ix,iy)=tmp(22*nz)
        m.stray_frac(ix,iy)=tmp(22*nz+1)
        m.ffactor(ix,iy)=tmp(22*nz+2)
        m.keep_el_p(ix,iy)=tmp(22*nz+3)
        m.keep_gas_p(ix,iy)=tmp(22*nz+4)
        m.keep_rho(ix,iy)=tmp(22*nz+5)
        m.keep_nH(ix,iy)=tmp(22*nz+6)
        m.keep_nHminus(ix,iy)=tmp(22*nz+7)
        m.keep_nHplus(ix,iy)=tmp(22*nz+8)
        m.keep_nH2(ix,iy)=tmp(22*nz+9)
        m.keep_nh2plus(ix,iy)=tmp(22*nz+10)
        m.abundance(ix,iy,0:91)=tmp(22*nz+12+1:22*nz+12+92)
        m.chrom_x(ix,iy)=tmp(22*nz+11)
        m.chrom_y(ix,iy)=tmp(22*nz+12)
     endfor
     ON_IOERROR,NULL
 endif

  if formatversion eq '1.6' then begin
     tmp=dblarr(13*nz+3)
     readu,inunit,tmp           ; Ignore first record (signature)
     ON_IOERROR,ers2
     for ipix=0l,npix-1 do begin
        if not eof(inunit) then readu,inunit,tmp
        ers2: m.z(0,ipix,*)=tmp(0*nz:(0+1)*nz-1)
        m.tau(0,ipix,*)=tmp(1*nz:(1+1)*nz-1)
        m.t(0,ipix,*)=tmp(2*nz:(2+1)*nz-1)
        m.gas_p(0,ipix,*)=tmp(3*nz:(3+1)*nz-1)
        m.rho(0,ipix,*)=tmp(4*nz:(4+1)*nz-1)
        m.el_p(0,ipix,*)=tmp(5*nz:(5+1)*nz-1)
        m.v_los(0,ipix,*)=tmp(6*nz:(6+1)*nz-1)
        m.v_mic(0,ipix,*)=tmp(7*nz:(7+1)*nz-1)
        m.b_los_z(0,ipix,*)=tmp(8*nz:(8+1)*nz-1)
        m.b_los_x(0,ipix,*)=tmp(9*nz:(9+1)*nz-1)
        m.b_los_y(0,ipix,*)=tmp(10*nz:(10+1)*nz-1)
;        m.b_local_inc(0,ipix,*)=tmp(11*nz:(11+1)*nz-1)
;        m.b_local_azi(0,ipix,*)=tmp(12*nz:(12+1)*nz-1)
        m.v_mac(0,ipix)=tmp(13*nz)
        m.stray_frac(0,ipix)=tmp(13*nz+1)
        m.ffactor(0,ipix)=tmp(13*nz+2)
        m.chrom_x[0,ipix]=0.
        m.chrom_y[0,ipix]=0.
     endfor
     ON_IOERROR,NULL
 endif

  free_lun,inunit

  if n_elements(chisq) gt 0 then begin
     filename=''
     if file_test('Chisq.dat') then filename='Chisq.dat'
     if (file_test('Chisq.dat_1')) then begin
        filename='Chisq.dat_1'
        i=1
        while (file_test('Chisq.dat_'+string(i,format='(i0)'))) do begin
           filename='Chisq.dat_'+string(i,format='(i0)')
           i=i+1
        endwhile
     endif
     chisq=fltarr(nx,ny)
     openr,iunit,'Chisq.dat_1',/get_lun,/swap_if_big_endian
     tmpchisq=0.d0
     for ix=0l,nx-1 do for iy=0l,ny-1 do begin
        if not eof(iunit) then readu,iunit,tmpchisq
        chisq(ix,iy)=tmpchisq
     endfor
     free_lun,inunit
  endif

  if keyword_set(outfile) then begin
     z=reform(m.z,nx,ny,nz)
     tau=reform(m.tau,nx,ny,nz)
     t=reform(m.t,nx,ny,nz)
     gas_p=reform(m.gas_p,nx,ny,nz)
     rho=reform(m.rho,nx,ny,nz)
     el_p=reform(m.el_p,nx,ny,nz)
     v_los=reform(m.v_los,nx,ny,nz)
     v_mic=reform(m.v_mic,nx,ny,nz)
     b_los_z=reform(m.b_los_z,nx,ny,nz)
     b_los_x=reform(m.b_los_x,nx,ny,nz)
     b_los_y=reform(m.b_los_y,nx,ny,nz)
     b_x=reform(m.b_x,nx,ny,nz)
     b_y=reform(m.b_y,nx,ny,nz)
     b_z=reform(m.b_z,nx,ny,nz)
     v_x=reform(m.v_x,nx,ny,nz)
     v_y=reform(m.v_y,nx,ny,nz)
     v_z=reform(m.v_z,nx,ny,nz)
     v_mac=reform(m.v_mac,nx,ny)
     stray_frac=reform(m.stray_frac,nx,ny)
     ffactor=reform(m.ffactor,nx,ny)
     save,file=outfile,z,tau,t,gas_p,rho,el_p,v_los,v_mic,b_los_z, $
          b_los_x,b_los_y,b_x,b_y,b_z,v_x,v_y,v_z, $
          v_mac,stray_frac,ffactor
  endif

; Read mask?
  mask=fltarr(nx,ny)
  ind=strpos(filename,'/',/reverse_search)
  if ind eq -1 then $
     path='./' $
  else $
     path=strmid(filename,0,ind+1)

  if (file_test(path+'maskinvert.dat_1')) then begin
     openr,inunit,path+'maskinvert.dat_1',/get_lun,/swap_if_big_endian
     tmp=0d0
     for ix=0l,nx-1 do for iy=0l,ny-1 do begin
        readu,inunit,tmp
        mask(ix,iy)=tmp
     endfor
     free_lun,inunit
  endif

  if (keyword_set(chrom)) then begin
     nx=(size(m.t))(1)
     ny=(size(m.t))(2)
     for ix=0,nx-1 do for iy=0,ny-1 do for ip=0,nz-1 do begin
        if m.tau[ix,iy,ip] lt m.chrom_x[ix,iy]-.1 then $
           m.t[ix,iy,ip]=m.t[ix,iy,ip]+m.chrom_y[ix,iy] $
        else if m.tau[ix,iy,ip] lt m.chrom_x[ix,iy] + 3 then $
           m.t[ix,iy,ip]=m.t[ix,iy,ip]+m.chrom_y[ix,iy]*exp(  -((m.tau[ix,iy,ip]-m.chrom_x[ix,iy])/(0.5))^2 )
     endfor
  endif
  
  return,m
end
