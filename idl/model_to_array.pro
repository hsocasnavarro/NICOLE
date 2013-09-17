; Function that takes a model either in the form of a structure or
; as a savefile (with named variables) and puts it in the form of a
; 13-column array (or optionally a 8-column array). If there is an
; array of models, keywords ix, iy (default to zero) point to the
; model to extract
; 

function model_to_array,struct=m,savefile=savefile,columns=columns,ix=ix, $
                        iy=iy,v_mac=v_mac,stray_frac=stray_frac,$
                        ffactor=ffactor

  if n_elements(ix) eq 0 then ix=0
  if n_elements(iy) eq 0 then iy=0
  if n_elements(columns) ne 0 then begin
     if columns ne 13 and columns ne 8 then begin
        print,'If columns keyword is present, it must be set to either 13 or 8'
        stop
     endif
  endif else columns=13

  if n_elements(m) ne 0 then begin
     nz=(size(m.t))(1)
     data=fltarr(columns,nz)
     if columns eq 13 then begin
        data(0,*)=m.z(ix,iy,*)
        data(1,*)=m.tau(ix,iy,*)
        data(2,*)=m.t(ix,iy,*)
        data(3,*)=m.gas_p(ix,iy,*)
        data(4,*)=m.rho(ix,iy,*)
        data(5,*)=m.el_p(ix,iy,*)
        data(6,*)=m.v_los(ix,iy,*)
        data(7,*)=m.v_mic(ix,iy,*)
        data(8,*)=m.b_long(ix,iy,*)
        data(9,*)=m.b_x(ix,iy,*)
        data(10,*)=m.b_y(ix,iy,*)
        data(11,*)=m.b_local_inc(ix,iy,*)
        data(12,*)=m.b_local_azi(ix,iy,*)
     endif else begin
        data(0,*)=m.tau(ix,iy,*)
        data(1,*)=m.t(ix,iy,*)
        data(2,*)=m.el_p(ix,iy,*)
        data(3,*)=m.v_mic(ix,iy,*)
        data(4,*)=sqrt(m.b_long(ix,iy,*)^2+m.b_x(ix,iy,*)^2+m.b_y(ix,iy,*)^2)
        data(5,*)=m.v_los(ix,iy,*)
        btan=sqrt(m.b_x(ix,iy,*)^2+m.b_y(ix,iy,*)^2)
        data(6,*)=abs(atan(btan,m.b_long(ix,iy,*)))
        data(7,*)=abs(atan(m.b_y(ix,iy,*),m.b_x(ix,iy,*)))
     endelse
     v_mac=m.v_mac
     stray_frac=m.stray_frac
     ffactor=m.ffactor
     return,data
  endif

  if n_elements(savefile) ne 0 then begin
     restore,savefile
     nz=(size(t))(2)
     data=fltarr(columns,nz)
     if columns eq 13 then begin
        data(0,*)=z(ix,iy,*)
        data(1,*)=tau(ix,iy,*)
        data(2,*)=t(ix,iy,*)
        data(3,*)=gas_p(ix,iy,*)
        data(4,*)=rho(ix,iy,*)
        data(5,*)=el_p(ix,iy,*)
        data(6,*)=v_los(ix,iy,*)
        data(7,*)=v_mic(ix,iy,*)
        data(8,*)=b_long(ix,iy,*)
        data(9,*)=b_x(ix,iy,*)
        data(10,*)=b_y(ix,iy,*)
        data(11,*)=b_local_inc(ix,iy,*)
        data(12,*)=b_local_azi(ix,iy,*)
     endif else begin
        data(0,*)=tau(ix,iy,*)
        data(1,*)=t(ix,iy,*)
        data(2,*)=el_p(ix,iy,*)
        data(3,*)=v_mic(ix,iy,*)
        data(4,*)=sqrt(b_long(ix,iy,*)^2+b_x(ix,iy,*)^2+b_y(ix,iy,*)^2)
        data(5,*)=v_los(ix,iy,*)
        btan=sqrt(b_x(ix,iy,*)^2+b_y(ix,iy,*)^2)
        data(6,*)=abs(atan(btan,b_long(ix,iy,*)))
        data(7,*)=abs(atan(b_y(ix,iy,*),b_x(ix,iy,*)))
     endelse
     return,data
  endif


end
