; Function to read a model atmosphere and cut a x,y range in it
;
;
; Inputs:
;   m_in: Model structure with the full model to cut (e.g., read using
;     read_model.pro)
;   xrange: 2-element array with the range of pixels to cut in the x-direction
;   yrange: 2-element array with the range of pixels to cut in the y-direction
;
; Outputs:
;  Function output is a model structure with the range selected


function cut_model,m_in,xrange,yrange
  nx0=(size(m_in.z))(1)
  ny0=(size(m_in.z))(2)
  nz0=(size(m_in.z))(3)
  nx=xrange(1)-xrange(0)+1
  ny=yrange(1)-yrange(0)+1
  nz=nz0
  indx=xrange(0)+indgen(nx)
  indy=yrange(0)+indgen(ny)
  indz=indgen(nz)
  m={z:fltarr(nx,ny,nz), tau:fltarr(nx,ny,nz),  t:fltarr(nx,ny,nz), $
     gas_p:fltarr(nx,ny,nz), rho:fltarr(nx,ny,nz), el_p:fltarr(nx,ny,nz), $
     v_los:fltarr(nx,ny,nz), v_mic:fltarr(nx,ny,nz), $
     b_los_z:fltarr(nx,ny,nz), b_los_x:fltarr(nx,ny,nz), $
     b_los_y:fltarr(nx,ny,nz), b_x:fltarr(nx,ny,nz), $
     b_y:fltarr(nx,ny,nz), b_z:fltarr(nx,ny,nz), $
     v_x:fltarr(nx,ny,nz), v_y:fltarr(nx,ny,nz), v_z:fltarr(nx,ny,nz), $
     v_mac: fltarr(nx,ny), stray_frac: fltarr(nx,ny),$
     ffactor: fltarr(nx,ny)}

  m.z=m_in.z(indx,indy,*)
  m.tau=m_in.tau(indx,indy,*)
  m.t=m_in.t(indx,indy,*)
  m.gas_p=m_in.gas_p(indx,indy,*)
  m.rho=m_in.rho(indx,indy,*)
  m.el_p=m_in.el_p(indx,indy,*)
  m.v_los=m_in.v_los(indx,indy,*)
  m.v_mic=m_in.v_mic(indx,indy,*)
  m.b_los_z=m_in.b_los_z(indx,indy,*)
  m.b_los_x=m_in.b_los_x(indx,indy,*)
  m.b_los_y=m_in.b_los_y(indx,indy,*)
  m.b_x=m_in.b_x(indx,indy,*)
  m.b_y=m_in.b_y(indx,indy,*)
  m.b_z=m_in.b_z(indx,indy,*)
  m.v_x=m_in.v_x(indx,indy,*)
  m.v_y=m_in.v_y(indx,indy,*)
  m.v_z=m_in.v_z(indx,indy,*)
  m.v_mac=m_in.v_mac(indx,indy,*)
  m.stray_frac=m_in.stray_frac(indx,indy,*)
  m.ffactor=m_in.ffactor(indx,indy,*) 

  return,m
end
