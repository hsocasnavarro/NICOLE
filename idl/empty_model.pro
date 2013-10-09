; Creates an IDL structure with an empty model

function empty_model,nx,ny,nz
  if n_elements(nx) eq 0 or n_elements(ny) eq 0 or n_elements(nz) eq 0 then begin
     print, 'Usage:'
     print, 'model=empty_model(nx,ny,nz)'
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
     keep_el_p: fltarr(nx,ny), keep_gas_p: fltarr(nx,ny), $
     keep_rho: fltarr(nx,ny), keep_nH: fltarr(nx,ny), $
     keep_nHminus: fltarr(nx,ny), keep_nHplus: fltarr(nx,ny), $
     keep_nH2: fltarr(nx,ny), keep_nh2plus: fltarr(nx,ny), $
     ffactor: fltarr(nx,ny), abundance: fltarr(nx,ny,92)}


  return, m

end
