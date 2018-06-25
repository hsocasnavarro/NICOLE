; This procedure creates a blank new model containing zeros in all
; arrays and variables except abundances, which are set to the 
; Grevesse & Sauval (1998) values

function new_model,nx,ny,nz


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
     ffactor: fltarr(nx,ny), abundance: fltarr(nx,ny,92), $
     chrom_x: fltarr(nx,ny), chrom_y: fltarr(nx,ny)}
  for ix=0,nx-1 do for iy=0,ny-1 do $
     m.abundance[ix,iy,*]=[ 12.00,10.93,1.10,1.40,2.55,8.52,7.92,8.83,4.56,8.08,6.33,7.58,$
                6.47,7.55,5.45,7.33,5.5,6.40,5.12,6.36,3.17,5.02,4.00,5.67,$
                5.39,7.50,4.92,6.25,4.21,4.60,2.88,3.41,-10,-10,-10,-10,$
                2.60,2.97,2.24,2.60,1.42,$
                1.92,-10,1.84,1.12,1.69,0.94,1.77,1.66,2.0,1.0,-10,-10,-10,-10,$
                2.13,1.17,1.58,0.71,1.50,-10,1.01,0.51,1.12,-0.1,1.14,0.26,$
                0.93,0.00,1.08,0.06,0.88,-10,1.11,-10,1.45,1.35,1.8,1.01,-10,$
                     0.9,1.95,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10]


  return, m

end

