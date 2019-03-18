; Take a model in NICOLE binary format and interpolate it to a new
; tau grid

pro interpolate_model,filein,fileout,newscale

if filein eq fileout then begin
   print,'Error. Input and output files must be different'
   stop
endif

nzout=n_elements(newscale)
modelin=read_model(filein)
nx=(size(modelin.t))(1)
ny=(size(modelin.t))(2)
nzin=(size(modelin.t))(3)

print,'nx,ny=',nx,ny

print,'nzin:',nzin,' nzout:',nzout

modelout=new_model(nx,ny,nzout)
for ix=0,nx-1 do for iy=0,ny-1 do begin
   for iz=0, nzout-1 do $
      modelout.tau[*,*,iz]=newscale[iz]
   tauin=modelin.tau[ix,iy,*]
   modelout.t[ix,iy,*]=interpol(modelin.t[ix,iy,*],tauin,newscale)
   modelout.v_mic[ix,iy,*]=interpol(modelin.v_mic[ix,iy,*],tauin,newscale)
   modelout.v_los[ix,iy,*]=interpol(modelin.v_los[ix,iy,*],tauin,newscale)
   modelout.b_los_z[ix,iy,*]=interpol(modelin.b_los_z[ix,iy,*],tauin,newscale)
   modelout.b_los_x[ix,iy,*]=interpol(modelin.b_los_x[ix,iy,*],tauin,newscale)
   modelout.b_los_y[ix,iy,*]=interpol(modelin.b_los_y[ix,iy,*],tauin,newscale)
   modelout.el_p[ix,iy,*]=exp(interpol(alog(modelin.el_p[ix,iy,*]),tauin,newscale))
   modelout.gas_p[ix,iy,*]=exp(interpol(alog(modelin.gas_p[ix,iy,*]),tauin,newscale))
   modelout.rho[ix,iy,*]=exp(interpol(alog(modelin.rho[ix,iy,*]),tauin,newscale))
   modelout.v_mac[ix,iy]=modelin.v_mac[ix,iy]
   modelout.stray_frac[ix,iy]=modelin.stray_frac[ix,iy]
   modelout.ffactor[ix,iy]=modelin.ffactor[ix,iy]
   modelout.abundance[ix,iy]=modelin.abundance[ix,iy]
   modelout.chrom_x[ix,iy]=modelin.chrom_x[ix,iy]
   modelout.chrom_y[ix,iy]=modelin.chrom_y[ix,iy]
endfor

idl_to_nicole,file=fileout,model=modelout

end

