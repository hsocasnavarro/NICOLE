; Replaces a column from a existing model structure m at coordinates ix,iy
; with column x2,y2 from model structure m2

pro model_column,m,ix,iy,m2,x2,y2

  m.tau(ix,iy,*)=m2.tau(x2,y2,*)
  m.z(ix,iy,*)=m2.z(x2,y2,*)
  m.t(ix,iy,*)=m2.t(x2,y2,*)
  m.v_los(ix,iy,*)=m2.v_los(x2,y2,*)
  m.el_p(ix,iy,*)=m2.el_p(x2,y2,*)
  m.gas_p(ix,iy,*)=m2.gas_p(x2,y2,*)
  m.rho(ix,iy,*)=m2.rho(x2,y2,*)
  m.v_mic(ix,iy,*)=m2.v_mic(x2,y2,*)
  m.b_los_x(ix,iy,*)=m2.b_los_x(x2,y2,*)
  m.b_los_y(ix,iy,*)=m2.b_los_y(x2,y2,*)
  m.b_los_z(ix,iy,*)=m2.b_los_z(x2,y2,*)
  m.b_x(ix,iy,*)=m2.b_x(x2,y2,*)
  m.b_y(ix,iy,*)=m2.b_y(x2,y2,*)
  m.b_z(ix,iy,*)=m2.b_z(x2,y2,*)
  m.v_x(ix,iy,*)=m2.v_x(x2,y2,*)
  m.v_y(ix,iy,*)=m2.v_y(x2,y2,*)
  m.v_z(ix,iy,*)=m2.v_z(x2,y2,*)
  m.nH(ix,iy,*)=m2.nH(x2,y2,*)
  m.nHplus(ix,iy,*)=m2.nHplus(x2,y2,*)
  m.nHminus(ix,iy,*)=m2.nHminus(x2,y2,*)
  m.nH2(ix,iy,*)=m2.nH2(x2,y2,*)
  m.nH2plus(ix,iy,*)=m2.nH2plus(x2,y2,*)
  m.v_mac(ix,iy)=m2.v_mac(x2,y2)
  m.stray_frac(ix,iy)=m2.stray_frac(x2,y2)
  m.chrom_x(ix,iy)=m2.chrom_x(x2,y2)
  m.chrom_y(ix,iy)=m2.chrom_y(x2,y2)
  m.ffactor(ix,iy)=m2.ffactor(x2,y2)
  m.abundance(ix,iy,*)=m2.abundance(x2,y2,*)
  m.keep_gas_p(ix,iy)=m2.keep_gas_p(x2,y2)
  m.keep_rho(ix,iy)=m2.keep_rho(x2,y2)
  m.keep_nH(ix,iy)=m2.keep_nH(x2,y2)
  m.keep_nHplus(ix,iy)=m2.keep_nHplus(x2,y2)
  m.keep_nHminus(ix,iy)=m2.keep_nHminus(x2,y2)
  m.keep_nH2(ix,iy)=m2.keep_nH2(x2,y2)
  m.keep_nH2plus(ix,iy)=m2.keep_nH2plus(x2,y2)

return

end

