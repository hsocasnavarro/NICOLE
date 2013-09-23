tmp=fltarr(95)
openr,1,'/d/navarro/atmospheric_models/hsra/hsra_gasp.txt' & readf,1,tmp & close,1
pg=tmp
openr,1,'/d/navarro/atmospheric_models/hsra/hsra_H_ion.txt' & readf,1,tmp & close,1
hion=tmp
openr,1,'/d/navarro/atmospheric_models/hsra/hsra_rho.txt' & readf,1,tmp & close,1
rho=tmp
openr,1,'/d/navarro/atmospheric_models/hsra/hsra_z.txt' & readf,1,tmp & close,1
z=tmp


end
