; Read debug files
; corenumber (defuault=0) is the number of core file to read (in case
; there are more than one in the current directory). Out output
; d, e and m are structures containing the debug information and the
; model atmosphere respectively. d has data stored in the core file
; whereas e has data in the various additional output files, such as
; Populations.dat or Cont_opacity.dat


pro read_debug,corenumber=corenumber,d,e,m

  if n_elements(corenumber) eq 0 then corenumber=1

  ndep=1
  nk=1
  ndata=1
  ntran=1
  tau5=fltarr(ndep)
  outputpop='F'
  outputcontop='F'
  outputnltesf='F'
  do_nlte='?'
  if file_test('debug.txt') then begin
     openr,unit,'debug.txt',/get_lun
     string=''
     readf,unit,string
     values=strsplit(string,/extract)
     outputpop=values(0)
     outputcontop=values(1)
     outputnltesf=values(2)
     readf,unit,ndep,ndata,ntran
     do_nlte=''
     readf,unit,do_nlte
     do_nlte=(strsplit(do_nlte,/extract))(0)
     if do_nlte eq 'T' then readf,unit,nk
     readf,unit,tau5
     free_lun,unit
     if ntran eq 0 then ntran=1
  endif

  maxnq=100
  e={mode:'',ndep:ndep,nk:nk,ndata:ndata,ntran:ntran,nlte:do_nlte,$
     nx:-1,ny:-1,ix:-1,iy:-1,irec:-1,icycle:-1, $
     tau5:tau5,n:fltarr(nk,ndep),nstar:fltarr(nk,ndep),$
     ngi:fltarr(ndep),ngj:fltarr(ndep),contop:fltarr(ndep),$
     contop5:fltarr(ndep),sf:fltarr(ntran,maxnq,ndep), $
     obs:fltarr(ndata)}

  tmp=dblarr(ndep)
  if outputpop eq 'T' then begin
     openu,unit,'Populations.dat',/swap_if_big,/get_lun
     if do_nlte eq 'T' then begin
        for i=1,nk do begin
           readu,unit,tmp
           e.n(i-1,*)=tmp
        endfor
        for i=1,nk do begin
           readu,unit,tmp
           e.nstar(i-1,*)=tmp
        endfor
     endif else begin
        readu,unit,tmp
        e.ngi=tmp
        readu,unit,tmp
        e.ngj=tmp
     endelse
     free_lun,unit
  endif

  if outputcontop eq 'T' then begin
     openu,unit,'Cont_opacity.dat',/swap_if_big,/get_lun
     tmp=dblarr(ndep)
     readu,unit,tmp
     e.contop=tmp
     readu,unit,tmp
     e.contop5=tmp
     free_lun,unit
  endif

  if outputpop eq 'T' then begin
     openu,unit,'NLTE_sf.dat',/swap_if_big,/get_lun
     readu,unit,tmp
     if max(tmp) gt  maxnq then begin
        print,"Can't read all wavelengths in source function file"
        print,"There is a transition with ",max(tmp)," wavelenghts."
        print,"Can only read up to maxnq=",maxnq
     endif
     e.sf=0.
     for i=1,ntran do begin
        for j=1, tmp[i] do begin
           readu,unit,tmp
           e.sf(i-1,j,*)=tmp
        endfor
     endfor
     free_lun,unit
  endif

  filename='core_'+string(corenumber,format='(i0)')+'.txt'
  if file_test(filename) then begin
     print,'Opening core file:'+filename
     openr,unit,filename,/get_lun
     str=''
     readf,unit,str
     if (str eq '  Inversion mode') then mode='i' else $
        if (str eq '  Synthesis mode') then mode='s' else $
           print,'Unknown mode'

     readf,unit,str
     readf,unit,ia,ib,ic,id
     ndep=ic & ndata=id & nx=ia & ny=ib
     readf,unit,str
     readf,unit,ia,ib,ic
     irec=ia & ix=ib & iy=ic
 
  d={mode:mode,ndep:ndep,nk:nk,ndata:ndata,ntran:ntran,nlte:do_nlte,$
     nx:nx,ny:ny,ix:ix,iy:iy,irec:irec,icycle:-1, $
     tau5:tau5,n:fltarr(nk,ndep),nstar:fltarr(nk,ndep),$
     ngi:fltarr(ndep),ngj:fltarr(ndep),contop:fltarr(ndep),$
     contop5:fltarr(ndep),sf:fltarr(ntran,ndep), $
     obs:fltarr(ndata)}
  d.nx=nx & d.ny=ny & d.ndep=ndep & d.ndata=ndata
  d.irec=irec & d.ix=ix & d.iy=iy

  m={z:fltarr(ndep), tau:fltarr(ndep),  t:fltarr(ndep), $
     gas_p:fltarr(ndep), rho:fltarr(ndep), el_p:fltarr(ndep), $
     v_los:fltarr(ndep), v_mic:fltarr(ndep), $
     b_los_z:fltarr(ndep), b_los_x:fltarr(ndep), $
     b_los_y:fltarr(ndep), b_x:fltarr(ndep), $
     b_y:fltarr(ndep), b_z:fltarr(ndep), $
     v_x:fltarr(ndep), v_y:fltarr(ndep), v_z:fltarr(ndep), $
     v_mac: 0., stray_frac: 0., ffactor: 0., keep_el_p:0., keep_gas_p: 0.,$
     keep_rho: 0., keep_nH: 0., keep_nHminus: 0., keep_nHplus: 0., $
     keep_nH2: 0., keep_nH2plus: 0., $
     abundance:fltarr(92), ane:fltarr(ndep),nh:fltarr(ndep), $
     nhplus:fltarr(ndep),nhminus:fltarr(ndep),nh2:fltarr(ndep), $
     nh2plus:fltarr(ndep)}
 
   if d.mode eq 'i' then begin
        readf,unit,str
        readf,unit,ia
        d.icycle=ia
        readf,unit,str
        tmp2=fltarr(d.ndata)
        readf,unit,tmp2
        d.obs=tmp2
     endif
     readf,unit,str
     tmp=fltarr(ndep)
     readf,unit,tmp
     m.z=tmp
     readf,unit,tmp
     m.tau=tmp
     readf,unit,tmp
     m.t=tmp
     readf,unit,tmp
     m.gas_p=tmp
     readf,unit,tmp
     m.rho=tmp
     readf,unit,tmp
     m.el_p=tmp
     readf,unit,tmp
     m.v_los=tmp
     readf,unit,tmp
     m.v_mic=tmp
     readf,unit,tmp
     m.b_los_z=tmp
     readf,unit,tmp
     m.b_los_x=tmp
     readf,unit,tmp
     m.b_los_y=tmp
     readf,unit,tmp
     m.b_x=tmp
     readf,unit,tmp
     m.b_y=tmp
     readf,unit,tmp
     m.b_z=tmp
     readf,unit,tmp
     m.v_x=tmp
     readf,unit,tmp
     m.v_y=tmp
     readf,unit,tmp
     m.v_z=tmp
     readf,unit,tmp
     m.ane=tmp
     readf,unit,tmp
     m.nh=tmp
     readf,unit,tmp
     m.nhplus=tmp
     readf,unit,tmp
     m.nhminus=tmp
     readf,unit,tmp
     m.nh2=tmp
     tmp=0.
     readf,unit,tmp
     m.v_mac=tmp
     readf,unit,tmp
     m.stray_frac=tmp
     readf,unit,tmp
     m.ffactor=tmp
     tmp=fltarr(92)
     readf,unit,tmp
     m.abundance=tmp
     free_lun,unit
  endif

end
