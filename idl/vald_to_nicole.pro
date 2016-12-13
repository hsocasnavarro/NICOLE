; This procedure takes a text file with the long-format output from
; the VAL-D database and extracts the spectral line information in 
; the format of NICOLE's LINES file

pro vald_to_nicole,filein=filein

openr,lun,filein,/get_lun
str=''
labels=['']
duplicated=['']
nlines=0
while (not eof(lun)) do begin
   readf,lun,str
   str=strtrim(str,2)
   str2=''
   if strmid(str,0,1) eq "'" then begin ; this is a line to parse
      readf,lun,str2 ; Clear following line
; Parse data
      nlines=nlines+1
      element=strmid(str,1,2)
      element=strtrim(element,2)
      ionstage=strmid(str,4,1)
      ionstage=strtrim(ionstage,2)
      wlength=strmid(str,7,16)
      wlength=strtrim(wlength,2)
      loggf=strmid(str,24,8)
      loggf=strtrim(loggf,2)
      pot=strmid(str,33,8)
      pot=strtrim(pot,2)
      jlow=strmid(str,42,5)
      jlow=strtrim(jlow,2)
      jup=strmid(str,57,5)
      jup=strtrim(jup,2)
      landelow=strmid(str,62,6)
      landelow=strtrim(landelow,2)
      landeup=strmid(str,69,6)
      landeup=strtrim(landeup,2)
      gamma_r=strmid(str,84,6)
      gamma_r=strtrim(gamma_r,2)
      gamma_stk=strmid(str,91,6)
      gamma_stk=strtrim(gamma_stk,2)
      gamma_vdw=strmid(str,98,6)
      gamma_vdw=strtrim(gamma_vdw,2) ;
; labels
      newlabel=element+ionstage+' '+strtrim(string(nint(wlength)),2)
      if (where(newlabel eq labels))(0) ne -1 then duplicated=[duplicated,newlabel]
         
      labels=[labels,newlabel]
; Print line section
      print,'['+newlabel+']'
      print,'Element= '+element
      print,'Ionization stage= '+ionstage
      print,'Wavelength= '+wlength
      print,'Excitation potential= '+pot
      print,'Log(gf)= '+loggf
      print,'Term (lower)= 5P'+jlow+' # Warning! Mock value'
      print,'Term (upper)= 5D'+jup+' # Warning! Mock value'
      if abs(gamma_r) gt 0.0001 or abs(gamma_stk) gt 0.0001 or abs(gamma_vdw) gt 0.0001 then begin
         print,'Collisions= 3 # Default '
         if abs(gamma_r) gt 0.0001 then print,'Gamma Radiative=',10^float(gamma_r)/1d8
         if abs(gamma_stk) gt 0.0001 then print,'Gamma Stark=', $
            10^float(gamma_stk)/1d8*1e12
         if abs(gamma_vdw) gt 0.0001 then print,'Gamma van der Waals=', $
            10^float(gamma_vdw)/1d8*1e16
      endif else begin
         print,'Collisions= Unsold # Default '
      endelse
      print,'Damping enhancement= 1.0 # Default '
      print,'Width= 2.0 # Default '
      print,'Mode= LTE # Default '
      print,''
   endif
endwhile

if n_elements(duplicated) gt 1 then begin
   print,''
   print ,'Warning! The following labels are duplicated'
   print,''
   for ii=1,n_elements(duplicated)-1 do print,'  '+duplicated(ii)
   print,''
endif
free_lun,lun


end
