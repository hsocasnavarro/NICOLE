; This procedure takes a text file with the long-format output from
; the VAL-D database and extracts the spectral line information in 
; the format of NICOLE's LINES file
;
; Using Extract ALL, short format with all boxes checked

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
; Parse data
      nlines=nlines+1
      strdata=strsplit(str,',',/extract)
      element=(strsplit(strdata[0],"' ",/extract))[0]
      ionstage=(strsplit(strdata[0],"' ",/extract))[1]
      wlength=(strdata[1])
      pot=(strdata[2])
      loggf=(strdata[3])
      gamma_r=(strdata[4])
      gamma_stk=(strdata[5])
      gamma_vdw=(strdata[6])
      lande=(strdata[7])
      jlow='0.'
      jup='0.'
; labels
      newlabel=element+ionstage+' '+strtrim(string(wlength,format='(f0.3)'),2)
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
