pro displayprof,mapin,i,q,u,v,si,sq,su,sv,_extra=extra
; NOTE: Uses tvframe.pro
;
; Shows a 2D image (e.g., a temperature map) and allows the user to 
; interactively move the mouse over it to display the Stokes profiles
; corresponding to each point in the image
; Right-click to exit
;
map=reform(mapin)
if (size(map))(0) ne 2 then begin
   print,'Map must be a 2D array'
   return
endif

if n_elements(i) gt 0 then begin
   if n_elements(q) eq 0 then q=i*0
   if n_elements(u) eq 0 then u=i*0
   if n_elements(v) eq 0 then v=i*0
endif

if n_elements(si) gt 0 then begin
   if n_elements(sq) eq 0 then sq=si*0
   if n_elements(su) eq 0 then su=si*0
   if n_elements(sv) eq 0 then sv=si*0
endif

window,/free
win1=!d.window
window,/free
win2=!d.window

oldpmulti=!p.multi

wset,win1
!p.multi=[0,1,1]
tvframe,map,_extra=extra

nx=(size(map))(1)
ny=(size(map))(2)

cursor,x,y,/nowait
oldx=0
oldy=0
while !mouse.button ne 4 do begin
   wset,win1
   !p.multi=[0,1,1]
   tvframe,map,_extra=extra,/noerase,/nodata

   cursor,x,y,/nowait,/norm
   x=(x-!x.s(0))/!x.s(1)
   y=(y-!y.s(0))/!y.s(1)
   x=long(x)
   y=long(y)

   if (x gt 0 and x lt nx-1 and y gt 0 and y lt ny-1 $
      and (x ne oldx or y ne oldy)) then begin
      wset,win2
      !p.multi=[0,2,2]
      plot,i(x,y,*)
      if n_elements(si) gt 0 then oplot,si(x,y,*),lin=2
      plot,q(x,y,*)
      if n_elements(sq) gt 0 then oplot,sq(x,y,*),lin=2
      plot,u(x,y,*)
      if n_elements(su) gt 0 then oplot,su(x,y,*),lin=2
      plot,v(x,y,*)
      if n_elements(sv) gt 0 then oplot,sv(x,y,*),lin=2
      !p.multi=[0,1,1]
      plot,indgen(10),xstyle=4,ystyle=4,/nodata,title=$
           '(x,y)=('+strtrim(string(x),2)+','+strtrim(string(y),2)+ $
           '). Pix value='+strtrim(string(map(x,y)),2),/noerase
      oldx=x
      oldy=y
   endif
endwhile

wdelete,win1
wdelete,win2
!p.multi=oldpmulti

end
