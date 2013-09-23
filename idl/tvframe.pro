function tvframe_minmax,array,NAN=nan
;+
; NAME:
;      MINMAX
; PURPOSE:
;      Return a 2 element array giving the minimum and maximum of a vector
; EXPLANATION:
;      Using MINMAX() is faster than doing a separate MAX and MIN.
;
; CALLING SEQUENCE:
;      value = minmax( array )
; INPUTS:
;      array - an IDL numeric scalar, vector or array.
;
; OUTPUTS:
;      value = a two element vector, 
;            value(0) = minimum value of array
;            value(1) = maximum value of array
;
; KEYWORDS:
;      NAN   - Set this keyword to cause the routine to check for occurrences of
;            the IEEE floating-point value NaN in the input data.  Elements with
;            the value NaN are treated as missing data.
;
; EXAMPLE:
;      Print the minimum and maximum of an image array, im
; 
;            IDL> print, minmax( im )
;
; PROCEDURE:
;      The MIN function is used with the MAX keyword
;
; REVISION HISTORY:
;      Written W. Landsman                January, 1990
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Added NaN keyword.      M. Buie       June 1998
;-
 On_error,2
 amin = min( array, MAX = amax, NAN=nan)
 return, [ amin, amax ]
 end
;+
; NAME:
;        TVFRAME
; PURPOSE:
;        Display an image with axis.
; CATEGORY:
;        Image dispaly
; CALLING SEQUENCE:
;        TVBOX,data
; INPUTS:
;        data = Two dimensional array of numerical type
; KEYWORD PARAMETERS:
;     XRANGE   : Array with at least 2 elements giving the range for
;                the annotation of the x-axis.
;                Only min(XRANGE) and max(XRANGE) are used.
;                Default is annotation with pixel indices.
;                (should not be confused with the standard plotting keword)
;     YRANGE   : Same as XRANGE but for y-axis.
;     POSITION : (output) position parameter for the axis in device units,
;                may be used in subsequent calls to PLOT or CONTOUR.
;                Example : TVBOX,a,/aspect,position=p
;                          CONTOUR,a,nlev=10,xsty=5,ysty=5,pos=p,/dev,/noeras
;                          (no keyword like /ASPECT exists for CONTOUR)
;     /sample  : If set, nearest neighbourhood method is used in
;                resizing the data array.
;                Default : averaging or linear interpolation.
;                (Corresponds to the SAMPLE-keyword of the REBIN function)
;     /center  : If set, tickmarks are placed at the center of each pixel.
;                ( example: tvframe,randomu(seed,20,20),/sam [,/center] )
;     /aspect  : If set, the aspect ratio of the data array is preserved
;                (quadratic pixel shapes).
;     /noscale : If set, data is not scaled. (much like TV compared to TVSCL)
;     /bar     : If set, an intensity bar is displayed.
;     /bbar    : Same as bar, only the image is not compressed horizontally
;                    to leave space for the bar
;     BRANGE   : Range for annotation of the intensity bar.
;                Default is the range of data.
;     BTITLE   : Title of the intensity bar.
;     BTICKS,BMINOR : Control the number of tickmarks for the intensity bar.
; Standard plotting keywords :
;      TITLE,SUBTITLE,XTITLE,YTITLE,TICKLEN,CHARSIZE,XSTYLE,YSTYLE,
;      XTICKS,YTICKS,XMINOR,YMINOR
;                Normally they are just passed through,
;                with the following exceptions :
;                the default value for TICKLEN is -.01 instead of .02 .
;                bit 0 of X/Y-STYLE is allways set (exact axis range).
;                bit 1 of X/Y-STYLE is  never  set (no extended axis range).
; OUTPUTS:
;        POSITION of the plot frame.
; COMMON BLOCKS:
;        none.
; SIDE EFFECTS:
;        none.
; RESTRICTIONS:
;        none.
; PROCEDURE:
;        The image is scaled to the size of the axis-box via the
;        RESIZE-function, and then displayed using TV or TVSCL.
; MODIFICATION HISTORY:
;        Written, A. Welz, Univ. Wuerzburg, Germany, Feb. 1991
;        Extended (/bar,xticks,..) by A.W. June 1992
;-
;
pro tvframe , ia    $
    , sample=sample , center=center , aspect=aspect , noscale=noscale    $
    , POSITION=POSITION    $
    , XRANGE=XRANGI , YRANGE=YRANGI   $
    , TITLE=TITLE , XTITLE=XTITLE , YTITLE=YTITLE , SUBTITLE=SUBTITLE   $
    , TICKLEN=TICKLEN , CHARSIZE=CHARSIZE   $
    , XTICKS=XTICKS , YTICKS=YTICKS , XMINOR=XMINOR , YMINOR=YMINOR $
    , XSTYLE=XSTYLI , YSTYLE=YSTYLI   $
    , bar=bar , bbar=bbar , BTITLE=BTITLE , BRANGE=BRANGE  $
    , BTICKS=BTICKS , BMINOR=BMINOR, D=D, SIZE=SIZE, XSIZE=XSIZE, YSIZE=YSIZE $
    , _EXTRA=EXTRA

a=ia
if (n_elements(D) ne 0) then begin
    bar=1
    aspect=1
    sample=1
endif
if (n_elements(BRANGE) ne 0) then begin
    a=a<max(BRANGE)>min(BRANGE)
;    if (min(a) gt min(BRANGE)) then begin
;        a(0,0)=min(BRANGE)
;        print,'Warning: Pixel (0,0) set to ',a(0,0),' for scaling'
;    endif
;    if (max(a) lt max(BRANGE)) then begin
;        a(0,1)=max(BRANGE)
;        print,'Warning: Pixel (0,1) set to ',a(0,1),' for scaling'
;    endif
endif
on_error,2
sa=size(reform(a)) & if sa(0) ne 2 then goto,errout
mina=min(a,max=maxa)
;
; set keyword parameters to default values if not present
if n_elements(XRANGI) eq 0 then XRANGI=[0,sa(1)-1]
if n_elements(YRANGI) eq 0 then YRANGI=[0,sa(2)-1]
if n_elements(   TITLE) eq 0 then    TITLE=''
if n_elements(SUBTITLE) eq 0 then SUBTITLE=''
if n_elements(  XTITLE) eq 0 then   XTITLE=''
if n_elements(  YTITLE) eq 0 then   YTITLE=''
if n_elements(TICKLEN) eq 0 then TICKLEN=-.01
if n_elements(XTICKS) eq 0 then XTICKS=0
if n_elements(XMINOR) eq 0 then XMINOR=0
if n_elements(YTICKS) eq 0 then YTICKS=0
if n_elements(YMINOR) eq 0 then YMINOR=0
if n_elements(CHARSIZE) eq 0 then CHARSIZE=1.0
XSTYLE=1  &  YSTYLE=1
if n_elements(XSTYLI) eq 1 then XSTYLE=( 1 or XSTYLI ) and 29
if n_elements(YSTYLI) eq 1 then YSTYLE=( 1 or YSTYLI ) and 29
if n_elements(BTITLE) eq 0 then BTITLE=''
if n_elements(BRANGE) eq 0 then BRANGE=float([mina,maxa])
if n_elements(BTICKS) eq 0 then BTICKS=0
if n_elements(BMINOR) eq 0 then BMINOR=0
;
XRANGE=float(tvframe_minmax(XRANGI))
YRANGE=float(tvframe_minmax(YRANGI))
;
if keyword_set(center) then begin
    xunit=0.5*(XRANGE(1)-XRANGE(0))/float(sa(1)-1)
    yunit=0.5*(YRANGE(1)-YRANGE(0))/float(sa(2)-1)
    XRANGE(0)=XRANGE(0)-xunit  &  XRANGE(1)=XRANGE(1)+xunit
    YRANGE(0)=YRANGE(0)-yunit  &  YRANGE(1)=YRANGE(1)+yunit
endif else begin
    xunit=(XRANGE(1)-XRANGE(0))/float(sa(1)-1)
    yunit=(YRANGE(1)-YRANGE(0))/float(sa(2)-1)
    XRANGE(1)=XRANGE(1)+xunit
    YRANGE(1)=YRANGE(1)+yunit
endelse
;
plot,XRANGE,YRANGE,/nodata,xstyle=xstyle or 4,ystyle=ystyle or 4    $
     ,  TITLE=' ',XTITLE=XTITLE,YTITLE=YTITLE    $
     ,  SUBTITLE=SUBTITLE      $
     ,  TICKLEN=TICKLEN,CHARSIZE=CHARSIZE   $
     ,  color=!p.background, _EXTRA=EXTRA
;
px = !x.window * !d.x_vsize     ;Position of frame in device units
py = !y.window * !d.y_vsize
sx = px(1)-px(0)                ;Size of frame in device units
sy = py(1)-py(0)
if keyword_set(bar) then sx = sx/1.25
if keyword_set(aspect) then begin
     f = float(sa(1))/sa(2)*sy/sx
     if f ge 1. then sy=sy/f else sx=sx*f
     sx=fix(sx)
endif

if keyword_set(size) then begin
    sx=sx*size
    sy=sy*size
endif
if keyword_set(xsize) then sx=sx*xsize
if keyword_set(ysize) then sy=sy*ysize

POSITION = [px(0),py(0),px(0)+sx,py(0)+sy]
if keyword_set(bar) or keyword_set(bbar) then begin
   bx    = fix(px(0)+sx*1.04)
   by    = fix(py(0))
   bsx   = fix(sx*0.08)
   bsy   = fix(sy)
   barpos= [bx,by,bx+bsx,by+bsy]
endif
;
mcol=( !D.N_COLORS - 1) > 0
;
if (!d.flags and 1) ne 0 then begin

      ;  scalable pixels

     if keyword_set(sample) then b=a else b=resize(a,256,256)
     if keyword_set(noscale) then begin
         tv, 0>b<mcol ,px(0),py(0),xsize=sx,ysize=sy,/device
     endif else begin
         maxb=max(b)
         minb=min(b)
         if (n_elements(BRANGE) ne 0) then begin
             if (min(b) gt min(BRANGE)) then minb=min(BRANGE) 
             if (max(b) lt max(BRANGE)) then maxb=max(BRANGE) 
         endif
         bb=(b-minb)/(maxb-minb)*255
         tv, bb ,px(0),py(0),xsize=sx,ysize=sy,/device
     endelse
     if keyword_set(bar) or keyword_set(bbar) then begin
        barim=findgen(1,256)/255*(maxa-mina)+mina
        if keyword_set(noscale) then begin
            tv, 0>barim<mcol ,bx,by,xsize=bsx,ysize=bsy,/device
        endif else begin
            tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device
        endelse
     endif

endif else begin

      ;  not scalable pixels

     if sx*sy gt 10e6 then begin
        print,' do you really want to allocate ',sx*sy   $
             ,' words of memory ? [y/n]'
        answer='n'   &  read,answer
        if answer ne 'y' then return
     endif
     if keyword_set(sample) then b=resize(a,sx,sy,/sample)  $
        else b=resize(a,sx,sy)
     if keyword_set(noscale) then begin
        tv, 0>b<mcol ,px(0),py(0),/device
     endif else begin
         maxb=max(b)
         minb=min(b)
         if (n_elements(BRANGE) ne 0) then begin
             if (min(b) gt min(BRANGE)) then minb=min(BRANGE)
             if (max(b) lt max(BRANGE)) then maxb=max(BRANGE)
         endif
         bb=(b-minb)/(maxb-minb)*255
         tv, bb ,px(0),py(0),/device
     endelse
     if keyword_set(bar) or keyword_set(bbar) then begin
        barim=findgen(1,bsy)/(bsy-1)*(maxa-mina)+mina
        barim=rebin(barim,bsx,bsy,/sample)
        if keyword_set(noscale) then begin
           tv, 0>barim<mcol ,bx,by,/device
        endif else begin
           tvscl, barim ,bx,by,/device
        endelse
     endif

endelse
;
if keyword_set(bar) or keyword_set(bbar) then begin
   plot,[0,1],BRANGE,/nodata,/noerase,pos=barpos,/device,xsty=5,ysty=5, $
        _EXTRA=EXTRA
   plots,[bx+bsx,bx,bx,bx+bsx],[by,by,by+bsy,by+bsy],/device
   axis,yaxis=1,bx+bsx,/device,yrange=BRANGE,ystyle=1   $
     ,  ytitle=BTITLE $
     ,  TICKLEN=-.15,CHARSIZE=CHARSIZE,YTICKS=BTICKS,YMINOR=BMINOR
endif
;
plot,XRANGE,YRANGE,/nodata,/noerase,xstyle=XSTYLE,ystyle=YSTYLE    $
     ,  POSITION=POSITION  ,/device    $
     ,  TITLE=TITLE,XTITLE=XTITLE,YTITLE=YTITLE    $
     ,  SUBTITLE=SUBTITLE      $
     ,  TICKLEN=TICKLEN,CHARSIZE=CHARSIZE   $
     ,  XTICKS=XTICKS,XMINOR=XMINOR , YTICKS=YTICKS,YMINOR=YMINOR $
     ,  _EXTRA=EXTRA
;
return
;
errout: print,' TVFRAME : data must be 2-dimensional !'
return
end
