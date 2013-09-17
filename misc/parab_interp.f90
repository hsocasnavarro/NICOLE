Subroutine Parab_interp(n, x, y, nn, xx, yy)
! This routine interpolates a vector x,y to another vector xx,yy
! Note: It doesn't check that xx is within x
!
  Implicit None
  Integer :: n, nn
  Real, dimension(n) :: x, y
  Real, dimension(3) :: x3, y3
  Real, dimension(nn) :: xx, yy
  Integer, dimension(1) :: imin1
  Integer :: ii, i1, i3
  Real :: pend
!
  If (n .eq. 1) then ! 1-point (absurd)
     yy=y
  Else if (n .eq. 2) then ! 2-point (linear interpolation)
     pend=(y(2)-y(1))/(x(2)-x(1))
     yy(1:nn)=pend*(xx(1:nn)-x(1)) + y(1)
  Else if (n .eq. 3) then ! 3-point (parabolic interpolation)
     Call ThreePt_Parab_interp(nn, x, y, xx, yy)
  Else ! More than 3 points (parabolic)
     Do ii=1, nn
        imin1=MinLoc(Abs( xx(ii)-x(:) ))
        i1=imin1(1)-1
        i3=imin1(1)+1
        If (i1 .le. 0) then
           i1=1
           i3=3
        Else if (i3 .gt. n) then
           i1=n-2
           i3=n
        End if
        x3(1:3)=x(i1:i3)
        y3(1:3)=y(i1:i3)
        Call ThreePt_parab_interp(1, x3, y3, xx(ii), yy(ii))
     End do
  End if
  Return
End Subroutine Parab_interp
! This routine performs a parabolic interpolation of the parabola that
! passes through 3 points x,y to a vector xx,yy
!
Subroutine ThreePt_Parab_interp(nn, x, y, xx, yy)
  Implicit None
  Integer :: nn
  Real, dimension (3) :: x, y
  Real, dimension (nn) :: xx, yy
  Real :: den, a, b, c, x1, x2, x3, y1, y2, y3
!
  x1=x(1)
  x2=x(2)
  x3=x(3)
  y1=y(1)
  y2=y(2)
  y3=y(3)
  den = (x1 - x2)*(x1 - x3)*(x2 - x3)
  a = (x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))/den
  b = (x3*x3*(y1-y2)+x1*x1*(y2-y3)+x2*x2*(y3-y1))/den
  c = (x1*x3*y2*(x3-x1)+x2*x2*(x3*y1-x1*y3)+ &
       x2*(x1*x1*y3-x3*x3*y1))/den
  yy = a*xx*xx + b*xx + c
  Return
End Subroutine ThreePt_Parab_interp
