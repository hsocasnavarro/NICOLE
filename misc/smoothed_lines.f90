! This routine performs an interpolation of the curve that
! passes through x,y to xx,yy. The interpolation is performed joining
! adjacent points by a straight line and then smoothing the result 
! by convolving it with a Gaussian whose width is a half of the distance
! between the two closest nodes.
!
Subroutine Smoothed_lines(n, nn, x, y, xx, yy)
  Implicit None
  Integer :: n, nn, inode, ipoint, jpoint, i
  Real, dimension (n) :: x, y
  Real, dimension (nn) :: xx, yy, yy2
  Real :: x1, x2, y1, y2, closest, expo, Norm
!
! Compute segments
!
  yy=0.
  Do inode=1, n-1
     x1=x(inode)
     x2=x(inode+1)
     y1=y(inode)
     y2=y(inode+1)
     Do ipoint=1, nn
        If (xx(ipoint) .ge. x1 .and. xx(ipoint) .le. x2) then
           yy(ipoint)=y1+(y2-y1)/(x2-x1)*(xx(ipoint)-x1)
        End If
     End do
  End do



! Do not smooth ! DEBUG
! For some reason the inversions don't work very well with this option
  Return




!
! Smooth the result
!     
!     Find the closest nodes
!
  closest=abs(x(1)-x(2))
  Do inode=2, n-1
     If (abs(x(inode)-x(inode+1)) .lt. closest) &
          closest=abs(x(inode)-x(inode+1))
  End do
!
! Convolve (note that the model is not necessarily equispaced)
!
  closest=closest/2.
  yy2=0.
  Do ipoint=1, nn
     Norm=0.
     Do jpoint=1,nn
        x1=xx(jpoint)-xx(ipoint)
        If (x1*x1/closest/closest .lt. 25.) then
           expo=exp(-x1*x1/(closest*closest))
           yy2(ipoint)=yy2(ipoint)+yy(jpoint)*expo
           Norm=Norm+expo
        End If
     End Do
     yy2(ipoint)=yy2(ipoint)/Norm
  End Do
  yy=yy2
!
! Done
!
  Return
End Subroutine Smoothed_lines
