! This routine performs a cubic interpolation of the curve that
! passes through x,y to xx,yy
!
Subroutine Cubic_interp(nn, x, y, xx, yy)
  Implicit None
  Integer :: nn, i
  Real, dimension (4) :: x, y
  Real, dimension (nn) :: xx, yy
  Real :: den, a, b, c, d, x1, x2, x3, x4, y1, y2, y3, y4
!
  x1=x(1)
  x2=x(2)
  x3=x(3)
  x4=x(4)
  y1=y(1)
  y2=y(2)
  y3=y(3)
  y4=y(4)
  den = (x1 - x2)*(x1 - x3)*(x2 - x3)*(x1 - x4)*(x2 - x4)*(x3 - x4)
  a=(x1*x4*(x1-x4)*(y2-y3) + x3*x3*(x4*(y1-y2)+x1*(y2-y4))+x2* &
       (x4*x4*(y1-y3)+x1*x1*(y3-y4)+x3*x3*(y4-y1))+ &
       x3*(x4*x4*(y2-y1)+x1*x1*(y4-y2))+ x2*x2* &
       (x4*(y3-y1)+x3*(y1-y4)+x1*(y4-y3)))/den

  b=(-x1*x4*(x1*x1-x4*x4)*(y2-y3)+x3*(x4*x4*x4*(y1-y2)+ &
       x1*x1*x1*(y2-y4))+x2*x2*x2*(x4*(y1-y3)+x1*(y3-y4)+ &
       x3*(y4-y1))+x3*x3*x3*(x4*(y2-y1)+x1*(y4-y2))+ &
       x2*(x4*x4*x4*(y3-y1)+x3*x3*x3*(y1-y1)+x1*x1*x1*(y4-y3)))/den

  c=(x1*x1*(x1-x4)*x4*x4*(y2-y3)+x3*x3*x3*(x4*x4*(y1-y2)+x1*x1*(y2-y4))+ &
       x2*x2*(x4*x4*x4*(y1-y3)+x1*x1*x1*(y3-y4)+x3*x3*x3*(y4-y1))+ &
       x3*x3*(x4*x4*x4*(y2-y1)+x1*x1*x1*(y4-y2))+ &
       x2*x2*x2*(x4*x4*(y3-y1)+x3*x3*(y1-y4)+x1*x1*(y4-y3)))/den

  d=(x1*(x1-x3)*x3*(x1-x4)*(x3-x4)*x4*y2+x2*x2*x2*(x1*(x1-x4)*x4*y3+ &
       x3*x3*(x1*y4-x4*y1)+x3*(x4*x4*y1-x1*x1*y4))+ &
       x2*(x1*x1*(x1-x4)*x4*x4*y3+x3*x3*x3*(x1*x1*y4-x4*x4*y1)+ &
       x3*x3*(x4*x4*x4*y1-x1*x1*x1*y4))+ &
       x2*x2*(x1*x4*(x4*x4-x1*x1)*y3+x3*x3*x3*(x4*y1-x1*y4)+ &
       x3*(x1*x1*x1*y4-x4*x4*x4*y1)))/den

  yy = a*(xx**3) + b*(xx**2) + c*xx + d

  Return
End Subroutine Cubic_interp
