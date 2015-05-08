MODULE bezier_math
  !
  ! Bezier-Splines utilities. These tools provide accurate Bezier Interpolants
  ! that minimize overshooting of the interpolated function. Additional tools
  ! to integrate a functions analytically, assuming a Bezier interpolant are
  ! also included (these are used to obtain tau_nu)
  ! 
  ! CONTAINS:
  ! cent_deriv     -> This subroutine computes centered derivatives
  ! bezier         -> Performs interpolation assuming a cuadratic bezier interpolant
  ! bezier3        -> Performs interpolation assuming a cubic Bezier interpolant
  ! bezier_cqintep -> Cumulative integral of a function, starting at element 1. Quadratic Bezier
  ! bezier_qintep  -> Total integral of a function. Quadratic Bezier.
  ! bezier_ccintep -> Cumulative integral of a function, starting at element 1. Cubic Bezier
  ! bezier_cintep  -> Total integral of a function. Cubic Bezier.
  !
  ! USAGE:
  ! call     cent_deriv(npoints_in, x_in, y_in, yprime)
  ! call         bezier(npoints_in, x_in, y_in, npoints_new, x_new, y_out)
  ! call        bezier3(npoints_in, x_in, y_in, npoints_new, x_new, y_out)
  ! call  bezier_qintep(npoints_in, x_in, y_in, result)
  ! call  bezier_cintep(npoints_in, x_in, y_in, result)
  ! call bezier_ccintep(npoints_in, x_in, y_in, cresult)
  ! call bezier_cqintep(npoints_in, x_in, y_in, cresult)
  !
  ! NOTE:
  !  The formal solvers alter the input arrays tau and abs matrix
  ! WHERE:
  ! npoints_in (int) is the number of elements of x_in, y_in, yprime, cresult (floats).
  ! npoints_new (int) is the number of elements of x_new, y_out (floats).
  ! result is a float. Notice that cresult is an array and result is not. 
  !
  ! REFERENCES: 
  ! Auer 2003 (Interpolants).
  ! Fritsch and Butland 1984 (Derivatives).
  !
  ! The analytical integral of the Bezier interpolants used in the integration and
  ! cumulative integrated was computed by the author.
  ! 
  ! Jaime de la Cruz Rodriguez (IFA-Uppsala University 2012)
  !
  ! MODIFICATIONS:
  ! 2012-10-04, JdlCR: Created!
  !
  ! -
  Implicit None
CONTAINS
  subroutine ctrapez(n,x,y,inte)
    Implicit none
    Integer :: n, np, k, k1
    real :: x(n), y(n), inte(n)
    real(8) :: tmp
    
    tmp = inte(1)
    do k = 2, n
       tmp = tmp + 0.5d0*(y(k) + y(k-1)) * (x(k) - x(k-1))
       inte(k) = tmp
    end do


    return
  end subroutine ctrapez
  subroutine trapez(n,x,y,inte)
    Implicit none
    Integer :: n, np, k, k1
    real :: x(n), y(n), inte(n)
    
    do k = 2, n
       inte(k) = 0.5d0*(y(k) + y(k-1)) * (x(k) - x(k-1))
    end do


    return
  end subroutine trapez
  subroutine cent_deriv(n,x,y,yp)
    Implicit None
    integer :: n, k
    Real :: x(n), y(n)
    real(8) :: der, der1, lambda, yp(n), dx , dx1
    
    do k = 2, n - 1
       dx = x(k) - x(k-1)
       der = (y(k) - y(k-1)) / dx
       
       dx1 = x(k+1) - x(k)
       der1 = (y(k+1) - y(k)) / dx1
       
       if(der*der1 .gt. 0.0d0) then 
          lambda = (1.d0 + dx1 / (dx1 + dx)) / 3.d0
          yp(k) = (der / (lambda * der1 + (1.d0 - lambda) * der)) * der1
       Else
          yp(k) = 0.0d0
       end if
    enddo
    
    yp(1) =  (y(1) - y(2)) / (x(1) - x(2))
    yp(n) = (y(n-1) - y(n)) / (x(n-1) - x(n))
    
    return
  end subroutine cent_deriv
  subroutine cent_deriv_d(n,x,y,yp)
    Implicit None
    integer :: n, k
    Real :: x(n), y(n)
    real(8) :: der, der1, lambda, yp(n), dx , dx1
    
    do k = 2, n - 1
       dx = x(k) ! x(k) - x(k-1)
       der = (y(k) - y(k-1)) / dx
       
       dx1 = x(k+1) ! x(k+1) - x(k)
       der1 = (y(k+1) - y(k)) / dx1
              
       if(der*der1 .gt. 0.0d0) then 
          lambda = (1.d0 + dx1 / (dx1 + dx)) / 3.d0
          yp(k) = (der / (lambda * der1 + (1.d0 - lambda) * der)) * der1
       Else
          yp(k) = 0.0d0
       end if
    enddo
    
    yp(1) =  (y(2) - y(1)) / (x(2)) ! (x(1) - x(2))
    yp(n) = (y(n) - y(n-1)) / (x(n)) ! (x(n-1) - x(n))
    
  

    return
  end subroutine cent_deriv_d
  
  subroutine bezier(n, x, y, np, xp, yp)
    Implicit None
    Integer :: n, np, k
    Real, dimension(n) :: x, y
    Real, dimension(np) :: xp, yp
    Real(8) :: cntrl, dx, yprime(n), lambda, u(np)
    
    
    !
    ! Compute derivatives 
    !
    cntrl = 0
    yprime = 0
    call cent_deriv(n, x, y, yprime)
    
    do k = 2, n
       dx =  x(k) - x(k-1)
       
       
       cntrl = 0.5d0 * (y(k) - 0.5d0*dx * yprime(k) + y(k-1) + 0.5d0*dx * yprime(k-1))
       if(cntrl .GT. max(y(k),y(k-1)) .OR. cntrl .LT. min(y(k),y(k-1))) cntrl = y(k-1) 

       where(xp .LT. x(k) .AND. xp .GE. x(k-1))
          u = (xp - x(k-1)) / dx
          yp = y(k-1) * (1.d0 - u)**2 + y(k) * u**2 + 2.d0 * cntrl * u * (1.d0 - u)
       End where
    end do
    
    !
    ! Points outside the domain
    !
    where(xp .GE. x(n))
       yp = y(n)
    end where
    where(xp .LE. x(1))
       yp = y(1)
    end where
  end subroutine bezier
  subroutine bezier3(n, x, y, np, xp, yp)
    Implicit None
    Integer :: n, np, k
    Real, dimension(n) :: x, y
    Real, dimension(np) :: xp, yp
    Real(8) :: c1, c2, yprime(n), dx, u(np), mmi, mma
    
    
    c1 = 0
    c2 = 0
    yprime = 0
    !
    ! Compute derivatives 
    !
    call cent_deriv(n, x, y, yprime)
    !
    do k = 2, n
       dx =  x(k) - x(k-1)
       
       c1 = y(k-1) + dx * yprime(k-1) / 3.0d0
       c2 = y(k) - dx * yprime(k) / 3.0d0

       mmi = min(y(k),y(k-1))
       mma = max(y(k),y(k-1))
       if(c1 .LT. mmi .OR. c1 .GT. mma) c1 = y(k-1)
       if(c2 .LT. mmi .OR. c2 .GT. mma) c2 = y(k)

       where(xp .LT. x(k) .AND. xp .GE. x(k-1))
          u = (xp - x(k-1)) / dx
          yp = y(k-1) * (1.d0 - u)**3 + y(k) * u**3 + &
               3.d0 * c1 * u * (1.d0 - u)**2 + 3.0d0 * c2 * u**2 * (1.0d0 - u)
       End where
       
    end do
    
    !
    ! Points outside the domain
    !
    where(xp .GE. x(n))
       yp = y(n)
    end where
    where(xp .LE. x(1))
       yp = y(1)
    end where
  end subroutine bezier3
  
  subroutine bezier_cqintep(n, x, y, inte)
    Implicit none
    Integer :: n, np, k, k1
    real :: x(n), y(n), inte(n)
    real(8) :: dx, dx2, dI, I, yprime(n), c, mma, mmi
    
    call cent_deriv(n, x, y, yprime)
    
    Do k = 1, n - 1
       k1 = k+1
       dx = x(k+1) - x(k)
       
       dx2 = dx * 0.5d0

       c = (y(k) + dx2 * yprime(k) + y(k1) - dx2 * yprime(k1)) * 0.5d0
     !  mmi = min(y(k), y(k1))
     !  mma = max(y(k), y(k1))
     !  if(c .gt. mma)  c = mma
     !  if(c .lt. mmi)  c = mmi

       inte(k1) = inte(k) +  dx / 3.0d0 * (y(k) + y(k1) + c)
    end do
    inte(n)=y(n)*(x(n)-x(n-1))
    
    return
  End subroutine bezier_cqintep
  
  subroutine bezier_qintep(n, x, y, inte)
    Implicit none
    Integer :: n, np, k, k1
    real :: x(n), y(n), inte(n)
    real(8) :: dx, dx2, dI, I, yprime(n), mmi, mma, c, c1
    
    call cent_deriv(n, x, y, yprime)
    
    inte(1)=y(1)*(x(2)-x(1))
    do k = 1, n - 1
       k1 = k+1
       dx = x(k1) - x(k)
       
       dx2 = dx * 0.5d0
       
       c = (y(k) + dx2 * yprime(k) + y(k+1) - dx2 * yprime(k+1)) * 0.5d0
   !    mmi = min(y(k), y(k1))
   !    mma = max(y(k), y(k1))
   !    if(c .gt. mma)  c = mma
   !    if(c .lt. mmi)  c = mmi

       inte(k1) = dx * (y(k) + y(k+1) + c) / 3.0d0
    end do
    inte(n)=y(n)*(x(n)-x(n-1))
    
    return
  End subroutine bezier_qintep
  
  subroutine bezier_ccintep(n, x, y, inte)
    Implicit none
    Integer :: n, np, k, k1
    real :: x(n), y(n), inte(n)
    real(8) :: dx, dx3, dI, I, yprime(n), c, c1, mmi, mma
    
    call cent_deriv(n, x, y, yprime)
    
    Do k = 1, n - 1
       k1 = k + 1
       dx = x(k+1) - x(k)
       
       dx3 = dx / 3.0d0
       
       c = y(k) + dx3 * yprime(k)
       c1 =  y(k1) - dx3 * yprime(k1)
       
  !     mmi = min(y(k),y(k1))
  !     mma = max(y(k),y(k1))
  !     if(c  .LT. mmi .OR.  c .GT. mma) c = y(k)
  !     if(c1 .LT. mmi .OR. c1 .GT. mma) c1 = y(k1)

       inte(k1) = inte(k) + dx * (y(k) + y(k1) + c + c1) * 0.25d0
    end do
    
    return
  End subroutine bezier_ccintep
  
  subroutine bezier_cintep(n, x, y, inte)
    Implicit none
    Integer :: n, np, k, k1
    real :: x(n), y(n), inte(n), c, c1
    real(8) :: dx, dx3, dI, I, yprime(n), mmi, mma
    
    call cent_deriv(n, x, y, yprime)
    
    
    do k = 1, n - 1
       k1 = k+1
       dx = x(k+1) - x(k)
       dx3 = dx / 3.0d0
       
       c = y(k) + dx3 * yprime(k)
       c1 =  y(k1) - dx3 * yprime(k1)
       
       
      ! mmi = min(y(k),y(k1))
      ! mma = max(y(k),y(k1))
      ! if(c  .LT. mmi .OR.  c .GT. mma) c = y(k)
      ! if(c1 .LT. mmi .OR. c1 .GT. mma) c1 = y(k1)

       inte(k1) = dx * (y(k) + y(k1) + c + c1) * 0.25d0
    end do
    
    return
  End subroutine bezier_cintep
END MODULE bezier_math

MODULE Bezier_solvers
  Implicit None
  !
  ! Bezier Solvers. Integrates the RTE for polarized light, assuming 
  ! Bezier-Splines interpolants. 
  !
  ! CONTAINS: full_cent_der    -> Aux. routine that computes derivatives of a 4xndep array.
  !           full_cent_der4   -> Aux. routine that computes derivatives of a 4x4xndep array.
  !           matinx8          -> Inverts a matrix
  !           delobezier       -> Integrates the RTE assuming a quadratic Bezier-Interpolant
  !           delobezier3      -> Integrates the RTE assuming a cubic Bezier-Interpolant
  ! 
  ! USAGE: the interface of this routines is identical to the other solvers in NICOLE. Same 
  !        input/output variables are expected.
  !
  ! REFERENCES: de la Cruz Rodriguez & Piskunov in prep.
  !             
  !
  ! Jaime de la Cruz Rodriguez (IFA-Uppsala University 2012)
  !
  ! Modifications:
  !            2012-10-04, JdlCR: Created!
  !
CONTAINS
  subroutine full_cent_der_d(n, x, y, yp, ubound, ibound)
    Implicit None
    Integer :: n, k, k1, ibound, ubound, kinit, kend
    Real :: x(n), yp(4,n), y(4, n)
    real(8) :: lambda, dx, dx1, der(4), der1(4)
    
    kinit = max(2, ubound)
    kend = min(n-1, ibound)

    yp = 0
    do k = kinit, kend
       dx = x(k)
       dx1 = x(k+1)
       
       der  = (y(:,k) - y(:,k-1)) / dx
       der1 = (y(:,k+1) - y(:,k)) / dx1
       
       lambda =  (1.0d0 + dx1 / (dx + dx1)) / 3.0d0
       where(der * der1 > 0.0d0)
          yp(:,k) = (der / (lambda * der1 + (1.0d0 - lambda) * der)) * der1
       end where
    end do
    
    if(ubound == 1) yp(:,1) = (y(:,1) - y(:,2)) / (-x(2)) ! (x(1) - x(2))
    if(ibound == n) yp(:,n) = (y(:,n-1) - y(:,n)) / (-x(n)) ! (x(n-1) - x(n))
    
  end subroutine full_cent_der_d
  subroutine full_cent_der(n, x, y, yp, ubound, ibound)
    Implicit None
    Integer :: n, k, k1, ibound, ubound, kinit, kend
    Real :: x(n), yp(4,n), y(4, n)
    real(8) :: lambda, dx, dx1, der(4), der1(4)
    
    kinit = max(2, ubound)
    kend = min(n-1, ibound)

    yp = 0
    do k = kinit, kend
       dx = x(k) - x(k-1)
       dx1 = x(k+1) - x(k)
       
       der  = (y(:,k) - y(:,k-1)) / dx
       der1 = (y(:,k+1) - y(:,k)) / dx1
       
       lambda =  (1.0d0 + dx1 / (dx + dx1)) / 3.0d0
       where(der * der1 > 0.0d0)
          yp(:,k) = (der / (lambda * der1 + (1.0d0 - lambda) * der)) * der1
       end where
    end do
    
    if(ubound == 1) yp(:,1) = (y(:,1) - y(:,2)) / (x(1) - x(2))
    if(ibound == n) yp(:,n) = (y(:,n-1) - y(:,n)) / (x(n-1) - x(n))
    
  end subroutine full_cent_der
  subroutine full_cent_der4(n, x, y, yp, ubound, ibound)
    Implicit None
    Integer :: n, k, k1, ibound, ubound, kinit, kend, i, j
    Real ::  y(n,4,4), yp(n,4,4), x(n)
    real(8) ::  lambda(n), dx, dx1, der, der1
    
    kinit = max(2, ubound)
    kend = min(n-1, ibound)

    yp = 0
    
    dx = x(kinit) - x(kinit-1)
    do k = kinit, kend
       dx1 = x(k+1) - x(k)
       lambda(k) =  (1.0d0 + dx1 / (dx + dx1)) / 3.0d0
       dx = dx1
    end do


    do j = 1,4
       do i = 1,4
          do k = kinit, kend
             
             der  = (y(k,i,j) - y(k-1,i,j)) / (x(k) - x(k-1))
             der1 = (y(k+1,i,j) - y(k,i,j)) / (x(k+1) - x(k))
             
             if(der * der1 .gt. 0.0d0) then
                yp(k,i,j) = (der / (lambda(k) * der1 + (1.0d0 - lambda(k)) * der)) * der1
             end if
          end do
       end do
    end do
    
    if(ubound == 1) yp(1,:,:) = (y(1,:,:) - y(2,:,:)) / (x(1) - x(2))
    if(ibound == n) yp(n,:,:) = (y(n-1,:,:) - y(n,:,:)) / (x(n-1) - x(n))
    
  end subroutine full_cent_der4
  subroutine full_cent_der4_d(n, x, y, yp, ubound, ibound)
    Implicit None
    Integer :: n, k, k1, ibound, ubound, kinit, kend, i, j
    Real ::  y(n,4,4), yp(n,4,4), x(n)
    real(8) ::  lambda(n), dx, dx1, der, der1
    
    kinit = max(2, ubound)
    kend = min(n-1, ibound)

    yp = 0
    
    dx = x(kinit) ! x(kinit) - x(kinit-1)
    do k = kinit, kend
       dx1 = x(k+1) ! x(k+1) - x(k)
       lambda(k) =  (1.0d0 + dx1 / (dx + dx1)) / 3.0d0
       dx = dx1
    end do


    do j = 1,4
       do i = 1,4
          do k = kinit, kend
             
             der  = (y(k,i,j) - y(k-1,i,j)) / x(k) ! (x(k) - x(k-1))
             der1 = (y(k+1,i,j) - y(k,i,j)) / x(k+1) ! (x(k+1) - x(k))
             
             if(der * der1 .gt. 0.0d0) then
                yp(k,i,j) = (der / (lambda(k) * der1 + (1.0d0 - lambda(k)) * der)) * der1
             end if
          end do
       end do
    end do
    
    if(ubound == 1) yp(1,:,:) = (y(1,:,:) - y(2,:,:)) / (-x(2)) ! (x(1) - x(2))
    if(ibound == n) yp(n,:,:) = (y(n-1,:,:) - y(n,:,:)) / (-x(n)) ! (x(n-1) - x(n))
    
  end subroutine full_cent_der4_d
  
  !
  ! Solvers
  !
  Subroutine delobezier3(npoints, dtau_nu, Ab, S, stokes)
    USE Bezier_math
    Implicit None
    Integer :: ipoint, npoints, ibound, ii, jj, k, k1, info, ubound
    Integer :: ipiv(4), myunit
    Real, dimension(npoints) :: tau_500, tau_nu, iab, S, dtau_nu
    Real, dimension(npoints, 4, 4) :: Ab, abp
    Real, dimension(4,npoints) :: sv, sp
    Real, dimension(4) :: stokes
    Real(8), dimension(4,4) :: ident,  A, tmpa, tmpb
    Real(8) :: dtau, dtau2, idtau3, dtau03, dtau3, dtau4
    Real(8) :: alpha, beta, gamma, mu, eps, ostokes(4)
    Real, Parameter :: Optically_thick = 100., Optically_thin=1.e-4

    !
    ! 4x4 Identity matrix
    !
    Data ident(:,1)/1.0d0,0.0d0,0.0d0,0.0d0/
    Data ident(:,2)/0.0d0,1.0d0,0.0d0,0.0d0/
    Data ident(:,3)/0.0d0,0.0d0,1.0d0,0.0d0/
    Data ident(:,4)/0.0d0,0.0d0,0.0d0,1.0d0/

    tau_nu(1)=0.
    do k=2, npoints
       tau_nu(k)=tau_nu(k-1)+dtau_nu(k)
    end do

    iab = ab(:,1,1)

    !
    ! Find the point for the boundary condition
    !
    ubound = 1
    ibound = npoints
    do k = 2, npoints-1
       If (Tau_Nu(k) .gt. (Optically_Thin) .and. Tau_Nu(k-1) .lt. (Optically_Thin)) &
            ubound=k-1
       If (Tau_Nu(k) .gt. (Optically_Thick) .and. Tau_Nu(k-1) .lt. (Optically_Thick)) &
            ibound=k
    End do
    if(ibound .lt. npoints) ibound = ibound + 1
    if(ubound .gt. 1) ubound = ubound - 1


    !
    ! Modified source vector
    !
    do k = max(1, ubound-1), min(npoints, ibound+1)
       do ii = 1,4
          Sv(ii,k) = ab(k,ii,1) * S(k) / iab(k)
       end do
    enddo

    !
    ! Convert K to K' (subtract diagonal and normalize to K(1,1))
    !
    do jj = 1, 4
       do ii = 1, 4
          do k = max(1, ubound-1), min(npoints, ibound+1)
             ab(k,ii,jj) = ab(k,ii,jj) / iab(k)
          enddo
       enddo
       ab(:,jj,jj) = 0.0d0 !ab(ii,ii,:) - 1.d0      
    enddo


    !
    ! Set boundary condition (I=S) at ibound
    !
    ostokes(2:4)=0.0d0
    ostokes(1) = S(ibound) ! Boundary cond.

    !
    ! Derivatives for Ab and Sv
    !
    call full_cent_der_d(npoints, dtau_nu, Sv, sp, ubound, ibound)
    call full_cent_der4_d(npoints, tau_nu, ab, abp, ubound, ibound)

    ibound = ibound - 1

    if(ubound .GE. ibound) then 
       print *, 'Error, ibound LT ubound'
       stop
    endif

    do k = ibound, ubound, -1
       k1 = k + 1

       !
       ! Compute dtau(k), ...
       ! 
       dtau = dtau_nu(k1) ! tau_nu(k1) - tau_nu(k)
       dtau03 = dtau / 3.0d0
       dtau2 = dtau * dtau
       dtau3 = dtau2 * dtau
       dtau4 = dtau2 * dtau2

       
       !
       ! Integration coeffs.
       !
       if(dtau >= 1.d-2) then ! Accurate integration
          eps = exp(-dtau)

          idtau3 = 1.0d0 / dtau3

          alpha = (-6.0d0 + 6.0d0 * dtau - 3.0d0*dtau2 + dtau3 + 6.0d0 * eps) * idtau3
          beta  = (6.d0 + (-6.0d0 - dtau*(6.0d0 + dtau*(3.d0 + dtau))) * eps) * idtau3
          gamma = 3.0d0 * (6.0d0 + (-4.0d0 + dtau)*dtau - 2.0d0 * (3.0d0+dtau)*eps) * idtau3
          mu = 3.0d0 * ( eps * (6.0d0 + dtau2 + 4.d0 * dtau) + 2.d0 * dtau - 6.d0) * idtau3
       else ! Expansion
          eps = 1.d0 - dtau + 0.5d0 * dtau2 - dtau3/6.d0 + dtau4/24.d0

          alpha = 0.25d0 * dtau - 0.05d0 * dtau2 + dtau3 / 120.d0  - dtau4 / 840.d0 
          beta  = 0.25d0 * dtau - 0.20d0 * dtau2 + dtau3 / 12.0d0  - dtau4 / 42.0d0 
          gamma = 0.25d0 * dtau - 0.10d0 * dtau2 + dtau3 * 0.025d0 - dtau4 / 210.d0 
          mu    = 0.25d0 * dtau - 0.15d0 * dtau2 + dtau3 * 0.05d0  - dtau4 / 84.d0 
       end if


       !
       ! Interval integration
       !
       tmpa =    (dtau03 * (matmul(Ab(k1,:,:), Ab(k1,:,:)) + abp(k1,:,:)  + Ab(k1,:,:))  - Ab(k1,:,:))
       tmpb =  - (dtau03 * (matmul(Ab(k,:,:),  Ab(k,:,:))  + abp(k,:,:)   + Ab(k,:,:))   + Ab(k,:,:))

       A =                  ident + alpha * Ab(k,:,:)  - gamma * tmpb ! I_b
       ostokes = matmul((eps*ident -  beta * Ab(k1,:,:) + mu    * tmpa), ostokes) + &! I_a
            matmul((beta *ident +    mu * (ident - dtau03*ab(k1,:,:))), Sv(:,k1)) + &! S_a
            matmul((alpha*ident + gamma * (ident + dtau03*Ab(k,:,:))) , Sv(:,k)) + &! S_b
            dtau03 * (gamma*sp(:,k) - mu*Sp(:,k1))! Sp_a and Sp_b


       ! Solve system AX = B with LaPack's routine. This
       ! seems to be numerically more stable than inverting the Matrix
       ! when it has singular values.  "stokes" is an input/output 
       ! array.
       !
      ! call dgesv(int(4), int(1), A, int(4), ipiv, ostokes, int(4), info)

       !
       ! Alternatively use matrix inversion
       !
       call matinx8(a)
       ostokes = matmul(a, ostokes)

    End Do

    stokes = ostokes

  End Subroutine delobezier3

  Subroutine delobezier(npoints, dtau_nu, Ab, S, stokes)
    USE Bezier_math
    Implicit None
    Integer :: ipoint, npoints, ibound, ii, jj, k, k1, info, ubound
    Integer :: ipiv(4)
    Real, dimension(npoints) :: dtau_nu, S, iab, tau_nu, tau_500
    Real, dimension(npoints,4, 4) :: Ab, abp
    Real, dimension(4,npoints) ::  sv, opac, sp
    Real, dimension(4) :: stokes
    Real(8), dimension(4,4) :: ident, A, tmpa, tmpb
    Real(8) :: dtau, dtau1, dtau2, idtau2, eps, dtau05
    Real(8) :: alpha, beta, gamma, lambda, ostokes(4)
    Real, Parameter :: Optically_thick = 100., Optically_thin=1.e-4

    !
    ! 4x4 Identity matrix
    !
    Data ident(:,1)/1.0d0,0.0d0,0.0d0,0.0d0/
    Data ident(:,2)/0.0d0,1.0d0,0.0d0,0.0d0/
    Data ident(:,3)/0.0d0,0.0d0,1.0d0,0.0d0/
    Data ident(:,4)/0.0d0,0.0d0,0.0d0,1.0d0/

    iab = ab(:,1,1)

    !
    ! Find the point for the boundary condition
    !
    ubound = 1
    ibound = npoints

    tau_nu(1)=0.
    do k=2, npoints
       tau_nu(k)=tau_nu(k-1)+dtau_nu(k)
    end do

    do k = 2, npoints-1
       If (Tau_Nu(k) .gt. (Optically_Thin) .and. Tau_Nu(k-1) .lt. (Optically_Thin)) &
            ubound=k-1
       If (Tau_Nu(k) .gt. (Optically_Thick) .and. Tau_Nu(k-1) .lt. (Optically_Thick)) &
            ibound=k
    End do
    if(ibound .lt. npoints) ibound = ibound + 1
    if(ubound .gt. 1) ubound = ubound - 1


    !
    ! Modified source vector
    !
    do k = max(1, ubound-1), min(npoints, ibound+1)
       do ii = 1,4
          Sv(ii,k) = ab(k,ii,1) * S(k) / iab(k)
       end do
    enddo

    !
    ! Convert K to K' (subtract diagonal and normalize to K(1,1))
    !
    do jj = 1, 4
       do ii = 1, 4
          do k = max(1, ubound-1), min(npoints, ibound+1)
             ab(k,ii,jj) = ab(k,ii,jj) / iab(k)
          enddo
       enddo
       ab(:,jj,jj) = 0.0d0 !ab(ii,ii,:) - 1.d0
    enddo


    !
    ! Set boundary condition (I=S) at ibound
    !
    ostokes(2:4)=0.0d0
    ostokes(1) = S(ibound) ! Boundary cond.

    !
    ! Derivatives of the source function and ab matrix
    !
    call full_cent_der_d(npoints, dtau_nu, Sv, sp, ubound, ibound)
    call full_cent_der4_d(npoints, dtau_nu, ab, abp, ubound, ibound)

    ibound = ibound - 1

    do k = ibound, ubound, -1
       !
       k1 = k + 1

       !
       ! Compute dtau(k), dtau(k-1) -> dtau, dtau1 in Nik's solver
       !
       dtau = dtau_nu(k1) ! tau_nu(k1) - tau_nu(k)
       dtau05 = dtau * 0.5d0
       dtau2 = dtau * dtau
       idtau2 = 1.0d0 / dtau2



       !
       ! Terms of the RT-Bezier integral
       !
       if(dtau < 1.e-2) then ! Analytical expansion
          eps = 1.d0 - dtau + 0.5d0 * dtau2 - dtau2*dtau/6.d0 + dtau2*dtau2/24.d0

          alpha = dtau/3.d0 - dtau2/12.d0   + dtau2*dtau/60.d0
          beta  = dtau/3.d0 - dtau2*0.25d0  + dtau2*dtau*0.1d0
          gamma = dtau/3.d0 - dtau2/6.0d0   + dtau2*dtau*0.05d0
       else                  ! Accurate calculation
          eps = exp(-dtau)
          alpha = (dtau2 - 2.0d0*dtau + 2.0d0 - 2.0d0*eps) * idtau2
          beta  = (2.0d0 - (2.0d0 + 2.0d0*dtau + dtau2)*eps) * idtau2
          gamma = (2.d0*dtau - 4.0d0 + (2.d0*dtau+4.0d0)*eps) * idtau2
       endif

       gamma = gamma * 0.5d0

       
       !
       ! RT Integration
       !
       tmpa =    (dtau05 * (matmul(Ab(k1,:,:), Ab(k1,:,:)) + abp(k1,:,:)  + Ab(k1,:,:))  - Ab(k1,:,:))
       tmpb =  - (dtau05 * (matmul(Ab(k,:,:),  Ab(k,:,:))  + abp(k,:,:)   + Ab(k,:,:))   + Ab(k,:,:))

       A =                   ident + alpha * Ab(k,:,:)  - gamma * tmpb ! I_b terms
       ostokes = matmul((eps*ident -  beta * Ab(k1,:,:) + gamma * tmpa), ostokes) + &! I_a terms
            matmul((beta *ident + gamma * (ident - dtau05*ab(k1,:,:))) , Sv(:,k1)) + &! S_a terms
            matmul((alpha*ident + gamma * (ident + dtau05*Ab(k,:,:)))  , Sv(:,k))  &! S_b terms
            +dtau05 * gamma * (Sp(:,k) - Sp(:,k1))! Sp_a and Sp_b terms


       !
       ! Solve system AX = B with LaPack's routine. This
       ! seems to be numerically more stable than inverting the Matrix
       ! when it has singular values.  "stokes" is an input/output 
       ! array.
       !
      ! call dgesv(int(4), int(1), A, int(4), ipiv, ostokes, int(4), info)

       !
       ! Alternatively use a explicit matrix inversion (uncomment lines)
       !
       call matinx8(A)
       ostokes = matmul(A,ostokes)

    End Do
    stokes = ostokes

  END Subroutine delobezier
  
  SUBROUTINE delobezier_scal(npoints, dtau_nu, Ab, S, stokes)
    Use bezier_math
    Implicit None
    Integer :: npoints, k, k1, ibound, ubound
    Real, dimension(npoints) :: dtau_nu, tau_nu, tau_500, S, iab
    Real, dimension(npoints,4,4) :: Ab
    Real, dimension(4) :: stokes
    Real(8) :: dtau, dtau1, dtau2, idtau2, eps, dtau05, I
    Real(8) :: alpha, beta, gamma, cntrl, sp(npoints)
    Real, Parameter :: Optically_thick = 100., Optically_thin=1.e-3


    iab = ab(:,1,1)

    tau_nu(1)=0.
    do k=2, npoints
       tau_nu(k)=tau_nu(k-1)+dtau_nu(k)
    end do

    !
    ! Find the point for the boundary condition
    !
    ubound = 1
    ibound = npoints

    do k = 2, npoints-1
       If (Tau_Nu(k) .gt. (Optically_Thin) .and. Tau_Nu(k-1) .lt. (Optically_Thin)) &
            ubound=k-1
       If (Tau_Nu(k) .gt. (Optically_Thick) .and. Tau_Nu(k-1) .lt. (Optically_Thick)) &
            ibound=k
    End do
    if(ibound .lt. npoints) ibound = ibound + 1
    if(ubound .gt. 1) ubound = ubound - 1

    !
    ! Source function derivatives
    !
    call cent_deriv(npoints, tau_nu, S, Sp)
    
    !
    ! Boundary condition
    !
    I = S(ibound)

    ibound = ibound - 1
    
    do k = ibound, ubound, -1
       !
       k1 = k + 1

       !
       ! Compute dtau(k), dtau(k-1) -> dtau, dtau1 in Nik's solver
       !
       dtau = dtau_nu(k1) ! tau_nu(k1) - tau_nu(k)
       dtau05 = dtau * 0.5d0
       dtau2 = dtau * dtau
       idtau2 = 1.0d0 / dtau2

       !
       ! Terms of the RT-Bezier integral
       !
       if(dtau < 1.e-2) then ! Analytical expansion
          eps = 1.d0 - dtau + 0.5d0 * dtau2 - dtau2*dtau/6.d0 + dtau2*dtau2/24.d0
          alpha = dtau/3.d0 - dtau2/12.d0   + dtau2*dtau/60.d0
          beta  = dtau/3.d0 - dtau2*0.25d0  + dtau2*dtau*0.1d0
          gamma = dtau/3.d0 - dtau2/6.0d0   + dtau2*dtau*0.05d0
       else                  ! Accurate calculation
          eps = exp(-dtau)
          alpha = (dtau2 - 2.0d0*dtau + 2.0d0 - 2.0d0*eps) * idtau2
          beta  = (2.0d0 - (2.0d0 + 2.0d0*dtau + dtau2)*eps) * idtau2
          gamma = (2.d0*dtau - 4.0d0 + (2.d0*dtau+4.0d0)*eps) * idtau2
       endif       

       !
       ! Integration
       !
       cntrl = 0.5d0 * (S(k)  + dtau05*Sp(k) + S(k1) - dtau05 * Sp(k1))
       if(cntrl .GT. max(s(k), s(k1)) .OR. cntrl .LT. min(s(k), s(k1))) cntrl = S(k)

       I = I * eps + alpha * S(k) + beta * S(k1) + gamma * cntrl
    end do

    stokes = 0
    stokes(1) = I
  end SUBROUTINE delobezier_scal
  Subroutine hermite2(npoints, dtau_nu, Ab, S, stokes)
    USE Bezier_math
    Implicit None
    Integer :: ipoint, npoints, ibound, ii, jj, k, k1, info, ubound
    Integer :: ipiv(4)
    Real, dimension(npoints) :: dtau_nu, S, iab, tau_nu, tau_500
    Real, dimension(npoints, 4, 4) :: Ab, abp
    Real, dimension(4,npoints) ::  sv, opac, sp
    Real, dimension(4) :: stokes, Jb(4), Ja(4), sJa(4), sJb(4),  sKa(4,4), sKb(4,4)
    Real(8), dimension(4,4) :: ident, A, tmpa, tmpb, k2, D
    Real(8) :: dtau, dtau1, dtau2, idtau2, eps, dtau05
    Real(8) :: alpha, beta, gamma, lambda, ostokes(4), E(4)
    Real, Parameter :: Optically_thick = 100., Optically_thin=1.e-4

    !
    ! 4x4 Identity matrix
    !
    Data ident(:,1)/1.0d0,0.0d0,0.0d0,0.0d0/
    Data ident(:,2)/0.0d0,1.0d0,0.0d0,0.0d0/
    Data ident(:,3)/0.0d0,0.0d0,1.0d0,0.0d0/
    Data ident(:,4)/0.0d0,0.0d0,0.0d0,1.0d0/

  
    tau_nu(1)=0.
    do k=2, npoints
       tau_nu(k)=tau_nu(k-1)+dtau_nu(k)
    end do

    iab = ab(:,1,1)

    !
    ! Find the point for the boundary condition
    !
    ubound = 1
    ibound = npoints

    do k = 2, npoints-1
       If (Tau_Nu(k) .gt. (Optically_Thin) .and. Tau_Nu(k-1) .lt. (Optically_Thin)) &
            ubound=k-1
       If (Tau_Nu(k) .gt. (Optically_Thick) .and. Tau_Nu(k-1) .lt. (Optically_Thick)) &
            ibound=k
    End do
    if(ibound .lt. npoints) ibound = ibound + 1
    if(ubound .gt. 1) ubound = ubound - 1


    !
    ! Source function
    !
    Sv(2:4,:) = 0
    Sv(1,:) = S

       do jj = 1,4
          do ii = 1, 4
             ab(:,ii,jj) = ab(:,ii,jj) / iab
          end do
       end do



    !
    ! Set boundary condition (I=S) at ibound
    !
    ostokes(2:4)=0.0d0
    ostokes(1) = S(ibound) ! Boundary cond.

    !
    ! Derivatives of the source function and ab matrix
    !
    call full_cent_der_d(npoints, dtau_nu, Sv, sp, ubound, ibound)
    call full_cent_der4_d(npoints, dtau_nu, ab, abp, ubound, ibound)

    ibound = ibound - 1

    do k = ibound, ubound, -1
       !
       k1 = k + 1
       !
       ! dtau(s)
       !
       dtau = dtau_nu(k1) ! (tau_nu(k1) - tau_nu(k))
       dtau2 = dtau * dtau
       dtau05 = dtau * 0.5d0

       !
       ! at points "a" (k+1)
       !
       sKa  = matmul(ab(k1,:,:),ab(k1,:,:)) - abp(k1,:,:)
       sJa = dtau2/12.0d0 * (matmul(ab(k1,:,:), sp(:,k1)) - matmul(ska,Sv(:,k1)))
       Ja = dtau05 * matmul(ab(k1,:,:), Sv(:,k1))
       

       !
       ! at points "b" (k)
       !
       sKb = matmul(ab(k,:,:),ab(k,:,:)) - abp(k,:,:)
       sJb = dtau2/12.0d0 * (matmul(ab(k,:,:), sp(:,k)) - matmul(sKb,Sv(:,k)))
       Jb = dtau05 * matmul(ab(k,:,:), Sv(:,k))

       !
       ! Compute intensity
       !
       D = ident + dtau05*ab(k,:,:) + dtau2/12.d0 * sKb
       call matinx8(D)
       
       E = matmul(D, (Jb + Ja + sJa - sJb))
       A = matmul(D, (ident - dtau05 * ab(k1,:,:) + dtau2/12.d0 * sKa))
       
       ostokes = matmul(A, ostokes) + E

    END do

    stokes = ostokes
  END Subroutine hermite2
  subroutine matinx8( a )
    implicit None
    real(8) :: absmax, a(4,4), b(4,4), det, fdeta, fabsmx
    integer i , j

    absmax = 0.0d0

    absmax = maxval(abs(a))

    if ( absmax .eq. 0. ) then 
       print *,'singularity problem. Zero or NaN matrix D in Hermite'

       print *,'a='
       do i = 1,4
          print *, abs(a(:,i))
       end do
       do i=1, 4
          do j=1, 4
             a(i,j)=0.
          end do
          a(i,i)=1.
       end do
    end if


    fabsmx = 1.d00 / absmax


    do i = 1 , 4
       do j = 1 , 4
          a ( i , j ) = a ( i , j ) * fabsmx
       end do
    end do



    b(1,1) = a(2,2) * a(3,3) * a(4,4) + a(2,3) * a(3,4) * a(4,2) &
         + a(2,4) * a(3,2) * a(4,3) - a(2,2) * a(3,4) * a(4,3) &
         - a(2,3) * a(3,2) * a(4,4) - a(2,4) * a(3,3) * a(4,2)
    b(2,1) = a(2,3) * a(3,1) * a(4,4) + a(2,4) * a(3,3) * a(4,1) &
         + a(2,1) * a(3,4) * a(4,3) - a(2,3) * a(3,4) * a(4,1) &
         - a(2,4) * a(3,1) * a(4,3) - a(2,1) * a(3,3) * a(4,4)
    b(3,1) = a(2,4) * a(3,1) * a(4,2) + a(2,1) * a(3,2) * a(4,4) &
         + a(2,2) * a(3,4) * a(4,1) - a(2,4) * a(3,2) * a(4,1) &
         - a(2,1) * a(3,4) * a(4,2) - a(2,2) * a(3,1) * a(4,4)
    b(4,1) = a(2,1) * a(3,3) * a(4,2) + a(2,2) * a(3,1) * a(4,3) &
         + a(2,3) * a(3,2) * a(4,1) - a(2,1) * a(3,2) * a(4,3) &
         - a(2,2) * a(3,3) * a(4,1) - a(2,3) * a(3,1) * a(4,2)
    b(1,2) = a(3,2) * a(4,4) * a(1,3) + a(3,3) * a(4,2) * a(1,4) &
         + a(3,4) * a(4,3) * a(1,2) - a(3,2) * a(4,3) * a(1,4) &
         - a(3,3) * a(4,4) * a(1,2) - a(3,4) * a(4,2) * a(1,3)
    b(2,2) = a(3,3) * a(4,4) * a(1,1) + a(3,4) * a(4,1) * a(1,3) &
         + a(3,1) * a(4,3) * a(1,4) - a(3,3) * a(4,1) * a(1,4) &
         - a(3,4) * a(4,3) * a(1,1) - a(3,1) * a(4,4) * a(1,3)
    b(3,2) = a(3,4) * a(4,2) * a(1,1) + a(3,1) * a(4,4) * a(1,2) &
         + a(3,2) * a(4,1) * a(1,4) - a(3,4) * a(4,1) * a(1,2) &
         - a(3,1) * a(4,2) * a(1,4) - a(3,2) * a(4,4) * a(1,1)
    b(4,2) = a(3,1) * a(4,2) * a(1,3) + a(3,2) * a(4,3) * a(1,1) &
         + a(3,3) * a(4,1) * a(1,2) - a(3,1) * a(4,3) * a(1,2) &
         - a(3,2) * a(4,1) * a(1,3) - a(3,3) * a(4,2) * a(1,1)
    b(1,3) = a(4,2) * a(1,3) * a(2,4) + a(4,3) * a(1,4) * a(2,2) &
         + a(4,4) * a(1,2) * a(2,3) - a(4,2) * a(1,4) * a(2,3) &
         - a(4,3) * a(1,2) * a(2,4) - a(4,4) * a(1,3) * a(2,2)
    b(2,3) = a(4,3) * a(1,1) * a(2,4) + a(4,4) * a(1,3) * a(2,1) &
         + a(4,1) * a(1,4) * a(2,3) - a(4,3) * a(1,4) * a(2,1) &
         - a(4,4) * a(1,1) * a(2,3) - a(4,1) * a(1,3) * a(2,4)
    b(3,3) = a(4,4) * a(1,1) * a(2,2) + a(4,1) * a(1,2) * a(2,4) &
         + a(4,2) * a(1,4) * a(2,1) - a(4,4) * a(1,2) * a(2,1) &
         - a(4,1) * a(1,4) * a(2,2) - a(4,2) * a(1,1) * a(2,4)
    b(4,3) = a(4,1) * a(1,3) * a(2,2) + a(4,2) * a(1,1) * a(2,3) &
         + a(4,3) * a(1,2) * a(2,1) - a(4,1) * a(1,2) * a(2,3) &
         - a(4,2) * a(1,3) * a(2,1) - a(4,3) * a(1,1) * a(2,2)
    b(1,4) = a(1,2) * a(2,4) * a(3,3) + a(1,3) * a(2,2) * a(3,4) &
         + a(1,4) * a(2,3) * a(3,2) - a(1,2) * a(2,3) * a(3,4) &
         - a(1,3) * a(2,4) * a(3,2) - a(1,4) * a(2,2) * a(3,3)
    b(2,4) = a(1,3) * a(2,4) * a(3,1) + a(1,4) * a(2,1) * a(3,3) &
         + a(1,1) * a(2,3) * a(3,4) - a(1,3) * a(2,1) * a(3,4) &
         - a(1,4) * a(2,3) * a(3,1) - a(1,1) * a(2,4) * a(3,3)
    b(3,4) = a(1,4) * a(2,2) * a(3,1) + a(1,1) * a(2,4) * a(3,2) &
         + a(1,2) * a(2,1) * a(3,4) - a(1,4) * a(2,1) * a(3,2) &
         - a(1,1) * a(2,2) * a(3,4) - a(1,2) * a(2,4) * a(3,1)
    b(4,4) = a(1,1) * a(2,2) * a(3,3) + a(1,2) * a(2,3) * a(3,1) &
         + a(1,3) * a(2,1) * a(3,2) - a(1,1) * a(2,3) * a(3,2) &
         - a(1,2) * a(2,1) * a(3,3) - a(1,3) * a(2,2) * a(3,1)

    det = a ( 1 , 1 ) * b ( 1 , 1 ) + a ( 1 , 2 ) * b ( 2 , 1 ) &
         + a ( 1 , 3 ) * b ( 3 , 1 ) + a ( 1 , 4 ) * b ( 4 , 1 )

    fdeta = fabsmx / det



    do  i = 1 , 4
       do  j = 1 , 4
          a ( i , j ) = b ( i , j ) * fdeta
       enddo
    enddo
    return
  end subroutine matinx8

  
END MODULE Bezier_solvers
