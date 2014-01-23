! HSN: Modifications to Andres' routines
! Background opacity and scattering coefficient in cm^2/cm^3
!
!  Include atomic_data module in general_routines_mod and in 
!    background_opacity_module to supply Mg abundances consistent with
!    rest of the code. Changed also get_ionpot_abundance_mg accordingly
!
!  In background opacity, include scattering term as an output argument.
!     Also, take input arguments as default real and convert internally
!     to work with kind=8  
module constants_mod
implicit none

real(kind=8), parameter :: PK = 1.3806503d-16, PC = 2.99792458d10
real(kind=8), parameter :: PH = 6.62606876d-27, PI = 3.14159265359d0
real(kind=8), parameter :: PME = 9.10938188d-28, PE = 4.8032d-10
real(kind=8), parameter :: CC = 2.99792458D10

end module constants_mod


module general_routines_mod
use constants_mod
use atomic_data, only: at_abund, at_ioniz1, at_ioniz2
implicit none
contains

! ---------------------------------------------------------
! Given x(:) and y(:) which tabulate a function and the derivative at the boundary points
! this function returns the second derivative of the spline at each point
! ---------------------------------------------------------
  subroutine splin1(x,y,yp1,ypn,y2)
    real(kind=8), INTENT(IN) :: x(:), y(:)
    real :: yp1, ypn
    real(kind=8), INTENT(INOUT) :: y2(size(x))
    integer :: n, i, k
    real(kind=8) :: p, qn, sig, un, u(size(x))

    n = size(x)
    
    if (yp1 > .99d30) then
       y2(1) = 0.d0
       u(1) = 0.d0
    else
       y2(1) = -0.5d0
       u(1) = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif

    do i = 2, n-1
       sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
       p = sig * y2(i-1)+2.d0
       y2(i) = (sig-1.d0)/p
       u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/&
            (x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    if (ypn > .99d30) then
       qn = 0.d0
       un = 0.d0
    else
       qn = 0.5d0
       un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
    
    do k = n-1, 1, -1
       y2(k) = y2(k)*y2(k+1)+u(k)
    enddo
    
  end subroutine splin1
  
! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function, returns the interpolation using
! splines of vector x(:) in y(:)
! ---------------------------------------------------------
  subroutine spline(xa,ya,x,y)
    real(kind=8), INTENT(INOUT) :: y(:)
    real(kind=8), INTENT(IN) :: xa(:), ya(:), x(:)
    real(kind=8) :: y2a(size(xa))
    integer :: n_x, n, i, k, khi, klo
    real(kind=8) :: a, b, h, extrap
    
    n = size(xa)
    n_x = size(x)
    call splin1(xa,ya,1.e30,1.e30,y2a)
    
    do i = 1, n_x
       ! Downward extrapolation 
       if (x(i) < xa(1)) then
          !					y(i) = ya(1)
          y(i) = ya(1) + (ya(1)-ya(2))/(xa(1)-xa(2)) * (xa(1) - x(i))
       else 
          
          ! Upward extrapolation
          if (x(i) > xa(n)) then
             !					y(i) = ya(n)
             y(i) = ya(n) + (ya(n)-ya(n-1)) / (xa(n)-xa(n-1)) * (x(i) - xa(n))
          else
             ! In range
             klo = 1
             khi = n
1            if(khi-klo > 1) then
                k = (khi+klo)/2
                if (xa(k) > x(i)) then
                   khi = k
                else
                   klo = k
                endif
                go to 1
             endif
             
             h = xa(khi)-xa(klo)
             
             if (h == 0.d0) then
                print *, 'bad xa input in spline'
                stop
             endif
             a = (xa(khi)-x(i))/h
             b = (x(i)-xa(klo))/h
             
             y(i) = a*ya(klo)+b*ya(khi)+((a**3.d0-a)*y2a(klo)+(b**3.d0-b)*y2a(khi))*(h**2.d0)/6.d0
             
          endif
       endif
    enddo
    
  end subroutine spline
  
! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function z(nxmny), returns the interpolation using
! spline interpolation in 2D
! ---------------------------------------------------------		
  function interpol_2d(xa,ya,z,x0,y0)
    real(kind=8) :: interpol_2d
    real(kind=8) :: xa(:), ya(:), z(:,:), x0, y0, res(1), xin(1)
    real(kind=8), allocatable :: temp(:)
    integer :: nx, ny, j
    nx = size(xa)
    ny = size(ya)
    allocate(temp(nx))
    do j = 1, nx
       xin(1) = y0
       call spline(ya,z(j,:),xin,res)
       temp(j) = res(1)
    enddo
    
    xin(1) = x0
    call spline(xa,temp,xin,res)
    deallocate(temp)
    
    interpol_2d = res(1)
    
  end function interpol_2d
  
!-----------------------------------------------------------------
! This subroutine calculates the atomic partition function for a Mg
!	t: array with the temperature
!	u1: partition function for the neutral element (also an array of the same size as t)
!	u2: partition function for the 1st ionization (also an array of the same size as t)
!	u3: partition function for the 2nd ionization (also an array of the same size as t)	
!     A. D. Wittmann, Goettingen (1975).(1974, solar phys. 35, 11)
!-----------------------------------------------------------------
  subroutine partition_atomic_mg(t, u1, u2, u3)
    real(kind=8) :: t(:), u1(size(t)), u2(size(t)), u3(size(t))
    real(kind=8), allocatable :: x(:), y(:), coeff(:)
    integer :: n, i
    
    n = size(t)
    allocate(x(n))
    allocate(y(n))
    
    x = alog(5040.e0/real(t))
    y = 1.d-3*t
    
    u1 = 1.d0+exp(-4.027262-x*(6.173172+x*(2.889176+x*(2.393895+.784131*x))))
    where(t > 8.d3)
       u1 = 2.757+t*(-7.8909e-4+t*7.4531e-8)
    endwhere
    u2 = 2.d0+exp(-7.721172-x*(7.600678+x*(1.966097+.212417*x)))
    where(t > 2.d4)
       u2 = 7.1041+t*(-1.0817e-3+t*4.7841d-8)
    endwhere
    u3 = 1.d0
    
    deallocate(x)
    deallocate(y)
    
  end subroutine partition_atomic_mg
  
!-----------------------------------------------------------------
! Calculates the ratio between two ionization stages using the Saha equation
!-----------------------------------------------------------------
  function saha(t, pe, ul ,uu, chi)
    real(kind=8) :: t(:), pe(:), ul(:), uu(:), chi, saha(size(t))
    real(kind=8) :: constant

    constant = 2.d0 * (2.d0*PI*PME)**(1.5d0) / PH**3 * PK**(2.5d0)
    saha = constant * uu / ul * t**2.5d0 * 10.d0**(-5040.d0*chi/t) / pe
    
  end function saha
  
!-----------------------------------------------------------------
! Return the ionization potentials and abundances of Mg
!-----------------------------------------------------------------
  subroutine get_ionpot_abundance_mg(ionization_pot1, ionization_pot2, abundance)
    real(kind=8) :: ionization_pot1, ionization_pot2, abundance
    integer :: iel
    !
    iel=12 ! Mg atomic number
    
    !		ionization_pot1 = 7.64400d0
    !		ionization_pot2 = 15.0300d0
    !		abundance = 10.d0**(7.53-12.d0)
    ionization_pot1 = At_ioniz1(iel)
    ionization_pot2 = At_ioniz2(iel)
    abundance = 10.d0**(At_abund(iel)-12.d0)
    
  end subroutine get_ionpot_abundance_mg
  
end module general_routines_mod

module asensio_background_opacity_module
  use constants_mod
  use general_routines_mod
  use atomic_data
  use debug_module
  implicit none
  Private
  Public :: asensio_background_opacity
contains
  
  !-----------------------------------------------------------------
  ! Calculates the negative hydrogen (H-) free-free continuum absorption coefficient per
  ! hydrogen atom (cm^2)
  !  REF: John 1989 A&A 193, 189
  !  INPUT:
  !		T : temperature in K (1400 < T < 100080)
  !		Pe: electron pressure
  !		lambda: wavelength in A (greater than 1880 A)
  !-----------------------------------------------------------------
  function hminus_ff(T, Pe, lambda_in)
    real(kind=8) :: hminus_ff, T, Pe, lambda_in
    real(kind=8), dimension(6) :: a1, b1, c1, d1, e1, f1, com1
    real(kind=8), dimension(4) :: a2, b2, c2, d2, e2, f2, com2
    real(kind=8) :: lambda, theta, result
    
    
    if (lambda_in .lt. 1800) then 
       hminus_ff = 0.
       return
    end if

    a1 = (/0.d0,2483.346d0,-3449.889d0,2200.04d0,-696.271d0,88.283d0/)
    b1 = (/0.d0,285.827d0,-1158.382d0,2427.719d0,-1841.4d0,444.517d0/)
    c1 = (/0.d0,-2054.291d0,8746.523d0,-13651.105d0,8624.97d0,-1863.864d0/)
    d1 = (/0.d0,2827.776d0,-11485.632d0,16755.524d0,-10051.53d0,2095.288d0/)
    e1 = (/0.d0,-1341.537d0,5303.609d0,-7510.494d0,4400.067d0,-901.788d0/)
    f1 = (/0.d0,208.952d0,-812.939d0,1132.738d0,-655.02d0,132.985d0/)
    a2 = (/518.1021d0,473.2636d0,-482.2089d0,115.5291d0/)
    b2 = (/-734.8666d0,1443.4137d0,-737.1616d0,169.6374d0/)
    c2 = (/1021.1775d0,-1977.3395d0,1096.8827d0,-245.649d0/)
    d2 = (/-479.0721d0,922.3575d0,-521.1341d0,114.243d0/)
    e2 = (/93.1373d0,-178.9275d0,101.7963d0,-21.9972d0/)
    f2 = (/-6.4285d0,12.36d0,-7.0571d0,1.5097d0/)
    
    ! Transform wavelength to microns		
    lambda = lambda_in / 1.d4
    theta = 5040.d0 / T
    
    result = 0.d0
    
    if (lambda < 0.3645) then
       com2 = a2*lambda**2 + b2 + c2/lambda + d2/lambda**2 + e2/lambda**3 + f2/lambda**4
       result = com2(1) * theta**1 + com2(2)*theta**1.5 + com2(3)*theta**2 + com2(4)*theta**2.5
    else
       com1 = a1*lambda**2 + b1 + c1/lambda + d1/lambda**2 + e1/lambda**3 + f1/lambda**4
       result = com1(1) * theta**1 + com1(2)*theta**1.5 + com1(3)*theta**2 + com1(4)*theta**2.5 + &
            com1(5)*theta**3 + com1(6)*theta**3.5
    endif

    hminus_ff = 1.d-29 * result * Pe 
    
  end function hminus_ff
  
  !-----------------------------------------------------------------
  ! Calculates the negative hydrogen (H-) free-free continuum absorption coefficient per
  ! hydrogen atom (cm^2)
  !  REF: John 1989 A&A 193, 189
  !  INPUT:
  !		T : temperature in K (1400 < T < 100080)
  !		Pe: electron pressure
  !		lambda: wavelength in A (greater than 1880 A)
  !-----------------------------------------------------------------
  function hminus_table(T, Pe, lambda_in)
    real(kind=8) :: hminus_table, T, Pe, lambda_in
    real(kind=8), dimension(9,17) :: obs
    real(kind=8) :: l_obs(17), theta_obs(9), theta, l, result
    
    l_obs = (/0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.65,1.8,2.0,2.5,3.0/)
    theta_obs = (/0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0/)
    
    obs(:,1) = (/1.676,3.135,9.077,22.41,50.0,104.0,205.6,390.8,720.0/)
    obs(:,2) = (/2.850,5.355,15.54,38.38,85.65,178.2,352.2,669.3,1233.0/)
    obs(:,3) = (/3.985,7.548,22.05,54.55,121.8,253.4,500.8,951.7,1753.0/)
    obs(:,4) = (/4.925,9.423,27.83,69.12,154.5,321.7,635.9,1208.0,2226.0/)
    obs(:,5) = (/5.594,10.81,32.32,80.78,181.1,377.4,76.4,1419.0,2614.0/)
    obs(:,6) = (/5.960,11.61,35.16,88.53,199.3,416.1,823.7,1566.0,2887.0/) 
    obs(:,7) = (/6.033,11.83,36.22,91.90,207.8,435.1,862.7,1642.0,3027.0/)
    obs(:,8) = (/5.832,11.48,35.49,90.72,206.2,433.1,860.3,1639.0,3024.0/)
    obs(:,9) = (/5.400,10.64,33.11,85.18,194.5,410.0,816.4,1558.0,2877.0/)
    obs(:,10) = (/4.017,7.819,24.29,62.94,144.9,307.5,615.3,1178.0,2181.0/)
    obs(:,11) = (/2.332,4.254,12.49,31.86,73.21,155.7,312.5,600.3,1114.0/)
    obs(:,12) = (/0.998,1.321,2.318,4.244,8.038,15.37,29.17,54.38,99.32/)
    obs(:,13) = (/0.915,1.103,1.456,1.785,2.086,2.362,2.638,2.894,3.137/)
    obs(:,14) = (/1.083,1.308,1.727,2.117,2.470,2.800,3.122,3.424,3.708/)
    obs(:,15) = (/1.335,1.610,2.124,2.610,3.029,3.466,3.862,4.237,4.584/)
    obs(:,16) = (/2.076,2.494,3.283,4.029,4.677,5.379,5.984,6.570,7.122/)
    obs(:,17) = (/2.966,3.562,4.685,5.720,6.694,7.605,8.463,9.276,10.06/)
    
    theta = 5040.d0 / T
    l = lambda_in / 1.d4
    
    result = interpol_2d(theta_obs,l_obs,obs,theta,l)
    
    hminus_table = 1.d-27 * result * Pe
    
  end function hminus_table
  
  !-----------------------------------------------------------------
  ! Calculates the negative hydrogen (H-) bound-free continuum absorption coefficient per
  ! hydrogen atom (cm^2)
  !  REF: John 1989 A&A 193, 189 (with a couple of mistakes)
  !  INPUT:
  !		T : temperature in K
  !		Pe: electron pressure
  !		lambda: wavelength in A (greater than 1250 A)
  !-----------------------------------------------------------------
  function hminus_bf(T, Pe, lambda_in)
    real(kind=8) :: hminus_bf, T, Pe, lambda_in
    real(kind=8) :: lambda, lambda0, alpha, cte, cc(6), sigma, com, result
    
    lambda0 = 1.6419d0
    alpha = PH*PC/PK * 1.d4
    cte = 0.75d-18
    cc = (/152.519d0,49.534d0,-118.858d0,92.536d0,-34.194d0,4.982d0/)
    
    ! Transform wavelength to microns		
    lambda = lambda_in / 1.d4
    
    if (lambda < lambda0) then
       com = 1.d0/lambda - 1.d0/lambda0
       sigma = cc(1) + cc(2)*com**0.5 + cc(3)*com**1 + cc(4)*com**1.5 + cc(5)*com**2 + cc(6)*com**2.5
       sigma = sigma * lambda**3 * com**1.5
       result = T**(-2.5) * exp(alpha/(T*lambda0)) * (1.d0 - exp(-alpha/(T*lambda)))
       result = result * cte * sigma
    else
       result = 0.d0
    endif

    hminus_bf = result * Pe
    
  end function hminus_bf
  
  !-----------------------------------------------------------------
  ! Calculates the continuum absorption coefficient per
  ! hydrogen atom (cm^2) due to Thomson scattering (scattering with free electrons)
  !  INPUT:
  !		Pe: electron pressure
  !		PH: neutral H atoms partial pressure
  !-----------------------------------------------------------------
  function thomson(Pel, PH)
    real(kind=8) :: thomson, Pel, PH
    real(kind=8) :: cte
    
    cte = 8.d0*PI/3.d0 * (PE/PC)**4 / PME**2 ! cte=6.65243527767e987e-25
    thomson = cte * (Pel/PH)

  end function thomson
  
  !-----------------------------------------------------------------
  ! Calculates the hydrogen (H) bound-free and free-free continuum absorption coefficient per
  ! hydrogen atom (cm^2)
  !  REF: Landi Degl'Innocenti 1976, A&ASS, 25, 379
  !  INPUT:
  !		T : temperature in K (T < 12000)
  !		lambda: wavelength in A 
  !-----------------------------------------------------------------
  function hydrogen(T, lambda_in)
    real(kind=8) :: hydrogen, T, lambda_in
    real(kind=8) :: r, c1, c2, cte, theta1, theta2, theta3, sum, gff
    integer :: i, n0
    
    r = 1.096776d-3
    c1 = 1.5777216d5
    c2 = 1.438668d8
    cte = 1.045d-26
    
    theta1 = c1 / T
    theta2 = c2 / (lambda_in * T)
    theta3 = 2.d0*theta1
    
    ! Lowest level which can be photoionized
    n0 = 1 + floor(sqrt(r*lambda_in))
    
    ! Sum over states that can be photoionized
    if (n0 <= 8) then
       sum = exp(theta1 / n0**2) / n0**3
       do i = n0+1, 8
          sum = sum + exp(theta1 / i**2) / i**3
       enddo
       sum = sum + (0.117d0 + exp(theta1/81.d0)) / theta3
    else
       sum = (0.117d0 + exp(theta1/n0**2)) / theta3
    endif
    
    ! Approximate the value of the Gaunt factor G_ff from Mihalas Eq (80) @ theta=1, x=0.5
    gff = (1.d0-exp(-theta2)) * exp(-theta1) * lambda_in**3
    
    hydrogen = cte * gff * sum

  end function hydrogen
  
  !-----------------------------------------------------------------
  ! Calculates the Rayleigh scattering on hydrogen (H) per
  ! hydrogen atom (cm^2)
  !  REF: Landi Degl'Innocenti 1976, A&ASS, 25, 379 (he quotes Dalgarno (1962) but the polynomial has a mistake)
  !  INPUT:
  !		lambda: wavelength in A 
  !-----------------------------------------------------------------
  function rayleigh_h(lambda_in)
    real(kind=8) :: rayleigh_h, lambda_in
    real(kind=8) :: c(3)
    
    c = (/5.799d-13,1.422d-6,2.784d0/)
   
    rayleigh_h = (c(1)+(c(2)+c(3)/lambda_in**2) / lambda_in**2) / lambda_in**4
    
  end function rayleigh_h
  
  !-----------------------------------------------------------------
  ! Calculates the Rayleigh scattering on molecular hydrogen (H2) per
  ! hydrogen atom (cm^2)
  !  REF: Landi Degl'Innocenti 1976, A&ASS, 25, 379 (he quotes Dalgarno & Williams (1962) ApJ 136, 960)
  !  INPUT:
  !		lambda: wavelength in A 
  !		PH: hydrogen partial pressure
  !		PH2: molecular hydrogen partial pressure
  !-----------------------------------------------------------------
  function rayleigh_h2(PH, PH2, lambda_in)
    real(kind=8) :: rayleigh_h2, PH, PH2, lambda_in
    real(kind=8) :: c(3)
    
    c = (/8.14d-13,1.28d-6,1.61d0/)
    
    rayleigh_h2 = (c(1)+(c(2)+c(3)/lambda_in**2) / lambda_in**2) / lambda_in**4 * (PH2/PH)
    
  end function rayleigh_h2
  
  !-----------------------------------------------------------------
  ! Calculates the bound-free and free-free absorption coefficient of neutral Mg
  ! per H particle (not per H atom, so remember to multiply by nH_tot
  !  INPUT:
  !		T : temperature (K)
  !		Pe: electron pressure
  !		lambda: wavelength in A (from 1900 to 9000 A)
  !-----------------------------------------------------------------
  function magnesium(T, Pe, lambda_in)
    real(kind=8) :: magnesium, T, Pe, lambda_in
    real(kind=8) :: lambda(29), cMg(29,9), temin, temax, testep, result
    real(kind=8), dimension(1) :: temper, u1, u2, u3, Pelec, n1overn0, n2overn1, n0overn
    real(kind=8) :: ionization_pot1, ionization_pot2, abundance
    integer :: nlambda, ix(1), iy
    
    temin = 4000.d0
    temax = 12000.d0
    testep = 1000.d0
    nlambda = 29
    
    lambda = (/1.900d3,2.200d3,2.513d3,2.514d3,2.900d3,3.300d3,3.756d3,3.757d3,4.000d3,4.400d3,&
         4.884d3,4.885d3,5.200d3,5.504d3,5.505d3,6.000d3,6.549d3,6.550d3,6.900d3,7.234d3,7.235d3,&
         7.260d3,7.291d3,7.292d3,7.700d3,8.113d3,8.114d3,8.500d3,9.000d3/)
    
    ! Table equispaced in temperature from 4000 K to 12000 K in steps of 1000 K			
    cmg(:,1) = (/2.587d-20,3.910d-20,5.626d-20,7.319d-23,9.530d-23,&
         1.024d-22,9.561d-23,9.338d-24,1.097d-23,1.400d-23,1.831d-23,&
         1.594d-23,1.856d-23,2.131d-23,2.124d-23,2.617d-23,3.237d-23,&
         1.962d-23,2.245d-23,2.537d-23,1.778d-23,1.795d-23,1.815d-23,&
         1.748d-24,1.944d-24,2.162d-24,1.868d-24,2.116d-24,2.466d-24/)
    
    cmg(:,2) = (/1.251d-19,1.892d-19,2.723d-19,9.667d-22,1.268d-21,1.389d-21,&
         1.349d-21,2.756d-22,3.232d-22,4.117d-22,5.365d-22,4.908d-22,&
         5.718d-22,6.571d-22,6.556d-22,8.088d-22,1.002d-21,6.424d-22,&
         7.352d-22,8.307d-22,5.936d-22,5.990d-22,6.057d-22,8.915d-23,&
         1.005d-22,1.131d-22,1.029d-22,1.166d-22,1.360d-22/)
    
    cmg(:,3) = (/3.583d-19,5.425d-19,7.810d-19,5.606d-21,7.419d-21,8.323d-21,&
         8.449d-21,2.683d-21,3.147d-21,4.007d-21,5.217d-21,4.887d-21,&
         5.698d-21,6.553d-21,6.541d-21,8.079d-21,1.002d-20,6.692d-21,&
         7.660d-21,8.656d-21,6.304d-21,6.361d-21,6.433d-21,1.280d-21,&
         1.455d-21,1.648d-21,1.539d-21,1.745d-21,2.037d-21/)
    
    cmg(:,4) = (/7.611d-19,1.154d-18,1.661d-18,2.029d-20,2.711d-20,3.114d-20,&
         3.296d-20,1.381d-20,1.621d-20,2.064d-20,2.687d-20,2.552d-20,&
         2.978d-20,3.427d-20,3.422d-20,4.231d-20,5.252d-20,3.622d-20,&
         4.147d-20,4.687d-20,3.476d-20,3.507d-20,3.547d-20,8.834d-21,&
         1.009d-20,1.147d-20,1.088d-20,1.235d-20,1.443d-20/)
    
    cmg(:,5) = (/1.342d-18,2.037d-18,2.934d-18,5.454d-20,7.355d-20,8.633d-20,&
         9.477d-20,4.763d-20,5.594d-20,7.131d-20,9.288d-20,8.900d-20,&
         1.039d-19,1.197d-19,1.195d-19,1.479d-19,1.838d-19,1.301d-19,&
         1.490d-19,1.685d-19,1.271d-19,1.283d-19,1.297d-19,3.838d-20,&
         4.397d-20,5.018d-20,4.808d-20,5.459d-20,6.383d-20/)
    
    cmg(:,6) = (/2.090d-18,3.177d-18,4.579d-18,1.198d-19,1.630d-19,1.950d-19,&
         2.207d-19,1.258d-19,1.478d-19,1.886d-19,2.458d-19,2.370d-19,&
         2.769d-19,3.191d-19,3.187d-19,3.950d-19,4.913d-19,3.556d-19,&
         4.074d-19,4.607d-19,3.530d-19,3.562d-19,3.603d-19,1.222d-19,&
         1.403d-19,1.605d-19,1.549d-19,1.760d-19,2.059d-19/)
    
    cmg(:,7) = (/2.990d-18,4.548d-18,6.558d-18,2.282d-19,3.127d-19,3.804d-19,&
         4.415d-19,2.753d-19,3.237d-19,4.134d-19,5.392d-19,5.222d-19,&
         6.107d-19,7.042d-19,7.035d-19,8.727d-19,1.087d-18,8.017d-19,&
         9.187d-19,1.039d-18,8.080d-19,8.154d-19,8.247d-19,3.122d-19,&
         3.549d-19,4.118d-19,3.995d-19,4.541d-19,5.317d-19/)
    
    cmg(:,8) = (/4.015d-18,6.115d-18,8.819d-18,3.909d-19,5.391d-19,6.651d-19,&
         7.882d-19,5.253d-19,6.183d-19,7.904d-19,1.032d-18,1.003d-18,&
         1.173d-18,1.354d-18,1.353d-18,1.680d-18,2.093d-18,1.571d-18,&
         1.801d-18,2.039d-18,1.606d-18,1.621d-18,1.639d-18,6.794d-19,&
         7.834d-19,8.991d-19,8.765d-19,9.959d-19,1.167d-18/)
    
    cmg(:,9) = (/5.146d-18,7.844d-18,1.132d-17,6.177d-19,8.563d-19,1.069d-18,&
         1.290d-18,9.045d-19,1.065d-18,1.363d-18,1.781d-18,1.736d-18,&
         2.033d-18,2.346d-18,2.345d-18,2.915d-18,3.636d-18,2.769d-18,&
         3.177d-18,3.597d-18,2.868d-18,2.895d-18,2.928d-18,1.309d-18,&
         1.512d-18,1.738d-18,1.698d-18,1.932d-18,2.265d-18/)
    
    if (lambda_in < lambda(1) .or. lambda_in > lambda(29)) then
       result = 0.d0
    else
       ix = minloc(abs(lambda_in-lambda))
       iy = nint((T-temin) / testep) + 1
       if (T <= temin) then
          iy = 1
       endif
       if (T >= temax) then
          iy = 9
       endif
       result = cmg(ix(1),iy)
       temper = T
       Pelec = Pe
       call partition_atomic_mg(temper, u1, u2, u3)
       call get_ionpot_abundance_mg(ionization_pot1, ionization_pot2, abundance)
       n1overn0 = saha(temper, Pelec, u1, u2, ionization_pot1)
       n2overn1 = saha(temper, Pelec, u2, u3, ionization_pot2)
       n0overn = 1.d0 / (1.d0 + n1overn0 + n2overn1 * n1overn0)
       result = result * abundance * n0overn(1)
    endif
    
    magnesium = result
    
  end function magnesium
  
  !-----------------------------------------------------------------
  ! Calculates the background opacity due to various species
  ! INPUT :
  !	T : temperature in K
  !	Pe : electron pressure
  !	PH, PH2 : hydrogen and H2 partial pressures in dyn*cm^-2
  !	lambda_in : wavelength in A
  ! OUTPUT : 
  !  Background opacity in cm^-1 including the following contributions:
  !   H- free-free
  !   H- bound-free
  !   H bound-free + free-free
  !   Thomson scattering
  !   Rayleign scattering by H and H2
  !   Mg bound-free + free-free
  !-----------------------------------------------------------------
  function asensio_background_opacity(T4, Pe4, Pg4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4, Scat)
    real :: asensio_background_opacity, T4, Pe4, Pg4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4, Scat
    real :: nu, chi_0, chi_e, eta, tmp
    real(kind=8) :: T, Pe, Pg, PH, PHminus, PHplus, PH2, PH2plus, lambda_in
    real(kind=8) :: factor, H_per_volume
    real(kind=8) :: contrib1, contrib2, contrib3, contrib4, contrib5, contrib6, contrib7, caca
    real, Parameter :: Min_Pe=1e-6, Max_Pe=1e6
    real, Parameter :: Min_Lambda=300., Max_Lambda=120000
    Logical, Save :: FirstTime=.True.
    Character (Len=512) :: String
    ! Convert default Real to kind=8	

    T=T4
    Pe=Pe4
    PH=PH4
    PHminus=PHminus4
    PHplus=PHplus4
    PH2=PH24
    PH2plus=PH2plus4
    lambda_in=lambda_in4
    !		

    If (lambda_in .lt. Min_Lambda) then
       Call Debug_log('In asensio_background_opacity. Lambda .lt. Min_Lambda parameter. Clipping it.',2)
       lambda_in=Min_Lambda
    End if
    If (lambda_in .gt. Max_Lambda) then
       Call Debug_log('In asensio_background_opacity. Lambda .gt. Max_Lambda parameter. Clipping it.',2)
       lambda_in=Max_Lambda
    End if
    If (Pe .lt. Min_Pe) then
       Call Debug_log('In asensio_background_opacity. Pe .lt. Min_Pe parameter. Clipping it.',2)
       Pe=Min_Pe
    End if
    If (Pe .gt. Max_Pe) then
       Call Debug_log('In asensio_background_opacity. Pe .gt. Max_Pe parameter. Clipping it.',2)
       Pe=Max_Pe
    End if


    If (PH .gt. 1.e-6*(PH2+PH2plus)) then ! There is atomic H, act normally
       factor = 1.d0 + (PHminus + PHplus + 2.d0*(PH2+PH2plus)) / PH
       contrib1 = hminus_ff(T, Pe, lambda_in)
       contrib2 = hminus_bf(T, Pe, lambda_in)
       contrib4 = thomson(Pe, PH)
       contrib5 = rayleigh_h(lambda_in)
       contrib6 = rayleigh_h2(PH, PH2, lambda_in)
       contrib3=0.
       contrib7=0.
    Else ! All H is in molecular form. Cant compute things per H atom
       factor=1.
       contrib1 = 0.
       contrib2 = 0.
       contrib4 = thomson(Pe, 1.)
       contrib5 = 0.
       contrib6 = rayleigh_h(lambda_in)
       contrib3 = 0.
       contrib7 = 0.
    End if

    If (lambda_in .gt. 4000) then ! For UV, neutral H and Mg are computed in the UV package
       contrib7 = magnesium(T, Pe, lambda_in) * factor
       If (PH .gt. 1.e-6*(PH2+PH2plus)) then ! There is atomic H, act normally
          contrib3 = hydrogen(T, lambda_in)
          ! Magnesium opacity is per total H atom. For this reason, we multiply by factor
          ! because nH*factor=nH_tot
       Else ! All H is in molecular form. Cant compute things per H atom
          contrib3 = 0.
       End if
    End if
    
    If (PH .gt. 1.e-6*(PH2+PH2plus)) then ! There is atomic H, act normally
       asensio_background_opacity = PH / (PK*T) * &
            (contrib1 + contrib2 + contrib3 + contrib4 + contrib5 + contrib6 + contrib7)
       Scat=(contrib4+contrib5+contrib6)*PH/(PK*T)
    Else ! All H is in molecular form. Cant compute things per H atom
       asensio_background_opacity = 1. / (PK*T) * &
            (contrib1 + contrib2 + contrib3 + contrib4 + contrib5 + contrib6 + contrib7)
       Scat=1. / (PK*T) * (contrib4+contrib5+contrib6)
    End if

    If (Scat .lt. 0) then
       Write (String,*) 'T, Pe, PH, PHminus, PHplus, PH2, PH2plus, lambda, c1, c2, c3=', &
            T4,Pe4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4, &
            contrib4,contrib5,contrib6
      Call Debug_log('In asensio_background_opacity. Scat .lt. 0. Clipping it.'//String,2)
       Scat = 0.
    End if
   
    If (asensio_background_opacity .lt. 0.) then
       Write (String,*) 'T, Pe, PH, PHminus, PHplus, PH2, PH2plus, lambda=', &
            T4,Pe4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4
      Call Debug_log('In asensio_background_opacity. Opac .lt. 0. Clipping it.'//String,2)
       asensio_background_opacity = 0.
    End if
   
    If (asensio_background_opacity .gt. 1.e30) then
       Write (String,*) 'T, Pe, PH, PHminus, PHplus, PH2, PH2plus, lambda=', &
            T4,Pe4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4
       Call Debug_log('In asensio_background_opacity. Opac .gt. 1.e30 Clipping it'//String,2)
       asensio_background_opacity = 1.e30
    End if
 
  end function asensio_background_opacity

end module asensio_background_opacity_module
