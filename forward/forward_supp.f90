!
Module forward_support
  Use Phys_constants
  Use Atomic_data
  Use Line_data_structure
  Use Eq_state

Contains
! This routine checks an optical depth scale to make sure there are no
! two consecutive points with the same value. If such an occurrence is found
! then the tau scale from that point to the bottom is shifted slightly by 1%
! of the degenerate value. The corrected scale is returned in the same vector
! It is assumed that the tau scale is sorted from top to bottom, so tau(1) is the
! top of the atmosphere
!
! If you change this, check also misc/z_to_tau.f90
  Subroutine Check_tau(tau)
    Implicit None
    Real, Dimension(:) :: tau
    Integer :: npoints, ipoint
    
    npoints=Size(tau)
    Do ipoint=npoints-1, 1, -1
       If (abs(tau(ipoint)-tau(ipoint+1)) .lt. .01*tau(ipoint) ) then
          tau(ipoint+1:npoints)=tau(ipoint+1:npoints)+.01*tau(ipoint)
       End if
    End do
    
  End Subroutine Check_tau

! This subroutine constructs the absorption matrix at all the depth-points
! from the array of elements Phi_x and Psi_x. The matrix is normalized to
! the continuum opacity at 5000 Angstroms (which enters through the array
! Cont).
!
Subroutine Abs_matrix(npoints, nwlengths, Phi_I, Phi_Q, Phi_U, Phi_V, Psi_Q, &
             Psi_U, Psi_V, Cont, Absorp) ! Construct absorption matrix
  Implicit None
  Integer :: npoints, nwlengths, iwave
  Real, Dimension (nwlengths, npoints, 4, 4) :: Absorp
  Real, Dimension (nwlengths, npoints) :: Phi_I, Phi_Q, Phi_U, Phi_V
  Real, Dimension (nwlengths, npoints) :: Psi_Q, Psi_U, Psi_V
  Real, Dimension (npoints) :: Cont
!
  Do iwave=1, nwlengths
     Phi_I(iwave, :)=Phi_I(iwave, :)/Cont(:)
     Phi_Q(iwave, :)=Phi_Q(iwave, :)/Cont(:)
     Phi_U(iwave, :)=Phi_U(iwave, :)/Cont(:)
     Phi_V(iwave, :)=Phi_V(iwave, :)/Cont(:)
     Psi_Q(iwave, :)=Psi_Q(iwave, :)/Cont(:)
     Psi_U(iwave, :)=Psi_U(iwave, :)/Cont(:)
     Psi_V(iwave, :)=Psi_V(iwave, :)/Cont(:)
  End do
!
  Absorp(1:nwlengths,1:npoints,1,1)=Phi_I(1:nwlengths,1:npoints)
  Absorp(1:nwlengths,1:npoints,2,2)=Phi_I(1:nwlengths,1:npoints)
  Absorp(1:nwlengths,1:npoints,3,3)=Phi_I(1:nwlengths,1:npoints)
  Absorp(1:nwlengths,1:npoints,4,4)=Phi_I(1:nwlengths,1:npoints)
!
  Absorp(1:nwlengths,1:npoints,1,2)=Phi_Q(1:nwlengths,1:npoints)
  Absorp(1:nwlengths,1:npoints,2,1)=Phi_Q(1:nwlengths,1:npoints)
!  
  Absorp(1:nwlengths,1:npoints,1,3)=Phi_U(1:nwlengths,1:npoints)
  Absorp(1:nwlengths,1:npoints,3,1)=Phi_U(1:nwlengths,1:npoints)
!  
  Absorp(1:nwlengths,1:npoints,1,4)=Phi_V(1:nwlengths,1:npoints)
  Absorp(1:nwlengths,1:npoints,4,1)=Phi_V(1:nwlengths,1:npoints)
!  
  Absorp(1:nwlengths,1:npoints,2,3)=Psi_V(1:nwlengths,1:npoints)
  Absorp(1:nwlengths,1:npoints,3,2)=-Psi_V(1:nwlengths,1:npoints)
!
  Absorp(1:nwlengths,1:npoints,2,4)=-Psi_U(1:nwlengths,1:npoints)
  Absorp(1:nwlengths,1:npoints,4,2)=Psi_U(1:nwlengths,1:npoints)
!
  Absorp(1:nwlengths,1:npoints,3,4)=Psi_Q(1:nwlengths,1:npoints)
  Absorp(1:nwlengths,1:npoints,4,3)=-Psi_Q(1:nwlengths,1:npoints)
!

  Return
End Subroutine Abs_matrix

! This subroutine computes the damping parameter for the absorption profile.
! Temp, El_P and Pg and in cgs units
! If GA and GQ are .gt. 0, they are used to compute the Stark and
! radiative broadenings in the same manner as MULTI.
! If GA and GQ are .le. 0 then an approximation is used (see Gray)
!
Subroutine damping(Line, Temp, El_p, Pg, Dldop, Damp, GA, GQ)
  Use Debug_module
  Implicit None
  Type (Line_data) :: Line
  Real :: Temp, El_p, Dldop, Damp, C6, ioniz, chi_l, a, b, metal
  Real :: gamma_r, gamma_vdw, gamma_s, Pg, dldopHz, GA, GQ
  Real :: Sigma, Alpha, X, GX, GAMMAF, GVW, K, A0, M0, H1FRC, HE1FRC, VBAR
  Real, Dimension(1) :: tmp1, tmp2, tmp3, tmp4, tmp5, T1, ElP1, Pg1
  Real, Dimension (10) :: Pp
  Logical, Save :: Warning
  Data Warning/.FALSE./
  PARAMETER (K=1.380658E-23,M0=1.660540E-27)
  PARAMETER (A0=5.29177249E-11) ! Bohr radius
!
  dldopHz=cc/Line%Wlength/Line%Wlength*Dldop/1.E-8 ! Dldop in Hz
!
! Radiative damping:
  If (GA .le. 0) then 
! Reference: Gray 1976, "The observation and analysis of stellar 
!  photospheres", pag 227, just after Eq. (11-19). This approximation
!  is poor, but the radiative damping is usually negligible anyway.
!
     gamma_r=0.22233/(Line%Wlength*Line%Wlength*1.E-16)
  Else
     gamma_r=GA ! Same as MULTI
  End if
  If (Line%collisions .eq. 3 .and. gamma_r .ge. 0) then ! Explicit Radiative coefficient
     gamma_r = Line%Gamma_Rad
     gamma_r = gamma_r * 1D8
  End if
!
! Van der Waals damping:
  If (Line%collisions .eq. 1 .or. (Line%Collisions .eq. 3 .and. Line%Gamma_vdW_16 .lt. 0) ) then 
     ! Reference: Gray (see above, pag 239, Eqs. (11-35) and (11-36))
     ! Formula used: 
     !    log (gamma_vdw) = 19.6 + 2./5.*log (C6(H)) + log Pg - 7./10.*log T
     If (Line%Ion_stage .eq. 1) then
        ioniz=At_ioniz1(Line%Atomic_number)
     else if (Line%Ion_stage .eq. 2) then
        ioniz=At_ioniz2(Line%Atomic_number)
     else
        If (.not. Warning) &
             Print *,'Ionization stage gt 2. Van der Waals broadening not considered'
        ioniz=-1
        Warning=.TRUE.
     end if
     gamma_vdw=0.
     If (ioniz .gt. -1) then 
        chi_l=1.24E4/Line%Wlength
        a=(ioniz-Line%Energy_low-chi_l)**2.
        b=(ioniz-Line%Energy_low)**2.
        C6=3.E-31*(1./a-1./b)

        If (C6 .lt. 0) then
           gamma_vdw=0.
        Else
           gamma_vdw=19.6 + 2./5.*log10(C6) + log10(Pg) - 7./10.*log10(Temp)
           gamma_vdw=10.**(gamma_vdw)
        End if

        gamma_vdw=gamma_vdw*Line%VDW_enh ! Empirical Van der Waals enhancement
     End if
  Else if (Line%collisions .eq. 2) then ! Barklem formula
                           ! Following
                           ! http://www.astro.uu.se/~barklem/howto.html
     Sigma=Line%Bark_sigma*A0*A0
     Alpha=Line%Bark_alpha
    !  Compute the Gamma function of X, this function is valid over the 
    !  range 1<X<2 ie. 0<ALPHA<2 which is always satisfied    
     X=2.-ALPHA*.5
     GX=X-1.0
     GAMMAF=1+(-.5748646+(.9512363+(-.6998588+(.4245549-.1010678*GX &
          )*GX)*GX)*GX)*GX
    ! Compute the halfwidth per unit perturber number density for this temp
     GVW=(4./PI)**(ALPHA*0.5)*GAMMAF*1.E4*SIGMA
     VBAR=SQRT(8.*K*Temp/PI/M0*(1./1.008+1./Line%Atomic_weight))
     GVW=GVW*((VBAR/1.E4)**(1.-ALPHA))
    ! Get H and He partial pressures
     T1=Temp
     ElP1=El_p
     Pg1=Pg
     Call compute_others_from_T_Pe_Pg(1,T1,Elp1, Pg1, tmp1,tmp2,tmp3,tmp4,tmp5)
    ! Fullwidth given H1FRC perturbers per cm^3 with approx He I broadening
    ! The factor of 2 comes from converting to the full width.
    ! The factor of 1.E6 converts from SI to cgs
     H1FRC=tmp1(1) ! nH=neutral H
     HE1FRC=0.1*H1FRC ! Approx neutral He by 0.10*neutral H

     GVW=GVW*(H1FRC+ 0.42*HE1FRC)*1.E6*2.
     gamma_vdw=GVW
     gamma_vdw=gamma_vdw*Line%VDW_enh ! Empirical Van der Waals enhancement
  Else If (Line%Collisions .ne. 3) then
     Print *,'Unknown collisional broadening in damping (forward_supp.f90)'
     Stop
  Endif
  If (Line%collisions .eq. 3) then ! Explicit vdW coefficient
    ! Get H partial pressures
     T1=Temp
     ElP1=El_p
     Pg1=Pg
     Call compute_others_from_T_Pe_Pg(1,T1,Elp1, Pg1, tmp1,tmp2,tmp3,tmp4,tmp5)
     If (Line%Gamma_vdW_16 .gt. 0) then
        gamma_vdw=Line%Gamma_vdW_16*(tmp1(1)/1D16)/(1D4**0.38)*(Temp**0.38)
        gamma_vdw = gamma_vdw * 1D8
     End if
     gamma_vdw=gamma_vdw*Line%VDW_enh ! Empirical Van der Waals enhancement!
  End if

! Stark damping
  If (GQ .le. 0) then ! debug
! Formula used: gamma_4=38.8*(v**(1./3.))*(C4**(2./3.))*N , from
!   Unsold 1955, "Physik der Sternatmospharen", 2nd ed., 
!   Springer-Verlag, pp. 326 and 331. According to Gray (ref above), 
!   pag 237, this is similar to 
!   log (gamma_4) = 19.4 + 2./3*log C4 + log Pe - 5./6.*log T
!   The value of log C4 is approximated by -14. (see Gray, Table 11-3)
     gamma_s = 19.4 + 2./3.*(-14.) + log10(El_p) - 5./6.*log10(Temp)
     gamma_s = 10.**(gamma_s)
  Else
!     gamma_s = GQ*El_p/bk/Temp ! Same as MULTI
     gamma_s = GQ*((El_p/bk/Temp)**0.6666666)*4.*Pi*0.426*0.6*(3**2-2**2) ! Same as MULTI for i=2 to i=3 transition
  End if
  If (Line%collisions .eq. 3 .and. gamma_s .ge. 0) then ! Explicit Stark coefficient
     gamma_s = Line%Gamma_Strk_12*(El_P/BK/Temp/1D12)/(1D4**0.17)*(Temp**0.17)
     gamma_s = gamma_s * 1D8
  End if
  
!  if (el_p .gt. 38 .and. el_p .lt. 80) then ! debug
!     print *,'gamma=',gamma_r,gamma_vdw,gamma_s
!     pause
!  endif


  Damp=(gamma_r+gamma_vdw+gamma_s)/4./Pi/dldopHz

  Return
End Subroutine damping

!
subroutine matinx ( a )

!	'exact' inversion of 4 X 4 matrix
  Use Debug_Module ! Optional
  implicit real ( a-h, o-z )

  dimension a ( 4 , 4 ) , b ( 4 , 4 )
!      dimension c(4,4)
  integer i , j

  absmax = 0.
  
  do i = 1 , 4
     do j = 1 , 4
        !          c(i,j)=a(i,j) 
        if ( abs ( a ( i , j ) ) .gt. absmax ) absmax = a ( i , j )
     end do
  end do

  if ( absmax .eq. 0. ) then 
     print *,'singularity problem. Zero or NaN matrix D in Hermite'
     print *,'a=',a
!     stop
     do i=1, 4
        do j=1, 4
           a(i,j)=0.
        end do
        a(i,i)=1.
     end do

!   Only if optional module debug is used
     Call Debug_Log('Singularity problem',1)

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
  !	if(abs(det).lt.1.e-2)then
  !		print*,det
  !                do i=1,4
  !                   do j=1,4
  !                       print*,'c(',i,',',j,')=',sngl(c(i,j))  
  !                   end do
  !                end do
  !	end if
  
  !      print*,'en matinx (4)',det
  
  fdeta = fabsmx / det
  
  !      print*,'en matinx (5)',fdeta
  
  
  do 5000 i = 1 , 4
     do 5001 j = 1 , 4
        a ( i , j ) = b ( i , j ) * fdeta
5001    continue
5000    continue
        !      	if(abs(det).lt.1.e-2)then
        !          do i=1,4
        !             do j=1,4
        !                sum=0.
        !                do l=1,4
        !                    sum=sum+a(i,l)*c(l,j)
        !                end do
        !                print*,'c*inversa(',i,',',j,')=',sum
        !             end do
        !           end do
        !         end if 
        
        !      print*,'en matinx (6)',fdeta
        
        return
end subroutine matinx

!______________________________________________________________________________
     subroutine matinx4 ( b )
!______________________________________________________________________________

       implicit real ( a-h, o-z )
       
       dimension b ( 4 , 4 ) 
       
       q=b(1,2)
       u=b(1,3)
       v=b(1,4)
       r=b(2,3)
       s=b(4,2)
       t=b(3,4)
       q2=q*q
       u2=u*u
       v2=v*v
       r2=r*r
       s2=s*s
       t2=t*t
       a=q*t+r*v+s*u
       a2=a*a
       taq=t*a+q
       qat=q*a-t
       sau=s*a+u
       uas=u*a-s
       rav=r*a+v
       var=v*a-r
       ur=u*r-s*v
       uv=u*v+r*s
       qr=q*r-t*v
       qs=q*s-t*u
       qu=q*u+t*s
       qv=q*v+t*r
       
       det=1.+t2+r2+s2-q2-u2-v2-a2
       
       b(1,1)=(1.+t2+r2+s2)/det
       b(2,1)=(ur-taq)/det
       b(3,1)=(-qr-sau)/det
       b(4,1)=(qs-rav)/det
       b(1,2)=(-ur-taq)/det
       b(2,2)=(1.+t2-u2-v2)/det
       b(3,2)=(qu-var)/det
       b(4,2)=(qv+uas)/det
       b(1,3)=(qr-sau)/det
       b(2,3)=(qu+var)/det
       b(3,3)=(1.+s2-q2-v2)/det
       b(4,3)=(uv-qat)/det
       b(1,4)=(-qs-rav)/det
       b(2,4)=(qv-uas)/det
       b(3,4)=(uv+qat)/det
       b(4,4)=(1.+r2-q2-u2)/det
       
       return
     end subroutine matinx4
!
! This routine performs a matrix product A.B=C. The dimension of A are nrA
! by ncA (number of rows and columns, respectively), while ncB is the number
! of columns in B.
!
Subroutine Multiply_matrix(nrA, ncA, ncB, A, B, C)
  Implicit None
  Integer :: nrA, ncA, ncB
  Real, Dimension (nrA, ncA) :: A
  Real, Dimension (ncA, ncB) :: B
  Real, Dimension (nrA, ncB) :: C
  Integer :: i, j, k
!
  Do i=1, nrA
     Do j=1, ncB
        C(i, j)=0.
        Do k=1, ncA
           C(i, j)=C(i, j)+A(i, k)*B(k, j)
        End do
     End do
  End do
  Return
End Subroutine Multiply_matrix


! This subroutine calculates the derivatives of a function y assuming a
! parabolic behavior between i-1, i and i+1 (being i the point where the
! derivative is to be evaluated). The derivatives at the boundaries are
! computed assuming a linear behavior. The vector x does not need to be
! equispaced.
!
Subroutine Parab_deriv(npoints, x, y, dy)
  Implicit None
  Integer :: npoints, i
  Real, Dimension (npoints) :: x, y, dy
  Real :: a, b, c, den
!
! Linear interpolation for the boundaries
!
  dy(1)=(y(2)-y(1))/(x(2)-x(1))
  dy(npoints)=(y(npoints)-y(npoints-1))/(x(npoints)-x(npoints-1))
!
! Parabolic interpolation for all the other points
!
  Do i=2,npoints-1
     den=(x(i)-x(i-1))*(x(i)-x(i+1))*(x(i-1)-x(i+1))
     a = (x(i+1)*(y(i-1)-y(i)) + x(i-1)*(y(i)-y(i+1)) + &
          x(i)*(y(i+1)-y(i-1)))/den
     b = (x(i+1)*x(i+1)*(y(i)-y(i-1)) + x(i)*x(i)*(y(i-1)-y(i+1)) + &
          x(i-1)*x(i-1)*(y(i+1)-y(i)))/den
!
!    y = a x^2 + bx + c => y' = 2 a x + b
!
     dy(i)=2.*a*x(i)+b
  End do
  Return
End Subroutine Parab_deriv

End Module forward_support
