! This subroutine takes a model and determines the electron pressure,
! gas pressure, density and the height scale by solving the hydrostatic
! equilibrium equation. The boundary condition for this equation is the
! electron pressure at the top of the atmosphere and is taken from the
! model Atmo.
!
Subroutine Hydrostatic(Params, Atmo)
  Use Param_structure
  Use Model_structure
  Use Eq_state
  Use Atomic_data
  Use Background_opacity_module, Only: Background_opacity
  Use Debug_module
  Use Profiling
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Atmo
  Type (Model), Save :: Saved_Atmo
  Integer :: ipoint, npoints, iatom, iters
  Integer, Dimension(1) :: imin
  Real, Parameter :: niters=30, Precision=1e-5
  Real, Dimension (Params%n_points) :: Kappa, Tau, temp
  Real, Dimension (10) :: Pp
  Real, Parameter :: nu500 = cc/(5000.*1.e-8), Min_temp = 2500., &
       Mu=12.566370614 ! Vacuum permeability (G^2 cm^3/erg)=4*Pi, P_m=B^2/2/Mu
  Real :: dtau, dif, n2P, Scat, chi_0, chi_e, eta
  Real :: Avmolweight, Asum, Wsum, Pg_Old, OldKappa, metal
  Logical :: Warning1, Warning2
  Logical, Save :: FirstTime=.True.
!
  Debug_warningflags(flag_hydrostatic)=0
  Debug_errorflags(flag_hydrostatic)=0
  Call time_routine('hydrostatic',.True.)
! Do we need to redo the calculation?
  If (FirstTime) then
     Allocate(Saved_Atmo%temp(Params%n_points))
     Allocate(Saved_Atmo%z_scale(Params%n_points))
     Allocate(Saved_Atmo%gas_p(Params%n_points))
     Allocate(Saved_Atmo%rho(Params%n_points))
     Allocate(Saved_Atmo%el_p(Params%n_points))
     Allocate(Saved_Atmo%ne(Params%n_points))
     Allocate(Saved_Atmo%nH(Params%n_points))
     Allocate(Saved_Atmo%nHplus(Params%n_points))
     Allocate(Saved_Atmo%nHminus(Params%n_points))
     Allocate(Saved_Atmo%nH2(Params%n_points))
     Allocate(Saved_Atmo%nH2plus(Params%n_points))
  Else
     If (MaxVal(Abs(Atmo%temp-Saved_Atmo%temp)) .lt. 10) then
        Atmo%z_scale=Saved_Atmo%z_scale
        Atmo%gas_p=Saved_Atmo%gas_p
        Atmo%rho=Saved_Atmo%rho
        Atmo%el_p=Saved_Atmo%el_p
        Atmo%ne=Saved_Atmo%ne
        Atmo%nH=Saved_Atmo%nH
        Atmo%nHplus=Saved_Atmo%nHplus
        Atmo%nHminus=Saved_Atmo%nHminus
        Atmo%nH2=Saved_Atmo%nH2
        Atmo%nH2plus=Saved_Atmo%nH2plus
        Call time_routine('hydrostatic',.False.)
        Return
     End if
  End if
!
  FirstTime=.False.
  npoints=Params%n_points
  Warning1=.FALSE.
  Warning2=.FALSE.
  Tau=10.**Atmo%ltau_500
!
! Calculate the average molecular weight of the particles
!
  Wsum=0.
  Asum=0.
  Do iatom=1, 92
     Wsum=Wsum+At_weight(iatom)*10**(At_abund(iatom)-12.)
     Asum=Asum+10**(At_abund(iatom)-12.)
  End do
  metal=At_abund(26)-7.5
!
! Boundary condition
!
  temp(1)=Atmo%Temp(1)
  If (temp(1) .lt. Min_temp) then ! Check temperature so that gasc won't crash
     Warning1=.TRUE.
     temp(1)=Min_temp
  End if

  !
  ! JdlCR: change the boundary condition from "preserve electron pressure"
  ! to "preserve gas pressure". It seems to solve the problem of dual
  ! solutions in the upper layers when running chromospheric inversions. 
  !
  !  Call Compute_Pg(1, temp(1), Atmo%El_p(1), Atmo%Gas_p(1))

  If (Params%Input_dens .ne. 'pel')  &
       Call Compute_Pe(1, temp(1), Atmo%Gas_p(1), Atmo%El_p(1))
  

  n2P=BK*temp(1)
  Call Compute_others_from_T_Pe_Pg(1,Temp(1), Atmo%El_p(1), Atmo%Gas_p(1), Atmo%nH(1), &
       Atmo%nHminus(1), Atmo%nHplus(1), Atmo%nH2(1), Atmo%nH2plus(1))
  Kappa(1)=Background_opacity(Temp(1), Atmo%El_p(1), Atmo%Gas_p(1), Atmo%nH(1)*n2P, &
       Atmo%nHminus(1)*n2P, Atmo%nHplus(1)*n2P, Atmo%nH2(1)*n2P, Atmo%nH2plus(1)*n2P,&
       5000., Scat)
  Avmolweight=Wsum/(Asum+ &
       Atmo%El_P(1)/Atmo%Gas_P(1))
  Atmo%Rho(1)=Atmo%Gas_p(1)*Avmolweight/Avog/bk/temp(1) ! Gas density
  Kappa(1)=Kappa(1)/Atmo%Rho(1) ! Convert to cm^2/g
  Atmo%Z_scale(1)=0. ! Boundary condition for the height scale
  Do ipoint=2, npoints
     dtau=10.**Atmo%ltau_500(ipoint) - 10.**Atmo%ltau_500(ipoint-1)
     temp(ipoint)=Atmo%Temp(ipoint)
     If (temp(ipoint) .lt. Min_temp) then ! Check temperature so that gasc won't crash
        temp(ipoint)=Min_temp
        Warning1=.TRUE.
     End if
! Initialization
     Kappa(ipoint)=Kappa(ipoint-1)
     dif=1
     iters=0
! Now iterate to find consistent kappa, gas_p, el_p
     Do while (iters .lt. niters .and. dif .gt. Precision)
!       Hydrostatic equilibrium equation.
        Atmo%Gas_p(ipoint)=Atmo%Gas_p(ipoint-1) + &
             Gravity*dtau/(.5*(Kappa(ipoint-1)+Kappa(ipoint))) ! + & 
!            (Atmo%B_str(ipoint-1)**2  - &
!             Atmo%B_str(ipoint)**2) /8./Pi ! Magnetic pressure term
        Call Compute_Pe(1, Temp(ipoint), Atmo%Gas_p(ipoint), Atmo%El_p(ipoint))
        OldKappa=Kappa(ipoint)
        Call Compute_others_from_T_Pe_Pg(1,Temp(ipoint), Atmo%El_p(ipoint), Atmo%Gas_p(ipoint), Atmo%nH(ipoint), &
             Atmo%nHminus(ipoint), Atmo%nHplus(ipoint), Atmo%nH2(ipoint), Atmo%nH2plus(ipoint))
        Kappa(ipoint)=Background_opacity(Temp(ipoint), Atmo%El_p(ipoint), Atmo%Gas_p(ipoint), Atmo%nH(ipoint)*n2P, &
             Atmo%nHminus(ipoint)*n2P, Atmo%nHplus(ipoint)*n2P, Atmo%nH2(ipoint)*n2P, Atmo%nH2plus(1)*n2P,&
             5000., Scat)
        Avmolweight=Wsum/(Asum+ &
             Atmo%El_p(ipoint)/Atmo%Gas_p(ipoint))
        Atmo%Rho(ipoint)=Atmo%Gas_p(ipoint)*Avmolweight/Avog/bk/temp(ipoint) ! Gas density
        Kappa(ipoint)=Kappa(ipoint)/Atmo%Rho(ipoint) ! Convert to cm^2/g
        Atmo%Z_scale(ipoint)=Atmo%Z_scale(ipoint-1) - &
             dtau/2./1.e5* &
             (1./(Kappa(ipoint)*Atmo%Rho(ipoint))+1./(Kappa(ipoint-1)*Atmo%Rho(ipoint-1)))
        dif=Abs(OldKappa-Kappa(ipoint))/(OldKappa+Kappa(ipoint))
        iters=iters+1
     End do
     Call Compute_Pe(1, Temp(ipoint), Atmo%Gas_p(ipoint), Atmo%El_p(ipoint))
  End Do ! ipoint
  If (dif .gt. Precision) Warning2=.TRUE.
  If (Warning1 .or. Warning2) then 
     If (Debug_level .ge. 1) then
        If (Debug_FileUnit .lt. 0) then
           Print *,'Error in hydrostatic.f90 writing debug info. This should not happen!!'
        Else
           If (Warning1) &
                Write (Debug_FileUnit,*) '  ** Clipped T in forward/hydrostatic'
           If (Warning2) &
                Write (Debug_FileUnit,*) '  ** Error in forward/hydrostatic'
           Write (Debug_FileUnit,*) '     Atmosphere that caused the error follows'
           Write (Debug_FileUnit,*) '     ------   Z   ------'
           Write (Debug_FileUnit,*) Atmo%Z_scale
           Write (Debug_FileUnit,*) '     ------   tau_500  ------'
           Write (Debug_FileUnit,*) Atmo%ltau_500
           Write (Debug_FileUnit,*) '     ------   T  ------'
           Write (Debug_FileUnit,*) Atmo%Temp
           Write (Debug_FileUnit,*) '     ------   Electron pressure  ------'
           Write (Debug_FileUnit,*) Atmo%El_p
           Write (Debug_FileUnit,*) '     ------   Gas pressure  ------'
           Write (Debug_FileUnit,*) Atmo%Gas_p
           Write (Debug_FileUnit,*) '     ------   Density  ------'
           Write (Debug_FileUnit,*) Atmo%Rho
           Call Debug_Log('Error in forward/hydrostatic',2)
           Debug_warningflags(flag_hydrostatic)=1
        End if
     End if
  End if

  imin=MinLoc(Abs(Atmo%ltau_500))
  Atmo%Z_scale(1:Params%n_points)=Atmo%Z_scale(1:Params%n_points) - &
       Atmo%Z_scale(imin(1))
  Atmo%ne(1:Params%n_points)=Atmo%el_p(1:Params%n_points)/ &
       bk/Atmo%Temp(1:Params%n_points)
!
  Saved_Atmo%temp=Atmo%temp
  Saved_Atmo%z_scale=Atmo%z_scale
  Saved_Atmo%gas_p=Atmo%gas_p
  Saved_Atmo%rho=Atmo%rho
  Saved_Atmo%el_p=Atmo%el_p
  Saved_Atmo%ne=Atmo%ne
  Saved_Atmo%nH=Atmo%nH
  Saved_Atmo%nHplus=Atmo%nHplus
  Saved_Atmo%nHminus=Atmo%nHminus
  Saved_Atmo%nH2=Atmo%nH2
  Saved_Atmo%nH2plus=Atmo%nH2plus
  Call time_routine('hydrostatic',.False.)

  Return
End Subroutine Hydrostatic
