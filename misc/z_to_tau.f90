! This routine takes a model in geometrical heights (with
! z(km), Gas_p(dyn/cm^2) and rho (g/cm^3) and fills in the vector
! log(tau_5000). Note that el_p, Gas_p and rho must be known before calling
! this routine, as well as nH, nHminus, nHplus and nH2
!
Subroutine Z_to_tau(Params, Atmo)
  Use Param_structure
  Use Eq_state
  Use Model_structure
  Use Atomic_data
  Use Background_opacity_module, Only: Background_opacity
  Use Debug_module
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Atmo
  Integer :: ipoint, idx, iatom
  Real :: Avmolweight, Asum, Wsum, metal
  Real :: T, Pe, Pg, Rho, dtau, Chi_0, Chi_e, Eta, H
  Real :: PH, PHminus, PHplus, PH2, PH2plus, Scat, n2P
  Real, Dimension(10) :: pp
  Real, Dimension(Params%n_points) :: Kappa
  Real, Parameter :: nu500 = cc/(5000.*1.e-8)
  Logical :: Warning
  Warning=.False.
!
! Check z scale for errors
!
  If (MaxVal(Abs(Atmo%Z_scale)) .lt. 1.e-3) then
     If (Params%Printout .ge. 0) then
        Print *,'Error. Something is wrong with the Z scale in the model'
        Print *,'It spans from ',MinVal(Atmo%Z_scale),' to ',MaxVal(Atmo%Z_scale)
        Stop
     End if
  End if
  If (MinVal(Atmo%Z_scale) .gt. 0 .or. MaxVal(Atmo%Z_scale) .lt. 200) then
     If (Params%Printout .ge. 1) then
        Print *,'Warning. There seems to be something weird with the Z scale in the model'
        Print *,'It spans from ',MinVal(Atmo%Z_scale),' to ',MaxVal(Atmo%Z_scale)
        Print *,'Proceeding anyway (hope you know what you are doing). If the results'
        Print *,'are not as expected, this might be the reason'
     End if
  End if
!
! Check el_p, gas_p and rho at z=0
!
  ipoint=-1
  Do idx=1, Params%n_points
     If (Abs(Atmo%Z_scale(idx))-MinVal(Abs(Atmo%Z_scale)) .lt. 1e-5) then
        ipoint=idx
     end if
  End do
  If (ipoint .eq. -1) then 
     Print *,'There seems to be something wrong with the Z-scale. Please, check it'
     Stop
  End if
  If (Atmo%El_p(ipoint) .lt. 1e-2*56.44 .or. Atmo%El_p(ipoint) .gt. 1e2*56.44) then
     If (Params%Printout .ge. 1) then
        Print *,'Warning. Electron pressure near z=0 is way off typical solar values'
        Print *,'Found (cgs):',Atmo%El_p(ipoint),'  HSRA has:',56.44
        Print *,'Proceeding anyway (hope you know what you are doing). If the results'
        Print *,'are not as expected, this might be the reason'
     End if
  End if
  If (Atmo%Gas_p(ipoint) .lt. 1e-2*1.13E5 .or. Atmo%Gas_p(ipoint) .gt. 1e2*1.31E5) then
     If (Params%Printout .ge. 1) then
        Print *,'Warning. Gas pressure near z=0 is way off typical solar values'
        Print *,'Found (cgs):',Atmo%Gas_p(ipoint),'  HSRA has:',1.31E5
        Print *,'Proceeding anyway (hope you know what you are doing). If the results'
        Print *,'are not as expected, this might be the reason'
     End if
  End if
  If (Atmo%Rho(ipoint) .lt. 1e-2*3.19E-7 .or. Atmo%Rho(ipoint) .gt. 1e2*3.19E-7) then
     If (Params%Printout .ge. 1) then
        Print *,'Warning. Density near z=0 is way off typical solar values'
        Print *,'Found (g/cm3):',Atmo%Rho(ipoint),'  HSRA has:',3.19E-7
        Print *,'Proceeding anyway (hope you know what you are doing). If the results'
        Print *,'are not as expected, this might be the reason'
     End if
  End if
!
! First depth-point
!
  n2P=BK*Atmo%Temp(1)
  metal=Atmo%Abundance(26)-7.5

  Call Compute_others_from_T_Pe_Pg(1, Atmo%Temp(1), &
       Atmo%El_p(1), Atmo%Gas_p(1), Atmo%nH(1), &
       Atmo%nHminus(1), Atmo%nHplus(1), &
       Atmo%nH2(1), Atmo%nH2plus(1))
  Kappa(1)=Background_opacity(Atmo%Temp(1), Atmo%El_p(1), Atmo%Gas_p(1), Atmo%nH(1)*n2P, &
       Atmo%nHminus(1)*n2P, Atmo%nHplus(1)*n2P, Atmo%nH2(1)*n2P, &
       Atmo%nH2plus(1)*n2P, 5000., Scat)
  Kappa(1)=Kappa(1)/Atmo%Rho(1) ! Convert to cm^2/g

!
! Boundary condition: tau=1e-9 at the top
  Atmo%ltau_500(1)=1.e-9
!
! Conversion at each depth-point
!
  Do ipoint=2, Params%n_points
     n2P=BK*Atmo%Temp(ipoint)
     Call Compute_others_from_T_Pe_Pg(1, Atmo%Temp(ipoint), &
          Atmo%El_p(ipoint), Atmo%Gas_p(ipoint), Atmo%nH(ipoint), &
          Atmo%nHminus(ipoint), Atmo%nHplus(ipoint), &
          Atmo%nH2(ipoint), Atmo%nH2plus(ipoint))
     Kappa(ipoint)=Background_opacity(Atmo%Temp(ipoint), Atmo%El_p(ipoint), Atmo%Gas_p(ipoint), &
          Atmo%nH(ipoint)*n2P, Atmo%nHminus(ipoint)*n2P, Atmo%nHplus(ipoint)*n2P,&
          Atmo%nH2(ipoint)*n2P, Atmo%nH2plus(ipoint)*n2P, 5000., Scat)
     Kappa(ipoint)=Kappa(ipoint)/Atmo%Rho(ipoint) ! Convert to cm^2/g
     dtau=.5*(Kappa(ipoint-1)*Atmo%Rho(ipoint-1) + &
          Kappa(ipoint)*Atmo%Rho(ipoint))* &
          (Atmo%z_scale(ipoint-1)-Atmo%z_scale(ipoint))*1.e5

     Atmo%ltau_500(ipoint)=Atmo%ltau_500(ipoint-1)+dtau
!     print *,'ipoint=',ipoint,'kappa,rho=',kappa(ipoint),atmo%rho(ipoint)
!     print *,'   tau0, dtau, tau=',atmo%ltau_500(ipoint)-dtau,dtau,atmo%ltau_500(ipoint)

  End do
!
! H is set to be the same for the first and second points.
!
  If (MinVal(Atmo%ltau_500) .lt. 0) then
      Atmo%ltau_500=Atmo%ltau_500-MinVal(Atmo%ltau_500)+1e-8
   endif
!
! Check that no two points have the same tau (similar to check_tau routine in forward.f90)
!
  Do ipoint=Params%n_points-1, 1, -1
     Do While (Atmo%ltau_500(ipoint) .gt.  Atmo%ltau_500(ipoint+1)-.01*Atmo%ltau_500(ipoint+1))
        Atmo%ltau_500(ipoint+1:Params%n_points)=Atmo%ltau_500(ipoint+1:Params%n_points)+.01*Atmo%ltau_500(ipoint)
     End Do
  End do
!
! Convert to Log(tau)
!
  Atmo%ltau_500=Log10(Atmo%ltau_500)

  Return
!
End Subroutine Z_to_tau
