! This routine takes a model in tau (log(tau_5000), T(k), El_p(dyn/cm^2)) 
! and fills in the variables (z(km), rho(g/cm^3), and Gas_p(dyn/cm^2)). Note
! that el_p, Gas_p and rho must be known before calling this routine, as well 
! as nH, nHminus, nHplus and nH2
!
Subroutine Tau_to_z(Params, Atmo)
  Use Param_structure
  Use Model_structure
  Use Background_opacity_module, Only: Background_opacity
  Use Eq_state, Only: Compute_others_from_T_Pe_Pg
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Atmo, Saved_atmo
  Integer :: idx, ipoint, i, j
  Integer, Dimension(1) :: imin
  Logical :: Error
  Real :: Scat, dtau, n2P, dz, metal
  Real, Dimension(Params%n_points) :: Kappa

  Saved_atmo=Atmo
  
!
! Check tau scale for errors
!
  If (MaxVal(Abs(Atmo%ltau_500)) .lt. 1.e-3) then
     If (Params%Printout .ge. 0) then
        Print *,'Error. Something is wrong with the log(tau) scale in the model'
        Print *,'It spans from ',MinVal(Atmo%ltau_500),' to ',MaxVal(Atmo%ltau_500)
        Stop
     End if
  End if
  If (MinVal(Atmo%ltau_500) .gt. -2 .or. MaxVal(Atmo%ltau_500) .lt. 0) then
     If (Params%Printout .ge. 1) then
        Print *,'Warning. There seems to be something weird with the tau scale in the model'
        Print *,'It spans from ',MinVal(Atmo%ltau_500),' to ',MaxVal(Atmo%ltau_500)
        Print *,'Proceeding anyway (hope you know what you are doing). If the results'
        Print *,'are not as expected, this might be the reason'
     End if
  End if
!
! Check el_p, gas_p and rho at z=0
!
  ipoint=-1
  Do idx=1, Params%n_points
     If (Abs(Atmo%ltau_500(idx))-MinVal(Abs(Atmo%ltau_500)) .lt. 1e-5) then
        ipoint=idx
     End if
  End do
  If (ipoint .eq. -1) then 
     Print *,'There seems to be something wrong with the tau-scale. Please, check it'
     Stop
  End if
  If (Atmo%El_p(ipoint) .lt. 1e-2*56.44 .or. Atmo%El_p(ipoint) .gt. 1e2*56.44) then
     If (Params%Printout .ge. 1) then
        Print *,'Warning. Electron pressure near tau=1 is way off typical solar values'
        Print *,'Found (cgs):',Atmo%El_p(ipoint),'  HSRA has:',56.44
        Print *,'Proceeding anyway (hope you know what you are doing). If the results'
        Print *,'are not as expected, this might be the reason'
     End if
  End if
  If (Atmo%Gas_p(ipoint) .lt. 1e-2*1.13E5 .or. Atmo%Gas_p(ipoint) .gt. 1e2*1.31E5) then
     If (Params%Printout .ge. 1) then
        Print *,'Warning. Gas pressure near tau=1 is way off typical solar values'
        Print *,'Found (cgs):',Atmo%Gas_p(ipoint),'  HSRA has:',1.31E5
        Print *,'Proceeding anyway (hope you know what you are doing). If the results'
        Print *,'are not as expected, this might be the reason'
     End if
  End if
  If (Atmo%Rho(ipoint) .lt. 1e-2*3.19E-7 .or. Atmo%Rho(ipoint) .gt. 1e2*3.19E-7) then
     If (Params%Printout .ge. 1) then
        Print *,'Warning. Density near tau=1 is way off typical solar values'
        Print *,'Found (g/cm3):',Atmo%Rho(ipoint),'  HSRA has:',3.19E-7
        Print *,'Proceeding anyway (hope you know what you are doing). If the results'
        Print *,'are not as expected, this might be the reason'
     End if
  End if
!
  Atmo%z_scale(1)=0.
  n2P=BK*Atmo%Temp(1)
  metal=Atmo%Abundance(26)-7.5


  Call Compute_others_from_T_Pe_Pg(Params%n_points, Atmo%Temp,  &
       Atmo%El_p, Atmo%Gas_p, Atmo%nH, Atmo%nHminus, Atmo%nHplus, &
       Atmo%nH2, Atmo%nH2plus)
  
  Kappa(1)=Background_opacity(Atmo%Temp(1), Atmo%El_p(1), Atmo%Gas_p(1), Atmo%nH(1)*n2P, &
       Atmo%nHminus(1)*n2P, Atmo%nHplus(1)*n2P, Atmo%nH2(1)*n2P, &
       Atmo%nH2plus(1)*n2P, 5000., Scat)
  Kappa(1)=Kappa(1)/Atmo%Rho(1)

  Do ipoint=2, Params%n_points
     dtau=10.**Atmo%ltau_500(ipoint) - 10.**Atmo%ltau_500(ipoint-1)
     n2P=BK*Atmo%Temp(ipoint)
     Kappa(ipoint)=Background_opacity(Atmo%Temp(ipoint), Atmo%El_p(ipoint), Atmo%Gas_p(ipoint), &
          Atmo%nH(ipoint)*n2P, Atmo%nHminus(ipoint)*n2P, Atmo%nHplus(ipoint)*n2P, &
          Atmo%nH2(ipoint)*n2P, Atmo%nH2plus(ipoint)*n2P, 5000., Scat)
     Kappa(ipoint)=Kappa(ipoint)/Atmo%Rho(ipoint) ! Convert to cm^2/g
     Atmo%Z_scale(ipoint)=Atmo%Z_scale(ipoint-1) - &
          dtau/2./1.e5* &
          (1./(Kappa(ipoint)*Atmo%Rho(ipoint))+1./(Kappa(ipoint-1)*Atmo%Rho(ipoint-1)))
  End do

  ! Make z=0 at point with min(abs(ltau_5000))
  imin=MinLoc(Abs(Atmo%ltau_500))
  Atmo%Z_scale(1:Params%n_points)=Atmo%Z_scale(1:Params%n_points) - &
       Atmo%Z_scale(imin(1))
! Refine with linear interpolation. Make z=0 at ltau_5000=0
  i=imin(1)
  j=imin(1)+1
  if (j .gt. Params%n_points) then
     j=Params%n_points
     i=i-1
     If (i .lt. 1) then
        Print *,'Error in z scale'
        Stop
     End if
  End if
  dz=-Atmo%ltau_500(i)*(Atmo%z_scale(j)-Atmo%z_scale(i))/(Atmo%ltau_500(j)-Atmo%ltau_500(i))
  Atmo%Z_scale(:)=Atmo%Z_scale(:)+dz
!
  Return
!
End Subroutine Tau_to_z
