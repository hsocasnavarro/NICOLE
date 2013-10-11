Module LTE
  Use Param_structure
  Use Model_structure
  Use Line_data_structure
  Use Atomic_data

Contains


Function Saha(npoints, T, Ne, U1, U2, E_ioniz)
!
! The Saha equation. This function returns the ratio of the atoms in the
! ionization stage 1 to those in the next ionization stage (e.g., neutrals
! over ions)
!
  Implicit none
  Integer :: npoints
  Real, dimension (npoints) :: T, Ne, U1, U2, Saha
  Double Precision, dimension (npoints) :: Sahad
  Real :: E_ioniz
!
  Sahad(1:npoints)=Ne(1:npoints)*2.07D-16*U1(1:npoints)/U2(1:npoints)* &
       (T(1:npoints)**(-1.5))*exp(E_ioniz/bk/T(1:npoints)*1.d0)
  Where(Sahad .gt. 1.e30) 
     Sahad=1.e30
  End Where

  Saha=Sahad
  Return
End Function Saha

Subroutine Saha123(npoints,iel,T,Ne,n0overn,n1overn,n2overn)
! Solve (using safety checks) the ionization equilibrium for the three lowest ionization
! stages of atomic species of atomic number iel
! Exception: If iel .eq. 1 (H atom), then the third ionization stage actually refers to
! the H- ion
! 
  Use Atomic_data
  Implicit none
  Integer :: iel, ip, npoints
  Real :: Ioniz
  Real, Parameter :: BK= 1.38066D-16, eV_to_cgs=1.602189D-12
  Real, dimension (npoints) :: T, Ne, U1, U2, U3, D1, D2, D3, n1overn, n2overn,  &
       n0overn,n3overn,nminusovern0,n1overn0,n2overn1
!
  Do ip=1, npoints
     Call Partition_f(iel, T(ip), U1(ip), U2(ip), U3(ip), D1(ip), D2(ip), D3(ip)) 
  End do
  Ioniz=At_ioniz1(iel)*eV_to_cgs
  n1overn0=1./saha(npoints, T, Ne, U1, U2, Ioniz)
  If (iel .ge. 2) then ! generic atom (say Fe)
     Ioniz=At_ioniz2(iel)*eV_to_cgs
     n2overn1=1./saha(npoints, T, Ne, U2, U3, Ioniz)
  Else ! H-
     U2(:)=1.0
     Ioniz=0.754*eV_to_cgs ! ionization energy for Hminus
     nminusovern0=saha(npoints, T, Ne, U2, U1, Ioniz) ! For Hminus
  End if
  Do ip=1, npoints
     If (iel .eq. 1) then ! H, H+ and H-
        If (nminusovern0(ip) .gt. 1e20) then ! It's all H-
           n2overn(ip)=1.
           n1overn(ip)=0.
           n0overn(ip)=0.
        Else if (n1overn0(ip) .gt. 1e20) then ! H is negligible (and so is H-), it's all H+
           n1overn(ip)=1.
           n2overn(ip)=0.
           n0overn(ip)=0.
        Else ! Consider them all
           n0overn(ip)=1./(1+nminusovern0(ip)+n1overn0(ip))
           n1overn(ip)=n1overn0(ip)*n0overn(ip)
           n2overn(ip)=nminusovern0(ip)*n0overn(ip)
        End if
     Else ! generic atom (say Fe)
        If (n2overn1(ip) .gt. 1e20) then ! It's all Fe++
           n2overn(ip)=1.
           n1overn(ip)=0.
           n0overn(ip)=0.
        Else if (n1overn0(ip) .gt. 1e20) then ! Fe is negligible, it's all Fe+ and Fe++
           n0overn(ip)=0.
           n2overn(ip)=1./(1+n2overn1(ip))
           n1overn(ip)=1.-n2overn(ip)
        Else ! Consider them all
           n0overn(ip)=1./(1.+n2overn1(ip)*n1overn0(ip)+n1overn0(ip))
           n1overn(ip)=n0overn(ip)*n1overn0(ip)
           n2overn(ip)=n1overn(ip)*n2overn1(ip)
        End if
     End if
     print *,'t,h=',t(60),n0overn(60),n1overn(60),n2overn(60)
     pause
  End do

End Subroutine Saha123

!****************************************************
!
      FUNCTION PLANCK(U,T)

!
!  CALCULATES PLANCK FUNCTION BNY AT FREQUENCY U, TEMP T
!
      DATA  EM/9.109534E-28/,UU/1.6605655E-24/
!
      X=HH*U/BK/T
!      PRINT *,'BK=',BK
      
!      PRINT *,'T=',T
!      PRINT *,'U=',U
!      PRINT *,'X=',X
      
      IF(X.LT.80.) THEN
        PLANCK=2.0*HH*U/CC*U/CC*U/(EXP(X)-1.0)
      ELSE
        PLANCK=2.0*HH*U/CC*U/CC*U*EXP(-X)
      ENDIF
!
      RETURN
    END FUNCTION PLANCK


Subroutine LTE_pop(Params, Line, Atmo, ng_i, ng_j)
!
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Atmo
  Type (Line_data) :: Line
  Real, dimension (Params%n_points) :: N_e, U0, U1, U2, U_ion
  Real, dimension (Params%n_points) :: N_ion, ng_i, ng_j, N0_N1, N1_N2, N_tot
  Integer :: index, npoints, idepth, iline, Atomic_number
  Real :: du0, du1, du2, E_ioniz1, E_ioniz2, E_exc_i, E_exc_j, Abund_mass
  Real :: nu
  Logical, Save :: Warning=.FALSE.
!
  npoints=Params%n_points
  Atomic_number=Line%Atomic_number
  Abund_mass=At_weight(Atomic_number)*(10.**(At_abund(Atomic_number)-12.))/ &
       Sum(At_weight(1:N_elements)*(10.**(At_abund(1:N_elements)-12.)))
  N_tot(1:npoints)=Abund_mass*Atmo%rho(1:npoints)/At_weight(Atomic_number)/mass_pr

  N_e(1:npoints)=Atmo%El_p(1:npoints)/bk/Atmo%Temp(1:npoints)*1.d0
  Do idepth=1, npoints
     Call Partition_f(Atomic_number, Atmo%Temp(idepth), &
          U0(idepth), U1(idepth), U2(idepth), du0, du1, du2)
  End do
  E_ioniz1=At_ioniz1(Atomic_number)*eV_to_cgs
  E_ioniz2=At_ioniz2(Atomic_number)*eV_to_cgs
!
! Saha ionization equations
!
  N0_N1(1:npoints)=Saha(npoints, Atmo%Temp, N_e, U0, U1, E_ioniz1)
  N1_N2(1:npoints)=Saha(npoints, Atmo%Temp, N_e, U1, U2, E_ioniz2)
  N_ion(1:npoints)=N_tot(1:npoints)/  &
       (1.+N0_N1(1:npoints)+1./N1_N2(1:npoints)) ! This is N_1
  U_ion(1:npoints)=U1(1:npoints)
  If (Line%Ion_stage .eq. 1) then ! N_0
     N_ion(1:npoints)=N_ion(1:npoints)*N0_N1(1:npoints)
     U_ion(1:npoints)=U0(1:npoints)
  Else If (Line%Ion_stage .ge. 3) then ! N_2
     N_ion(1:npoints)=N_ion(1:npoints)/N1_N2(1:npoints)
     U_ion(1:npoints)=U2(1:npoints)
  End if
  If (Line%Ion_stage .gt. 3 .and. .not. Warning) then
     Warning=.TRUE.
     If (Params%Printout .ge. 1) Print *,'Warning. Ionization stage gt 3. Assuming 3 in LTE pop'
  End if

!
! Boltzmann excitation equations
!
  nu=cc/(Line%Wlength*1e-8) ! s^-1
  E_exc_i=Line%Energy_low*eV_to_cgs
  E_exc_j=E_exc_i+hh*nu
  ng_i(1:npoints)=N_ion(1:npoints)*exp(-E_exc_i/bk/Atmo%Temp(1:npoints))/ &
       U_ion(1:npoints)
  ng_j(1:npoints)=N_ion(1:npoints)*exp(-E_exc_j/bk/Atmo%Temp(1:npoints))/ &
       U_ion(1:npoints)
  Return
End Subroutine LTE_pop
!

End Module LTE

