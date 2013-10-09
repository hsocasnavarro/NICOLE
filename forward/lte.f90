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

