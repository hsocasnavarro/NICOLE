Module UV_opacity_TOPbase
! UV opacities using data from The Opacity Project (TOPbase)
! Only neutrals and first ion, photoionization processes
! Considers most elements up to Z=26 (Fe), many levels
! Valid only in the wavelength interval from 500A to 4000A
!
! Returns opacity cross-section (cm2) per H nucleus
!
! (Hector Socas-Navarro thanks Carlos Allende Prieto, Manuel Bautista
!   and Sultana Nahar for nicely formatted FeI and FeII data)
!
  Use UV_opacity_TOPbase_data ! This is where the binary data are stored
  Use Atomic_data, Only: Partition_f, At_ioniz1, At_abund
  Use Phys_constants
  Use LTE, only: Saha
  Implicit None

Contains

  Real Function UVopacity_TOPbase(T, Pe, Lambda, Scat)
    Real, Parameter :: megabarn_to_cm2=1.e-18
    Real :: T, Pe, Lambda, Scat, Opacity
    Real, Dimension(1) :: TT, U1, U2, U3, dU1, dU2, dU3, n0overn1, Ne
    Real :: Eioniz, Excit, Ab, novernh, Opac
    Integer :: ilam, iz, ilev
    Logical :: FirstTime=.True.

    If (FirstTime) then ! Read data from file the first time
       Call UVopacity_TOPbase_init
       FirstTime=.False.
    End if

    TT(1)=T
    Ne(1)=Pe/BK/T
    ilam=Lambda-500.+.5
    
    If (ilam .lt. 1) ilam=1
    If (ilam .gt. 3500) ilam=3500

    Opac=0.0
    Do iz=1,nz ! Loop in elements
       If (ilam .le. Maxval(nlambda(iz,:))) then ! Is at least one level relevant?
          Ab=10**(At_abund(zs(iz))-12.)
          Call Partition_f(zs(iz),TT(1),U1(1),U2(1),U3(1),dU1(1),dU2(1),dU3(1))
          Eioniz=At_ioniz1(zs(iz))
          Eioniz=Eioniz*ev_to_cgs
          n0overn1=Saha(1, TT, Ne, U1, U2, Eioniz)
          Do ilev=1, nlevels(iz) ! Levels to consider
             If (ilam .le. nlambda(iz,ilev)) then ! Is in the wavelength range?
                Excit=Energy(iz,ilev)
                If (ionstage(iz,ilev) .eq. 1) then ! It's a neutral
                   novernh=Ab/(1.+1./n0overn1(1))
                Else ! It's singly ionized
                   novernh=Ab/(1.+n0overn1(1))
                End if
                novernh=novernh*exp(-Excit/BK/T)/U1(1) ! Boltzmann
                Opac=Opac+novernh*XSec(iz,ilev,ilam)*megabarn_to_cm2
             End if
          End do
       End if
    End do

    UVOpacity_TOPbase=Opac
    Scat=0.
    Return
  End Function UVopacity_TOPbase

End Module UV_opacity_TOPbase
