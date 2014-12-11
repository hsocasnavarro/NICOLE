module background_opacity_module
  Use Phys_constants
  Use Wittmann_opac_module
  Use Profiling
  Use Asensio_background_opacity_module, Only: Asensio_background_opacity
  Use UV_opacity_TOPbase, Only: UVopacity_TOPbase
  Use UV_opacity_DM, Only: UVopacity_DM
  Use Eq_state
  Implicit None
!  Integer, Parameter :: maxnfudge=1000
  Integer :: Opacity_package=1, Opacity_Package_UV=1 ! Default, will be overwritten in main program
  Logical, Dimension(92) :: elneglectopac
  Integer :: nfudge=0 ! These are defined in forward.f90
  Real :: Fudge_start(1000), Fudge_end(1000), Fudge(1000)

Contains
  Function Background_opacity(T4, Pe4, Pg4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4, Scat)
    Real :: T4, Pe4, Pg4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4    
    Real :: Scat, Scat2, Background_opacity, TotH, llambda_in4, Opac
    Real :: nu, chi_0, chi_e, eta, num, den, Rho
    Integer :: i

    Call Time_routine('background_opac',.True.)

    TotH=ph4+phminus4+phplus4+2*ph24+2*ph2plus4
    If (Opacity_Package .eq. 1) then ! Use Wittmann'
       Opac=Wittmann_opac(T4, Pe4, Pg4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4, Scat, elneglectopac)
       Background_opacity=Opac
       ! If lambda .le. 4000 use specific module for UV
       if (lambda_in4 .le. 4000) then
          If (Opacity_Package_UV .eq. 1) then ! TOPbase
             Background_opacity=Opac+ &
                  UVopacity_TOPbase(T4, Pe4, lambda_in4, Scat2, elneglectopac)*TotH/BK/T4
             Scat=Scat+Scat2*TotH/BK/T4
          Else If (Opacity_Package_UV .eq. 2) then ! Dragon-Mutschlecner (1980)
             Background_opacity=Opac+ &
                  UVopacity_DM(T4, Pe4, lambda_in4, Scat2, elneglectopac)*TotH/BK/T4
             Scat=Scat+Scat2*TotH/BK/T4
          Else
             Print *,'Unknow UV Opacity Package'
             Stop
          End if
       End if

    Else if (Opacity_Package .eq. 2) then ! Use Asensio's
       Opac=Asensio_background_opacity(T4, Pe4, Pg4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4, Scat, elneglectopac)
       Background_opacity=Opac
       ! If lambda .le. 4000 use specific module for UV
       if (lambda_in4 .le. 4000) then
          If (Opacity_Package_UV .eq. 1) then ! TOPbase
             Background_opacity=Opac+ &
                  UVopacity_TOPbase(T4, Pe4, lambda_in4, Scat2, elneglectopac)*TotH/BK/T4
             Scat=Scat+Scat2*TotH/BK/T4
          Else If (Opacity_Package_UV .eq. 2) then ! Dragon-Mutschlecner (1980)
             Background_opacity=Opac+ &
                  UVopacity_DM(T4, Pe4, lambda_in4, Scat2, elneglectopac)*TotH/BK/T4
             Scat=Scat+Scat2*TotH/BK/T4
          Else
             Print *,'Unknow UV Opacity Package'
             Stop
          End if
       End if
       
    Else if (Opacity_Package .eq. 3) then ! Use SOPAS
       TotH=ph4+phminus4+phplus4+2*ph24+2*ph2plus4
       nu=cc/(lambda_in4*1e-8) ! s^-1
       Call Sopas(1, 2, nu, TotH, Pe4, T4, Pg4, chi_0, chi_e, eta, elneglectopac)
       Background_opacity=chi_0+chi_e
       Scat=chi_e
    endif

! Apply fudge factors?
    Do i=1, nfudge
       If (lambda_in4 .ge. Fudge_start(i) .and. lambda_in4 .le. Fudge_end(i)) then
          Background_opacity=Background_opacity*Fudge(i)
          Scat=Scat*Fudge(i)
       End if
    End do

    Call Time_routine('background_opac',.False.)
    Return
  End Function Background_opacity

End module background_opacity_module
