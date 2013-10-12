module background_opacity_module
  Use Phys_constants
  Use Wittmann_opac_module
  Use Asensio_background_opacity_module, Only: Asensio_background_opacity
  Implicit None
  Integer :: Opacity_package=1 ! Default, will be overwritten in main program

Contains
  Function Background_opacity(T4, Pe4, Pg4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4, Scat)
    Real :: T4, Pe4, Pg4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4    
    Real :: Scat, Background_opacity, TotH
    Real :: nu, chi_0, chi_e, eta, num, den, Rho
    Integer :: i

    If (Opacity_Package .eq. 1) then ! Use Wittmann'
       Background_opacity=Wittmann_opac(T4, Pe4, Pg4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4, Scat)
    Else if (Opacity_Package .eq. 2) then ! Use Asensio's
       TotH=ph4+phminus4+phplus4+2*ph24+2*ph2plus4
       Background_opacity=Asensio_background_opacity(T4, Pe4, Pg4, PH4, PHminus4, PHplus4, PH24, PH2plus4, lambda_in4, Scat)
!       Return
    Else if (Opacity_Package .eq. 3) then ! Use SOPAS
       TotH=ph4+phminus4+phplus4+2*ph24+2*ph2plus4
       nu=cc/(lambda_in4*1e-8) ! s^-1
       Call Sopas(1, 2, nu, TotH, Pe4, T4, Pg4, chi_0, chi_e, eta)
       Background_opacity=chi_0+chi_e
       Scat=chi_e
!       Return
    endif
  End Function Background_opacity

End module background_opacity_module
