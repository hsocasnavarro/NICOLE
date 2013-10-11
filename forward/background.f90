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
       num=0.
       den=0.
       Do i=1,n_elements
          num=num+10**(at_abund(i)-12.)*At_weight(i)
          den=den+10**(at_abund(i)-12.)
       End do
       TotH=ph4+phminus4+phplus4+2*ph24+2*ph2plus4
       den=den+Pe4/TotH
       Rho=1.66e-24*TotH/T4/BK*num/den
!       print *,'rho=',rho,num,den
!       print *,'back=',background_opacity,background_opacity/rho
!       print *,'sc=',scat,scat/rho
       scat=scat*6.023e23/num
       Background_opacity=Background_opacity/Rho
       background_opacity=1.e-13
!       Return
    Else if (Opacity_Package .eq. 2) then ! Use Asensio's
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

    if (int(T4) .eq. 4690 .or. int(T4) .eq. 5300) then
       print *,opacity_package,lambda_in4,T4,background_opacity,scat
!       pause
    endif

  End Function Background_opacity

End module background_opacity_module
