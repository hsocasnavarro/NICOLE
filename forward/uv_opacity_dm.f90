Module UV_opacity_DM
! UV opacities using data from the paper of Dragon &
! Mutschlecner (1980)
!
! Returns opacity cross-section (cm2) per H nucleus
!
!
  Use UV_opacity_DM_data ! This is where the binary data are stored
  Use Atomic_data, Only: Partition_f, At_ioniz1, At_abund
  Use Phys_constants
  Use LTE, only: Saha
  Implicit None

Contains

  Real Function UVopacity_DM(T, Pe, Lambda, Scat, ignore)
    Real :: T, Pe, Lambda, Scat, Opacity
    Real, Dimension(1) :: TT, U1, U2, U3, dU1, dU2, dU3, n0overn1, Ne
    Real :: Eioniz, Excit, Ab, novernh, Opac, nu, s, A, nul, alpha, alphal
    Integer :: ilev, z, n0
    Real :: Theta, Sum, gff, helium, hydrogen
    Logical, Dimension(92) :: ignore

    TT(1)=T
    Ne(1)=Pe/BK/T
    nu=cc/(Lambda*1e-8)
    
    Opac=0.0
    Do ilev=1, ncontrib ! Loop in levels
       If (lambda .le. DM_lambdal(ilev)) then
          z=iel(ilev)
          Ab=10**(At_abund(z)-12.)
          Call Partition_f(z,TT(1),U1(1),U2(1),U3(1),dU1(1),dU2(1),dU3(1))
          Eioniz=At_ioniz1(z)
          Eioniz=Eioniz*ev_to_cgs
          n0overn1=Saha(1, TT, Ne, U1, U2, Eioniz)
          Excit=DM_Energy(ilev)*ev_to_cgs ! cgs
          novernh=Ab/(1.+1./n0overn1(1)) ! They're all neutral
          novernh=novernh*exp(-Excit/BK/T)/U1(1) ! Boltzmann
          A=DM_A(ilev)
          nul=DM_nul(ilev)*1e15 ! Hz
          s=DM_s(ilev)
          alphal=DM_alphal(ilev)*1e-18 ! cm 2
          alpha=alphal*( A*((nul/nu)**s) + (1.-A)*((nul/nu)**(s+1.)) )
          If (ignore(z)) alpha=0.
          Opac=Opac+novernh*alpha
       End if
    End do

! Add H and He

! Neutral H photoionization
    Call Partition_f(1,TT(1),U1(1),U2(1),U3(1),dU1(1),dU2(1),dU3(1))
    Eioniz=At_ioniz1(1)
    Eioniz=Eioniz*ev_to_cgs
    Theta=1.5777216E5/T
    n0 = 1 + floor(sqrt(1.096776d-3*lambda))
    If (n0 .le. 8) then
       Sum=Exp(Theta/n0**2)/n0**3
       Do ilev=n0+1, 8
          Sum=Sum+Exp(Theta/ilev**2)/ ilev**3
       End do
       Sum=Sum+(0.117+Exp(Theta/81.))/(2.*Theta)
    Else
       Sum=(0.117 + Exp(Theta/n0**2))/Theta
    End if
    gff = (1.-exp(-1.438668E8/(lambda*T) ))*exp(-Theta) * lambda**3
    hydrogen=1.045E-26*gff*Sum
    n0overn1=Saha(1, TT, Ne, U1, U2, Eioniz)
    hydrogen=hydrogen*(1./(1.+1./n0overn1(1))) ! Convert from cm2 per neutral H to cm2 per H particle
    If (ignore(1)) hydrogen=0.
!    
    Helium=0. ! debug
    If (ignore(2)) Helium=0.
!
    Opac=Opac+hydrogen+helium

    Scat=0.
    UVOpacity_DM=Opac
    Return
  End Function UVopacity_DM

End Module UV_OPACITY_DM
