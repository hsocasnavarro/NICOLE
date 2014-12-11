! This module solves the forward problem, i.e., the synthesis of the
! spectral lines assuming LTE.
!
!
Module Forward_module
  Use Forward_support
  Use Param_structure
  Use Model_structure
  Use Eq_state
  Use Zeeman_splitting
  Use LTE
  Use NLTE_module
  Use Atomic_data
  Use Gauss_quad
  Use Line_data_structure
  Use Bezier_math
  Use Bezier_solvers
  Use Background_opacity_module
  Use nr
  Use File_operations
  Use Profiling
  Use Debug_module
  Character (len=25) :: Forward_ver='NICOLE Forward v3.6'
!
Contains
! CONHSRA calculates the continuum in the HSRA model by means of
! a polynomial fit. The maximum error is .125%	

real function conhsra(x)

  implicit real (a-h,o-z)
  real c1(11),c2(11),c3(11),c4(11),c5(11),c6(11),c7(11)
  data c1/-4.906054765549e13,1.684734544039e11,1.507254517567e7 &
       ,-7561.242976546,0.,0.,0.,0.,0.,0.,0./
  data c2/-4.4650822755e14,6.1319780351059e11,-9.350928003805e7 &
       ,0.,0.,0.,0.,0.,0.,0.,0./
  data c3/-1.025961e15,1.3172859e12,-3.873465e8 &
       ,46486.541,-2.049,0.,0.,0.,0.,0.,0./
  data c4/4.861821e15,-2.2589885e12,4.3764376e8 &
       ,-39279.61444,1.34388,0.,0.,0.,0.,0.,0./
!	data c5/4.353082e17,-1.261996e14,5.31125e9,908890.899,43.24361 &
  !             7,-2.0462285d-2,1.0018226758d-6,0.,0.,0.,0./
  data c5/1.758394e15,-3.293986e11,1.6782617e7,0.,0.,0. &
       ,0.,0.,0.,0.,0./
  data c6/1.61455557e16,-6.544209e12,1.0159316e9 &
       ,-70695.58136,1.852022,0.,0.,0.,0.,0.,0./
  data c7/7.97805136e14,-1.16906597e11,5.315222e6 &
       ,-4.57327954,-3.473452d-3,0.,0.,0.,0.,0.,0./

  if(x.lt.3644.15)then
     conhsra=c1(1)+x*(c1(2)+x*(c1(3)+x*(c1(4)+x*c1(5))))
  else if(x.lt.3750.)then
     conhsra=c2(1)+x*(c2(2)+x*(c2(3)+x*(c2(4)+x*c2(5))))
  else if(x.lt.6250.)then
     conhsra=c3(1)+x*(c3(2)+x*(c3(3)+x*(c3(4)+x*c3(5))))
  else if(x.lt.8300.)then
     conhsra=c4(1)+x*(c4(2)+x*(c4(3)+x*(c4(4)+x*c4(5))))
  else if(x.lt.8850.)then
     conhsra=c5(1)+x*(c5(2)+x*(c5(3)+x*(c5(4)+x*c5(5))))
  else if(x.lt.10000.)then
     conhsra=c6(1)+x*(c6(2)+x*(c6(3)+x*(c6(4)+x*c6(5))))
  else
     conhsra=c7(1)+x*(c7(2)+x*(c7(3)+x*(c7(4)+x*c7(5))))
  end if

  if (x .gt. 18000.) then
     print *,'Sorry. conhsra doesnt work over 1.8 microns'
     print *,'Contact the author to fix this or use a different'
     print *,'normalization scheme'
     Stop
  End if

  return
end function conhsra

! This routine adds a fraction of stray light a given synthetic profile.
!
Subroutine Add_stray_light(Params, Atmo, Region, Profile)
  Implicit None
  Type (Parameters) :: Params
  Type (model) :: Atmo
  Type (Region_data), dimension (Params%n_regions) :: Region
  Real, Dimension (Params%n_data) :: Profile
  Integer :: ireg, iwave, idata
  Real :: ff
!
  idata=1
  Do ireg=1, Params%n_regions
     ff=1.-Atmo%stray
     Do iwave=1, Region(ireg)%nwavelengths
        Profile(idata:idata+3) = ff*Profile(idata:idata+3) + &
             (1.-ff)*Params%Stray_prof(idata:idata+3)
        idata=idata+4
     End do
  End do
!
  Return
End Subroutine Add_stray_light
! This routine convolves the profile with a Gaussian whose width is the
! macroturbulent velocity.
!
Subroutine Convolve_profile(Params, Region, Atmo, Profile)
  Implicit None
  Type (Parameters) :: Params
  Type (Region_data), dimension (Params%n_regions) :: Region
  Type (model) :: Atmo
  Real, Dimension (Params%n_data) :: Profile, gauss, data
  Real, Dimension (Params%n_data*2) :: ans, respns
  Integer :: i, n, m, mmax, n2, iregion, istokes, iwave, istart, iunit, indreg1
  Real :: x, Wave, Maximum
  Logical, Save :: FirstTime=.True., ExistsInstrumProf=.False.
  Real, Dimension(:,:), Pointer, Save :: Instrum_prof
!  Integer, Dimension(:), Pointer, Save :: mm
!  Logical :: Exists
!
!  If (FirstTime) then 
!     FirstTime=.False.
     ! Read instrumental profile from file?
!     Inquire (File='Instrumental_profile.dat', Exist=Exists)
!     If (Params%IProfExists) then
  ExistsInstrumProf=.False.
  If (MaxVal(Abs(Params%IProf)) .gt. 0.0001) ExistsInstrumProf=.True.
!        If (Params%Printout .ge. 1) &
!           Print *,'Reading Instrumental_profile.dat'
!        Allocate(mm(Params%n_regions))
!        ! File format is: 
!        !  For each region:
!        !  First line -> number of wavelength points. It must be odd and
!        !     cannot be greater than the number of wavelengths in the grid 
!        !     file
!        !  Each line -> Instrument profile value, centered at first point
!        !     then from m/2+2 to m is the left part of the profile
!        Call Open_file(iunit,'Instrumental_profile.dat')
!        mmax=0
!        Do n=1, Params%n_regions ! Find region with most wavelengths
!           Read (iunit,*) m
!           mm(n)=m
!           If (m .gt. mmax) mmax=m
!           Do i=1, m
!              Read (iunit,*) x
!           End do
!        End do
!        Allocate(Instrum_prof(Params%n_regions,mmax))
!        Close(iunit)
!        Call Open_file(iunit,'Instrumental_profile.dat')
!        Do n=1, Params%n_regions
!           Read (iunit,*) m
!           Do i=1, m
!              Read (iunit,*) Instrum_prof(n, i)
!           End do
!           If (Mod(m,2) .ne. 1) then ! If even number of points, make it odd
!              mm(n)=mm(n)-1
!              Instrum_prof(n,m)=0.
!           End if
!        End do
!        Close (iunit)
!     End if
!  End if
!
  If (Atmo%v_mac .lt. 1.e-6 .and. .not. ExistsInstrumProf) Return
!
  istart=1 ! Index to the first point of the line
  indreg1=1
  Do iregion=1, Params%n_regions
     If (ExistsInstrumProf) then
        ! Use instrumental profile from file
        m=Params%mm(iregion)
        if (m/2 .ne. m/2.) then
           Print *,'Error. Instrumental profile for region ',iregion
           Print *,'should be an odd number. It is ',m
           Stop
        End if
        gauss(1:m)=Params%IProf(indreg1:indreg1+m-1)
        indreg1=indreg1+m
     Else
        ! Construct the gaussian
        Wave=Region(iregion)%First_wlength + &
             Region(iregion)%nwavelengths/2.*Region(iregion)%Wave_step
        m=5.*(Atmo%v_mac/cc*Wave/Region(iregion)%Wave_step*Region(iregion)%Macro_enh)
        m=(m/2.)
        m=m*2+1 ! Make sure m is an odd integer
        gauss(1)=1. ! To prevent the case where v_mac .eq. 0.
        If (m .gt. 1) then
           Do i=1, m/2+1 ! Positive x
              x=(i-1)*Region(iregion)%Wave_step/ &
                   (Atmo%v_mac*Region(iregion)%Macro_enh/cc*Wave)
              gauss(i)=Exp(-x*x)
           End do
        End if
        Do i=m/2+2, m ! Negative x
           gauss(i)=gauss(m-i+1)
        End do
     End if
     gauss(1:m)=gauss(1:m)/Sum(gauss(1:m)) ! Normalization
!
! Construct the profile to convolve (one for each Stokes profiles)
!
     n=Region(iregion)%nwavelengths ! Number of points
     If (m .gt. n) then
        Print *,'Error in Convolve_profile. m .gt. n !!!'
        Return
     Else
! Find nearest power of 2
        i=1
        Do while (2.**i .lt. (n+m))
           i=i+1
        End do
        n2=2.**i
        If (n2 .gt. Params%n_data) then
           Print *,'Number of wlengths=',Params%n_data/4,'. n2=',n2
           Print *,'Cannot convolve'
           Return
        End if
! 
        Do istokes=1, 4
           Do iwave=1, n
              data(iwave)=Profile(istart+(iwave-1)*4+istokes-1)
           End do
           Maximum=MaxVal(Abs(data(1:n)))
           If (Maximum .lt. 1e-5) Maximum=1.

           data(n+1:n+1+m/2)=data(n)
           data(n+1+m/2+1:n2)=data(1)
           data(1:n2)=data(1:n2)/Maximum
!
! Perform the convolution
!
           respns(:)=0.
           respns(1:m)=gauss(1:m)
           ans(1:n2)=convlv(data(1:n2), respns(1:m), 1)
           ans(1:n2)=ans(1:n2)*Maximum
!
! Store the result back in Profile
!
           Do iwave=1, n
              Profile(istart+(iwave-1)*4+istokes-1)=ans(iwave)
           End do
        End do ! Next Stokes parameter
     End if
     istart=istart + 4*n ! Update istart to point to the next region
  End do ! Next spectral region
!
! Done!
!
  Return
End Subroutine Convolve_profile

! This subroutine takes one of the density variables (electron pressure, 
! electron density, gas pressure or gas density) in c.g.s units and fills
! the other two using the chemical equilibrium routines. If electron density
! is used, it is first converted to electron pressure. The input parameter 
! input_dens must be set to either (lower case) 'pel','nel','pgas' or 'dens' 
! to specify which one to take.
!
! This routine also reverses the model if necessary to make sure that
! tau is monotonically increasing (i.e., atmosphere starts at the top)

Subroutine Fill_densities(Params, Input_dens, Atmo)
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Atmo, Atmo_in
  Character (len = 4) :: Input_dens
  Real, Dimension (10) :: Pp
  Integer, Parameter :: niters=50
  Real, Parameter :: Precision = 1e-5 ! 0.01
  Integer :: idepth, iters, iatom
  Real :: Pel, Pg, Wsum, Asum, Avmolweight, dif, temp
  Logical :: Warning
  Real, Parameter :: Min_Pe=1e-4
!
  Call Time_routine('Fill densities',.true.)
!
!
  Atmo_in=Atmo
!
! Calculate the average molecular weight of the particles
!
  Wsum=0.
  Asum=0.
  Do iatom=1, 92
     Wsum=Wsum+At_weight(iatom)*10**(At_abund(iatom)-12.)
     Asum=Asum+10**(At_abund(iatom)-12.)
  End do
!
  If (Input_dens .eq. 'nel' .or. Input_dens .eq. 'pel') then
!     el_p has already been computed by the Python script run_nicole.py
     Atmo%ne=Atmo%el_p/bk/Atmo%temp
     Call Compute_Pg(Params%n_points, Atmo%Temp, Atmo%El_p, Atmo%Gas_p)
     Do idepth=1, Params%n_points
        Avmolweight=Wsum/(Asum+ &
             Atmo%El_p(idepth)/Atmo%Gas_p(idepth))
!             Atmo%ne(idepth)/(Atmo%nH(idepth)+Atmo%nHplus(idepth)+ &
!             Atmo%nHminus(idepth)+Atmo%nH2(idepth))) ! Average molecular weight.
        Atmo%Rho(idepth)=Atmo%Gas_p(idepth)*Avmolweight/Avog/bk/Atmo%Temp(idepth)
     End do
  Else if (Input_dens .eq. 'pgas') then
     Call Compute_Pe(Params%n_points, Atmo%Temp, Atmo%Gas_p, Atmo%El_p)
     Atmo%ne=Atmo%el_p/bk/Atmo%temp
     Do idepth=1, Params%n_points
        Avmolweight=Wsum/(Asum+ &
             Atmo%El_p(idepth)/Atmo%Gas_p(idepth))
!             Atmo%ne(idepth)/(Atmo%nH(idepth)+Atmo%nHplus(idepth)+ &
!             Atmo%nHminus(idepth)+Atmo%nH2(idepth))) ! Average molecular weight.
        Atmo%Rho(idepth)=Atmo%Gas_p(idepth)*Avmolweight/Avog/bk/Atmo%Temp(idepth)
     End do
  Else if (Input_dens .eq. 'dens') then
     Do idepth=1, Params%n_points ! Initial guess to Gas_P without electrons
        Avmolweight=Wsum/Asum ! Average molecular weight.
        Atmo%Gas_p(idepth)=Atmo%Rho(idepth)/Avmolweight*Avog*bk*Atmo%Temp(idepth)
     End do
     Do iters=1,3
        Call Compute_Pe(Params%n_points, Atmo%Temp, Atmo%Gas_p, Atmo%El_p)
        Atmo%ne=Atmo%el_p/bk/Atmo%temp
        Do idepth=1, Params%n_points
           Avmolweight=Wsum/(Asum+ &
                Atmo%El_p(idepth)/Atmo%Gas_p(idepth))
   !             Atmo%ne(idepth)/(Atmo%nH(idepth)+Atmo%nHplus(idepth)+ &
   !             Atmo%nHminus(idepth)+Atmo%nH2(idepth))) ! Average molecular weight.
           Atmo%Gas_p(idepth)=Atmo%Rho(idepth)/Avmolweight*Avog*bk*Atmo%Temp(idepth)
        End do
     End do
  Else
     Print *,'Error: Input density must be one of the following: Pel, Nel, Dens or Pgas'
     Print *,'File NICOLE.input. Found:',Input_dens
     Stop
  End if
! Check atmosphere for sanity
  Do idepth=1, Params%n_points
     If (Atmo%El_p(idepth) .lt. Min_Pe) then
        Call Debug_log('In fill_densities. Pe .lt. Min_Pe parameter. Clipping it',2)
        Atmo%El_p(idepth)=Min_Pe
     End if
  End do
!
  Call Reset_densities(Atmo, Atmo_in) ! Reset densities in keep switches
!

!
  Call Time_routine('Fill densities',.false.)
  Return

End Subroutine Fill_densities
! This routine calls a suitable formal solution routine to compute the
! emergent Stokes vector.
! If formal_solut == 0 , it is chosen automatically (Bezier or WPM)
!                    == 1 , it uses cubic Bezier splines
!                    == 2 , it uses DELO Bezier
!                    == 3 , it uses Bezier scalar (Q=U=V=0)
!                    == 4 , it uses Hermite
!                    == 5 , it uses WPM
!                    == 6 , it uses DELO
!                    == 7 , it uses DELOPAR (DELO with parabolic interpolat)
!                    == 8 , it uses Short Characteristics (Q=U=V=0)
!
! By Hermite I mean the Hermitian solution method of Bellot Rubio, Ruiz
!   Cobo and Collados 1998, ApJ, 506, 805
! By WPM I mean the Weakly Polarizing Media approximation of Sanchez
!   Almeida and Trujillo Bueno 1999, ApJ, in press.
!
Subroutine Formal_solution(npoints, formal_solut, ltau_500, Absorp_in, &
     Source_f, Stokes, ichoice)
  Implicit None
  Integer :: formal_choice, ipoint, npoints, formal_solut
  Real, Dimension (npoints) :: ltau_500, Source_f, dtau_nu, tau_500
  Real, Dimension (npoints,4,4) :: Absorp_in, Absorp
  Real, Dimension (4) :: Stokes
  Logical :: suitable, CheckNaN
  Real, Parameter :: Threshold = 10.
  Integer :: ichoice, itry, k, i, j
!
  itry=1
  formal_choice=formal_solut
  If (formal_choice .eq. 0) then ! Choose one
!    See if WPM is suitable for this case
     suitable=.TRUE.
     ipoint=1
     Do while (ipoint .le. npoints .and. suitable)
        suitable = Absorp(ipoint,1,1) .gt. Threshold*( &
             Sum(Abs(Absorp(ipoint,1,2:4))))
        ipoint=ipoint+1
     End Do
     If (suitable) then ! WPM
        formal_choice=5 !  WPM
     else ! Bezier
        formal_choice=1 ! Bezier 
     End if
  End if
!
! Compute monochromatic delta tau scale
!
  tau_500=10.d0**ltau_500
  dtau_nu(1)=1e-10
  call bezier_qintep(npoints, tau_500, Absorp_in(:,1,1), dtau_nu)
!
  Stokes(:)=0.
  Do While (itry .eq. 1 .or. (CheckNaN(Stokes(1)) .and. itry .eq. 2) )
!
     Absorp=Absorp_in
     If (formal_choice .eq. 1) then
        Call delobezier3(npoints, dtau_nu, Absorp, Source_f, Stokes)
        ichoice=1
     Else If (formal_choice .eq. 2) then
        Call delobezier(npoints, dtau_nu, Absorp, Source_f, Stokes)
        ichoice=2
     Else If (formal_choice .eq. 3) then
        Call delobezier_scal(npoints, dtau_nu, Absorp, Source_f, Stokes)
        ichoice=3
     Else If (formal_choice .eq. 4) then
        Call Hermite(npoints, dtau_nu, Absorp, Source_f, Stokes)
        !     Call myold_Hermite(npoints, ltau_500, Absorp, Source_f, Stokes)
        ichoice=4
     Else if (formal_choice .eq. 5) then
        Call WPM(npoints, dtau_nu, Absorp, Source_f, Stokes)
        ichoice=5
     Else if (formal_choice .eq. 6) then
        Call Delolin(npoints, dtau_nu, Absorp, Source_f, Stokes)
        ichoice=6
     Else if (formal_choice .eq. 7) then
        Call Delopar(npoints, dtau_nu, Absorp, Source_f, Stokes)
        ichoice=7
     Else if (formal_choice .eq. 8) then
        Call SC_solver(npoints, dtau_nu, Absorp, Source_f, Stokes)
        ichoice=8
     End if
     !
     If (CheckNaN(Stokes(1))) then ! Something went wrong
        print *,'something went wrong',stokes
        If (itry .eq. 2) then ! Bad result, couldn't fix it
           print *,'cant fix it'
           Stokes(:)=1e10
           Call Debug_log('NaNs in formal solution',1)
        Else ! Let's try to fix it
           print *,'trying to fix it...'
           Do k=1, npoints
              If (dtau_nu(k) .lt. 1e-8) then
!                 print *,'k,dtau_nu=',k,dtau_nu(k)
                 dtau_nu(k)=1e-8
              endif
              If (Absorp_in(k,1,1) .lt. 0) then
!                 print *,'k,i,j,absorp=',k,i,j,absorp(k,i,j)
                 Absorp_in(k,1,1)=1.e-12
              End if
           End do
        End if
     End if
     itry=itry+1
  End do
End Subroutine Formal_solution

!
!**********************************************************************
! This subroutine computes the Zeeman-split profiles
!
Subroutine Profiles(Params, Line, nwlengths, nw, IndexWave, Wave, Atmo, Damp, Dldop, &
             Line_op, Phi_I, Phi_Q, Phi_U, Phi_V, Psi_Q, Psi_U, Psi_V, &
             Cont_op_5000, Incomplete)
  Use Param_structure
  Use Model_structure
  Use Line_data_structure
  Use Profiling
  Implicit None
  Type (Parameters) :: Params
  Type (Line_data) :: Line
  Type (Model) :: Atmo
  Integer :: ind, npoints, nr, nl, np, nwlengths, nw
  Real, Parameter :: ZeroFieldThreshold=1.E-6
  Real, Dimension(nwlengths) :: Wave, DWave
  Real, Dimension(nw) :: vv
  Integer, Dimension(nwlengths) :: IndexWave
  Integer, Dimension(2) :: ji, jf, mult
  Character, Dimension(2) :: desig
  Real, Dimension (transitions) :: dlp, dll, dlr, sp, sl, sr
  Real, Dimension (Params%n_points) :: sg, cg, s2g, c2g, sf, cf, B_str, inc, azi
  Real, Dimension (Params%n_points) :: Damp, Dldop, Line_op, Cont_op_5000
  Real, Dimension (nwlengths, Params%n_points) :: Phi_I, Phi_Q, Phi_U, Phi_V, &
       Psi_Q, Psi_U, Psi_V
  Real :: t13, t14, dldopcm, dldopcms, a, v, mag
  Double Precision :: etar, vetar, getar, ettar, ettvr, esar, vesar, gesar
  Double Precision :: essar, essvr, ettmr, essmr
  Double Precision :: etal, vetal, getal, ettal, ettvl, esal, vesal, gesal
  Double Precision :: essal, essvl, ettml, essml
  Double Precision :: etap, vetap, getap, ettap, ettvp, esap, vesap, gesap
  Double Precision :: essap, essvp, ettmp, essmp
  Double Precision :: tm, tn, sm, sn, tpm, spm
  Double Precision, Dimension(nwlengths) :: vtm, vtn, vsm, vsn, vtpm, vspm
  Double Precision, Dimension(nwlengths) :: vec_etar, vec_esar, vec_etal, &
       vec_esal, vec_etap, vec_esap
  Real :: dlo, wc
  Real, Dimension (2) :: tam
  Integer :: iwave
  Logical :: Incomplete, Exists
  Logical, Save :: FirstTime
! For Hyperfine Structure case
  Integer, Save :: Unit, nHFlines, ilin
  Real :: MagField
  Logical, Save :: Constant_B
  Integer, Dimension (:), Allocatable, Save :: HFIndex
  Real, Dimension(:), Allocatable, Save :: Alow, Blow, Aup, Bup, SpinI
! Special case definitions for 15652
  Real, Dimension (9) :: dlp15652, dll15652, dlr15652, &
       sp15652, sl15652 ,sr15652
  Data dlp15652/-4.4000149E-2,-3.3000112E-02,-2.2000074E-02,-1.1000037E-02, &
       0.00,1.1000037E-02,2.2000074E-02,3.3000112E-02,4.4000149E-02/
  Data dll15652/1.466000,1.477000,1.488000,1.499000, &
    1.510000,1.521000,1.532000,1.543000,1.554000/
  Data dlr15652/-1.554000,-1.543000,-1.532000,-1.521000, &
   -1.510000,-1.499000,-1.488000,-1.477000,-1.466000/
  Data sp15652/5.4545455E-02,9.6969694E-02,0.1272727,0.1454545, &
   0.1515152,0.1454545,0.1272727,9.6969694E-02,5.4545455E-02/
  Data sl15652/6.0606059E-03,1.8181818E-02,3.6363635E-02,6.0606062E-02, &
   9.0909094E-02,0.1272727,0.1696970,0.2181818,0.2727273/
  Data sr15652/0.2727273,0.2181818,0.1696970,0.1272727, &
   9.0909094E-02, 6.0606062E-02, 3.6363635E-02, 1.8181818E-02, 6.0606059E-03/
  Data FirstTime/.TRUE./
!
  Call Time_routine('profiles',.True.)
  npoints=Params%n_points
  Do ind=1, npoints
     If (Abs(Atmo%B_x(ind)) .lt. ZeroFieldThreshold) &
          Atmo%B_x(ind)=ZeroFieldThreshold
     If (Abs(Atmo%B_y(ind)) .lt. ZeroFieldThreshold) &
          Atmo%B_y(ind)=ZeroFieldThreshold
     If (Abs(Atmo%B_long(ind)) .lt. ZeroFieldThreshold) &
          Atmo%B_long(ind)=ZeroFieldThreshold
  End do
  B_str=Sqrt(Atmo%B_long**2+Atmo%B_x**2+Atmo%B_y**2)
  inc=atan2(Sqrt(Atmo%B_x**2+Atmo%B_y**2),Atmo%B_long)/Pi*180.
  azi=atan2(Atmo%B_y, Atmo%B_x)/Pi*180.
  If (FirstTime) then
     Constant_B=.FALSE.
     If (Line%Hyperfine) then
        If (MaxVal(B_str)-MinVal(B_str) .lt. 10) Constant_B=.TRUE.
     End if
!     FirstTime=.FALSE.
  End if
!
  Incomplete=.FALSE.
  Phi_I(:,:)=0.
  Phi_Q(:,:)=0.
  Phi_U(:,:)=0.
  Phi_V(:,:)=0.
  Psi_Q(:,:)=0.
  Psi_U(:,:)=0.
  Psi_V(:,:)=0.
  t13=0. ! Derivative of log (Dldop) with temperature
  t14=0. ! Derivative of log (Dldop) with microturbulence
  ji(1)=int(Line%J_low)
  ji(2)=int(Line%J_up)
  jf(1)=int(10*(Line%J_low-ji(1)))
  jf(2)=int(10*(Line%J_up-ji(2)))
  tam(1)=Line%J_low
  tam(2)=Line%J_up
  mult(1)=Line%Mult_low
  mult(2)=Line%Mult_up
  desig(1)=Line%Desig_low
  desig(2)=Line%Desig_up
  dlo=4.6686E-5*(Line%Wlength**2.)*1.E-16 ! Lambda_B in mA
! Special case for FeI 15652 A ?
  If (Line%Wlength .ge. 15652.7 .and. Line%Wlength .le. 15652.9 &
       .and. Line%Atomic_number .eq. 26 .and. Line%Ion_stage .eq. 1) then 
     If (FirstTime .and. Params%Printout .ge. 1) &
          Print *,'Special case for FeI 15652'
     np=9
     nl=9
     nr=9
     dlp(1:np)=dlp15652(1:np)*dlo
     dll(1:nl)=dll15652(1:nl)*dlo
     dlr(1:nr)=dlr15652(1:nr)*dlo
     sp(1:np)=sp15652(1:np)
     sl(1:nl)=sl15652(1:nl)
     sr(1:nr)=sr15652(1:nr)
  Else
     Call Zeeman(transitions, mult, desig, tam, ji, jf, dlo, np, nl, nr, &
          dlp, dll, dlr, sp, sl, sr)
     If (np .gt. transitions .or. nl .gt. transitions .or. &
          nr .gt. transitions) then
        Print *,'Too many transitions (',np,',', nl,',', nr,')'
        Print *,'Profiles has been prepared to work with ',transitions, &
             '. STOP in Profile'
        Stop
     End if
  End if
  sg(1:npoints)=Sin(inc(1:npoints)*pi/180.)
  cg(1:npoints)=Cos(inc(1:npoints)*pi/180.)
  s2g(1:npoints)=sg(1:npoints)*sg(1:npoints)
  c2g(1:npoints)=1.+cg(1:npoints)*cg(1:npoints)
  sf(1:npoints)=Sin(2.*azi(1:npoints)*pi/180.)
  cf(1:npoints)=Cos(2.*azi(1:npoints)*pi/180.)
  wc=Line%Wlength*1.E-8/cc ! Inverse of line frequency (in s)
  Do iwave=1, nwlengths
     DWave(iwave)=Wave(iwave)-Line%Wlength
  End do
  If (MinVal(Abs(DWave(:))) .lt. Line%Width) then
     Do ind=npoints, 1, -1
        a=Damp(ind)
        dldopcm=Dldop(ind)*1.E-8 ! Dldop in cm
        dldopcms=Dldop(ind)/Line%Wlength*cc ! Dldop in cm/s
        vv(1:nw)=DWave(IndexWave(1:nw))/Dldop(ind) - Atmo%v_los(ind)/dldopcms
        mag=B_str(ind)/dldopcm
        If (Line%Hyperfine) then
           If (.not. Constant_B .or. ind .eq. npoints) then
              Call Zeeman_hyperfine(line, B_str(ind), np, nl, nr, &
                   dlp, dll, dlr, sp, sl, sr)
           End if
           mag=1.d-8*line%Wlength**2/dldop(ind)
        End if
        If (Params%speed .eq. 0) then ! Use old point-by-point voigt
           Do iwave=1, nw
              v=vv(iwave)
              call mvoigt(nr,dlr,sr,a,v,mag,t13,t14,wc,dldopcm,etar,vetar, &
                   getar,ettar,ettvr,esar,vesar,gesar,essar,essvr, &
                   ettmr,essmr)
              call mvoigt(nl,dll,sl,a,v,mag,t13,t14,wc,dldopcm,etal,vetal, &
                   getal,ettal,ettvl,esal,vesal,gesal,essal,essvl, &
                   ettml,essml)
              call mvoigt(np,dlp,sp,a,v,mag,t13,t14,wc,dldopcm,etap,vetap, &
                   getap,ettap,ettvp,esap,vesap,gesap,essap,essvp, &
                   ettmp,essmp)
              tm=.5*(etar+etal)
              tn=.5*(etar-etal)
              sm=.5*(esar+esal)
              sn=.5*(esar-esal)
              tpm=.5*(etap-tm)
              spm=.5*(esap-sm)
              Phi_I(IndexWave(iwave),ind)=Line_op(ind)*0.5*(etap*s2g(ind)+tm*c2g(ind))
              Phi_Q(IndexWave(iwave),ind)=Line_op(ind)*tpm*s2g(ind)*cf(ind)
              Phi_U(IndexWave(iwave),ind)=Line_op(ind)*tpm*s2g(ind)*sf(ind)
              Phi_V(IndexWave(iwave),ind)=Line_op(ind)*tn*cg(ind)
              Psi_Q(IndexWave(iwave),ind)=Line_op(ind)*spm*s2g(ind)*cf(ind)
              Psi_U(IndexWave(iwave),ind)=Line_op(ind)*spm*s2g(ind)*sf(ind)
              Psi_V(IndexWave(iwave),ind)=Line_op(ind)*sn*cg(ind)
           End do
        Else ! Use new vector voigt function
           Call mvoigt2(nr,nw,dlr,sr,a,vv,mag,wc,dldopcm,vec_etar,vec_esar)
           Call mvoigt2(nl,nw,dll,sl,a,vv,mag,wc,dldopcm,vec_etal,vec_esal)
           Call mvoigt2(np,nw,dlp,sp,a,vv,mag,wc,dldopcm,vec_etap,vec_esap)
           vtm(1:nw)=.5*(vec_etar(1:nw)+vec_etal(1:nw))
           vtn(1:nw)=.5*(vec_etar(1:nw)-vec_etal(1:nw))
           vsm(1:nw)=.5*(vec_esar(1:nw)+vec_esal(1:nw))
           vsn(1:nw)=.5*(vec_esar(1:nw)-vec_esal(1:nw))
           vtpm(1:nw)=.5*(vec_etap(1:nw)-vtm(1:nw))
           vspm(1:nw)=.5*(vec_esap(1:nw)-vsm(1:nw))
           Phi_I(IndexWave(1:nw),ind)=Line_op(ind)*0.5*(vec_etap(1:nw)*s2g(ind)+vtm(1:nw)*c2g(ind))
           Phi_Q(IndexWave(1:nw),ind)=Line_op(ind)*vtpm(1:nw)*s2g(ind)*cf(ind)
           Phi_U(IndexWave(1:nw),ind)=Line_op(ind)*vtpm(1:nw)*s2g(ind)*sf(ind)
           Phi_V(IndexWave(1:nw),ind)=Line_op(ind)*vtn(1:nw)*cg(ind)
           Psi_Q(IndexWave(1:nw),ind)=Line_op(ind)*vspm(1:nw)*s2g(ind)*cf(ind)
           Psi_U(IndexWave(1:nw),ind)=Line_op(ind)*vspm(1:nw)*s2g(ind)*sf(ind)
           Psi_V(IndexWave(1:nw),ind)=Line_op(ind)*vsn(1:nw)*cg(ind)
        End if
        If (Phi_I(1, ind) .gt. 1.e-2*Cont_op_5000(ind)) &
             Incomplete=.TRUE.
        If (Phi_I(nwlengths, ind) .gt. 1.e-2*Cont_op_5000(ind)) &
             Incomplete=.TRUE.
     End do ! Do ind=npoints, 1, -1
  End if ! If (MaxVal(Abs(Dwave(:))) .gt. Line%Width)
  Call Time_routine('profiles',.False.)
  Return
End Subroutine Profiles

! This subroutine solves the scalar (unpolarized) radiative transfer
! equation using the Hermitian solution of Bellot Rubio, Ruiz Cobo and 
! Collados 1998, ApJ, 506, 805.
!
! The derivatives of the source function is computed assuming a 
! parabolic form.
!
Subroutine Scalar_hermite(npoints, tau_nu, S, Intensity)
  Implicit None
  Integer :: npoints, ibound, ipoint
  Real, Dimension (npoints) :: tau_nu, S, Intensity
  Real, Dimension (npoints) :: S_p
  Real :: A, E, J, J_down, Cal_J, Cal_J_down, D, dtau
  Real, Parameter :: Optically_thick = 100, Optically_thin=1.e-2
  Logical, Save :: Warning_thin=.TRUE., Warning_thick=.TRUE.
  logical :: checknan
!
! Calculate derivatives of S
!
  Call Parab_deriv(npoints, tau_nu, S, S_p)
!
! Find the point for the boundary condition
!
  ibound=npoints-1
  Do while (tau_nu(ibound) .gt. Optically_thick .and. ibound .gt. 3)
     ibound=ibound-1
  End do
  If (tau_nu(ibound) .lt. 10. .and. Warning_thick) then
     Print *,'Warning! Line not optically thick at the bottom'
     Print *,'Tau_nu at the bottom = ',tau_nu(ibound)
     Print *,'This warning will not be shown any more during this run'
     Warning_thick=.FALSE.
  End if
!
! Set boundary condition (I=S + dS/d(tau_nu)) at ibound
!
  Intensity(ibound:npoints)=S(ibound:npoints)!+S_p(ibound:npoints)
!
! Ray propagation
!
  ipoint=ibound-1
  Do while (ipoint .gt. 1 .and. tau_nu(ipoint) .gt. Optically_thin) ! Loop
     dtau=tau_nu(ipoint+1)-tau_nu(ipoint)
!
! Down quantities (point "i" in the paper)
!
     J_down=dtau/2.*S(ipoint+1)    
     Cal_J_down=dtau*dtau/12.*(S_p(ipoint+1) - S(ipoint+1))
!
! Up quantities (point "i+1" in the paper)
!
     J=dtau/2.*S(ipoint)
     Cal_J=dtau*dtau/12.*(S_p(ipoint) - S(ipoint))
     D=1. + dtau/2. + dtau*dtau/12.
     A=(1. - dtau/2. + dtau*dtau/12.) / D
     E=(J + J_down + Cal_J_down - Cal_J) / D
     Intensity(ipoint) = A*Intensity(ipoint+1) + E
     ipoint=ipoint-1
  End do
 ! End loop at the top of the atmo or where tau_nu < Optically_thin
  If (tau_nu(1) .gt. 1.e-2 .and. Warning_thin) then
     Print *,'Warning! Line not optically thin in the surface'
     Print *,'Tau_nu at the surface = ',tau_nu(1)
     Print *,'This warning will not be shown any more during this run'
     Warning_thin=.FALSE.
  End if
  Intensity(1:ipoint)=Intensity(ipoint+1)
  Return
End Subroutine Scalar_hermite
! This function computes the HSRA continuum at the specified wavelength (in 
! Angstroms) and heliocentric angle (cos [mu]). The source function and optical
! depth scales are tabulated (output from LILIA v4.2) and the formal solution 
! is carried out using the scalar Hermitian method. Output is intensity in 
! cgs units.

Real function syn_hsra(lambda, theta)
  Use Background_opacity_module, Only: Background_opacity
  Use Eq_state, Only: compute_others_from_T_Pe_Pg
  Implicit None
  Real :: lambda, theta, metal
  Real :: Wave_cm, nu, nu500, n2P, Scat
  Real, Parameter :: cc = 2.99792458D10, hh = 6.62618D-27, bk = 1.38066D-16
  Integer :: idepth
  Integer, Parameter :: npoints=44
  Real, Dimension (npoints) :: Temp, el_p, ltau_500, ltau_nu, gas_p
  Real, Dimension (npoints) :: Source_f, Cont_op, cont_op_500, &
       tau_nu, tau_500, Intensity
  Real, Dimension (npoints) :: nH, nHminus, nHplus, nH2, nH2plus
  Real, Dimension (10) :: Pp
! HSRA model around ltau_500=1
  Data ltau_500/  -2.900000     ,   -2.800000     ,   -2.700000     ,   -2.600000     , &
   -2.500000     ,   -2.400000     ,   -2.300000     ,   -2.200000     , &
   -2.100000     ,   -2.000000     ,   -1.900000     ,   -1.800000     , &
   -1.700000     ,   -1.600000     ,   -1.500000     ,   -1.400000     , &
   -1.300000     ,   -1.200000     ,   -1.100000     ,   -1.000000     , &
  -0.9000000     ,  -0.8000000     ,  -0.7000000     ,  -0.6000000     , &
  -0.5000000     ,  -0.4000000     ,  -0.3000000     ,  -0.2000000     , &
  -0.1000000     ,    0.000000     ,   0.1000000     ,   0.2000000     , &
   0.3000000     ,   0.4000000     ,   0.5000000     ,   0.6000000     , &
   0.7000000     ,   0.8000000     ,   0.9000000     ,    1.000000     , &
    1.100000     ,    1.200000     ,    1.300000     ,    1.400000     /
  Data Temp/   4405.000     ,    4430.000     ,    4460.000     ,    4490.000     , &
    4525.000     ,    4550.000     ,    4575.000     ,    4600.000     , &
    4630.000     ,    4660.000     ,    4690.000     ,    4720.000     , &
    4750.000     ,    4790.000     ,    4840.000     ,    4895.000     , &
    4950.000     ,    5010.000     ,    5080.000     ,    5160.000     , &
    5240.000     ,    5330.000     ,    5430.000     ,    5540.000     , &
    5650.000     ,    5765.000     ,    5890.000     ,    6035.000     , &
    6200.000     ,    6390.000     ,    6610.000     ,    6860.000     , &
    7140.900     ,    7440.000     ,    7750.000     ,    8030.000     , &
    8290.000     ,    8520.000     ,    8710.000     ,    8880.000     , &
    9050.000     ,    9220.000     ,    9390.000     ,    9560.000  /
  Data el_p/  0.2727410     ,   0.3088174     ,   0.3554601     ,   0.4084350     , &
   0.4706356     ,   0.5313036     ,   0.6006021     ,   0.6796458     , &
   0.7814314     ,   0.8967699     ,    1.016245     ,    1.165931     , &
    1.320946     ,    1.527207     ,    1.776870     ,    2.072720     , &
    2.415603     ,    2.824712     ,    3.329539     ,    3.961881     , &
    4.727736     ,    5.719847     ,    7.043319     ,    8.864668     , &
    11.17472     ,    14.48614     ,    19.20516     ,    26.47779     , &
    37.89124     ,    56.54964     ,    88.19545     ,    142.3658     , &
    235.7640     ,    388.8344     ,    629.1764     ,    944.5691     , &
    1347.037     ,    1814.897     ,    2300.430     ,    2825.199     , &
    3445.770     ,    4175.081     ,    5028.727     ,    6020.554  /
  Data gas_p/  3765.968     ,    4254.446     ,    4876.757     ,    5580.309     , &
    6374.859     ,    7184.000     ,    8110.760     ,    9170.515     , &
    10516.81     ,    12036.57     ,    13571.12     ,    15534.97     , &
    17518.44     ,    20044.50     ,    22887.84     ,    26098.81     , &
    29734.65     ,    33845.06     ,    38483.77     ,    43714.83     , &
    49620.61     ,    56240.48     ,    63611.58     ,    71732.42     , &
    79563.52     ,    89252.80     ,    99432.09     ,    109863.3     , &
    120199.5     ,    130097.9     ,    139246.6     ,    147431.7     , &
    154576.7     ,    160726.7     ,    166037.4     ,    170733.3     , &
    175078.0     ,    179243.7     ,    183409.1     ,    187728.1     , &
    192268.1     ,    197049.3     ,    202219.2     ,    207666.5  /  
  Data nH/  5.6194241E+15 , 6.3123542E+15 , 7.1868574E+15 , 8.1685392E+15 ,&
   9.2592580E+15 , 1.0376913E+16 , 1.1651208E+16 , 1.3101552E+16,&
   1.4927041E+16 , 1.6973615E+16 , 1.9014572E+16 , 2.1626889E+16,&
   2.4233335E+16 , 2.7495303E+16 , 3.1070867E+16 , 3.5031665E+16,&
   3.9468203E+16 , 4.4386224E+16 , 4.9774994E+16 , 5.5665680E+16,&
   6.2222553E+16 , 6.9334632E+16 , 7.6979841E+16 , 8.5085836E+16,&
   9.2538910E+16 , 1.0173731E+17 , 1.1093311E+17 , 1.1962041E+17,&
   1.2737926E+17 , 1.3374387E+17 , 1.3833491E+17 , 1.4103766E+17,&
   1.4189576E+17 , 1.4134895E+17 , 1.3977811E+17 , 1.3821148E+17,&
   1.3665632E+17 , 1.3542692E+17 , 1.3485021E+17 , 1.3465291E+17,&
   1.3447756E+17 , 1.3431769E+17 , 1.3425305E+17 , 1.3418188E+17 /
  Data nHplus/  2.4705971E+09 , 3.0427884E+09 , 3.8894730E+09 , 4.9554954E+09,&
   6.5225533E+09 , 7.9511496E+09 , 9.6766044E+09 , 1.1756768E+10 ,&
   1.4787671E+10 , 1.8543409E+10 , 2.3131349E+10 , 2.8854268E+10 ,&
   3.5807687E+10 , 4.7355167E+10 , 6.6335642E+10 , 9.5125291E+10 ,&
   1.3528567E+11 , 1.9640495E+11 , 2.9854938E+11 , 4.7222243E+11 ,&
   7.3312934E+11 , 1.1714653E+12 , 1.9084469E+12 , 3.1375323E+12 ,&
   4.9500604E+12 , 7.7060170E+12 , 1.1952516E+13 , 1.8907509E+13 ,&
   3.0178338E+13 , 4.8791047E+13 , 8.0094622E+13 , 1.3246594E+14 ,&
   2.1984871E+14 , 3.5763283E+14 , 5.6527863E+14 , 8.2739250E+14 ,&
   1.1504122E+15 , 1.5145412E+15 , 1.8829566E+15 , 2.2727187E+15 ,&
   2.7243826E+15 , 3.2446717E+15 , 3.8418700E+15 , 4.5223778E+15 /
  Data nHminus/  6507388.   ,    8069774.  ,    1.0261648E+07 , 1.3007725E+07,&
   1.6414706E+07 , 2.0267806E+07 , 2.5110990E+07 , 3.1196458E+07 ,&
   3.9717520E+07 , 5.0385260E+07 , 6.2197508E+07 , 7.8943320E+07 ,&
   9.7500912E+07 , 1.2334002E+08 , 1.5505795E+08 , 1.9426640E+08 ,&
   2.4317294E+08 , 3.0379898E+08 , 3.7864086E+08 , 4.7177434E+08 ,&
   5.9003590E+08 , 7.4105338E+08 , 9.3826496E+08 , 1.2022258E+09 ,&
   1.5215749E+09 , 1.9991647E+09 , 2.6521277E+09 , 3.5799480E+09 ,&
   4.9065283E+09 , 6.8368241E+09 , 9.6835512E+09 , 1.3843898E+10 ,&
   1.9853416E+10 , 2.8048089E+10 , 3.8720676E+10 , 5.0664645E+10 ,&
   6.3905403E+10 , 7.7657268E+10 , 9.0949157E+10 , 1.0454604E+11 ,&
   1.1961636E+11 , 1.3625697E+11 , 1.5471411E+11 , 1.7496125E+11 /
  Data nH2/  2.0795303E+12 , 2.4565381E+12 , 2.9337308E+12 , 3.5061314E+12 ,&
   4.1329657E+12 , 4.8697862E+12 , 5.7788276E+12 , 6.8729999E+12 ,&
   8.3156958E+12 , 9.9739899E+12 , 1.1681122E+13 , 1.4111733E+13 ,&
   1.6529361E+13 , 1.9460996E+13 , 2.2251781E+13 , 2.5128151E+13 ,&
   2.8428456E+13 , 3.1779828E+13 , 3.4718800E+13 , 3.7187492E+13 ,&
   3.9976381E+13 , 4.2163618E+13 , 4.3630392E+13 , 4.4308850E+13 ,&
   4.3885095E+13 , 4.4396129E+13 , 4.3866061E+13 , 4.1571655E+13 ,&
   3.7799361E+13 , 3.2801676E+13 , 2.7082621E+13 , 2.1437214E+13 ,&
   1.6350040E+13 , 1.2313220E+13 , 9.2649030E+12 , 7.2831457E+12 ,&
   5.8975275E+12 , 4.9536125E+12 , 4.3455430E+12 , 3.9019838E+12 ,&
   3.5198224E+12 , 3.1886218E+12 , 2.9036980E+12 , 2.6534952E+12 /
!
! Compute S and continuum opacities for HSRA
  nu=cc/(lambda*1e-8)! s^-1
  Wave_cm=cc/nu
  Source_f(1:npoints)=2.*hh/Wave_cm*cc*cc/(Wave_cm**4)/  & ! LTE
       (exp(hh*nu/bk/Temp(1:npoints))-1.) ! Planck function (c.g.s.)
  nu500=cc/(5000*1e-8)! s^-1
!
!  Call Compute_others_from_T_Pe_Pg(npoints, Temp, El_p, Gas_p, &
!       nH, nHminus, nHplus, nH2, nH2plus)     
  Do idepth=1, npoints
     n2P=BK*Temp(idepth)
     Cont_op(idepth)=Background_opacity(Temp(idepth), El_p(idepth), Gas_p(idepth), nH(idepth)*n2P, &
          nHminus(idepth)*n2P, nHplus(idepth)*n2P, nH2(idepth)*n2P, 0., lambda, Scat)
     Cont_op_500(idepth)=Background_opacity(Temp(idepth), El_p(idepth), Gas_p(idepth), nH(idepth)*n2P, &
          nHminus(idepth)*n2P, nHplus(idepth)*n2P, nH2(idepth)*n2P, 0., 5000., Scat)
  End do
  Cont_op=Cont_op/Cont_op_500
!
! Construct the tau_nu (monochromatic optical depth)
!
  tau_500(1:npoints)=10**(ltau_500(1:npoints)-alog10(theta))
  tau_nu(1)=tau_500(1)*Cont_op(1)
  Do idepth=2, npoints
     tau_nu(idepth)=tau_nu(idepth-1)+ &
          (tau_500(idepth)-tau_500(idepth-1))* &
          0.5*(Cont_op(idepth)+Cont_op(idepth-1))
  End do
!
  Call Scalar_hermite(npoints, tau_nu, Source_f, Intensity)
!
  syn_hsra=Intensity(1)
  Return 
!
End function syn_hsra
! This subroutine solves the Stokes radiative transfer equation as 4 
! independent scalar equations, using the Weakly Polarizing Media
! approximation of Sanchez Almeida and Trujillo Bueno 1999, ApJ, in
! press.
!
Subroutine WPM(npoints, dtau_nu, Ab, Source_f, Emergent_stokes)
  Implicit None
  Integer :: npoints, ipoint
  Real, Dimension (npoints) :: ltau_500, Source_f, dtau_nu
  Real, Dimension (npoints) :: Opac, tau_500, Stokes_component
  Real, Dimension (npoints,4,4) :: Ab, Ab2
  Real, Dimension (4) :: Emergent_stokes, stk
!
! Solve for Stokes I without polarization
!
!  Opac(1:npoints)=Ab(1:npoints,1,1)
!  tau_nu(1)=tau_500(1)*Opac(1)
!  Do ipoint=2, npoints
!     tau_nu(ipoint)=tau_nu(ipoint-1)+ &
!          (tau_500(ipoint)-tau_500(ipoint-1))* &
!          0.5*(Opac(ipoint)+Opac(ipoint-1))
!  End do
!  Call Check_tau(tau_nu)
  Call delobezier_scal(npoints, dtau_nu, Ab, Source_f, stk)
!  Call Scalar_hermite(npoints, tau_nu, Source_f, Stokes_component)
  Emergent_stokes(1)=stk(1)
!
! Now solve for I + V and subtract the I computed above
!
  Ab2(1:npoints,1,1)=Ab(1:npoints,1,1) + Ab(1:npoints,1,4)
!  tau_nu(1)=tau_500(1)*Opac(1)
!  Do ipoint=2, npoints
!     tau_nu(ipoint)=tau_nu(ipoint-1)+ &
!          (tau_500(ipoint)-tau_500(ipoint-1))* &
!          0.5*(Opac(ipoint)+Opac(ipoint-1))
!  End do 
!  Call Check_tau(tau_nu)
  Call delobezier_scal(npoints, dtau_nu, Ab2, Source_f, stk)
!  Call Scalar_hermite(npoints, tau_nu, Source_f, Stokes_component)
  Emergent_stokes(4)=stk(1) - Emergent_stokes(1) ! I+V - I
!
! Now solve for I + Q and subtract the I computed above
!
  Ab2(1:npoints,1,1)=Ab(1:npoints,1,1) + Ab(1:npoints,1,2)
!  Opac(1:npoints)=Ab(1:npoints,1,1) + Ab(1:npoints,1,2)
!  tau_nu(1)=tau_500(1)*Opac(1)
!  Do ipoint=2, npoints
!     tau_nu(ipoint)=tau_nu(ipoint-1)+ &
!          (tau_500(ipoint)-tau_500(ipoint-1))* &
!          0.5*(Opac(ipoint)+Opac(ipoint-1))
!  End do 
!  Call Check_tau(tau_nu)
  Call delobezier_scal(npoints, dtau_nu, Ab2, Source_f, stk)
!  Call Scalar_hermite(npoints, tau_nu, Source_f, Stokes_component)
  Emergent_stokes(2)=stk(1) - Emergent_stokes(1) ! I+Q - I
!
! Now solve for I + U and subtract the I computed above
!
  Ab2(1:npoints,1,1)=Ab(1:npoints,1,1) + Ab(1:npoints,1,3)
!  Opac(1:npoints)=Ab(1:npoints,1,1) + Ab(1:npoints,1,3)
!  tau_nu(1)=tau_500(1)*Opac(1)
!  Do ipoint=2, npoints
!     tau_nu(ipoint)=tau_nu(ipoint-1)+ &
!          (tau_500(ipoint)-tau_500(ipoint-1))* &
!          0.5*(Opac(ipoint)+Opac(ipoint-1))
!  End do 
!  Call Check_tau(tau_nu)
  Call delobezier_scal(npoints, dtau_nu, Ab2, Source_f, stk)
!  Call Scalar_hermite(npoints, tau_nu, Source_f, Stokes_component)
  Emergent_stokes(3)=stk(1) - Emergent_stokes(1) ! I+U - I
! Done!
  Return
End Subroutine WPM

! This subroutine solves the Stokes radiative transfer equation using the
! Hermitian solution of Bellot Rubio, Ruiz Cobo and Collados 1998, 
! ApJ, 506, 805.
! The derivatives of the source function and the absorption matrix are
! computed assuming a parabolic form for them.
!
  Subroutine hermite(npoints, dtau_nu, Ab, S, stokes)
    USE Bezier_math
    Implicit None
    Integer :: ipoint, npoints, ibound, ii, jj, k, k1, info, ubound
    Integer :: ipiv(4)
    Real, dimension(npoints) :: tau_nu, S, iab, dtau_nu, tau_500
    Real, dimension(npoints, 4, 4) :: Ab, abp
    Real, dimension(4,npoints) ::  sv, opac, sp
    Real, dimension(4) :: stokes, Jb(4), Ja(4), sJa(4), sJb(4),  sKa(4,4), sKb(4,4)
    Real (kind=8), dimension(4,4) :: ident, A, tmpa, tmpb, k2, D
    Real (kind=8) :: dtau, dtau1, dtau2, idtau2, eps, dtau05, tt
    Real (kind=8) :: alpha, beta, gamma, lambda, ostokes(4), E(4)
    Real, Parameter :: Optically_thick = 100., Optically_thin=1.e-4

    !
    ! 4x4 Identity matrix
    !
    Data ident(:,1)/1.0d0,0.0d0,0.0d0,0.0d0/
    Data ident(:,2)/0.0d0,1.0d0,0.0d0,0.0d0/
    Data ident(:,3)/0.0d0,0.0d0,1.0d0,0.0d0/
    Data ident(:,4)/0.0d0,0.0d0,0.0d0,1.0d0/

    iab = ab(:,1,1)


    !
    ! Find the point for the boundary condition
    !
    ubound = 1
    ibound = npoints
    tt=dtau_nu(1)
    do k = 2, npoints-1
       If (tt+dtau_nu(k) .ge. (Optically_Thin) .and. tt .lt. (Optically_Thin)) &
            ubound=k-1
       If (tt+dtau_nu(k) .ge. (Optically_Thick) .and. tt .lt. (Optically_Thick)) &
            ibound=k
       tt=tt+dtau_nu(k)
    End do
    if(ibound .lt. npoints) ibound = ibound + 1
    if(ubound .gt. 1) ubound = ubound - 1


    !
    ! Source function
    !
    Sv(2:4,:) = 0
    Sv(1,:) = S

       do jj = 1,4
          do ii = 1, 4
             ab(:,ii,jj) = ab(:,ii,jj) / iab
          end do
       end do



    !
    ! Set boundary condition (I=S) at ibound
    !
    ostokes(2:4)=0.0d0
    ostokes(1) = S(ibound) ! Boundary cond.

    !
    ! Derivatives of the source function and ab matrix
    !
    tau_nu(1)=0.
    do k=2,npoints
       tau_nu(k)=tau_nu(k-1)+dtau_nu(k-1)
    end do
    call full_cent_der_d(npoints, dtau_nu, Sv, sp, ubound, ibound)
    call full_cent_der4_d(npoints, dtau_nu, ab, abp, ubound, ibound)

    ibound = ibound - 1

    do k = ibound, ubound, -1
       !
       k1 = k + 1
       !
       ! dtau(s)
       !
       dtau = dtau_nu(k1) ! (tau_nu(k1) - tau_nu(k))
       dtau2 = dtau * dtau
       dtau05 = dtau * 0.5d0

       !
       ! at points "a" (k+1)
       !
       sKa  = matmul(ab(k1,:,:),ab(k1,:,:)) - abp(k1,:,:)
       sJa = dtau2/12.0d0 * (matmul(ab(k1,:,:), sp(:,k1)) - matmul(ska,Sv(:,k1)))
       Ja = dtau05 * matmul(ab(k1,:,:), Sv(:,k1))
       

       !
       ! at points "b" (k)
       !
       sKb = matmul(ab(k,:,:),ab(k,:,:)) - abp(k,:,:)
       sJb = dtau2/12.0d0 * (matmul(ab(k,:,:), sp(:,k)) - matmul(sKb,Sv(:,k)))
       Jb = dtau05 * matmul(ab(k,:,:), Sv(:,k))

       !
       ! Compute intensity
       !
       D = ident + dtau05*ab(k,:,:) + dtau2/12.d0 * sKb
       call matinx8(D)
       
       E = matmul(D, (Jb + Ja + sJa - sJb))
       A = matmul(D, (ident - dtau05 * ab(k1,:,:) + dtau2/12.d0 * sKa))
       
       ostokes = matmul(A, ostokes) + E

    END do

    stokes = ostokes
  END Subroutine hermite


Subroutine myold_Hermite(npoints, ltau_500, Ab, S, Emergent_stokes)
  Implicit None
  Integer :: ipoint, npoints, ibound, jj, kk
  Real, Dimension (npoints) :: ltau_500, S, S_p, Tmp, tau_nu, &
       tau_500
  Real, Dimension (npoints, 4) :: Stokes
  Real, Dimension (npoints, 4, 4) :: Ab, Ab_p
  Real, Dimension (4, 4) :: K2, A, Cal_K, D, Id, Tmp44
  Real, Dimension (4) :: Emergent_stokes, Tmp4
  Real, Dimension (4) :: J, Cal_J, E, J_down, Cal_J_down
  Real :: dtau
  Real, Parameter :: Optically_thick = 100., Optically_thin=1.e-4
  Logical, Save :: Warning_thin=.TRUE., Warning_thick=.TRUE.
  logical :: checknan
!
! First construct the tau_nu (monochromatic optical depth) and
! tau_5000 scales
!
  tau_500(1:npoints)=10**ltau_500(1:npoints)
  tau_nu(1)=(10**ltau_500(1))*Ab(1, 1, 1)
  Do ipoint=2, npoints
     tau_nu(ipoint)=tau_nu(ipoint-1)+ &
          (tau_500(ipoint)-tau_500(ipoint-1))* &
          0.5*(Ab(ipoint, 1, 1)+Ab(ipoint-1,1, 1))
  End do
!  Call Check_tau(tau_nu)
!
! Calculate the derivatives of Source funct. and Absorpt. mat.
!
  Call Parab_deriv(npoints, tau_500, Ab( 1:npoints,1,1), Tmp)
  Ab_p(1:npoints,1,1)=Tmp
  Call Parab_deriv(npoints, tau_500, Ab( 1:npoints,1,2), Tmp)
  Ab_p(1:npoints,1,2)=Tmp
  Call Parab_deriv(npoints, tau_500, Ab( 1:npoints,1,3), Tmp)
  Ab_p(1:npoints,1,3)=Tmp
  Call Parab_deriv(npoints, tau_500, Ab( 1:npoints,1,4), Tmp)
  Ab_p(1:npoints,1,4)=Tmp
  Call Parab_deriv(npoints, tau_500, Ab( 1:npoints,2,3), Tmp)
  Ab_p(1:npoints,2,3)=Tmp
  Call Parab_deriv(npoints, tau_500, Ab( 1:npoints,2,4), Tmp)
  Ab_p(1:npoints,2,4)=Tmp
  Call Parab_deriv(npoints, tau_500, Ab( 1:npoints,3,4), Tmp)
  Ab_p(1:npoints,3,4)=Tmp

! Invoke symmetry properties to construct the rest of the abs. matrix derivat.
  Ab_p( 1:npoints,2,2)=Ab_p( 1:npoints,1,1)
  Ab_p( 1:npoints,3,3)=Ab_p( 1:npoints,1,1)
  Ab_p( 1:npoints,4,4)=Ab_p( 1:npoints,1,1)
  Ab_p( 1:npoints,2,1)=Ab_p( 1:npoints,1,2)
  Ab_p( 1:npoints,3,1)=Ab_p( 1:npoints,1,3)
  Ab_p( 1:npoints,4,1)=Ab_p( 1:npoints,1,4)
  Ab_p( 1:npoints,3,2)=-Ab_p(1:npoints,2,3)
  Ab_p( 1:npoints,4,2)=-Ab_p( 1:npoints,2,4)
  Ab_p( 1:npoints,4,3)=-Ab_p( 1:npoints,3,4)
! Derivatives of source function (1st component only)
  Call Parab_deriv(npoints, tau_500, S, S_p)
!
! Find the point for the boundary condition
!
  ibound=npoints-1
  Do while (tau_nu(ibound) .gt. Optically_thick .and. ibound .gt. 3)
     ibound=ibound-1
  End do
  If (tau_nu(ibound) .lt. 10. .and. Warning_thick) then
     Print *,'Warning! Line not optically thick at the bottom'
     Print *,'Tau_nu at the bottom = ',tau_nu(ibound)
     Print *,'This warning will not be shown any more during this run'
     Warning_thick=.FALSE.
  End if
!
! Set boundary condition (I=S + dS/dtau) at ibound
!
  Stokes(ibound, 1)=S(ibound) + 0.*S_p(ibound)/ & ! Boundary cond.
       Ab(ibound,1, 1) ! Note that S_p is with respect to
                        ! tau_500 and we need it with tau_nu
  Stokes(ibound,2:4)=0.
!
! Ray propagation
!
  Id(1:4, 1:4)=0.
  Id(1,1)=1.
  Id(2,2)=1.
  Id(3,3)=1.
  Id(4,4)=1. ! Define identity matrix
  ipoint=ibound-1
  Do while (ipoint .gt. 1 .and. tau_nu(ipoint) .gt. Optically_thin) ! Loop
     dtau=tau_500(ipoint+1)-tau_500(ipoint)
!
! Down quantities (point "i" in the paper)
!
     J_down(1:4)=dtau/2.*Ab(ipoint+1,1:4, 1)*S(ipoint+1)
     Call Multiply_matrix(4, 4, 4, Ab(ipoint+1,1:4, 1:4), &
          Ab(ipoint+1,1:4, 1:4), K2) ! K^2
     Cal_J_down(1:4)=dtau*dtau/12.*(Ab(ipoint+1,1:4, 1)*S_p(ipoint+1)+ &
          Ab_p(ipoint+1,1:4, 1)*S(ipoint+1) - K2(1:4, 1)*S(ipoint+1))
!
! Up quantities (point "i+1" in the paper)
!
     J(1:4)=dtau/2.*Ab(ipoint,1:4, 1)*S(ipoint)
     Call Multiply_matrix(4, 4, 4, Ab(ipoint+1,1:4, 1:4), &
          Ab(ipoint,1:4, 1:4), K2) ! K^2
     Cal_J(1:4)=dtau*dtau/12.*(Ab(ipoint,1:4, 1)*S_p(ipoint)+ &
          Ab_p(ipoint,1:4, 1)*S(ipoint) - K2(1:4, 1)*S(ipoint))
     Cal_K(1:4, 1:4)=K2(1:4, 1:4) - Ab_p(ipoint,1:4, 1:4)
     D(1:4, 1:4)=Id(1:4, 1:4) + dtau/2.*Ab(ipoint,1:4, 1:4) + &
          dtau*dtau/12.*Cal_K(1:4, 1:4)
     Call matinx(D) ! Now the matrix D is actually the D^-1 in the paper
     Tmp44(1:4, 1:4)=Id(1:4, 1:4) - dtau/2.*Ab(ipoint+1, 1:4, 1:4) + &
          dtau*dtau/12.*Cal_K(1:4, 1:4)
     Call Multiply_matrix(4, 4, 4, D, Tmp44, A)
     Tmp4(1:4)=J(1:4)+J_down(1:4)+Cal_J_down(1:4)-Cal_J(1:4)
     Call Multiply_matrix(4, 4, 1, D, Tmp4, E)
     Call Multiply_matrix(4, 4, 1, A, Stokes(ipoint+1,1:4), Tmp4)
     Stokes(ipoint,1:4)=Tmp4(1:4) + E(1:4) ! I(i+1)=A.I(i)+E
     ipoint=ipoint-1
  End do ! End loop at the top of the atmo or where tau_nu < Optically_thin
  If (tau_nu(1) .gt. 1.e-2 .and. Warning_thin) then
     Print *,'Warning! Line not optically thin in the surface'
     Print *,'Tau_nu at the surface = ',tau_nu(1)
     Print *,'This warning will not be shown any more during this run'
     Warning_thin=.FALSE.
  End if
  Emergent_Stokes(1:4)=Stokes(ipoint,1:4)
  Return
End Subroutine myold_Hermite

Subroutine SC_solver(npoints, dtau_nu, Ab, S, Emergent_stokes)
! *** NOTE: For intensity ONLY. Returns Q=U=V=0.
  Implicit None
  Integer :: npoints, ibound, ipoint
  Real, Dimension (npoints) :: ltau_500, S, dtau_nu
  Real, Dimension (npoints) :: Opac, tau_nu, Stokes_component
  Real, Dimension (npoints,4, 4) :: Ab
  Real, Dimension (4) :: Emergent_stokes
  Real, Dimension (npoints) :: Intensity
  Real, Parameter :: Optically_thick = 100., Optically_thin=1.e-4
  Real :: T, Ex, Ex1, E0, E1, E2, D2, D3, D4, &
       ALF, BET, GAM, DELTAI, SDELTA, DELTAL
  Double Precision :: DTM, DTP, EXPDTP, EXPDTM
  Logical, Save :: Warning_thin=.TRUE., Warning_thick=.TRUE.
  Integer :: k
!
  Opac(1:npoints)=Ab(1:npoints,1,1)
!
! Find the point for the boundary condition
!
  tau_nu(1)=0.
  do k=2,npoints
     tau_nu(k)=tau_nu(k-1)+dtau_nu(k-1)
  end do
  ibound=npoints-1
  Do while (tau_nu(ibound) .gt. Optically_thick .and. ibound .gt. 3)
     ibound=ibound-1
  End do
  If (tau_nu(ibound) .lt. 10. .and. Warning_thick) then
     Print *,'Warning! Line not optically thick at the bottom'
     Print *,'Tau_nu at the bottom = ',tau_nu(ibound)
     Print *,'This warning will not be shown any more during this run'
     Warning_thick=.FALSE.
  End if
!
! Set boundary condition (I=S) at ibound
!
  Intensity(ibound)=S(ibound) ! Boundary cond.
  ipoint=ibound-1
  Do while (ipoint .gt. 1 .and. tau_nu(ipoint) .gt. Optically_thin) ! Loop
     dtm=dtau_nu(ipoint+1)
     dtp=dtau_nu(ipoint)
     EXPDTM=Exp(-DTM)
     E0=1.-EXPDTM
     E1=DTM-E0
     E2=DTM*DTM-2.*E1
     ALF=E0+(E2-(DTP+2*DTM)*E1)/(DTM*(DTM+DTP))
     BET=((DTM+DTP)*E1-E2)/DTM/DTP
     GAM=(E2-DTM*E1)/(DTP*(DTM+DTP))
     DELTAI=ALF*S(ipoint+1)+BET*S(ipoint)+GAM*S(ipoint-1)
     Intensity(ipoint)=Intensity(ipoint+1)*EXPDTM+DELTAI
     ipoint=ipoint-1
  End Do
  ipoint=ipoint+1
  if (ipoint .gt. npoints) ipoint=npoints
  Emergent_stokes(1)=Intensity(ipoint)
!
End Subroutine SC_solver
!
Subroutine Delopar(npoints, ltau_500, Ab_matrix, Source, Emergent_Stokes)
  Implicit None
  Integer :: npoints, ibound, ipoint, k, km, kp, i, j
  Real, Dimension (npoints) :: ltau_500, Source, tau_nu, dtau_nu, tau_500
  Real, Dimension (npoints,4) :: Opac
  Real, Dimension (npoints,4,4) :: Ab_matrix
  Real, Dimension (npoints,4) :: Stokes, Source_vector
  Real, Dimension (4) :: Emergent_Stokes
  Real :: chim, chi0, chip, dtp, dtm, exu
  Real :: psim, psi0, psip, psim_lin, psi0_lin
  Real, Dimension (4) :: sm, s0, sp
  Real, Dimension (4,4) :: mat1, mat2, mat3, identity
  Real, Parameter :: Optically_thick = 100., Optically_thin=1.e-4
  Logical, Save :: Warning_thin=.TRUE., Warning_thick=.TRUE.
!
! First construct the tau_5000 depth scale
!
  tau_500(1:npoints)=10**ltau_500(1:npoints)
!
  Opac(1:npoints,1)=Ab_matrix(1:npoints,1,1)
  tau_nu(1)=tau_500(1)*Opac(1,1)
  Do ipoint=2, npoints
     tau_nu(ipoint)=tau_nu(ipoint-1)+ &
          (tau_500(ipoint)-tau_500(ipoint-1))* &
          0.5*(Opac(ipoint,1)+Opac(ipoint-1,1))
  End do   
!  Call Check_tau(tau_nu)
!
  dtau_nu(2:npoints)=tau_nu(2:npoints)-tau_nu(1:npoints-1)
  dtau_nu(1)=dtau_nu(2)
!
! Find the point for the boundary condition
!
  ibound=npoints-1
  Do while (tau_nu(ibound) .gt. Optically_thick .and. ibound .gt. 3)
     ibound=ibound-1
  End do
  If (tau_nu(ibound) .lt. 10. .and. Warning_thick) then
     Print *,'Warning! Line not optically thick at the bottom'
     Print *,'Tau_nu at the bottom = ',tau_nu(ibound)
     Print *,'This warning will not be shown any more during this run'
     Warning_thick=.FALSE.
  End if
!
! Set boundary condition (I=S) at ibound
!
  Stokes(ibound:npoints,2:4)=0.
  Stokes(ibound:npoints,1)=Source(ibound:npoints) ! Boundary cond.
!
! Identity matrix
  identity(:,:)=0.
  Do i=1,4
     identity(i,i)=1.
  End do
!
  do i = 1, 4
     Opac(:,i)=Ab_matrix(:,1,i)
  End do
! Transform K into K* and then into K'
  do i = 1, 4
     do j = 1, 4
        ab_matrix(:,i,j) = ab_matrix(:,i,j) / opac(:,1)
     enddo
     ab_matrix(:,i,i) = ab_matrix(:,i,i) - 1.d0
     source_vector(:,i) = opac(:,i) * Source / opac(:,1)
  enddo
!
  k=ibound-1
  Do while (k .gt. 1 .and. tau_nu(k) .gt. Optically_thin) ! Loop
! Parabolic short-characteristics
     If (k .ne. 1) then
        km=k+1
        kp=k-1
        sm = source_vector(km,:)
        s0 = source_vector(k,:)
        sp = source_vector(kp,:)
        dtm = tau_nu(km)-tau_nu(k)
        dtp = tau_nu(k)-tau_nu(kp)
     else
! Linear short-characteristics
        km = k + 1
        sm = source_vector(km,:)
        s0 = source_vector(k,:)
        sp = 0.d0
        dtm = tau_nu(km)-tau_nu(k)
        dtp = 0.
     endif
!
     if (dtm >= 1.d-4) then
        exu = exp(-dtm)
     else
        exu = 1.d0 - dtm + 0.5d0 * dtm**2.d0
     endif
     call lin_sc(dtm,psim_lin,psi0_lin)
     mat1 = exu * identity - psim_lin*ab_matrix(km,:,:)
     mat2 = identity + psi0_lin * ab_matrix(k,:,:)
     call matinx(mat2)    

     if (k .ne. 1) then
        call par_sc(dtm,dtp,psim,psi0,psip)
        Stokes(k,:) = matmul(mat2,matmul(mat1,Stokes(km,:)) + psim*sm + psi0*s0 + psip*sp)
     else
        call lin_sc(dtm,psim,psi0)
        Stokes(k,:) = matmul(mat2,matmul(mat1,Stokes(km,:)) + psim*sm + psi0*s0)
     endif
     k=k-1
  enddo ! End while (k .gt. 1 .and ...)
!
  Do i=1,4
     Stokes(1:k,i)=Stokes(k+1,i)
  End do
  Emergent_stokes=Stokes(1,:)
!
  Return
  contains
    subroutine lin_sc(dtm, psim, psi0)
      real :: short_car_linea
      real, INTENT(IN) :: dtm
      real, INTENT(INOUT) :: psim, psi0
      real :: exu, u0, u1, c0, cm, d2
      
      if (dtm >= 1.d-4) then
         exu=exp(-dtm)
         u0=1.d0-exu
         u1=dtm-1.d0+exu
         
         c0=u1/dtm
         cm=u0-c0
      else
         d2=dtm**2.d0
         c0=(dtm/2.d0)-(d2/6.d0)
         cm=(dtm/2.d0)-(d2/3.d0)
      endif
      psi0 = c0
      psim = cm
      
    end subroutine lin_sc
    
    subroutine par_sc(dtm, dtp, psim, psi0, psip)
      real :: short_car_parab
      real, INTENT(IN) :: dtm, dtp
      real, INTENT(INOUT) :: psim, psi0, psip
      real :: exu, u0, u1 ,u2, d2, d3, d4
      
      if (dtm >= 1.d-4) then
         exu=exp(-dtm)
         u0=1.d0-exu
         u1=dtm-1.d0+exu
         u2=dtm**2-2.d0*dtm+2.d0-2.d0*exu
      else
         d2=dtm**2
         d3=dtm**3
         d4=dtm**4
         u0=dtm-(d2/2.d0)
         u1=(d2/2.d0)-(d3/6.d0)
         u2=(d3/3.d0)-(d4/12.d0)
      endif
      
      if (dtm * dtp /= 0.d0 .and. dtm**2 /= 0.d0 .and. dtp**2 /= 0.d0) then
         psim=(u2-u1*(dtp+2.d0*dtm))/(dtm*(dtm+dtp))+u0
         psi0=(u1*(dtm+dtp)-u2)/(dtm*dtp)
         psip=(u2-dtm*u1)/(dtp*(dtm+dtp))
      else
         psim = 0.d0
         psi0 = 0.d0
         psip = 0.d0
      endif
      
    end subroutine par_sc

End Subroutine Delopar

Subroutine Delolin(npoints, ltau_500, Ab_matrix, Source, Emergent_Stokes)
  Implicit None
  Integer :: npoints, ibound, ipoint, k, km, kp, i, j
  Real, Dimension (npoints) :: ltau_500, Source, tau_nu, dtau_nu, tau_500
  Real, Dimension (npoints,4) :: Opac
  Real, Dimension (npoints,4, 4) :: Ab_matrix
  Real, Dimension (npoints,4) :: Stokes, Source_vector
  Real, Dimension (4) :: Emergent_Stokes
  Real :: chim, chi0, dtm, exu
  Real :: psim, psi0, psip, psim_lin, psi0_lin
  Real, Dimension (4) :: sm, s0
  Real, Dimension (4,4) :: mat1, mat2, mat3, identity
  Real, Parameter :: Optically_thick = 100., Optically_thin=1.e-4
  Logical, Save :: Warning_thin=.TRUE., Warning_thick=.TRUE.
!
! First construct the tau_5000 depth scale
!
  tau_500(1:npoints)=10**ltau_500(1:npoints)
!
  Opac(1:npoints,1)=Ab_matrix(1:npoints,1,1)
  tau_nu(1)=tau_500(1)*Opac(1,1)
  Do ipoint=2, npoints
     tau_nu(ipoint)=tau_nu(ipoint-1)+ &
          (tau_500(ipoint)-tau_500(ipoint-1))* &
          0.5*(Opac(ipoint,1)+Opac(ipoint-1,1))
  End do   
!  Call Check_tau(tau_nu)
!
  dtau_nu(2:npoints)=tau_nu(2:npoints)-tau_nu(1:npoints-1)
  dtau_nu(1)=dtau_nu(2)
!
! Find the point for the boundary condition
!
  ibound=npoints-1
  Do while (tau_nu(ibound) .gt. Optically_thick .and. ibound .gt. 3)
     ibound=ibound-1
  End do
  If (tau_nu(ibound) .lt. 10. .and. Warning_thick) then
     Print *,'Warning! Line not optically thick at the bottom'
     Print *,'Tau_nu at the bottom = ',tau_nu(ibound)
     Print *,'This warning will not be shown any more during this run'
     Warning_thick=.FALSE.
  End if
!
! Set boundary condition (I=S) at ibound
!
  Stokes(ibound:npoints,2:4)=0.
  Stokes(ibound:npoints,1)=Source(ibound:npoints) ! Boundary cond.
!
! Identity matrix
  identity(:,:)=0.
  Do i=1,4
     identity(i,i)=1.
  End do
!
  do i = 1, 4
     Opac(:,i)=Ab_matrix(:,1,i)
  End do
! Transform K into K* and then into K'
  do i = 1, 4
     do j = 1, 4
        ab_matrix(:,i,j) = ab_matrix(:,i,j) / opac(:,1)
     enddo
     ab_matrix(:,i,i) = ab_matrix(:,i,i) - 1.d0
     source_vector(:,i) = opac(:,i) * Source / opac(:,1)
  enddo
!
  k=ibound-1
  Do while (k .gt. 1 .and. tau_nu(k) .gt. Optically_thin) ! Loop
! Linear short-characteristics
     km = k + 1
     sm = source_vector(km,:)
     s0 = source_vector(k,:)
     dtm = tau_nu(km)-tau_nu(k)
!
     if (dtm >= 1.d-4) then
        exu = exp(-dtm)
     else
        exu = 1.d0 - dtm + 0.5d0 * dtm**2.d0
     endif
     call lin_sc(dtm,psim_lin,psi0_lin)
     mat1 = exu * identity - psim_lin*ab_matrix(km,:,:)
     mat2 = identity + psi0_lin * ab_matrix(k,:,:)
     call matinx(mat2)    

     call lin_sc(dtm,psim,psi0)
     Stokes(k,:) = matmul(mat2,matmul(mat1,Stokes(km,:)) + psim*sm + psi0*s0)

     k=k-1
  enddo ! End while (k .gt. 1 .and ...)
!
  Do i=1,4
     Stokes(1:k,i)=Stokes(k+1,i)
  End do
  Emergent_stokes=Stokes(1,:)
!
  Return
  contains
    subroutine lin_sc(dtm, psim, psi0)
      real :: short_car_linea
      real, INTENT(IN) :: dtm
      real, INTENT(INOUT) :: psim, psi0
      real :: exu, u0, u1, c0, cm, d2
      
      if (dtm >= 1.d-4) then
         exu=exp(-dtm)
         u0=1.d0-exu
         u1=dtm-1.d0+exu
         
         c0=u1/dtm
         cm=u0-c0
      else
         d2=dtm**2.d0
         c0=(dtm/2.d0)-(d2/6.d0)
         cm=(dtm/2.d0)-(d2/3.d0)
      endif
      psi0 = c0
      psim = cm
      
    end subroutine lin_sc
    
End Subroutine Delolin

! For ipol_dscale
Function smooth(n, var, w) result(res)
  Implicit none
  Integer :: n, w, i, j, p0, p1, p
  Real (kind=8) :: tmp
  Real, dimension(n) :: var, res
  Real, dimension(:), allocatable :: pvar
  !
  ! Check that window is not larger than array!
  !
  if(w.GT.n) w = n
  p = w / 2
  
  allocate(pvar(n+p*2))
  !
  pvar(p+1:p+n) = var
  pvar(1:p+1) = var(2)
  pvar(p+n:n+2*p) = var(n)

  !
  res(:) = 0.0
  res(n) = var(n)
  !
  ! aind(1) is 0. Thus start at i = 2
  ! 
  tmp = sum(pvar(2:w+1)) / dble(w)

  res(2) = tmp
  do i = 3, n-1     
     tmp = tmp + ( pvar(i+w-1)-pvar(i-1)) / dble(w)
     res(i) = tmp
  enddo
  !
  deallocate(pvar)
  !
  return
END Function smooth
!
FUNCTION intepf(np, x, y, nxp, xp) Result(yp)
  Implicit None
  Integer :: np, nxp
  Real :: x(np), y(np), xp(nxp), yp(nxp)
  !
  Real (kind=8) :: lp1(np-1), lp2(np-1), fp1(np-1), fp2(np-1)
  Real (kind=8) :: xpi, xpi1, l1, l2, a0, a1, b0, b1
  Integer :: ii, i, i1, jj, np1
  !
  np1 = np - 1
  !
  ! Interpolation coefficients
  !
  lp1(1:np1) = 1.d0 / (x(1:np1) - x(2:np))
  lp2(1:np1) = 1.d0 / (x(2:np) - x(1:np1))
  !
  fp1(2:np1) = (y(3:np) - y(1:np-2)) / (x(3:np) - x(1:np-2))
  fp1(1) = (y(2) - y(1)) / (x(2) - x(1))
  !
  fp2(np1) = (y(np) - y(np1)) / (x(np) - x(np1))
  fp2(1:np-2) = fp1(2:np-1)
  !
  ! XP dependent part
  !
  do i = 1, np
     !
     i1 = i - 1
     if(i1 .LT. 1) i1 = 1
     !
     do jj = 1, nxp
        !
        if((xp(jj) .LE. x(i)) .AND. (xp(jj) .GE. x(i1))) then 
           xpi = xp(jj) - x(i1)
           xpi1 = xp(jj) - x(i)
           !
           l1 = (xpi1 * lp1(i1))**2
           l2 = (xpi  * lp2(i1))**2
           !
           yp(jj) = y(i1) * (1.0d0 - (2.0d0 * lp1(i1) * xpi)) * l1 + &
                y(i) * (1.0d0 - (2.0d0 * lp2(i1) * xpi1)) * l2 + &
                fp2(i1) * xpi1 * l2 + fp1(i1) * xpi * l1
        end if
        !
     end do
     !
  end do
  !
  ! out of bounds values? -> Linear extrapolation
  !
  if((maxval(xp) .GT. maxval(x)) .OR. (minval(xp) .LT. minval(x))) then
     !
     a0 = fp1(1)
     b0 = y(1) - a0 * x(1)
     !
     a1 = (y(np) - y(np1)) / (x(np) - x(np1))
     b1 = y(np) - a1 * x(np)
     !
     do ii = 1, nxp
        if(xp(ii) .LT. x(1)) yp(ii) = a0 * xp(ii) + b0
        if(xp(ii) .GT. x(np)) yp(ii) = a1 * xp(ii) + b1
     end do
  end if
  !
  return
END FUNCTION intepf

!
! For debugging purposes:
!subroutine printvars(atmo, k)
!  Use Model_structure 
!  Implicit none
!  Type (Model), target :: Atmo
!  integer  :: k
!  print *, 'z', atmo%z_scale(k)
!  print *, 'tau', atmo%ltau_500(k)
!  print *, 'temp', atmo%temp(k)
!  print *, 'rho', atmo%rho(k)
!  print *, 'gas_p', atmo%gas_p(k)
!  print *, 'el_p', atmo%el_p(k)
!  print *, 'v_los', atmo%v_los(k)
!  print *, 'b_long', atmo%b_long(k)
!  print *, 'b_x', atmo%b_x(k)
!  print *, 'b_y', atmo%b_y(k)
!  print *, 'ne', atmo%ne(k)
!  print *, 'nH', atmo%nH(k)
!  print *, 'nhplus', atmo%nhplus(k)
!  print *, 'nhminus', atmo%nhminus(k)
!  print *, 'nh2', atmo%nh2(k)
!  print *, ' '
!END subroutine printvars


SUBROUTINE ipol_dscale(params, atmo, Line)
! Subroutine ipol_dscale: Adapted from Mats Carlsson MULTI_3D.
!      Computes an optimal depth scale for RT purposes based on gradients
!      in T, rho and tau_500. The atmospheric model is interpolated to the 
!      new grid.
!
! ---------------------------------------------------------------------------
! Modifications: 
!      2011-03-06: 
!              JdlCR: added routine to the source code.
!
!      2011-03-14:
!              JdlCR: departure coefficients are also interpolated
!                     onto the new grid of depth points.
!   
!      2012-07-05:
!              JdlCR: Hermitian interpolation can rarely produce 
!                     overshooting. Changed to Bezier-spline instead.
!----------------------------------------------------------------------------
!
  Use Param_structure
  Use Line_data_structure
  Use Model_structure 
  Implicit None
  Type (Parameters), target :: Params
  Type (Model), target :: Atmo
  Type (Line_data), dimension (Params%n_lines) :: Line
  integer :: k, k0, k1, nn, kk, ii, jj, sm2
  Integer, pointer :: ndep
  Integer, parameter :: sm = 5
  Real, parameter :: ttop = 50000.0, taumax = 2.0
  Real, dimension(params%n_points), target :: aind2, xx, zin, tin
  Real :: tdiv, taudiv, rdiv
  Real, dimension(:), pointer ::  x
  Real, pointer :: xp
  Real, dimension(params%n_points) :: f
  Logical :: go
  integer, save ::ncall = 0 
  !
  ndep => params%n_points
  !
  ! Compute bounds
  !
  do k = 1, ndep
     zin(k) = atmo%z_scale(k)
     tin(k) = atmo%ltau_500(k)
     xx(k) = k - 1.0
  end do
  !
  k0 = 1
  k1 = 1
  !
  do k = 2, ndep
     if((atmo%temp(k).GT.ttop).AND.(k0.EQ.(k-1))) k0 = k
     if(atmo%ltau_500(k).LT.taumax) k1 = k
  end do
  !
  !k0 = k0 + 1
  !
  if(k0.GE.k1) then 
     do k = 1, ndep
        print *, atmo%ltau_500(k), atmo%temp(k)
     end do
     RETURN
  end if
  !
  ! Compute gradients -> scale factor from multi_3d -> 1. / log10(1.1) = 24.1588
  !
  aind2(:) = 0.0d0
  !
  do k = (k0 + 1), k1
     tdiv =   abs(log10(atmo%temp(k))     - log10(atmo%temp(k-1))) * 24.158858
     rdiv =   abs(log10(atmo%rho(k))      - log10(atmo%rho(k-1))) * 24.158858
     taudiv =    (atmo%ltau_500(k) - atmo%ltau_500(k-1)) * 10.0
     aind2(k) = aind2(k-1) + max(max(tdiv,rdiv), taudiv)
  end do
  !
  ! Normalize index range
  !
  do k = 1, ndep
     aind2(k) = aind2(k) / aind2(k1) * (ndep - 1.0)
  end do
  !
  ! Smooth gradients?
  aind2(k0:k1) = smooth((k1 - k0 + 1), aind2(k0:k1), sm)
  !
  ! New height scale and interpolate atmosphere 
  !
  nn = k1 - k0 + 1
  do k = 1, nn
     f(k) = atmo%z_scale(k+k0-1)
  end do
  x => aind2(k0:k1)
  !
  ! Height
  Call bezier3(nn,x, f, ndep, xx, atmo%z_scale) 
  !
  ! Tau
  do k = 1, nn
     f(k) = atmo%ltau_500(k+k0-1)
  end do
  Call bezier3(nn,x, f, ndep, xx, Atmo%ltau_500)
 
  !
  ! Now we interpolate in tau (not aind2)
  nullify(x)
  x => tin(k0:k1)
  do k = 1, ndep
     xx(k) = atmo%ltau_500(k)
  end do
  !
  ! Temperature
  do k = 1, nn
     f(k) = log(atmo%temp(k+k0-1))
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%temp)
  atmo%temp(1:ndep) = exp(atmo%temp(1:ndep))
  !
  ! Check temperature
  do k = 1, ndep
     xp => xx(k)
     if(atmo%temp(k).LT.1.0d3) then
        kk = 2
        go = .TRUE.
        do while((kk.LE.nn).AND.(go))
           if(((xp - x(kk-1)) * &
                (xp - x(kK))).LE.0.0d0) then
              !
              atmo%temp(k) = f(kk-1) + (xp - x(kk-1)) / (x(kk) - x(kk-1)) * &
                   (f(kk) - f(kk-1))
              atmo%temp(k) = exp(atmo%temp(k))
              !
              go = .FALSE.
           end if
           kk = kk + 1 
        end do
     end if
     nullify(xp)
  end do
  !
  ! Electron pressure
  do k = 1, nn
     f(k) = log(atmo%el_p(k+k0-1))
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%el_p)
  atmo%el_p(1:ndep) = exp(atmo%el_p(1:ndep))
  !
  ! Gas pressure
  do k = 1, nn
     f(k) = log(atmo%gas_p(k+k0-1))
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%gas_p)
  atmo%gas_p(1:ndep) = exp(atmo%gas_p(1:ndep))
  !
  ! Gas density
  do k = 1, nn
     f(k) = log(atmo%rho(k+k0-1))
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%rho)
  atmo%rho(1:ndep) = exp(atmo%rho(1:ndep))
  !
  ! nH
  do k = 1, nn
     f(k) = log(atmo%nH(k+k0-1))
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%nH)
  atmo%nH(1:ndep) = exp(atmo%nH(1:ndep))
  !
  ! nHminus
  do k = 1, nn
     f(k) = log(atmo%nHminus(k+k0-1))
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%nHminus)
  atmo%nHminus(1:ndep) = exp(atmo%nHminus(1:ndep))
  !
  ! nHplus
  do k = 1, nn
     f(k) = log(atmo%nHplus(k+k0-1))
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%nHplus)
  atmo%nHplus(1:ndep) = exp(atmo%nHplus(1:ndep))
  !
  ! nH2
  do k = 1, nn
     f(k) = log(atmo%nH2(k+k0-1))
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%nH2)
  atmo%nH2(1:ndep) = exp(atmo%nH2(1:ndep))
  !
  ! nH2plus
  do k = 1, nn
     f(k) = log(atmo%nH2plus(k+k0-1))
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%nH2plus)
  atmo%nH2plus(1:ndep) = exp(atmo%nH2plus(1:ndep))
  !
  ! L.O.S. velocity
  do k = 1, nn
     f(k) = atmo%v_los(k+k0-1) * 1.e-5
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%v_los)
  atmo%v_los(1:ndep) =  atmo%v_los(1:ndep)* 1.0e5
  !
  ! Microturbulence
  do k = 1, nn
     f(k) = atmo%v_mic(k+k0-1) * 1.e-5
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%v_mic)
  atmo%v_mic(1:ndep) =  atmo%v_mic * 1.e5
  !
  ! b_long
  do k = 1, nn
     f(k) = atmo%b_long(k+k0-1)
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%b_long)
  !
  ! b_x
  do k = 1, nn
     f(k) = atmo%b_x(k+k0-1)
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%b_x)
  !
  ! b_y
  do k = 1, nn
     f(k) = atmo%b_y(k+k0-1)
  end do
  Call bezier3(nn,x, f, ndep, xx, atmo%b_y)
  ! 
  ! Departure coefficients 
  !
  do ii = 1, Params%n_lines
     !
     ! b_low
     do k = 1, nn
        f(k) = log(Line(ii)%b_low(k+k0-1))
     end do
     Call bezier3(nn,x, f, ndep, xx, line(ii)%b_low)
     line(ii)%b_low(1:ndep) = exp(line(ii)%b_low)
     !
     ! b_up
     do k = 1, nn
        f(k) = log(Line(ii)%b_up(k+k0-1))
     end do
     Call bezier3(nn,x, f, ndep, xx,line(ii)%b_up)
     line(ii)%b_up(1:ndep) = exp(line(ii)%b_up)
  end do
  !
  if(associated(ndep)) nullify(ndep)
  if(associated(x)) nullify(x)
  !
END SUBROUTINE ipol_dscale
!
! Note: Filling factor of component 2 (ff2) is always assumed to be ff2=1-ff1
!
Subroutine Forward(Params, Line, Region, Atmo, Syn_profile, Hydro)
  Implicit None
  Type (Parameters) :: Params
  Type (Line_data), dimension (Params%n_lines) :: Line
  Type (Region_data), dimension (Params%n_regions) :: Region
  Type (Model_2comp) :: Atmo
  Real, Dimension (Params%n_data) :: Syn_profile, syn2
  Logical :: Hydro
!
  Call Forward_1comp(Params, Line, Region, Atmo%Comp1, Syn_profile, Hydro)
!
  If (Params%TwoComp .and. Atmo%Comp1%ffactor .lt. 0.99) then
     Call Forward_1comp(Params,Line,Region,Atmo%Comp2,syn2,Hydro)
     Syn_profile=Atmo%Comp1%ffactor*Syn_profile + (1.-Atmo%Comp1%ffactor)*syn2
  End if
!
! Use macroturbulence and stray light factor from the first component
  If (Params%Skip_lambda .eq. 1) &
       Call Convolve_profile(Params, Region, Atmo%Comp1, Syn_profile) ! Macroturbulence
  If (Atmo%Comp1%stray .gt. 0.001) &
       Call Add_stray_light(Params, Atmo%Comp1, Region, Syn_profile) ! Add stray light &
!
  Return 
!
End Subroutine Forward
!
Subroutine Forward_1comp(Params, Line, Region, Atmo_in, Syn_profile, Hydro)
  Implicit None
  Type (Parameters) :: Params
  Type (Line_data), dimension (Params%n_lines) :: Line
  Type (Region_data), dimension (Params%n_regions) :: Region
  Type (NLTE_Atom), Save :: Atom
  Type (NLTE_input), Save :: NLTEInput
  Type (NLTE_variables), Save :: NLTE
  Type (Model) :: Atmo, Atmo_in, Saved, Atmo_pre
  Integer, Parameter :: NQuad=3 ! Number of points in quadrature (2 to 5)
  Logical :: Hydro, Incomplete, End
  Real, Dimension (Params%n_data) :: Syn_profile
  Real, Dimension (Params%n_points) :: Source_f, ltau_500_mu, term
  Real, Dimension (Params%n_points) :: Cont_op_5000, Cont_op, Cont_op_5000_2
  Real, Dimension (Params%n_points) :: ng_i, ng_j
  Real, Dimension (Params%n_lines, Params%n_points) :: Dldop, Damp, Line_op
  Real, Dimension (Params%n_points,4,4) :: Absorp_height
  Real, Dimension (:,:,:,:), Allocatable :: Absorp, TotAbsorp
  Real, Dimension (:,:), Allocatable :: Phi_I, Phi_Q, Phi_U, Phi_V
  Real, Dimension (:,:), Allocatable :: Psi_Q, Psi_U, Psi_V, TotEmis
  Real, Dimension (:), Allocatable :: Wave
  Integer, Dimension (:), Allocatable :: IndexWave
  Real, Dimension (10) :: Pp
  Real, Dimension (4) :: Stokes
  Real, Dimension (3) :: x, y
  Real, Dimension (1) :: imin1, DWave1
  Real, Dimension (NQuad) :: XMU, WMU
  Real :: reference_cont, nu, Wave_cm, Av_mol_weight
  Real :: mu, DWave, last_update, metal
  Real :: PH, PHminus, PHplus, PH2, PH2plus, n2P, Scat
  Integer :: npoints, idepth, iline, iwave, nwlengths, idata, ichoice, iunit, ntrans
  Integer :: iregion, i, j, nw, istart, iend, NMu, imu, itran, imin, irec, iostat
  Integer, Parameter :: nformalsolutions=8
  Integer, Dimension(nformalsolutions) :: nformal
  Logical, Save :: do_NLTE
  Logical :: NLTE_done
  Logical :: CheckNaN
  Character (len=256) :: String
  Real, Parameter :: Min_Pe=1e-4
!
  Call Time_routine('forward',.True.)

  npoints=Params%n_points  
  Atmo=Atmo_in ! Model assign operation

  At_abund(1:N_elements)=Atmo%Abundance(1:N_elements)
  If (At_abund(1) .lt. 11.99 .or. At_abund(1) .gt. 12.01) then
     Print *,"Something's wrong with the abundance scale. Aborting (in forward)"
     Stop
  End if
!  First re-interpolate tau grid according to de la Cruz Rodriguez
!   (adapterd from Carlsson'¡s MULTI_3D). Changes are not propagated outside
  If (Params%Reinterpolate .gt. 0) & 
       Call ipol_dscale(Params, Atmo, Line)  
! Check atmosphere for sanity
  Do idepth=1, npoints
     If (Atmo%El_p(idepth) .lt. Min_Pe) then
        Call Debug_log('In forward. Pe .lt. Min_Pe parameter. Clipping it',2)
        Atmo%El_p(idepth)=Min_Pe
     End if
  End do
!
  Atmo_pre=Atmo
!
  Cont_op_5000_2=-10 ! Initialize
  Debug_warningflags(flag_forward)=0
  Debug_errorflags(flag_forward)=0
  NLTE_done=.FALSE.
!
  Params%def_abund=0 ! Force to not use ANN. Need to revise this

  Do iline=1, Params%n_lines
     ! Complete line information
     Line(iline)%Atomic_weight=At_weight(Line(iline)%Atomic_number)
     If (Line(iline)%Atomic_number .gt. N_elements) then
        Print *,'Error with element identification in forward.f90'
        Stop
     End if
  End do
!
  Do iregion=1, Params%n_regions ! Set background opacity enhancement factors
     nfudge=nfudge+1
     If (nfudge .gt. 1000) then
        Print *,'maxnfudge too low in forward/background.f90. Need at least:',Params%n_regions
        Stop
     End if
     Fudge_start(nfudge)=Region(iregion)%First_wlength
     Fudge_end(nfudge)=Region(iregion)%First_wlength+ &
          (Region(iregion)%nwavelengths-1)*Region(iregion)%Wave_step
     Fudge(nfudge)=Region(iregion)%Opacity_enh
  End do

  ! Need to do NLTE calculation or are we in pure LTE mode?
  do_NLTE=.False.
  Do iline=1, Params%n_lines
     If (Line(iline)%NLTEtransition .ne. 0) do_NLTE=.True.
  End do
  !
  If (do_NLTE) then
     NLTEInput%Hydro=Hydro
     NLTE%Atmo=Atmo
     Call NLTE_init(Params, NLTEinput, NLTE, Atom, Atmo)
     Call Read_NLTE_lines(Params, NLTEInput, NLTE, Atom, Line)
  End if
  If (.not. Allocated(Line(1)%b_low) ) then
     Do iline=1, Params%n_lines
        Allocate (Line(iline)%b_low(Params%n_points))
        Line(iline)%b_low(:)=1.
        Allocate (Line(iline)%b_up(Params%n_points))
        Line(iline)%b_up(:)=1.
     End do
  End if
  If (do_NLTE) then
     NLTEInput%Hydro=Hydro
     NLTE%Atmo=Atmo
     Call NLTE_init(Params, NLTEinput, NLTE, Atom, Atmo)
  End if

! If flux computation is required, set Gaussian weights
  If (Params%heliocentric .ne. 0) then
     NMu=1
     XMU(1)=Params%heliocentric
     WMU(1)=1.
  Else
     If (NQuad .gt. 5) then
        Print *,'Too many quadrature points'
        Stop
     End if
     NMu=NQuad
     XMU(1:NQuad)=XMU_Gauss_Quad(NQuad,1:Nquad)
     WMU(1:NQuad)=WMU_Gauss_Quad(NQuad,1:Nquad)
  End if
! Now let's put the model in hydrostatic equilibrium and fill in
! the columns Gas_P, El_p, Rho and Z_scale of the model.
!
  Saved=Atmo
  If (Hydro) then     
     Call Hydrostatic(Params, Atmo)
  Else ! Trust El_P, rho and Gas_P and fill the nH, nHplus, nHminus, 
       !       nH2 columns of the model

!     Call Compute_others_from_T_Pe_Pg(Params%n_points, Atmo%Temp, Atmo%El_p, Atmo%Gas_p, &
!          Atmo%nH, Atmo%nHminus, Atmo%nHplus, Atmo%nH2, Atmo%nH2plus)
     Atmo%ne=Atmo%el_p/bk/Atmo%temp
     If (Atmo%Keep_El_p .gt. 0.9) Atmo%El_p=Saved%El_p
     If (Atmo%Keep_Gas_p .gt. 0.9) Atmo%Gas_p=Saved%Gas_p
     If (Atmo%Keep_Rho .gt. 0.9) Atmo%Rho=Saved%Rho
     If (Atmo%Keep_nH .gt. 0.9) Atmo%nH=Saved%nH
     If (Atmo%Keep_nHminus .gt. 0.9) Atmo%nHminus=Saved%nHminus
     If (Atmo%Keep_nHplus .gt. 0.9) Atmo%nHplus=Saved%nHplus
     If (Atmo%Keep_nH2 .gt. 0.9) Atmo%nH2=Saved%nH2
     If (Atmo%Keep_nH2plus .gt. 0.9) Atmo%nH2plus=Saved%nH2plus
!
  End if
!

  metal=at_abund(26)-7.5
  Call compute_others_from_T_Pe_Pg(Params%n_points,Atmo%Temp, Atmo%el_p, Atmo%Gas_p,&
       Atmo%nH,Atmo%nHminus,Atmo%nHplus,Atmo%nH2,Atmo%nH2plus)
    Call Reset_densities(Atmo, Saved)
!

  nformal(1:nformalsolutions)=0. ! Initialize the f. s. counter
  idata=1 ! Index for Syn_profile
  Syn_profile(:)=-1e10
!
! Calculate line- and wavelength-independent data.
!
  Do idepth=1, npoints !$$ PARALLEL LOOP START 
     n2P=BK*Atmo%Temp(idepth)
     Cont_op_5000(idepth)=Background_opacity(Atmo%Temp(idepth), Atmo%El_p(idepth), Atmo%Gas_p(idepth), &
     Atmo%nH(idepth)*n2P, Atmo%nHminus(idepth)*n2P, Atmo%nHplus(idepth)*n2P, Atmo%nH2(idepth)*n2P, &
          Atmo%nH2plus(idepth)*n2P, 5000., Scat)
  End do !$$ PARALLEL LOOP ENDe
  Cont_op_5000=Cont_op_5000/Atmo%Rho ! Convert to cm^2/g  
!
! Loop in spectral lines. Compute line-dependent, wave-independent data.
!
  Do iline=1,Params%n_lines !$$ PARALLEL LOOP START 
     Dldop(iline, 1:npoints)=Line(iline)%Wlength/cc* & ! Doppler width in A
          Sqrt( 2.*bk*Atmo%Temp(1:npoints)/&
          At_weight(Line(iline)%Atomic_number)/mass_pr+ &
          Atmo%v_mic(1:npoints)**2. + Line(iline)%extra_vmic**2 )
     If (Line(iline)%NLTEtransition .eq. 0) then ! LTE line
        Do idepth=1, npoints
           Call Damping(Line(iline), Atmo%Temp(idepth), Atmo%El_p(idepth), &
                Atmo%Gas_p(idepth), Dldop(iline, idepth), &
                Damp(iline, idepth), -1.0, -1.0) ! Compute damping
        End do
        Call LTE_pop(Params, Line(iline), Atmo, ng_i, ng_j)
        If (Line(iline)%DepCoefMode .eq. 2) then ! Contents of b_low and b_up 
           ng_i(:)=Line(iline)%b_low         ! are actually n/g
           ng_j(:)=Line(iline)%b_up          ! rather than dep coefs
        End if
     Else ! NLTE
        itran=Line(iline)%NLTEtransition
        Do idepth=1, npoints
           Call Damping(Line(iline), Atmo%Temp(idepth), Atmo%El_p(idepth), &
                Atmo%Gas_p(idepth), Dldop(iline, idepth), &
                Damp(iline, idepth), Atom%GA(itran), Atom%GQ(itran)) 
        End do
        NLTE%Atmo=Atmo
        If (NLTEInput%VelFree) NLTE%Atmo%v_los=0.
        If (NLTE_done .eqv. .FALSE.) then
!           Call NLTE_init(Params, NLTEinput, NLTE, Atom, Atmo)
           NLTE%linear=Params%NLTE_Linear
           NLTEInput%UseColSwitch=0
           Call SolveStat(NLTE, NLTEInput, Atom)
           If (Debug_errorflags(flag_NLTE) .ge. 1) then ! Try again with different init
              If (Params%printout .ge. 3) Print *,'NLTE iteration did not converge. Trying again'
              NLTE%linear=1
              Call SolveStat(NLTE, NLTEInput, Atom)
           End if
           If (Debug_errorflags(flag_NLTE) .ge. 1) then ! Try again with different init
              If (Params%printout .ge. 3) Print *,'NLTE iteration did not converge. Trying again'
              NLTE%N=NLTE%NStar
              NLTEInput%NGacc=.False.
              NLTEInput%MaxIters=300
              NLTEInput%IStart=0
              Call SolveStat(NLTE, NLTEInput, Atom)
              ! Reset values
              NLTEInput%NGacc=.True.
              NLTEInput%MaxIters=500
              NLTEInput%IStart=1
           End if
           If (Debug_errorflags(flag_NLTE) .ge. 1) then ! Try again with different colswitch
              If (Params%printout .ge. 3) Print *,'Again NLTE iteration did not converge. One more try...'
              NLTE%N=NLTE%NStar
              NLTEInput%NGacc=.False.
              NLTEInput%MaxIters=300
              NLTEInput%UseColSwitch=1
              NLTEInput%IStart=0
              Call SolveStat(NLTE, NLTEInput, Atom)
              ! Reset values
              NLTEInput%NGacc=.True.
              NLTEInput%MaxIters=500
              NLTEInput%UseColSwitch=0
              NLTEInput%IStart=1
           End if
           If (Debug_errorflags(flag_NLTE) .ge. 1) then 
              If (Params%printout .ge. 3) Print *,'NLTE calculation: Giving up'
           End if
           NLTE_done=.TRUE.
        End if
        If (Debug_errorflags(flag_NLTE) .ge. 1) then
           Syn_profile(:)=1e10
           Debug_errorflags(flag_forward)=1
        End if
        i=Atom%i(itran)
        j=Atom%j(itran)
        If (Line(iline)%DepCoefMode .eq. 2) then   ! Contents of b_low and b_up 
           NLTE%N(i,:)=Line(iline)%b_low/Atom%g(i) ! are actually n/g
           NLTE%N(j,:)=Line(iline)%b_up/Atom%g(j)  ! rather than dep coefs
        End if
        ng_i(:)=NLTE%N(i,:)/Atom%g(i)*Line(iline)%NLTE_nl_ratio
        ng_j(:)=NLTE%N(j,:)/Atom%g(j)*Line(iline)%NLTE_nu_ratio
        Line(iline)%NLTESource_f(1:npoints,1:Line(iline)%NLTEgridsize)= &
             NLTE%Source_f(1:npoints,1:Atom%NQ(itran),itran)
     End if
     Line_op(iline, 1:npoints)=4.9947E-21*(10.**Line(iline)%loggf)* & 
          (ng_i(1:npoints)-ng_j(1:npoints))/Dldop(iline, 1:npoints)* &
          (Line(iline)%Wlength**2) ! Still needs the profile
     Line_op(iline, 1:npoints)=Line_op(iline, 1:npoints) / &
          Atmo%Rho(1:npoints) ! cm^2/cm^3 to cm^2/g
  End do ! End line loop

!
! Start loop in spectral regions and wavelengths.
!
  Do iregion=1, Params%n_regions
     last_update=-10 ! Make sure that opacities will be computed the first time
     nwlengths=Region(iregion)%nwavelengths
     If (Params%Reference_cont .eq. 0) then ! Do not normalize
        reference_cont=1.e15 ! For numerical precision. Will remove this factor
                             ! at the end
     Else If (Params%Reference_cont .eq. 1) then
!        reference_cont=conhsra(Region(iregion)%First_wlength + &
!             Region(iregion)%nwavelengths * &
!             Region(iregion)%Wave_step/2.) ! HSRA continuum at middlepoint
        reference_cont=syn_hsra(Region(iregion)%First_wlength + &
             Region(iregion)%nwavelengths * &
             Region(iregion)%Wave_step/2., 1.0)! HSRA continuum
     Else If (Params%Reference_cont .eq. 2) then
!        reference_cont=conhsra(5000.) ! HSRA continuum at 500nm
        reference_cont=syn_hsra(5000., 1.0)! HSRA continuum
     Else If (Params%Reference_cont .eq. 3) then
      reference_cont=syn_hsra(Region(iregion)%First_wlength + &
             Region(iregion)%nwavelengths * &
             Region(iregion)%Wave_step/2., Params%heliocentric)! HSRA continuum
     Else If (Params%Reference_cont .eq. 4) then
        reference_cont=1. ! Normalize to local continuum
     Else
        Print *,'Unknown spectrum normalization option'
        Stop
     End if

! Compute wavelength profiles
     If (Allocated(TotAbsorp)) then
        Deallocate (TotAbsorp)
        Deallocate (TotEmis)
        Deallocate (Absorp)
        Deallocate (Phi_I)
        Deallocate (Phi_Q)
        Deallocate (Phi_U)
        Deallocate (Phi_V)
        Deallocate (Psi_Q)
        Deallocate (Psi_U)
        Deallocate (Psi_V)
        Deallocate (Wave)
     End if
     Allocate (TotAbsorp(nwlengths, Params%n_points, 4, 4))
     Allocate (TotEmis(nwlengths, Params%n_points))
     Allocate (Absorp(nwlengths, Params%n_points, 4, 4))
     Allocate (Phi_I(nwlengths, Params%n_points))
     Allocate (Phi_Q(nwlengths, Params%n_points))
     Allocate (Phi_U(nwlengths, Params%n_points))
     Allocate (Phi_V(nwlengths, Params%n_points))
     Allocate (Psi_Q(nwlengths, Params%n_points))
     Allocate (Psi_U(nwlengths, Params%n_points))
     Allocate (Psi_V(nwlengths, Params%n_points))
     Allocate (Wave(nwlengths))
     TotAbsorp(:,:,:,:)=0.
     TotEmis(:,:)=0.
     Do iwave=1, nwlengths
        Wave(iwave)=Region(iregion)%First_wlength + &
             (iwave-1)*Region(iregion)%Wave_step ! Lambda in A
     End do

     Do iline=1, Params%n_lines ! Loop in lines to blend		
        If (MaxVal(Line_op(iline,:)/Cont_op_5000(:)) .gt. 1.e-1) then ! Consider line			
           Incomplete=.TRUE.
           Do While (Incomplete)
              istart= (Line(iline)%Wlength-Line(iline)%Width - Wave(1))/ &
                   Region(iregion)%Wave_step -1
              iend= (Line(iline)%Wlength+Line(iline)%Width - Wave(1))/ &
                   Region(iregion)%Wave_step +1
              If (istart .gt. nwlengths) istart=nwlengths-2
              If (istart .lt. 1) istart=1
              If (iend .lt. 1) iend=3
              If (iend .gt. nwlengths) iend=nwlengths
              If (iend .le. istart+2) then
                 istart=istart-2
                 iend=istart+2
              End if
              If (istart .lt. 1) istart=1
              If (iend .gt. nwlengths) iend=nwlengths
              nw=0
              Do iwave=istart, iend, Params%Skip_lambda
                 nw=nw+1
              End do
              If (Allocated(IndexWave)) Deallocate(IndexWave)
              Allocate (IndexWave(nw))
              nw=0
              Do iwave=istart, iend, Params%Skip_lambda
                 nw=nw+1
                 IndexWave(nw)=iwave
              End do
              Call Profiles(Params, Line(iline), nwlengths, nw, IndexWave, Wave(:), Atmo, &
                   Damp(iline,:), Dldop(iline,:), Line_op(iline,:), &
                   Phi_I(:,:), Phi_Q(:,:), &
                   Phi_U(:,:), Phi_V(:,:), Psi_Q(:,:), &
                   Psi_U(:,:), Psi_V(:,:), &
                   Cont_op_5000(:), Incomplete) ! Compute profiles
              If (istart .eq. 1 .and. iend .eq. nwlengths) Incomplete=.FALSE.
              If (Line(iline)%Width .gt. 5.) Incomplete=.FALSE.
              If (Incomplete) Line(iline)%Width=Line(iline)%Width*2.
           End do ! While incomplete
           Call Abs_matrix(npoints, nwlengths,  Phi_I(:,:), Phi_Q(:,:), &
                Phi_U(:,:), Phi_V(:,:), &
                Psi_Q(:,:), Psi_U(:,:), Psi_V(:,:), &
                Cont_op_5000, Absorp(:,:,:,:)) ! Build normalized absorp. m.
           Do iwave=1, nwlengths, Params%skip_lambda
              DWave=Wave(iwave)-Line(iline)%Wlength
              If (Line(iline)%NLTEtransition .ge. 1) then
                 If (Dwave .lt. MinVal(Line(iline)%NLTEgrid)) DWave=Abs(DWave)
              Endif
              nu=cc/(Wave(iwave)*1e-8) ! s^-1
              Wave_cm=cc/nu
              If (Line(iline)%NLTEtransition .gt. 0) then
                 If (DWave .ge. MinVal(Line(iline)%NLTEgrid) .and. &
                   DWave .le. MaxVal(Line(iline)%NLTEgrid)) then ! Interpolate NLTE source function
                    Do idepth=1, Params%n_points
                       Dwave1(1)=Dwave
                       Call bezier3(Line(iline)%NLTEgridsize,&
                            Line(iline)%NLTEgrid(:), &
                            Line(iline)%NLTESource_f(idepth,:),1, &
                            DWave1,Source_f(idepth))
                    End do
! S has units of 2*hh*nu^3/cc^2, we need  2*hh*c^2/(lambda^5). 
                    Source_f(:)=Source_f(:)*(nu/cc)*nu
                    Term=Absorp(iwave,:,1,1)*Source_f(:)
                    If (Line(iline)%DepCoefMode .eq. 1) &
                         Term=Term*Line(iline)%b_up ! Apply departure coeff
                    TotEmis(iwave,:)=TotEmis(iwave,:)+Term(:)
                 Else ! LTE
                    Term(:)=Absorp(iwave,:,1,1)* &
                         2.*hh/Wave_cm*cc*cc/(Wave_cm**4)/  & 
                         (exp(hh*nu/bk/Atmo%Temp(:))-1.) ! Planck f (cgs)
                    If(Line(iline)%DepCoefMode .eq. 2) then
                       Term(:) = Absorp(iwave,:,1,1)* &
                            2.*hh/Wave_cm*cc*cc/(Wave_cm**4)/  & 
                            ((ng_i(:) / ng_j(:)) - 1.) ! From populations read
                    End if
                    If (Line(iline)%DepCoefMode .eq. 1) &
                         Term=Term*Line(iline)%b_up ! Apply departure coeff
                    TotEmis(iwave,:)=TotEmis(iwave,:)+Term(:)
                 End if
              Else ! LTE
                 Term(:)=Absorp(iwave,:,1,1)* &
                      2.*hh/Wave_cm*cc*cc/(Wave_cm**4)/  & 
                      (exp(hh*nu/bk/Atmo%Temp(:))-1.) ! Planck f (cgs)
                 If(Line(iline)%DepCoefMode .eq. 2) then
                    Term(:) = Absorp(iwave,:,1,1)* &
                         2.*hh/Wave_cm*cc*cc/(Wave_cm**4)/  & 
                         ((ng_i(:) / ng_j(:)) - 1.) ! From populations read
                 End if
                 If (Line(iline)%DepCoefMode .eq. 1) &
                      Term=Term*Line(iline)%b_up ! Apply departure coeff
                 TotEmis(iwave,:)=TotEmis(iwave,:)+Term(:)
              End if ! NLTE transition

              If (Debug_level .ge. 1) then
                 if (minval(absorp(iwave,:,1,1)) .lt. 0) then
                    Call Debug_Log('Negative absorption in routine forward',2)
                    Do idepth=1,npoints
                       If (Absorp(iwave,idepth,1,1) .lt. 0) then
                          Write (String,*) 'idepth, Ab, T, el_p, gas_p, ng_i, ng_j=', &
      idepth,absorp(iwave,idepth,1,1),atmo%temp(idepth),atmo%el_p(idepth),&
      atmo%gas_p(idepth),ng_i(idepth),ng_j(idepth),iwave 
                          Call Debug_Log('Negative absorption in routine forward',1)
                       End if
                    End do
                 End if
              End if
           End do ! Do in iwave
  
           If (Line(iline)%DepCoefMode .eq. 1) then
              Do idepth=1, Params%n_points ! Apply depature coefficients
                 Absorp(:,idepth,:,:)=Absorp(:,idepth,:,:)* &
                      Line(iline)%b_low(idepth)
              End do
           End if

           TotAbsorp(:,:,:,:)=TotAbsorp(:,:,:,:)+ &
                Absorp(:,:,:,:) ! Total absorption of all blends (normalized)
        End if ! Consider line
     End do ! End loop in lines to blend
!
     Do iwave=1, nwlengths, Params%skip_lambda !$$ PARALLEL LOOP START
        If (abs(last_update-Wave(iwave)) .gt. Params%Update_opac) then !Update background opacities
           nu=cc/(Wave(iwave)*1e-8) ! s^-1
           Wave_cm=cc/nu
           Call Compute_others_from_T_Pe_Pg(Params%n_points, Atmo%Temp, Atmo%El_p, Atmo%Gas_p, &
                Atmo%nH, Atmo%nHminus, Atmo%nHplus, Atmo%nH2, Atmo%nH2plus)     
           Call Reset_densities(Atmo, Atmo_pre)
           Do idepth=1, npoints
              n2P=BK*Atmo%Temp(idepth)
              Cont_op(idepth)=Background_opacity(Atmo%Temp(idepth), Atmo%El_p(idepth), &
                   Atmo%Gas_p(idepth), Atmo%nH(idepth)*n2P, Atmo%nHminus(idepth)*n2P, &
                   Atmo%nHplus(idepth)*n2P, Atmo%nH2(idepth)*n2P, Atmo%nH2plus(idepth)*n2p, &
                   Wave(iwave), Scat)
           End do
           If (Cont_op_5000_2(1) .lt. 0) then
              Do idepth=1, npoints
                 Cont_op_5000_2(idepth)=Background_opacity(Atmo%Temp(idepth), &
                      Atmo%El_p(idepth), Atmo%Gas_p(idepth),Atmo%nH(idepth)*n2P, Atmo%nHminus(idepth)*n2P, &
                      Atmo%nHplus(idepth)*n2P, Atmo%nH2(idepth)*n2P, Atmo%nH2plus(idepth)*n2P, 5000., Scat)
              End do
              Cont_op_5000_2=Cont_op_5000_2/Atmo%rho ! Convert to cm2/g
           End if
           Cont_op=Cont_op/Atmo%rho ! Convert to cm2/g
!           Cont_op=Cont_op/Cont_op_5000_2
           Last_update=Wave(iwave)
        End if
 ! Add continuum opacity and emission at this wlength
        Do i=1, 4
           TotAbsorp(iwave, 1:Params%n_points, i, i)= &
                TotAbsorp(iwave, 1:Params%n_points, i, i)+ &
                Cont_op(1:Params%n_points)/ & 
                Cont_op_5000(1:Params%n_points)
        End do
 !
        TotEmis(iwave,:)=TotEmis(iwave,:)+Cont_op(:)/Cont_op_5000(:)* &
              2.*hh/Wave_cm*cc*cc/(Wave_cm**4)/  & 
                      (exp(hh*nu/bk/Atmo%Temp(1:npoints))-1.) ! Planck f (cgs)
        Source_f=TotEmis(iwave,:)/TotAbsorp(iwave,:,1,1)
!
        Syn_profile(idata:idata+3)=0. 
        Absorp_height(:,:,:)=TotAbsorp(iwave,:,:,:)

        Do imu=1, NMu ! Loop in angles to compute flux
           mu=XMU(imu) ! mu=Cos. of heliocentric angle
           ltau_500_mu=Atmo%ltau_500 - Log10(mu) ! Optical depth along l.o.s.
           Call time_routine('formalsolution',.True.)
           Call formal_solution(Params%n_points, Params%formal_solution, &
                ltau_500_mu, Absorp_height, Source_f, &
                Stokes, ichoice) ! Formal solution
           Call time_routine('formalsolution',.False.)
           If (Params%reference_cont .eq. 4) then ! Normalize to local cont
              If (iwave .eq. 1) then ! First point
                 If (imu .eq. 1) reference_cont=0.
                 reference_cont=reference_cont+WMu(imu)*Stokes(1)

              End if
           End if
           Stokes(1)=Stokes(1)+Region(iregion)%Bias ! Add spectrally-flat bias
           Syn_profile(idata:idata+3)=Syn_profile(idata:idata+3) + &
                WMu(imu)*Stokes(1:4)
           If (CheckNaN(Sum(Syn_profile(idata:idata+3)))) then
              If (Debug_level .ge. 1) then
                 Write (Debug_FileUnit,*) '  ** Error in forward: Returning NaNs'
                 Write (Debug_FileUnit,*) '     Region, wavelength, mu:'
                 Write (Debug_FileUnit,*) iregion,iwave,imu
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
                 Write (Debug_FileUnit,*) '     ------   Background Opacity  ------'
                 Write (Debug_FileUnit,*) Cont_op
                 Write (Debug_FileUnit,*) '     ------   TotAbsorp (1,1)  ------'
                 Write (Debug_FileUnit,*) TotAbsorp(iwave,:,1,1)
                 Write (Debug_FileUnit,*) '     ------   TotAbsorp (1,4)  ------'
                 Write (Debug_FileUnit,*) TotAbsorp(iwave,:,1,4)
                 Write (Debug_FileUnit,*) '     ------   TotAbsorp (2,3)  ------'
                 Write (Debug_FileUnit,*) TotAbsorp(iwave,:,2,3)
                 Write (Debug_FileUnit,*) '     ------   Source_f  ------'
                 Write (Debug_FileUnit,*) Source_f
                 Write (Debug_FileUnit,*) '     ------   NLTE?  ------'
                 Write (Debug_FileUnit,*) NLTE_done
                 If (NLTE_done) then
                    Write (Debug_FileUnit,*) '     -----  N_low ----- '
                    Write (Debug_FileUnit,*) NLTE%N(i,:)
                    Write (Debug_FileUnit,*) '     -----  N_up  ----- '
                    Write (Debug_FileUnit,*) NLTE%N(j,:)
                 End if
                 Call Debug_Log('Error in Forward, discarding results',1)
                 Syn_profile(:)=1E10 ! Discard this result
              End if
           End if
        End do

        If (ichoice .gt. nformalsolutions) then
           Print *,'Error in formal solution. Unknown formal solution method'
           Stop
        End if
        nformal(ichoice)=nformal(ichoice)+1 ! How many solutions of each type
        Syn_profile(idata:idata+3)=Syn_profile(idata:idata+3)/reference_cont
        idata=idata+4*Params%Skip_lambda ! Update the data index
     End do !$$ PARALLEL LOOP END


     Deallocate (TotAbsorp)
     Deallocate (TotEmis)
     Deallocate (Absorp)
     Deallocate (Phi_I)
     Deallocate (Phi_Q)
     Deallocate (Phi_U)
     Deallocate (Phi_V)
     Deallocate (Psi_Q)
     Deallocate (Psi_U)
     Deallocate (Psi_V)
     Deallocate (Wave)
     If (Allocated(IndexWave)) Deallocate (IndexWave)
  End do !$$ PARALLEL LOOP END

  If (Params%printout .ge. 4 .and.  Params%formal_solution .eq. 0) then
     Print *,'Number of Stokes formal solutions using Delobezier: ',nformal(1)
     Print *,'Number of Stokes formal solutions using WPM: ',nformal(4)
  End if

  If (Params%Reference_cont .eq. 0) & ! Remove normalization factor
       Syn_profile=Syn_profile*1.e15 
!
  If (Debug_errorflags(flag_NLTE) .ge. 1) Debug_errorflags(flag_forward)=1
  If (Debug_warningflags(flag_NLTE) .ge. 1) Debug_warningflags(flag_forward)=1
! Additional outputs?
  If (Debug_OutputPop .or. Debug_OutputContOp .or. Debug_OutputNLTEsf) then
     Call Open_file(iunit,'debug.txt')
     Write (iunit,*) Debug_OutputPop,Debug_OutputContOp,Debug_OutputNLTEsf
     Write (iunit,*) Params%n_points, Params%n_data,ATOM%NLIN+ATOM%NCNT
     Write (iunit,*) do_NLTE
     If (do_NLTE) then
        Write (iunit,*) ATOM%NK
     End if
     Write (iunit,*) Atmo%ltau_500
     Close (iunit)
  End if
  If (Debug_OutputPop) then
     Call Open_file_direct (iunit,'Populations.dat',RealBytes*Params%n_points)
     If (do_NLTE) then
        irec=1
        Do i=1, ATOM%NK
           Call Write_direct(LittleEndian,iunit,irec,NLTE%N(i,:),Params%n_points,iostat)
           irec=irec+1
        End do
        Do i=1, ATOM%NK
           Call Write_direct(LittleEndian,iunit,irec,NLTE%NStar(i,:),Params%n_points,iostat)
           irec=irec+1
        End do
     Else
        Call Write_direct(LittleEndian,iunit,1,ng_i,Params%n_points,iostat)
        Call Write_direct(LittleEndian,iunit,2,ng_j,Params%n_points,iostat)
     End if
     Close (iunit)
  End if
  If (Debug_OutputContOp) then ! en cm^2/g
     Call Open_file_direct (iunit,'Cont_opacity.dat',RealBytes*Params%n_points)
     Call Write_direct(LittleEndian,iunit,1,Cont_op(:),Params%n_points,iostat)
     Call Write_direct(LittleEndian,iunit,2,Cont_op_5000(:)*Atmo%Rho,Params%n_points,iostat) ! cm^2/cm^3
     Close (iunit)
  End if
  If (Debug_OutputNLTEsf .and. do_NLTE) then
     Call Open_file_direct (iunit,'NLTE_sf.dat',RealBytes*Params%n_points)
     irec=1
     Do itran=1, ATOM%NLIN+ATOM%NCNT
        Do j=1, Atom%NQ(itran)
           Call Write_direct(LittleEndian,iunit,irec,NLTE%Source_f(:,j,itran),Params%n_points,iostat)
           irec=irec+1
        End do
     End do
     Close (iunit)
  End if
!
  Call Time_routine('forward',.False.)
  Return
!
End Subroutine Forward_1comp


Subroutine Reset_densities(Atmo, Atmo_pre)
! If the user has enabled the keep switches, reset densities to their
! original values (i.e., those in Atmo_pre model)
!
  Type (Model) :: Atmo, Atmo_pre
!

  If (Atmo%Keep_El_p .gt. 0.9) Atmo%El_p=Atmo_pre%El_p
  If (Atmo%Keep_Gas_p .gt. 0.9) Atmo%Gas_p=Atmo_pre%Gas_p
  If (Atmo%Keep_Rho .gt. 0.9) Atmo%Rho=Atmo_pre%Rho
  If (Atmo%Keep_nH .gt. 0.9) Atmo%nH=Atmo_pre%nH
  If (Atmo%Keep_nHminus .gt. 0.9) Atmo%nHminus=Atmo_pre%nHminus
  If (Atmo%Keep_nHplus .gt. 0.9) Atmo%nHplus=Atmo_pre%nHplus
  If (Atmo%Keep_nH2 .gt. 0.9) Atmo%nH2=Atmo_pre%nH2
  If (Atmo%Keep_nH2plus .gt. 0.9) Atmo%nH2plus=Atmo_pre%nH2plus
  Return
!
  End Subroutine Reset_densities

End Module Forward_module
