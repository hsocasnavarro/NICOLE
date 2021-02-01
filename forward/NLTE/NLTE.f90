! Requires
!   Model_structure
!   Params_structure
!     Phys_constants
!   Atomic_data
!   Line_data_structure
!   Eq_state
!   LTE
!   Background_opacity_module, Only: Background_opacity
!   Gauss_quad
!   Profiling
!
!   Plus all other modules in this file, namely:
!      NLTE_Atom_structure
!      NLTE_inp
!      VirtualFiles
!      NLTE_Vars

Module NLTE_Atom_structure
!
! Variable names are reminiscent of those used in Carlsson's MULTI:
! However, some MULTI features are not supported here for simplicity:
!    IW .eq. 1 is not supported
!    ITRAD .eq. 4 is not supported
  Integer, Parameter :: MaxNFreqs=3000, MaxNTERMS=100, MaxNLIN=1000
  Integer, Parameter :: MaxNColStrings=5000
!
  Type NLTE_Atom
!
! Variable definitions:
! element=2-character string with the element name (e.g., 'Ca')
! Z=Nuclear charge (atomic number)
! ion_stage=ionization stage of ground level (e.g., 2 for CaII)
! NK=Number of energy levels
! NLIN=Number of bound-bound transitions
! NCONT=Number of bound-free transitions in detail
! NFIX=Number of transitions treated with radiation temperatures
! Abund and AWgt=These values will not be read from the model atom file. 
!        They will be replaced with the contents of forward/atomic_data.f90
! 
     Character (Len=2) :: element
     Integer :: ion_stage, Z
     Integer :: NK, NLIN, NCNT, NFIX
     Real :: Abund, AWgt
!
! Use pointers for allocatable variables
!
! NK-sized vectors:
! cm_1=Energy level in cm^-1
! eV=Same in eV
! g=Statistical weight
! lvl_label=String with level label
! ion=Ionization stage the level belongs to
     Real, Dimension(:), pointer :: cm_1, eV, g
     Character (Len=256), Dimension(:), pointer :: lvl_label
     Integer, Dimension(:), pointer :: ion
!
! NLIN-sized vectors (bound-bound line transitions)
! Q0=Core-to-wing transition (in Doppler units, defined in NLTE_input)).
!       Wavelengths are spaced 
!       linearly in the core and logarithmically in the wings
! GA, GW, GS=Damping constants
! A=Einsteing A coefficient
     Real, Dimension(:), pointer :: Q0, GA, GW, GQ, A
!
! (NLIN+NCNT)-sized vectors (b-b and b-f radiative transitions):
! i=lower level index
! j=upper level index
! NQ=Number of pts for frequency quadrature
! QMAX=Furthest frequency (in Doppler units)
! f=Oscillator strength
! Alamb=Vacuum wavelength (in Angstroms)
     Integer, Dimension(:), pointer :: i, j, NQ
     Real, Dimension(:), pointer :: f, QMAX, Alamb
!
! NFIX-sized vectors (b-f fixed radiative transitions with radiation temp):
! iFIX=lower level index
! jFIX=upper level index
! ITrad=Photospheric-Chromospheric switch used in MULTI (can be 1, 2 or 3)
! IPho=1-continuum, 0-line
! A0=Cross-section at the limit (cm^-2)
! TRad=Radiation temperature (K)
! f=Oscillator strength
     Integer, Dimension(:), pointer :: iFIX, jFIX, ITRad, IPho
     Real, Dimension(:), pointer :: A0, TRad
! 
! Other arrays:
! Q(1:NQ(1:NLIN+NCNT),1:NLIN+NCNT): Wavelength grid (in Doppler units)
!          Note: If read from file ATOM via QMAX.lt.0 switch, then
!          units are in Angstroms, converted to Doppler in freqquad
! FRQ(0:NQ(1:NLIN+NCNT),1:NLIN+NCNT): Freq grid for continuum trans (s^-1)
! WQ(1:NQ(1:NLIN+NCNT),1:NLIN+NCNT): Wavelength quadrature weights
! ALPHAC(1:NQ(1:NLIN+NCNT),1:NLIN+NCNT): Absorption profile
! WMU
! B(1:NK, 1:NK)=Einstein B coefficients (note the different dimensions wrt A)
! krad(1:NK,1:NK)=Index to transition between each pair of levels
     Real, Dimension(:,:), pointer :: Q, WQ, AlphaC, B, FRQ
     Integer, Dimension(:,:), pointer :: krad
!
! Information for collisional rates
! This is stored as strings that will be decoded by the corresponding routine
!
! ColRoutine=6-character array with the name of the routine to use.
!             Note: Currently, only CA2COL and GENCOL are implemented
! MaxNColStrings=Parameter setting the maximum number of lines with col data
! NColStrings=Actual of lines with collisional data
! ColStr(1:NColStrings)=Character*256 strings containing collisional data
!
     Character (Len=6) :: ColRoutine
     Integer :: NColStrings
     Character (Len=500), Dimension(MaxNColStrings) :: ColStr
! For multiple terms
     Real, Dimension(MaxNTERMS) :: DETERM, WTERM, FTERM, GATERM, GWTERM, GQTERM
     Integer, Dimension(MaxNLIN) :: NTERM
     Integer, Dimension(MaxNTERMS,MaxNLIN) :: KTERM
  End Type NLTE_Atom
!
End Module NLTE_Atom_structure
!
Module NLTE_inp
  Type NLTE_input
!
! Elim2=When to stop iterating (e.g., 1e-3)
! Elim1=When to stop updating the rate matrix E (e.g., same as Elim2)
! QNORM=Typical Doppler width (in km/s) used to define frequency 
!                                units (e.g., 8)
! NMU=Number of points to use in angular quadrature (e.g., 3)
! VelFree=Whether or not to use the Velocity-free approximation
! CPER=Artificial enhancment of collisional rates. Should be 1. Setting it 
!          to a large value results in LTE populations. May be useful for
!          debugging or investigating issues.
! UseColSwitch=If 1, use Colisional-Radiative Switching (Default 0)
! VirtualFiles=Large arrays containing opacities, profile and radiation field
!          information may be stored in memory (virtual files) or disk.
!          Using virtual files is much faster, but some large computations
!          may demand disk files.
! NGacc=Whether or not to use NG acceleration (e.g., .TRUE.)
! OptThin, OptThick=Minimum and maximum optical thickness for the
!          radiative transfer. Use 0 and Inf to do the transfer over
!          the entire atmo. Otherwise the atmo will be assumed to be
!          transparent above OptThin and use the diffussion approximation
!          below OptThick
! Linear=Use Linear interpolation in the NLTE formal solution (more stable
!          but less accurate)
! MaxIters=Maximum number of allowed iterations in solvestat (e.g.,200)
! IStart=As in MULTI, initialization for the populations: 0=Start with
!          zero radiation field and solve the statistical equilibrium eqs;
!          1=Start with LTE populations
! ISum=As in MULTI, number of equation to replace with particle conservation 
!          equation. If zero, pick the level with the largest population.
!          I've found that 1 usually yields best results.
! Verbose=Verbosity level
! Write=Write pop.dat and ltepop.dat files containing atomic populations
! Hydro=Put model in hydrostatic equilibrium to compute gas pressure and
!          density from electron pressure. Otherwise, trust the gas_p and
!          rho supplied as input.
! NumLambdaIters=Number of Lambda iterations to be performed before the
!          accelerated iteration
! NLTE_Formal_solution: 1 for Bezier, 2 for Short Characteristics
! NLTE_LtePop: 1 -> Use normal routine; 2-> Use MULTI routine
!
     Real :: elim1, elim2, QNORM, CPER, OptThin, OptThick
     Integer :: NMU, MaxIters, Verbose, ISum, IStart, UseColSwitch, &
          NumLambdaIters
     Integer :: NLTE_formal_solution
     Logical :: VelFree, VirtualFiles, NGacc, Write, Hydro
  End Type NLTE_input
End Module NLTE_inp
!
Module VirtualFiles
  Use File_Operations
  Use Param_structure, only : Realbytes ! Just to get the RealBytes parameter
  Integer :: ipointopc, ipointphi ! Virtual file pointers. Initialized in
                                  !   OpenVirtualFile or RewindVirtualFile
  Integer :: NDEP, NRAD, MQ, NMU ! Initialized in NLTE_init
  Integer :: lJNY, lPHI, lOPC ! Unit numbers
  Real, Dimension(:), Allocatable :: DJNY,DPHI,DXCONT,DSCAT,DSC ! Storage
!
  Contains
    Subroutine OpenVirtualFile(Filename, Virtual)
      Implicit None
      Character (Len=*) :: Filename
      Character (Len=3) :: File
      Logical :: Virtual
      Integer :: IRECW, iunit
!                      
      File=Filename(1:3)
      IRECW=NDEP*NRAD*MQ
      If (File .ne. 'JNY' .and. File .ne. 'OPC' .and. File .ne. 'PHI') then
         Print *,'Trying to open wrong file:',Filename
         Print *,'in forward/NLTE/NLTE.f90, routine openvirtualfiles'
         Stop
      End if
      If (Virtual) then ! Memory storage
         If (File .eq. 'JNY') then
            lJNY=-10
            If (.not. Allocated(DJNY)) &
                 Allocate (DJNY(NDEP*NRAD*MQ))
         Else if (File .eq. 'OPC') then
            lOPC=-20
            ipointopc=1
            If (.not. allocated(DXCONT)) then
               Allocate (DXCONT(NDEP*NRAD*MQ))
               Allocate (DSCAT(NDEP*NRAD*MQ))
               ALLOCATE (DSC(NDEP*NRAD*MQ))
            End if
         Else if (File .eq. 'PHI') then
            lPHI=-30
            ipointphi=1
            If (.not. allocated(DPHI)) &
            Allocate (DPHI(NDEP*NRAD*MQ*NMU))
         End if
         Return
      Else
         If (File .eq. 'JNY') then ! Unformatted direct
            Call Open_file_direct(iunit, File, IRECW)
            lJNY=iunit
         Else if (File .eq. 'PHI') then ! Formatted sequential
            Call Open_file(iunit, File)
            lPHI=iunit
         Else if (File .eq. 'OPC') then ! Formatted sequential
            Call Open_file(iunit, File)
            lOPC=iunit
         End if
      End if
      Return
    End Subroutine OpenVirtualFile
!
    Subroutine CloseVirtualFile(iunit)
      Implicit None
      Character (len=3) :: File
      Integer :: iunit
      If (File .ne. 'JNY' .and. File .ne. 'OPC' .and. File .ne. 'PHI') then
         Print *,'Trying to close wrong file:',File
         Print *,'in forward/NLTE/NLTE.f90, routine closevirtualfile'
         Stop
      End if
      If (File .eq. 'JNY') then
         iunit=lJNY
      Else if (File .eq. 'OPC') then
         iunit=lOPC
      Else if (File .eq. 'PHI') then
         iunit=lPHI
      End if
      If (iunit .gt. 0) then ! Disk file
         Call Close_file (iunit)
      Else ! Memory file
         If (iunit .eq. -10 .and. Allocated(DJNY)) Deallocate (DJNY)
         If (iunit .eq. -20) then
            ipointopc=1
            If (Allocated(DXCONT)) then
               Deallocate (DXCONT)
               Deallocate (DSCAT)
               Deallocate (DSC)
            End if
         End If
         If (iunit .eq. -30) Then
            ipointphi=1
            If (Allocated(DPHI)) &
                 Deallocate (DPHI)
         End if
      End if
      Return
    End Subroutine CloseVirtualFile
!
    Subroutine RewindVirtualFile(File)
      Implicit None
      Character (len=3) :: File
      Integer :: iunit
      If (File .ne. 'JNY' .and. File .ne. 'OPC' .and. File .ne. 'PHI') then
         Print *,'Trying to rewind wrong file:',File
         Print *,'in forward/NLTE/NLTE.f90, routine rewindvirtualfile'
         Stop
      End if
      If (File .eq. 'JNY') then
         iunit=lJNY
      Else if (File .eq. 'OPC') then
         iunit=lOPC
      Else if (File .eq. 'PHI') then
         iunit=lPHI
      End if
      If (iunit .gt. 0) then
         Rewind iunit
      Else
         If (iunit .eq. -20) ipointopc=1
         If (iunit .eq. -30) ipointphi=1
      End if
      Return
    End Subroutine RewindVirtualFile
!
    Subroutine ReadX(XCONT, SCAT, SC)
      Implicit None
      Integer :: initial, final, ND, iunit
      Real, Dimension(*) :: XCONT, SCAT, SC

      iunit=lOPC
      If (iunit .gt. 0) then
         Read (iunit, *) XCONT(1:NDEP), SCAT(1:NDEP), SC(1:NDEP)
      Else
         initial=(ipointopc-1)*NDEP+1
         final=ipointopc*NDEP
         XCONT(1:NDEP)=DXCONT(initial:final)
         SCAT(1:NDEP)=DSCAT(initial:final)
         SC(1:NDEP)=DSC(initial:final)
         ipointopc=ipointopc+1
      End if
      Return
    End Subroutine ReadX
!
    Subroutine WriteX(XCONT, SCAT, SC)
      Implicit None
      Integer :: initial, final, ND, iunit
      Real, Dimension(*) :: XCONT, SCAT, SC

      iunit=lOPC
      If (iunit .gt. 0) then
         Write (iunit, *) XCONT(1:NDEP), SCAT(1:NDEP), SC(1:NDEP)
      Else
         initial=(ipointopc-1)*NDEP+1
         final=ipointopc*NDEP
         DXCONT(initial:final)=XCONT(1:NDEP)
         DSCAT(initial:final)=SCAT(1:NDEP)
         DSC(initial:final)=SC(1:NDEP)
         ipointopc=ipointopc+1
      End if
      Return
    End Subroutine WriteX
!
    Subroutine ReadJ(irec, JNY)
      Implicit None
      Integer :: initial, final, ND, irec, iunit
      Real, Dimension(*) :: JNY

      iunit=lJNY
      If (iunit .gt. 0) then
         Read (iunit, REC=irec) JNY(1:NDEP)
      Else
         initial=(irec-1)*NDEP+1
         final=irec*NDEP
         JNY(1:NDEP)=DJNY(initial:final)
      End if
      Return
    End Subroutine ReadJ
!
    Subroutine WriteJ(irec, JNY)
      Implicit None
      Integer :: initial, final, ND, irec, iunit
      Real, Dimension(*) :: JNY(1:NDEP)

      iunit=lJNY
      If (iunit .gt. 0) then
         Write (iunit, REC=irec) JNY
      Else
         initial=(irec-1)*NDEP+1
         final=irec*NDEP
         DJNY(initial:final)=JNY(1:NDEP)
      End if
      Return
    End Subroutine WriteJ
!
    Subroutine ReadP(PHI)
      Implicit None
      Integer :: initial, final, ND, iunit
      Real, Dimension(*) :: PHI

      iunit=lPHI
      If (iunit .gt. 0) then
         Read (iunit, *) PHI(1:NDEP)
      Else
         initial=(ipointphi-1)*NDEP+1
         final=ipointphi*NDEP
         PHI(1:NDEP)=DPHI(initial:final)
         ipointphi=ipointphi+1
      End if
      Return
    End Subroutine ReadP
!
    Subroutine WriteP(PHI)
      Implicit None
      Integer :: initial, final, ND, iunit
      Real, Dimension(*) :: PHI

      iunit=lPHI
      If (iunit .gt. 0) then
         Write (iunit, *) PHI(1:NDEP)
      Else
         initial=(ipointphi-1)*NDEP+1
         final=ipointphi*NDEP
         DPHI(initial:final)=PHI(1:NDEP)
         ipointphi=ipointphi+1
      End if
      Return
    End Subroutine WriteP
!
End Module VirtualFiles
!
Module NLTE_vars
!
! Structure containing physical variables needed for the NLTE computation
!
  Use Model_structure
  Type NLTE_variables
!
! NDEP is the number of depth-points in the NLTE grid
! Nstar(NK,NDEP) and N(NK,NDEP) are LTE and NLTE population number 
! densities in cm^-3
!
! W, F and C are the rate matrices (NK,NK,NDEP): F->Fixed radiative trans,
!            C->Collisional rates. 
!
! BPlanck: Planck function B_nu
! Xnorm(1:NDEP)=Background continuum opacity (g/cm^2)
! XMu(1:NMU), WMU(1:NMU)=Gaussian integration weights
!
! Wphi=1./Integral_{angle, freq} (Phi_{mu,nu}) d_mu d_nu (see voigprofs)
! Sl(1:NDEP,1:NLIN)=Line source function (cgs)
! Source_f(1:NDEP, 1:NQ, 1:NLIN+NCNT)=Total source function (cgs)
!
! Error is set to .true. if something has gone wrong in the NLTE calculation
!  (e.g., non-convergence, NaN results, negative populations or other 
!  oddities) 
     Integer :: NDEP, MQ, NMU, NRAD, Linear
     Real, Dimension(:,:), pointer :: Nstar, N, BPlanck, WPhi, Sl
     Real, Dimension(:,:,:), pointer :: W, F, C, Source_f
     Real, Dimension(:,:), pointer :: Source_l_f
     Real, Dimension(:), pointer :: XMu, WMu, Xnorm
     Real :: OptThin, OptThick
     Logical :: Error=.False.
     Type (Model) :: Atmo
!
  End Type NLTE_variables
!
End Module NLTE_vars
!
! And now the main NLTE module
!
Module NLTE_module
  Use NLTE_Atom_structure
  Use NLTE_inp
  Use VirtualFiles
  Use NLTE_vars
  Use Param_structure
  Use Line_data_structure
  Use Forward_support
  Use Gauss_Quad
  Use LTE
  Use Eq_state
  Use Forward_support
  Use Background_opacity_module, Only: Background_opacity
  Use Zeeman_splitting
  Use Bezier_math
  Use Debug_module
  Use File_Operations
  Use Profiling
  Character (len = 25) :: NLTE_ver='NICOLE NLTE v1.0'
Contains
!
! Initialize NLTE computation
! The model is set in hydrostatic equilibrium in order to fill in the 
! required arrays (such as ne, nH, rho, etc)
!
Subroutine NLTE_init(Params, NLTEinput, NLTE, Atom, Atmo, ForceInit)
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Atmo
  Type (NLTE_Atom) :: Atom
  Type (NLTE_input) :: NLTEInput
  Type (NLTE_variables) :: NLTE
  Integer :: ielem, idepth
  Real, Dimension (:), Allocatable :: tmp1
  Logical, Save :: FirstTime=.TRUE.
  Logical :: Error, ForceInit
  Real, Dimension (10) :: Pp
  !
  If (ForceInit) FirstTime=.True.
  !
  If (FirstTime) then ! All intializations
     Call Read_atom(Atom)
     Call Find_atomic_number
     Atom%Abund=At_abund(Atom%Z)
     Atom%AWgt=At_weight(Atom%Z)
!
     Call Set_hardwired_values ! Need to externalize hardwired values
!
     Call NLTE_allocations ! Allocate NLTE atmospheric arrays
     Call FreqQuad(Atom, NLTEInput) ! Frequency quadrature
     Call AngQuad(NLTE, NLTEInput) ! Angular quadrature
!  For Virtual Files
     NDEP=NLTE%NDEP
     NRAD=Atom%NLIN+Atom%NCNT
     MQ=MaxVal(Atom%NQ(1:NRAD))
     NMU=NLTEInput%NMu
     NLTE%MQ=MQ
     NLTE%NMU=NMU
     NLTE%NRAD=NRAD
     Call OpenVirtualFile('OPC',NLTEInput%VirtualFiles)
     Call OpenVirtualFile('PHI',NLTEInput%VirtualFiles)
     Call OpenVirtualFile('JNY',NLTEInput%VirtualFiles)
  End if

! From here on may be needed for repeated calls with different atmospheres


  If (NLTEInput%Hydro) then
     Call Hydrostatic(Params, NLTE%Atmo) ! To fill in some arrays in Atmo
  Else ! Trust El_P, rho and Gas_P and fill the nH, nHplus, nHminus, 
     !       nH2 columns of the model

     NLTE%Atmo%ne=NLTE%Atmo%el_p/BK/NLTE%Atmo%Temp
     Call Compute_others_from_T_Pe_Pg(NLTE%NDEP, NLTE%Atmo%Temp, &
          NLTE%Atmo%El_p, NLTE%Atmo%Gas_p, &
          NLTE%Atmo%nH, NLTE%Atmo%nHminus, NLTE%Atmo%nHplus, NLTE%Atmo%nH2, &
          NLTE%Atmo%nH2plus)
  End if

  If (Params%NLTE_ltepop .eq. 1) then  ! Compute LTE populations
     Call LTE_pop_2(NLTE, Atom)
  Else
     Call LTE_pop_3(NLTE, Atom) ! routine from MULTI
  End if

  If (FirstTime) & ! initialize LTE only first time or each inversion step
       NLTE%N(:,:)=NLTE%NStar(:,:) ! Populations initial guess
  FirstTime=.FALSE.
  Call Collisions(NLTE, Atom) ! Collisional rates
  NLTE%C=NLTE%C*NLTEInput%CPER ! Enhance collisions if CPER .gt. 1
  Call FixedRad(NLTE, Atom) ! Continuum transitions with radiation temp
!
  Call BackgroundOpac(NLTE, NLTEInput, Atom) ! Compute background opacities
  Call VoigtProfs(NLTE, NLTEInput, Atom) ! Compute Voigt profiles

!
  Return
!
    contains
      Subroutine NLTE_allocations
        If (.not. Associated(NLTE%C)) then
           Allocate (NLTE%C(Atom%NK, Atom%NK, NLTE%NDEP))
           NLTE%C(:,:,:)=0.
           Allocate (NLTE%F(Atom%NK, Atom%NK, NLTE%NDEP))
           NLTE%F(:,:,:)=0.
           Allocate (NLTE%W(Atom%NK, Atom%NK, NLTE%NDEP))
           NLTE%W(:,:,:)=0.
           Allocate (NLTE%Nstar(Atom%NK, NLTE%NDEP))
           Allocate (NLTE%N(Atom%NK, NLTE%NDEP))
           Allocate (NLTE%Xnorm(NLTE%NDEP))
           Allocate (NLTE%BPlanck(NLTE%NDEP, Atom%NLIN+Atom%NCNT))
           Allocate (NLTE%XMu(NLTEInput%NMU))
           Allocate (NLTE%WMu(NLTEInput%NMU))
           Allocate (NLTE%WPhi(NLTE%NDEP,Atom%NLIN+Atom%NCNT))
           Allocate (NLTE%Sl(NLTE%NDEP,Atom%NLIN+Atom%NCNT))
           Allocate (NLTE%Source_f(NLTE%NDEP,MaxNFreqs,Atom%NLIN+Atom%NCNT))
           Allocate (NLTE%Source_l_f(NLTE%NDEP,Atom%NLIN+Atom%NCNT))
        End if
!
        Return
      End Subroutine NLTE_allocations
!
      Subroutine Find_atomic_number
        Character (len=256) :: str1, str2
        str1=Atom%element
        Call Tolower(str1)
        Do ielem=1, N_elements
           str2=Atom_char(ielem)
           Call Tolower(str2)
           If (str1 .eq. str2) Atom%Z=ielem
        End do
        Return
      End Subroutine Find_atomic_number
!
      Subroutine Set_hardwired_values
! These should
! be consistent with those in routine forward after the
! reset values comment
        NLTE%NDEP=Params%n_points
        NLTEInput%Elim1=Params%NLTE_Elim1
           ! debug. Doesn't work if elim1 .ne. elim2
        NLTEInput%Elim2=NLTEInput%Elim1 !   Need to look into it
        NLTEInput%IStart=Params%NLTE_Istart
        NLTEInput%CPER=Params%NLTE_CPER
        NLTEInput%UseColSwitch=Params%NLTE_UseColSwitch
        NLTEInput%NMU=Params%NLTE_NMU
        NLTEInput%QNORM=Params%NLTE_QNORM
        NLTEInput%VelFree=Params%NLTE_VelFree
        NLTEInput%VirtualFiles=.True.
        NLTEInput%NGacc=Params%NLTE_NGAcc
        NLTE%Linear=Params%NLTE_Linear
        NLTE%OptThin=Params%NLTE_OptThin
        NLTE%OptThick=Params%NLTE_OptThick
        NLTEInput%MaxIters=Params%NLTE_MaxIters
        NLTEInput%NLTE_formal_solution=Params%NLTE_Formal_solution
        NLTEInput%Verbose=Params%printout
        NLTEInput%ISum=Params%NLTE_Isum
        NLTEInput%NumLambdaIters=Params%NLTE_NumLambdaIters
        NLTEInput%Write=.False.
        NLTE%Atmo=Atmo ! Use the same grid
        If (NLTEInput%VelFree) NLTE%Atmo%v_los=0. ! Velocity free approximation
        Return
      End Subroutine Set_hardwired_values

End Subroutine NLTE_init
!
Subroutine Allocate_atom_arrays(Atom)

  Implicit None
  Type (NLTE_Atom) :: Atom
  Integer :: status
  !
  If (.not. Associated(Atom%cm_1)) then
     Allocate (Atom%cm_1(Atom%NK),stat=status)
     Allocate (Atom%eV(Atom%NK),stat=status)
     Allocate (Atom%g(Atom%NK),stat=status)
     Allocate (Atom%lvl_label(Atom%NK),stat=status)
     Allocate (Atom%ion(Atom%NK),stat=status)
     !
     Allocate (Atom%GA(Atom%NLIN),stat=status)
     Allocate (Atom%GW(Atom%NLIN),stat=status)
     Allocate (Atom%GQ(Atom%NLIN),stat=status)
     !
     Allocate (Atom%A(Atom%NLIN+Atom%NCNT),stat=status)
     Allocate (Atom%i(Atom%NLIN+Atom%NCNT),stat=status)
     Allocate (Atom%j(Atom%NLIN+Atom%NCNT),stat=status)
     Allocate (Atom%NQ(Atom%NLIN+Atom%NCNT),stat=status)
     Allocate (Atom%Q0(Atom%NLIN+Atom%NCNT),stat=status)
     Allocate (Atom%QMAX(Atom%NLIN+Atom%NCNT),stat=status)
     Allocate (Atom%f(Atom%NLIN+Atom%NCNT),stat=status)
     Allocate (Atom%Alamb(Atom%NLIN+Atom%NCNT),stat=status)
     !
     Allocate (Atom%iFIX(Atom%NFIX),stat=status)
     Allocate (Atom%jFIX(Atom%NFIX),stat=status)
     Allocate (Atom%ITRad(Atom%NFIX),stat=status)
     Allocate (Atom%IPho(Atom%NFIX),stat=status)
     Allocate (Atom%A0(Atom%NFIX),stat=status)
     Allocate (Atom%TRad(Atom%NFIX),stat=status)
     !
     Allocate (Atom%krad(Atom%NK,Atom%NK),stat=status)
     Allocate (Atom%B(Atom%NK,Atom%NK),stat=status)
     Allocate (Atom%FRQ(0:MaxNFreqs,Atom%NLIN+Atom%NCNT),stat=status)
     Allocate (Atom%Q(MaxNFreqs,Atom%NLIN+Atom%NCNT),stat=status)
     Allocate (Atom%WQ(MaxNFreqs,Atom%NLIN+Atom%NCNT),stat=status)
     Allocate (Atom%AlphaC(MaxNFreqs,Atom%NLIN+Atom%NCNT),stat=status)
  End if
!
  Return
!
End Subroutine Allocate_atom_arrays
!
! This subroutine reads the atomic file data from the ATOM file
! The format is similar to that of MULTI, except that comments
! are marked by the ! symbol, instead of the *
! There are some limitations as to the features supported here
! Check forward/NLTE/NLTE_atom_struct.f90 for details on these 
! limitations.
!
Subroutine Read_atom(Atom)
  Implicit None
  Type (NLTE_Atom) :: Atom
  Character (len=256) :: String
  Character (len=254) :: char_stage
  Integer :: fileunit, i, iw, lvlswap, i0, inu, itrm, KTRM
  Logical :: End, Error
  Logical, Save :: WarningIssued=.False.
!
  Call Open_file(fileunit, 'ATOM')
  Call Read_next_nocomment(fileunit, String, End)
  Call Tolower(String)
  String=AdjustL(String)
  Atom%element=String(1:2)
  char_stage=String(3:256)
  Atom%ion_stage=1 ! Default (if nothing else is found)
  If (SCAN(String,'0123456789ivx') .ne. 0) &
       Call Roman_to_Int(char_stage,Atom%ion_stage,Error) ! In case it's Roman
100  Call Read_next_nocomment(fileunit, String, End)
!  Read (String,*) Atom%Abund, Atom%AWgt ! These values will be overriden
  Call Read_next_nocomment(fileunit, String, End)
  Read (String,*) Atom%NK, Atom%NLIN, Atom%NCNT, Atom%NFIX
  Call Allocate_Atom_arrays(Atom)
  Do i=1, Atom%NK ! Read level information
     Call Read_next_nocomment(fileunit, String, End)
     Read (String,*) Atom%cm_1(i), Atom%g(i), Atom%lvl_label(i), Atom%ion(i)
  End Do
  Atom%eV(:)=Atom%cm_1(:)*cc*hh/ee
  KTRM=0
  Do i=1, Atom%NLIN ! Read b-b transition information
     Call Read_next_nocomment(fileunit, String, End)
     Read (String,*) Atom%j(i), Atom%i(i), Atom%f(i), Atom%NQ(i), &
          Atom%QMAX(i), Atom%Q0(i), iw, Atom%GA(i), Atom%GW(i), Atom%GQ(i)
     If (Atom%NQ(i) .gt. MaxNFreqs) then
        Print *,'Line transition ',i,' has too many frquencies:',Atom%NQ(i)
        Print *,'MaxNFreqs=',MaxNFreqs,' in NLTE_Atom_structure'
        Stop
     End if
     If (Atom%QMax(i) .lt. 0) then
        Read (fileunit, *) (Atom%Q(inu, i),inu=1,Atom%NQ(i))
     End if
     If (iw .ge. 2) then
        If (.not. WarningIssued) then
           Print *,'Warning. Reading multiple term transition in model atom, line transition:',i
           Print *,"This doesn't work exactly as in MULTI because the term structure is ignored"
           Print *,' for the damping calculation. Should be viewed as an approximation'
           WarningIssued=.True.
        End if
        
        Do itrm=1, iw
           KTRM=KTRM+1
           If (KTRM .gt. MaxNTERMS) then
              Print *,'NTERM too large ',KTRM,' compared to maximum dimension'
              Print *,'in NLTE_Atom_structure MaxNTERMS=',MaxNTERMS
              Stop
           End if
           Call Read_next_nocomment(fileunit, String, End)
           Atom%KTERM(itrm,i)=KTRM
           Read (String,*) Atom%DETERM(KTRM),Atom%WTERM(KTRM),Atom%FTERM(KTRM), &
                Atom%GATERM(KTRM),Atom%GWTERM(KTRM),Atom%GQTERM(KTRM)
           
        End Do
     End if
     If (Atom%Q0(i) .lt. 0 .or. iw .eq. 1) then
        Print *,'MULTI feature not supported in model atom'
        Print *,'Line ransition ',i
        Stop
     End if
     If (Atom%j(i) .lt. Atom%i(i)) then ! Swap levels if i .gt. j
        lvlswap=Atom%j(i)
        Atom%j(i)=Atom%i(i)
        Atom%i(i)=lvlswap
     End if
     Atom%krad(Atom%i(i),Atom%j(i))=i
     Atom%krad(Atom%j(i),Atom%i(i))=i
     Atom%Alamb(i)=hh*cc/ee*1.e8/(Atom%eV(Atom%j(i))-Atom%eV(Atom%i(i)))
     Atom%A(i)=Atom%f(i)*6.671e15*Atom%g(Atom%i(i))/Atom%g(Atom%j(i))/&
          Atom%Alamb(i)/Atom%Alamb(i)
     Atom%B(Atom%j(i),Atom%i(i))=Atom%Alamb(i)**3/(2.*hh*cc*1.e24)*Atom%A(i)
     Atom%B(Atom%i(i),Atom%j(i))=Atom%g(Atom%j(i))/Atom%g(Atom%i(i))* &
          Atom%B(Atom%j(i),Atom%i(i))
  End Do
  Do i=1, Atom%NCNT ! Read b-f transition information
     i0=Atom%NLIN+i
     Call Read_next_nocomment(fileunit, String, End)
     Read (String,*) Atom%j(i0), Atom%i(i0), Atom%f(i0), Atom%NQ(i0), &
          Atom%QMAX(i0)
     If (Atom%NQ(i0) .gt. MaxNFreqs) then
        Print *,'Continuum transition ',i,' has too many frquencies'
        Print *,'MaxNFreqs=',MaxNFreqs,' in NLTE_Atom_structure'
        Stop
     End if
     If (Atom%QMAX(i0) .lt. 0) then
        Do inu=1, Atom%NQ(i0)
           Call Read_next_nocomment(fileunit, String, End)
           Read (String,*) Atom%Q(inu, i0), Atom%AlphaC(inu, i0)
        End do
        Do inu=2, Atom%NQ(i0)
           If (Atom%Q(inu, i0) .gt. Atom%Q(inu-1, i0)) then
              Print *,'Continuum transition: ',i
              Print *,'Wavelengths not decreasing'
              Stop
           End if
        End do
     End if
     If (Atom%j(i0) .lt. Atom%i(i0)) then ! Swap levels if i .gt. j
        lvlswap=Atom%j(i0)
        Atom%j(i0)=Atom%i(i0)
        Atom%i(i0)=lvlswap
     End if
     Atom%krad(Atom%i(i0),Atom%j(i0))=i0
     Atom%krad(Atom%j(i0),Atom%i(i0))=i0
     Atom%Alamb(i0)=hh*cc/ee*1.e8/(Atom%eV(Atom%j(i0))-Atom%eV(Atom%i(i0)))
     Atom%A(i0)=Atom%f(i0)*6.671e15*Atom%g(Atom%i(i0))/Atom%g(Atom%j(i0))/&
          Atom%Alamb(i0)/Atom%Alamb(i0)
  End do
  Do i=1, Atom%NFIX ! Read fixed (radiation temp) b-f transition information
     Call Read_next_nocomment(fileunit, String, End)
     Read (String,*) Atom%jFix(i), Atom%iFix(i), Atom%IPho(i), Atom%A0(i), &
          Atom%TRad(i), Atom%ITRad(i)
     If (Atom%jFix(i) .lt. Atom%iFix(i)) then ! Swap levels if i .gt. j
        lvlswap=Atom%jFix(i)
        Atom%jFix(i)=Atom%iFix(i)
        Atom%iFix(i)=lvlswap
     End if
     If (Atom%ITRad(i) .gt. 3) then
        Print *,'MULTI feature not supported in model atom'
        Print *,'Fixed ransition ',i
        Stop
     End if
  End Do
  If (MaxVal(Atom%NQ(:)) .gt. MaxNFreqs) then 
     Print *,'MaxNFreqs is defined in NLTE_Atom_structure as ',MaxNFreqs
     Print *,'Should be at least ',MaxVal(Atom%NQ(:))
     Stop
  End if
! Now store information for calculation of collisional rates
  Call Read_next_nocomment(fileunit, String, End)
  String=Trim(AdjustL(String))
  Atom%ColRoutine=String(1:6)
  Atom%NColStrings=0
  Do While (.not. End)
     Call Read_next_nocomment(fileunit, String, End)
     If (.not. End) then 
        Atom%NColStrings=Atom%NColStrings+1
        If (Atom%NColStrings .gt. MaxNColStrings) then 
           print *,'In NLTE_Atom_structure:'
           print *,'MaxNColStrings=',MaxNColStrings,' (too low)'
           Stop
        End if
        Atom%ColStr(Atom%NColStrings)=String
     End if
  End Do
  Call Close_file (fileunit)
  Return
!
End Subroutine Read_atom
!
Subroutine BackgroundOpac(NLTE, NLTEInput, Atom)
  Implicit None
  Type (NLTE_variables) :: NLTE
  Type (NLTE_input) :: NLTEInput
  Type (NLTE_Atom) :: Atom
  Type (Model) :: Atmo
  Real, Dimension (NLTE%NDEP) :: XCONT, SC, SCAT
  Real, Dimension (1) :: Dummy
  Real, Dimension (:), Allocatable, Save :: Cont_op_5000_2
  Logical, Save :: FirstTime=.True.
  Real :: Freq, P_el, T_el, Theta, P_g, n2P, ScatCoef, metal, lambda
  Integer :: ihtot, idepth, itran, inu, i0
!
! Calculate standard opacities at 500nm (will be used for normalization)
!

  If (FirstTime) then
     If (.not. Allocated(Cont_op_5000_2)) &
          Allocate(Cont_op_5000_2(NLTE%NDEP))
     Cont_op_5000_2(:)=-10
     FirstTime=.False.
  End if
  metal=NLTE%Atmo%Abundance(26)-7.5
  Atmo=NLTE%Atmo
  Do idepth=1, NLTE%NDEP
     n2P=BK*NLTE%Atmo%Temp(idepth)
     NLTE%Xnorm(idepth)=Background_opacity(Atmo%Temp(idepth), Atmo%El_p(idepth), &
          Atmo%Gas_p(idepth), Atmo%nH(idepth)*n2P, Atmo%nHminus(idepth)*n2P, Atmo%nHplus(idepth)*n2P, &
          Atmo%nH2(idepth)*n2P, Atmo%nH2plus(idepth)*n2P, 5000., Dummy(1))
  End do
!
! Calculate background opacities at central wavelength of each line transition
!
  Call RewindVirtualFile('OPC')
  Do itran=1, Atom%NLIN
     freq=2.99792458e18/Atom%Alamb(itran)
     Do idepth=1, NLTE%NDEP
        n2P=BK*Atmo%Temp(idepth)
        XCONT(idepth)=Background_opacity(Atmo%Temp(idepth), Atmo%El_p(idepth), Atmo%Gas_p(idepth),Atmo%nH(idepth)*n2P, &
             Atmo%nHminus(idepth)*n2P, Atmo%nHplus(idepth)*n2P, Atmo%nH2(idepth)*n2P, &
             Atmo%nH2plus(idepth)*n2P, Atom%Alamb(itran), Scat(idepth))

        If (Cont_op_5000_2(idepth) .lt. 0) then

           Cont_op_5000_2(idepth)=Background_opacity(Atmo%Temp(idepth), Atmo%El_p(idepth),&
                Atmo%Gas_p(idepth),Atmo%nH(idepth)*n2P,Atmo%nHminus(idepth)*n2P, &
                Atmo%nHplus(idepth)*n2P, Atmo%nH2(idepth)*n2P, Atmo%nH2plus(idepth)*n2P, &
                5000., Dummy(1))
        End if
!        XCONT(idepth)=XCONT(idepth)/Cont_op_5000_2(idepth)

        NLTE%BPlanck(idepth,itran)=Planck(freq,NLTE%Atmo%Temp(idepth))
        SC(idepth)=(XCONT(idepth)-Scat(idepth))/XCONT(idepth)* &
             NLTE%BPlanck(idepth,itran)
        SCAT(idepth)=SCAT(idepth)/XCONT(idepth)
        XCONT(idepth)=XCONT(idepth)/NLTE%Xnorm(idepth)
        if (xcont(idepth) .lt. 0) &
             Call Debug_Log('Error!! Negative NLTE continuum opaicity',1)
     End do
     Call WriteX(XCONT, SCAT, SC)
  End do
!
! Calculate background opacities at each wavelength of continuum transitions
!
  Do itran=Atom%NLIN+1, Atom%NLIN+Atom%NCNT
     Do inu=1, Atom%NQ(itran)
        freq=Atom%FRQ(inu, itran)
        lambda=2.99792458e18/freq
        Do idepth=1, NLTE%NDEP
           n2P=BK*Atmo%Temp(idepth)
           XCONT(idepth)=Background_opacity(Atmo%Temp(idepth), Atmo%El_p(idepth), Atmo%Gas_p(idepth), &
                Atmo%nH(idepth)*n2P,Atmo%nHminus(idepth)*n2P, Atmo%nHplus(idepth)*n2P, &
                Atmo%nH2(idepth)*n2P, Atmo%nH2plus(idepth)*n2P, lambda, Dummy(1))
!           XCONT(idepth)=XCONT(idepth)/Cont_op_5000_2(idepth)
           NLTE%BPlanck(idepth,itran)=Planck(freq,NLTE%Atmo%Temp(idepth))
           SC(idepth)=(XCONT(idepth)-Dummy(1))/XCONT(idepth)* &
                NLTE%BPlanck(idepth,itran)
           SCAT(idepth)=Dummy(1)/XCONT(idepth)
           XCONT(idepth)=XCONT(idepth)/NLTE%XNorm(idepth)
        End do
        Call WriteX(XCONT, SCAT, SC)
     End do
  End do
!
  Return
!
End Subroutine BackgroundOpac
! This subroutine takes a structure Atom and computes the collisional rates
!
Subroutine Collisions(NLTE, Atom)
  Implicit None
  Type (NLTE_Atom) :: Atom
  Type (NLTE_variables) :: NLTE
!
  If (Atom%ColRoutine .eq. 'CA2COL') then
     Call Ca2Col(Atom, NLTE)
  Else if (Atom%ColRoutine .eq. 'GENCOL') then
     Call GenCol(Atom, NLTE)
  Else if (Atom%ColRoutine .eq. 'HCOL') then
     Call HCol(Atom, NLTE)
  Else
     Print *,'Collisional routine not supported:',Atom%ColRoutine
     Stop
  End if

End Subroutine Collisions

Subroutine Ca2Col(Atom, NLTE)
! Adapted from MULTI
 Use Param_structure
 Use NLTE_Atom_structure
 Use NLTE_vars
 Implicit None
 Type (NLTE_Atom) :: Atom
 Type (NLTE_variables) :: NLTE
 Real, Dimension(:), Allocatable :: CE, CI
 Real, Dimension(:,:,:), Allocatable :: C
 Integer :: i, j, idepth, imax
 Real :: DC
!
 Allocate (C(Atom%NK, Atom%NK, NLTE%NDEP))
 Allocate (CE(Atom%NK-1))
 Allocate (CI(Atom%NK-1))
 C(:,:,:)=0.
!
! Start with de-excitation rates first
!
!
! Bound-bound collisions with electrons
!
 Do j=2, Atom%NK-1
    imax=j-1
    Read (Atom%ColStr(j-1),*) CE(1:imax)
    Do i=1, j-1
       Do idepth=1, NLTE%NDEP
          C(j,i,idepth)=CE(i)*Sqrt(5000./NLTE%Atmo%Temp(idepth))*NLTE%Atmo%ne(idepth)
          C(i,j,idepth)=C(j,i,idepth)*NLTE%Nstar(j,idepth)/NLTE%Nstar(i,idepth)
       End do
    End do
 End do
!
! Collisions with neutral hydrogen
! (Only for Ca 6-level atom, we consider leves 2-3 and 4-5)
!
 If (Atom%NK .eq. 6) then
    Do idepth=1, NLTE%NDEP
       DC=1.5E-9*(NLTE%Atmo%Temp(idepth)/1000.)**0.33*NLTE%Atmo%nH(idepth)
       C(4,5,idepth)=C(4,5,idepth)+DC
       C(2,3,idepth)=C(2,3,idepth)+DC
       C(5,4,idepth)=C(5,4,idepth)+DC*NLTE%Nstar(4,idepth)/NLTE%Nstar(5,idepth)
       C(3,2,idepth)=C(3,2,idepth)+DC*NLTE%Nstar(2,idepth)/NLTE%Nstar(3,idepth)
    End do
 Else
    Print *,'Collisions with H ignored for this atom'
 End if
!
! Bound-free rates
!
 If (Atom%ion(6) .eq. Atom%ion(1)+1) then
    Read (Atom%ColStr(Atom%NK-1),*) CI(1:Atom%NK-1)
    Do i=1, Atom%NK-1
       Do idepth=1,NLTE%NDEP
          C(i,Atom%NK,idepth)=CI(i)*NLTE%Atmo%ne(idepth)*SQRT(NLTE%Atmo%Temp(idepth))*EXP(-ee/bk/NLTE%Atmo%Temp(idepth)* &
               (Atom%eV(Atom%NK)-Atom%eV(i)))
          C(Atom%NK,i,idepth)=C(i,Atom%NK,idepth)*NLTE%Nstar(i,idepth)/NLTE%Nstar(Atom%NK,idepth)
       End do
    End do
 Else
    Print *,'forward/NLTE/collisions: This atom does not have the expected structure...'
    Stop
 Endif
!
 NLTE%C=C
 Deallocate (C,CE,CI)
!
 Return
End Subroutine Ca2Col

Subroutine GenCol(Atom, NLTE)
! Adapted from MULTI
! Note: Some features have been removed in order to port that routine here!!
!
 Use Param_structure
 Use NLTE_Atom_structure
 Use NLTE_vars
 Type (NLTE_Atom) :: Atom
 Type (NLTE_variables) :: NLTE
!  general routine for computing collisional rates
!  format is either given in terms of a collision strength omega
!         or                          a ce(te) rate
!
!  omegas are used for charged ions in general, since they remain roughl
!  constant with temperature.  rate is prop. to 1/sqrt(te) * exp(delta-e
!
!  ce(te) values are used largely for neutrals for same reason.
!  rate is prop. to sqrt(te) * exp(delta-e)
!
!  ci(te) values are used for ionization rates: these have a
!         sqrt(te) * exp(delta-e) dependence as for ce(te)
! 
!  originally coded by p.g. judge,  april 2nd, 1987
!  modified by p.g. judge, may 21st 1987
!    derivatives with respect to ne and t also returned
!    file read only once and values stored
!    reading of data performed by rcol
!
!:
!: gencl  07-01-29  new routine: (mats carlsson)
!:        gencol split into rcoll (reads data) and gencl (interpolation)
!:        to enable multiple calls to gencol without rereading atom file
!:        especially important for _3d version
!:        based on gencol for radyn
!:

      PARAMETER (MTGRD = 50, MSHELL = 5, MID=3000)
!                                                                       
!  LOCAL VARIABLES: UP TO MTGRD POINTS IN TEMPERATURE GRID ALLOWED      
!  THIS CAN BE INCREASED                                                
!                                                                       
      PARAMETER (x100=500.d0)
      DIMENSION TGRID (MTGRD) , CGRID (MTGRD),work(mtgrd,6)                  
      Real, Parameter :: EK=ee/bk
      Real, Dimension(:,:,:), Allocatable :: C
      Real, Dimension(:), Allocatable :: TEMP, PE, NE, CEA, CT
      Real, Dimension(:,:), Allocatable :: NH
      Integer :: unit
      Character (Len=256) :: ColString
      Character (Len=80) :: text
      Character (Len=20) :: Atomid
!  ALPRD IS TRUE IF INITIA HAS BEEN CALLED TO CALCULATE ALPHA COLLISIONS
      LOGICAL :: ALPRD, ctneg, FirstTime=.True.
      DATA ALPRD / .FALSE. /                                             
      SAVE ALPRD    
! From Block CGENCL
      Dimension ntmp(mid),il(mid),ih(mid),ielc(mid)
      Dimension cgrd(4,mtgrd*2,mid),tgrd(mtgrd*2,mid)
      Dimension key(mid)
      character (len=8) :: key
      Dimension coeff(mtgrd,mid),cdi(5,mshell,mid)
      Dimension ncoeff(mid)

!
 Allocate (C(Atom%NK, Atom%NK, NLTE%NDEP))
 Allocate (CEA(NLTE%NDEP))
 Allocate (CT(NLTE%NDEP))

 C(:,:,:)=0.
 CEA(:)=0.
 CT(:)=0.


 NK=Atom%NK
 NDEP=NLTE%NDEP
 Allocate (Temp(NDEP))
 Allocate (NE(NDEP))
 Allocate (PE(NDEP))
 Allocate (NH(6,NDEP))
 Temp(:)=NLTE%Atmo%Temp(:)
 Pe(:)=NLTE%Atmo%El_p(:)
 NE(:)=NLTE%Atmo%El_p(:)/bk/Temp(:)
 NH(:,:)=0. 

 NH(1,:)=NLTE%Atmo%nH(:) ! Neutral H
 NH(6,:)=NLTE%Atmo%nHplus(:) ! Ionized H

 Atomid=ATOM%element
! Call Open_file(unit,'__scratch.dat')
! Do ipcount=1, Atom%NColStrings
!    Write (unit,'(A256)') Atom%ColStr(ipcount)
! End do
! Call Close_file (unit)
! Open (unit,File='__scratch.dat')
 unit=-1
 ipcount=1

 If (FirstTime) then ! MULTI routine rcoll
!
!  assume initially that data are given independent of temperature:
!  i.e. that ntemp is originally set = 1
!
      ntemp=1
      nid=0
      iel=1
!
!  read the keyword 'key'
!
      do 1800 i=1,mid
!
!  read atom file and store values in arrays
!  key,ntmp,tgrd,il,ih
!  put in possibility of fixed rate
!  rate goes from il to ih, not necessarily from low to high
!  (il and ih not converted in gencol to ilo and ihi as for other
!  rates)
!
!        read(unit,'(a)',end=1998) text
        If (ipcount .gt. Atom%NColStrings) Goto 1998
        ColString=Atom%Colstr(ipcount)
        ipcount=ipcount+1
        Read(ColString,'(a)') text
        call lcase(text)
        k0=1
        call getwrd(text,k0,k1,k2)

        if(text(k1:k2).eq.'ar85-che+') then
          key(i)='ar85-chp'
        else if(k2-k1.le.7) then
          key(i)='       '
          key(i)(1:k2-k1+1)=text(k1:k2)
        else
          Print *,'rcoll: key more than 8 characters'
          Stop
        endif
        if (key(i) .eq. 'end     ') then
          goto 1999
        else if (key(i) .eq. 'temp    ') then
!          read(unit,*) ntemp,(tgrid(it),it=1,min(ntemp,mtgrd))
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
           Read(ColString,*) ntemp,(tgrid(it),it=1,min(ntemp,mtgrd))
          if(ntemp .gt. mtgrd) then
            Print *,' rcoll: ntemp.gt.mtgrd'
            Stop
          endif
          goto 1800
        else if (key(i)(1:3) .eq. 'ohm' &
        .or.     key(i) .eq. 'ce      ' &
        .or.     key(i) .eq. 'cp      ' &
        .or.     key(i) .eq. 'ch      ' &
        .or.     key(i) .eq. 'chi     ' &
        .or.     key(i) .eq. 'ci      ' &
        .or.     key(i) .eq. 'ch0     ' &
        .or.     key(i) .eq. 'ch+     ' &
        .or.     key(i) .eq. 'calp    ' &
        .or.     key(i) .eq. 'reco    ' &
        .or.     key(i) .eq. 'fixrat  ') then
          if(key(i).eq.'calp    ') then
            print *,'rcoll: calp not implemented in _3d'
            Stop
          endif
!          read(unit,*) il(i),ih(i),(cgrid(it),it=1,ntemp)
          ColString=Atom%Colstr(ipcount)
          ipcount=ipcount+1
          Read(ColString,*) il(i),ih(i),(cgrid(it),it=1,ntemp)
! 
          if(key(i) .eq. 'ch+     ' .or. key(i) .eq. 'ch0     ') then
            do jt=1,ntemp
              cgrid(jt)=log10(cgrid(jt))
            end do
          endif
! 
!  23-feb-1994 p.g.judge modifications end:
! 
!
!  set up spline coefficients
!
          if(ntemp.ge.4) then
            m=4
            x55=5.5d0
            call tautsp(tgrid,cgrid,ntemp,x55,work,tgrd(1,i), &
              cgrd(1,1,i),ntmp(i),m,iflag)
          else
!            if(ntemp.le.2) then
!              write(ljoblo,*) 'linear interpolation: ',key(i) &
!               ,  ' between levels', il(i), ih(i)
!            else
!              write(ljoblo,*) 'splin interpolation: ',key(i) &
!              ,  ' between levels', il(i), ih(i)
!            endif
            do j=1,ntemp
              tgrd(j,i)=tgrid(j)
              cgrd(1,j,i)=cgrid(j)
            enddo
            ntmp(i)=ntemp
          endif
!
!  additional key-words from philip judge
!
        else if (key(i) .eq. 'semi    ') then
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
           read(ColString,*) il(i),ih(i),coeff(1,i)
        else if (key(i) .eq. 'ltdr    ') then
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
          read(ColString,*) il(i),ih(i),(coeff(j,i),j=1,5)
        else if (key(i) .eq. 'corona  ') then
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
          read(ColString,*) il(i),ih(i)
        else if (key(i) .eq. 'ar85-rr ') then
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
           read(ColString,*) il(i),ih(i),(coeff(j,i),j=1,2)
        else if (key(i) .eq. 'ar85-cdi') then
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
           read(ColString,*) il(i),ih(i),ncoeff(i)
          if(ncoeff(i) .gt. mshell) then
             Print *,'gencol: ncdi .gt. mshell'
             Stop
          endif
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
          read(ColString,*) ((cdi(j,l,i),j=1,5),l=1,ncoeff(i))
        else if (key(i) .eq. 'ar85-cea') then
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
          read(ColString,*) il(i),ih(i),coeff(1,i)
        else if (key(i)(1:7) .eq. 'ar85-ch') then
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
          read(ColString,*) il(i),ih(i)
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
          read(ColString,*) (coeff(j,i),j=1,6)
        else if (key(i) .eq. 'shull82 ') then
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
          read(ColString,*) il(i),ih(i), &
            (coeff(j,i),j=1,8)
        else if (key(i) .eq. 'burgess ') then
           ColString=Atom%Colstr(ipcount)
           ipcount=ipcount+1
          read(ColString,*) il(i),ih(i),coeff(1,i)
!
        else if(key(i)(1:1).ne.' ')then
          Print *,'rcoll:  unknown keyword '//key(i)// &
                  ' in atomic data'
          Stop
        endif
  1800 continue
!
      Print *,'rcoll: nid.gt.mid'
      Stop
  1998 continue
      Print *,' rcoll: keyword ''end'' not found in collision file'
      Stop
  1999 continue
      nid=i-1

 End if ! End multi routine rcoll ! FirstTime
!
! Start multi gencl routine
!
      ctneg=.false.
!
!  assume initially that data are given independent of temperature:
!  i.e. that ntemp is originally set = 1
!
      ntemp=1
!
!  read the keyword 'key'
!
      do 800 i=1,nid
!
        if (key(i) .eq. 'end     ') then
          goto 999
        else if (key(i) .eq. 'temp    ') then
          goto 800
        else if (      key(i)(1:3) .ne. 'ohm' &
        .and.          key(i) .ne. 'ce      ' &
        .and.          key(i) .ne. 'cp      ' &
        .and.          key(i) .ne. 'ch      ' &
        .and.          key(i) .ne. 'calp    ' &
        .and.          key(i) .ne. 'reco    ' &
        .and.          key(i) .ne. 'ch0     ' &
        .and.          key(i) .ne. 'ch+     ' &
        .and.          key(i) .ne. 'ci      ' &
        .and.          key(i) .ne. 'chi     ' &
        .and.          key(i) .ne. 'fixrat  ' &
!
        .and.          key(i) .ne. 'semi    ' &
        .and.          key(i) .ne. 'ltdr    ' &
        .and.          key(i) .ne. 'corona  ' &
        .and.          key(i) .ne. 'shull82 ' &
        .and.          key(i) .ne. 'burgess ' &
        .and.          key(i)(1:4) .ne. 'ar85'  )then
          goto 800
        else if ( key(i)(1:3) .eq. 'ohm' &
        .or.      key(i) .eq. 'ce      ' &
        .or.      key(i) .eq. 'cp      ' &
        .or.      key(i) .eq. 'ch      ' &
        .or.      key(i) .eq. 'calp    ' &
        .or.      key(i) .eq. 'reco    ' &
        .or.      key(i) .eq. 'ch0     ' &
        .or.      key(i) .eq. 'ch+     ' &
        .or.      key(i) .eq. 'ci      ' &
        .or.      key(i) .eq. 'chi     ' &
        .or.      key(i) .eq. 'fixrat  ') then
!
!  make interpolation for all depths
!
          do k=1,ndep
            tgp=temp(k)
            tgp=min(tgp,tgrd(ntmp(i),i))
            tgp=max(tgp,tgrd(1,i))
            if(ntmp(i).ge.4) then
              ct(k)=ppvalu(tgrd(1,i),cgrd(1,1,i),ntmp(i),4,tgp,0, &
                           mflag)
            else
              do j=1,ntmp(i)
                tgrid(j)=tgrd(j,i)
                cgrid(j)=cgrd(1,j,i)
              enddo
              ct(k)=splin(temp(k),tgrid,cgrid,ntmp(i),ntmp(i))
            endif
            if(key(i).eq.'ch+     '.or.key(i).eq.'ch0     ') then
              ct(k)=10.**ct(k)
            endif
            if(ct(k).lt.0.0d0) then

! the following condition has been introduced to force the ct value
! to a minimum level (zero) in case of collisions with hi (a la drawin)
! with very low (but positive) cross-sections that may cause
! the spline to be slightly negative   

              ilo=min( il(i), ih(i) )
              ihi=max( il(i), ih(i) )
!              write(ljoblo,48) key(i),ilo,ihi,temp(k)
!   48         format('gencol: ct negative;key,ilo,ihi,t=',a,2i4,f15.0)
              if(key(i).eq.'chi     '.or.key(i).eq.'ch      ') then
                ct(k) = 0.
!                write(ljoblo,'(a)') 'gencol: ct set to zero'
              else
                ctneg=.true.
              endif
            endif
          end do
        endif
!
!  identify the upper and lower levels:
!
        ilo=min( il(i), ih(i) )
        ihi=max( il(i), ih(i) )
!
!  check to see that the levels are between 1 and nk
!
        if(ilo .lt. 1 .or. ihi .gt. nk) then
          print *,' gencol:  level index outside range (1,nk)'
          stop
        endif
!
!  omegas are given (+ve ions)
!
        if(key(i)(1:3) .eq. 'ohm') then
          do  k=1,ndep
            cdn = 8.63e-06 * ct(k) * ne(k) / ( atom%g(ihi)*sqrt(temp(k)) )
            cup = cdn * nlte%nstar(ihi,k) / nlte%nstar(ilo,k)
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
            c(ilo,ihi,k) = cup + c(ilo,ihi,k)
          end do

!
!  ce values are given (neutrals)
!
        else if (key(i) .eq. 'ce      ') then
          do k=1,ndep
            cdn = ne(k) * ct(k) * atom%g(ilo) * sqrt(temp(k)) / atom%g(ihi)
            cup= cdn * nlte%nstar(ihi,k) / nlte%nstar(ilo,k)
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
            c(ilo,ihi,k) = cup + c(ilo,ihi,k)
          end do
!
!  cp values are given (b-b collisions with protons)
!
        else if (key(i) .eq. 'cp      ') then
          do k=1,ndep
            cdn = nh(6,k) * ct(k)
            cup= cdn * nlte%nstar(ihi,k) / nlte%nstar(ilo,k)
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
            c(ilo,ihi,k) = cup + c(ilo,ihi,k)
          end do
!
!  ch values are given (collisions with neutral hydrogen)
!
        else if (key(i) .eq. 'ch      ') then
          do k=1,ndep
            cdn = nh(1,k) * ct(k)
            cup= cdn * nlte%nstar(ihi,k) / nlte%nstar(ilo,k)
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
            c(ilo,ihi,k) = cup + c(ilo,ihi,k)
          end do
!
!  added by aegp 27/06/01
!  chi values are given (b-f collisions with neutral hydrogen)
!
        else if (key(i) .eq. 'chi     ') then
          do k=1,ndep
            cup = nh(1,k) * ct(k)
            cdn= cup * nlte%nstar(ilo,k) / nlte%nstar(ihi,k) 
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
            c(ilo,ihi,k) = cup + c(ilo,ihi,k)
          end do
!
!  ci values are given
!
        else if (key(i) .eq. 'ci      ') then
          do k=1,ndep
            dekt=  (atom%ev(ihi)-atom%ev(ilo)) * ek / temp(k)
            cup = ne(k) * ct(k) * exp(-dekt) * sqrt(temp(k))
            cdn = cup * nlte%nstar(ilo,k) / nlte%nstar(ihi,k)
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
            c(ilo,ihi,k) = cup + c(ilo,ihi,k)
          end do
!
!  ch or ch+ values given:  charge transfer collisions with 
!            hydrogen atoms/ions  nb:  ordering important here
!  example input:
!
!  ch0
!  1  2  1.e-9 1.e-9 1.e-9
!
!  this assumes that the rate 1.e-9*nh(1,k) will be added to the
!  collision rate from level 1 to level 2
! 
        else if (key(i) .eq. 'ch0     ') then
          do k=1,ndep
            c(il(i),ih(i),k)= nh(1,k) * ct(k)  + c(il(i),ih(i),k)
          end do
        else if (key(i) .eq. 'ch+     ') then
          do k=1,ndep
            c(il(i),ih(i),k)= nh(6,k) * ct(k)  + c(il(i),ih(i),k)
          end do
!
!  calp values are given (collisions with alpha particles)
!  atom has to be he
!  initia is called for first rate to read populations
!  from rstrt
!
        else if(key(i).eq.'calp    ') then
          if(atomid(1:2).ne.'he') then
            print *,'gencol: atom has to be he for calpha'
            stop
          endif
          if(.not.alprd) then
             Print *,'Call to rdalp not implemented'
             Stop
          endif
          do k=1,ndep
            cdn = nlte%n(nk,k) * ct(k) 
            cup= cdn * nlte%nstar(ihi,k) / nlte%nstar(ilo,k)
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
            c(ilo,ihi,k) = cup + c(ilo,ihi,k)
          end do
!
!  recombination coefficients are given
!  c(j-i)=ne*ct
!
        else if(key(i).eq.'reco    ') then
          do k=1,ndep
            c(ihi,ilo,k)= ne(k)*ct(k) + c(ihi,ilo,k)
          end do
!
!  fixrat values are given
!  fixed rates from il to ih (not necessarily from low to high)
!  no reverse rate
!
        else if (key(i) .eq. 'fixrat  ') then
          do k=1,ndep
            c(il(i),ih(i),k)= c(il(i),ih(i),k) + ct(k)
          end do
        else if (key(i) .eq. 'semi    ') then
          fab=coeff(1,i)
          eupcm=atom%ev(ihi)*ee/cc/hh
          elocm=atom%ev(ilo)*ee/cc/hh
          fem = fab*atom%g(ilo)/atom%g(ihi)
          do k=1,ndep
            ct(k)=semic(atom%ion(ilo),eupcm,elocm,fem,temp(k),iii)
            cdn=ct(k)*ne(k)
            cup = cdn * nlte%nstar(ihi,k) / nlte%nstar(ilo,k)
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
            c(ilo,ihi,k) = cup + c(ilo,ihi,k)
          end do
!
!  low-temperature dielectronic recombination coeffs. are given
!  formula (9) of nussbaumer and storey (a+a paper ii)
!
        else if (key(i) .eq. 'ltdr    ') then
          altdr=coeff(1,i)
          bltdr=coeff(2,i)
          cltdr=coeff(3,i)
          dltdr=coeff(4,i)
          fltdr=coeff(5,i)
          do k=1,ndep
            t4ltdr=temp(k)/1.e4
            cdn = 1.d-12*ne(k)* (altdr/t4ltdr + bltdr + &
                   cltdr*t4ltdr +dltdr*t4ltdr*t4ltdr) * &
                   exp(-fltdr/t4ltdr)/t4ltdr/sqrt(t4ltdr)
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
          end do
        else if (key(i) .eq. 'corona  ') then
          do k=1,ndep
            call coronr(atomid(1:3),atom%ion(ilo),temp(k),pe(k),recrat,rirat)
            c(ihi,ilo,k) =  c(ihi,ilo,k) + recrat*ne(k)
            c(ilo,ihi,k) =  c(ilo,ihi,k) + rirat*ne(k)
          end do
        else if(key(i)(1:7) .eq. 'ar85-rr') then
          aradar=coeff(1,i)
          etaar=coeff(2,i)
          do k=1,ndep
            cdn= aradar*(temp(k)/1.e4)**etaar
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
          end do
        else if(key(i)(1:8) .eq. 'ar85-cdi') then
          do k=1,ndep
            cup=0.d0
            do j=1,ncoeff(i)
              xj=cdi(1,j,i)*ee/bk/temp(k)
              fac=exp(-xj)*sqrt(xj)
              fxj=fac*( cdi(2,j,i)+cdi(3,j,i)*(1.d0+xj) &
                 +(cdi(4,j,i)-xj*(cdi(2,j,i)+cdi(3,j,i) &
                       *(2.d0+xj)))*fone(xj)+ &
               cdi(5,j,i)*xj*ftwo(xj) )
              fac=6.69d-07/cdi(1,j,i)/sqrt(cdi(1,j,i))
!
              cup = cup+fac*fxj
            end do
            if(cup .lt. -1.d-30) then
               print *,'gencol: cdi-cup .lt.0'
               stop
            endif
            c(ilo,ihi,k) = cup*ne(k) + c(ilo,ihi,k)
            c(ihi,ilo,k) = cup*ne(k)*nlte%nstar(ilo,k)/nlte%nstar(ihi,k) &
               + c(ihi,ilo,k)
          end do
!
!  addition p. judge 24-jan-1994: end
!
!
!  22-feb-1994 p.g.judge modifications start:
!
        else if(key(i) .eq. 'ar85-cea') then
          ceafak=coeff(1,i)
          call ar85cea(ilo,ihi,cea,NDEP,Temp,ATOM)
          do k=1,ndep
            cea(k)=cea(k)*ceafak
            c(ilo,ihi,k) = cea(k)*ne(k) + c(ilo,ihi,k)
!
!  this is incorrect since population of upper level is the upper level
!  of the doubly excited state, **, i.e. this should be multiplied by
!  g(**)/g(ihi) exp(e(**) - e(ihi)).  a better approximation would seem to
!  be 0. times this to avoid problems. hence comment this line out
!
!          c(ihi,ilo,k) = cea(k)*ne(k)*nstar(ilo,k)/nstar(ihi,k)
!     *        + c(ihi,ilo,k)
!
!  22-feb-1994 p.g.judge modifications end:
!
          end do
        else if(key(i)(1:7) .eq. 'ar85-ch') then
          ar85t1=coeff(1,i)
          ar85t2=coeff(2,i)
          ar85a=coeff(3,i)
          ar85b=coeff(4,i)
          ar85c=coeff(5,i)
          ar85d=coeff(6,i)
! 
!   ar85-ch    charge transfer recombination with neutral hydrogen
! 
          if(key(i)(1:8).eq.'ar85-ch ') then
            do k=1,ndep
              if(temp(k).ge.ar85t1 .and. temp(k).le.ar85t2) then 
                t4=temp(k)/1.d4
                cup = ar85a * 1.d-9 * t4**ar85b * (1. + ar85c  &
                     * exp(ar85d*t4))*nh(1,k)           
                c(il(i),ih(i),k) = cup + c(il(i),ih(i),k)
              endif
            end do
! 
!   ar85-ch+   charge transfer with ionized hydrogen
! 
          else if(key(i)(1:8).eq.'ar85-ch+') then
            do k=1,ndep
              if(temp(k) .ge. ar85t1 .and. temp(k) .le. ar85t2) then 
                t4=temp(k)/1.d4
                cup = ar85a * 1.d-9 * t4**ar85b * exp(-ar85c*t4)  &
                  * exp(-ar85d*ee/bk/temp(k))*nh(6,k) 
                c(il(i),ih(i),k) = cup + c(il(i),ih(i),k)
              endif
            end do
!  
!   ar85-che   charge transfer with neutral helium
! 
          else if(key(i)(1:8) .eq. 'ar85-che') then
            do k=1,ndep
              if(temp(k) .ge. ar85t1 .and. temp(k) .le. ar85t2) then 
                t4=temp(k)/1.d4
                hpop=nh(1,k)+nh(2,k)+nh(3,k)+nh(4,k)+ &
                 nh(5,k)+nh(6,k)
                call hepop(temp(k),hpop,he1,he2,he3)
                cup = ar85a * 1.d-9 * t4**ar85b * (1. + ar85c  &
                     * exp(ar85d*t4))*he1 
                c(il(i),ih(i),k) = cup + c(il(i),ih(i),k)
              endif
            end do
! 
!   ar85-chp  charge transfer with ionized helium
! 
          else if(key(i)(1:8) .eq. 'ar85-chp') then 
            do k=1,ndep
              if(temp(k) .ge. ar85t1 .and. temp(k) .le. ar85t2) then 
                t4=temp(k)/1.d4
                hpop=nh(1,k)+nh(2,k)+nh(3,k)+nh(4,k)+ &
                 nh(5,k)+nh(6,k)
                call hepop(temp(k),hpop,he1,he2,he3)
                cup = ar85a * 1.d-9 * t4**ar85b * exp(-ar85c*t4)  &
                     * exp(-ar85d*ee/bk/temp(k))*he2 
                c(il(i),ih(i),k) = cup + c(il(i),ih(i),k)
              endif
            end do
          endif
        else if(key(i)(1:7) .eq. 'shull82') then
          acolsh=coeff(1,i)
          tcolsh=coeff(2,i)
          aradsh=coeff(3,i)
          xradsh=coeff(4,i)
          adish=coeff(5,i)
          bdish=coeff(6,i)
          t0sh=coeff(7,i)
          t1sh=coeff(8,i)
!
!  november 9, 1994, modifications p. judge, begin
!  idea is to reduce dielectronic recombination rate
!  owing to collisional ionization of rydberg states
!  following summers 1974, appleton lab rept.
!
          summrs=1.0d0
          do k=1,ndep
! to remove density dependence in dielectronic rate comment out
! the next line
            summrs=summers(ilo,ihi,ne(k),ATOM)
            cdn=aradsh*(temp(k)/1.d4)**(-xradsh)  +  summrs* &
              adish /temp(k)/sqrt(temp(k)) * exp(-t0sh/temp(k))*   &
              (1.d0+bdish * (exp(- t1sh/temp(k))))
!
!  november 9, 1994, modifications p. judge, end
!
            cdn=cdn*ne(k)
            cup=acolsh * sqrt(temp(k)) * exp( -tcolsh / temp(k))  &
               / (1.d0 + 0.1d0 * temp(k) / tcolsh) 
            cup=cup*ne(k)
!*
!*  3-body recombination (high density limit)
!*
            cdn=cdn+cup*nlte%nstar(ilo,k)/nlte%nstar(ihi,k)
!*
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
            c(ilo,ihi,k) = cup + c(ilo,ihi,k)
          end do
!
        else if(key(i) .eq. 'burgess ') then
          de=(atom%ev(ihi)-atom%ev(ilo))
          zc=atom%ion(ilo)-1
          betab = 0.25d0*(sqrt((100.d0*zc +91)/(4.d0*zc+3))-5.d0)
          cbar=2.3d0
          do k=1,ndep
            dekt=de*ee/bk/temp(k)
            ddekt=-dekt/temp(k)
            dekt=min(x100,dekt)
            dekti=1.d0/dekt
            wlog=log(1.d0 + dekti)
            wb = wlog**(betab/(1.d0+dekti))
!
            exin1=expint(1,dekt,dum)
            dexin1=-exp(-dekt)/dekt*ddekt
            qb = 2.1715d-08*cbar*(13.6d0/de)**1.5d0 * sqrt(dekt) * &
               exin1*wb
!  gencol 95-05-22 vh actelc=1 to be replaced later by
!                     routine that computes number
!                     of active electrons from term information
            actelc=1.d0
            cup=qb*actelc*ne(k)
            cdn=cup*nlte%nstar(ilo,k)/nlte%nstar(ihi,k)
            c(ilo,ihi,k) = cup + c(ilo,ihi,k)
            c(ihi,ilo,k) = cdn + c(ihi,ilo,k)
          end do
        else if(key(i)(1:6).eq.'radrat')then
           crat=coeff(1,i)
           do k=1,ndep
             c(il(i),ih(i),k)=crat+c(il(i),ih(i),k)
           end do
        endif
  800 continue
!
  999 continue

!        Close (unit,Status='DELETE')
      If (ctneg) then
         print *,'gencol: ct negative'
         Stop
      End if
!
      NLTE%C=C
!
      Deallocate (C, CEA, CT, Temp, NE, PE, NH)
      return
!
End Subroutine GenCol
!
SUBROUTINE LJUST (TEXT)                                            
!                                                                       
!  LEFT JUSTIFIES THE STRING TEXT                                       
!                                                                       
      CHARACTER ( * ) TEXT                                               
!                                                                       
      L = LEN (TEXT)                                                     
      DO 100 J = 1, L                                                    
         IF (TEXT (J:J) .NE.' ') GOTO 200                                
  100 END DO                                                             
  200 CONTINUE                                                           
      DO 300 I = 1, L                                                    
         IF (J.LE.L) THEN                                                
            TEXT (I:I) = TEXT (J:J)                                      
         ELSE                                                            
            TEXT (I:I) = ' '                                             
         ENDIF                                                           
         J = J + 1                                                       
  300 END DO                                                             
!                                                                       
      RETURN                                                             
    END SUBROUTINE LJUST
!                                                                       
!***********************************************************************
!                                                                       
      FUNCTION FONE (X)                                                  
!                                                                       
! NEW ROUTINE 24-JAN-1994                                               
! CALCULATES F1(X) NEEDED FOR COLLISIONAL RATES OF ARNAUD AND ROTHENFLUG
! MODIFIED BY P. JUDGE 15-MAR-1994                                      
! FOR LARGE VALUES OF X, USE ASYMPTOTIC LIMIT                           
!                                                                       
      DATA SPLIT / 50.0 /                                                
      IF (X.LE.SPLIT) THEN                                               
         FONE = EXP (X) * EXPINT (1, X, EX)                              
      ELSE                                                               
         FONE = 1. / X                                                   
      ENDIF                                                              
      RETURN                                                             
      END FUNCTION FONE                             
!
!*********************************************************************
!
      SUBROUTINE AR85CEA(ILO,IHI,CEA,NDEP,Temp,ATOM)
        Use Phys_constants
!
!  NEW ROUTINE FOR COMPUTING COLLISIONAL AUTOIONIZATION 
!  RATES USING FORMALISM AND FORMULAE FROM ARNAUD AND ROTHENFLUG 1985
!
!
! AR85CEA 94-02-22  NEW ROUTINE: (PHILIP JUDGE)
!        
      Type (NLTE_Atom) :: Atom
      REAL, DIMENSION(NDEP) :: CEA, Temp
      CHARACTER (len=2) :: CSEQ, ARELEM
      LOGICAL :: DEBUG, Defined
      INTEGER :: NDEP
      DATA DEBUG /.FALSE./

      Defined=.False.

!  INITIALIZE OUTPUT TO ZERO
      DO 2 K=1,NDEP
        CEA(K)=0.0
    2 CONTINUE
!
!
!  FIND ELEMENT
!
      ARELEM = Atom%Element
      Call ToUpper(ARELEM)
      IZ=IATOMN(ARELEM)   
      ZZ=0.0 + IZ
      IF (IZ .LT. 1 .OR. IZ .GT. 92) THEN 
        Print *,'AR85-CEA:  ATOMIC NUMBER = ',IZ, ' FOR ',ARELEM
        Print *,'AR85-CEA:  NO AUTOIONIZATION INCLUDED'
        Stop
      ELSE IF(DEBUG) THEN
        Print *,'AR85-CEA: ',ARELEM, ' CHARGE ',Atom%ION(ILO)-1
        Stop
      ENDIF
!
!  FIND ISO-ELECTRONIC SEQUENCE
!
      ICHRGE=Atom%ION(ILO)-1
      ISOSEQ=IZ-ICHRGE
      CSEQ=ATOMNM(ISOSEQ)  
      DO 1 K=1, NDEP
        BKT=BK*TEMP(K)/EE
!
!*********************************************************************
! LITHIUM SEQUENCE
!*********************************************************************
!
        IF (CSEQ .EQ. 'LI') THEN 
           Defined=.True.
          CEAB=1. / (1. + 2.E-4*ZZ*ZZ*ZZ)
          CEAZEFF=(ZZ-1.3)
          CEAIEA=13.6*((ZZ-0.835)*(ZZ-0.835) - 0.25*(ZZ-1.62)*&
           (ZZ-1.62))
          Y=CEAIEA/BKT
          F1Y=FONE(Y)
!PREPRINT          CUP= 8.0E+10 * CEAB /CEAZEFF/CEAZEFF/SQRT(BKT) 
!PREPRINT     *   * EXP(-Y)*(2.22*F1Y+0.67*(1.-Y*F1Y)+0.49*Y*F1Y+1.2*(Y-F1Y))
          CUP=1.60E-07*1.2*CEAB/CEAZEFF/CEAZEFF/SQRT(BKT)*EXP(-Y)* &
           (2.22*F1Y+0.67*(1.-Y*F1Y)+0.49*Y*F1Y+1.2*Y*(1.-Y*F1Y))
!
!  TWO SPECIAL CASES:
!
! C IV - APP A AR85
          IF(ARELEM .EQ. 'C ') CUP = CUP*0.6
! N V  - APP A AR85
          IF(ARELEM .EQ. 'N ') CUP = CUP*0.8
! O VI - APP A AR85
          IF(ARELEM .EQ. 'O ') CUP = CUP*1.25  
!
!*********************************************************************
! SODIUM SEQUENCE
!*********************************************************************
!
        ELSE IF (CSEQ .EQ. 'NA') THEN 
           Defined=.True.
          IF (IZ .LE. 16) THEN 
            CEAA=2.8E-17*(ZZ-11.)**(-0.7)
            CEAIEA=26.*(ZZ-10.)
            Y=CEAIEA/BKT
            F1Y=FONE(Y)
            CUP= 6.69E+7 * CEAA *CEAIEA/SQRT(BKT) * EXP(-Y)*(1. - Y*F1Y)
          ELSE IF (IZ .GE. 18 .AND.  IZ .LE. 28) THEN 
            CEAA=1.3E-14*(ZZ-10.)**(-3.73)
            CEAIEA=11.*(ZZ-10.)*SQRT(ZZ-10.)
            Y=CEAIEA/BKT
            F1Y=FONE(Y)
            CUP= 6.69E+7 * CEAA *CEAIEA/SQRT(BKT) * EXP(-Y)  &
               *(1. - 0.5*(Y -Y*Y + Y*Y*Y*F1Y))
          ELSE 
            RETURN
          ENDIF
        ENDIF
!
!*********************************************************************
! MAGNESIUM-SULFUR SEQUENCES
!*********************************************************************
!
        IF(CSEQ .EQ. 'MG' .OR. CSEQ .EQ. 'AL'  &
         .OR. CSEQ .EQ. 'SI' .OR.   &
         CSEQ .EQ. 'P ' .OR. CSEQ .EQ. 'S ') THEN 
          IF(CSEQ .EQ. 'MG') CEAIEA=10.3*(ZZ-10.)**1.52
          IF(CSEQ .EQ. 'AL') CEAIEA=18.0*(ZZ-11.)**1.33
          IF(CSEQ .EQ. 'SI') CEAIEA=18.4*(ZZ-12.)**1.36
          IF(CSEQ .EQ. 'P ' ) CEAIEA=23.7*(ZZ-13.)**1.29
          IF(CSEQ .EQ. 'S ' ) CEAIEA=40.1*(ZZ-14.)**1.1
          CEAA=4.0E-13/ZZ/ZZ/CEAIEA
          Y=CEAIEA/BKT
          F1Y=FONE(Y)
           Defined=.True.
          CUP= 6.69E+7 * CEAA *CEAIEA/SQRT(BKT) * EXP(-Y)   &
            *(1. - 0.5*(Y -Y*Y + Y*Y*Y*F1Y))
        ENDIF
!
!
!*********************************************************************
!  SPECIAL CASES
!*********************************************************************
!
!
        IF(ARELEM .EQ. 'CA' .AND.  (ICHRGE .EQ. 0)) THEN 
          CEAA = 9.8E-17
          CEAIEA= 25.
          CEAB=1.12         
          Y=CEAIEA/BKT
          F1Y=FONE(Y)
           Defined=.True.
          CUP=6.69E+7*CEAA *CEAIEA/SQRT(BKT) * EXP(-Y)*(1. + CEAB*F1Y)
!          WRITE(LJOBLO,*)'AR85-CEA SPECIAL CASE ',ARELEM,  &
!           ' ION ICHRGE ',ICHRGE
        ELSE IF(ARELEM .EQ. 'CA' .AND.  (ICHRGE .EQ. 1)) THEN 
          CEAA = 6.0E-17
          CEAIEA= 25.
          CEAB=1.12         
          Y=CEAIEA/BKT
          F1Y=FONE(Y)
           Defined=.True.
          CUP=6.69E+7*CEAA *CEAIEA/SQRT(BKT) * EXP(-Y)*(1. + CEAB*F1Y)
        ELSE IF(ARELEM .EQ. 'FE' .AND.  (ICHRGE .EQ. 3)) THEN 
          CEAA = 1.8E-17
          CEAIEA= 60.
          CEAB=1.0         
          Y=CEAIEA/BKT
          F1Y=FONE(Y)
           Defined=.True.
          CUP=6.69E+7*CEAA *CEAIEA/SQRT(BKT) * EXP(-Y)*(1. + CEAB*F1Y)
        ELSE IF(ARELEM .EQ. 'FE' .AND.  (ICHRGE .EQ. 4)) THEN 
         CEAA = 5.0E-17
         CEAIEA= 73.
         CEAB=1.0         
         Y=CEAIEA/BKT
         F1Y=FONE(Y)
         Defined=.True.
         CUP=6.69E+7*CEAA *CEAIEA/SQRT(BKT) * EXP(-Y)*(1. + CEAB*F1Y)
        ENDIF
        If (.not. Defined) then
           Print *,'There has been some error computing collisional rates in routine AR85CEA'
           Print *,'ARELEM=',ARELEM
           Stop
        End if
        CEA(K)=CUP
    1 CONTINUE

      RETURN
      End Subroutine AR85CEA
                                      
      FUNCTION SUMMERS(ILO,IHI,EDENS,ATOM)
!
!  NEW ROUTINE FOR COMPUTING COLLISIONAL REDUCTION OF DIELECTRONIC 
!  RECOMBINATION RATE BY ELECTRON COLLISIONS FOLLOWING SUMMERS 1974
!  Inputs:
!    ILO  - index of lower level of transition
!    IHI  - index of upper level of transition
!    EDENS - electron density in /cm3
! SUMMERS 94-11-09  NEW ROUTINE: (PHILIP JUDGE)
!        
      Type (NLTE_Atom) :: Atom
      CHARACTER (len=2) ::   CSEQ, ARELEM
      LOGICAL DEBUG
      DATA DEBUG /.FALSE./
      SUMMERS=1.0
!
!
!  FIND ELEMENT
!
      ARELEM = ATOM%element(1:2)
      IZ=IATOMN(ARELEM)   
      ZZ=0.0 + IZ
      IF (IZ .LT. 1 .OR. IZ .GT. 92) THEN 
        Print *,'SUMMERS:  ATOMIC NUMBER = ',IZ, ' FOR ',ARELEM
        Print *,'SUMMERS:  NO DENSITY EFFECTS IN DIEL RATE'
        RETURN
      ELSE IF(DEBUG) THEN
!        WRITE(LJOBLO,*)'SUMMERS: ',ARELEM, ' CHARGE ',ION(ILO)
      ENDIF
! 
!  FIND ISO-ELECTRONIC SEQUENCE OF RECOMBINING ION
!
      ICHRGE=Atom%ION(ILO)
      ZZ=FLOAT(ICHRGE)
      ISOSEQ=IZ-ICHRGE
      CSEQ=ATOMNM(ISOSEQ)  
!
!*********************************************************************
! LITHIUM, SODIUM, POTASSIUM SEQUENCES
!*********************************************************************
!
      RHO0=2000.0
      IF (CSEQ .EQ. 'LI' .OR. CSEQ .EQ. 'NA' .OR.  &
         CSEQ .EQ. 'K') RHO0=30.0
      RHOQ=EDENS/ZZ/ZZ/ZZ/ZZ/ZZ/ZZ/ZZ
      SUMMERS=1. / (1.+ RHOQ/RHO0)**0.14
      IF(DEBUG) THEN
!         WRITE(LJOBLO,*)'SUMMERS: ',LABEL(IHI),' -> ',  &
!          LABEL(ILO), ' SEQUENCE =', CSEQ
!         WRITE(LJOBLO,*)'SUMMERS: DENSITY FACTOR = ',SUMMERS,  &
!            ' NE = ', EDENS, 'RHO0=',RHO0
      ENDIF
      RETURN
    END FUNCTION SUMMERS

!                                                                       
!********************************************************************** 
!                                                                       
      FUNCTION FTWO (X)                                                  
!                                                                       
! NEW ROUTINE 24-JAN-1994  P. JUDGE                                     
! CALCULATES F2(X) NEEDED FOR COLLISIONAL RATES OF ARNAUD AND ROTHENFLUG
! IF ARGUMENT X IS < BREAK, THEN USE THE EQNS 5 AND 6 OF HUMMER (1983, J
! BREAK GIVEN BY HUMMER AS 4.0                                          
!                                                                       
! NOTE THAT THE SUGGESTED POLYNOMIAL EXPANSION OF ARNAUD AND ROTHENFLUG 
! FAILS FOR ARGUMENTS < BREAK- SEE HUMMER'S ORIGINAL PAPER.             
! IT IS UNCLEAR IF THE TABULATED IONIZATION FRACTIONS OF ARNAUD AND ROTH
! CONTAIN THIS ERROR.                                                   
!                                                                       
      DIMENSION P (15) , Q (15)                                          
      DATA BREAK / 4.0 /                                                 
      DATA P / 1.0, 2.1658E+02, 2.0336E+04, 1.0911E+06, 3.7114E+07,     &
      8.3963E+08, 1.2889E+10, 1.3449E+11, 9.4002E+11, 4.2571E+12,       &
      1.1743E+13, 1.7549E+13, 1.0806E+13, 4.9776E+11, 0.0000 / , Q /    &
      1.0, 2.1958E+02, 2.0984E+04, 1.1517E+06, 4.0349E+07, 9.4900E+08,  &
      1.5345E+10, 1.7182E+11, 1.3249E+12, 6.9071E+12, 2.3531E+13,       &
      4.9432E+13, 5.7760E+13, 3.0225E+13, 3.3641E+12 /                  
      IF (X.GE.BREAK) THEN                                               
         PX = P (1)                                                      
         XFACT = 1.0                                                     
         DO 1 I = 2, 15                                                  
            XFACT = XFACT / X                                            
            PX = PX + P (I) * XFACT                                      
    1    END DO                                                          
!        PX = P(1) + XI*(P(2) + (XI*P(3) +(XI*P(4) + (XI*P(5) + (XI*P(6)
!     *   + (XI*P(7) + (XI*P(8) +(XI*P(9) + (XI*P(10) + (XI*P(11) +     
!     *   (XI*P(12) + (XI*P(13) +(XI*P(14)+XI*P(15))))))))))))))        
!        QX = Q(1) + XI*(Q(2) + (XI*Q(3) +(XI*Q(4) + (XI*Q(5) + (XI*Q(6)
!     *   + (XI*Q(7) + (XI*Q(8) +(XI*Q(9) + (XI*Q(10) + (XI*Q(11) +     
!     *   (XI*Q(12) + (XI*Q(13) +(XI*Q(14)+XI*Q(15))))))))))))))        
         QX = Q (1)                                                      
         XFACT = 1.0                                                     
         DO 2 I = 2, 15                                                  
            XFACT = XFACT / X                                            
            QX = QX + Q (I) * XFACT                                      
    2    END DO                                                          
         FTWO = PX / QX / X / X                                          
      ELSE                                                               
!                                                                       
!  HUMMER'S EQUNS 5 AND 6                                               
!  GAMMA IS EULER'S CONSTANT (ABRAMOVICZ+STEGUN)                        
!                                                                       
         GAMMA = 0.5772156649                                            
         F0X = PI * PI / 12.                                             
         TERM = 1.0                                                      
         COUNT = 0.0                                                     
         FACT = 1.0                                                      
         XFACT = 1.0                                                     
         DO WHILE (ABS (TERM / F0X) .GT.1.E-8)                           
         COUNT = COUNT + 1.0                                             
         FACT = FACT * COUNT                                             
         XFACT = XFACT * ( - X)                                          
         TERM = XFACT / COUNT / COUNT / FACT                             
         F0X = F0X + TERM                                                
         IF (COUNT.GT.100.) then
            print *,'FTWO: TOO MANY ITERATIONS'
            Stop
         End if
         ENDDO                                                           
         FTWO = EXP (X) * ( (LOG (X) + GAMMA) * (LOG (X) + GAMMA)       &
         * 0.5 + F0X)                                                   
      ENDIF                                                              
      RETURN                                                             
      END FUNCTION FTWO                             
!                                                                       
! ****************************************************************      
!                                                                       
      FUNCTION SPLIN (T, X, Y, N, NNN)                                   
!                                                                       
!  SPLINE INTERPOLATION OF T IN A GIVEN TABLE OF POINTS(X(I),Y(I)).     
!                                                                       
! Y(I) IS A FUNCTION OF X(I) - X(N)>X(N-1)>....X(1) -OR- X(1)> X(2)     
! > ....> X(N) - WHEN X IS NORMALISED X(MIN) = 0, XMAX=1,               
! NORM(ABS(X(I+1)-X(I))) < 1.0E-8,(I=1,N-1). THE ARRAY X IS DIVIDED     
! INTO OVERLAPPING INTERVALS. EACH INTERVAL CONSISTS OF NNN POINTS,     
! THE OVERLAP IS INT(LOG(NNN)+3) POINTS. ( 3<=NNN<=100,N).              
! EACH TIME AFTER A CALL WHEN T CHANGES OF INTERVAL OR AFTER A CALL     
! WITH ANOTHER SET OF POINTS, A NEW CUBIC SPLINE OF NNN POINTS HAS      
! TO BE COMPUTED, THEREFORE THE COMPUTING TIME DEPENDS STRONGLY OF      
! THE SUCCESSION OF CALLS.                                              
! FOR EXAMPLE: SUCCESSIVE CALLS, WHERE T IS RANDOMLY CHOSEN WILL COST   
! A LOT OF COMPUTING TIME - (T OFTEN CHANGES OF INTERVAL) --->          
! TAKE NNN THEN AS LARGE AS POSSIBLE. SUCCESSIVE CALLS WHERE T CHANGES  
! SLIGHTLY COST LESS COMPUTING TIME, SO IT WILL BE CLEAR THAT YOUR      
! CHOICE OF NNN DEPENDS HEAVILY UPON ITS WAY OF USEAGE.                 
!                                                                       
! COMPUTING TIME:   FOR THE COMPUTATION OF SPLINE COEFFICIENTS          
!        THIS ROUTINE NEEDS ABOUT 21*N MULTIPLICATIONS OR DIVISIONS     
!        FOR THE INTERPOLATION BY MEANS OF THE SPLINE COEFFICIENTS      
!        THE ROUTINE NEEDS 9 MULTIPLICATIONS OR DIVISIONS.              
!                                                                       
! EXTRAPOLATION: IF T IS NOT IN [X(1),X(N)] T IS LINEAR EXTRAPOLATED    
!                                                                       
! A PART OF THE ALGORITM FOR THE COMPUTATION OF A NATURAL SPLINE        
! IS TAKEN FROM T.N.E. GREVILLE, THE LINEAR SYSTEM IS SOLVED BY         
! THE METHOD OF SUCCESSIVE OVERRELAXATION. (ERROR = 1.0E-6)             
!                                                                       
! REF: T.N.E. GREVILLE, MATHEMATICS RESEARCH CENTER, U.S. ARMY,         
!        UNIVERSITY OF WISCONSIN. MATHEMATICAL METHODS.                 
!                                                                       
! THIS ROUTINE WAS IMPLEMENTED BY E.B.J. VAN DER ZALM, STERRENWACHT     
! UTRECHT (MAY 7TH 1981). ADAPTED FOR VAX\11 BY PAUL KUIN AT OXFORD.    
! WRITTEN IN REAL*4 BY P.G. JUDGE, OXFORD                               
! MODIFIED SO THAT IF THE ARGUMENT IS OUTSIDE THE RANGE OF X-VALUES     
! THEN THE FIRST (OR LAST) Y-VALUE IS TAKEN                             
!                                                                       
! INPUT:         T        ARGUMENT   (REAL)                             
!                X        ARRAY OF ARGUMENTS        (REAL)              
!                Y        ARRAY OF FUNCTION VALUES (REAL)               
!                N        LENGTH OF THE ARRAYS X,Y                      
!                        IF N = 2 ---> LINEAR INTERPOLATION             
!                NNN        NUMBER OF POINTS OF THE CUBIC SPLINE        
!                        IF NNN > N   ---> NNN = N                      
!                        IF NNN > 100 ---> NNN = 100                    
!                        IF NNN = 2 LINEAR INTERPOLATION ADOPTED.       
!                                                                       
!                        IF NNN IS ZERO NNN=7  (ALTERED BY P JUDGE)     
!                       (USED TO BE IF NNN OMITTED )                    
! OUTPUT:        SPLIN INTERPOLATED VALUE AT T                          
!                                                                       
      DIMENSION X (N) , Y (N) , S2 (100) , S3 (100) , DELY (100) , H (  &
      100) , H2 (100) , DELSQY (100) , C (100) , B (100) , XX (100) ,   &
      YY (100)                                                          
      LOGICAL FOUND                                              
      CHARACTER(4) INIT                                                  
      DATA INIT / 'NOTO' /                                               
      SAVE                                                               
!                                                                       
!  NN = NUMBER OF POINTS OF THE SPLINE INTERVAL                         
!                                                                       
      NN = MIN (7, N)                                                    
      IF (NNN.EQ.0) GOTO 10                                              
! PGJ        IF (%LOC(NNN).EQ.0) GO TO 10                               
      NN = NNN                                                           
      NN = MIN (N, NNN, 100)                                             
!                                                                       
!           IF T IN [X(K)-ERROR,X(K)+ERROR] T = X(K)                    
!                                                                       
   10 ERROR = 1.0E-12                                                    
      IF (T.NE.0) ERROR = ABS (1.0E-12 * T)                              
      FOUND = SEARCH (X, T, 1, N, K, ERROR)                              
      IF (.NOT.FOUND) GOTO 20                                            
      SPLIN = Y (K)                                                      
      RETURN                                                             
   20 IF (N.NE.2.AND.NN.NE.2) GOTO 30                                    
!                                                                       
!           IF N=2 OR NNN=2 LINEAR INTERPOLATION IS ADOPTED             
!                                                                       
      SPLIN = (Y (K + 1) - Y (K) ) * (T - X (K) ) / (X (K + 1) - X (K) )&
      + Y (K)                                                           
      RETURN                                                             
   30 IF (N.NE.1) GOTO 40                                                
      SPLIN = Y (1)                                                      
      RETURN                                                             
!                                                                       
!           INTERVAL OF THE SPLINE IS COMPUTED                          
!                                                                       
   40 IOVER = LOG10 (FLOAT (NN) ) * 3                                    
      INT = (K - IOVER) / (NN - 2 * IOVER) + 1                           
      I1 = (INT - 1) * (NN - 2 * IOVER)                                  
      I2 = I1 + NN - 1                                                   
      IF (I1.GE.1) GOTO 50                                               
      I2 = NN                                                            
      I1 = 1                                                             
   50 IF (I2.LE.N) GOTO 60                                               
      I2 = N                                                             
      I1 = N - NN + 1                                                    
   60 ISHIFT = I1 - 1                                                    
      N1 = NN - 1                                                        
!                                                                       
!           INITIALIZATION CHECK                                        
!                                                                       
      IF (INIT.NE.'OKAY') GOTO 80                                        
!                                                                       
!           INTERVAL CHECK                                              
!                                                                       
      DO 70 I = 1, NN                                                    
         IF (XX (I) .NE.X (I + ISHIFT) ) GOTO 90                         
         IF (YY (I) .NE.Y (I + ISHIFT) ) GOTO 90                         
   70 END DO                                                             
      GOTO 200                                                           
   80 INIT = 'OKAY'                                                      
!                                                                       
!           START OF THE COMPUTATION OF SPLINE COEFFICIENTS             
!                                                                       
   90 RMIN = 1.0E37                                                      
      RMAX = - 1.0E37                                                    
      DO 100 I = 1, NN                                                   
         IF (Y (I + ISHIFT) .LT.RMIN) RMIN = Y (I + ISHIFT)              
         IF (Y (I + ISHIFT) .GT.RMAX) RMAX = Y (I + ISHIFT)              
  100 END DO                                                             
!                                                                       
!           COMPUTATION OF THE NORM FACTORS XNORM AND YNORM             
!                                                                       
      XNORM = 1 / (X (NN + ISHIFT) - X (1 + ISHIFT) )                    
      YNORM = RMAX - RMIN                                                
      IF (YNORM.EQ.0) YNORM = RMAX                                       
      IF (YNORM.EQ.0) YNORM = 1                                          
      YNORM = 1 / YNORM                                                  
      DO 110 I = 1, N1                                                   
         H (I) = (X (I + 1 + ISHIFT) - X (I + ISHIFT) ) * XNORM          
         DELY (I) = ( (Y (I + 1 + ISHIFT) - Y (I + ISHIFT) ) / H (I) )  &
         * YNORM                                                        
  110 END DO                                                             
      DO 120 I = 2, N1                                                   
         H2 (I) = H (I - 1) + H (I)                                      
         B (I) = 0.5 * H (I - 1) / H2 (I)                                
         DELSQY (I) = (DELY (I) - DELY (I - 1) ) / H2 (I)                
         S2 (I) = 2. * DELSQY (I)                                        
         C (I) = 3. * DELSQY (I)                                         
  120 END DO                                                             
      S2 (1) = 0                                                         
      S2 (N1 + 1) = 0                                                    
!                                                                       
!        SOLUTION OF THE LINEAR SYSTEM OF THE SPLINE COEFFICIENTS       
!        BY SUCCESSIVE OVERRELAXATION.                                  
!        CONSTANTS: EPS = ERROR CRITERION IN THE ITERATIVE SOLUTION.    
!                   OMEGA = RELAXATION COEFFICIENT.                     
!                                                                       
      EPSLN = 1.0E-8                                                     
      OMEGA = 1.0717968                                                  
  130 ETA = 0.                                                           
      DO 160 I = 2, N1                                                   
         W = (C (I) - B (I) * S2 (I - 1) - (0.5 - B (I) ) * S2 (I + 1)  &
         - S2 (I) ) * OMEGA                                             
         IF ((ABS (W) - ETA) .LT. 0) GOTO 150
         IF ((ABS (W) - ETA) .EQ. 0) GOTO 150
         IF ((ABS (W) - ETA) .GT. 0) GOTO 140
  140    ETA = ABS (W)                                                   
  150    S2 (I) = S2 (I) + W                                             
  160 END DO                                                             
      IF ((ETA - EPSLN) .LT. 0) GOTO 170
      IF ((ETA - EPSLN) .EQ. 0) GOTO 130
      IF ((ETA - EPSLN) .GT. 0) GOTO 130
  170 DO 180 I = 1, N1                                                   
         S3 (I) = (S2 (I + 1) - S2 (I) ) / H (I)                         
  180 END DO                                                             
!                                                                       
!        X AND Y STORED IN XX AND YY FOR INTERVAL CHECK                 
!                                                                       
      DO 190 I = 1, NN                                                   
         XX (I) = X (I + ISHIFT)                                         
         YY (I) = Y (I + ISHIFT)                                         
  190 END DO                                                             
  200 I = K - I1 + 1                                                     
      HT1 = (T - X (I + ISHIFT) ) * XNORM                                
      HT2 = (T - X (I + ISHIFT + 1) ) * XNORM                            
      IF ( (T - X (1) ) * XNORM.GT.0) GOTO 210                           
!                                                                       
!        EXTRAPOLATION   T < X(1)                                       
!                                                                       
      SPLIN = Y (1)                                                      
      RETURN                                                             
  210 IF ( (T - X (N) ) * XNORM.GT.0) GOTO 220                           
!                                                                       
!        INTERPOLATION BY MEANS OF THE SPLINE COEFFICIENTS              
!                                                                       
      SS2 = S2 (I) + HT1 * S3 (I)                                        
      PROD = HT1 * HT2                                                   
      DELSQS = (S2 (I) + S2 (I + 1) + SS2) / 6.                          
      SPLIN = Y (I + ISHIFT) + (HT1 * DELY (I) + PROD * DELSQS) / YNORM  
      RETURN                                                             
!                                                                       
!        EXTRAPOLATION   T > X(N)                                       
!                                                                       
  220 SPLIN = Y (N)                                                      
      RETURN                                                             
      END FUNCTION SPLIN            
!                                                                       
!***********************************************************************
!                                                                       
      FUNCTION SEMIC (Z, EUP, ELO, FEM, TEMP, IFLAG)                     
!                                                                       
!  RETURNS COLLISION RATE DOWNWARDS USING SEMI-EMPIRICAL GF'S           
!                                                                       
!  P.G. JUDGE, NOVEMBER 1987.                                           
!                                                                       
!  INPUT:                                                               
!         Z      - CHARGE OF ION+1 (E.G. FOR C II Z=2)                  
!         EUP    - ENERGY OF UPPER LEVEL IN CM-1                        
!         ELO    - ENERGY OF LOWER LEVEL IN CM-1                        
!         FEM    - EMISSION OSCILLATOR STRENGTH ( = GF/G(UPPER))        
!         TEMP   - TEMPERATURE IN DEGREES KELVIN                        
!                                                                       
!  OUTPUT:                                                              
!         IFLAG  -  = 1 FOR VAN REGEMORTER APPROXIMATION                
!         IFLAG  -  = 2 FOR SHEVELKO       APPROXIMATION                
!         SEMIC-  DOWNWARD COLLISION RATE IN CM3 /S                     
!                                                                       
      INTEGER Z                                                          
      DATA C1, RY / 1.43882, 109737.312 /                                
!                                                                       
! INITIALIZE                                                            
!                                                                       
      DELTE = EUP - ELO                                                  
      BETA = DELTE * C1 / TEMP                                           
      SEMIC = 0.0                                                        
!                                                                       
!  IFLAG=1 => VAN REGEMORTER APPROXN.                                   
!                                                                       
      IF (BETA.GE.0.01) THEN                                             
         IFLAG = 1                                                       
         SEMIC = 3.2E-7 * FEM * (RY / DELTE) **1.5 * SQRT (BETA)        &
         * PSEMI (Z, BETA)                                              
!                                                                       
!  IFLAG=2 => SHEVELKO APPROX                                           
!                                                                       
      ELSEIF (BETA.LT.0.01) THEN                                         
         IFLAG = 2                                                       
         SEMIC = 1.74E-07 * FEM * EXP (BETA) / SQRT (BETA * DELTE / RY) &
         * LOG (ELO / DELTE * SQRT (RY / C1 / TEMP) )                   
      ENDIF                                                              
!                                                                       
      RETURN                                                             
      END FUNCTION SEMIC   
!                                                                       
! **********************************************************************
!                                                                       
      FUNCTION PSEMI (Z, B)                                              
!                                                                       
!  THERMAL P-FUNCTION FOR SEMIC EMPIRICAL COLLISION RATES               
!  REFERENCE: SOBEL'MAN - 'ATOMIC PHYSICS II.  '                        
!                                                                       
!:                                                                      
!: PSEMI  94-08-29  MODIFICATIONS: (MATS CARLSSON)                      
!:        MULTIPLE RETURN CHANGED TO SINGLE RETURN TO AVOID             
!:        COMPILER BUG IN DIGITAL FORTRAN                               
!:                                                                      
      PARAMETER (N1 = 10)                                                
      DIMENSION BETREF (N1) , PNREF (N1) , PCREF (N1)                    
      INTEGER Z                                                          
      DATA BETREF / - 2.0, - 1.699, - 1.398, - 1.0, - 0.699, - 0.398,   &
      0.0, 0.301, 0.602, 1.0 /                                          
      DATA PNREF / 1.160, 0.956, 0.758, 0.493, 0.331, 0.209, 0.100,     &
      0.063, 0.040, 0.023 /                                             
      DATA PCREF / 1.160, 0.977, 0.788, 0.554, 0.403, 0.290, 0.214,     &
      0.201, 0.200, 0.200 /                                             
!
      PSEMI=0.
!                                                                       
!  LIMIT OF HIGH TEMPERATURE: KT >> E                                   
!                                                                       
      IF (B.LT.0.01) THEN                                                
         PSEMI = 0.27566 * (0.577 + LOG (B) )                            
!                                                                       
!  INTERMEDIATE TEMPS (MOST IMPORTANT FOR EQM PLASMAS)                  
!  SPLINE INTERPOLATION ONTO LOGB GRID:                                 
!                                                                       
      ELSEIF (B.GT.0.01.AND.B.LT.10.0) THEN                              
         BB = LOG10 (B)                                                  
         IF (Z.EQ.1) PSEMI = SPLIN (BB, BETREF, PNREF, N1, N1)           
         IF (Z.GT.1) PSEMI = SPLIN (BB, BETREF, PCREF, N1, N1)           
!                                                                       
!  LIMIT OF LOW TEMPERATURE: KT << E                                    
!                                                                       
      ELSEIF (B.GT.10.0) THEN                                            
         IF (Z.EQ.1) PSEMI = 0.066 / SQRT (B)                            
         IF (Z.GT.1) PSEMI = 0.200                                       
      ENDIF                                                              
      RETURN                                                             
      END FUNCTION PSEMI        
!
      !                                                                       
! **********************************************************************
!                                                                       
      LOGICAL FUNCTION SEARCH (ARR, X, IA, IB, K, ERR)                   
!                                                                       
! * SEARCHES THE POINT X WITHIN AN ERROR BOUND -ERR- IN ARRAY ARR      *
! * IF THAT POINT IS NOT FOUND, IT GIVES THE                           *
! * INTERVAL IN THE ARRAY, WHERE THAT POINT FITS.                      *
! * THE INTERVAL K IF X IN [ARR(K),ARR(K+1)), INTERVAL 1 IF X IN       *
! * (-INF,ARR(2)), INTERVAL N-1 IF X IN [ARR(N-1),INF)-(ARR INCREASING)*
! * INTERVAL 1 IF X IN (INF,ARR(2)),INTERVAL N-1 IF X IN [ARR(N-1),    *
! * -INF)-(ARR DECREASING)                                             *
! * IN THE WORST CASE A SEARCH COSTS 2*LOG2(N) CYCLES.                 *
! * THIS ROUTINE IS VERY ECONOMIC, WHEN IN A SEQUENCE OF SEARCHES,     *
! * THE DIFFERENCES BETWEEN THE POINTS WHICH ARE TO BE SEARCHED ARE    *
! * SMALL.                                                             *
! * FOR EXAMPLE: A SEQUENCE OF SEARCHES WHICH GOES ALONG THE ARRAY,    *
! * OR A SEQUENCE OF SEARCHES WHICH REMAINS IN A SMALL AREA OF         *
! * THE ARRAY.                                                         *
! *                                                                    *
! * INPUT:     - ARR         ARRAY WITH A NON DECREASING OR A NON      *
! *                          INCREASING FUCNTION.                      *
! *            - IA,IB       LOWER AND UPPER LIMIT OF THE ARRAY WITHIN *
! *                          IS SEARCHED.                              *
! *            - X           VALUE WHICH YOU WANT TO SEARCH            *
! *            - K           LAST INDEX WHICH WAS FOUND                *
! *                          SEARCH                                    *
! *            - ERR         ERROR BOUND                               *
! *                                                                    *
! * OUTPUT:    - K           INDEX OF THE POINT OR INTERVAL OF THE     *
! *                          ARRAY WHICH IS FOUND.                     *
! *            - SEARCH      TRUE - X IS FOUND WITHIN THE ERRORBOUND   *
! *                                  -ERR-                             *
! *                          FALSE - X IS NOT FOUND.                   *
! *                                                                    *
! * SEND IN BY E.V.D.ZALM, UTRECHT, STERREWACHT, MARCH 24TH 1981       *
! **********************************************************************
      DIMENSION ARR (IB)                                                 
!                                                                       
      ONE = 1.0                                                          
      SEARCH = .TRUE.                                                    
      IF (ARR (IB) - ARR (IA) .LT.0.0) GOTO 100                          
    5 IF (X.LT.ARR (IA) .OR.X.GT.ARR (IB) ) GOTO 25                      
      IF (K.LT.IA.OR.K.GT.IB - 1) GOTO 25                                
      L = K                                                              
      IP = SIGN (ONE, X - ARR (K) )                                      
      IF (IP.LT.0.) GOTO 20                                              
   10 L = MIN (K + IP, IB)                                               
      IF (ARR (L) .GE.X) GOTO 30                                         
      IP = IP * 2                                                        
      K = L                                                              
      GOTO 10                                                            
   20 K = MAX (L + IP, IA)                                               
      IF (ARR (K) .LE.X) GOTO 30                                         
      IP = IP * 2                                                        
      L = K                                                              
      GOTO 20                                                            
   25 K = IA                                                             
      L = IB                                                             
   30 IF ( (L - K) .LE.1) GOTO 50                                        
      I = INT ( (K + L) / 2.)                                            
      IF (ARR (I) .GE.X) GOTO 40                                         
      K = I                                                              
      GOTO 30                                                            
   40 L = I                                                              
      GOTO 30                                                            
  100 IF (X.LT.ARR (IB) .OR.X.GT.ARR (IA) ) GOTO 250                     
      IF (K.LT.IA.OR.K.GT.IB - 1) GOTO 250                               
      L = K                                                              
      IP = SIGN (ONE, ARR (K) - X)                                       
      IF (IP.LT.0.) GOTO 200                                             
  110 L = MIN (K + IP, IB)                                               
      IF (ARR (L) .LE.X) GOTO 300                                        
      IP = IP * 2                                                        
      K = L                                                              
      GOTO 110                                                           
  200 K = MAX (L + IP, IA)                                               
      IF (ARR (K) .GE.X) GOTO 300                                        
      IP = IP * 2                                                        
      L = K                                                              
      GOTO 200                                                           
  250 K = IA                                                             
      L = IB                                                             
  300 IF ( (L - K) .LE.1) GOTO 50                                        
      I = INT ( (K + L) / 2.)                                            
      IF (ARR (I) .LE.X) GOTO 400                                        
      K = I                                                              
      GOTO 300                                                           
  400 L = I                                                              
      GOTO 300                                                           
   50 IF (ABS (X - ARR (K) ) .LE.ERR) RETURN                             
      IF (ABS (X - ARR (L) ) .GT.ERR) GOTO 60                            
      K = L                                                              
      RETURN                                                             
   60 SEARCH = .FALSE.                                                   
      RETURN                                                             
      END FUNCTION SEARCH            
!
!                                                                       
!*********************************************************************  
!                                                                       
      INTEGER FUNCTION IATOMN (STRING)                                   
!                                                                       
!  ATOMN 94-02-22  NEW ROUTINE: (PHILIP JUDGE)                          
!  GIVES ATOMIC NUMBER OF ARELEMENT IF STRING IS A STRING CONTAINING    
!  THE NAME OF THE ARELEMENT, E.G. IF INPUT  IS 'H ' IT WILL RETURN 1   
!                                                                       
      PARAMETER (NDATA = 28)                                             
      CHARACTER(2) ARELEM (NDATA), ARELEM2(NDATA)
      CHARACTER ( * ) STRING                                             
      DATA ARELEM / 'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F '&
     &, 'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA'&
     &, 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI' /                
      DATA ARELEM2 / 'h ', 'he', 'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f '&
     &, 'ne', 'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar', 'k ', 'ca'&
     &, 'sc', 'ti', 'v ', 'cr', 'mn', 'fe', 'co', 'ni' /                
      IATOMN = - 1                                                       
      DO 1 I = 1, NDATA                                                  
         IF (STRING.EQ.ARELEM (I) .or. String.EQ. ARELEM2(I)) THEN                                 
            IATOMN = I                                                   
            RETURN                                                       
         ENDIF                                                           
    1 END DO                                                             
      END FUNCTION IATOMN                           
!                                                                       
!*********************************************************************  
!                                                                       
      CHARACTER(2) FUNCTION ATOMNM (I)                                   
!                                                                       
!  ATOMNM 94-02-22  NEW ROUTINE: (PHILIP JUDGE)                         
!  GIVES ATOMIC NAME OF ARELEMENT IF I IS AN INTEGER CONTAINING         
!  E.G. IF INPUT  IS 'H ' IT WILL RETURN 1, 'HE', IT WILL RETURN 2, ETC.
!                                                                       
      PARAMETER (NDATA = 28)                                             
      CHARACTER(2) ARELEM (NDATA)                                        
      DATA ARELEM / 'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F '&
     &, 'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA'&
     &, 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI' /                
      IF (I.LE.NDATA.AND.I.GE.1) THEN                                    
         ATOMNM = ARELEM (I)                                             
      ELSE                                                               
      ATOMNM = '  '                                                      
      ENDIF                                                              
      RETURN                                                             
      END FUNCTION ATOMNM              
!                                                                       
!*********************************************************************  
!                                                                       
      FUNCTION EXPINT (N, X, EX)                                         
!                                                                       
!  COMPUTES THE N'TH EXPONENTIAL INTEGRAL OF X                          
!  INPUT - X,  INDEPENDENT VARIABLE (-100. .LE. X .LE. +100.)           
!         N,  ORDER OF DESIRED EXPONENTIAL INTEGRAL (1 .LE. N .LE. 8)   
!  OUTPUT - EXPINT,  THE DESIRED RESULT                                 
!           EX,  EXPF(-X)                                               
!  NOTE   RETURNS WITH E1(0)=0, (NOT INFINITY).                         
!  BASED ON THE SHARE ROUTINE NUEXPI, WRITTEN BY J. W. COOLEY,          
!  COURANT INSTITUTE OF MATHEMATICAL SCIENCES, NEW YORK UNIVERSITY      
!  OBTAINED FROM RUDOLF LOESER                                          
!-----GENERAL COMPILATION OF 1 AUGUST 1967.                             
!                                                                       
      DIMENSION TAB (20) , XINT (7)                                      
      DATA XINT / 1., 2., 3., 4., 5., 6., 7. /                           
      DATA TAB / .2707662555, .2131473101, .1746297218, .1477309984,    &
      .1280843565, .1131470205, .1014028126, .0919145454, .0840790292,  &
      .0774922515, .0718735405, .0670215610, .0627878642, .0590604044,  &
      .0557529077, .0527977953, .0501413386, .0477402600, .0455592945,  &
      .0435694088 /                                                     
      DATA XSAVE / 0. /                                                  
!                                                                       
      U = X                                                              
      IF (U .LT. 0) GOTO 603
      IF (U .EQ. 0) GOTO 602
      IF (U .GT. 0) GOTO 603
  602 EX = 1.                                                            
      IF ((N - 1) .LT. 0) GOTO 800
      IF ((N - 1) .EQ. 0) GOTO 800
      IF ((N - 1) .GT. 0) GOTO 801
  800 EXPINT = 0.                                                        
      GOTO 777                                                           
  801 EXPINT = 1. / XINT (N - 1)                                         
      GOTO 777                                                           
  603 IF ((U - XSAVE) .LT. 0) GOTO 604
      IF ((U - XSAVE) .EQ. 0) GOTO 503
      IF ((U - XSAVE) .GT. 0) GOTO 604
  604 XSAVE = U                                                          
      XM = - U                                                           
      EMX = EXP (XM)                                                     
!                                                                       
!  SELECT METHOD FOR COMPUTING EI(XM)                                   
!                                                                       
      IF ((XM - 24.5) .LT. 0) GOTO 501
      IF ((XM - 24.5) .EQ. 0) GOTO 400
      IF ((XM - 24.5) .GT. 0) GOTO 400
  501 IF ((XM - 5.) .LT. 0) GOTO 502
      IF ((XM - 5.) .EQ. 0) GOTO 300
      IF ((XM - 5.) .GT. 0) GOTO 300
  502 IF ((XM + 1.) .LT. 0) GOTO 100
      IF ((XM + 1.) .EQ. 0) GOTO 200
      IF ((XM + 1.) .GT. 0) GOTO 200
  503 EISAVE = - ARG                                                     
      EXSAVE = EMX                                                       
!                                                                       
!  NOW RECURSE TO HIGHER ORDERS                                         
!                                                                       
      IF ((N - 1) .LT. 0) GOTO 507
      IF ((N - 1) .EQ. 0) GOTO 507
      IF ((N - 1) .GT. 0) GOTO 505
  505 DO 506 I = 2, N                                                    
         EISAVE = (U * EISAVE-EXSAVE) / ( - XINT (I - 1) )               
  506 END DO                                                             
  507 EXPINT = EISAVE                                                    
      EX = EXSAVE                                                        
  777 RETURN                                                             
!                                                                       
!  EI(XM) FOR XM .LT. -1.0                                              
!  HASTINGS POLYNOMIAL APPROXIMATION                                    
!                                                                       
  100 ARG = ( ( ( ( ( (U + 8.573328740) * U + 18.05901697) * U +        &
      8.634760893) * U + .2677737343) / XM) * EMX) / ( ( ( (U +         &
      9.573322345) * U + 25.63295615) * U + 21.09965308) * U +          &
      3.958496923)                                                      
      GOTO 503                                                           
!     EI(XM) FOR -1. .LE. XM .LT. 5.0                                   
!     POWER SERIES EXPANSION ABOUT ZERO                                 
  200 ARG = LOG (ABS (XM) )                                              
      ARG = ( ( ( ( ( ( ( ( ( ( ( ( ( ( ( (.41159050E-14 * XM +         &
      .71745406E-13) * XM + .76404637E-12) * XM + .11395905E-10)        &
      * XM + .17540077E-9) * XM + .23002666E-8) * XM + .27536018E-7)    &
      * XM + .30588626E-6) * XM + .31003842E-5) * XM + .28346991E-4)    &
      * XM + .23148057E-3) * XM + .0016666574) * XM + .010416668)       &
      * XM + .055555572) * XM + .25) * XM + .99999999) * XM + .57721566)&
      + ARG                                                             
      GOTO 503                                                           
!                                                                       
!  EI(XM) FOR 5.0 .LE. XM .LT. 24.5                                     
!  TABLE LOOK-UP AND INTERPOLATION                                      
!                                                                       
  300 I = XM + .5                                                        
      XZERO = I                                                          
      DELTA = XZERO - XM                                                 
      ARG = TAB (I - 4)                                                  
      IF (DELTA .LT. 0) GOTO 303
      IF (DELTA .EQ. 0) GOTO 305
      IF (DELTA .GT. 0) GOTO 303
  303 Y = ARG                                                            
      DELTAX = DELTA / XZERO                                             
      POWER = 1. / DELTAX                                                
      DO I = 1, 7                                                    
         POWER = POWER * DELTAX                                          
         Y = ( (Y - POWER / XZERO) * DELTA) / XINT (I)                   
         ARG = ARG + Y                                                   
         IF ((ABS (Y / ARG) - 1.E-8) .LT. 0) GOTO 305
         IF ((ABS (Y / ARG) - 1.E-8) .EQ. 0) GOTO 304
         IF ((ABS (Y / ARG) - 1.E-8) .GT. 0) GOTO 304
      END DO
  304 CONTINUE                                                             
  305 ARG = EMX * ARG                                                    
      GOTO 503                                                           
!     EI(XM) FOR 24.5 .LE. XM                                           
!     TRUNCATED CONTINUED FRACTION                                      
  400 ARG = ( ( ( (XM - 15.) * XM + 58.) * XM - 50.) * EMX) / ( ( ( (XM &
      - 16.) * XM + 72.) * XM - 96.) * XM + 24.)                        
      GOTO 503                                                           
      END FUNCTION EXPINT                

Subroutine HCol(Atom, NLTE)
! Adapted from MULTI
 Use Param_structure
 Use NLTE_Atom_structure
 Use NLTE_vars
 Implicit None
 Type (NLTE_Atom) :: Atom
 Type (NLTE_variables) :: NLTE
 Real, Dimension(:,:,:), Allocatable :: C
 Integer :: i, j, idepth, k
 Double Precision T

 Allocate (C(Atom%NK, Atom%NK, NLTE%NDEP))
 C(:,:,:)=0.

 Do idepth=1, NLTE%NDEP
    T=NLTE%Atmo%Temp(idepth)
    Do j=2, Atom%NK-1       
       Do i=1, j-1
          C(i,j,idepth)=NLTE%Atmo%ne(idepth)*H1BB(i,j,T)
          C(j,i,idepth)=C(i,j,idepth)*NLTE%NStar(i,idepth)/NLTE%NStar(j,idepth)
       End do
    End do
    Do i=1, Atom%NK-1
       C(i, Atom%NK, idepth)=NLTE%Atmo%ne(idepth)*H1BF(i,T)
       C(Atom%NK, i, idepth)=C(i, Atom%NK, idepth)*NLTE%NStar(i,idepth)/ &
            NLTE%NStar(Atom%NK,idepth)
    End do
 End do

 NLTE%C=C
 Deallocate (C)
!
End Subroutine HCol

FUNCTION H1BB (I,J,T)
!
!  DETERMINES H COLLISIONAL RATE COEFFICIENTS UP TO N=50
!
! CONVERSION FACTORS:
!
!      CM**-1 TO EV    HH * CC / EE
!      CM**-1 TO K     HH * CC / BK
!
!      K     TO HZ    BK / HH
!      CM**-1 TO HZ    CC
!
! CONSTANTS:
!
!      EE   ELECTRON CHARGE      1.602189E-12         
!      HH   PLANCK CONSTANT      (ERG S)
!      CC   VELOCITY OF LIGHT    (CM S**-1)
!      EM   ELECTRON REST MASS   (G)
!      UU   ATOMIC MASS UNIT     (G)
!      BK   BOLTZMANN CONSTANT   (ERG K**-1)
!      PI   3.14159265359
!
!
! INPUT:
!
!   I        LOWER LEVEL
!   J        UPPER LEVEL
!   T        TEMPERATURE (K)
!
! OUTPUT:
!
!   H1BB    UPWARDS COLLISIONAL RATE                  I TO J
!
!: H1BB   90-02-28  MODIFICATIONS: (MARTIN J. STIFT)
!:        REWRITTEN.  UP TO 50 BOUND LEVELS + ONE CONTINUUM LEVEL
!:
      implicit none
      INTEGER I,J,K
      DOUBLE PRECISION &
              T,S,X, &
              AI,AMN,AJ,BM,BMN,DELMN,DERG,ENG, &
              DION,FIJ,GAMMN,GMX,G0,G1,G2, &
              EION,EV(50),EE,HH,CC,EM,UU,BK,PI,TEV, &
              B1S2P(9),B1S2S(9),RH1,RH2,H1BB
!
! POLYNOMIAL FITS FOR COLLISIONAL RATE COEFFICIENTS (LYMAN ALPHA)
!
! H(1S) TO H (2P)
!
      DATA B1S2P / -2.814949375869D+1, &
                    1.009828023274D+1, &
                   -4.771961915818D+0, &
                    1.467805963618D+0, &
                   -2.979799374553D-1, &
                    3.861631407174D-2, &
                   -3.051685780771D-3, &
                    1.335472720988D-4, &
                   -2.476088392502D-6/
!
! H(1S) TO H(2S)
!
      DATA B1S2S / -2.833259375256D+1, &
                    9.587356325603D+0, &
                   -4.833579851041D+0, &
                    1.415863373520D+0, &
                   -2.537887918825D-1, &
                    2.800713977946D-2, &
                   -1.871408172571D-3, &
                    6.986668318407D-5, &
                   -1.123758504195D-6/
!
! HYDROGEN MEAN ENERGY LEVELS TAKEN FROM BASHKIN & STONER (1975), PP. 2-3
!
      DATA EV /      0.000, &
                 82259.047, &
                 97492.279, &
                102823.882, &
                105291.646, &
                106632.157, &
                107440.441, &
                107965.049, &
                108324.718, &
                41 * 0.000/
!
! HYDROGEN IONIZATION LEVEL TAKEN FROM BASHKIN & STONER (1975), P. 3
!
      DATA EION / 109678.764/
!
! MISCELLANOUS CONSTANTS
!
      DATA EE / 1.602189D-12/,&
           HH / 6.626176D-27/,&
           CC / 2.99792458D10/,&
           EM / 9.109534D-28/,&
           UU / 1.6605655D-24/,&
           BK / 1.380662D-16/,&
           PI / 3.14159265359/
!
! J MUST BE GREATER THAN I
!
      IF (J.LE.I) then
         print *,'H1BB: ERROR IN ORDER OF ENERGY LEVELS'
         stop
      endif
!
      AI = FLOAT(I)
      AJ = FLOAT(J)
!
! DETERMINE ENERGY LEVELS UP TO N=50
!
      DO 1 K = 10,50
         EV(K) = EION * (1. - 1. / FLOAT(K*K))
    1 CONTINUE
!
! EV, EION, ENG  IN   CM**-1
!
! ENG   ENERGY ASSOCIATED TO TEMPERATURE
!
      ENG = T * BK / HH / CC
!
! ENERGY DIFFERENCE AND IONISATION ENERGY FROM LOWER LEVEL
!
      DERG = EV(J) - EV(I)
      DION = EION  - EV(I)
!
! JOHNSON (1972) FORMULA FOR OSCILLATOR STRENGTH
! SEE JANEV ET AL. (1987), PP. 315-316
!
      X = 1. - (AI/AJ) * (AI/AJ)
!
      IF (I.EQ.1) THEN
         G0 =  1.1330
         G1 = -0.4059
         G2 =  0.07014
      ELSE IF (I.EQ.2) THEN
         G0 =  1.0785
         G1 = -0.2319
         G2 =  0.02947
      ELSE IF (I.GE.3) THEN
         G0 =   0.9935 + (0.2328 - 0.1296 / AI) / AI
         G1 = -(0.6282 - (0.5598 - 0.5299 / AI) / AI ) / AI
         G2 =  (0.3887 - (1.1810 - 1.4700 / AI) / AI ) / (AI * AI)
      ENDIF
!
      GMX = G0 + (G1 + G2 / X) / X
!
! OSCILLATOR STRENGTH - AGREES WITH VALUES IN JANEV ET AL. (1987), P. 315
!
      FIJ = 1.96028052 * GMX * (AI/AJ) / (AJ * AJ * X * X * X)
!
! APN, BP AND BPN OF VRIENS AND SMEETS (1980), EQNS. (11), (13) AND (12)
!
      AMN = 2. * EION * FIJ / DERG
!
      BM = ( 1.4 * LOG(AI) - 0.7 - (0.51 - (1.16 - 0.55 / AI) / &
           AI ) / AI ) / AI
!
      BMN = 4. * (EION/DERG) * (EION/DERG) / (AJ * AJ * AJ) * &
            (1. + (4./3.) * (DION/DERG) + &
             BM * (DION/DERG) * (DION/DERG) )
!
! DELTA_PN AND GAMMA_PN OF VRIENS AND SMEETS (1980), EQNS. (18) AND (19)
!
      S = ABS(AI - AJ)
      DELMN = EXP(-BMN/AMN) + 0.06 * (S/AI) * (S/AI) / AJ
      GAMMN = EION * (3. + 11. * (S/AI)*(S/AI)) * &
              LOG(1. + AI * AI * AI * ENG / EION) / &
              (6. + 1.6 * AJ * S + 0.3 / (S * S) +  &
               0.8 * AJ * SQRT(AJ/S) * ABS(S-0.6))
!
! UPWARD COLLISIONAL RATE  (I TO J) : 
! VRIENS AND SMEETS (1980), EQN. (17), EXCEPT LYMAN ALPHA
! LYMAN ALPHA ACCORDING TO JANEV ET AL. (1987) PP. 18-21 AND 257
!
! TEV IS IN EV IN POLYNOMIAL EXPRESSIONS
!
      IF (I.EQ.1 .AND. J.EQ.2) THEN
!
         TEV = LOG(T * BK / EE)
!
         RH1 = B1S2P(1) + TEV * (B1S2P(2) + TEV * (B1S2P(3) + &
                          TEV * (B1S2P(4) + TEV * (B1S2P(5) + &
                          TEV * (B1S2P(6) + TEV * (B1S2P(7) + &
                          TEV * (B1S2P(8) + TEV *  B1S2P(9) )))))))
         RH2 = B1S2S(1) + TEV * (B1S2S(2) + TEV * (B1S2S(3) + &
                          TEV * (B1S2S(4) + TEV * (B1S2S(5) + &
                          TEV * (B1S2S(6) + TEV * (B1S2S(7) + &
                          TEV * (B1S2S(8) + TEV *  B1S2S(9) )))))))
!
         H1BB = EXP(RH1) + EXP(RH2)
!
      ELSE
!
         H1BB = 1.6D-7 * SQRT (EE / HH / CC) * &
                 (AMN * LOG(0.3D0 * ENG / EION + DELMN) + BMN) * &
                 SQRT(ENG) * EXP(-DERG/ENG) / (ENG + GAMMN) 
!
      ENDIF
!
      RETURN
      END Function H1BB
!
!
      FUNCTION H1BF(I,T)
!
!  DETERMINES H COLLISIONAL IONIZATION COEFFICIENTS UP TO N=50
!
! CONVERSION FACTORS:
!
!      CM**-1 TO EV    HH * CC / EE
!      CM**-1 TO K     HH * CC / BK
!
!      K     TO HZ    BK / HH
!      CM**-1 TO HZ    CC
!
! CONSTANTS:
!
!      EE   ELECTRON CHARGE      1.602189E-12         
!      HH   PLANCK CONSTANT      (ERG S)
!      CC   VELOCITY OF LIGHT    (CM S**-1)
!      EM   ELECTRON REST MASS   (G)
!      UU   ATOMIC MASS UNIT     (G)
!      BK   BOLTZMANN CONSTANT   (ERG K**-1)
!      PI   3.14159265359
!
!
! INPUT:
!
!   I        LOWER LEVEL
!   T        TEMPERATURE (K)
!
! OUTPUT:
!
!   H1BF      COLLISIONAL IONIZATION RATE               I TO K
!
!...................................................................
!
!: H1BF   90-02-28  MODIFICATIONS: (MARTIN J. STIFT)
!:        REWRITTEN.  UP TO 50 BOUND LEVELS + ONE CONTINUUM LEVEL
!:
      Implicit None
      INTEGER I,K
      DOUBLE PRECISION &
              T,RATIK,AI,ENG,DION, &
              EION,EV(50),EE,HH,CC,EM,UU,BK,PI,TEV, &
              B1SI(9),B2SI(9),H1BF
!
! H1S TO H+
!
      DATA B1SI /  -3.271396786375D+1, &
                    1.353655609057D+1, &
                   -5.739328757388D+0, &
                    1.563154982022D+0, &
                   -2.877056004391D-1, &
                    3.482559773737D-2, &
                   -2.631976175590D-3, &
                    1.119543953861D-4, &
                   -2.039149852002D-6/
!
! H2S TO H+
!
      DATA B2SI /  -1.973476726029D+1, &
                    3.992702671457D+0, &
                   -1.773436308973D+0, &
                    5.331949621358D-1, &
                   -1.181042453190D-1, &
                    1.763136575032D-2, &
                   -1.616005335321D-3, &
                    8.093908992682D-5, &
                   -1.686664454913D-6/
!
! HYDROGEN MEAN ENERGY LEVELS TAKEN FROM BASHKIN & STONER (1975), PP. 2-3
!
      DATA EV /      0.000,&
                 82259.047,&
                 97492.279,&
                102823.882,&
                105291.646,&
                106632.157,&
                107440.441,&
                107965.049,&
                108324.718,&
                41 * 0.000/
!
! HYDROGEN IONIZATION LEVEL TAKEN FROM BASHKIN & STONER (1975), P. 3
!
      DATA EION / 109678.764/
!
! MISCELLANOUS CONSTANTS
!
      DATA EE / 1.602189E-12/,&
           HH / 6.626176E-27/,&
           CC / 2.99792458E10/,&
           EM / 9.109534E-28/,&
           UU / 1.6605655E-24/,&
           BK / 1.380662E-16/,&
           PI / 3.14159265359/
!
      AI = FLOAT(I)
!
! DETERMINE ENERGY LEVELS UP TO N=50
!
      DO 1 K = 10,50
         EV(K) = EION * (1. - 1. / FLOAT(K*K))
    1 CONTINUE
!
! DION AND ENG  IN   CM**-1
!
! ENG   ENERGY ASSOCIATED TO TEMPERATURE
! DION  IONISATION ENERGY FROM LOWER LEVEL
!
      ENG  = T * BK / HH / CC
      DION = EION  - EV(I)
!
! COLLISIONAL IONIZATION RATE  (I TO K) : VRIENS AND SMEETS (198), EQN. (8)
!
      IF (I.EQ.1) THEN
!
         TEV = LOG(T * BK / EE)
         RATIK = B1SI(1) + TEV * (B1SI(2) + TEV * (B1SI(3) + &
                           TEV * (B1SI(4) + TEV * (B1SI(5) + &
                           TEV * (B1SI(6) + TEV * (B1SI(7) + &
                           TEV * (B1SI(8) + TEV *  B1SI(9) )))))))
         H1BF = EXP(RATIK)
!
      ELSE IF (I.EQ.2) THEN
!
         TEV = LOG(T * BK / EE)
         RATIK = B2SI(1) + TEV * (B2SI(2) + TEV * (B2SI(3) + &
                           TEV * (B2SI(4) + TEV * (B2SI(5) + &
                           TEV * (B2SI(6) + TEV * (B2SI(7) + &
                           TEV * (B2SI(8) + TEV *  B2SI(9) )))))))
         H1BF = EXP(RATIK)
!
      ELSE
!
         H1BF = 9.56D-6 * EXP(-DION/ENG) / ENG / SQRT(ENG) * &
           SQRT(EE / HH / CC) * (EE / HH / CC) / &
          ( (DION/ENG)**2.33 + 4.38*(DION/ENG)**1.72 + 1.32*(DION/ENG) )
!
      ENDIF
!
      RETURN
      END Function H1bf

!
!
      Subroutine ENEQ(NK, NDEP, A4, B4, NEWMAT)                           
!
! Adapted from MULTI
!                                                                       
!  SOLVES THE EQUATION  SYSTEM  A*X=B.                                  
!  WHEN NEWMAT=TRUE, THE SYSTEM IS REARRANGED INTO U*X=L*B, WHERE U     
!  IS UPPER AND L IS LOWER TRIANGULAR. THESE ARE THEN REUSED IN LATER   
!  CALLS WITH NEWMAT=FALSE AND NEW RIGHT HAND SIDES B. THE SOLUTION     
!  VECTOR IS RETURNED IN B. NO PIVOTING, I.E. THE MATRIX A IS ASSUMED   
!  TO HAVE NONZERO DIAGONAL ELEMENTS.                                   
!                                                                       
!  CODED BY: A. NORDLUND (FEB-1979)                                     
!                                                                       
!  THIS IS A MODIFIED VERSION OF EQSYST WHICH TESTS FOR ZERO ELEMENTS   
!  BELOW THE DIAGONAL AND ALSO STOPS AT THE LAST NON-ZERO ELEMENT ABOVE 
!  THE DIAGONAL. CONSIDERABLE SAVINGS ARE OBTAINED FOR LOOSE MATRICES.  
!                                                                       
!  THIS IS A COLUMN ORIENTED VERSION (M. CARLSSON JAN-1986)             
!  TEMPORARY SCALARS BL, ALL, ALM AND BK ARE USED TO SHOW THE COMPILER  
!  THAT THERE IS NO VECTOR DEPENDENCY IN THE INNERMOST DO-LOOP          
!                                                                       
!: ENEQ   90-06-08  NEW ROUTINE: (MARTIN J STIFT)                       
!:        LIKE EQSYST BUT FOR A BLOCK DIAGONAL GRAND MATRIX             
!:                                                                      
        Implicit None
        Integer :: NK, NDEP, N, MD, L, K, M
        Real, Dimension(NK, NK, NDEP) :: A4
        Real, Dimension(NK, NDEP) :: B4
        Real (kind=8), Dimension(NK, NK, NDEP) :: A
        Real (kind=8), Dimension(NK, NDEP) :: B
        Integer, Dimension(:,:), Allocatable, Save :: LASTN
        Logical :: NEWMAT
        Logical, SAVE :: Warning=.False.
        Real (kind=8) :: ALL, BL, ALM, BK
!             
        Call time_routine('eneq',.True.)
        If (.not. Allocated(LASTN)) Allocate(LASTN(NK,NDEP))
        A=A4
        B=B4
!                                

        N = NK                                                             
!       
        Do MD=1, NDEP
!                                                                       
         If (NEWMAT) then
!                                                                       
!     FIND THE LAST NON-ZERO ELEMENT IN EACH COLUMN                     
!                                                                       
            Do L=1, N                                               
               Do K=N, L+1, -1                                   
                  If (A(K, L, MD) .ne. 0.0) GOTO 20                      
               End do
               K=L                                                     
   20          LASTN(L, MD)=K                                         
            End do
!                                                                       
!  COLUMN LOOP: ELIMINATE ELEMENTS BELOW THE DIAGONAL IN COLUMN L.      
!                                                                       
            Do L=1, N - 1                                           
!                                                                       
!  STORE -A(K,L)/A(L,L) IN ELEMENT A(K,L)                               
!  MULTIPLY RIGHT HAND SIDE WITH -A(K,L)/A(L,L)                         
!                    
               ALL = A (L, L, MD)  
               If (ALL .eq. 0) then
                  If (.not. Warning) Then
                     Print *,'Warning in ENEQ, diagonal element of A is zero'
                     Print *,'(trying to recover...)'
                     Print *,'Point:',MD
!                     Warning=.True.
                  End if
                  ALL=1.
                  A(L, L, MD)=ALL
               End if
                                      
               BL = B (L, MD)                                            
               Do K = L + 1, LASTN (L, MD)                            
                  A (K, L, MD) = - A (K, L, MD) / ALL                    
                  B (K, MD) = B (K, MD) + A (K, L, MD) * BL              
               End do
!                                                                       
!  ADD FRACTION -A(K,L)/A(L,L) OF ROW L TO ROW K.                       
!  IN EACH COLUMN GO THROUGH ALL ROWS                                   
!                                                                       
               Do M = L + 1, N                                        
                  If (A (L, M, MD) .ne. 0.0) then
                     ALM = A (L, M, MD)                                  
                     LASTN (M, MD)=MAX(LASTN (L, MD) , LASTN (M, MD))
                     Do K = L + 1, LASTN (L, MD)                      
                        A (K, M, MD) = A (K, M, MD) + A (K, L, MD)      &
                        * ALM                                           
                     End do
                  End if
               End do
!                                                                       
            End do ! Do in L
!                                                                       
         Else ! if .not. NEWMAT
!                                                                       
            Do L = 1, N - 1                                           
               BL = B (L, MD)                                            
               Do K = L + 1, LASTN (L, MD)                            
                  B (K, MD) = B (K, MD) + A (K, L, MD) * BL              
               End do
            End do
!                                                                       
         End if ! End If NEWMAT
!                                                                       
! BACKSUBSTITUTE                                                        
!                                                                       
         Do K = N, 1, - 1                                            
            BK = B (K, MD)                                               
            Do L = K + 1, N                                          
               BK = BK - A (K, L, MD) * B (L, MD)                        
            End do
            B (K, MD) = BK / A (K, K, MD)                                
         End do
!                                                                       
      End do ! Do in MD
      B4=B
      A4=A
!                                                                       
      Call time_routine('eneq',.False.)
      Return
    End Subroutine ENEQ
! This subroutine takes a structure Atom and computes the fixed
! radiative rates (Radiation temperatures).
!
Subroutine FixedRad(NLTE, Atom)
  Implicit None
  Type (NLTE_Atom) :: Atom
  Type (NLTE_variables) :: NLTE
  Integer :: K, i, j, itran
  Real :: T, TR, Delt, XXI, CLAM, CA0, RIJK, RJIK, STIMD, EX, X2, &
       HKTT
!
  NLTE%F(:,:,:)=0.
  Do itran=1, Atom%NFIX
     j=Atom%JFIX(itran)
     i=Atom%IFIX(itran)
     XXI=Atom%eV(J)-Atom%eV(I)
     CLAM=hh*cc/ee/XXI
     Do K=1, NLTE%NDEP
        T=NLTE%Atmo%Temp(K)
        DELT=-1.0
        If (K .ne. NLTE%NDEP) DELT=NLTE%Atmo%Temp(K)-NLTE%Atmo%Temp(K+1)
        TR=T
        If (Atom%ITRad(itran) .eq. 2) then
           If (DELT .ge. 0.) TR=Atom%TRad(itran)
           If (DELT .lt. 0. .and. Atom%TRad(itran) .gt. T) &
                TR=Atom%TRad(itran)
        Else if (Atom%ITRad(itran) .eq. 3 .and. Atom%TRad(itran) .lt. T &
             .and. DELT .ge. 0.) then
           TR=Atom%TRad(itran)
        Endif
!
        If (Atom%IPho(itran) .eq. 1) then
           CA0=8.*Pi*cc/CLAM/CLAM/CLAM*Atom%A0(itran)
           RIJK=CA0*Rint(ee/bk/TR*XXI)
           RJIK=RIJK
           If (TR .ne. T) RJIK=CA0*RINT(ee/bk/T*XXI)
           HKTT=ee/bk*XXI/(T*TR)
           STIMD=RINT1(T,TR,HKTT)*CA0
           RJIK=(RJIK+STIMD)*NLTE%NSTAR(i,K)/NLTE%NSTAR(j,K)
        Else
           EX=EXP(-ee/bk/TR*XXI)
           X2=(EE*XXI/HH)**2
           RIJK=7.421e-22*Atom%A0(itran)*X2*EX/(1.-EX)
           RJIK=RIJK
           If (TR .ne. T) then
              EX=EXP(-ee/bk/T*XXI)
              RJIK=7.421e-22*Atom%A0(itran)*X2*EX/(1.-EX)
           End if
           RJIK=RJIK*Atom%g(i)/(EX*Atom%g(j))
        End if
        NLTE%F(i,j,K)=NLTE%F(i,j,K)+RIJK
        NLTE%F(j,i,K)=NLTE%F(j,i,K)+RJIK
     End do
  End do
  Return
!
End Subroutine FixedRad

Function Exint1(X)
!
!  CALCULATES THE EXPONENTIAL INTEGRAL FUNCTION E1(X)
!
!  REFERENCE: HANDBOOK OF MATHEMATICAL FUNCTIONS,
!  NBS APPLIED MATHEMATICS SERIES 55 (1964) (ED. ABRAMOVITZ
!  AND STEGUN), EQUATIONS 5.1.53 AND 5.1.54.
!
  Implicit None
  Real, Dimension(6) :: A
  Real, Dimension(4) :: B
  Real :: Exint1, S, Y, X
  Integer :: i
  Data A/-0.57721566,0.99999193,-0.24991055,0.05519968, &
       -0.00976004,0.00107857/
  Data B/2.334733,0.250621,3.330657,1.681534/
!
  If (X .gt. 80.0) then
     Exint1=0.
  Else if (X .gt. 1.0) then
     S=X*X+B(1)*X+B(2)
     S=S/(X*X+B(3)*X+B(4))
     Exint1=S*EXP(-X)/X
  Else
     S=0.0
     Y=1.0
     Do i=1,6
        S=S+A(I)*Y
        Y=Y*X
     End do
     Exint1=S-LOG(X)
  End if
  Return
End Function Exint1

Function Rint(X)
!
!  EVALUATES THE INTEGRAL NEEDED IN CALCULATING THE FIXED RAD. RATES
!  BY FINDING THE SUM OF THE EXPONENTIAL INTEGRALS.
!
  Implicit None
  Real :: Sum, X, XX, DSUM, Rint
  Integer :: I
!
  SUM=Exint1(X)
  If (SUM .eq. 0.0) then
     Rint=0.
     Return
  Else
     Do I=2,999
        XX=I*X
        DSUM=Exint1(XX)
        SUM=SUM+DSUM
        If (DSUM/SUM .le. 1.E-5) then
           Rint=SUM
           Return
        End if
     End do
  End if
  Rint=SUM
  Return
End Function Rint

Function Rint1(TE, TR, HKTT)
  Implicit None
  Real :: TE, TR, HKTT, SUM, XX, DSUM, Rint1
  Integer :: i
!
  SUM=Exint1((TE+TR)*HKTT)
  If (SUM .eq. 0.0) then
     Rint1=SUM
     Return
  Else
     Do i=2, 999
        XX=(i*TE+TR)*HKTT
        DSUM=EXINT1(XX)
        SUM=SUM+DSUM
        If (DSUM/SUM .le. 1.e-5) then
           Rint1=SUM
           Return
        End if
     End do
  End if
  Rint1=SUM
  Return
End Function Rint1
!
Subroutine FormalSolution(NLTE, imu, inu, itran, iformal, X, S, RNu, P, LStMuNu, UseLinear)
! Given the opacity and source function, this routine solves the
! scalar radiative transfer equation. It returns P (the average
! incoming and outgoing intensity) and LStMuNu (the monochromatic local 
! Lambda-operator.
!
! X(1:NDEP)=Total opacity/Continuum opacity at 500nm
! S(1:NDEP)=Source function (arbitrary units)
! RNu(1:NDEP)=Ratio of line to total opacity (needed to compute LStMuNu)
! IncidentInt=Incident Intensity (same units as S)
! P(1:NDEP)=Output (same units as S)
! LStMuNu(1:NDEP)=Output (dimensionless)
!
  Implicit None
  Type (NLTE_variables) :: NLTE
  Real :: OptThin, OptThick
  Integer :: imu, inu, itran, idepth, BoundUp, BoundLow, iformal, ii
  Real, Dimension(0:NLTE%NDEP) :: Iminus, Iplus
  Real, Dimension(NLTE%NDEP) :: X, S, P, LStMuNu, Tau_500, Tau_Nu, Dtau_Nu,&
       RNu
  Real :: IncidentInt, CMu, T
  Logical :: CheckNaN, UseLinear(NLTE%NDEP)
  Real :: u
!
  OptThin=NLTE%OptThin
  OptThick=NLTE%OptThick
  Iminus(:)=0.
  Iplus(:)=0.
  BoundUp=1
  BoundLow=NLTE%NDEP
!
! Compute monochromatic optical depths
  CMu=1./NLTE%XMu(imu)
  Tau_500(:)=10**(NLTE%Atmo%ltau_500)
  T=Tau_500(1)*X(1)*CMu
  Tau_Nu(1)=T
  Dtau_Nu(1)=T

!  Dtau_Nu(2:NLTE%NDEP)=.5*(X(2:NLTE%NDEP)+X(1:NLTE%NDEP-1))* &
!       (Tau_500(2:NLTE%NDEP)-Tau_500(1:NLTE%NDEP-1))*CMu

!  call bezier_qintep(NLTE%ndep, tau_500, X, dtau_nu) 
!  dtau_nu = dtau_nu * cmu

  call get_dtau

  Do idepth=2, NLTE%NDEP
     Tau_Nu(idepth)=Tau_Nu(idepth-1)+DTau_Nu(idepth)
  End do
  ! First point of the grid
  tau_nu(1) = exp(2.0d0 * log(tau_nu(2)) - log(tau_nu(3)))
  dtau_nu(1) = tau_nu(2) - tau_nu(1)

  Do idepth=2, NLTE%NDEP
     If (Tau_Nu(idepth) .gt. OptThin .and. Tau_Nu(idepth-1) .lt. OptThin) &
          BoundUp=idepth-1
     If (Tau_Nu(idepth) .gt. OptThick .and. Tau_Nu(idepth-1) .lt. OptThick) &
          BoundLow=idepth
  End do

  BoundLow=BoundLow+1
  If (BoundLow .gt. NLTE%NDEP) BoundLow=NLTE%NDEP
  BoundUp=BoundUp-1
  If (BoundUp .lt. 1) BoundUp=1

!  Call Check_tau(Tau_Nu)

  If (iformal .eq. 1) then
     Call Bezier_NLTE
 Else
     Call ShortCharacteristics
  End If
!
  Return
!
  contains

 subroutine bezier_nlte
    
      Implicit none
      integer :: k0, k1, dk, k
      real(8) :: rlu(nlte%ndep), ex, ex1, t, sgrad, dtau, der, der1
      logical :: stopme
      
      !
      ! Init variables
      !
      iplus(:)=0.
      iminus(:)=0.
      p(:)=0.
      Rlu(1:boundup)=0.
      Rlu(boundlow:NLTE%ndep)=1.
      Rlu(boundup:boundlow) = 1.0d0 - Rnu(boundup:boundlow)
      !
      LStMuNu = 0.d0
      !
      T = Tau_Nu(BoundUp)
      
      If (T .lt. 1e-3) then ! Optically thin
         Ex1=T*(1.d0-T*(0.5d0 - T * ((1.0d0 / 6.0d0) - T / 24.d0)))
      Else
         if (T .lt. 100) then
            Ex1 = 1.0d0 - Exp(-T)
         Else
            Ex1 = 1.0d0
         End if
      End if
      !
      Ex = 1.0d0 - Ex1
      
      !
      ! Fill-in the transparent part
      !
      Iminus(0)=0.0d0
      Iminus(1:BoundUp)= Iminus(0)!S(1:BoundUp)*(T-0.5d0*T**2)
      
      !
      ! Boundary conditions at the bottom
      !
         do k = boundlow, NLTE%ndep
!            If (k .gt. 1 .and. k .lt. NLTE%NDEP) then
               ! sgrad=(S(k+1)-S(k-1)) &
               !      /(DTau_Nu(k+1)+DTau_Nu(k))
            
               !
!               der = (s(k) - s(k-1)) / dtau_nu(k)
!               der1 = (s(k + 1) - s(k)) / dtau_nu(k+1)
!               If(der*der1 .lt. 0.0d0) then
!                  sgrad = 0
!               Else
!                  !
!                  sgrad = (der * dtau_nu(k+1) + der1 * dtau_nu(k)) / (dtau_nu(k+1) + dtau_nu(k))
!               endif
!            Else 
!               sgrad=(S(k)-S(k-1))/DTau_nu(k)
!            End if
            Iplus(k) = S(k) ! + sgrad * NLTE%XMu(imu)
            Iminus(k) = Iplus(k)
         end do

      If (BoundLow .ge. 2) then ! If all is not optically thick
         
         ! Upwards integration
         k0 = boundlow
         k1 = boundup
         dk = -1
         if (UseLinear(k0)) then ! debug . Entire ray is linear
            call linear_ray(k0, k1, dk)
         Else
            call bezier_ray(k0, k1, dk)
         End if
         
         ! Downwards integration
         k0 = boundup
         k1 = boundlow
         dk = 1
         if (UseLinear(k0)) then ! debug . Entire ray is linear
            call linear_ray(k0, k1, dk)
         Else
            call bezier_ray(k0, k1, dk)
         End if
      end if
      
      !
      ! Compute contribution to J and normalize lambda operator
      !
      Iplus(0:Boundup) = Iplus(boundup)
      !
      do k = 1, NLTE%ndep
         P(k) = (Iplus(k) + Iminus(k)) * 0.5d0
         LStMuNu(k) = LStMuNu(k) * 0.5d0
      end do
      
      !
      do k = boundlow, NLTE%ndep-1
         LStMuNu(k) = 1.0d0
      end do
      LStMuNu(BoundLow:NLTE%NDEP)=1.d0-RNU(BoundLow:NLTE%NDEP)
      LStMuNu(NLTE%ndep) = 0.5d0*(1.0d0-RNU(NLTE%NDEP))
      Iplus(0) = 2.0d0 * (Ex * P(1) + S(1) * 0.5d0 * Ex1**2)
      
      do k = boundup, boundlow
       !  print *, rlu(k),  lstmunu(k)
         lstmunu(k) = min(1.0 , lstmunu(k))
         lstmunu(k) = max(1.e-12, lstmunu(k))
      end do
      !
 end subroutine bezier_nlte
    subroutine get_dtau
      Implicit none
      integer :: k, ku, kd
      real(8) :: xu, xd, x0, dx, dx1, der, der1
      real(8) :: ctu, ctd, lambda, cmu03, dt, dt1
      real(8) :: xp(nlte%ndep)
      
      dtau_nu(2:nlte%ndep) = tau_500(2:nlte%ndep)-tau_500(1:nlte%ndep-1)
      der = (x(2) - x(1)) / dtau_nu(2)
      xp(1) = der

      cmu03 = 1.d0 / (3.d0 * NLTE%XMu(imu))

      do k = 2, nlte%ndep - 1
         ku = k-1
         kd = k+1
         dt = dtau_nu(k)
         dt1 = dtau_nu(kd)
         der1 = (x(kd) - x(k)) / dt1
         if(der*der1 .gt. 0.0d0) then
            lambda = (1.0+dt1/(dt+dt1))/3.0d0
            xp(k) = DER/(LAMBDA*DER1+(1.0d0-LAMBDA)*DER)*DER1
            !xp(k) = (der * dt1 + der1*dt) / (dt+dt1)
         else
            xp(k) = 0.0d0
         end if
         
         der = der1
      end do
      xp(nlte%ndep) = der
      
      do k = 2,nlte%ndep
         ku = k-1
         kd = k+1
         dt = dtau_nu(k)
         xu = x(ku)
         x0 = x(k)
         !
         if(k .eq. nlte%ndep) then
            ctu =  xu + 0.5d0 * dt * xp(ku)
         elseif(k .eq. 2) then
            ctu = x0 - 0.5d0 * dt * xp(k)
         else
            ctu = ( (x0 - 0.5d0 * dt * xp(k))  +  (xu + 0.5d0 * dt * xp(ku)) ) * 0.5d0
         endif

         !if(ctu .gt. max(x0, xu) .or. (ctu .lt. min(x0,xu)) .or. (ctu .lt. 0.d0)) ctu = xu

         dtau_nu(k) = dt * (ctu + x0 + xu) * cmu03
      enddo
      dtau_nu(1) = 0.d0
    end subroutine get_dtau
  subroutine linear_ray(k0, k1, dk)
      
      Implicit None
      
      integer :: k, kd, ku, dk, k0, k1, kend
      real(8) :: u0, u1, u2
      real(8) :: dt, dt1, dtau, dt2, dt3
      real(8) :: c0, cu, eps, cmu05 ,rlu, dI
      real(8) :: xu, xd, x0, dx, dx1, der, der1, der_save
      real(8) :: ctu, ctd, xp, lambda, cmu03

       ! Ray integration
      !
      cmu05 = 0.5d0  / NLTE%XMu(imu)
      cmu03 = 1.0d0 / (3.0d0 * NLTE%XMu(imu))
      !
      !call get_dtau(tau_500, x, dt, k0+dk, dk)
      !dt = dt / NLTE%XMu(imu)

      dt1 = cmu05 * (x(k0+dk) + x(k0)) * abs(Tau_500(k0+dk) - Tau_500(k0))

      do k = k0 + dk, k1, dk
         ku = k - dk
        ! kd = k + dk

         if(dk .eq. -1) then
            dt = dtau_nu(k+1)
         !   dt1 = dtau_nu(k)
         else
            dt = dtau_nu(k)
        !    dt1 = dtau_nu(k+1)
         end if

         dt2 = dt**2
         dt3 = dt2 * dt
         ! Introduce U{0,1} to compute the interpolation coeff.
         if(dt .gt. 1e-2) then 
            if(dt .lt. 100) then
               eps = exp(-dt)
            else
               eps = 0.0d0
            end if
            u0 = 1.0d0 - eps
            u1 = dt - u0
         else
            dt3 = dt2 * dt
            dt2 = dt * dt
            eps = 1.0d0 - dt + dt2 * 0.5d0 - dt3 / 6.0d0 !+ dt3*dt/24.0d0
            u1 = 0.5d0 * dt2 - dt3 / 6.0d0 !- (dt3*dt) / 24.0d0
            u0 = dt - u1
         endif
         
         ! Interpolation coeff.
         c0 = u1 / dt
         cu = u0 - c0
         dI = cu * s(ku) + c0 * s(k) 
         
         if(dk .eq. 1) then
            iminus(k) = iminus(ku) * eps + dI
         else
            iplus(k) = iplus(ku) * eps + dI
         endif

         lstmunu(k) = LStMuNu(k) + c0 * (1.d0-rnu(k))
         
         ! Reuse value 
      end do

    end subroutine linear_ray
    subroutine bezier_ray(k0, k1, dk)
      
      Implicit None
      
      integer :: k, kd, ku, dk, k0, k1, kend
      real(8) :: u0, u1, u2
      real(8) :: dt, dt1, dtau, dt2, dt3, rlu, dt33
      real(8) :: S0, Su, Sd, sp, ome
      real(8) :: cp, mmi, mma, dI
      real(8) :: c0, cu, cd, eps, cmu05, cp1
      real(8) ::  xu, xd, x0, dx, dx1, der, der1, der_save
      real(8) :: ctu, ctd, xp, lambda, cmu03, dlst
      logical :: bound

      ! Reach boundary of the array? Do linear integration there

      if((k1 .eq. nlte%ndep) .OR. (k1 .eq. 1)) then
         bound = .true.
         kend = k1 - dk
      Else
         bound = .false.
         kend = k1 
      end if

      !
      ! Ray integration
      !
!      cmu05 = 0.5d0  / NLTE%XMu(imu)
      cmu03 = 1.d0 / (3.0d0 * NLTE%XMu(imu))
      !
      do k = k0 + dk, kend, dk
         ku = k - dk
         kd = k + dk
         
         ! Source function values
         S0 = s(k)
         Su = s(ku)
         Sd = s(kd)
         
         if(dk .eq. -1) then
            dt = dtau_nu(k+1)
            dt1 = dtau_nu(k)
         else
            dt = dtau_nu(k)
            dt1 = dtau_nu(k+1)
         end if

         dt2 = dt**2    
         
       
         !
         ! Introduce formulation based on the control point
         !
         if(dt .ge. 5.d-3) then
            if(dt .gt. 100) then 
               eps = 0.0d0
            else
               eps = exp(-dt)
            endif
            c0 = 1.0d0 - 2.0d0 * (dt + eps - 1.0d0) / dt2
            cu = (2.0d0 - (2.0d0 + 2.0d0*dt + dt2)*eps) / dt2
            cd = 2.d0*(dt - 2.0d0 + (dt + 2.0d0)*eps) / dt2
           ! ome = 1.d0 + 2.d0 * (eps * (1.d0 - dt) - 1.d0) / dt2
         Else
            dt3 = dt * dt2
            eps = 1.d0 - dt + 0.5d0 * dt2 - dt3/6.d0 + dt2*dt2/24.d0 - dt3*dt2/120.d0 + dt3*dt3/720.d0

            c0 = -dt * (-120 + 30*dt - 6*dt2 + dt3 + 2*dt2*dt2) / 360.d0
            
            cu = -dt *(-60 + 45*dt - 18*dt2 + 5*dt3) / 180.d0
            
            cd = dt*(120 + (-2 + dt)*dt*(30 + dt*(6 + 5*dt))) / 360.d0
         endif
         
         der = (S0-Su) / dt
         der1= (Sd-S0) / dt1


         if(der*der1 .le. 0.0d0) then
            !Max/Min in the source function
            cp = S0
         Else
            !sp = (dt*der1 + dt1 * der) / (dt + dt1)
            !double lambda = (1.0 + dzd / (dzd + dzu)) / 3.0;
            ! dki = (deu / (lambda * ded + (1.0 - lambda) * deu)) * ded;
            lambda = (1.0d0 + dt1 / (dt + dt1)) / 3.0d0;
            sp = (der / (lambda * der1 + (1.0 - lambda) * der)) * der1;
            cp = S0 - 0.5d0 * dt * sp
            !
           ! if((cp .le. max(S0,Su)) .AND. (cp .ge. min(S0,Su))) then  ! Normal case
           !    cp1 = S0 + 0.5d0 * dt1 * sp 
           !    if(cp1 .gt. max(S0,Sd) .or. (cp1 .lt. min(S0,Sd))) then ! Overshooting in downwind interval
           !       !
           !       sp = 2.d0 * (Sd - S0)/dt1
           !       cp = S0 - dt * 0.5d0 * sp
           !    endif
           ! Else ! Overshooting in the upwind interval
           !    cp = Su
           ! endif
         endif

         dI = S0 * c0 + Su * cu + cp * cd

         ! Integrate and define lambda operator
         if(dk .eq. 1) then
            iminus(k) = iminus(ku) * eps + dI
         else
            iplus(k) = iplus(ku) * eps + dI
         endif
         
         ome = c0 + cd

         LStMuNu(k) = LStMuNu(k) + ome * (1.d0-rnu(k))
      end do

      ! Do last point of the array? Use linear
      if(bound) call linear_ray(k1-dk, k1, dk)

    end subroutine bezier_ray
!
    Subroutine ShortCharacteristics
!
! FORMAL SOLVER USING SHORT CARACTERISTICS (SEE GORDON, OLSON, KUNASZ
! 1987). 
!
      Implicit None
      Integer :: step, K, KINIT, KEND
      Real :: T, Ex, Ex1, E0, E1, E2, D2, D3, D4, &
           ALF, BET, GAM, DELTAI, SDELTA, DELTAL, sgrad
      Double Precision :: DTM, DTP, EXPDTP, EXPDTM
      Logical, Save :: Warning=.False.
!
      LStMuNu(1:NLTE%NDEP)=0.
!
! It could happen that DTP was not defined and nicole crashed.
! Better to initialize it...
!
      DTP=0.
!
      T=Tau_Nu(BoundUp)
      If (T .lt. 1e-3) then ! Optically thin
         Ex1=T*(1.-T*(.5-T*(0.1666667-T*0.041666667)))
      Else
         if (T .lt. 20) then
            Ex1=1.-Exp(-T)
         Else
            Ex1=1.
         End if
      End if
      Ex=1.-Ex1
      Iminus(0)=0.
!      Iminus(1:BoundUp)=S(1:BoundUp)*(T-.5*T*T)
      Iminus(1:BoundUp)=Iminus(0)
      If (Iminus(1) .lt. 0) Iminus(1)=0.
! Bottom boundary condition: S+dS/dTau_Nu
      Do k=BoundLow, NLTE%NDEP
!         If (k .ge. 2 .and. k .le. NLTE%NDEP-1) then
!            sgrad=(S(k+1)-S(k-1)) &
!                 /(DTau_Nu(k+1)+DTau_Nu(k))
!         Else if (k .le. NLTE%NDEP-1) then
!            sgrad=(S(k+1)-S(k))/DTau_nu(k+1)
!         Else if (k .ge. 2) then
!            sgrad=(S(k)-S(k-1))/DTau_nu(k)
!         End if
         Iplus(K)=S(K) !+sgrad ! Max(sgrad,-S(K)*.9)
         Iminus(K)=Iplus(k)
      End do
      If (BoundLow .ge. 2) then ! Do the transfer if not all optically thick
         Do step=-1,1,2 ! FIRST INWARDS AND THEN OUTWARDS
            If (step .eq. 1) then ! IF INWARDS
               KINIT=BoundUp
               KEND=BoundLow
            Else                ! IF OUTWARDS
               KINIT=BoundLow
               KEND=BoundUp
            End if
            Do K=KINIT+step, KEND-step, step
               IF (step .eq. 1) then
                  DTM=DTau_Nu(K)
                  DTP=DTau_Nu(K+1)
               Else
                  DTM=DTau_Nu(K+1)
                  DTP=DTau_Nu(K)
               End if
               If (DTM .ge. 1e-3) then
                  EXPDTM=Exp(-DTM)
                  E0=1.-EXPDTM
                  E1=DTM-E0
                  E2=DTM*DTM-2.*E1
               Else
                  EXPDTM=1.-DTM+(DTM*DTM)/2. ! TAYLOR EXPANSION
                  D2=DTM*DTM ! TO IMPROVE PRECISION AT SMALL OPTICAL DEPTHS
                  D3=DTM*D2
                  D4=DTM*D3
                  E0=DTM-(D2/2.)
                  E1=(D2/2.)-(D3/6.)
                  E2=(D3/3.)-(D4/12.)
               End if
               ALF=E0+(E2-(DTP+2*DTM)*E1)/(DTM*(DTM+DTP))
               BET=((DTM+DTP)*E1-E2)/DTM/DTP
               GAM=(E2-DTM*E1)/(DTP*(DTM+DTP))
               DELTAI=ALF*S(K-STEP)+BET*S(K)+GAM*S(K+STEP)

               SDELTA=1.-RNu(K)
               DELTAL=BET*SDELTA
               If (step .eq. 1) then
                  Iminus(K)=Iminus(K-1)*EXPDTM+DELTAI
! Valgrind complains about Iminus(2) here
                  IF ( Iminus(K) .lt. 0 .or. DELTAL .gt. 1. .or. UseLinear(K))  then
! COMPUTE AGAIN USING LINEAR INTERPOLATION
!                     print *,'using linear interp 1 at ',k
!                     print *,'dtau=',dtm,dtp
!                     print *,'s=',s(k-1:k+1)
!                     print *,'deltai,deltal=',deltai,deltal,expdtm
                     ALF=E0-E1/DTM
                     BET=E1/DTM
                     DELTAI=ALF*S(K-step)+BET*S(K)
                     Iminus(K)=Iminus(K-1)*EXPDTM+DELTAI
                     GAM=0.
                     LStMuNu(K)=LStMuNu(K)+BET*SDELTA
!!$                     If (Iminus(K) .lt. 0 .or. LStMuNu(K) .lt. 0 .or. &
!!$                          LStMuNu(K) .gt. 2) then
!!$                        Print *,'K=',K
!!$                        Print *,'Iminus=',Iminus(K-1:K)
!!$                        Print *,'LStMuNu=',LStMuNu(K-1:K)
!!$                        Print *,'S(K-1:K)=',S(K-1:K)
!!$                        Print *,'dtm=',dtm
!!$                        If (IMinus(K) .lt. 0) Stop
!!$                     End if
                  Else
                     LStMuNu(K)=LStMuNu(K)+DELTAL
                  End if
               Else
                  Iplus(K)=Iplus(K+1)*EXPDTM+DELTAI

                  IF (Iplus(K) .lt. 0 .or. DELTAL .gt. -1. .or. UseLinear(K) ) then
! COMPUTE AGAIN USING LINEAR INTERPOLATION
!                     print *,'using linear interp 2 at ',k
!                     print *,'dtau=',dtm,dtp
!                     print *,'s=',s(k-1:k+1)
!                     print *,'deltai,deltal=',deltai,deltal,expdtm
!                     stop
                     ALF=E0-E1/DTM
                     BET=E1/DTM
                     DELTAI=ALF*S(K-STEP)+BET*S(K)
                     Iplus(K)=Iplus(K+1)*EXPDTM+DELTAI
                     GAM=0.
                     LStMuNu(K)=LStMuNu(K)+BET*SDELTA
!!$                     If (Iplus(K) .lt. 0 .or. LStMuNu(K) .lt. 0 .or. &
!!$                          LStMuNu(K) .gt. 2) then
!!$                        Print *,'K=',K
!!$                        Print *,'IPlus=',IPlus(K:K+1)
!!$                        Print *,'LStMuNu=',LStMuNu(K:K+1)
!!$                        Print *,'S(K:K+1)=',S(K:K+1)
!!$                        Print *,'dtm=',dtm
!!$                        If (IPlus(K) .lt. 0) Stop
!!$                     End if
                  Else
                     LStMuNu(K)=LStMuNu(K)+DELTAL
                  End if
               End if
            End do ! Do in K
! USE LINEAR INTERPOLATION FOR THE LAST POINT
            IF (DTP .ge. 1e-3) then ! DELTA-TAU BETWEEN THE LAST TWO POINTS
               EXPDTP=Exp(-DTP)        ! IS DTP
               E0=1.-EXPDTP 
               E1=DTP-E0      
            Else
               EXPDTP=1.-DTP+(DTP*DTP)/2. ! TAYLOR EXPANSION 
               E0=DTP-DTP*DTP/2.          ! TO IMPROVE PRECISION AT SMALL
               E1=DTP*DTP/2.              ! OPTICAL DEPTHS
            End if
            K=KEND
            ALF=E0-E1/DTP
            BET=E1/DTP
            DELTAI=ALF*S(K-STEP)+BET*S(K)
            SDELTA=1.-RNu(K)
            LStMuNu(K)=LStMuNu(K)+BET*SDELTA
            IF (Step .eq. 1) then
               Iminus(K)=Iminus(K-1)*EXPDTP+DELTAI
            Else
               Iplus(K)=Iplus(K+1)*EXPDTP+DELTAI
            End if
         End do ! Do in step
      End if ! End check if all is optically thick

!      print *,'iplus top=',iplus(1:4)
!      print *,'lst top=',lstmunu(1:4)
!      print *,'iminus bottom=',iminus(boundlow-3:boundlow)
!      print *,'lst bottom=',lstmunu(boundlow-3:boundlow)
!      stop
      If (boundup .ge. 2) &
           Iplus(1:BoundUp)=Iplus(BoundUp)
      P(1:NLTE%NDEP)=(Iplus(1:NLTE%NDEP)+Iminus(1:NLTE%NDEP))/2.
      Iplus(0)=2.*(Ex*P(1)+S(1)*0.5*Ex1**2)
      LStMuNu(1:NLTE%NDEP)=LStMuNu(1:NLTE%NDEP)/2.
      LStMuNu(BoundLow:NLTE%NDEP)=1.-RNU(BoundLow:NLTE%NDEP)
      LStMuNu(NLTE%NDEP)=.5*(1.-RNU(NLTE%NDEP))
      If (MinVal(P) .lt. 0. .or. &
           MinVal(LStMuNu) .lt. 0. .or. &
           MaxVal(LStMuNu) .gt. 1.) then
         If (.not. Warning) then
            Print *,'WARNING: Negative I or LStar <0 or >1'
            Warning=.True.
         End if
         Do K=1, NLTE%NDEP
            If (P(K) .lt. 0) then
               P(K)=1.e-10
            endif
            If (LStMuNu(K) .lt. 0.) LStMuNu(K)=1.e-10
            If (LStMuNu(K) .gt. 1.) LStMuNu(K)=0.999
         End do
      End if


!      print *,'SCiplus=',iplus(:)
!      print *,'SCiminus=',iminus(:)
!      print *,'SClambda=',lstmunu(:)
!      stop
      Return

    End Subroutine ShortCharacteristics

End Subroutine FormalSolution
!
! Uses Gaussian quadrature for angular integrations
!
Subroutine AngQuad(NLTE,NLTEInput)
  Implicit None
  Type (NLTE_variables) :: NLTE
  Type (NLTE_input) :: NLTEInput
!
  If (NLTEInput%NMU .gt. 5) then
     Print *,'Too many angles in angular quadrature:',NLTEInput%NMU
     Print *,'This version supports only up to NMU=5'
     Print *,'Should be easy to modify gaussquad and'
     Print *,'angquad to support more'
     Stop
  End if
  NLTE%XMu(1:NLTEInput%NMU)=XMU_Gauss_Quad(NLTEInput%NMU,1:NLTEInput%NMU)
  NLTE%WMu(1:NLTEInput%NMU)=WMU_Gauss_Quad(NLTEInput%NMU,1:NLTEInput%NMU)
!
End Subroutine AngQuad
!
! This routine computes quadrature points and weights for line and
! continuum transitions. In the case of continuum transitions, Atom%AlphaC
! is also filled assuming a nu^-3 dependence (unless Q0 or QMAX .lt. 0, in 
! which case it is explicitly given in the file ATOM).
!
! For cont transitions, QMAX is the shortest wavelength considered (in Angs)
! For line transitions, QMAX is the farthest wavelength (in Doppler units)
!
Subroutine FreqQuad(Atom, NLTEInput)
  Implicit None
  Type (NLTE_Atom) :: Atom
  Type (NLTE_input) :: NLTEInput
  Real, Dimension (0:MaxNFreqs) :: X
  Integer :: itran, ix, ncore, nwing
  Real :: DX
!
! Line transitions
!
  Do itran=1, Atom%NLIN
     If (Atom%QMax(itran) .ge. 0) then ! If Qmax < 0, Q has been read from file
        ncore=(Atom%NQ(itran)-1)/2 ! One point is for line center (Q=0)
        nwing=Atom%NQ(itran)-ncore
        if (Atom%QMax(itran) .eq. Atom%Q0(itran)) then ! All is equispaced
           ncore=Atom%NQ(itran)
           nwing=0
        End if
!
! Symmetric profile, consider only 1 side
!
! core
        Atom%Q(1,itran)=0. ! line center
        If (ncore .gt. 1) then
           DX=Atom%Q0(itran)/(ncore-1)
           Do ix=1, ncore-1
              Atom%Q(ix+1,itran)=ix*DX
           End do
        End if
! wing
        if (nwing .gt. 0) then
           DX=( Log10(Atom%QMAX(itran))-Log10(Atom%Q0(itran)) )/ & 
                (nwing)
           Do ix=1, nwing
              Atom%Q(ncore+ix,itran)=10.**(Log10(Atom%Q0(itran))+(ix)*DX)
           End do
        End if
     End if ! End if in Qmax .ge. 0
!
! Asymmetric profile, consider both sides
!
     If (.not. NLTEInput%VelFree) then
        X(1:Atom%NQ(itran))=Atom%Q(1:Atom%NQ(itran),itran)
        Do ix=1, Atom%NQ(itran)/2
           Atom%Q(ix,itran)=-X(Atom%NQ(itran)+2-ix*2)
           Atom%Q(Atom%NQ(itran)/2+ix,itran)=X(ix*2)
        End do
     End if
!
! Weights
!
     Atom%WQ(1,itran)=0.5*(Atom%Q(2,itran)-Atom%Q(1,itran))
     Atom%WQ(Atom%NQ(itran),itran)=0.5*(Atom%Q(Atom%NQ(itran),itran)- &
          Atom%Q(Atom%NQ(itran)-1,itran))
     Do ix=2, Atom%NQ(itran)-1
        Atom%WQ(ix,itran)=0.5*(Atom%Q(ix+1,itran)-Atom%Q(ix-1,itran))
     End do
     If (NLTEInput%VelFree) Atom%WQ(2:,itran)=2.*Atom%WQ(2:,itran)
  End do
!
! Continuum transitions
!
 Do itran=Atom%NLIN+1, Atom%NLIN+Atom%NCNT
    X(0)=cc/Atom%Alamb(itran)*1.e8 ! Freq in Hz
    If (atom%QMAX(itran) .ge. 0.) then ! Compute frequency grid
       X(Atom%NQ(itran))=cc/Atom%QMAX(itran)*1.e8 ! Hz
       DX=(X(Atom%NQ(itran)) - X(0))/(Atom%NQ(itran)-1.0)
       Do ix=1, Atom%NQ(itran)
          X(ix)=X(0)+DX*(ix-1)
       End do
       Do ix=1, Atom%NQ(itran)
          Atom%Q(ix,itran)=cc*1e-5/NLTEInput%QNORM*(X(ix)/X(0)-1.)
          Atom%AlphaC(ix,itran)=Atom%f(itran)*((X(0)/X(ix))**3)
       End do
       Atom%FRQ(0:Atom%NQ(itran),itran)=X(0:Atom%NQ(itran))
    Else ! Grid was read from ATOM file
       Atom%FRQ(0,itran)=cc/Atom%ALAMB(itran)*1.E8
       Atom%FRQ(1:Atom%NQ(itran),itran)=cc/Atom%Q(1:Atom%NQ(itran),itran)*1.e8
       Atom%Q(1:Atom%NQ(itran),itran)=cc*1.e-5/NlteInput%QNORM* &
            (Atom%FRQ(1:Atom%NQ(itran),itran)/Atom%FRQ(0,itran)-1.)
    End if
    Atom%WQ(1,itran)=0.5*(Atom%Q(2,itran)-Atom%Q(1,itran))
    Atom%WQ(Atom%NQ(itran),itran)=0.5*(Atom%Q(Atom%NQ(itran),itran)- &
         Atom%Q(Atom%NQ(itran)-1,itran))
    Do ix=2, Atom%NQ(itran)-1
       Atom%WQ(ix,itran)=0.5*(Atom%Q(ix+1,itran)-Atom%Q(ix-1,itran))
    End do
!
 End do
!

End Subroutine FreqQuad

SUBROUTINE NG(AN,NK,NDEP,NGON, Printout)
!
!  PERFORMS NG ACCELERATION
!  SEE L.H. AUER , P. 101
!  IN KALKOFEN, ED., "NUMERICAL RADIATIVE TRANSFER", 
!  CAMBRIDGE UNIVERSITY PRESS 1987
!
!: NG     90-09-11  NEW ROUTINE: (MARTIN J. STIFT, MATS CARLSSON)
!:
!
!
      Real, DIMENSION(NK,NDEP) :: AN
      LOGICAL :: NGON, Printout
!
      Real, DIMENSION(3,3) :: AA
      Real, Dimension(3) :: BB
      Real, Dimension(:,:,:), Allocatable, Save :: YS
      Integer, Save :: ICALL
      DATA ICALL/0/
!
!  IF NGON=.FALSE., RESET ITERATION COUNT AND RETURN
!
      If (.not. Allocated(YS)) Allocate(YS(5,NK,NDEP))
      IF(.NOT.NGON) THEN
        ICALL=0
        RETURN
      ENDIF
      ICALL=ICALL+1
!
! STORE ITERATED LEVEL POPULATIONS FOR NG ACCELERATION
!
      IS0 = MOD(ICALL-1,5) + 1
      DO 910 I = 1,NK
        DO 900 K = 1,NDEP
          YS(IS0,I,K) = (AN(I,K))
  900   END DO
  910 END DO
!
      IF (IS0.EQ.5) THEN
!
         If (Printout) &
              Print *,' NG acceleration'

        DO 960 I = 1,NK
!
          DO 930 K1 = 1,3
            DO 920 K2 = 1,3
              AA(K1,K2) = 0.
  920       END DO
            BB(K1) = 0.
  930     END DO
!
          DO 940 K = 1,NDEP
            WT = 1. / YS(5,I,K)**2
            D0 = YS(5,I,K) - YS(4,I,K)
            D1 = D0 - YS(4,I,K) + YS(3,I,K)
            D2 = D0 - YS(3,I,K) + YS(2,I,K)
            D3 = D0 - YS(2,I,K) + YS(1,I,K)
            AA(1,1) = AA(1,1) + WT * D1 * D1
            AA(1,2) = AA(1,2) + WT * D1 * D2
            AA(1,3) = AA(1,3) + WT * D1 * D3
            AA(2,2) = AA(2,2) + WT * D2 * D2
            AA(2,3) = AA(2,3) + WT * D2 * D3
            AA(3,3) = AA(3,3) + WT * D3 * D3
            AA(2,1) = AA(1,2)
            AA(3,1) = AA(1,3)
            AA(3,2) = AA(2,3)
            BB(1)   = BB(1)   + WT * D0 * D1
            BB(2)   = BB(2)   + WT * D0 * D2
            BB(3)   = BB(3)   + WT * D0 * D3
  940     END DO
!
          CALL LINEQ (AA,BB,3,3)

          DO 950 K = 1,NDEP
            AN(I,K) = (1. - BB(1) - BB(2) - BB(3)) * YS(5,I,K) + &
                     BB(1) * YS(4,I,K) +  &
                     BB(2) * YS(3,I,K) +  &
                     BB(3) * YS(2,I,K)
!            AN(I,K)=Exp(AN(I,K))
  950     END DO
!
  960   END DO
ENDIF
!
RETURN
End SUBROUTINE NG

SUBROUTINE LINEQ(A,B,N,M)
!
!  FINDS SOLUTION OF SYSTEM OF LINEAR EQUATIONS
!  WITH GAUSSIAN ELIMINATION WITH PIVOTING
!
!: LINEQ  90-06-05  NEW ROUTINE: (MARTIN J STIFT)
!:
!:        09-05-01  MODIFICATIONS: (MATS CARLSSON)
!:        MADE ROUTINE GENERAL AND NOT SPECIALIZED TO 3X3 MATRICES
!:
!     INCLUDE 'PREC'
!
      PARAMETER (MDIM=10000)
      DIMENSION A(M,M),B(M),C(MDIM),ICOL(MDIM)
!
      IF(N.GT.MDIM) THEN
         PRINT *,'LINEQ: N.GT.MDIM'
         STOP
      ENDIF
!
! INITIALIZE COLUMN COUNT AND STARTING POINTS
!
      DO 10 I = 1,N
        ICOL(I) = I
   10 CONTINUE
!
      IBEG = 1
      JBEG = 1
!
! DETERMINE PIVOT
!
   20 AMA = ABS(A(IBEG,JBEG))
      IMA = IBEG
      JMA = JBEG
      DO 30 I = IBEG,N
        DO 30 J = JBEG,N
          IF(ABS(A(I,J)).LE.AMA) GOTO 30
          AMA = ABS(A(I,J))
          IMA = I
          JMA = J
   30 CONTINUE
!
! ORDER MATRIX DEPENDING ON PIVOT
!
      DO 40 I = 1,N
        TEMP = A(I,JMA)
        A(I,JMA) = A(I,JBEG)
        A(I,JBEG) = TEMP
   40 CONTINUE
!
      DO 50 J = JBEG,N
        TEMP = A(IMA,J)
        A(IMA,J) = A(IBEG,J)
        A(IBEG,J) = TEMP
   50 CONTINUE
!
      TEMP = B(IMA)
      B(IMA) = B(IBEG)
      B(IBEG) = TEMP
      IT = ICOL(JBEG)
      ICOL(JBEG) = ICOL(JMA)
      ICOL(JMA)= IT
      IBEG = IBEG + 1
      JBEG = JBEG + 1
!
! ELIMINATE
!
      IMIN = IBEG - 1
      JMIN = JBEG - 1
      DO 70 I = IBEG,N
        QUOT = A(I,JMIN) / A(IMIN,JMIN)
        DO 60 J = JBEG,N
          A(I,J) = A(I,J) - QUOT * A(IMIN,J)
   60   CONTINUE
        B(I) = B(I) - QUOT * B(IMIN)
   70 CONTINUE
      IF (IBEG.LT.N) GOTO 20
!
! DETERMINE COEFFICIENTS
!
      B(N) = B(N) / A(N,N)
      N1 = N - 1
      DO 90 I = N1,1,-1
        I1 = I + 1
        DO 80 J = N,I1,-1
          B(I) = B(I) - A(I,J) * B(J)
   80   CONTINUE
        B(I) = B(I) / A(I,I)
   90 CONTINUE
!
! REORDER COEFFICIENTS
!
      DO 100 I = 1,N
        C(ICOL(I))=B(I)
  100 CONTINUE
!
      DO 110 I = 1,N
        B(I) = C(I)
  110 CONTINUE
!
      RETURN
END SUBROUTINE LINEQ
!
! This subroutine reads the file NLTE_lines with information to match
! the lines in the LINES file to the NLTE transitions in the ATOM file.
! The format is:
!      iline     itran      ratio_l     ratio_u
! where iline is the index of a line in the LINES file, itran is the
! corresponding transition in the ATOM file and ratio_l is the ratio
! of lower level populations of the line divided by the level population
! in the ATOM file. Same for ratio_u and the upper level.
!
! The information is then used to fill in the corresponding fields in the 
! Line structure
!
Subroutine Read_NLTE_lines(Params, NLTEInput, NLTE, Atom, Line)
  Implicit None
  Type (Parameters) :: Params
  Type (Line_data), Dimension (Params%n_lines) :: Line
  Type (NLTE_Atom) :: Atom
  Type (NLTE_input) :: NLTEInput
  Type (NLTE_variables) :: NLTE
  Integer :: i2, iline, ix
  Real :: ratiol, ratiou
  Character (len=256) :: String
!
! Initialize
!
  Line(:)%NLTEgridsize=1
  Do iline=1, Params%n_lines
     If (Allocated(Line(iline)%NLTEgrid)) &
          Deallocate(Line(iline)%NLTEgrid)
     Allocate (Line(iline)%NLTEgrid(1))
     Line(iline)%NLTEgrid(1)=0.
     If (Allocated(Line(iline)%NLTESource_f)) &
          Deallocate(Line(iline)%NLTESource_f)
     If (Allocated(Line(iline)%NLTESource_l_f)) &
          Deallocate(Line(iline)%NLTESource_l_f)
     Allocate (Line(iline)%NLTESource_f(NLTE%NDEP,1))
     Allocate (Line(iline)%NLTESource_l_f(NLTE%NDEP))
     Line(iline)%NLTESource_f(NLTE%NDEP,1)=0.
     Line(iline)%NLTESource_l_f(NLTE%NDEP)=0.
     i2=Line(iline)%NLTEtransition
     If (i2 .ge. 1) then
        Line(iline)%NLTEgridsize=Atom%NQ(i2)
        If (Allocated(Line(iline)%NLTEgrid)) &
             Deallocate(Line(iline)%NLTEgrid)
        Allocate(Line(iline)%NLTEgrid(Atom%NQ(i2)))
        If (Allocated(Line(iline)%NLTESource_f)) &
             Deallocate(Line(iline)%NLTESource_f)
        If (Allocated(Line(iline)%NLTESource_l_f)) &
             Deallocate(Line(iline)%NLTESource_l_f)
        Allocate(Line(iline)%NLTESource_f(NLTE%NDEP,Atom%NQ(i2)))
        Allocate(Line(iline)%NLTESource_l_f(NLTE%NDEP))
        Line(iline)%NLTEgrid=Atom%Q(1:Atom%NQ(i2),i2)*NLTEInput%QNORM* &
             1.e5/cc*Atom%Alamb(i2) ! NLTE grid in A from line center
        Line(iline)%NLTESource_f(:,:)=NLTE%Source_f(:,1:Atom%NQ(i2),i2) ! cgs units
        Line(iline)%NLTESource_l_f(:)=NLTE%Source_l_f(:,i2) ! cgs units
     End if
  End do
!
  Return
End Subroutine Read_NLTE_lines
! Compute Voigt profiles Phi and integrated profile WPhi for each transition
! Store the resulting Phi in virtual file PHI and WPhi in the NLTE structure
!
Subroutine VoigtProfs(NLTE, NLTEInput, Atom)
  Implicit None
  Type (NLTE_Atom) :: Atom
  Type (NLTE_variables) :: NLTE
  Type (NLTE_input) :: NLTEInput
  Type (Line_data) :: Line
  Real, Dimension (NLTE%NDEP) :: Phi
  Real :: WQMU, xlam, DlDop, DlDopcms, DlDopkms, Damp, H, F, dqterm
  Integer :: itran, imu, inu, idepth, nwave, itrm, ktrm
!
  Call RewindVirtualFile('PHI')
  NLTE%WPhi(:,:)=0.
  !
  Do itran=1, Atom%NLIN
     Line%Wlength=Atom%Alamb(itran)
     Line%ion_stage=Atom%ion(Atom%i(itran))
     Line%Atomic_number=Atom%Z
     Line%Energy_low=Atom%eV(Atom%i(itran))
     Line%VDW_enh=Atom%GW(itran)
     Line%collisions=1
     Do inu=1, Atom%NQ(itran)
        Do imu=1, NLTEInput%NMU
           WQMU=Atom%WQ(inu, itran)*NLTE%WMu(imu)
           Do idepth=1, NLTE%NDEP
              DlDop=Atom%Alamb(itran)/cc* &
                   Sqrt(2.*bk*NLTE%Atmo%Temp(idepth)/Atom%AWgt/mass_pr+ &
                   NLTE%Atmo%v_mic(idepth)**2.)  ! Angstroms ! debug
              DlDopcms=DlDop/Atom%Alamb(itran)*cc ! DlDop in cm/s
              DlDopkms=DlDopcms*1e-5 ! DlDop in km/s
              If (.not. NLTEInput%VelFree .or. imu .eq. 1) then
                 Call Damping(Line, NLTE%Atmo%Temp(idepth), &
                      NLTE%Atmo%El_p(idepth), NLTE%Atmo%Gas_P(idepth), &
                      DlDop, Damp, Atom%GA(itran), Atom%GQ(itran))
                 If (Atom%Nterm(itran) .le. 1) then
                    xlam=Atom%Q(inu,itran)*NLTEInput%QNORM/DlDopkms
                    xlam=xlam-NLTE%XMu(imu)*NLTE%Atmo%v_los(idepth)/&
                         DlDopcms ! Wlength in Doppler units (dimensionless)
                    H=voigt(Damp, xlam, 0)
                    Phi(idepth)=H/(DlDopkms/NLTEInput%QNORM*SqrtPi)
                 Else
                    Phi(idepth)=0.
                    Do itrm=1,Atom%Nterm(itran)
                       KTRM=Atom%KTERM(itrm,itran)
                       dqterm=Atom%DETERM(ktrm)*Atom%alamb(itran)*1.e-13*cc/NLTEInput%QNORM
                       xlam=(Atom%Q(inu,itran)-dqterm)*NLTEInput%QNORM/DlDopkms
                       xlam=xlam-NLTE%XMu(imu)*NLTE%Atmo%v_los(idepth)/&
                            DlDopcms ! Wlength in Doppler units (dimensionless)
                       Phi(idepth)=Phi(idepth)+Atom%Wterm(ktrm)*H
                    End do
                    Phi(idepth)=H/(DlDopkms/NLTEInput%QNORM*SqrtPi)
                 End if ! Atom%Nterm .le. 1
              End if
              NLTE%WPhi(idepth,itran)=NLTE%WPhi(idepth,itran)+ &
                   WQMU*Phi(idepth)
           End do ! End idepth do
           If (.not. NLTEInput%VelFree .or. imu .eq. 1) &
                Call WriteP(PHI)
        End do ! End imu do
     End do ! End inu do
  End do ! End itran do
!
  NLTE%Wphi(:,1:ATOM%NLIN)=1./NLTE%Wphi(:,1:ATOM%NLIN)

!
  Call RewindVirtualFile('PHI')
  Return
!
End Subroutine VoigtProfs
!
! Compute LTE level populations using information in the ATOM file
! Note: This routine requires a valid Atom%ne array, so a call to
! hydrostatic before this is a good idea.
!
Subroutine LTE_pop_2(NLTE, Atom)
!
  Implicit None
  Type (NLTE_Atom) :: Atom
  Type (NLTE_variables) :: NLTE
  Real, dimension (NLTE%NDEP) :: N_e, U0, U1, U2, U_ion
  Real, dimension (NLTE%NDEP) :: N_ion, N0_N1, N1_N2, N_tot, N0, N1, N2
  Integer :: index, npoints, idepth, iline, Atomic_number, ilev
  Real :: du0, du1, du2, E_ioniz1, E_ioniz2, E_exc_i, E_exc_j, Abund_mass
  Real :: E_ion, E_exc, E_base
  Logical, Save :: Warning=.FALSE.
! Interface for the function Saha
!  Interface
!     Function Saha(npoints, T, Ne, U1, U2, E_ioniz)
!       Real, dimension (npoints) :: T, Ne, U1, U2, Saha
!       Real :: E_ioniz
!       Integer :: npoints
!     End Function Saha
!  End interface
!
  npoints=NLTE%NDEP
  Atomic_number=Atom%Z
  Abund_mass=Atom%AWgt*(10.**(Atom%Abund-12.))/ &
       Sum(At_weight(1:N_elements)*(10.**(At_abund(1:N_elements)-12.)))
  N_tot(1:npoints)=Abund_mass*NLTE%Atmo%rho(1:npoints)/Atom%Awgt/mass_pr
  N_e(1:npoints)=NLTE%Atmo%ne
!
  Do idepth=1, npoints
     Call Partition_f(Atomic_number, NLTE%Atmo%Temp(idepth), &
          U0(idepth), U1(idepth), U2(idepth), du0, du1, du2)
  End do
  E_ioniz1=At_ioniz1(Atomic_number)*eV_to_cgs
  E_ioniz2=At_ioniz2(Atomic_number)*eV_to_cgs
!
! Saha ionization equations
!
  N0_N1(1:npoints)=Saha(npoints, NLTE%Atmo%Temp, N_e, U0, U1, E_ioniz1)
  N1_N2(1:npoints)=Saha(npoints, NLTE%Atmo%Temp, N_e, U1, U2, E_ioniz2)
  N1(1:npoints)=N_tot(1:npoints)/  &
       (1.+N0_N1(1:npoints)+1./N1_N2(1:npoints)) ! N_1 (cm^-3)
  N0(1:npoints)=N0_N1(1:npoints)*N1(1:npoints) ! N_0 (cm^-3)
  N2(1:npoints)=N1(1:npoints)/N1_N2(1:npoints) ! N_2 (cm^-3)
!

  Do ilev=1, Atom%NK
     If (Atom%ion(ilev) .eq. 1) then
        N_ion(1:npoints)=N0(1:npoints)
        U_ion(1:npoints)=U0(1:npoints)
        E_ion=0
     Else if (Atom%ion(ilev) .eq. 2) then
        N_ion(1:npoints)=N1(1:npoints)
        U_ion(1:npoints)=U1(1:npoints)
        E_ion=E_ioniz1
     Else if (Atom%ion(ilev) .ge. 3) then
        N_ion(1:npoints)=N2(1:npoints)
        U_ion(1:npoints)=U2(1:npoints)
        E_ion=E_ioniz1+E_ioniz2
     End if
     E_base=0.
     If (Atom%ion_stage .eq. 2) then
        E_base=At_ioniz1(Atomic_number)
     Else if (Atom%ion_stage .eq. 3) then
        E_base=At_ioniz2(Atomic_number)
     End if
     If (Atom%ion(ilev) .gt. 3 .and. .not. Warning) then
        Print *,'Ionization stage gt 3 in LTE_pop_2'
        Stop
     End if
     E_exc=( (Atom%eV(ilev)+E_base)*eV_to_cgs - E_ion)
     If (E_exc .lt. 0) then 
        If (E_exc .gt. -1e-18) then ! Numerical accuracy
           E_exc=0.
        Else
           Print *,'Negative excitation potential in LTE_pop_2'
           Stop
        End if
     End if
     NLTE%Nstar(ilev,1:npoints)=N_ion(1:npoints)*Atom%g(ilev)* & ! Boltzmann
          exp(-E_exc/bk/NLTE%Atmo%Temp(1:npoints))/U_ion(1:npoints)
  End do

! Check that the sum of all levels doesn't exceed totn
  Do idepth=1, NLTE%NDEP
     If (Sum(NLTE%Nstar(:,idepth))/N_Tot(idepth) .gt. 1.1) then
        Print *,'Error in NLTE module, routine ltepop2'
        Print *,'At depth point:',idepth
        Print *,'The total population of all levels is ', &
             Sum(NLTE%Nstar(:,idepth))
        Print *,'According to abundance data it should be .lt. ',N_Tot(idepth)
        Print *,'This might indicate an error in the G values in the ATOM file'
        Print *,'Trying to recover...'
        If (idepth .ge. 2) then 
           NLTE%Nstar(:,idepth)=NLTE%Nstar(:,idepth-1)/ &
                N_Tot(idepth-1)*N_Tot(idepth)
        Else
           NLTE%Nstar(1,idepth)=N_Tot(idepth)
           NLTE%Nstar(2:,idepth)=0.
        End if
     End if
  End do
!

  Return
End Subroutine LTE_pop_2

Subroutine LTE_pop_3(NLTE, Atom)
!*==LTEPOP.spg  processed by SPAG 6.70Dc at 16:16 on 22 May 2013
      Use Phys_constants
      IMPLICIT NONE
      Type (NLTE_Atom) :: Atom
      Type (NLTE_variables) :: NLTE
!*--LTEPOP4
!*** Start of declarations inserted by SPAG
      REAL ccon , conl , dxi , ek, em , evi ,   &
         & sumn , t , tns , tnsl , TOTN
      INTEGER i , k , k1 , k2 , l , ljoblo , ndxi 

!*** End of declarations inserted by SPAG
!
!  CALCULATES LTE POPULATIONS
!:
!: LTEPOP 87-04-07  MODIFICATIONS: (PHILIP JUDGE)
!:        A DANGER SIGN IS OUTPUT WHEN ZERO POPULATIONS ARE FOUND.
!:
!:        88-02-01  MODIFICATIONS: (MATS CARLSSON)
!:        PROVISION FOR MOLECULES ADDED
!:
!:        88-07-01  MODIFICATIONS: (MATS CARLSSON)
!:        DEBYE LOWERING OF IONIZATION POTENTIAL TAKEN INTO ACCOUNT
!:
!:        92-09-10  MODIFICATIONS: (MATS CARLSSON)
!:        DANGER SIGN LIMIT CHANGED FROM 1.E-37 TO 0.0
!:
!:        98-09-04  MODIFICATIONS: (MATS CARLSSON)
!:        QCALC REPLACED BY QPART FOR CO, CH, CN
!:
!:        00-09-12  MODIFICATIONS: (ANA)
!:        ADDED OH MOLECULE
!:
!
!  PGJ ADDED FOLLOWING LINE
      CHARACTER (len=4) :: elemid
!
      LOGICAL molec
      integer :: npoints, atomic_number
      real :: abund_mass, n_e, n_tot
      DIMENSION tns(atom%nk), ndxi(20), totn(nlte%ndep), n_e(nlte%ndep), n_tot(nlte%ndep)
!*
!* 88-07-01 MATS CARLSSON
!* THE ENERGIES OF IONIZATION ARE REDUCED BY ION*DXI, FOLLOWING BASCHEK ET
!* AL., ABH. HAMB. VIII, 26 EQ. (10).
!* SINCE ALL EV REFERS TO THE LOWEST LEVEL, IT IS NECESSARY TO SUBTRACT
!* FROM THE ENERGY OF A LEVEL ALL DXI EXPERIENCED BY LOWER IONIZATION
!* STAGES. THIS NUMBER IS STORED IN NDXI(ION)
!*

      em=9.109534D-28
      npoints=NLTE%NDEP
      Atomic_number=Atom%Z
      Abund_mass=Atom%AWgt*(10.**(Atom%Abund-12.))/ &
           Sum(At_weight(1:N_elements)*(10.**(At_abund(1:N_elements)-12.)))
      N_tot(1:npoints)=Abund_mass*NLTE%Atmo%rho(1:npoints)/Atom%Awgt/mass_pr
      N_e(1:npoints)=NLTE%Atmo%ne
      totn(:)=n_tot(:)

      ndxi(1) = 0
      DO i = 2 , 20
         ndxi(i) = ndxi(i-1) + (i-1)
      ENDDO
      DO i = 1 , 20
         ndxi(i) = ndxi(i) - ndxi(ATOM%ION(1))
      ENDDO
!
      elemid = ATOM%element
      ek=ee/bk
      molec = elemid.EQ.'ch' .OR. elemid.EQ.'co' .OR.                   &
            & elemid.EQ.'cn' .OR. elemid.EQ.'oh'
      ccon = 0.5*(hh/SQRT(2.*pi*em)/SQRT(bk))**3
      DO k = 1 , npoints
         t = nlte%atmo%TEMP(k)
         conl = LOG(ccon*nlte%atmo%NE(k)) - 1.5*LOG(t)
         dxi = 4.98E-4*5040./t*SQRT(nlte%atmo%NE(k)*bk*t)
         sumn = 1.
         DO i = 2 , atom%nk
            evi = atom%EV(i) - dxi*ndxi(ATOM%ION(i))
            tnsl = LOG(Atom%G(i)) - LOG(Atom%G(1)) - ek/t*evi
            IF ( ATOM%ION(i).GT.ATOM%ION(1) ) THEN
               l = ATOM%ION(i) - ATOM%ION(1)
               tnsl = tnsl - FLOAT(l)*conl
            ENDIF
            tns(i) = EXP(tnsl)
            sumn = sumn + tns(i)
         ENDDO
         IF ( .NOT.molec ) THEN
            NLTE%NSTAR(1,k) = TOTN(k)/sumn
         ELSE
            print *,'molecules not implemented in lte_pop_3 (NLTE.f90)'
            stop
!            NLTE%NSTAR(1,k) = TOTN(k)*Atom%G(1)*EXP(-ek/t*atom%EV(1))/QPART(elemid,k)
         ENDIF
         DO i = 2 , atom%nk
            NLTE%NSTAR(i,k) = tns(i)*NLTE%NSTAR(1,k)
         ENDDO
!
!  SCALE TOTAL ABUNDANCES IF TREATED MODEL IS A MOLECULE
!
         IF ( molec ) TOTN(k) = sumn*NLTE%NSTAR(1,k)


      ENDDO

End Subroutine LTE_pop_3

!  
! Solve statistical equilibrium equations and compute NLTE popluations
! Use NLTE%N as starting guess (probably initialized to NLTE%NStar) and
! return the result in the same array
!
Subroutine SolveStat(NLTE, NLTEInput, Atom)
  Implicit None
  Type (NLTE_Atom) :: Atom
  Type (NLTE_variables) :: NLTE
  Type (NLTE_variables), Save :: Saved_NLTE 
  Type (NLTE_input) :: NLTEInput
  Real, Dimension (Atom%NK, NLTE%NDEP) :: E, NOld, Save_N
  Real, Dimension (Atom%NK, Atom%NK, NLTE%NDEP) :: LU
  Real, Dimension (NLTE%NDEP) :: Phi, Xcont, Scat, Sc, Jnu, Gij, Z, &
       LStar, PhiJ, SBck, Alpha, X, Rnu, S, P, &
       LStMuNu, Jeff
  Real, Dimension (0:NLTE%NDEP) :: iminus, iplus
  Integer, Dimension (NLTE%NDEP) :: largepop
  Logical, Dimension (NLTE%NDEP) :: Depth_converged
  logical :: checknan
  Real :: popmax, TotN, RelChg, Gijk, hnu3c2, hc4p, WQMu, IncidentInt, XC
  Real :: CRSW = 1.0
  Integer :: i, j, iter, itran, idepth, imu, inu, irec, unit
  Character (Len=15) :: String
  Character (Len=256) :: Message
  Logical :: Converged, Newmat, Cont, resetNG
  Logical, Save :: FirstTime=.True.
  Logical, Dimension (NLTE%NDEP) :: UseLinear
! Do we need to redo the whole calculation?
  Call time_routine('solvestat',.True.)
  UseLinear(:)=.False.
  If (NLTE%Linear .eq. 1) UseLinear(:)=.True.
  Debug_errorflags(flag_NLTE)=0
  Debug_warningflags(flag_NLTE)=0

  If (FirstTime) then
     Call Allocate_model(Saved_NLTE%Atmo, NLTE%NDEP)
     Saved_NLTE%Atmo%temp(:)=0.
     Saved_NLTE%Atmo%v_los(:)=0.
     Saved_NLTE%Atmo%v_mic(:)=0.
     Allocate(Saved_NLTE%N(Atom%NK,NLTE%NDEP))
     Saved_NLTE%N(:,:)=0.
     Allocate(Saved_NLTE%Source_f(NLTE%NDEP,MaxNFreqs,Atom%NLIN+Atom%NCNT))
     Saved_NLTE%Source_f(:,:,:)=0.
     Allocate(Saved_NLTE%Source_l_f(NLTE%NDEP,Atom%NLIN+Atom%NCNT))
     Saved_NLTE%Source_l_f(:,:)=0.
  Else
     If (MaxVal(Abs(NLTE%atmo%temp-Saved_NLTE%atmo%temp)) .lt. 1e-5 &
          .and. ( MaxVal(Abs(NLTE%atmo%v_los-Saved_NLTE%atmo%v_los)) .lt. 1e-5 &
                   .or. NLTEInput%VelFree ) &
          .and. MaxVal(Abs(NLTE%atmo%v_mic-Saved_NLTE%atmo%v_mic)) .lt. 1e-5 &
          .and. .not. NLTE%Error ) & 
     then 
        If (NLTEInput%Verbose .ge. 3) &
             Print *,'Using previous populations'
        Call time_routine('solvestat',.False.)
        NLTE%N(:,:)=Saved_NLTE%N(:,:)
        NLTE%Source_f(:,:,:)=Saved_NLTE%Source_f(:,:,:)
        NLTE%Source_l_f(:,:)=Saved_NLTE%Source_l_f(:,:)
        NLTE%Error=Saved_NLTE%Error
        NLTE%atmo=Saved_NLTE%atmo
        Return
     End if
  End if

!
! First calculation of radiation field, needed to initialize Jnu for Sbck
  Relchg=1e10
  Call Radiation
  FirstTime=.False.
  Converged=.False.
  Depth_converged(:)=.False.
  Newmat=.True.
  iter=0

  If (NLTEInput%MaxIters .lt. 1) then
     Converged=.True.
  End if

  Do while (iter .lt. NLTEInput%MaxIters .and. .not. Converged)
     iter=iter+1
     NLTE%Error=.False.
 ! Initialize arrays
     E(:,:)=0.
     If (Newmat) then
        NLTE%W(:,:,:)=0.
     End if
! Rewind virtual files
     Call RewindVirtualFile('PHI')
     Call RewindVirtualFile('OPC')
     irec=0
! Compute terms for line transitions
     Do itran=1, Atom%NLIN+Atom%NCNT
        i=Atom%i(itran)
        j=Atom%j(itran)
        LStar(:)=0.
        PhiJ(:)=0.
! Source functions, opacities and cross-sections
        Cont=.False.
        If (itran .gt. Atom%NLIN) then
           Cont=.TRUE.
           NLTE%WPHI(:,itran)=1.0
        Else
           Gijk=Atom%g(i)/Atom%g(j)
           hnu3c2=Atom%A(itran)/Atom%B(j,i) ! It's actually 2*h*nu^3/c^2
           Z(:)=NLTE%N(i,:)-Gijk*NLTE%N(j,:)
           Z(:)=Max(Z(:),Z(:)*0.) ! Guarantee positivity (no masers)
           Gij(:)=Gijk
           NLTE%Sl(:,itran)=hnu3c2*NLTE%N(j,:)*Gijk/(Z(:)+1e-20) ! 1e-20 is the
                                                     ! minimum line opac allowed
        End if
        Do inu=1, Atom%NQ(itran)
           If (Cont .or. inu .eq. 1) then
              Call ReadX(Xcont, Scat, Sc)
              LStar(:)=0.
           End if

           hc4p=hh*cc*1.e-5/NLTEInput%QNORM/4./Pi
           irec=irec+1
           Call ReadJ(irec, Jnu)
           SBck(:)=Sc(:)+Scat(:)*Jnu(:)
           If (Debug_level .ge. 1) then
              If (checknan(Sum(SBck))) then
                 Do idepth=1, NLTE%NDEP
                    If (CheckNan(SBck(idepth))) then
                       Write (Message,*) 'SBck .lt. 0 in solvestat. ip, itran, S, Sc=', &
                            idepth,itran,S(idepth),Scat(idepth)
                       Call Debug_Log(Message,1)
                       Debug_errorflags(flag_NLTE)=1
                    End if
                 End do
              End if
           End if
           Jnu(:)=0. ! Reset Jnu to recompute at this frequency
           If (Cont) then
              Do idepth=1, NLTE%NDEP ! Do not add twice the same continuum
                 XC=XCont(idepth)
                 XC=MAX(1.e-3*XCont(idepth),XCont(idepth)- &
                      Z(idepth)*Alpha(idepth)/NLTE%Xnorm(idepth))
                 XCont(idepth)=XC
              End do
              hnu3c2=2.*hh*Atom%FRQ(inu,itran)/cc*Atom%FRQ(inu,itran)/cc * &
                   Atom%FRQ(inu,itran)
              Gij(:)=NLTE%NStar(i,:)/NLTE%NStar(j,:)*Exp(-hh* &
                   Atom%FRQ(inu,itran)/bk/NLTE%Atmo%Temp(:))
              Alpha(:)=Atom%AlphaC(inu,itran)
              Z(:)=NLTE%N(i,:)-Gij(:)*NLTE%N(j,:)
              Z(:)=Max(Z(:),Z(:)*0.) ! Guarantee positivity (no masers)
              NLTE%Sl(:,itran)=hnu3c2*NLTE%N(j,:)*Gij(:)/(Z(:)+1e-20) ! 1e-20 is
                                                  !the minimum line opac allowed
           End if
           Do imu=1, NLTEInput%NMU
              WQMu=Atom%WQ(inu,itran)*NLTE%WMu(imu)
              If (.not. Cont .and. (.not. NLTEInput%VelFree .or. imu .eq. 1)) then
                 Call ReadP(Phi)
                 Alpha(:)=Atom%B(i,j)*Phi(:)*hc4p*NLTE%WPhi(:,itran)
              End if
              If (.not. NLTEInput%VelFree .or. imu .eq. 1) then
                 X(:)=Z(:)*Alpha(:)/NLTE%Xnorm(:)+Xcont(:)
!                 Rnu(:)=Xcont(:)/X(:)
                 Rnu(:)=Xcont(:)/(Z(:)*Alpha(:)/NLTE%Xnorm(:)+Xcont(:))
                 S(:)=(1.-Rnu(:))*NLTE%Sl(:,itran)+Rnu(:)*SBck(:)

              End if
              IncidentInt=0. ! No incident radiation
! Compute monochromatic radiation field and local operator
              If (imu .eq. 1) NLTE%Source_f(:,inu,itran)=S(:)             
!              If (imu .eq. 1) NLTE%Source_f(:,inu,itran)=NLTE%Sl(:,itran)
              If (imu .eq. 1) NLTE%Source_l_f(:,itran)=NLTE%Sl(:,itran)    

              Call time_routine('NLTE_formalsolution',.True.)
              Call FormalSolution(NLTE, imu, inu, itran, NLTEInput%NLTE_formal_solution, &
                   X, S, RNu, P, LStMuNu, UseLinear)
              Call time_routine('NLTE_formalsolution',.False.)
              If (iter .le. NLTEInput%NumLambdaIters) LStMuNu(:)=0.
              If (Debug_level .ge. 1) then
                 If (checknan(Sum(P))) then
                    Do idepth=1, NLTE%NDEP
                       If (CheckNan(P(idepth))) then
                          Write (Message,*) 'P is NaN in solvestat. ip, itran, S, Scat, P=',&
                               idepth,itran,S(idepth),Scat(idepth)
                          Call Debug_Log(Message,2)
                          Write (Message,*) 'P is NaN in solvestat. X, Sl, SBck, S, Rnu=',&
                               X(idepth),NLTE%SL(idepth,itran),SBck(idepth),S(idepth),RNu(idepth)
                          Call Debug_Log(Message,1)
                          Debug_errorflags(flag_NLTE)=1
                       End if
                    End do
                 End if
              End if
! Initialize with populations from zero radiation field
!   (need to initialize to LTE first) 
              If (NLTEInput%IStart .eq. 0 .and. iter .eq. 1) then
                 P(:)=0.
                 LStMuNu(:)=0.
              End if
!
              Jnu(:)=Jnu(:)+NLTE%WMu(imu)*P(:)

              PhiJ(:)=PhiJ(:)+WQMu*Phi(:)*P(:)*NLTE%WPhi(:,itran)
              If (Cont) then
                 LStar(:)=LStar(:)+NLTE%WMu(imu)*LStMuNu(:)
              Else
                 LStar(:)=LStar(:)+WQMu*LStMuNu(:)*Phi(:)*NLTE%WPhi(:,itran)
              End if
           End do ! Do in imu
!
           If (MaxVal(LStar(:)) .gt. 1.) then
              Do idepth=1, NLTE%NDEP
!                 If (lstar(idepth) .gt. 1.000) then
!                    print *,'lstar .gt. 1',lstar(idepth)
!                    pause
!                 endif
                 LStar(idepth)=1.
              End do
           End if
           If (MinVal(LStar(:)) .lt. 0) then
              Do idepth=1, NLTE%NDEP
                 If (lstar(idepth) .lt. 0) then
                    print *,'lstar .lt. 0'
                 endif
                 If (LStar(idepth) .lt. 0) LStar(idepth)=0.
              End do
           End if

           If (Cont) then ! b-f transitions              
              Jeff(:)=Jnu(:)-LStar(:)*NLTE%Sl(:,itran)*1.00 ! debug
              If (MinVal(Jeff) .lt. 0) then 
                 Do idepth=1, NLTE%NDEP
                    If (Jeff(idepth) .lt. 0) Jeff(idepth)=0.
!                    If (Lstar(idepth) .gt. 1.D0 .or. Lstar(idepth) .lt. 0.D0) print *,'Lst(',idepth,')=',Lstar(idepth)
!                    If (Jeff(idepth) .lt. 0) then
!                       print *,'i,jnu,lst,sl=',idepth,jnu(idepth),lstar(idepth),nlte%sl(idepth,itran)
!                       print *,'jeff .lt. 0'
!                       stop
!                    endif
!                    If (Jeff(idepth) .lt. 0) Jeff(idepth)=Jnu(idepth)
                 End do
              End if

              ! Term with n_u (changed sign)
              NLTE%W(i,j,:)=NLTE%W(i,j,:)+Atom%WQ(inu,itran)/hc4p/ &
                   Atom%FRQ(inu,itran)*Atom%FRQ(0,itran)*Gij(:)*Alpha(:)*&
                   (hnu3c2*(1.-LStar(:)*1.00)+Jeff(:)) ! debug
              ! Term with n_l
              NLTE%W(j,i,:)=NLTE%W(j,i,:)+Atom%WQ(inu,itran)/hc4p/ &
                   Atom%Frq(inu,itran)*Atom%FRQ(0,itran)*Alpha(:)*Jeff(:)
              if (minval(nlte%w(i,j,:)) .lt. 0) then
                 do idepth=1,nlte%ndep
                    if (nlte%w(i,j,idepth) .lt. 0) then
                       print *,'cccccccc'
                       print *,'i,idepth=',i,idepth,' n=',nlte%n(i,idepth)
                       print *,'w=',nlte%w(i,:,idepth)
                       print *,'lst=',lstar(idepth)
                       print *,'jeff=',jeff(idepth)
                       stop
                    endif

                 end do
              endif

           End if
           Call WriteJ(irec,Jnu)
        End do ! Do in inu

        If (.not. Cont) then ! b-b transitions
           Jeff(:)=Phij(:)-LStar(:)*NLTE%Sl(:,itran)
           if (MinVal(Jeff(:)) .lt. 0) then
              Do idepth=1, NLTE%NDEP
                 If (Jeff(idepth) .lt. 0) then
!                    print *,'Jeff .lt. 0'
                    If (Jeff(idepth) .le. 0.1*PhiJ(idepth)) then
                       Jeff(idepth)=0.
                    Else
                       Jeff(idepth)=PhiJ(idepth)
                    End if
                 End if
              End do
           End if
           ! Term with n_u (changed sign)
           NLTE%W(i,j,:)=NLTE%W(i,j,:)+(Atom%A(itran)*(1.-LStar(:))+ &
                Atom%B(j,i)*Jeff(:))
           ! Term with n_l
           NLTE%W(j,i,:)=NLTE%W(j,i,:)+Atom%B(i,j)*Jeff(:)
        End if ! If .not. Cont
     End do ! Do in itran

! Collisional-radiative switching
     CRSW=1.0
     If (NLTEInput%UseColSwitch .ge. 1) then
        CRSW=1.+ 10**(4.-iter/100.)
        If (NLTEInput%Verbose .ge. 3) &
             Print *,'  Collisional switch=',CRSW
     End if
! Add collisional and fixed rates (C and F have inverted indeces wrt W):
!  W is rate into level i from level j, while C is rate from level i to j

     Do i=1, Atom%NK
        Do j=1, Atom%NK
           NLTE%W(i,j,:)=NLTE%W(i,j,:)+NLTE%C(j,i,:)*CRSW+NLTE%F(j,i,:)
        End do
     End do

! The sum over the i-th column gives the total rate out of the i-th level
     Do i=1, Atom%NK
        Do j=1, i-1
           NLTE%W(i,i,:)=NLTE%W(i,i,:)-NLTE%W(j,i,:)
        End do
        Do j=i+1, Atom%NK
           NLTE%W(i,i,:)=NLTE%W(i,i,:)-NLTE%W(j,i,:)
        End do
     End do
!
     largepop(:)=NLTEInput%ISum
     If (NLTEInput%ISum .eq. 0) then
! Find out which level has the largest population at each depth-point
        Do idepth=1, NLTE%NDEP
           popmax=NLTE%N(1,idepth)
           largepop(idepth)=1
           Do i=2, Atom%NK
              If (NLTE%N(i,idepth) .gt. popmax) then
                 popmax=NLTE%N(i,idepth)
                 largepop(idepth)=i
              End if
           End do
        End do
     End if

! Particle conservation equation
     E(:,:)=0.
     Do idepth=1,NLTE%NDEP
        NLTE%W(largepop(idepth),:,idepth)=1.
!        E(largepop(idepth),idepth)=Sum(NLTE%N(:,idepth))
        E(largepop(idepth),idepth)=Sum(NLTE%NStar(:,idepth))
     End do
!
! Ok. Matrices are now ready. Update populations
!

     NOld(:,:)=NLTE%N(:,:)
     Do idepth=1,NLTE%NDEP
        If (.not. Depth_converged(idepth)) then 
!
!        Solve for dn/n_old
           E(:,idepth)=E(:,idepth)-MatMul(NLTE%W(:,:,idepth), &
                NLTE%N(:,idepth))
           Do j=1, Atom%NK
              NLTE%W(:,j,idepth)=NLTE%W(:,j,idepth)*NLTE%N(j,idepth)
           End do
           If (NewMat) then
              Call EneQ(Atom%NK,1,NLTE%W(:,:,idepth),E(:,idepth),.TRUE.)
           Else
              Call EneQ(Atom%NK,1,LU(:,:,idepth),E(:,idepth),.FALSE.)
           End if
           NLTE%N(:,idepth)=NOLD(:,idepth)*(1.+E(:,idepth))
!
!       Solve for n_new/n_old
!           E(:,idepth)=E(:,idepth)
!           Do j=1, Atom%NK
!              NLTE%W(:,j,idepth)=NLTE%W(:,j,idepth)*NLTE%N(j,idepth)
!           End do
!           If (NewMat) then
!              Call EneQ(Atom%NK,1,NLTE%W(:,:,idepth),E(:,idepth),.TRUE.)
!           Else
!              Call EneQ(Atom%NK,1,LU(:,:,idepth),E(:,idepth),.FALSE.)
!           End if
!           NLTE%N(:,idepth)=E(:,idepth)*NLTE%N(:,idepth)
!
!       Solve for n_new 
!           E(:,idepth)=E(:,idepth)
!           Do j=1, Atom%NK
!              NLTE%W(:,j,idepth)=NLTE%W(:,j,idepth)
!           End do
!           If (NewMat) then
!              Call EneQ(Atom%NK,1,NLTE%W(:,:,idepth),E(:,idepth),.TRUE.)
!           Else
!              Call EneQ(Atom%NK,1,LU(:,:,idepth),E(:,idepth),.FALSE.)
!           End if
!           NLTE%N(:,idepth)=E(:,idepth)
!
        End if
     End do
!
     If (CheckNaN(Sum(NLTE%N(:,:)))) then
        If (NLTEInput%Verbose .ge. 2) &
             Print *,'Warning! NaN populations. Resetting to LTE!'
        Do idepth=1, NDEP
           Do i=1, ATOM%NK
              If (CheckNaN(NLTE%N(i,idepth))) &
                   NLTE%N(i,idepth)=NLTE%NStar(i,idepth)
           End do
        End do
        NLTE%Error=.True.
     End if
     If (MinVal(NLTE%N(:,:)).LE.0) then 
        If (MinVal(NLTE%N(:,:)).LE.0) then
           If (NLTEInput%Verbose .ge. 2) &
                Print *,'Warning! Negative populations!'
           NLTE%Error=.True.
           NLTE%N(:,:)=ABS(NLTE%N(:,:))
        End if
        If (MinVal(NLTE%N(:,:)).EQ.0) then
           If (NLTEInput%Verbose .ge. 2) &
                Print *,'Warning! Zero populations!'
           Do i=1,Atom%NK
              Do idepth=1, NLTE%NDEP
                 NLTE%Error=.True.
                 If (NLTE%N(i,idepth) .eq. 0) &
                      NLTE%N(i,idepth)=1.E-10
              End do
           End do
        End if
     End if

     RelChg=MaxVal(Abs((NOLD(:,:)-NLTE%N(:,:))/NLTE%N(:,:)))

     If (CheckNaN(RelChg)) then
        RelChg=0.
        NLTE%Error=.True.
     End if

     Do idepth=1, NLTE%NDEP
        If (MaxVal(Abs((NOLD(:,idepth)-NLTE%N(:,idepth))/ &
             NLTE%N(:,idepth))) .lt. NLTEInput%elim2) then 
            Depth_converged(idepth)=.True. 
!!           This doesn't work because W contains the populations
!!           LU(:,:,idepth)=NLTE%W(:,:,idepth)
        End if
     End do

! NG acceleration
     Save_N=NLTE%N
     If (iter .eq. 1) & ! Reset NG acceleration
          Call NG(NLTE%N, ATOM%NK, NLTE%NDEP, .false., .false.)

     If (iter .gt. 5 .and. NLTEInput%NGacc) then
        Call NG(NLTE%N, ATOM%NK, NLTE%NDEP, .true., NLTEInput%Verbose .ge. 4)
     End if
     If (minval(nlte%n) .lt. 0 .or. checknan(Sum(nlte%n))) then
        If (NLTEInput%Verbose .ge. 3) &
             print *,'Negative or NaN pop after NG!!' 
        NLTE%N=Save_N ! Revert to previous and reset NG
        Call NG(NLTE%N, ATOM%NK, NLTE%NDEP, .false., .false.)
!        NLTE%Error=.True.
     end if
!
     If (NLTEInput%Verbose .ge. 4) &
          Print *,'POPULATIONS RELATIVE CHANGE: ', RelChg
!     If (RelChg .lt. NLTEInput%elim1 .and. NewMat) THEN
!        NewMat=.FALSE.
!        LU(:,:,:)=NLTE%W(:,:,:)
!     End if
! Stop?
     If (RelChg .lt. NLTEInput%elim2) Converged=.TRUE.
!
  End do ! Do while (iterations)
  If (NLTEInput%Verbose .ge. 3) Print *,'NLTE Conv   ',RelChg,iter
  If (NLTEInput%Verbose .ge. 3 .or. (NLTEInput%Verbose .ge. 2 .and. .not. Converged)) &
       Print *,'Populations maximum relative change: ', RelChg
!

  If (.not. Converged) then
     NLTE%Error=.True.
     If (NLTEInput%Verbose .ge. 4) then
        Do idepth=1, NDEP
           If (idepth .eq. 1) &
                Print *,'Conergence analysis:'
           Print *,'i=',idepth,' tau=',NLTE%Atmo%ltau_500(idepth),' relchg=',&
                MaxVal(Abs((NOLD(:,idepth)-NLTE%N(:,idepth))/ &
             NLTE%N(:,idepth))),' level=',MaxLoc(Abs((NOLD(:,idepth)-NLTE%N(:,idepth))/ &
             NLTE%N(:,idepth)))
        End do
     End if
  End if
!
  Do idepth=1, NLTE%NDEP
     If (CheckNaN(NLTE%N(1, idepth))) then
        Debug_errorflags(flag_NLTE)=1
        If (Debug_level .ge. 1) then
           Write (Debug_FileUnit,*) '  ** Error in forward/solvestat'
           Write (Debug_FileUnit,*) '     Atmosphere that caused the error follows'
           Write (Debug_FileUnit,*) '     ------   Z   ------'
           Write (Debug_FileUnit,*) NLTE%Atmo%Z_scale
           Write (Debug_FileUnit,*) '     ------   tau_500  ------'
           Write (Debug_FileUnit,*) NLTE%Atmo%ltau_500
           Write (Debug_FileUnit,*) '     ------   T  ------'
           Write (Debug_FileUnit,*) NLTE%Atmo%Temp
           Write (Debug_FileUnit,*) '     ------   Electron pressure  ------'
           Write (Debug_FileUnit,*) NLTE%Atmo%El_p
           Write (Debug_FileUnit,*) '     ------   Gas pressure  ------'
           Write (Debug_FileUnit,*) NLTE%Atmo%Gas_p
           Write (Debug_FileUnit,*) '     ------   Density  ------'
           Write (Debug_FileUnit,*) NLTE%Atmo%Rho
        End if
        Write (String,*) idepth
        Call Debug_Log('Error!! NLTE calculation returned NaN for populations at i='//Trim(Adjustl(String)),1)
     End if
  End do

  If (RelChg .lt. NLTEInput%elim2 .and. .not. NLTE%Error) then
     Saved_NLTE%N(1:Atom%NK,1:NLTE%NDEP)=NLTE%N(1:Atom%NK,1:NLTE%NDEP)
     Saved_NLTE%Source_f(:,:,:)=NLTE%Source_f(:,:,:)
     Saved_NLTE%Source_l_f(:,:)=NLTE%Source_l_f(:,:)
     Saved_NLTE%Error=NLTE%Error
     !Saved_NLTE%atmo=NLTE%atmo
     Call Model_assign(Saved_NLTE%atmo,NLTE%atmo)
     Saved_NLTE%atmo%temp=NLTE%atmo%temp
     Saved_NLTE%atmo%v_los=NLTE%atmo%v_los
     Saved_NLTE%atmo%v_mic=NLTE%atmo%v_mic
  Else
     Saved_NLTE%atmo%temp(:)=0. ! Make sure we don't reuse these populations
     NLTE%Error=.True.
  End if
!
  If (NLTE%Error) Debug_errorflags(flag_NLTE)=1
!
  If (NLTEInput%Write) then
     Call Open_file_unform(unit, 'ltepop.dat')
     Write (unit) NLTE%NStar(:,:)
     Call Close_file (unit)
     Call Open_file_unform(unit, 'pop.dat')
     Write (unit) NLTE%N(:,:)
     Call Close_file (unit)
  End if
!
  Call time_routine('solvestat',.False.)
  Return
!
  contains
    Subroutine Radiation

      integer :: itran, inu
! Computes Jnu for all transitions, freqs and angles and store it.
! This routine is analogous to the computation in solvestat, except
! that it doesn't compute LStar, Jeff, etc
! Rewind virtual files
      Call RewindVirtualFile('PHI')
      Call RewindVirtualFile('OPC')
      irec=0
! Compute terms for line transitions
      Do itran=1, Atom%NLIN+Atom%NCNT
         i=Atom%i(itran)
         j=Atom%j(itran)
! Source functions, opacities and cross-sections
         Cont=.False.
         If (itran .gt. Atom%NLIN) then
            Cont=.TRUE.
            NLTE%WPHI(:,itran)=1.0
         Else
            Gijk=Atom%g(i)/Atom%g(j)
            hnu3c2=Atom%A(itran)/Atom%B(j,i)
            Z(:)=NLTE%N(i,:)-Gijk*NLTE%N(j,:)
            Z(:)=Max(Z(:),Z(:)*0.) ! Guarantee positivity (no masers)
            Gij(:)=Gijk
            NLTE%Sl(:,itran)=hnu3c2*NLTE%N(j,:)*Gijk/(Z(:)+1e-20) ! 1e-20 is the
                                                     ! minimum line opac allowed
         End if
         Do inu=1, Atom%NQ(itran)
            If (Cont .or. inu .eq. 1) Call ReadX(Xcont, Scat, Sc)
            hc4p=hh*cc*1.e-5/NLTEInput%QNORM/4./Pi
            irec=irec+1
            SBck(:)=Sc(:)
            Jnu(:)=0.
            If (Cont) then
               Do idepth=1, NLTE%NDEP ! Do not add twice the same continuum
                  XC=XCont(idepth)
                  XC=MAX(1.e-3*XCont(idepth),XCont(idepth)- &
                       Z(idepth)*Alpha(idepth)/NLTE%Xnorm(idepth))
                  XCont(idepth)=XC
               End do
               hnu3c2=2.*hh*Atom%FRQ(inu,itran)/cc*Atom%FRQ(inu,itran)/cc * &
                    Atom%FRQ(inu,itran)
               Gij(:)=NLTE%NStar(i,:)/NLTE%NStar(j,:)*Exp(-hh* &
                    Atom%FRQ(inu,itran)/bk/NLTE%Atmo%Temp(:))
               Alpha(:)=Atom%AlphaC(inu,itran)
               Z(:)=NLTE%N(i,:)-Gij(:)*NLTE%N(j,:)
               Z(:)=Max(Z(:),Z(:)*0.) ! Guarantee positivity (no masers)
               NLTE%Sl(:,itran)=hnu3c2*NLTE%N(j,:)*Gij(:)/(Z(:)+1e-20) ! 1e-20
                                             ! is the minimum line opac allowed
            End if
            NLTE%Source_f(:,inu,itran)=S(:)
            Do imu=1, NLTEInput%NMU
               If (.not. Cont .and. (.not. NLTEInput%VelFree .or. imu .eq. 1)) then
                  Call ReadP(Phi)
                  Alpha(:)=Atom%B(i,j)*Phi(:)*hc4p*NLTE%WPhi(:,itran)
               End if
               If (.not. NLTEInput%VelFree .or. imu .eq. 1) then
                  X(:)=Z(:)*Alpha(:)/NLTE%Xnorm(:)+Xcont(:)
                  If (Debug_level .ge. 1) then
                     If (CheckNan(Sum(X)) .or. MinVal(X) .lt. 0) then
                        Do idepth=1, NLTE%NDEP
                           If (CheckNan(X(idepth)) .or. X(idepth) .lt. 0) then
                              Write (Message,*) 'X is NAN or .lt. 0. ip, X, Z, Alpha,XNorm,XCont=', &
                                   idepth,X(idepth),Z(idepth),Alpha(idepth),NLTE%XNorm(idepth),XCont(idepth)
                              Call Debug_Log(Message,1)
                              Debug_errorflags(flag_NLTE)=1
                           End if
                        End do
                     End if
                  End if
!                  Rnu(:)=Xcont(:)/X(:)  
                  Rnu(:)=Xcont(:)/(Z(:)*Alpha(:)/NLTE%Xnorm(:)+Xcont(:))
                  S(:)=(1.-Rnu(:))*NLTE%Sl(:,itran)+Rnu(:)*SBck(:)
               End if
               Iminus(0)=0. ! No incident radiation
! Compute monochromatic radiation field and local operator
               Call FormalSolution(NLTE, imu, inu, itran, NLTEInput%NLTE_formal_solution, &
               X, S, RNu, P, LStMuNu, UseLinear)
               If (Debug_level .ge. 1) then
                  If (checknan(Sum(P))) then
                     Do idepth=1, NLTE%NDEP
                        If (CheckNan(P(idepth))) then
                           Write (Message,*) 'P .lt. 0 in NLTE/radiation. ip, Sc, Scat, Jnu=', &
                                idepth,itran,inu,imu
                           Call Debug_Log(Message,2)
                           Write (Message,*) 'P .lt. 0 in NLTE/radiation. X, Sl, SBck, S, Rnu=',&
                                X(idepth),NLTE%SL(idepth,itran),SBck(idepth),S(idepth),RNu(idepth)
                           Call Debug_Log(Message,1)
                           Debug_errorflags(flag_NLTE)=1
                        End if
                     End do
                  End if
               End if
               !
               Jnu(:)=Jnu(:)+NLTE%WMu(imu)*P(:)
               If (Debug_level .ge. 1) then
                  If (checknan(Sum(Jnu))) then
                    Do idepth=1, NLTE%NDEP
                       If (CheckNan(Jnu(idepth))) then
                          Write (Message,*) 'Jnu .lt. 0 in NLTE/radiation. ip, itran,inu,imu,P=',idepth,itran,inu,imu,P(idepth)
                          Call Debug_Log(Message,1)
                          Debug_errorflags(flag_NLTE)=1
                       End if
                    End do                    
                 End if
              End if
            End do ! Do in imu
            Call WriteJ(irec,Jnu)
         End do ! Do in inu
     End do ! Do in itran
     Return
!
   End Subroutine Radiation
!
 End Subroutine SolveStat

!
! MULTI Collisional routines
!

      subroutine lcase(text)
!
!  converts a string to all lower case
!
      character (len=*) text
!
      l=len(text)
      do 100 i=1,l
        ic=ichar(text(i:i))
        if(ic.ge.65 .and. ic.le.90) text(i:i)=char(ic+32)
  100 continue
!
      return
    end subroutine lcase
!
      SUBROUTINE GETWRD(TEXT,K0,K1,K2)
!
!  FINDS NEXT WORD IN TEXT FROM INDEX K0. NEXT WORD IS TEXT(K1:K2)
!  THE NEXT WORD STARTS AT THE FIRST ALPHANUMERIC CHARACTER AT K0
!  OR AFTER. IT ENDS WITH THE LAST ALPHANUMERIC CHARACTER IN A ROW
!  FROM THE START
!
      INTEGER MSEPAR
      PARAMETER (MSEPAR=7)
      CHARACTER (LEN=*) TEXT
      CHARACTER SEPAR(MSEPAR)
      INTEGER K0,K1,K2,I,J
      DATA SEPAR/' ','(',')','=','*','/',','/
!
      K1=0
      DO 400 I=K0,LEN(TEXT)
        IF(K1.EQ.0) THEN
          DO 100 J=1,MSEPAR
            IF(TEXT(I:I).EQ.SEPAR(J)) GOTO 200
  100     CONTINUE
          K1=I
!
!  NOT START OF WORD
!
  200     CONTINUE
        ELSE
          DO 300 J=1,MSEPAR
            IF(TEXT(I:I).EQ.SEPAR(J)) GOTO 500
  300     CONTINUE
        ENDIF
  400 CONTINUE
!
!  NO NEW WORD. RETURN K1=K2=0
!
      K1=0
      K2=0
      GOTO 999
!
!  NEW WORD IN TEXT(K1:I-1)
!
  500 CONTINUE
      K2=I-1
! 
  999 CONTINUE
      RETURN
    End subroutine getwrd

      SUBROUTINE TAUTSP ( TAU, GTAU, NTAU, GAMMA, S, &
                          BREAK, COEF, L, K, IFLAG )
!  FROM  * A PRACTICAL GUIDE TO SPLINES *  BY C. DE BOOR
! CONSTRUCTS CUBIC SPLINE INTERPOLANT TO GIVEN DATA
!         TAU(I), GTAU(I), I=1,...,NTAU.
!  IF  GAMMA .GT. 0., ADDITIONAL KNOTS ARE INTRODUCED WHERE NEEDED TO
!  MAKE THE INTERPOLANT MORE FLEXIBLE LOCALLY. THIS AVOIDS EXTRANEOUS
!  INFLECTION POINTS TYPICAL OF CUBIC SPLINE INTERPOLATION AT KNOTS TO
!  RAPIDLY CHANGING DATA.
!
!  PARAMETERS
!            INPUT
!  TAU      SEQUENCE OF DATA POINTS. MUST BE STRICTLY INCREASING.
!  GTAU     CORRESPONDING SEQUENCE OF FUNCTION VALUES.
!  NTAU     NUMBER OF DATA POINTS. MUST BE AT LEAST  4 .
!  GAMMA    INDICATES WHETHER ADDITIONAL FLEXIBILITY IS DESIRED.
!          = 0., NO ADDITIONAL KNOTS
!          IN (0.,3.), UNDER CERTAIN CONDITIONS ON THE GIVEN DATA AT
!                POINTS I-1, I, I+1, AND I+2, A KNOT IS ADDED IN THE
!                I-TH INTERVAL, I=2,...,NTAU-2. SEE DESCRIPTION OF METH-
!                OD BELOW. THE INTERPOLANT GETS ROUNDED WITH INCREASING
!                GAMMA. A VALUE OF  2.5  FOR GAMMA IS TYPICAL.
!          IN (3.,6.), SAME , EXCEPT THAT KNOTS MIGHT ALSO BE ADDED IN
!                INTERVALS IN WHICH AN INFLECTION POINT WOULD BE PERMIT-
!                TED.  A VALUE OF  5.5  FOR GAMMA IS TYPICAL.
!            OUTPUT
!  BREAK, COEF, L, K  GIVE THE PP-REPRESENTATION OF THE INTERPOLANT.
!          SPECIFICALLY, FOR BREAK(I) .LE. X .LE. BREAK(I+1), THE
!        INTERPOLANT HAS THE FORM
!  F(X) = COEF(1,I) +DX(COEF(2,I) +(DX/2)(COEF(3,I) +(DX/3)COEF(4,I)))
!        WITH  DX = X - BREAK(I) AND I=1,...,L .
!  IFLAG   = 1, OK
!          = 2, INPUT WAS INCORRECT. A PRINTOUT SPECIFYING THE MISTAKE
!            WAS MADE.
!            WORKSPACE
!  S     IS REQUIRED, OF SIZE (NTAU,6). THE INDIVIDUAL COLUMNS OF THIS
!        ARRAY CONTAIN THE FOLLOWING QUANTITIES MENTIONED IN THE WRITE-
!        UP AND BELOW.
!     S(.,1) = DTAU = TAU(.+1) - TAU
!     S(.,2) = DIAG = DIAGONAL IN LINEAR SYSTEM
!     S(.,3) = U = UPPER DIAGONAL IN LINEAR SYSTEM
!     S(.,4) = R = RIGHT SIDE FOR LINEAR SYSTEM (INITIALLY)
!            = FSECND = SOLUTION OF LINEAR SYSTEM , NAMELY THE SECOND
!                       DERIVATIVES OF INTERPOLANT AT  TAU
!     S(.,5) = Z = INDICATOR OF ADDITIONAL KNOTS
!     S(.,6) = 1/HSECND(1,X) WITH X = Z OR = 1-Z. SEE BELOW.
!
!  ------  M E T H O D  ------
!  ON THE I-TH INTERVAL, (TAU(I), TAU(I+1)), THE INTERPOLANT IS OF THE
!  FORM
!  (*)  F(U(X)) = A + B*U + C*H(U,Z) + D*H(1-U,1-Z) ,
!  WITH  U = U(X) = (X - TAU(I))/DTAU(I). HERE,
!       Z = Z(I) = ADDG(I+1)/(ADDG(I) + ADDG(I+1))
!  (= .5, IN CASE THE DENOMINATOR VANISHES). WITH
!       ADDG(J) = ABS(DDG(J)), DDG(J) = DG(J+1) - DG(J),
!       DG(J) = DIVDIF(J) = (GTAU(J+1) - GTAU(J))/DTAU(J)
!  AND
!       H(U,Z) = ALPHA*U**3 + (1 - ALPHA)*(MAX(((U-ZETA)/(1-ZETA)),0)**3
!  WITH
!       ALPHA(Z) = (1-GAMMA/3)/ZETA
!       ZETA(Z) = 1 - GAMMA*MIN((1 - Z), 1/3)
!  THUS, FOR 1/3 .LE. Z .LE. 2/3,  F  IS JUST A CUBIC POLYNOMIAL ON
!  THE INTERVAL I. OTHERWISE, IT HAS ONE ADDITIONAL KNOT, AT
!         TAU(I) + ZETA*DTAU(I) .
!  AS  Z  APPROACHES  1, H(.,Z) HAS AN INCREASINGLY SHARP BEND  NEAR 1,
!  THUS ALLOWING  F  TO TURN RAPIDLY NEAR THE ADDITIONAL KNOT.
!     IN TERMS OF F(J) = GTAU(J) AND
!       FSECND(J) = 2.DERIVATIVE OF  F  AT  TAU(J),
!  THE COEFFICIENTS FOR (*) ARE GIVEN AS
!       A = F(I) - D
!       B = (F(I+1) - F(I)) - (C - D)
!       C = FSECND(I+1)*DTAU(I)**2/HSECND(1,Z)
!       D = FSECND(I)*DTAU(I)**2/HSECND(1,1-Z)
!  HENCE CAN BE COMPUTED ONCE FSECND(I),I=1,...,NTAU, IS FIXED.
!   F  IS AUTOMATICALLY CONTINUOUS AND HAS A CONTINUOUS SECOND DERIVAT-
!  IVE (EXCEPT WHEN Z = 0 OR 1 FOR SOME I). WE DETERMINE FSCND(.) FROM
!  THE REQUIREMENT THAT ALSO THE FIRST DERIVATIVE OF  F  BE CONTIN-
!  UOUS. IN ADDITION, WE REQUIRE THAT THE THIRD DERIVATIVE BE CONTINUOUS
!  ACROSS  TAU(2) AND ACROSS  TAU(NTAU-1) . THIS LEADS TO A STRICTLY
!  DIAGONALLY DOMINANT TRIDIAGONAL LINEAR SYSTEM FOR THE FSECND(I)
!  WHICH WE SOLVE BY GAUSS ELIMINATION WITHOUT PIVOTING.
!
!:
!: TAUTSP 89-08-30 MATS CARLSSON
!:        TEST OF ONEMZT CHANGED FROM 0. TO 1.E-11 TO AVOID DIVISION
!:        BY UNDERFLOWED EXPRESSION (ONEMZT**3)
!:
!:        89-09-05 MATS CARLSSON
!:        DIMENSION OF BREAK AND COEF CHANGED FROM 1 TO * TO MAKE 
!:        POSSIBLE CHECK OF INDEX OUT OF BOUNDS
!:
      INTEGER IFLAG,K,L,NTAU,   I,METHOD,NTAUM1
      DIMENSION BREAK(*),COEF(4,*),GTAU(NTAU),S(NTAU,6),TAU(NTAU)
      REAL, SAVE :: ONE
      DATA ONE/1.0/
!
!
!  THERE MUST BE AT LEAST  4  INTERPOLATION POINTS.
      IF (NTAU .GE. 4)                  GO TO 5
      PRINT *,'(8H NTAU =', NTAU,'  SHOULD BE .GE. 4'
      GO TO 999
!
!CONSTRUCT DELTA TAU AND FIRST AND SECOND (DIVIDED) DIFFERENCES OF DATA
!
    5 NTAUM1 = NTAU - 1
      DO I=1,NTAUM1
         S(I,1) = TAU(I+1) - TAU(I)
         IF (S(I,1) .GT. 0.)            GO TO 6
         PRINT *,'7H POINT ',I,TAU(I),' AND THE NEXT,2E15.6,15H ARE DISORDERED'
         GO TO 999
    6    S(I+1,4) = (GTAU(I+1)-GTAU(I))/S(I,1)
      END DO
      DO I=2,NTAUM1
    7    S(I,4) = S(I+1,4) - S(I,4)
      END DO
!
! CONSTRUCT SYSTEM OF EQUATIONS FOR SECOND DERIVATIVES AT  TAU. AT EACH
!  INTERIOR DATA POINT, THERE IS ONE CONTINUITY EQUATION, AT THE FIRST
!  AND THE LAST INTERIOR DATA POINT THERE IS AN ADDITIONAL ONE FOR A
!  TOTAL OF  NTAU  EQUATIONS IN  NTAU  UNKNOWNS.
!
      I = 2
      S(2,2) = S(1,1)/3.        
      SIXTH = 1./6.
      METHOD = 2
      GAM = GAMMA
      IF (GAM .LE. 0.)   METHOD = 1
      IF ( GAM .LE. 3.)                 GO TO 9
      METHOD = 3
      GAM = GAM - 3.
    9 ONEMG3 = 1. - GAM/3.
!                 ------ LOOP OVER I ------
   10 CONTINUE
!          CONSTRUCT Z(I) AND ZETA(I)
      Z = .5
      IF (METHOD .EQ. 1) GOTO 19
      IF (METHOD .EQ. 2) GOTO 11
      IF (METHOD .EQ. 3) GOTO 12

   11 IF (S(I,4)*S(I+1,4) .LT. 0.)      GO TO 19
   12 TEMP = ABS(S(I+1,4))
      DENOM = ABS(S(I,4)) + TEMP
      IF (DENOM .EQ. 0.)                GO TO 19
      Z = TEMP/DENOM
      IF (ABS(Z - .5) .LE. SIXTH)  Z = .5
   19 S(I,5) = Z
!   ******SET UP PART OF THE I-TH EQUATION WHICH DEPENDS ON
!         THE I-TH INTERVAL
      IF (Z - .5 .LT. 0) GOTO 21
      IF (Z - .5 .EQ. 0) GOTO 22
      IF (Z - .5 .GT. 0) GOTO 23

   21 ZETA = GAM*Z
      ONEMZT = 1. - ZETA
      ZT2 = ZETA**2
      ALPHA = ONE
      IF (ONEMG3/ONEMZT .LT. ONE) ALPHA=ONEMG3/ONEMZT
      FACTOR = ZETA/(ALPHA*(ZT2-1.) + 1.)
      S(I,6) = ZETA*FACTOR/6.
      S(I,2) = S(I,2) + S(I,1)*((1.-ALPHA*ONEMZT)*FACTOR/2. - S(I,6))
!     IF Z = 0 AND THE PREVIOUS Z = 1, THEN D(I) = 0. SINCE THEN
!     ALSO U(I-1) = L(I+1) = 0, ITS VALUE DOES NOT MATTER. RESET
!     D(I) = 1 TO INSURE NONZERO PIVOT IN ELIMINATION.
      IF (S(I,2) .LE. 0.) S(I,2) = 1.
      S(I,3) = S(I,1)/6.
                                        GO TO 25
   22 S(I,2) = S(I,2) + S(I,1)/3.
      S(I,3) = S(I,1)/6.
                                        GO TO 25
   23 ONEMZT = GAM*(1. - Z)
      ZETA = 1. - ONEMZT
      ALPHA = ONE
      IF (ONEMG3/ZETA .LT. ONE) ALPHA=ONEMG3/ZETA
      FACTOR = ONEMZT/(1. - ALPHA*ZETA*(1.+ONEMZT))
      S(I,6) = ONEMZT*FACTOR/6.
      S(I,2) = S(I,2) + S(I,1)/3.
      S(I,3) = S(I,6)*S(I,1)
   25 IF (I .GT. 2)                     GO TO 30
      S(1,5) = .5
!  ******THE FIRST TWO EQUATIONS ENFORCE CONTINUITY OF THE FIRST AND OF
!        THE THIRD DERIVATIVE ACROSS TAU(2).
      S(1,2) = S(1,1)/6.
      S(1,3) = S(2,2)
      ENTRY3 = S(2,3)

      IF (Z - .5 .LT. 0) GOTO 26
      IF (Z - .5 .EQ. 0) GOTO 27
      IF (Z - .5 .GT. 0) GOTO 28

   26 FACTR2 = ZETA*(ALPHA*(ZT2-1.) + 1.)/(ALPHA*(ZETA*ZT2-1.)+1.)
      RATIO = FACTR2*S(2,1)/S(1,2)
      S(2,2) = FACTR2*S(2,1) + S(1,1)
      S(2,3) = -FACTR2*S(1,1)
                                        GO TO 29
   27 RATIO = S(2,1)/S(1,2)
      S(2,2) = S(2,1) + S(1,1)
      S(2,3) = -S(1,1)
                                        GO TO 29
   28 RATIO = S(2,1)/S(1,2)
      S(2,2) = S(2,1) + S(1,1)
      S(2,3) = -S(1,1)*6.*ALPHA*S(2,6)
!       AT THIS POINT, THE FIRST TWO EQUATIONS READ
!              DIAG(1)*X1 + U(1)*X2 + ENTRY3*X3 = R(2)
!       -RATIO*DIAG(1)*X1 + DIAG(2)*X2 + U(2)*X3 = 0.
!       ELIMINATE FIRST UNKNOWN FROM SECOND EQUATION
   29 S(2,2) = RATIO*S(1,3) + S(2,2)
      S(2,3) = RATIO*ENTRY3 + S(2,3)
      S(1,4) = S(2,4)
      S(2,4) = RATIO*S(1,4)
                                        GO TO 35
   30 CONTINUE
!  ******THE I-TH EQUATION ENFORCES CONTINUITY OF THE FIRST DERIVATIVE
!        ACROSS TAU(I). IT HAS BEEN SET UP IN STATEMENTS 35 UP TO 40
!        AND 21 UP TO 25 AND READS NOW
!         -RATIO*DIAG(I-1)*XI-1 + DIAG(I)*XI + U(I)*XI+1 = R(I) .
!        ELIMINATE (I-1)ST UNKNOWN FROM THIS EQUATION
      S(I,2) = RATIO*S(I-1,3) + S(I,2)
      S(I,4) = RATIO*S(I-1,4) + S(I,4)
!
!  ******SET UP THE PART OF THE NEXT EQUATION WHICH DEPENDS ON THE
!        I-TH INTERVAL.
   35 IF (Z - .5 .LT. 0) GOTO 36
      IF (Z - .5 .EQ. 0) GOTO 37 
      IF (Z - .5 .GT. 0) GOTO 38

   36 RATIO = -S(I,6)*S(I,1)/S(I,2)
      S(I+1,2) = S(I,1)/3.
                                        GO TO 40
   37 RATIO = -(S(I,1)/6.)/S(I,2)
      S(I+1,2) = S(I,1)/3.
                                        GO TO 40
   38 RATIO = -(S(I,1)/6.)/S(I,2)
      S(I+1,2) = S(I,1)*((1.-ZETA*ALPHA)*FACTOR/2. - S(I,6))
!         ------  END OF I LOOP ------
   40 I = I+1
      IF (I .LT. NTAUM1)                GO TO 10
      S(I,5) = .5
!
!        ------  LAST TWO EQUATIONS  ------
!  THE LAST TWO EQUATIONS ENFORCE CONTINUITY OF THIRD DERIVATIVE AND
!  OF FIRST DERIVATIVE ACROSS  TAU(NTAU-1).
      ENTRY = RATIO*S(I-1,3) + S(I,2) + S(I,1)/3.
      S(I+1,2) = S(I,1)/6.
      S(I+1,4) = RATIO*S(I-1,4) + S(I,4)
      IF (Z - .5 .LT. 0) GOTO 41
      IF (Z - .5 .EQ. 0) GOTO 42 
      IF (Z - .5 .GT. 0) GOTO 43

   41 RATIO = S(I,1)*6.*S(I-1,6)*ALPHA/S(I-1,2)
      S(I,2) = RATIO*S(I-1,3) + S(I,1) + S(I-1,1)
      S(I,3) = -S(I-1,1)
                                        GO TO 45
   42 RATIO = S(I,1)/S(I-1,2)
      S(I,2) = RATIO*S(I-1,3) + S(I,1) + S(I-1,1)
      S(I,3) = -S(I-1,1)
                                        GO TO 45
   43 FACTR2 = ONEMZT*(ALPHA*(ONEMZT**2-1.)+1.)/ &
                     (ALPHA*(ONEMZT**3-1.)+1.)
      RATIO = FACTR2*S(I,1)/S(I-1,2)
      S(I,2) = RATIO*S(I-1,3) + FACTR2*S(I-1,1) + S(I,1)
      S(I,3) = -FACTR2*S(I-1,1)
!     AT THIS POINT, THE LAST TWO EQUATIONS READ
!             DIAG(I)*XI + U(I)*XI+1 = R(I)
!      -RATIO*DIAG(I)*XI + DIAG(I+1)*XI+1 = R(I+1)
!     ELIMINATE XI FROM LAST EQUATION
   45 S(I,4) = RATIO*S(I-1,4)
      RATIO = -ENTRY/S(I,2)
      S(I+1,2) = RATIO*S(I,3) + S(I+1,2)
      S(I+1,4) = RATIO*S(I,4) + S(I+1,4)
! 
!        ------ BACK SUBSTITUTION ------
!
      S(NTAU,4) = S(NTAU,4)/S(NTAU,2)
   50    S(I,4) = (S(I,4) - S(I,3)*S(I+1,4))/S(I,2)
         I = I - 1
         IF (I .GT. 1)                  GO TO 50
      S(1,4) = (S(1,4)-S(1,3)*S(2,4)-ENTRY3*S(3,4))/S(1,2)
!
!        ------ CONSTRUCT POLYNOMIAL PIECES ------
!
      BREAK(1) = TAU(1)
      L = 1
      DO I=1,NTAUM1
         COEF(1,L) = GTAU(I)
         COEF(3,L) = S(I,4)
         DIVDIF = (GTAU(I+1)-GTAU(I))/S(I,1)
         Z = S(I,5)

         IF (Z - .5 .LT. 0) GOTO 61
         IF (Z - .5 .EQ. 0) GOTO 62 
         IF (Z - .5 .GT. 0) GOTO 63

   61    IF (Z .EQ. 0.)                 GO TO 65
         ZETA = GAM*Z
         ONEMZT = 1. - ZETA
         C = S(I+1,4)/6.
         D = S(I,4)*S(I,6)
         L = L + 1
         DEL = ZETA*S(I,1)
         BREAK(L) = TAU(I) + DEL
         ZT2 = ZETA**2
         ALPHA = ONE
         IF (ONEMG3/ONEMZT .LT. ONE) ALPHA=ONEMG3/ONEMZT
         FACTOR = ONEMZT**2*ALPHA
         COEF(1,L) = GTAU(I) + DIVDIF*DEL &
                   + S(I,1)**2*(D*ONEMZT*(FACTOR-1.)+C*ZETA*(ZT2-1.))
         COEF(2,L) = DIVDIF + S(I,1)*(D*(1.-3.*FACTOR)+C*(3.*ZT2-1.))
         COEF(3,L) = 6.*(D*ALPHA*ONEMZT + C*ZETA)
         COEF(4,L) = 6.*(C - D*ALPHA)/S(I,1)
         COEF(4,L-1) = COEF(4,L) - 6.*D*(1.-ALPHA)/(DEL*ZT2)
         COEF(2,L-1) = COEF(2,L) - DEL*(COEF(3,L) -(DEL/2.)*COEF(4,L-1))
                                        GO TO 68
   62    COEF(2,L) = DIVDIF - S(I,1)*(2.*S(I,4) + S(I+1,4))/6.
         COEF(4,L) = (S(I+1,4)-S(I,4))/S(I,1)
                                        GO TO 68
   63    ONEMZT = GAM*(1. - Z)
         IF (ABS(ONEMZT) .LT. 1.E-11)            GO TO 65
         ZETA = 1. - ONEMZT
         ALPHA = ONE
         IF (ONEMG3/ZETA .LT. ONE) ALPHA=ONEMG3/ZETA
         C = S(I+1,4)*S(I,6)
         D = S(I,4)/6.
         DEL = ZETA*S(I,1)
         BREAK(L+1) = TAU(I) + DEL
         COEF(2,L) = DIVDIF - S(I,1)*(2.*D + C)
         COEF(4,L) = 6.*(C*ALPHA - D)/S(I,1)
         L = L + 1
         COEF(4,L) = COEF(4,L-1) + 6.*(1.-ALPHA)*C/(S(I,1)*ONEMZT**3)
         COEF(3,L) = COEF(3,L-1) + DEL*COEF(4,L-1)
         COEF(2,L) = COEF(2,L-1)+DEL*(COEF(3,L-1)+(DEL/2.)*COEF(4,L-1))
         COEF(1,L) = COEF(1,L-1)+DEL*(COEF(2,L-1)+(DEL/2.)*(COEF(3,L-1) &
                        +(DEL/3.)*COEF(4,L-1)))
                                        GO TO 68
   65    COEF(2,L) = DIVDIF
         COEF(3,L) = 0.
         COEF(4,L) = 0.
   68    L = L + 1
   70    BREAK(L) = TAU(I+1)
      END DO
!*      L = L - 1   DELETION OF THIS STATEMENT MAKES EXTRAPOL WARNING POSSIBLE
      K = 4
      IFLAG = 1
                                        RETURN
  999 IFLAG = 2
                                        RETURN
      END Subroutine tautsp
!*********************************************************************  
      SUBROUTINE INTERV ( XT, LXT, X, LEFT, MFLAG ) 
!  FROM  * A PRACTICAL GUIDE TO SPLINES *  BY C. DE BOOR                
!OMPUTES  LEFT = MAX( I , 1 .LE. I .LE. LXT  .AND.  XT(I) .LE. X )  .   
!                                                                       
!******  I N P U T  ******                                              
!  XT.....A REAL*8 SEQUENCE, OF LENGTH  LXT , ASSUMED TO BE NONDECREASIN
!  LXT.....NUMBER OF TERMS IN THE SEQUENCE  XT .                        
!  X.....THE POINT WHOSE LOCATION WITH RESPECT TO THE SEQUENCE  XT  IS  
!        TO BE DETERMINED.                                              
!                                                                       
!******  O U T P U T  ******                                            
!  LEFT, MFLAG.....BOTH INTEGERS, WHOSE VALUE IS                        
!                                                                       
!   1     -1      IF               X .LT.  XT(1)                        
!   I      0      IF   XT(I)  .LE. X .LT. XT(I+1)                       
!  LXT     1      IF  XT(LXT) .LE. X                                    
!                                                                       
!        IN PARTICULAR,  MFLAG = 0 IS THE 'USUAL' CASE.  MFLAG .NE. 0   
!        INDICATES THAT  X  LIES OUTSIDE THE HALFOPEN INTERVAL          
!        XT(1) .LE. Y .LT. XT(LXT) . THE ASYMMETRIC TREATMENT OF THE    
!        INTERVAL IS DUE TO THE DECISION TO MAKE ALL PP FUNCTIONS CONT- 
!        INUOUS FROM THE RIGHT.                                         
!                                                                       
!******  M E T H O D  ******                                            
!  THE PROGRAM IS DESIGNED TO BE EFFICIENT IN THE COMMON SITUATION THAT 
!  IT IS CALLED REPEATEDLY, WITH  X  TAKEN FROM AN INCREASING OR DECREA-
!  SING SEQUENCE. THIS WILL HAPPEN, E.G., WHEN A PP FUNCTION IS TO BE   
!  GRAPHED. THE FIRST GUESS FOR  LEFT  IS THEREFORE TAKEN TO BE THE VAL-
!  UE RETURNED AT THE PREVIOUS CALL AND STORED IN THE  L O C A L  VARIA-
!  BLE  ILO . A FIRST CHECK ASCERTAINS THAT  ILO .LT. LXT (THIS IS NEC- 
!  ESSARY SINCE THE PRESENT CALL MAY HAVE NOTHING TO DO WITH THE PREVI- 
!  OUS CALL). THEN, IF  XT(ILO) .LE. X .LT. XT(ILO+1), WE SET  LEFT =   
!  ILO  AND ARE DONE AFTER JUST THREE COMPARISONS.                      
!     OTHERWISE, WE REPEATEDLY DOUBLE THE DIFFERENCE  ISTEP = IHI - ILO 
!  WHILE ALSO MOVING  ILO  AND  IHI  IN THE DIRECTION OF  X , UNTIL     
!                      XT(ILO) .LE. X .LT. XT(IHI) ,                    
!  AFTER WHICH WE USE BISECTION TO GET, IN ADDITION, ILO+1 = IHI .      
!  LEFT = ILO  IS THEN RETURNED.                                        
!                                                                       
!:                                                                      
!: INTERV 88-08-30 MATS CARLSSON                                        
!:        STOP IF MORE THAN 20 EXTRAPOLATIONS                           
!:                                                                      
!:        88-10-26 MATS CARLSSON                                        
!:        SAVE STATEMENT INSERTED, EXTRAPOLATION CHECK CHANGED          
!:        TO WORK FOR THE CASE X=XT                                     
!:                                                                      
!:        94-01-06 MATS CARLSSON                                        
!:        IEXTR IN COMMON TO ENABLE ZEROING OUTSIDE ROUTINE             
!:        INITIALIZATION WITH DATA STATEMENT IN VIOLATION OF F77 STANDAR
!:                                                                      
!:        94-07-20 MATS CARLSSON                                        
!:        EXTRAPOLATION RESULTS IN ZERO VALUE AND WARNING MESSAGE       
!:        ONLY 10 WARNINGS PRINTED BUT PROGRAM CONTINUES                
!:                                                                      
!:        95-08-16 MATS CARLSSON                                        
!:        INITIALIZATION OF IEXTR MOVED TO BLOCK DATA SUBPROGRAM        
!:                                                                      
!:        95-11-27 MATS CARLSSON                                        
!:        MFLAG=1 ONLY SET IF X.GT.XT                                   
!:                                                                      
      COMMON /CINTV/ IEXTR 
      INTEGER LEFT,LXT,MFLAG,   IHI,ILO,ISTEP,MIDDLE 
      DIMENSION XT(LXT) 
      DATA ILO /1/ 
      SAVE ILO 
      IHI = ILO + 1 
      IF (IHI .LT. LXT)                 GO TO 20 
         IF (X .GE. XT(LXT))            GO TO 110 
         IF (LXT .LE. 1)                GO TO 90 
         ILO = LXT - 1 
         IHI = LXT 
!                                                                       
   20 IF (X .GE. XT(IHI))               GO TO 40 
      IF (X .GE. XT(ILO))               GO TO 100 
!                                                                       
!              **** NOW X .LT. XT(ILO) . DECREASE  ILO  TO CAPTURE  X . 
      ISTEP = 1 
   31    IHI = ILO 
         ILO = IHI - ISTEP 
         IF (ILO .LE. 1)                GO TO 35 
         IF (X .GE. XT(ILO))            GO TO 50 
         ISTEP = ISTEP*2 
                                        GO TO 31 
   35 ILO = 1 
      IF (X .LT. XT(1))                 GO TO 90 
                                        GO TO 50 
!              **** NOW X .GE. XT(IHI) . INCREASE  IHI  TO CAPTURE  X . 
   40 ISTEP = 1 
   41    ILO = IHI 
         IHI = ILO + ISTEP 
         IF (IHI .GE. LXT)              GO TO 45 
         IF (X .LT. XT(IHI))            GO TO 50 
         ISTEP = ISTEP*2 
                                        GO TO 41 
   45 IF (X .GE. XT(LXT))               GO TO 110 
      IHI = LXT 
!                                                                       
!           **** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL. 
   50 MIDDLE = (ILO + IHI)/2 
      IF (MIDDLE .EQ. ILO)              GO TO 100 
!     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1 .       
      IF (X .LT. XT(MIDDLE))            GO TO 53 
         ILO = MIDDLE 
                                        GO TO 50 
   53    IHI = MIDDLE 
                                        GO TO 50 
!**** SET OUTPUT AND RETURN.                                            
   90 MFLAG = -1 
      LEFT = 1 
      IF(IEXTR.LE.10) THEN 
        WRITE(*,92) X,XT(1) 
   92   FORMAT(' INTERV: X OUTSIDE XT, VALUE SET TO ZERO, X,XT=',       &
     &   1P,2E14.6)                                                     
      ENDIF 
      IEXTR=IEXTR+1 
                                        RETURN 
  100 MFLAG = 0 
      LEFT = ILO 
                                        RETURN 
  110 CONTINUE 
      LEFT = LXT-1 
      MFLAG = 0 
      IF(X.GT.XT(LXT)) THEN 
        MFLAG = 1 
        IF(IEXTR.LE.10) THEN 
          WRITE(*,92) X,XT(LXT) 
        ENDIF 
        IEXTR=IEXTR+1 
      ENDIF 
      RETURN 
     END SUBROUTINE INTERV
!*********************************************************************  
      FUNCTION PPVALU (BREAK, COEF, L, K, X, JDERIV ,MFLAG) 
!                                                                       
!: PPVALU 94-07-20 MATS CARLSSON                                        
!:        SETS VALUE TO 0 IF EXTRAPOLATION                              
!:                                                                      
!  FROM  * A PRACTICAL GUIDE TO SPLINES *  BY C. DE BOOR                
!ALLS  INTERV                                                           
!ALCULATES VALUE AT  X  OF  JDERIV-TH DERIVATIVE OF PP FCT FROM PP-REPR 
!                                                                       
!******  I N P U T  ******                                              
!  BREAK, COEF, L, K.....FORMS THE PP-REPRESENTATION OF THE FUNCTION  F 
!        TO BE EVALUATED. SPECIFICALLY, THE J-TH DERIVATIVE OF  F  IS   
!        GIVEN BY                                                       
!                                                                       
!     (D**J)F(X) = COEF(J+1,I) + H*(COEF(J+2,I) + H*( ... (COEF(K-1,I) +
!                             + H*COEF(K,I)/(K-J-1))/(K-J-2) ... )/2)/1 
!                                                                       
!        WITH  H = X - BREAK(I),  AND                                   
!                                                                       
!       I  =  MAX( 1 , MAX( J ,  BREAK(J) .LE. X , 1 .LE. J .LE. L ) ). 
!                                                                       
!  X.....THE POINT AT WHICH TO EVALUATE.                                
!  JDERIV.....INTEGER GIVING THE ORDER OF THE DERIVATIVE TO BE EVALUAT- 
!        ED.  A S S U M E D  TO BE ZERO OR POSITIVE.                    
!                                                                       
!******  O U T P U T  ******                                            
!  PPVALU.....THE VALUE OF THE (JDERIV)-TH DERIVATIVE OF  F  AT  X.     
!                                                                       
!******  M E T H O D  ******                                            
!     THE INTERVAL INDEX  I , APPROPRIATE FOR  X , IS FOUND THROUGH A   
!  CALL TO  INTERV . THE FORMULA ABOVE FOR THE  JDERIV-TH DERIVATIVE    
!  OF  F  IS THEN EVALUATED (BY NESTED MULTIPLICATION).                 
!                                                                       
!: PPVALU 95-11-06 MATS CARLSSON                                        
!:        RETURNS VALUE OF FLAG FROM INTERV                             
!:                                                                      
      INTEGER JDERIV,K,L,   I,M,MFLAG 
      DIMENSION BREAK(L),COEF(K,L) 
      PPVALU = 0. 
      FMMJDR = K - JDERIV 
!              DERIVATIVES OF ORDER  K  OR HIGHER ARE IDENTICALLY ZERO. 
      IF (FMMJDR .LE. 0.)               GO TO 99 
!                                                                       
!              FIND INDEX  I  OF LARGEST BREAKPOINT TO THE LEFT OF  X . 
      CALL INTERV ( BREAK, L, X, I, MFLAG ) 
!                                                                       
      IF(MFLAG.NE.0) THEN 
        PPVALU=0.0 
        RETURN 
      ENDIF 
!                                                                       
!      EVALUATE  JDERIV-TH DERIVATIVE OF  I-TH POLYNOMIAL PIECE AT  X . 
      H = X - BREAK(I) 
      M = K 
    9    PPVALU = (PPVALU/FMMJDR)*H + COEF(M,I) 
         M = M - 1 
         FMMJDR = FMMJDR - 1. 
         IF (FMMJDR .GT. 0.)            GO TO 9 
   99                                   RETURN 
       END FUNCTION PPVALU

!                                                                       
! **********************************************************************
!                                                                       
      SUBROUTINE CORONR(CEL,ION,TEMP,PRESS,RRATE,RIRATE) 
!                                                                       
!  CORONAL ION BALANCE CALCULATION                                      
!  NEW ROUTINE FOR USE WITH V2.0 OPACITY PACKAGE AND GENCOL             
!                                                                       
!  INPUT:  ELEMENT IDENTIFIER    CEL (CHARACTER)                        
!          ION STAGE IDENTIFIER  ION (INTEGER)                          
!          TEMP   (ELECTRON TEMPERATURE IN K       (RL)                 
!          PRESS  (ELECTRON PRESSURE IN DYNE/CM2)  (RL)                 
!                                                                       
!  OUTPUT: RRATE    RECOMBINATION RATE (CM3/S)                          
!          RIRATE   IONIZATION RATE    (CM3/S)                          
!                                                                       
!  EXAMPLE:                                                             
!          CEL='C  ', ION=1, WILL RETURN                                
!          THE RECOMBINATION RATE PER C II ATOM PER ELECTRON            
!          AND THE IONIZATION RATE PER C I ATOM PER ELECTRON            
!                                                                       
!  NOTES:                                                               
!          COMPUTATIONS FROM TABLES OF SHULL & VAN STEENBERG,           
!          AP J SUPPL  48 P95                                           
!          NO PE-DEPENDENCE, THIS COULD BE ADDED WHEN COMPUTATIONS BECOM
!          AVAILABLE                                                    
!          IF THE ION FRACTION IS NOT FOUND THEN THE ROUTINE WILL STOP  
!                                                                       
!  89-03-28 NEW ROUTINE (PHILIP JUDGE)                                  
!                                                                       
!: CORONR 92-08-10  MODIFICATIONS: (MATS CARLSSON)                      
!:        INTEGER ARRAYS WITH ELEMENT NAME CHANGED TO CHARACTER ARRAYS  
!:                                                                      
!:        92-09-09  MODIFICATIONS: (MATS CARLSSON)                      
!:        INTEGER IEL IN WRITE STATEMENT CHANGED TO CHARACTER           
!:                                                                      
!:        95-09-06  MODIFICATIONS: (MATS CARLSSON)                      
!:        COMMON BLOCK CLU CHANGED TO INCLUDE FILE                      
!:                                                                      
!                                                                       
      PARAMETER (MCOR=30) 
      INTEGER NCOR,IONCOR(MCOR) 
      DOUBLE PRECISION ARADC(MCOR),XRADC(MCOR),ADIC(MCOR),T0C(MCOR),    &
     & BDIC(MCOR),T1C(MCOR),ACOLC(MCOR),TCOLC(MCOR),TECOR,DECOR         
      COMMON /CCORON/ TECOR,DECOR,                                      &
     & ARADC,XRADC,ADIC,T0C,BDIC,T1C,ACOLC,TCOLC,NCOR,IONCOR            
      CHARACTER (LEN=3) :: CCOR 
      COMMON /CCCOR/ CCOR(MCOR) 
      CHARACTER (LEN=*) :: CEL 
!                                                                       
!  IDENTIFY THE ELEMENT AND IONIZATION STAGE TO BE COMPUTED             
!                                                                       
      DO 100 N=1,NCOR 
        IF( CEL .EQ. CCOR(N) .AND. ION .EQ. IONCOR(N) )GOTO 99 
  100 END DO 
      Print *,' CORONR: SPECIES',CEL,ION,' NOT FOUND IN ABSDAT FILE'
      Stop
!                                                                       
!  ELEMENT HAS BEEN COMPUTED- NOW RETURN THE RATES                      
!                                                                       
   99 RRATE=ARADC(N)*(TEMP/1.E4)**(-XRADC(N))  +                        &
     &  ADIC(N) * (TEMP**(-1.5)) * EXP(-T0C(N)/TEMP)*                   &
     &  (1.E0+BDIC(N) * (EXP(- T1C(N)/TEMP)))                           
      RIRATE=ACOLC(N) * SQRT(TEMP) * EXP( -TCOLC(N) / TEMP)             &
     &    / (1.E0 + 0.1 * TEMP / TCOLC(N))                              
      RETURN 
      END Subroutine Coronr           
!                                                                       
!*********************************************************************  
!                                                                       
      SUBROUTINE HEPOP(T,TOTH,HE1,HE2,HE3) 
!                                                                       
!  RATIOS OF HELIUM IONIZATION FRACTIONS FROM TABLE IV OF ARNAUD AND ROT
!  ABHE IS ABUNDANCE OF HELIUM - ASSUMED = 0.1                          
!                                                                       
      PARAMETER(NT=19) 
      DIMENSION TT(NT),ONE(NT),TWO(NT) 
      DATA TT /3.50, 4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60,          &
     &  4.70, 4.80, 4.90,                                               &
     &  5.00, 5.10, 5.20, 5.30, 5.40, 5.50, 5.60, 5.70/                 
      DATA ONE /0.0, 0.0 ,0.0 ,0.0 ,0.0 ,-0.07,-0.51,-1.33,-2.07,       &
     & -2.63,-3.20,-3.94,-4.67,-5.32,-5.90,-6.42,                       &
     & -6.90,-7.33,-7.73/                                               
      DATA TWO /-20.0,-9.05,-6.10,-3.75,-2.12,-0.84,-0.16,              &
     & -0.02,-0.01,-0.05,-0.34,-0.96,-1.60,-2.16,-2.63,                 &
     & -3.03,-3.38,-3.68,-3.94/                                         
      ABHE=0.1 
      TLOG=LOG10(T) 
      F1=10.**(SPLIN(TLOG,TT,ONE,NT,NT)) 
      F2=10.**(SPLIN(TLOG,TT,TWO,NT,NT)) 
      HE1=TOTH*ABHE*F1 
      HE2=TOTH*ABHE*F2 
      ALFA=1.-F1-F2 
      IF(ALFA .LT. 0.) THEN
         ALFA=1.E-30 
      End if
      ALFA=TOTH*ABHE*ALFA 
      RETURN 
      End Subroutine HePop
      

End Module NLTE_module
