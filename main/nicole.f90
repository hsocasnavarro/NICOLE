!                      N I C O L E   v 24.01
!       Non-LTE Inversion COde based on the Lorien Engine
!         By Hector Socas-Navarro, Jaime de la Cruz and
!                     Andres Asensio Ramos
!
!
! Please send your comments and bug reports to hsocas@iac.es
!
!
Program Nicole
  Use Param_structure
  Use Line_data_structure
  Use NICOLE_inp
  Use Compex_module
  Use Model_structure
  Use Lorien_module
  Use Forward_module
  Use Nodes_information
  Use Atomic_data
  Use File_operations
  Use debug_module
  Use Profiling
  Use Background_opacity_module, Only: Opacity_Package, Opacity_Package_UV, elneglectopac
  Use Eq_state, Only: eqstate_switch_others, eqstate_switch, eqstate_pe_consistency
  Implicit None
  Type (Parameters) :: Params
  Type (Line_data), dimension (:), allocatable :: Line
  Type (Region_data), dimension (:), allocatable :: Region
  Type (Model_2comp) :: Guess_model, Errors, Best_model, Start_Model
  Type (NICOLE_input) :: Input
  Type (Nodes_info) :: Nodes
  Integer, Parameter :: k18 = selected_int_kind(18)
  Real, Dimension (:), Allocatable :: Obs_profile, Syn_profile, Sigma, Sigma2, &
       TmpModel, TmpModel2, Tmp3, KeepVars
  Real, Dimension(100) :: TmpNodeloc
  Real :: Chisq, BestChisq
  Logical, Dimension (:), allocatable :: Brute_force
  Integer, Dimension(1) :: Locs
  Real, Dimension (:), allocatable :: MaskInvert, ChisqinFile
  Integer(Kind=8), Dimension (:), allocatable :: IntRecord
  Integer(Kind=4), Dimension (2) :: Integer4
  Integer :: status, nwlengths, iline, iregion, icycle
  Integer :: datainunit, modeloutunit, modeloutunit2, modeloutuniterr
  Integer :: profoutunit, i_inv, chisqunit
  Integer :: modelinunit, modelinunit2, strayinunit, headerunit, maskunit
  Integer :: tmpunit, depunit, depoutunit, percentunit, discard
  Integer :: ix, iy, nPix_x, nPix_y, ix1, iy1, ix2, iy2, irec0, icycle0
  Integer :: irec, iost1, Restart, ind
  Integer :: ivar, i0, i1, IProfUnit, logfileunit, perc, perc_old
  Character (Len=256) :: Filename
  Character (Len=10) :: Date, Time, Cyclesuffix
  Logical :: CheckNaN, FileExists, DepCoefExists, DepCoefOut, DidInvert, NodeLocExists
  Logical :: Completed
  Integer :: myrank, nprocs, SizeModel, isize, islave
  Integer ::  TaskComing, serialmode=0
  Integer, allocatable :: inverted_host(:), host_busy(:), pixel_to_do(:)
  Integer :: ierr, NotifyWhenFinished, request, ioldnl=-1, ioldnr=-1
  !
  Interface ! Needs to be defined so that DataIn can be left assumed-shape
     Subroutine Write_direct(CorrectEndian, iunit, irec, Datain, sizedata, iost)
       Logical, Intent(In) :: Correctendian
       integer, Intent(In) :: iunit, irec, sizedata, iost
       Real, Dimension(:), Intent(In) :: DataIn ! Assumed shape allows to 
                                                ! interpret real as size-1 array
     End Subroutine Write_direct
  End Interface
  !
  External :: CheckNaN
  myrank=0
  nprocs=1
  !
  If(myrank.eq.0) then
     Allocate(inverted_host(nprocs))
     inverted_host(1:nprocs)=0
     If (nprocs .eq. 1) serialmode=1
     Allocate(host_busy(0:nprocs-1+serialmode))
     host_busy(:)=0 
     host_busy(0)=1 ! Master is always busy
  end if
  Allocate (KeepVars(8))

  If (myrank .eq. 0) then
     Print *,''
     Print *,''
     Print *,'*************** N I C O L E   v 24.01 ******************'
     Print *,''
     Print *,'Lorien version: ',Lorien_ver
     Print *,'Forward version: ',Forward_ver
     Print *,'Compex version: ',Compex_ver
     Print *,''
     Print *,'********************************************************'
     Print *,''
     Print *,' This is the serial build'
     Print *,''
  EndIf
  Call Time_routine('main',.True.)
  !
  ! Determine number of cycles
  !
  Inquire(File='__input.dat_1',Exist=FileExists)
  If (FileExists) then 
     Call Open_file(headerunit,'__input.dat_1')
  Else
     Call Open_file(headerunit,'__input.dat')
  End if
  Read (headerunit, *) Params%ncycles
  If (Params%ncycles .gt. 1) then
     Print *,'Sorry. ncycles .gt. 1 is not supported any more'
     Print *,'Please break your problem down into single-cycle runs'
     Stop
     ! This is because there is a problem with changing the number of
     ! free parameters in the saved arrays (dydx_saved, dchidx_saved and
     ! d2chid2x_saved) in compute_dchisq_dx (lorien/lorien.f90)
  End if
  Read (headerunit, *) icycle0
  Call Close_File (headerunit)
  !
  ! Loop in cycles
  !
  Do icycle=icycle0, Params%ncycles
  ! {{{
     If (Params%Printout .ge. 1 .and. Params%ncycles .gt. 1) &
          Print *,'Inversion Cycle:',icycle

  !
  ! Read input file and set various parameters
  !
     ! {{{

  ! Initialize file units to -1
  datainunit=-1 ; modeloutunit=-1 ; modeloutuniterr=-1 ; profoutunit=-1
  iprofunit=-1 ; modeloutunit2=-1
  modelinunit=-1 ; strayinunit=-1 ; headerunit=-1 ; chisqunit=-1 ; maskunit=-1
  tmpunit=-1 ; depunit=-1; depoutunit=-1; modelinunit2=-1
  Write(Cyclesuffix,'("_",i0)') icycle
  ! 
  ! Header file is read by all processes
  Inquire(File='__input.dat'//Cyclesuffix,Exist=FileExists)
  If (FileExists) then 
     Call Open_file(headerunit,'__input.dat'//Cyclesuffix)
  Else
     Call Open_file(headerunit,'__input.dat') 
  End if
  Read (headerunit, *) discard
  Read (headerunit, *) discard
  Read (headerunit, *) nPix_x
  Read (headerunit, *) nPix_y
  Read (headerunit, *) Params%n_points
  Read (headerunit, *) ix1,iy1
  Read (headerunit, *) ix2,iy2
  Read (headerunit, *) irec0
  Read (headerunit, *) Input%mode
  Read (headerunit, '(A)') Input%model_in_file
  Read (headerunit, '(A)') Input%model_in_file_2
  Read (headerunit, '(A)') Input%obs_profile_file
  Read (headerunit, '(A)') Input%syn_profile_file
  Read (headerunit, '(A)') Input%model_out_file
  Read (headerunit, '(A)') Input%model_out_file_2
  Read (headerunit, *) Params%formal_solution, Params%formal_boundary_cond
  Read (headerunit, *) Params%heliocentric
  Read (headerunit, *) Params%printout
  Read (headerunit, *) Params%max_inv_iters
  Read (headerunit, '(A)') Input%Stray_profile_file
  Read (headerunit, *) Params%noise
  Read (headerunit, *) Input%maxinv
  Read (headerunit, *) Input%acceptchisq
  Read (headerunit, *) Params%speed
  Read (headerunit, *) Params%always_compute_deriv,  Params%cent_der
  Read (headerunit, *) Gravity
  Read (headerunit, *) Params%Regularization
  Read (headerunit, *) Params%update_opac
  Read (headerunit, *) Params%Negligible_opacity
  Read (headerunit, *) Params%reference_cont
  Read (headerunit, '(A)') Input%input_dens
  Read (headerunit, *) KeepVars(1:8)
  Read (headerunit, '(A)') Input%hscale
  Read (headerunit, *) Opacity_package, Opacity_Package_UV
  Read (headerunit, *) i0
  elneglectopac(:)=.False.
  Do ivar=1, i0
     Read (headerunit, *) i1
     elneglectopac(i1)=.True.
  End do
  Read (headerunit, *) eqstate_switch, eqstate_switch_others, eqstate_pe_consistency
  Read (headerunit, *) Input%set_hydro, Input%set_nH, Restart
  Read (headerunit, *) Input%depcoef_mode, Input%write_depcoef
  if (ioldnl .ge. 0) ioldnl=Params%n_lines
  Read (headerunit, *) Params%n_lines
  if (ioldnr .ge. 0) ioldnr=Params%n_regions
  Read (headerunit, *) Params%n_regions
  Call Set_params(Input, Params)
  Params%nPix=nPix_x*nPix_y
  Params%Input_dens=Input%Input_dens
  !
  ! Allocate memory for line and region arrays
  !
  If (ioldnl .ge. 0) then ! It is not the first cycle
     If (ioldnl .ne. Params%n_lines .or. ioldnr .ne. Params%n_regions) then !
        ! Reset dimensions for Line and Region arrays below
        Deallocate(Line)
        Deallocate(Region)
     End if
  End if
  ioldnl=Params%n_lines
  ioldnr=Params%n_regions
  If (.not. Allocated(Line)) then
     Allocate(Line(Params%n_lines), stat=status)
     Allocate(Region(Params%n_regions), stat=status)
  End if
  !
  ! Read line data
  !
  Do iregion=1, Params%n_regions
     Read (headerunit, *) Region(iregion)%First_wlength, &
          Region(iregion)%Wave_step, Region(iregion)%nwavelengths, &
          Region(iregion)%Macro_enh, Region(iregion)%layer, &
          Region(iregion)%Opacity_enh, Region(iregion)%Obs_additive, &
          Region(iregion)%Obs_multiplicative, Region(iregion)%Obs_gauss_sigma
  End do
  Do iline=1, Params%n_lines
     Read (headerunit, *) Line(iline)%Width, &
          Line(iline)%Elem, Line(iline)%Ion_stage, Line(iline)%Wlength, &
          Line(iline)%VDW_enh, Line(iline)%Energy_low, Line(iline)%loggf, &
          Line(iline)%Desig_low, Line(iline)%Mult_low, Line(iline)%J_low, &
          Line(iline)%Desig_up, Line(iline)%Mult_up, Line(iline)%J_up, &
          Line(iline)%collisions, &
          Line(iline)%Bark_sigma, Line(iline)%Bark_alpha, &
          Line(iline)%Gamma_rad, Line(iline)%Gamma_Strk_12, Line(iline)%Gamma_vdW_16, &
          Line(iline)%NLTEtransition, Line(iline)%NLTE_nl_ratio, &
          Line(iline)%NLTE_nu_ratio, Line(iline)%Hyperfine, &
          Line(iline)%HF_Alow, Line(iline)%HF_Blow, &
          Line(iline)%HF_Aup, Line(iline)%HF_Bup, Line(iline)%SpinI, &
          Line(iline)%Extra_vmic
     If (Input%depcoef_mode .eq. 2) then
        If (.not. Input%write_depcoef) then
           Line(iline)%DepCoefMode=2
        Else
           Line(iline)%DepCoefMode=3
        Endif
     End if
!
     Line(iline)%Atomic_number=1
     Do While (Atom_char(Line(iline)%Atomic_number) .ne. Line(iline)%Elem &
          .and. Line(iline)%Atomic_number .le. N_elements)
        Line(iline)%Atomic_number=Line(iline)%Atomic_number+1
     End Do
  End do
  ! Read abundances
  Read (headerunit, *) Starting_at_abund
  ! Nodes
  Read (headerunit, *) nodes%n_nodes_t
  Read (headerunit, *) nodes%n_nodes_v
  Read (headerunit, *) nodes%n_nodes_mic
  Read (headerunit, *) nodes%n_nodes_blong
  Read (headerunit, *) nodes%n_nodes_bx
  Read (headerunit, *) nodes%n_nodes_by
  Read (headerunit, *) nodes%n_nodes_chrom_x
  Read (headerunit, *) nodes%n_nodes_chrom_y
  Read (headerunit, *) nodes%n_nodes_stray
  Read (headerunit, *) nodes%n_nodes_mac
  Read (headerunit, *) nodes%n_nodes_ffactor
  Read (headerunit, *) nodes%n_nodes_ab
  ! Read abundances to invert
  If (nodes%n_nodes_ab .gt. 0) then
     If (Associated(nodes%i_nodes_ab)) Deallocate(nodes%i_nodes_ab)
     Allocate(nodes%i_nodes_ab(nodes%n_nodes_ab))
     Read (headerunit, *) nodes%i_nodes_ab(1:nodes%n_nodes_ab)
  End if
  Read (headerunit, *) nodes%n_nodes_t2
  Read (headerunit, *) nodes%n_nodes_v2
  Read (headerunit, *) nodes%n_nodes_mic2
  Read (headerunit, *) nodes%n_nodes_blong2
  Read (headerunit, *) nodes%n_nodes_bx2
  Read (headerunit, *) nodes%n_nodes_by2
  Read (headerunit, *) nodes%n_nodes_chrom_x2
  Read (headerunit, *) nodes%n_nodes_chrom_y2
  Read (headerunit, *) nodes%n_nodes_ab2
  ! Read abundances to invert
  If (nodes%n_nodes_ab2 .gt. 0) then
     If (Associated(nodes%i_nodes_ab2)) Deallocate(nodes%i_nodes_ab2)
     Allocate(nodes%i_nodes_ab2(nodes%n_nodes_ab2))
     Read (headerunit, *) nodes%i_nodes_ab2(1:nodes%n_nodes_ab2)
  End if
  ! Debug mode? (Note: Only Debug_level works in MPI mode)
  Read (headerunit, *) Debug_level, Params%Reinterpolate, Debug_outputpop, &
       Debug_outputcontop, Debug_outputNLTEsf
  ! NLTE Parameters
  Read (headerunit, *) Params%NLTE_Elim1
  Read (headerunit, *) Params%NLTE_ISUM, Params%NLTE_IStart
  Read (headerunit, *) Params%NLTE_NMU, Params%NLTE_UseColSwitch, Params%NLTE_Formal_Solution
  Read (headerunit, *) Params%NLTE_QNorm, Params%NLTE_CPER
  Read (headerunit, *) Params%NLTE_VelFree, Params%NLTE_NGAcc, Params%NLTE_NumLambdaIters
  Read (headerunit, *) Params%NLTE_OptThin, Params%NLTE_OptThick, Params%NLTE_Linear
  Read (headerunit, *) Params%NLTE_MaxIters
  Read (headerunit, *) Params%NLTE_ltepop

  ! }}}
  !
  ! Find out number of data and allocate memory for arrays
  !
  Params%n_data=0
  Do iregion=1, Params%n_regions
     nwlengths=Region(iregion)%nwavelengths
     Params%n_data=Params%n_data+nwlengths
  End do
  Params%n_data=Params%n_data*4 ! 4 Stokes parameters for each wlength
  If (Allocated(TmpModel)) Deallocate(TmpModel)
  If (Allocated(TmpModel2)) Deallocate(TmpModel2)
  Allocate(TmpModel(nvarsdepth*Params%n_points+nvarssingle))
  Allocate(TmpModel2(nvarsdepth*Params%n_points+nvarssingle))
  TmpModel2(:)=0.
  If (.not. Allocated(Syn_profile)) then
     Allocate(Syn_profile(Params%n_data), stat=status) 
     Allocate(Params%Stray_prof(Params%n_data), stat=status)
     Allocate(Params%IProf(Params%n_data/4))
     Allocate (ChisqinFile(Params%nPix))
     Allocate (MaskInvert(Params%nPix))
     Allocate(Sigma(Params%n_data), Stat=status)
     Allocate(Sigma2(Params%n_data), Stat=status)
     Do iline=1, Params%n_lines
        Allocate(Line(iline)%b_low(Params%n_points))
        Allocate(Line(iline)%b_up(Params%n_points))
     End do
     Allocate(Obs_profile(Params%n_data), Stat=status)
     Allocate(pixel_to_do(Params%nPix))
     Allocate (Tmp3(2*Params%n_lines*Params%n_points)) ! For depcoef
  End if
  Syn_profile(:)=0.
  Params%Stray_prof(:)=0.
  Params%IProf(:)=0.
  Pixel_to_do(:)=1
  Call Allocate_model_2comp(Guess_model, Params%n_points)
  If (n_model_elements .lt. N_elements) then
     Print *,'n_model_elements=',n_model_elements,', N_elements=',N_elements
     Stop
  End if
  Call Allocate_model_2comp(Best_model, Params%n_points)
  Call Allocate_model_2comp(Start_model, Params%n_points)
  !
  ! Some initializations
  !
  If (Starting_at_abund(1) .gt. 0) then ! Starting abunds
     ! Copy starting abunds to model if it doesn't have its own abundances
     ! (flagged by whether abundance(H)=12)
     If (Abs(Guess_model%Comp1%abundance(1)-12.0) .gt. 0.001) &
          Guess_model%Comp1%Abundance(1:N_elements)=Starting_At_Abund(1:N_elements) 
     If (Abs(Guess_model%Comp2%abundance(1)-12.0) .gt. 0.001) &
          Guess_model%Comp2%Abundance(1:N_elements)=Starting_At_Abund(1:N_elements) 
  Endif
  If (KeepVars(1) .gt. -.01) Guess_model%comp1%Keep_El_p=KeepVars(1)
  If (KeepVars(2) .gt. -.01)   Guess_model%comp1%Keep_Gas_p=KeepVars(2)
  If (KeepVars(3) .gt. -.01)   Guess_model%comp1%Keep_Rho=KeepVars(3)
  If (KeepVars(4) .gt. -.01)  Guess_model%comp1%Keep_nH=KeepVars(4)
  If (KeepVars(5) .gt. -.01)  Guess_model%comp1%Keep_nHminus=KeepVars(5)
  If (KeepVars(6) .gt. -.01)  Guess_model%comp1%Keep_nHplus=KeepVars(6)
  If (KeepVars(7) .gt. -.01)   Guess_model%comp1%Keep_nH2=KeepVars(7)
  If (KeepVars(8) .gt. -.01)  Guess_model%comp1%Keep_nH2plus=KeepVars(8)
  If (KeepVars(1) .gt. -.01) Guess_model%comp2%Keep_El_p=KeepVars(1)
  If (KeepVars(2) .gt. -.01)   Guess_model%comp2%Keep_Gas_p=KeepVars(2)
  If (KeepVars(3) .gt. -.01)   Guess_model%comp2%Keep_Rho=KeepVars(3)
  If (KeepVars(4) .gt. -.01)  Guess_model%comp2%Keep_nH=KeepVars(4)
  If (KeepVars(5) .gt. -.01)  Guess_model%comp2%Keep_nHminus=KeepVars(5)
  If (KeepVars(6) .gt. -.01)  Guess_model%comp2%Keep_nHplus=KeepVars(6)
  If (KeepVars(7) .gt. -.01)   Guess_model%comp2%Keep_nH2=KeepVars(7)
  If (KeepVars(8) .gt. -.01)  Guess_model%comp2%Keep_nH2plus=KeepVars(8)
  ChisqinFile(:)=1.e11
  MaskInvert(:)=1.
  Do iline=1, Params%n_lines
     Line(iline)%b_low(:)=1.
     Line(iline)%b_up(:)=1.
  End Do
  !
  ! Read files that both master and slaves will read
  !
  ! {{{

  ! Read maskinvert (if it exists)
  Inquire(File='maskinvert.dat', Exist=FileExists)
  If (FileExists) then
     Call Open_file_direct(maskunit,'maskinvert.dat', &
          RealBytes)
  Else
     Inquire(File='maskinvert.dat'//Cyclesuffix, Exist=FileExists)
     If (FileExists) then
        Call Open_file_direct(maskunit,'maskinvert.dat'//Cyclesuffix, &
             RealBytes)
     End if
  End if
  If (FileExists) then
     irec=0
     Do ix=1, Params%nPix
        irec=irec+1
        Call Read_direct(LittleEndian, &
             maskunit, irec, MaskInvert(ix), 1, iost1)
        If (iost1 .ne. 0) MaskInvert(ix)=1.
        If (MaskInvert(ix) .lt. 0.1) Pixel_to_do(ix)=0
     End do
     Call Close_File(maskunit)
     MaskInvert(Params%npix)=1 ! Always do last point no matter what 
  End if
  ! If we're restarting a previous run, read the Chisq from 
  !  previous run (if it exists)
  If (Restart .eq. 1) then
     Inquire(File='Chisq.dat', Exist=FileExists)
     If (FileExists) then
        Call Open_file_direct(chisqunit,'Chisq.dat', &
             RealBytes)
     Else
        Inquire(File='Chisq.dat'//Cyclesuffix, Exist=FileExists)
        If (FileExists) then
           Call Open_file_direct(chisqunit,'Chisq.dat'//Cyclesuffix, &
                RealBytes)
        End if
     End if
     If (FileExists) then
        irec=0
        Do ix=1, Params%nPix
           irec=irec+1
           Call Read_direct(LittleEndian, &
                chisqunit, irec, ChisqinFile(ix), 1, iost1)
           If (iost1 .ne. 0) ChisqinFile(ix)=1.E10
           If (ChisqinFile(ix) .lt. Input%acceptchisq .and. Restart .eq. 1) Pixel_to_do(irec)=0
        End do
        Call Close_File(chisqunit)
     End if
  End if
  ! node locations (if it exists)
  UserNodeLocations(:,:)=-1
  Inquire (File='nodelocations.dat'//CycleSuffix, Exist=NodeLocExists)
  If (NodeLocExists) then
     Call Open_file_direct(tmpunit,'nodelocations.dat'//CycleSuffix,RealBytes)
     irec=1
     Do i0=1,12
        Call Read_direct(LittleEndian, tmpunit, irec, TmpNodeLoc(1), &
             1, iost1)
        UserNodeLocations(i0,1)=TmpNodeLoc(1)
        irec=irec+1
        If (UserNodeLocations(i0,1) .gt. 0) then
           Do i1=2,100
              Call Read_direct(LittleEndian, tmpunit, irec, TmpNodeLoc(i1), &
                   1, iost1)
              UserNodeLocations(i0,i1)=TmpNodeLoc(i1)
              irec=irec+1
           End do
        End if
     End do
     Call Close_File(tmpunit)
  End if
  ! depcoef (if it exists)
  Inquire (File='depcoef.dat', Exist=DepCoefExists) ! File with departure coefficients 
  If (.not. DepCoefExists) then
     Do iline=1, Params%n_lines
        Line(iline)%DepCoefMode=0
     End do
  End if
  ! Instrumental profile exists?
  Inquire (File='Instrumental_profile.dat', Exist=Params%IProfExists)
  If (Params%IProfExists) then
     Call Open_file_direct(IProfUnit,'Instrumental_profile.dat', RealBytes*Params%n_data/4)
  Else
     Inquire (File='Instrumental_profile.dat'//CycleSuffix, Exist=Params%IProfExists)
     If (Params%IProfExists) &
     Call Open_file_direct(IProfUnit,'Instrumental_profile.dat'//CycleSuffix, RealBytes*Params%n_data/4)
  End if
  If (IProfUnit .gt. 0) then
     If (.not. Allocated(Params%mm)) &
          Allocate(Params%mm(Params%n_regions))
     Call Read_Direct(LittleEndian, IProfUnit, 1, Params%IProf, Params%n_data/4, iost1)
     Do iregion=1, Params%n_regions
        Params%mm(iregion)=Params%IProf(iregion)
     End do
     If (Params%Printout .ge. 2 .and. myrank .eq. 0) &
          Print *,'Reading Instrumental Profile'
     Call Close_File(IProfUnit)
  Endif
  Params%WeightsExist=.FALSE.
  Inquire (File='Weights.pro', Exist=Params%WeightsExist)
  If (Params%WeightsExist) Params%WeightsFilename='Weights.pro'
  If (.not. Params%WeightsExist) then 
     Inquire (File='Weights.pro'//CycleSuffix, Exist=Params%WeightsExist)
     If (Params%WeightsExist) Params%WeightsFilename='Weights.pro'//CycleSuffix
  End if
  If(Params%WeightsExist .and. Input%mode .eq. 'i') then
     call read_weights(Params,Sigma)
     If(Params%Printout .ge. 1 .and. myrank .eq. 0) &
          write(*,*) 'Reading Weights file:',Params%WeightsFilename
  End if
  If (Trim(Input%model_in_file_2) .ne. '') Params%TwoComp=.True.


  ! }}}
  ! 
  ! Master: Open input files
  !
  If (myrank .eq. 0) then
     ! {{{

     If (Restart .eq. -1) then
        Inquire (File=Trim(Input%Model_out_file),Exist=FileExists)
        If (FileExists) then
           Call Open_File (tmpunit,Trim(Input%Model_out_file))
           Close (tmpunit, STATUS='DELETE')
        End if
        Inquire (File=Trim(Input%Model_out_file)//'.err',Exist=FileExists)
        If (FileExists) then
           Call Open_File (tmpunit,Trim(Input%Model_out_file)//'.err')
           Close (tmpunit, STATUS='DELETE')
        End if
        Inquire (File=Trim(Input%Syn_profile_file),Exist=FileExists)
        If (FileExists) then
           Call Open_File (tmpunit,Trim(Input%Syn_profile_file))
           Close (tmpunit, STATUS='DELETE')
        End if
     End if

     ! If this is not a Restart run, check if output files already exist
     ! In that case issue a warning as unwanted side effects might occur
     If (Restart .eq. 0) then
        Inquire(File=Trim(Input%Model_out_file),Exist=FileExists)
        If (FileExists) then
           Print *,'WARNING!! Outputfile already exists:',Trim(Input%Model_out_file)
           Print *,'This could result in a mixture of old and new results in that file.'
           Print *,'Hopefully you know what you are doing. Proceeding anyway'
        End if
        Inquire(File=Trim(Input%Model_out_file)//'.err',Exist=FileExists)
        If (FileExists) then
           Print *,'WARNING!! Outputfile already exists:',Trim(Input%Model_out_file)//'.err'
           Print *,'This could result in a mixture of old and new results in that file.'
           Print *,'Hopefully you know what you are doing. Proceeding anyway'
        End if
        Inquire(File=Trim(Input%Syn_profile_file),Exist=FileExists)
        If (FileExists) then
           Print *,'WARNING!! Outputfile already exists:',Trim(Input%Syn_profile_file)
           Print *,'This could result in a mixture of old and new results in that file.'
           Print *,'Hopefully you know what you are doing. Proceeding anyway'
        End if
     End if
     ! Open model in unit (leave it open, we'll read it as we go)
     Call Open_file_direct(modelinunit, Trim(Input%model_in_file), &
          RealBytes*(nvarsdepth*Params%n_points+nvarssingle))
     If (Trim(Input%model_in_file_2) .ne. '') &
        Call Open_file_direct(modelinunit2, Trim(Input%model_in_file_2), &
             RealBytes*(nvarsdepth*Params%n_points+nvarssingle))
     ! Open profile out unit (leave it open, we'll write it as we go)
     Call Open_file_direct(profoutunit, Trim(Input%Syn_profile_file), &
       RealBytes*Params%n_data)
     Allocate (IntRecord(Params%n_data)) ! Write signature in first record
     IntRecord(:)=0
     IntRecord(1)=3328834590979877230_k18
     IntRecord(2)=2314885530823713331_k18
! Para leerlos como cadenas en Python:
! import struct
! int_record_1 = 4049129056382445934
! int_record_2 = 2314885530823504944
! bytes_record_1 = struct.pack('<q', int_record_1)
! bytes_record_2 = struct.pack('<q', int_record_2)
! string_representation = (bytes_record_1 + bytes_record_2).decode('utf-8')
! print(string_representation)          
     If (LittleEndian) then
        Integer4(1)=nPix_x
        Integer4(2)=nPix_y
     Else
        Integer4(1)=nPix_y
        Integer4(2)=nPix_x
     End if
     IntRecord(3)=Transfer(Integer4,IntRecord(3))
     IntRecord(4)=Params%n_data/4
     Call Write_direct_int(LittleEndian,profoutunit,1,IntRecord,Params%n_data,iost1)
     Deallocate(IntRecord)
     If (Input%write_depcoef) &
          Call Open_file_direct(depoutunit, 'depcoef_out.dat'//Cyclesuffix, &
          RealBytes*2*Params%n_lines*Params%n_points)
     ! depcoef (if it exists)
     Inquire (File='depcoef.dat', Exist=DepCoefExists) ! File with departure coefficients 
     If (DepCoefExists) then
        !    Open depcoef file (leave it open, we'll read it as we go)
        Call Open_file_direct(depunit, 'depcoef.dat', &
             RealBytes*2*Params%n_lines*Params%n_points)
        If (Params%Printout .ge. 2 .and. myrank .eq. 0) &
             Print *,'Reading depcoef.dat'
     Else
        Inquire (File='depcoef.dat'//Cyclesuffix, Exist=DepCoefExists)        
        If (DepCoefExists) then
           Call Open_file_direct(depunit, 'depcoef.dat'//Cyclesuffix, &
                RealBytes*2*Params%n_lines*Params%n_points)
        End if
     End if
     ! Open depcoef file unit (leave it open, we'll read it as we go)
     If (.not. DepCoefExists) then
        Do iline=1, Params%n_lines
           Line(iline)%DepCoefMode=0
        End do
     End if
     ! Open instrumental profile file (if it exists)
     Inquire (File='Instrumental_profile.dat', Exist=Params%IProfExists)
     Params%Iprof(:)=0.
     If (Params%IProfExists) then
        Call Open_file_direct(IProfUnit,'Instrumental_profile.dat', RealBytes*Params%n_data/4)
     Else
        Inquire (File='Instrumental_profile.dat'//Cyclesuffix, Exist=Params%IProfExists)
        If (Params%IProfExists) &
             Call Open_file_direct(IProfUnit,'Instrumental_profile.dat'//CycleSuffix, RealBytes*Params%n_data/4)
           
     End if
     ! Open stray light profile unit (leave it open, we'll read it as we go)
     Params%Stray_prof(1: Params%n_data)=0.
     If (Input%Stray_profile_file .ne. '') &
          Call Open_file_direct(strayinunit, Input%Stray_profile_file, &
          RealBytes*Params%n_data)
     If (Input%Model_out_file .ne. '') then
        Call Open_file_direct(modeloutunit, Trim(Input%Model_out_file), &
             RealBytes*(nvarsdepth*Params%n_points+nvarssingle))
        Allocate (IntRecord(nvarsdepth*Params%n_points+nvarssingle)) ! Write signature 
        IntRecord(:)=0
        IntRecord(1)=4049129056382445934_k18
        IntRecord(2)=2314885530823504944_k18
! Para leerlos como cadenas en Python:
! import struct
! int_record_1 = 4049129056382445934
! int_record_2 = 2314885530823504944
! bytes_record_1 = struct.pack('<q', int_record_1)
! bytes_record_2 = struct.pack('<q', int_record_2)
! string_representation = (bytes_record_1 + bytes_record_2).decode('utf-8')
! print(string_representation)
        If (LittleEndian) then
           Integer4(1)=nPix_x
           Integer4(2)=nPix_y
        Else
           Integer4(1)=nPix_y
           Integer4(2)=nPix_x
        End if
        IntRecord(3)=Transfer(Integer4,IntRecord(3))
        IntRecord(4)=Params%n_points
        Call Write_direct_int(LittleEndian,modeloutunit,1,IntRecord, &
             nvarsdepth*Params%n_points+nvarssingle,iost1)
        If (Params%TwoComp) then
           Call Open_file_direct(modeloutunit2, Trim(Input%Model_out_file_2), &
                RealBytes*(nvarsdepth*Params%n_points+nvarssingle))
           Call Write_direct_int(LittleEndian,modeloutunit2,1,IntRecord, &
                nvarsdepth*Params%n_points+nvarssingle,iost1)
        End if
        Deallocate(IntRecord)
     End if
     ! Open files needed only in inversion mode
     If (Input%mode .eq. 'i' .or. Input%mode .eq. 'c') then
        ! {{{

        ! Open other files needed in inversion mode and leave open
        Call Open_file_direct(chisqunit,'Chisq.dat'//Cyclesuffix, &
             RealBytes) ! To write
        Call Open_file_direct(datainunit, Input%Obs_profile_file, &
             RealBytes*Params%n_data)
        Call Open_file_direct(modeloutuniterr, Trim(Input%Model_out_file)//'.err', &
             RealBytes*(nvarsdepth*Params%n_points+nvarssingle))
        Allocate (IntRecord(nvarsdepth*Params%n_points+nvarssingle)) ! Write signature 
        IntRecord(:)=0
        IntRecord(1)=4049129056382445934_k18
        IntRecord(2)=2314885530823504944_k18
        If (LittleEndian) then
           Integer4(1)=nPix_x
           Integer4(2)=nPix_y
        Else
           Integer4(1)=nPix_y
           Integer4(2)=nPix_x
        End if
        IntRecord(3)=Transfer(Integer4,IntRecord(3))
        IntRecord(4)=Params%n_points
        Call Write_direct_int(LittleEndian,modeloutuniterr,1,IntRecord, &
             nvarsdepth*Params%n_points+nvarssingle,iost1)
        Deallocate(IntRecord)

        ! }}}
     End if ! end mode eq i

     ! }}}
  End If ! myrank .eq. 0
  !
  !
  ! Loop in points to synthesize/invert (pixels)
  ! Start operations: Master / Slave approach
  !
  If (myrank .eq. 0) then ! Master block. Loop to send jobs to slaves and wait for response
     ! {{{

     Call open_file (logfileunit,'logfile'//Cyclesuffix)
     Call open_file (percentunit,'percent_done'//Cyclesuffix)
     perc=0
     perc_old=-1
     ! Loop while there are pixels to do or slaves running
     Do While (MaxVal(Pixel_to_do) .eq. 1 .or. MaxVal(host_busy(1:)) .eq. 1)
        ! If there's an available slave, send a pixel. Otherwise, skip and listen
        If (MaxVal(Pixel_to_do) .eq. 1 .and. MinVal(host_busy) .eq. 0) then 
           ! Send a pixel to the first available slave
           ! {{{

           Locs=MaxLoc(Pixel_to_do) ! First pixel still pending to do
           irec=Locs(1)
           ix=(irec-1)/nPix_y +1
           iy=Modulo(irec-1,nPix_y) +1
           
           If (MinVal(host_busy) .eq. 1) then
              Print *,'No slave processes available. Something is wrong!'
              Stop
           End if
           Locs=MinLoc(host_busy) ! Identify first idle slave
           islave=Locs(1) - 1 ! -1 to translate indexing from 1 -> n to 0 -> n-1
           If (serialmode .eq. 1) islave=0
           host_busy(islave)=1 ! Mark slave as busy
           !
           ! Send data to slave
           !
           if(Input%mode.eq.'i') Call Read_direct(LittleEndian,datainunit,irec+1, &
                Obs_profile, Params%n_data, iost1)
           If (Input%Stray_profile_file .ne. '') then
              Call Read_direct(LittleEndian,strayinunit,irec+1, &
                   Params%Stray_prof, Params%n_data, iost1)
           End if
           If (IprofUnit .gt. 0) then
              Call Read_direct(LittleEndian,IprofUnit,irec+1, &
                   Params%IProf, Params%n_data/4, iost1)
              If (iost1 .gt. 0) then
                 Print *,'Error reading Instrumental_profile.dat',iost1
                 Stop
              End if
           End if
           Call Read_direct(LittleEndian,modelinunit,irec+1,TmpModel, & 
                nvarsdepth*Params%n_points+nvarssingle,iost1)           
              If (iost1 .gt. 0) then
                 Print *,'Error reading model',iost1
                 Stop
              End if
           If (Params%TwoComp) &
           Call Read_direct(LittleEndian,modelinunit2,irec+1,TmpModel2, & 
                nvarsdepth*Params%n_points+nvarssingle,iost1)           
              If (iost1 .gt. 0) then
                 Print *,'Error reading model 2',iost1
                 Stop
              End if
           SizeModel=Size(TmpModel)
           If (Starting_at_abund(1) .gt. 0) then
              If (Abs(TmpModel(SizeModel-N_elements+1)-12.0) .gt. 0.001) & 
                TmpModel(SizeModel-N_elements+1:SizeModel)=Starting_At_Abund(1:N_elements)
           End if
           If (Params%TwoComp) then
              If (Starting_at_abund(1) .gt. 0) then
                 If (Abs(TmpModel2(SizeModel-N_elements+1)-12.0) .gt. 0.001) & 
                      TmpModel2(SizeModel-N_elements+1:SizeModel)=Starting_At_Abund(1:N_elements)
              End if
           End if
           Call Record_to_model_2comp(Params%n_points, Guess_model, TmpModel, TmpModel2, KeepVars, 1)
           If (.not. Params%TwoComp) then
              Call model_assign(guess_model%comp2,guess_model%comp1)
              Guess_model%Comp1%ffactor=1.0
           End if
           TaskComing=1 ! Notify slave that a new task is coming

           If (DepCoefExists) then
              Call Read_direct(LittleEndian,depunit,irec,Tmp3, 2*Params%n_lines*Params%n_points,iost1)
              isize=Size(Tmp3)
           End if
           If (Input%write_depcoef) then
              
           End if
           If (nprocs .eq. 1) then ! If serial mode do the job here
              Call Do_Task
              Completed=.True. ! No need to wait in serial mode
           End if
           Pixel_to_do(irec)=0 ! Don't send this pixel to another slave
           Call Date_and_Time(Date, Time)
           If(nprocs.gt.1) then
              Write (logfileunit,'(A,I0,A,I0,A,I0,A,I0)') Time(1:2)//':'//Time(3:4) &
                   //':'//Time(5:6)//':Point:',irec,'(x=',ix,',y=',iy,') assigned to processor:',islave
           Else
              Write (logfileunit,*) Time(1:2)//':'//Time(3:4) &
                   //':'//Time(5:6)//':Started calculations for point:', irec,'(x=',ix,',y=',iy,')'
           End if

           ! }}}
        End if ! If (MaxVal(Pixels_to_do) .eq. 1 .and. MinVal(host_busy) .eq. 0)
        !
        ! Listen to see if anyone has finished
        !
        Do islave=1, nprocs-1+serialmode ! Check if any slave has finished its task
           If (Completed) then ! This slave has finished. Collect results and write to disk
              If (Input%Mode .eq. 'i' .and. DidInvert) then
                 Call Write_direct(LittleEndian, profoutunit, &
                      irec+1, Syn_Profile, Params%n_data, iost1)
                 Call Write_direct(LittleEndian, modeloutunit, &
                      irec+1, TmpModel, SizeModel, iost1)
                 If (Params%TwoComp) &
                      Call Write_direct(LittleEndian, modeloutunit2, &
                      irec+1, TmpModel2, SizeModel, iost1)
                 Call Write_direct(LittleEndian, chisqunit, &
                      irec, (/Chisq/), 1, iost1)
              End if
              If (Input%Mode .eq. 's' .and. DidInvert) then
                 Call Write_direct(LittleEndian, profoutunit, &
                      irec+1, Syn_Profile, Params%n_data, iost1)
                 If (modeloutunit .gt. 0) &
                      Call Write_direct(LittleEndian, modeloutunit, &
                      irec+1, TmpModel, SizeModel, iost1)
                 If (modeloutunit2 .gt. 0) &
                      Call Write_direct(LittleEndian, modeloutunit2, &
                      irec+1, TmpModel2, SizeModel, iost1)
              End if
              If (Input%Mode .eq. 'c' .and. Didinvert) then
                 Call Write_direct(LittleEndian, modeloutunit, &
                      irec+1, TmpModel, SizeModel, iost1)
                 If (Params%TwoComp) &
                      Call Write_direct(LittleEndian, modeloutunit2, &
                      irec+1, TmpModel2, SizeModel, iost1)
              End if
              If (Input%write_depcoef) then
                 Call Write_direct(LittleEndian, depoutunit, irec, Tmp3, 2*Params%n_lines*Params%n_points, iost1)
              End if
              ! Update percent file
              perc=nint(real(irec)*100/real(Params%nPix))
              If (perc .gt. perc_old) then
                 Call Date_and_Time(Date, Time)
                 write (percentunit,*) perc,'% done'
                 perc_old=perc
              End if
              Call MyFlush(percentunit)
              ! This slave is now free
              host_busy(islave)=0 
           End if ! If (Completed)
        End do ! Do islave=1, nprocs
        If (profoutunit .gt. 0) Call MyFlush(profoutunit)
        If (modeloutunit .gt. 0) Call MyFlush(modeloutunit)
        If (modeloutunit2 .gt. 0) Call MyFlush(modeloutunit2)
        If (chisqunit .gt. 0) Call MyFlush(chisqunit)
     End do ! While (MaxVal(Pixel_to_do) .eq. 1)
     !
     ! Notify slaves that all tasks are done so they don't need to wait for more tasks
     TaskComing=0
     !
     Write(logfileunit,*)' '
     Write(logfileunit,*)'---PROCESS-STATS---'
     Do islave=1,nprocs
        write(logfileunit,'(A,i0,A,I5,A,I7,A,F6.2,A)') 'host: ',&
             islave,', rank: ',islave-1,', pixels processed: ',&
             inverted_host(islave),' (',100*real(inverted_host(islave))/real(Params%nPix),'%)'
     end do
     !
     ! Close files in preparation for next cycle
     !
     Call Close_File(logfileunit)
     Call Close_File(percentunit)
     Call Close_File(modelinunit)
     Call Close_File(modelinunit2)
     Call Close_File(profoutunit)
     Call Close_File(depoutunit)
     If (DepCoefExists) Call Close_File(depunit)
     If (Input%Stray_profile_file .ne. '') Call Close_File(strayinunit)
     If (IProfUnit .gt. 0) then
        Call Close_File(IProfUnit)
        IProfUnit=-1
     End if
     Call Close_File(modeloutunit)
     Call Close_File(modeloutunit2)
     Call Close_File(chisqunit)
     Call Close_File(datainunit)
     Call Close_File(modeloutuniterr)

     ! }}}

  End if ! If (myrank .eq. 0)

  If (myrank .ge. 1) then ! Slave block. Receive job from master
     ! {{{

     Do while (.True.)
        ! Is there a task coming?
        If (TaskComing .eq. 0) Exit ! Master says no more tasks. Break loop
        ! If yes, collect data from master
        ! Got data. Do calculations now
        Call Do_Task ! Routine included with Contains, has access to all variables here
        !
        ! Done. Notify master by sending an asynchronous message
        ! Send results back to master
     End do ! Do while (.True)

     ! }}}
  End if ! If (myrank .ge. 1)
  !
  ! Next cycle. Free memory for possibly different new arrays
  !
  Deallocate(Syn_profile)
  Deallocate(Params%Stray_prof)
  Deallocate(Params%IProf)
  Deallocate (ChisqinFile)
  Deallocate (MaskInvert)
  Deallocate(Sigma)
  Deallocate(Sigma2)
  Do iline=1, Params%n_lines
     Deallocate(Line(iline)%b_low)
     Deallocate(Line(iline)%b_up)
  End do
  Deallocate(Obs_profile)
  Deallocate(pixel_to_do)
  Deallocate (Tmp3) ! For depcoef

  ! }}}
  End do ! Do in cycles
  !
  Call Time_routine('main',.False.)
  Call Time_log

  If (myrank .eq. 0) print *, "DONE"
  Stop

Contains
  Subroutine Do_Task
    Call Record_to_Model_2comp(Params%n_points, Guess_model, TmpModel, TmpModel2, KeepVars, 1) ! Extract model from array
    ix=(irec-1)/nPix_y +1
    iy=Modulo(irec-1,nPix_y) +1  
    BestChisq=ChisqinFile(irec)
    Best_model=Guess_model
    Start_model=Guess_model
    If (Pixel_to_do(irec) .eq. 1 .and. ix .ge. ix1 .and. ix .le. ix2 &
         .and. iy .ge. iy1 .and. iy .le. iy2 &
         .and. irec .ge. irec0 ) then 

       If (debug_level .ge. 1) then ! Open core file and write synthesis info
          ! {{{

          Write(Filename,'("core_",i0,".txt")') myrank
          Call Open_file(Debug_FileUnit, Filename)
          Write (Debug_FileUnit,*) ' Any mode'
          Write (Debug_FileUnit,*) ' nPix_x, nPix_y, nz, ndata'
          Write (Debug_FileUnit,*) nPix_x, nPix_y, Params%n_points, Params%n_data
          Write (Debug_FileUnit,*) ' irec, ix, iy'
          Write (Debug_FileUnit,*) irec, ix, iy ! irec starts with 1
          Write (Debug_FileUnit,*) ' Model. Component 1'
          Write (Debug_FileUnit,*) Guess_model%comp1%z_scale
          Write (Debug_FileUnit,*) Guess_model%comp1%ltau_500
          Write (Debug_FileUnit,*) Guess_model%comp1%temp
          Write (Debug_FileUnit,*) Guess_model%comp1%gas_p
          Write (Debug_FileUnit,*) Guess_model%comp1%rho
          Write (Debug_FileUnit,*) Guess_model%comp1%el_p
          Write (Debug_FileUnit,*) Guess_model%comp1%v_los
          Write (Debug_FileUnit,*) Guess_model%comp1%v_mic
          Write (Debug_FileUnit,*) Guess_model%comp1%b_long
          Write (Debug_FileUnit,*) Guess_model%comp1%b_x
          Write (Debug_FileUnit,*) Guess_model%comp1%b_y
          Write (Debug_FileUnit,*) Guess_model%comp1%bl_x
          Write (Debug_FileUnit,*) Guess_model%comp1%bl_y
          Write (Debug_FileUnit,*) Guess_model%comp1%bl_z
          Write (Debug_FileUnit,*) Guess_model%comp1%vl_x
          Write (Debug_FileUnit,*) Guess_model%comp1%vl_y
          Write (Debug_FileUnit,*) Guess_model%comp1%vl_z
          Write (Debug_FileUnit,*) Guess_model%comp1%ne
          Write (Debug_FileUnit,*) Guess_model%comp1%nH
          Write (Debug_FileUnit,*) Guess_model%comp1%nHplus
          Write (Debug_FileUnit,*) Guess_model%comp1%nHminus
          Write (Debug_FileUnit,*) Guess_model%comp1%nH2
          Write (Debug_FileUnit,*) Guess_model%comp1%v_mac
          Write (Debug_FileUnit,*) Guess_model%comp1%stray
          Write (Debug_FileUnit,*) Guess_model%comp1%ffactor
          Write (Debug_FileUnit,*) Guess_model%comp1%Abundance
          Write (Debug_FileUnit,*) ' Model. Component 2'
          Write (Debug_FileUnit,*) Guess_model%comp2%z_scale
          Write (Debug_FileUnit,*) Guess_model%comp2%ltau_500
          Write (Debug_FileUnit,*) Guess_model%comp2%temp
          Write (Debug_FileUnit,*) Guess_model%comp2%gas_p
          Write (Debug_FileUnit,*) Guess_model%comp2%rho
          Write (Debug_FileUnit,*) Guess_model%comp2%el_p
          Write (Debug_FileUnit,*) Guess_model%comp2%v_los
          Write (Debug_FileUnit,*) Guess_model%comp2%v_mic
          Write (Debug_FileUnit,*) Guess_model%comp2%b_long
          Write (Debug_FileUnit,*) Guess_model%comp2%b_x
          Write (Debug_FileUnit,*) Guess_model%comp2%b_y
          Write (Debug_FileUnit,*) Guess_model%comp2%bl_x
          Write (Debug_FileUnit,*) Guess_model%comp2%bl_y
          Write (Debug_FileUnit,*) Guess_model%comp2%bl_z
          Write (Debug_FileUnit,*) Guess_model%comp2%vl_x
          Write (Debug_FileUnit,*) Guess_model%comp2%vl_y
          Write (Debug_FileUnit,*) Guess_model%comp2%vl_z
          Write (Debug_FileUnit,*) Guess_model%comp2%ne
          Write (Debug_FileUnit,*) Guess_model%comp2%nH
          Write (Debug_FileUnit,*) Guess_model%comp2%nHplus
          Write (Debug_FileUnit,*) Guess_model%comp2%nHminus
          Write (Debug_FileUnit,*) Guess_model%comp2%nH2
          Write (Debug_FileUnit,*) Guess_model%comp2%v_mac
          Write (Debug_FileUnit,*) Guess_model%comp2%stray
          Write (Debug_FileUnit,*) Guess_model%comp2%ffactor
          Write (Debug_FileUnit,*) Guess_model%comp2%Abundance

          ! }}}
       End if
       Best_model=Guess_model
       Start_model=Guess_model
       Call Fill_densities(Params, Input%input_dens, Guess_model%Comp1)
       If (Params%TwoComp) &
            Call Fill_densities(Params, Input%input_dens, Guess_model%Comp2)
       If (Input%set_hydro) then
          Call Hydrostatic(Params, Guess_model%Comp1)
          Call Fill_densities(Params, Input%input_dens, Guess_model%Comp1)
          If (Params%TwoComp) then
             Call Hydrostatic(Params, Guess_model%Comp2)
             Call Fill_densities(Params, Input%input_dens, Guess_model%Comp2)
          End if
       End if
       If (Input%hscale .eq. 'z') then
          Call z_to_tau(Params, Guess_model%Comp1)
          If (Params%TwoComp) & 
               Call z_to_tau(Params, Guess_model%Comp2)
       Else If (Input%hscale .eq. 't') then
          Call tau_to_z(Params, Guess_model%Comp1)
          If (Params%TwoComp) & 
               Call tau_to_z(Params, Guess_model%Comp2)
       Else
          Print *,'Unknown hscale value'
          Stop
       End if

       If (debug_level .ge. 1) then
          Call MyFlush(Debug_FileUnit)
          Close (Debug_FileUnit,Status='DELETE')
          Debug_FileUnit=-1
       End if

       If (Input%mode .eq. 'c') then ! Convert mode
          Call Record_to_model_2comp(Params%n_points, Guess_model, TmpModel, TmpModel2,KeepVars,-1)
          DidInvert=.True.
          Return
       End if

    End if
    ! External departure coefficients?
    If (DepCoefExists) then ! Unpack from Tmp3
       Do iline=1, Params%n_lines
          i0=(iline-1)*2*Params%n_points +1
          i1=i0+Params%n_points-1
          Line(iline)%b_low(:)=Tmp3(i0:i1)
          i0=i1+1
          i1=i0+Params%n_points-1
          Line(iline)%b_up(:)=Tmp3(i0:i1)
       End do
    End if
    i_inv=1
    If (Input%mode .eq. 'i') then ! Call Lorien to do inversion
       ! {{{

       DidInvert=.False.
       If (Params%Reinterpolate .gt. 0) then
          Print *,'You should not use the Optimize Grid option in inversion mode'
          Print *,'Aborting...'
          Stop
       End if
       Do while ( i_inv .le. Input%maxinv .and. ( BestChisq .gt. Input%acceptchisq &
            .and. Pixel_to_do(irec) .eq. 1 &
            .and. ix .ge. ix1 .and. ix .le. ix2 &
            .and. iy .ge. iy1 .and. iy .le. iy2 &
            .and. irec .ge. irec0 )  )
          If (i_inv .gt. 1) then
             Call Randomize_Model(Params, Nodes, Start_Model%Comp1, Guess_Model%Comp1)
             Call Hydrostatic(Params, Guess_Model%Comp1)
             Call Fill_densities(Params, Input%input_dens, Guess_model%Comp1)
             If (Params%TwoComp) then
                Call Randomize_Model(Params, Nodes, Start_Model%Comp2, Guess_Model%Comp2)
                Call Hydrostatic(Params, Guess_Model%Comp2)
                Call Fill_densities(Params, Input%input_dens, Guess_model%Comp2)
             End if
          End if
          If (Params%Printout .ge. 1) &
               Print *,' Inversion try:',i_inv
          Chisq=1e10
          If (debug_level .ge. 1) then ! Open core file and write inversion info
             ! {{{

             Write(Filename,'("core_",i0,".txt")') myrank
             Call Open_file(Debug_FileUnit, Filename)
             Write (Debug_FileUnit,*) ' Inversion mode'
             Write (Debug_FileUnit,*) ' nPix_x, nPix_y, nz, ndata'
             Write (Debug_FileUnit,*) nPix_x, nPix_y, Params%n_points, Params%n_data
             Write (Debug_FileUnit,*) ' irec, ix, iy'
             Write (Debug_FileUnit,*) irec, ix, iy ! irec starts with 1
             Write (Debug_FileUnit,*) ' icycle'
             Write (Debug_FileUnit,*) icycle ! irec starts with 1
             Write (Debug_FileUnit,*) ' Observed profile'
             Write (Debug_FileUnit,*) Obs_profile
             Write (Debug_FileUnit,*) ' Guess model. Component 1'
             Write (Debug_FileUnit,*) Guess_model%comp1%z_scale
             Write (Debug_FileUnit,*) Guess_model%comp1%ltau_500
             Write (Debug_FileUnit,*) Guess_model%comp1%temp
             Write (Debug_FileUnit,*) Guess_model%comp1%gas_p
             Write (Debug_FileUnit,*) Guess_model%comp1%rho
             Write (Debug_FileUnit,*) Guess_model%comp1%el_p
             Write (Debug_FileUnit,*) Guess_model%comp1%v_los
             Write (Debug_FileUnit,*) Guess_model%comp1%v_mic
             Write (Debug_FileUnit,*) Guess_model%comp1%b_long
             Write (Debug_FileUnit,*) Guess_model%comp1%b_x
             Write (Debug_FileUnit,*) Guess_model%comp1%b_y
             Write (Debug_FileUnit,*) Guess_model%comp1%bl_x
             Write (Debug_FileUnit,*) Guess_model%comp1%bl_y
             Write (Debug_FileUnit,*) Guess_model%comp1%bl_z
             Write (Debug_FileUnit,*) Guess_model%comp1%vl_x
             Write (Debug_FileUnit,*) Guess_model%comp1%vl_y
             Write (Debug_FileUnit,*) Guess_model%comp1%vl_z
             Write (Debug_FileUnit,*) Guess_model%comp1%ne
             Write (Debug_FileUnit,*) Guess_model%comp1%nH
             Write (Debug_FileUnit,*) Guess_model%comp1%nHplus
             Write (Debug_FileUnit,*) Guess_model%comp1%nHminus
             Write (Debug_FileUnit,*) Guess_model%comp1%nH2
             Write (Debug_FileUnit,*) Guess_model%comp1%v_mac
             Write (Debug_FileUnit,*) Guess_model%comp1%stray
             Write (Debug_FileUnit,*) Guess_model%comp1%ffactor
             Write (Debug_FileUnit,*) Guess_model%comp1%Abundance
             Write (Debug_FileUnit,*) ' Guess model. Component 2'
             Write (Debug_FileUnit,*) Guess_model%comp2%z_scale
             Write (Debug_FileUnit,*) Guess_model%comp2%ltau_500
             Write (Debug_FileUnit,*) Guess_model%comp2%temp
             Write (Debug_FileUnit,*) Guess_model%comp2%gas_p
             Write (Debug_FileUnit,*) Guess_model%comp2%rho
             Write (Debug_FileUnit,*) Guess_model%comp2%el_p
             Write (Debug_FileUnit,*) Guess_model%comp2%v_los
             Write (Debug_FileUnit,*) Guess_model%comp2%v_mic
             Write (Debug_FileUnit,*) Guess_model%comp2%b_long
             Write (Debug_FileUnit,*) Guess_model%comp2%b_x
             Write (Debug_FileUnit,*) Guess_model%comp2%b_y
             Write (Debug_FileUnit,*) Guess_model%comp2%bl_x
             Write (Debug_FileUnit,*) Guess_model%comp2%bl_y
             Write (Debug_FileUnit,*) Guess_model%comp2%bl_z
             Write (Debug_FileUnit,*) Guess_model%comp2%vl_x
             Write (Debug_FileUnit,*) Guess_model%comp2%vl_y
             Write (Debug_FileUnit,*) Guess_model%comp2%vl_z
             Write (Debug_FileUnit,*) Guess_model%comp2%ne
             Write (Debug_FileUnit,*) Guess_model%comp2%nH
             Write (Debug_FileUnit,*) Guess_model%comp2%nHplus
             Write (Debug_FileUnit,*) Guess_model%comp2%nHminus
             Write (Debug_FileUnit,*) Guess_model%comp2%nH2
             Write (Debug_FileUnit,*) Guess_model%comp2%v_mac
             Write (Debug_FileUnit,*) Guess_model%comp2%stray
             Write (Debug_FileUnit,*) Guess_model%comp2%ffactor
             Write (Debug_FileUnit,*) Guess_model%comp2%Abundance

             ! }}}
          End if
          Call Select_number_of_nodes(Params, Guess_model, Nodes, icycle)          
          Allocate(Brute_force(Params%n_free_parameters), Stat=status)
          Nodes%Reference_model=Guess_model ! Model_assign operation
          If (.not.Params%WeightsExist) then
             Call Compute_weights(Params, Obs_profile, Sigma, icycle)
          end if
          !                       
          If (status .gt. 0) then
             Print *,'Out of memory while allocating Brute_force'
             Stop
          End if
          Brute_force(1:Params%n_free_parameters)=.TRUE.
          Sigma2(:)=Sigma
          Call Lorien(Params, Nodes, Line, Region, Guess_model, Obs_profile, &
               Sigma, Brute_force, Errors, Chisq)
          Deallocate (Brute_force)                             
          If (CheckNAN(Chisq)) Chisq=1.E10 
          If (Chisq .le. BestChisq) then
             DidInvert=.True.
             Best_Model=Guess_Model
             BestChisq=Chisq
          End if
          If (Params%Printout .ge. 1) &
               Print *,' Chisq=',Chisq,'. Best so far=',BestChisq
          i_inv=i_inv+1
          Chisq=BestChisq
          Guess_Model=Best_Model
          Call tau_to_z(Params, Guess_model%Comp1)
          If (Params%TwoComp) & 
               Call tau_to_z(Params, Guess_model%Comp2)          
          Call Forward(Params, Line, Region, Guess_model, Syn_profile, &
               .TRUE.)
          Call Record_to_model_2comp(Params%n_points, Guess_model, TmpModel, TmpModel2,KeepVars, -1)
          If (debug_level .ge. 1) then
             Call MyFlush(Debug_FileUnit)
             Close (Debug_FileUnit,Status='DELETE')
             Debug_FileUnit=-1
          End if
       End Do ! End while

       ! }}}
    End if ! input%mode .eq. 'i'
!
    If (Input%mode .eq. 's') then
       Syn_profile(:)=0.
       DidInvert=.False.
       If ( Pixel_to_do(irec) .eq. 1 .and. ix .ge. ix1 .and. ix .le. ix2 &
            .and. iy .ge. iy1 .and. iy .le. iy2 &
            .and. irec .ge. irec0 ) then
          If (debug_level .ge. 1) then ! Open core file and write synthesis info
             ! {{{

             Write(Filename,'("core_",i0,".txt")') myrank
             Call Open_file(Debug_FileUnit, Filename)
             Write (Debug_FileUnit,*) ' Synthesis mode'
             Write (Debug_FileUnit,*) ' nPix_x, nPix_y, nz, ndata'
             Write (Debug_FileUnit,*) nPix_x, nPix_y, Params%n_points, Params%n_data
             Write (Debug_FileUnit,*) ' irec, ix, iy'
             Write (Debug_FileUnit,*) irec, ix, iy ! irec starts with 1
             Write (Debug_FileUnit,*) ' Model. Component 1'
             Write (Debug_FileUnit,*) Guess_model%comp1%z_scale
             Write (Debug_FileUnit,*) Guess_model%comp1%ltau_500
             Write (Debug_FileUnit,*) Guess_model%comp1%temp
             Write (Debug_FileUnit,*) Guess_model%comp1%gas_p
             Write (Debug_FileUnit,*) Guess_model%comp1%rho
             Write (Debug_FileUnit,*) Guess_model%comp1%el_p
             Write (Debug_FileUnit,*) Guess_model%comp1%v_los
             Write (Debug_FileUnit,*) Guess_model%comp1%v_mic
             Write (Debug_FileUnit,*) Guess_model%comp1%b_long
             Write (Debug_FileUnit,*) Guess_model%comp1%b_x
             Write (Debug_FileUnit,*) Guess_model%comp1%b_y
             Write (Debug_FileUnit,*) Guess_model%comp1%bl_x
             Write (Debug_FileUnit,*) Guess_model%comp1%bl_y
             Write (Debug_FileUnit,*) Guess_model%comp1%bl_z
             Write (Debug_FileUnit,*) Guess_model%comp1%vl_x
             Write (Debug_FileUnit,*) Guess_model%comp1%vl_y
             Write (Debug_FileUnit,*) Guess_model%comp1%vl_z
             Write (Debug_FileUnit,*) Guess_model%comp1%ne
             Write (Debug_FileUnit,*) Guess_model%comp1%nH
             Write (Debug_FileUnit,*) Guess_model%comp1%nHplus
             Write (Debug_FileUnit,*) Guess_model%comp1%nHminus
             Write (Debug_FileUnit,*) Guess_model%comp1%nH2
             Write (Debug_FileUnit,*) Guess_model%comp1%v_mac
             Write (Debug_FileUnit,*) Guess_model%comp1%stray
             Write (Debug_FileUnit,*) Guess_model%comp1%ffactor
             Write (Debug_FileUnit,*) Guess_model%comp1%Abundance
             Write (Debug_FileUnit,*) ' Model. Component 2'
             Write (Debug_FileUnit,*) Guess_model%comp2%z_scale
             Write (Debug_FileUnit,*) Guess_model%comp2%ltau_500
             Write (Debug_FileUnit,*) Guess_model%comp2%temp
             Write (Debug_FileUnit,*) Guess_model%comp2%gas_p
             Write (Debug_FileUnit,*) Guess_model%comp2%rho
             Write (Debug_FileUnit,*) Guess_model%comp2%el_p
             Write (Debug_FileUnit,*) Guess_model%comp2%v_los
             Write (Debug_FileUnit,*) Guess_model%comp2%v_mic
             Write (Debug_FileUnit,*) Guess_model%comp2%b_long
             Write (Debug_FileUnit,*) Guess_model%comp2%b_x
             Write (Debug_FileUnit,*) Guess_model%comp2%b_y
             Write (Debug_FileUnit,*) Guess_model%comp2%bl_x
             Write (Debug_FileUnit,*) Guess_model%comp2%bl_y
             Write (Debug_FileUnit,*) Guess_model%comp2%bl_z
             Write (Debug_FileUnit,*) Guess_model%comp2%vl_x
             Write (Debug_FileUnit,*) Guess_model%comp2%vl_y
             Write (Debug_FileUnit,*) Guess_model%comp2%vl_z
             Write (Debug_FileUnit,*) Guess_model%comp2%ne
             Write (Debug_FileUnit,*) Guess_model%comp2%nH
             Write (Debug_FileUnit,*) Guess_model%comp2%nHplus
             Write (Debug_FileUnit,*) Guess_model%comp2%nHminus
             Write (Debug_FileUnit,*) Guess_model%comp2%nH2
             Write (Debug_FileUnit,*) Guess_model%comp2%v_mac
             Write (Debug_FileUnit,*) Guess_model%comp2%stray
             Write (Debug_FileUnit,*) Guess_model%comp2%ffactor
             Write (Debug_FileUnit,*) Guess_model%comp2%Abundance

             ! }}}
          End if
          Call Forward(Params, Line, Region, Guess_model, Syn_profile, Input%set_hydro)
	  Call Fill_densities(Params, Input%input_dens, Guess_model%Comp1)
          If (Params%TwoComp) &
               Call Fill_densities(Params, Input%input_dens, Guess_model%Comp2)
          Call Record_to_Model_2comp(Params%n_points, Guess_model, TmpModel, TmpModel2, KeepVars,-1) 
          If (debug_level .ge. 1) then
             Write (Debug_FileUnit,*) ' Profile:'
             Write (Debug_FileUnit,*) Syn_profile
             Call MyFlush(Debug_FileUnit)
             Close (Debug_FileUnit,Status='DELETE')
             Debug_FileUnit=-1
          End if
          DidInvert=.True.
       End if ! If ix ge ix1 and ix le ix2 & ...
    End if ! mode .eq. 's'
    Do iline=1, Params%n_lines
       i0=(iline-1)*2*Params%n_points +1
       i1=i0+Params%n_points-1
       Tmp3(i0:i1)=Line(iline)%b_low(:)
       i0=i1+1
       i1=i0+Params%n_points-1
       Tmp3(i0:i1)=Line(iline)%b_up(:)
    End do

    If (Params%Printout .ge. 1) &
         Print *,'Point ',irec,' of ',Params%nPix,' done by process ',myrank

    Return

  End Subroutine Do_Task


End Program Nicole

Subroutine getmyrank(myrank)

  Integer :: myrank, status


  myrank=0
  Return

End Subroutine getmyrank
