! SUBROUTINE zeeman_hyperfine
! 
! PURPOSE: computes the hyperfine splitting and strength of a line for a 
!          given B_field_strength. 
!          Adapted from zeeman_profile_hyperfine (Asensio-Ramos).
!
! OUTPUT:  The results are stored in the variables nPi, nLeft, nRight, SplitPi, 
!          SplitLeft, SplitRight, Strpi, StrLeft, StrRight. 
!          The splitting is given in cm^-1 units. 
!
! KNOWN ISSUES: the output variables are dimensioned using the parameter 
!               "transition" (param_structure). nPi, nLeft, nRight can be larger for some 
!               lines. This could be solved by first checking the number of components
!               and then allocating the arrays, but in that case the subroutine ZEEMAN 
!               needs to be adapted too.
!               
! Adapted by Jaime de la Cruz Rodriguez (IFA-UU 2011).
!
! ---
! MODIFICATION history:
!
! ---
!
Subroutine Zeeman_hyperfine(linea, Bfield4,  nPi, nLeft, nRight, &
     SplitPi, SplitLeft, SplitRight, Strpi, StrLeft, StrRight)
  Use Line_data_structure
  Use hyperfine
  Use Param_structure
  Use Profiling
  Implicit none
  Type (Line_data) :: Linea
  Real :: Bfield4
  integer :: nfl_up, f2m_up, nfl_low, f2m_low, nPi, nRight, nLeft, loop
  Logical, save :: firsttime = .TRUE.
  integer, allocatable :: nflev_up(:), nflev_low(:)
  real(kind=8), allocatable :: energy_up(:,:), c_up(:,:,:), energy_low(:,:), c_low(:,:,:)
  real, dimension(transitions) :: SplitPi, SplitRight, SplitLeft, StrPi, StrRight, StrLeft
  real(kind=8) :: s_low, s_up, l_low, l_up
  real(kind=8) :: Jlow, Jup, SpinI, Aup, Alow, Bup, Blow, Bfield
  !
  ! Init factorials and convert multiplicity and L to doubles.
  !
  Call Time_routine('zeeman_hyperfine',.True.)
  If(firsttime) then
     firsttime = .FALSE.
     call factrl
  endif
  Bfield=Bfield4
  Jup=Linea%J_up
  Jlow=Linea%J_low
  SpinI=Linea%SpinI
  Aup=Linea%HF_Aup
  Alow=Linea%HF_Alow
  Bup=Linea%HF_Bup
  Blow=Linea%HF_Blow
  !
  ! Get L and S values 
  !
  S_low = (linea%mult_low - 1.0d0) * 0.5d0
  S_up  = (linea%mult_up - 1.0d0) * 0.5d0
  L_low = Lchar2double(linea%desig_low)
  L_up = Lchar2double(linea%desig_up)
  !
  ! Calculate the size of the energy and eigenvectors arrays
  !
  Call Time_routine('hyperfine_size',.True.)
  call hyper_size(Jup, SpinI, nfl_up, f2m_up)
  call hyper_size(Jlow, SpinI, nfl_low, f2m_low)
  Call Time_routine('hyperfine_size',.False.)
  !
  ! Allocate and init vars
  !
  allocate(nflev_up(0:2*f2m_up))
  allocate(energy_up(0:2*f2m_up,0:nfl_up-1))
  allocate(c_up(0:2*f2m_up,0:nfl_up-1,0:f2m_up))
  nflev_up(:) = 0
  energy_up(:,:) = 0.d0
  c_up(:,:,:) = 0.d0
  !
  allocate(nflev_low(0:2*f2m_low))
  allocate(energy_low(0:2*f2m_low,0:nfl_low-1))
  allocate(c_low(0:2*f2m_low,0:nfl_low-1,0:f2m_low))
  nflev_low(:) = 0
  energy_low(:,:) = 0.d0
  c_low(:,:,:) = 0.d0
  !
  ! Diagonalize the hyperfine+magnetic Hamiltonian of the upper and lower levels
  !
  Call Time_routine('hyperfine_diagonalize',.True.)
  call hyper_diagonalize(L_up, S_up, Jup, SpinI, Aup, &
       Bup, Bfield, nflev_up, energy_up, c_up, f2m_up, nfl_up)
  call hyper_diagonalize(L_low, S_low, Jlow, SpinI, Alow, &
       Blow, Bfield, nflev_low, energy_low, c_low, f2m_low, nfl_low)
  Call Time_routine('hyperfine_diagonalize',.False.)
  !
  ! Calculate the number of Zeeman components for the transition
  !
  call hyper_components(Jup, Jlow, SpinI, f2m_up, &
       f2m_low, nflev_up, nflev_low, nPi, nLeft, nRight)
  !
  ! Check that npi, nRight, nLeft < transitions (defined in param_structure)
  !
  if((nPi.GT.transitions).OR.(nRight.GT.transitions).OR.(nLeft.GT.transitions)) then
     print *, 'zeeman_hyperfine : Error -> too many transitions'
     print *, 'zeeman_hyperfine : Recompile with greater transition parameter!'
     print *, '   Max number of components allowed:', transitions
     print *, '   Number of SigmaRight components:', nRight
     print *, '   Number of SigmaLeft components:', nLeft
     print *, '   Number of Pi components:', nPi
     Stop
  end if
  !
  ! Compute the Zeeman pattern, returning the splitting and strength of each component
  !
  call hyper_pattern(SpinI, L_up, S_up, Jup, Aup, Bup, &
       L_low, S_low, Jlow, Alow,Blow, Bfield, nflev_up, &
       energy_up, c_up, f2m_up, nfl_up, nflev_low, energy_low, c_low, f2m_low, &
       nfl_low, SplitPi, SplitRight, SplitLeft, StrPi, StrRight, StrLeft)
  !
  ! Clean-up
  !
  deallocate(nflev_low)
  deallocate(energy_low)
  deallocate(c_low)
  !
  deallocate(nflev_up)
  deallocate(energy_up)
  deallocate(c_up) 
  Call Time_routine('zeeman_hyperfine',.False.)
  !
END Subroutine Zeeman_hyperfine

