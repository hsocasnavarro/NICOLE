Module Param_structure
  Use Phys_constants
!
! Hardware-dependent parameters
!
! Recl parameter for binary read input/output when reading/writing a 64-bit
! real number
  Integer, Parameter :: RealBytes=8 ! DO NOT MODIFY THIS LINE!!!
! Machine endianness
!  True=Little Endian, e.g. i386
!  False=Big endian, e.g. Motorola, PowerPC
!  Logical, Parameter :: LittleEndian=.True. ! DO NOT MODIFY THIS LINE
  Logical, Parameter :: LittleEndian=.True. ! DO NOT MODIFY THIS LINE
!
! Run-dependent parameters
!
  Type Parameters
     Integer :: n_data, n_free_parameters, n_variables
     Integer :: max_inv_iters
     Integer :: n_points, n_lines, n_regions, ncycles
     Integer :: formal_solution, formal_boundary_cond, printout
     Integer :: nPix, def_abund, speed, always_compute_deriv, cent_der
     Integer :: recompute_deriv, reference_cont
     Integer :: Skip_lambda, Reinterpolate
     Real :: SVD_threshold, heliocentric, helioazim
     Real :: Noise, Update_opac, Negligible_opacity
     Real :: Regularization, Regul_temp, Regul_stray, Regul_mic
     Logical :: IProfExists, WeightsExist
     Logical :: TwoComp=.False.
     Integer :: mmax
     Integer, Dimension(:), Allocatable :: mm
     Real, dimension(:), allocatable :: IProf
     Real, Dimension (:), pointer :: Stray_prof ! For the stray light
     Character (len=20) :: WeightsFilename 
     character (len = 4) :: hscale
     Real :: NLTE_Elim1=1e-3, NLTE_Elim2=1e-3
     Real ::  NLTE_QNORM, NLTE_CPER, NLTE_OptThin, NLTE_OptThick
     Integer :: NLTE_linear
     Integer :: NLTE_NMU, NLTE_MaxIters, NLTE_Verbose, NLTE_ISum, &
          NLTE_ltepop
     Integer :: NLTE_IStart, NLTE_UseColSwitch, &
          NLTE_NumLambdaIters
     Integer :: NLTE_formal_solution
     Logical :: NLTE_VelFree, NLTE_NGacc
     character (len = 4) :: input_dens     
  End Type Parameters
! Maximum allowed Zeeman components per spectral line
  Integer, Parameter :: transitions=5000
End Module Param_structure
