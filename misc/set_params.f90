! Subroutine to fill in the Params structure
!
Subroutine Set_params(Input, Params)
!
  Use NICOLE_inp
  Use Param_structure
  Implicit None
  Type (NICOLE_input) :: Input
  Type (Parameters) :: Params
  Integer File_unit
  Logical Exists
  Character (len = 256) :: Line
!
  If (Input%mode .eq. 'i') Params%NLTE_Elim1=1e-4
  Params%helioazim=Input%helioazim
  Params%SVD_threshold=1.e-5 ! This threshold is hardwired
!  Params%Regularization=1. ! Read from input file
  Params%Regul_temp=.5 ! Hardwired
  Params%Regul_stray=.5 ! Hardwired
  Params%ncycles=2 ! Default. May be changed by nodes.dat file
  Params%Skip_lambda=1 ! By default, use full wavelength grid

  Return
End Subroutine Set_params
