SUBROUTINE read_weights(Params,Sigma)
 Use Param_structure
  Implicit None
  Type (Parameters) :: Params
  Real, dimension(Params%n_data) :: Sigma
!
! Copied from the original compute_weights
!
  Call Read_profile(Params, Params%WeightsFilename, Sigma)
!
! Normalize weights so that their minimum is the noise
!
  Sigma(1:Params%n_data)=Sigma(1:Params%n_data)/MinVal(Sigma(1:Params%n_data))* &
       Params%Noise
!
! Discard values smaller than -10
!
!  Where (Profile .le. -10.) Sigma=1.e10 
!
END SUBROUTINE read_weights
  
