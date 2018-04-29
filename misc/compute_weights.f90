! This subroutine computes the weights for the chi^2 function so that all
! the Stokes parameters have the same amplitude. If icycle .eq. 1 then I
! is given Norm_I times more weight than Q, U, V.
!
Subroutine Compute_weights(Params, Profile, Sigma, icycle)
  Use Param_structure
  Implicit None
  Type (Parameters) :: Params
  Real, dimension(Params%n_data) :: Profile, AbsProfile, Sigma
  Real :: Imax, Qmax, Umax, Vmax
  Real :: Norm_I, MinAmpli
  Integer :: ind, status, icycle, iunit, Get_lun
  Logical :: Exists
  Logical, Save :: FirstTime=.True.
!
!
!
  Sigma=0. ! Initialization of Sigma.
!!$!
!!$! Is there a weights.pro file?
!!$!
!!$  Inquire (File='Weights.pro', Exist=Exists)
!!$  If (Exists) then
!!$     Call Read_profile(Params, 'Weights.pro', Sigma)
!!$!
!!$! Normalize weights so that their minimum is the noise
!!$!
!!$     Sigma(1:Params%n_data)=Sigma(1:Params%n_data)/MinVal(Sigma(1:Params%n_data))* &
!!$          Params%Noise
!!$!
!!$! Discard values smaller than -10
!!$!
!!$     Where (Profile .le. -10.) Sigma=1.e10 
!!$     Return
!!$  End if
!
! If there's no file, use standard weights except where observed profile
! is lt -10 (these points will be discarded)
!
  Norm_I=10. ! Weight for I is Norm_I times Q, U, V. It may be used to
!  If (icycle .eq. 1) Norm_I=1. ! give I more weight in the first iterations.
!
! Find the largest values of the different Stokes parameters.
!
  AbsProfile=Abs(Profile)
  Imax=0.
  Qmax=0.
  Umax=0.
  Vmax=0.
  ind=1
  Do while (ind .lt. Params%n_data)
     If (AbsProfile(ind) .gt. Imax .and. Profile(ind) .gt. -10) &
          Imax=AbsProfile(ind)
     If (AbsProfile(ind+1) .gt. Qmax .and. Profile(ind+1) .gt. -10) &
          Qmax=AbsProfile(ind+1)
     If (AbsProfile(ind+2) .gt. Umax .and. Profile(ind+2) .gt. -10) &
          Umax=AbsProfile(ind+2)
     If (AbsProfile(ind+3) .gt. Vmax .and. Profile(ind+3) .gt. -10) &
          Vmax=AbsProfile(ind+3)
     ind=ind+4
  End do
  If (Imax .lt. Params%Noise*100) Imax=Params%Noise*100
  If (Qmax .lt. Params%Noise*100) Qmax=Params%Noise*100
  If (Umax .lt. Params%Noise*100) Umax=Params%Noise*100
  If (Vmax .lt. Params%Noise*100) Vmax=Params%Noise*100
  If (Imax .gt. 1e7) Imax=1e6
!
! Fill in the Sigma array with the proper weights
!
  Do ind=1, Params%n_data-3, 4
     Sigma(ind)=Imax*.01
     Sigma(ind+1)=Qmax*.01
     Sigma(ind+2)=Umax*.01
     Sigma(ind+3)=Vmax*.01
  End do
!
! Normalize weights so that their minimum is the noise
!
  Sigma(1:Params%n_data)=Sigma(1:Params%n_data)/MinVal(Sigma(1:Params%n_data))* &
       Params%Noise
!
! Discard values smaller than -10
!
  Where (Profile .le. -10.) Sigma=1.e10 
  If (FirstTime) then
     FirstTime=.False.
     iunit=Get_lun()
     Open(Unit=iunit, File='Weights_used.pro')
     Do ind=1, Params%n_data-3, 4
        Write (iunit, *) 0., Sigma(ind), Sigma(ind+1), Sigma(ind+2), Sigma(ind+3)
     End do
     Close (iunit)
  End if
  Return
End Subroutine Compute_weights
