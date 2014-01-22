!        L O R I E N: the LOvely Reusable Inversion ENgine
!              (c) 1999 Hector Socas-Navarro
!
!     LORIEN is an inversion kernel consisting of a library of routines,
! that can be used to easily develop inversion codes for a variety of 
! purposes. Check the included documentation for instructions.
!
! Please, send bug reports to navarro@hao.ucar.edu
!
! Version control:
!
! v4.5:
!  * Encapsulated in a module
!
! v4.1:
!  * Changed compute_dydx to compute_dchisq_dx. Regularization terms are now
! incorporated in compute_chisq and considered in compute_dchisq_dx (but only
! in the first derivative)
!
! v3.1:
!  * Added check while computing derivatives to see if some parameters are
! getting out of range with the perturbation. If they are, then the 
! perturbation is subtracted instead of added.
!
! v2.01:
!  * Added history of lambda to stop useless iterations
!
! v2.0:
!  * Regularization. Deviations of the guess model from input model are 
! added to the Chisq
!  * Several IMPORTANT bug fixes
!  * Removed model smoothing (it was interfering with derivatives)
!
! v1.5:
!  * The structure (array) Region is taken as argument and passed to forward,
! exactly as with Line. 
!  * A NWChisq variable contains the Chisq without weights (although discarded
! points are still taken into account). Subroutine compute_chisq now returns
! both
!
! v1.4:
!  * Lorien now returns the final Chisq as an output argument
!  * Corrected several issues regarding the computation of the error bars
!  * A new user-provided routine expand_errors is called by compute_trial_model
! because sometimes the vector of errors in the free parameters needs to be
! expanded in a different way.
!  * Added stopping criterion. No more than 75 iterations allowed.
!  * Added call to the (user-provided) check_boundaries routine, which 
! should check that the current guess model does not exceed some preset
! values.
!  * The model_operations module is now moved to the model_struct.f90 file
! in the compex/ directory.
!  * Changed the shape and dimensions of the array Derivatives. Note that
! it should be accordingly changend in the routines Forward and Compute_dydx
!  * Added input argument Nodes, with information about nodes required by
! the routines Compress, Compress_deriv and Expand. The structure Nodes is 
! passed to these routines through Compute_dydx and Compute_trial_model.
!  * Fixed the declaration of the input variable Line (it must be an array),
! and passed it to Compute_dydx which needs it, in order to pass it to 
! Forward.
!  * Lambda, the diagonal element, is initialized to 10 instead of 1
!  * Changed the rejection criterion from Chisq .lt. Old_chisq to 
! Chisq .lt. 0.99*Old_chisq.
!
! v1.3:
!  * The structure line, which contains the required line data for solving 
! the forward problem is taken as input and passed to the forward routine.
!
! v1.2:
!  * Error bars are now given for the last accepted model, instead of just
! the last (usually rejected) model as in v1.1. Thanks to Rafael Manso for
! this one.
!
! v1.1: 
!  * Fixed two minor bugs present in v1.0
!  * Error bars now computed using the diagonal of the covariance matrix,
! instead of simply the inverse of the diagonal elements of the Hessian matrx.
!
Module Lorien_module
  Use Param_structure
  Use Compex_module
  Use Forward_module
  Use Model_structure
  Use Line_data_structure
  Use Nodes_information
  Use Profiling
  Use Debug_module
  Character (len=25) :: Lorien_ver='LORIEN Version 4.2'
!
Contains
  Subroutine Lorien(Params, Nodes, Line, Region, Guess_model, Obs_profile, &
     Sigma, Brute_force, Errors, Last_accepted_nwchisq)
  Implicit None
  Type (Parameters) :: Params
  Type (Line_data), dimension(Params%n_lines) :: Line
  Type (Region_data), dimension(Params%n_regions) :: Region
  Type (model_2comp) :: Guess_model, Trial_model, Errors, Trial_errors
  Type (Nodes_info) :: Nodes
  Real, Dimension (Params%n_data) :: Obs_profile, Syn_profile, Sigma, Trial_profile
  Logical, Dimension (Params%n_free_parameters) :: Brute_force
  Real, Dimension (Params%n_data, Params%n_variables, &
       Params%n_points) :: Derivatives
  Real, Dimension (Params%n_data, Params%n_free_parameters) :: Dydx
  Real, Dimension (Params%n_free_parameters) :: DChiDx
  Real, Dimension (Params%n_free_parameters, Params%n_free_parameters) :: &
       D2ChiD2x
  Real, Dimension (Params%max_inv_iters) :: HistLambda, HistChisq
  Real :: Lambda, Chisq, NWChisq, Last_accepted_chisq, Last_accepted_nwchisq
  Real :: Regul, Pertur
  Integer :: n_failed, n_iter
  Logical :: StopIter, Adjusted, CheckNaN
!
  Debug_errorflags(flag_lorien)=0
  Debug_warningflags(flag_lorien)=0
!
! Start with a coarse wavelength grid
  Adjusted=.False.
  Params%Skip_lambda=1
!
  Pertur=0.1
  Params%recompute_deriv=1
  Call Forward(Params, Line, Region, Guess_model, Syn_profile, .TRUE.)
  If (Debug_errorflags(flag_forward) .ge. 1) then
     Print *,'This guess model produces errors in the synthesis'
     Call Debug_Log('Lorien: This guess model produces errors in the synthesis',2)
     Last_accepted_nwchisq=1e10
     Return
  End if
  Call Compute_dchisq_dx(Params, Line, Region, Nodes, Brute_force, Pertur, Guess_model, &
       Obs_Profile, Syn_profile, Sigma, Derivatives, Dydx, DChiDx, D2ChiD2x)
  Call Compute_chisq(Params, Obs_profile, Syn_profile, Sigma, &
       Nodes, Guess_model, Chisq, NWChisq, Regul)
  If (Debug_errorflags(flag_dchisq_dx) .ge. 1) then 
     Print *,'This guess model produces errors in the derivatives'
     Call Debug_Log('Lorien: This guess model produces errors in the derivatives',2)
     Last_accepted_nwchisq=1e10
     Return
  End if
  HistChisq(1)=Chisq
  Last_accepted_chisq=Chisq
  Last_accepted_nwchisq=NWChisq
  Lambda=1e-3
  Call Printout(0, Lambda, Chisq, NWChisq, Params%Regularization*Regul, Params, &
       Line, Region, Guess_model, Guess_model, Syn_profile, .TRUE.)
  n_failed=0
  n_iter=1
  StopIter=.FALSE.
  Do while ( .not. StopIter )
     Call Compute_trial_model(Params, Nodes, Guess_model, Lambda, &
          Obs_profile, Syn_profile, NWChisq, Sigma, Dydx, Trial_model, &
          Trial_errors, Chisq, DChiDx, D2ChiD2x)
     If (n_iter .eq. 1) Errors=Trial_errors ! For the first iteration
     Call Check_boundaries(Params, Nodes, Guess_model%Comp1, Trial_model%Comp1)
     If (Params%TwoComp) &
          Call Check_boundaries(Params, Nodes, Guess_model%Comp2, Trial_model%Comp2)
     Call Forward(Params, Line, Region, Trial_model, Trial_profile,.TRUE.)
     Call Compute_chisq(Params, Obs_profile, Trial_profile, Sigma, &
          Nodes, Trial_model, Chisq, NWChisq, Regul)
     If (Debug_errorflags(flag_forward) .ge. 1) Chisq=1e10
     HistChisq(n_iter)=Chisq
     HistLambda(n_iter)=Lambda
     If (n_iter .ge. 6) then
        If (Minval(HistLambda(n_iter-5:n_iter)) .ge. 1) StopIter=.TRUE.
     End if
! Check that chisq is valid, otherwise reject model
     If (CheckNaN(Chisq)) &
          Chisq=1e10
! Adjust wavelength grid
     Call Adjust_wave_grid(Params, Chisq, Adjusted)
     If (Adjusted) then
        Call Forward(Params, Line, Region, Trial_model, Trial_profile,.FALSE.)
        Params%Recompute_deriv=1
     End if
! Refine derivatives
     If (n_iter .ge. 5) Pertur=0.01
!
     If (Chisq .lt. Last_accepted_chisq .or. Adjusted) then
        Last_accepted_chisq=Chisq
        Last_accepted_nwchisq=NWChisq
        If (.not. Adjusted) Params%recompute_deriv=0 
        If (Lambda .gt. 1.e-8) Lambda=Lambda/10.
        Guess_model=Trial_model 
        Errors=Trial_errors
        Syn_profile=Trial_profile
        Call Printout(n_iter, Lambda, Chisq, NWChisq, Params%regularization*Regul, &
             Params, Line, Region,Guess_model, Errors, Syn_profile, .TRUE.)
        If (n_iter .gt. 12) then
           If (Abs(MinVal(HistChisq(n_iter-10:n_iter-6))-Chisq) .lt. &
                0.01*Chisq) then
              If (Params%recompute_deriv .eq. 1) then 
                 StopIter=.TRUE.
                 Print *,'Convergence too slow. Giving up this cycle.'
              Else
                 Params%recompute_deriv=1
              End if
           End if
        End if
        Call Compute_dchisq_dx(Params, Line, Region, Nodes, Brute_force, &
             Pertur, Guess_model, Obs_profile, Syn_profile, Sigma, Derivatives, Dydx, &
             DchiDx, D2ChiD2x)
        n_failed=0
     Else
        Lambda=Lambda*10
        If (Lambda .lt. 1.e-4) Lambda=1.e-4
        If (Lambda .gt. 1.) Lambda=1.
        Call Printout(n_iter, Lambda, Chisq, NWChisq, Params%regularization*Regul, &
             Params, Line, Region,Trial_model, Trial_errors, Syn_profile, .False.)
        If (Pertur .lt. 0.09) n_failed=n_failed+1
        If (n_failed .eq. 1 .and. Params%recompute_deriv .eq. 0) then
           Params%recompute_deriv=1
           Call Compute_dchisq_dx(Params, Line, Region, Nodes, Brute_force, &
                Pertur, Guess_model, Obs_profile, Syn_profile, Sigma, Derivatives, Dydx, &
                DchiDx, D2ChiD2x)
        End if
     End if
     If (n_failed .gt. 3 .and. Pertur .lt. 0.09) StopIter=.True.
     If (n_iter .ge. Params%max_inv_iters) StopIter=.True.
     n_iter=n_iter+1
  End do
  Params%Skip_lambda=1 
  Call Deallocate_model_2comp(Trial_model)
  Call Deallocate_model_2comp(Trial_errors)
  Call Deallocate_model_2comp(Errors)
  Return
!
End Subroutine Lorien
! This subroutine evaluates Chisq for the current guess model
!
Subroutine Compute_chisq(Params, Obs_profile, Syn_profile, Sigma, &
     Nodes, Guess_model, Chisq, NWChisq, Deviation)
  Implicit None
  Type (Parameters) :: Params
  Type (model_2comp) :: Guess_model
  Type (Nodes_info) :: Nodes
  Real, Dimension (Params%n_data) :: Obs_profile, Syn_profile, Sigma
  Real :: Chisq, NWChisq, Sig, Deviation
  Integer :: i, ndata
!
!  Sig=Params%Noise
  Chisq=0.
  NWChisq=0.
  ndata=0
  Do i=1,Params%n_data
     If (Sigma(i) .lt. 1.e8 .and. Syn_profile(i) .gt. -1e9) then
        Chisq=Chisq+((Obs_profile(i)-Syn_profile(i))**2)/(Sigma(i)**2)
!        NWChisq=NWChisq + &
!             ((Obs_profile(i)-Syn_profile(i))**2)/(Sig**2)
        ndata=ndata+1
     End if
  End do
  Chisq=Chisq/ndata
  NWChisq=Chisq
!
! Add regularization
!
  Call Regul_term(Params, Nodes, Guess_model, Deviation)
!  Chisq=Chisq*(1.+Params%Regularization*Deviation)
  Chisq=Chisq+Params%Regularization*Deviation
!
  Return
!
End subroutine Compute_chisq
! This subroutine fills in the arrays dchidx and d2chid2x with the derivatives 
! of the merit function chi-squred with respect to the free parameters.
!
Subroutine Compute_dchisq_dx(Params, Line, Region, Nodes, Brute_force, &
     Pertur, Guess_model, Obs_Profile, Syn_Profile, Sigma, &
     Derivatives, Dydx, DChiDx, D2ChiD2x)
  Use Debug_module
  Implicit None
  Type (Parameters) :: Params
  Type (Line_data), dimension(Params%n_lines) :: Line
  Type (Region_data), dimension(Params%n_regions) :: Region
  Type (model_2comp) :: Guess_model, Pert_atmo
  Type (Nodes_info) :: Nodes
  Logical, Dimension (Params%n_free_parameters) :: Brute_force
  Real, Dimension (Params%n_data, Params%n_variables, &
       Params%n_points) :: Derivatives, &
       Derivatives2
  Real, Dimension (Params%n_data, Params%n_free_parameters) :: Dydx
  Real, Dimension (Params%n_free_parameters) :: DChiDx, DRegulDx
  Real, Dimension (Params%n_free_parameters,Params%n_free_parameters):: D2ChiD2x
  Real, Dimension (:), Allocatable, SAVE :: DChiDx_saved
  Real, Dimension (:,:), Allocatable, SAVE :: Dydx_saved, D2ChiD2x_saved
  Real, Dimension (Params%n_data) :: Obs_Profile, Syn_Profile, Pert_profile, &
       Sigma
  Real, Dimension (Params%n_free_parameters) :: X
  Real :: OldX, Chisq0, Chisq, NWChisq, Sigma2, Pertur, Regul, Pert_Regul
  Integer :: i_param, j_param, i_data, ndata
  Logical :: ForwError
!
  Call Time_routine('compute_dchisq_dx',.True.)
  Debug_errorflags(flag_dchisq_dx)=0
  Debug_warningflags(flag_dchisq_dx)=0
  If (.not. Allocated(Dydx_saved)) then
     If (Params%n_free_parameters .gt. 100) then
        Print *,'Need to redimension Dydx_saved in compute_dydx.f90'
        print *,'Number of free parameters:',Params%n_free_parameters
        Stop
     End if
     Allocate(Dydx_saved(Params%n_data, 100))
     Allocate(DChiDx_saved(Params%n_data))
     Allocate(D2ChiD2X_saved(Params%n_data, Params%n_data))
  End if
  DyDx(:,:)=0.
  If (Params%always_compute_deriv .eq. 1 .or. &
       Params%recompute_deriv .eq. 1) then
     Nodes%Reference_model=Guess_model ! Model_assign operation
     Call Compress(Params, Nodes, Guess_model, X)
     Call Compress_deriv(Params, Nodes, Derivatives, Dydx)
     Call Compute_chisq(Params, Obs_profile, Syn_profile, Sigma, &
          Nodes, Guess_model, Chisq0, NWChisq, Regul)
     Do i_param=1,Params%n_free_parameters
        If (params%printout .ge. 3) print *,'Computing deriv of parameter:',i_param
        If (Brute_force(i_param)) then
           OldX=X(i_param)
           If (X(i_param)+Pertur .lt. X_max(i_param)) then
              X(i_param)=X(i_param)+Pertur ! X is dimensionless and ~1
           Else
              X(i_param)=X(i_param)-Pertur 
           End if
           Call Expand(Params, Nodes, X, Pert_atmo)
           Call Check_boundaries(Params, Nodes, Guess_model%Comp1, Pert_atmo%Comp1)
           If (Params%TwoComp) & 
                Call Check_boundaries(Params, Nodes, Guess_model%Comp2, Pert_atmo%Comp2)
           Call Forward(Params, Line, Region, Pert_atmo, Pert_profile, &
                .TRUE.)
           If (Debug_errorflags(flag_forward) .ge. 1) &
                Debug_errorflags(flag_dchisq_dx)=1
           Call Compute_chisq(Params, Obs_profile, Pert_profile, Sigma, &
                Nodes, Pert_atmo, Chisq, NWChisq, Pert_regul)
           DChiDx(i_param)=(Chisq-Chisq0)/(X(i_param)-OldX)
           Do i_data=1,Params%n_data
              Dydx(i_data, i_param)=(Pert_profile(i_data)-Syn_Profile(i_data))/ &
                   (X(i_param)-OldX)
           End do
           DRegulDx(i_param)=Params%regularization* &
                (Pert_Regul-Regul)/(X(i_param)-OldX)
           X(i_param)=OldX
        Else
           Print *,'Stop in compute_dydx'
           Stop
        End if
     End do
! D2chid2x neglects additional contributions to chisq other than 
! weighted profile differences, according to Eq 15.5.7
     Do i_param=1,Params%n_free_parameters
        Do j_param=1,Params%n_free_parameters
           D2chid2x(i_param, j_param)=0.
           ndata=0
           Do i_data=1,Params%n_data
              If (Sigma(i_data) .lt. 1.e8 .and. Syn_profile(i_data) .gt. -1e9) then
                 ndata=ndata+1
                 Sigma2=Sigma(i_data)*Sigma(i_data)
                 D2chid2x(i_param, j_param)=D2chid2x(i_param, j_param) &
                      + 2.*Dydx(i_data, i_param)* & 
                      Dydx(i_data, j_param)/Sigma2
              End if
           End do
           D2chid2x(i_param, j_param)=D2chid2x(i_param, j_param)/ndata
           D2chid2x(i_param, j_param)=D2chid2x(i_param, j_param)+ &
                DRegulDx(i_param)*DRegulDx(j_param)
        End do
     End do
     Dydx_saved(1:Params%n_data,1:Params%n_free_parameters) = &
          Dydx(1:Params%n_data,1:Params%n_free_parameters)
     DchiDx_saved(1:Params%n_free_parameters) = &
          DchiDx(1:Params%n_free_parameters)
     D2chiD2x_saved(1:Params%n_free_parameters, 1:Params%n_free_parameters) = &
          D2chiD2x(1:Params%n_free_parameters, 1:Params%n_free_parameters)
  Else
     If (Params%printout .ge. 3) print *,'restoring derivatives'
     Dydx(1:Params%n_data,1:Params%n_free_parameters) = &
          Dydx_saved(1:Params%n_data,1:Params%n_free_parameters)
     DchiDx(1:Params%n_free_parameters) = &
          DchiDx_saved(1:Params%n_free_parameters)
     D2chiD2x(1:Params%n_free_parameters, 1:Params%n_free_parameters) = &
          D2chiD2x_saved(1:Params%n_free_parameters, 1:Params%n_free_parameters)
  End if
  Call Time_routine('compute_dchisq_dx',.False.)
  Return
!
End Subroutine Compute_dchisq_dx
! This routines computes the corrections to the current model.
!
Subroutine Compute_trial_model(Params, Nodes, Guess_model, Lambda, &
     Obs_profile, Syn_profile, NWChisq, Sigma, Dydx, Trial_model, Errors, &
     WChisq, DChiDx, D2ChiD2x)
  Implicit None
  Type (Parameters) :: Params
  Type (model_2comp) :: Guess_model, Trial_model, Errors
  Type (Nodes_info) :: Nodes
  Real, Dimension (Params%n_data, Params%n_free_parameters) :: Dydx
  Real, Dimension (Params%n_free_parameters) :: X, X_trial, DeltaX, Dchidx,&
       Beta_orig, Beta, X_errors, DRegulDx
  Real, Dimension (Params%n_free_parameters, Params%n_free_parameters) :: &
       Alpha_orig, Alpha, D2chid2x, D2RegulD2x
  Real, Dimension (Params%n_data) :: Obs_profile, Syn_profile, Sigma
  Real :: Lambda, Sigma2, Chisq, NWChisq, WChisq, Regul
  Integer :: i_param, j_param, i_data, ndata
  Logical, Dimension (Params%n_free_parameters) :: Zeroed
!
  Nodes%Reference_model=Guess_model ! Model_assign operation
  Call Compress(Params, Nodes, Guess_model, X)
  Chisq=WChisq
!
  Do i_param=1,Params%n_free_parameters
     Beta_orig(i_param)=-.5*Dchidx(i_param)
     Do j_param=1,Params%n_free_parameters
        Alpha_orig(i_param, j_param)=.5*D2chid2x(i_param, j_param)
     End do
  End do
!
  Beta=Beta_orig
  Alpha=Alpha_orig
!
! Alpha has the diagonal multiplied by the lambda factor and it's used
! for the calculation of the new model. Alpha_orig is used for the
! calculation of the errors.
!
  Do i_param=1,Params%n_free_parameters
     Alpha(i_param, i_param)=Alpha(i_param, i_param)*(1.+Lambda)
  End do
  Call SVD_solve(Params%n_free_parameters, Params%SVD_threshold, &
       Alpha, Beta, DeltaX, Zeroed)
  Do i_param=1, Params%n_free_parameters
     X_trial(i_param)=X(i_param)+DeltaX(i_param)
  End do
! Call SVD again for the errors.
! Beta was destroyed in the call to SVD_solve, so we use Beta_orig now
!!$  Call SVD_solve(Params%n_free_parameters, 1.e-3, &
!!$       Alpha_orig, Beta_orig, DeltaX, Zeroed)
  Do i_param=1, Params%n_free_parameters
!     Use (A^-1)_ii
!     The diagonal of the covariance matrix is returned in Beta_orig
!
!     Use the 1./A_ii instead of (A^-1)_ii. Note that in this case, the
!     call to SVD_solve above must be commented out.

     If (Alpha_orig(i_param, i_param) .gt. 1.E-8) then
        Beta_orig(i_param)=1./Alpha_orig(i_param, i_param)
     Else
        Beta_orig(i_param)=1000
     End If
     X_errors(i_param)=Sqrt(Chisq*Beta_orig(i_param))
  End do
!
  Call Expand(Params, Nodes, X_trial, Trial_model)
!
  Call Expand_errors(Params, Nodes, Trial_model%Comp1, X_errors, Errors%Comp1)
  Call Check_boundaries(Params, Nodes, Guess_model%comp1, Trial_model%Comp1)
  If (Params%TwoComp) then
     Call Check_boundaries(Params, Nodes, Guess_model%comp2, Trial_model%Comp2)
     Call Expand_errors(Params, Nodes, Trial_model%Comp2, X_errors, Errors%Comp2)
  End if
  Return
!
End Subroutine Compute_trial_model
!
Subroutine Regul_term(Params, Nod, Atmo_2comp, Deviation)
  Implicit None
  Type (Parameters) :: Params
  Type (model_2comp) :: Atmo_2comp
  Type (model) :: Atmo, Ref
  Type (Nodes_info) :: Nod
  Real :: Deviation, Mean, Norm, sigma, x1, expo
  Integer :: i_param, inode, ifree, nodes, idepth, idepth2
  Real, Dimension(Params%n_free_parameters) :: X, tau
  Real, Dimension(Params%n_free_parameters, Params%n_free_parameters) :: &
       PRD2, D2RegulD2x
  Real, Dimension(Params%n_points) :: xx, y, y2, kernel
  Real, Dimension(10) :: Reguls, Regul_weights

!
  Reguls(:)=0.
  Regul_weights(:)=1.
  ! Component 1
  Atmo=Atmo_2comp%Comp1
  ! 
  ! Velocity. Square deviation from the mean
  !
  Mean=Sum(Atmo%v_los)/Params%n_points
  y=(Atmo%v_los-Mean)/1e5
  Reguls(1)=0.
  Do idepth=2, params%n_points
     Reguls(1)=Reguls(1)+Abs( y(idepth)*(Atmo%ltau_500(idepth)-Atmo%ltau_500(idepth-1)) )
  End do
  Reguls(1)=Sqrt(Reguls(1))
  Regul_weights(1)=10.
  !
  ! Temperature increase
  !
  Do idepth=2, params%n_points
     if (atmo%temp(idepth) .lt. atmo%temp(idepth-1)) &
          Reguls(2)=Reguls(2)+(atmo%temp(idepth-1)-atmo%temp(idepth))* &
          (atmo%ltau_500(idepth)-atmo%ltau_500(idepth-1))/500.
  end do
  Regul_weights(2)=1.
  !
  ! Filling factor
  !
  Reguls(3)=1.-Atmo%ffactor
  Regul_weights(3)=1.
  !
  ! Smooth T
  !
  y=Atmo%Temp
  xx=Atmo%ltau_500
  y2(:)=0.
  sigma=1. ! Half-width of smoothing Gaussian (in ltau units)
  Do idepth=1, Params%n_points
     Norm=0.
     Do idepth2=1, Params%n_points
        x1=xx(idepth2)-xx(idepth)
        If (x1*x1/sigma/sigma .lt. 25.) then ! closer than 5-sigma
           expo=exp(-x1*x1/sigma/sigma)
           y2(idepth)=y2(idepth)+y(idepth2)*expo
           Norm=Norm+expo
        End if
     End do
     y2(idepth)=y2(idepth)/Norm
  End do
  y=Atmo%Temp-y2
  Reguls(4)=Sqrt(Sum( y**2 ) )*(atmo%ltau_500(Params%n_points)-atmo%ltau_500(1))/500.
  Regul_weights(4)=0.01
  !
  ! v_mic
  !
  Reguls(5)=Sum( .1*(10**(Atmo%ltau_500(2:Params%n_points)- &
       Atmo%ltau_500(1:Params%n_points-1)))*Atmo%v_mic(1:Params%n_points-1)*1e-5 ) ! v_mic in km/s
  Regul_weights(5)=0.1
  ! 
  ! B. Square deviation from the mean
  !
  Reguls(6)=0.

  Mean=Sum(Atmo%B_long)/Params%n_points
  y=(Atmo%B_long-Mean)
  Do idepth=2, params%n_points
     Reguls(6)=Reguls(6)+Abs( y(idepth)*(Atmo%ltau_500(idepth)-Atmo%ltau_500(idepth-1)) )
  End do

  Mean=Sum(Atmo%B_x)/Params%n_points
  y=(Atmo%B_x-Mean)
  Do idepth=2, params%n_points
     Reguls(6)=Reguls(6)+Abs( y(idepth)*(Atmo%ltau_500(idepth)-Atmo%ltau_500(idepth-1)) )
  End do

  Mean=Sum(Atmo%B_y)/Params%n_points
  y=(Atmo%B_y-Mean)
  Do idepth=2, params%n_points
     Reguls(6)=Reguls(6)+Abs( y(idepth)*(Atmo%ltau_500(idepth)-Atmo%ltau_500(idepth-1)) )
  End do

  Reguls(6)=Sqrt(Reguls(6))
  Regul_weights(6)=10.

  !
  !
  ! Component 2
  !
  !
  Atmo=Atmo_2comp%Comp2
  ! 
  ! Velocity. Square deviation from the mean
  !
  Mean=Sum(Atmo%v_los)/Params%n_points
  y=(Atmo%v_los-Mean)/1e5
  Do idepth=2, params%n_points
     Reguls(1)=Reguls(1)+ Abs( y(idepth)*(Atmo%ltau_500(idepth)-Atmo%ltau_500(idepth-1)) )
  End do
  !
  ! Temperature increase
  !
  Do idepth=2, params%n_points
     if (atmo%temp(idepth) .lt. atmo%temp(idepth-1)) &
          Reguls(2)=Reguls(2)+(atmo%temp(idepth-1)-atmo%temp(idepth))* &
          (atmo%ltau_500(idepth)-atmo%ltau_500(idepth-1))/500.
  end do
  !
  ! Smooth T
  !
  y=Atmo%Temp
  xx=Atmo%ltau_500
  y2(:)=0.
  sigma=1. ! Half-width of smoothing Gaussian (in ltau units)
  Do idepth=1, Params%n_points
     Norm=0.
     Do idepth2=1, Params%n_points
        x1=xx(idepth2)-xx(idepth)
        If (x1*x1/sigma/sigma .lt. 25.) then ! closer than 5-sigma
           expo=exp(-x1*x1/sigma/sigma)
           y2(idepth)=y2(idepth)+y(idepth2)*expo
           Norm=Norm+expo
        End if
     End do
     y2(idepth)=y2(idepth)/Norm
  End do
  y=Atmo%Temp-y2
  Reguls(4)=Reguls(4) + Sqrt(Sum( y**2 ) )*(atmo%ltau_500(Params%n_points)-atmo%ltau_500(1))/500.
  !
  ! v_mic
  !
  Reguls(5)=Reguls(5) + Sum( .1*(10**(Atmo%ltau_500(2:Params%n_points)- &
       Atmo%ltau_500(1:Params%n_points-1)))*Atmo%v_mic(1:Params%n_points-1)*1e-5 )
  ! 
  ! B. Square deviation from the mean
  !

  Mean=Sum(Atmo%B_long)/Params%n_points
  y=(Atmo%B_long-Mean)
  Do idepth=2, params%n_points
     Reguls(6)=Reguls(6)+Abs( y(idepth)*(Atmo%ltau_500(idepth)-Atmo%ltau_500(idepth-1)) )
  End do

  Mean=Sum(Atmo%B_x)/Params%n_points
  y=(Atmo%B_x-Mean)
  Do idepth=2, params%n_points
     Reguls(6)=Reguls(6)+Abs( y(idepth)*(Atmo%ltau_500(idepth)-Atmo%ltau_500(idepth-1)) )
  End do

  Mean=Sum(Atmo%B_y)/Params%n_points
  y=(Atmo%B_y-Mean)
  Do idepth=2, params%n_points
     Reguls(6)=Reguls(6)+Abs( y(idepth)*(Atmo%ltau_500(idepth)-Atmo%ltau_500(idepth-1)) )
  End do

  Reguls(6)=Sqrt(Reguls(6))
  Regul_weights(6)=10.


  Deviation=Sum(Reguls*Regul_weights)
!  print *,'reguls=',reguls(1:4)*regul_weights(1:4)


  Return
!
End Subroutine Regul_term

Subroutine Regul_variable(nn, X, Y, n, PR, PRD, PRD2)
  Implicit None
  Integer :: nn, n, i, j, k, l, inode, nodes
  Real, Dimension(nn) :: X, Y, PRD
  Real, Dimension(nn, nn) :: PRD2
  Real :: PR, Y0, X0, Mean

  nodes=n
  If (n .le. 1) then ! No regularization is necessary
     PR=0.
     PRD(:)=0.
     PRD2(:,:)=0.
     Return
  End if
  If (n .le. 3) then ! Penalize deviations from a constant (use Sigma^2)
     Mean=Sum(Y(1:nodes))/nodes
     Do i=1, nodes
        PR=(Y(i)-Mean)**2/(nodes-1)
        PRD(i)=2./(nodes-1)*( (Y(i)-Mean) - Sum(Y(1:nodes)-Mean)/nodes )
        Do j=1, nodes
           PrD2(i, j)=-2./nodes/(nodes-1)
        End do
        PRD2(i, i)=2./nodes ! Diagonal terms are different
     End do
  End if
  If (n .ge. 4) then ! Piece-wise penalize deviations from straight line
     PR=0.
     PRD(:)=0.
     PRD2(:,:)=0.
! First point
     X0=(X(1)-X(2))/(X(3)-X(2))
     Y0=(Y(3)-Y(2))*X0+Y(2)
     PR=PR+(Y(1)-Y0)**2
     PRD(1)=PRD(1)+2.*(Y(1)-Y0)
     PRD(2)=PRD(2)+2.*(Y(1)-Y0)*(1*X0 - 1)
     PRD(3)=PRD(3)+2.*(Y(1)-Y0)*(-1*X0)
     PRD2(1,1)=PRD2(1,1)+2.
     PRD2(1,2)=PRD2(1,2)+2.*(1.*X0 - 1)
     PRD2(1,3)=PRD2(1,3)+2.*(-1*X0)
     PRD2(2,1)=PRD2(2,1)+2.*(1)*(1*X0 - 1)
     PRD2(2,2)=PRD2(2,2)+2.*(1*X0 - 1)*(1*X0 - 1)
     PRD2(2,3)=PRD2(2,3)+2.*(-1*X0)*(1*X0 - 1)
     PRD2(3,1)=PRD2(3,1)+2.*(1)*(-1*X0)
     PRD2(3,2)=PRD2(3,2)+2.*(1*X0 - 1)*(-1*X0)
     PRD2(3,3)=PRD2(3,3)+2.*(-1*X0)*(-1*X0)
! Last point
     X0=(X(nodes)-X(nodes-1))/(X(nodes-2)-X(nodes-1))
     Y0=(Y(nodes-2)-Y(nodes-1))*X0+Y(nodes-1)
     PR=PR+(Y(nodes)-Y0)**2
     PRD(nodes)=PRD(nodes)+2.*(Y(nodes)-Y0)
     PRD(nodes-1)=PRD(nodes-1)+2.*(Y(nodes)-Y0)*(1*X0 - 1)
     PRD(nodes-2)=PRD(nodes-2)+2.*(Y(nodes)-Y0)*(-1*x0)
     PRD2(nodes,nodes)=PRD2(nodes,nodes)+2.
     PRD2(nodes,nodes-1)=PRD2(nodes,nodes-1)+2.*(1.*X0 - 1)
     PRD2(nodes,nodes-2)=PRD2(nodes,nodes-2)+2.*(-1*X0)
     PRD2(nodes-1,nodes)=PRD2(nodes-1,nodes)+2.*(1)*(1*X0 - 1)
     PRD2(nodes-1,nodes-1)=PRD2(nodes-1,nodes-1)+2.*(1*X0 - 1)*(1*X0 - 1)
     PRD2(nodes-1,nodes-2)=PRD2(nodes-1,nodes-2)+2.*(-1*X0)*(1*X0 - 1)
     PRD2(nodes-2,nodes)=PRD2(nodes-2,nodes)+2.*(1)*(-1*X0)
     PRD2(nodes-2,nodes-1)=PRD2(nodes-2,nodes-1)+2.*(1*X0 - 1)*(-1*X0)
     PRD2(nodes-2,nodes-2)=PRD2(nodes-2,nodes-2)+2.*(-1*X0)*(-1*X0)
! Intermediate points
     Do i=2, nodes-1
        X0=(X(i)-X(i-1))/(X(i+1)-X(i-1))
        Y0=(Y(i+1)-Y(i-1))*X0+Y(i-1)
        PR=PR+(Y(i)-Y0)**2
        PRD(i)=PRD(i)+2.*(Y(i)-Y0)
        PRD(i-1)=PRD(i-1)+2.*(Y(i)-Y0)*(1*X0 - 1)
        PRD(i+1)=PRD(i+1)+2.*(Y(i)-Y0)*(-1*X0)
        PRD2(i,i)=PRD2(i,i)+2.
        PRD2(i,i-1)=PRD2(i,i-1)+2.*(1*X0 - 1)
        PRD2(i,i+1)=PRD2(i,i+1)+2.*(-1*X0)
        PRD2(i-1,i)=PRD2(i-1,i)+2.*(1*X0 - 1)
        PRD2(i-1,i-1)=PRD2(i-1,i-1)+2.*(1*X0-1)*(1*X0 - 1)
        PRD2(i-1,i+1)=PRD2(i-1,i+1)+2.*(-1*X0)*(1*X0 - 1)
        PRD2(i+1,i)=PRD2(i+1,i)+2.*(-1*X0)
        PRD2(i+1,i-1)=PRD2(i+1,i-1)+2.*(1*X0-1)*(-1*X0)
        PRD2(i+1,i+1)=PRD2(i+1,i+1)+2.*(-1*X0)*(-1*X0)
     End do
  End if

     
  Return
End Subroutine Regul_variable

End Module Lorien_module
