Module Compex_module
  Character (len=25) :: Compex_ver='NICOLE Compex v3.5'

Contains
! This subroutine performs a compression of the atmosphere onto the
! free parameters vector X. The free parameters are the difference
! between the guess model and the reference model stored in 
! Nodes%Reference_model (it's the starting model for the inversion),
! normalized to the values in the Nodes_information module.
!
Subroutine Compress(Params, Nodes, Guess_model_2comp, X)
  Use Param_structure
  Use Model_structure
  Use Nodes_information
  Implicit None
  Type (Parameters) :: Params
  Type (Model_2comp) :: Guess_model_2comp
  Type (Model) :: Guess_model, Ref
  Type (Nodes_info) :: Nodes
  Real, Dimension (Params%n_free_parameters) :: X, Node_values
  Real, Dimension (Params%n_points) :: y
  Integer :: ifree, npoints, idx
!
  If (Params%n_free_parameters .gt. 100) then
     print *,'Free parameters: ',Params%n_free_parameters
     print *,'Need to change definition of X_max in nodes_info.f90'
     stop
  End if
  npoints=Params%n_points
  ifree=1
!
! Component 1
!
  Guess_model=Guess_model_2comp%Comp1
  Ref=Nodes%Reference_model%Comp1
! Temperature
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%temp, &
       Nodes%n_nodes_t, Nodes%i_nodes_t, Ref%temp, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_t ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_t-1)=Node_values(1:Nodes%n_nodes_t)
  X_max(ifree:ifree+Nodes%n_nodes_t-1)=Max_t/Norm_t
  ifree=ifree+Nodes%n_nodes_t
! Velocity
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%v_los, &
       Nodes%n_nodes_v, Nodes%i_nodes_v, Ref%v_los, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_v ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_v-1)=Node_values(1:Nodes%n_nodes_v)
  X_max(ifree:ifree+Nodes%n_nodes_v-1)=Max_v/Norm_v
  ifree=ifree+Nodes%n_nodes_v
! Microturbulence
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%v_mic, &
       Nodes%n_nodes_mic, Nodes%i_nodes_mic, Ref%v_mic, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_mic ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_mic-1)=Node_values(1:Nodes%n_nodes_mic)
  X_max(ifree:ifree+Nodes%n_nodes_mic-1)=Max_mic/Norm_mic
  ifree=ifree+Nodes%n_nodes_mic
! Magnetic strength
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%b_long, &
       Nodes%n_nodes_blong, Nodes%i_nodes_blong, Ref%b_long, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_blong ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_blong-1)=Node_values(1:Nodes%n_nodes_blong)
  X_max(ifree:ifree+Nodes%n_nodes_blong-1)=Max_blong/Norm_blong
  ifree=ifree+Nodes%n_nodes_blong
! Field
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%b_x, &
       Nodes%n_nodes_bx, Nodes%i_nodes_bx, Ref%b_x, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_bx ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_bx-1)=Node_values(1:Nodes%n_nodes_bx)
  X_max(ifree:ifree+Nodes%n_nodes_bx-1)=Max_bx/Norm_bx
  ifree=ifree+Nodes%n_nodes_bx
! Field
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%b_y, &
       Nodes%n_nodes_by, Nodes%i_nodes_by, Ref%b_y, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_by ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_by-1)=Node_values(1:Nodes%n_nodes_by)
  X_max(ifree:ifree+Nodes%n_nodes_by-1)=Max_by/Norm_by
  ifree=ifree+Nodes%n_nodes_by
! Macroturbulence
  If (Nodes%n_nodes_mac .eq. 1) then
     X(ifree)=(Guess_model%v_mac-Ref%v_mac)/Norm_mac
     X_max(ifree)=Max_mac/Norm_mac
     ifree=ifree+1
  End if
! Percentage of stray light
  If (Nodes%n_nodes_stray .eq. 1) then
     X(ifree)=(Guess_model%stray-Ref%stray)/Norm_stray
     X_max(ifree)=(Max_stray-Ref%stray)/Norm_stray
     ifree=ifree+1
  End if
! Field ffactor
  If (Nodes%n_nodes_ffactor .eq. 1) then
     X(ifree)=(Guess_model%ffactor-Ref%ffactor)/ &
          Norm_ffactor
     X_max(ifree)=(Max_ffactor-Ref%ffactor)/Norm_ffactor
     ifree=ifree+1
  End if
! Abundances
  If (Nodes%n_nodes_ab .gt. 0) then
     Do idx=1, Nodes%n_nodes_ab
        Node_values(idx)=Guess_model%Abundance(Nodes%i_nodes_ab(idx)) - &
             Ref%Abundance(Nodes%i_nodes_ab(idx))
        X(ifree)=Node_values(idx)/Norm_ab
        X_max(ifree)=(Max_ab-Ref%Abundance(Nodes%i_nodes_ab(idx)))/Norm_ab
        ifree=ifree+1
     End do
  End if
!
  If (.not. Params%TwoComp) Return
!
! Component 2
!
  Guess_model=Guess_model_2comp%Comp2
  Ref=Nodes%Reference_model%Comp2
! Temperature
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%temp, &
       Nodes%n_nodes_t2, Nodes%i_nodes_t2, Ref%temp, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_t ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_t2-1)=Node_values(1:Nodes%n_nodes_t2)
  X_max(ifree:ifree+Nodes%n_nodes_t2-1)=Max_t/Norm_t
  ifree=ifree+Nodes%n_nodes_t2
! Velocity
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%v_los, &
       Nodes%n_nodes_v2, Nodes%i_nodes_v2, Ref%v_los, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_v ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_v2-1)=Node_values(1:Nodes%n_nodes_v2)
  X_max(ifree:ifree+Nodes%n_nodes_v2-1)=Max_v/Norm_v
  ifree=ifree+Nodes%n_nodes_v2
! Microturbulence
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%v_mic, &
       Nodes%n_nodes_mic2, Nodes%i_nodes_mic2, Ref%v_mic, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_mic ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_mic2-1)=Node_values(1:Nodes%n_nodes_mic2)
  X_max(ifree:ifree+Nodes%n_nodes_mic2-1)=Max_mic/Norm_mic
  ifree=ifree+Nodes%n_nodes_mic2
! Magnetic strength
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%b_long, &
       Nodes%n_nodes_blong2, Nodes%i_nodes_blong2, Ref%b_long, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_blong ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_blong2-1)=Node_values(1:Nodes%n_nodes_blong2)
  X_max(ifree:ifree+Nodes%n_nodes_blong2-1)=Max_blong/Norm_blong
  ifree=ifree+Nodes%n_nodes_blong2
! Field
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%b_x, &
       Nodes%n_nodes_bx2, Nodes%i_nodes_bx2, Ref%b_x, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_bx ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_bx2-1)=Node_values(1:Nodes%n_nodes_bx2)
  X_max(ifree:ifree+Nodes%n_nodes_bx2-1)=Max_bx/Norm_bx
  ifree=ifree+Nodes%n_nodes_bx2
! Field
  Call Compress_variable(npoints, Guess_model%ltau_500, Guess_model%b_y, &
       Nodes%n_nodes_by2, Nodes%i_nodes_by2, Ref%b_y, &
       Node_values) ! Extract the values at the nodes
  Node_values=Node_values/Norm_by ! Normalize these values.
  X(ifree:ifree+Nodes%n_nodes_by2-1)=Node_values(1:Nodes%n_nodes_by2)
  X_max(ifree:ifree+Nodes%n_nodes_by2-1)=Max_by/Norm_by
  ifree=ifree+Nodes%n_nodes_by2
! Abundances
  If (Nodes%n_nodes_ab2 .gt. 0) then
     Do idx=1, Nodes%n_nodes_ab2
        Node_values(idx)=Guess_model%Abundance(Nodes%i_nodes_ab2(idx)) - &
             Ref%Abundance(Nodes%i_nodes_ab2(idx))
        X(ifree)=Node_values(idx)/Norm_ab
        X_max(ifree)=(Max_ab-Ref%Abundance(Nodes%i_nodes_ab2(idx)))/Norm_ab
        ifree=ifree+1
     End do
  End if
!
! Done!
!
  Return
!
End Subroutine Compress
! This subroutine performs the expansion of the atmosphere from the
! free parameters vector X. The result is returned in New_model.
!
Subroutine Expand(Params, Nodes, X, New_model_2comp)
  Use Param_structure
  Use Model_structure
  Use Nodes_information
  Implicit None
  Type (Parameters) :: Params
  Type (Model_2comp) :: New_model_2comp
  Type (Model) :: New_model, Ref
  Type (Nodes_info) :: Nodes
  Real, Dimension (Params%n_free_parameters) :: X, Node_values
  Real, Dimension (Params%n_points) :: y, inc, azi
  Integer :: ifree, npoints, ind, idx
  Logical :: CheckNaN, Error
!
  npoints=Params%n_points
  ifree=1
!
! Component 1
!
  New_model=Nodes%Reference_model%Comp1
  Ref=Nodes%Reference_model%Comp1
! Temperature
  If (Nodes%n_nodes_t .gt. 0) then
     Node_values(1:Nodes%n_nodes_t)=X(ifree:ifree+Nodes%n_nodes_t-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_t, Nodes%i_nodes_t, y)
     New_model%temp=Ref%temp+y*Norm_t
     ifree=ifree+Nodes%n_nodes_t
  End if
! Velocity
  If (Nodes%n_nodes_v .gt. 0) then
     Node_values(1:Nodes%n_nodes_v)=X(ifree:ifree+Nodes%n_nodes_v-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_v, Nodes%i_nodes_v, y)
     New_model%v_los=Ref%v_los+y*Norm_v
     ifree=ifree+Nodes%n_nodes_v
  End if
! Microturbulence
  If (Nodes%n_nodes_mic .gt. 0) then
     Node_values(1:Nodes%n_nodes_mic)=X(ifree:ifree+Nodes%n_nodes_mic-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_mic, Nodes%i_nodes_mic, y)
     New_model%v_mic=Ref%v_mic+y*Norm_mic
     ifree=ifree+Nodes%n_nodes_mic
  End if
! B long
  If (Nodes%n_nodes_blong .gt. 0) then
     Node_values(1:Nodes%n_nodes_blong)=X(ifree:ifree+Nodes%n_nodes_blong-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_blong, Nodes%i_nodes_blong, y)
     New_model%b_long=Ref%b_long+y*Norm_blong
     ifree=ifree+Nodes%n_nodes_blong
  End if
! Bx
  If (Nodes%n_nodes_bx .gt. 0) then
     Node_values(1:Nodes%n_nodes_bx)=X(ifree:ifree+Nodes%n_nodes_bx-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_bx, Nodes%i_nodes_bx, y)
     New_model%b_x=Ref%b_x+y*Norm_bx
     ifree=ifree+Nodes%n_nodes_bx
  End if
! By
  If (Nodes%n_nodes_by .gt. 0) then
     Node_values(1:Nodes%n_nodes_by)=X(ifree:ifree+Nodes%n_nodes_by-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_by, Nodes%i_nodes_by, y)
     New_model%b_y=Ref%b_y+y*Norm_by
     ifree=ifree+Nodes%n_nodes_by
  End if
! Force Bx to be always positive by choosing that disambiguation direction
  Do ind=1, Params%n_points
     If (New_model%b_x(ind) .lt. 0) then
        New_model%b_x(ind)=-New_model%b_x(ind)
        New_model%b_y(ind)=-New_model%b_y(ind)
     End if
  End do
! Macroturbulence
  If (Nodes%n_nodes_mac .eq. 1) then
     New_model%v_mac=Ref%v_mac+X(ifree)*Norm_mac
     ifree=ifree+1
  End if
! Percentage of stray light
  If (Nodes%n_nodes_stray .eq. 1) then
     New_model%stray=Ref%stray+X(ifree)*Norm_stray
     ifree=ifree+1
  End if
! Field ffactor
  If (Nodes%n_nodes_ffactor .eq. 1) then
     New_model%ffactor=Ref%ffactor+X(ifree)*Norm_ffactor
     ifree=ifree+1
  End if
! Abundances
  If (Nodes%n_nodes_ab .gt. 0) then
     Do idx=1, Nodes%n_nodes_ab
        New_model%abundance(Nodes%i_nodes_ab(idx))=Ref%Abundance(Nodes%i_nodes_ab(idx)) + &
             X(ifree)*Norm_ab
        ifree=ifree+1
     End do
  End if
!
  New_model_2comp%Comp1=New_model
  Call Hydrostatic(Params, New_model_2comp%Comp1, Error)
  If (.not. Params%TwoComp) then
     New_model_2comp%Comp2=New_model_2comp%Comp1
     Return
  End if
!
! Component 2
!
  New_model=Nodes%Reference_model%Comp2
  Ref=Nodes%Reference_model%Comp2
! Temperature
  If (Nodes%n_nodes_t2 .gt. 0) then
     Node_values(1:Nodes%n_nodes_t2)=X(ifree:ifree+Nodes%n_nodes_t2-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_t2, Nodes%i_nodes_t2, y)
     New_model%temp=Ref%temp+y*Norm_t
     ifree=ifree+Nodes%n_nodes_t2
  End if
! Velocity
  If (Nodes%n_nodes_v2 .gt. 0) then
     Node_values(1:Nodes%n_nodes_v2)=X(ifree:ifree+Nodes%n_nodes_v2-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_v2, Nodes%i_nodes_v2, y)
     New_model%v_los=Ref%v_los+y*Norm_v
     ifree=ifree+Nodes%n_nodes_v2
  End if
! Microturbulence
  If (Nodes%n_nodes_mic2 .gt. 0) then
     Node_values(1:Nodes%n_nodes_mic2)=X(ifree:ifree+Nodes%n_nodes_mic2-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_mic2, Nodes%i_nodes_mic2, y)
     New_model%v_mic=Ref%v_mic+y*Norm_mic
     ifree=ifree+Nodes%n_nodes_mic2
  End if
! B long
  If (Nodes%n_nodes_blong2 .gt. 0) then
     Node_values(1:Nodes%n_nodes_blong2)=X(ifree:ifree+Nodes%n_nodes_blong2-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_blong2, Nodes%i_nodes_blong2, y)
     New_model%b_long=Ref%b_long+y*Norm_blong
     ifree=ifree+Nodes%n_nodes_blong2
  End if
! Bx
  If (Nodes%n_nodes_bx2 .gt. 0) then
     Node_values(1:Nodes%n_nodes_bx2)=X(ifree:ifree+Nodes%n_nodes_bx2-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_bx2, Nodes%i_nodes_bx2, y)
     New_model%b_x=Ref%b_x+y*Norm_bx
     ifree=ifree+Nodes%n_nodes_bx2
  End if
! By
  If (Nodes%n_nodes_by2 .gt. 0) then
     Node_values(1:Nodes%n_nodes_by2)=X(ifree:ifree+Nodes%n_nodes_by2-1)
     Call Expand_variable(npoints, New_model%ltau_500, Node_values, &
          Nodes%n_nodes_by2, Nodes%i_nodes_by2, y)
     New_model%b_y=Ref%b_y+y*Norm_by
     ifree=ifree+Nodes%n_nodes_by2
  End if
! Force Bx to be always positive by choosing that disambiguation direction
  Do ind=1, Params%n_points
     If (New_model%b_x(ind) .lt. 0) then
        New_model%b_x(ind)=-New_model%b_x(ind)
        New_model%b_y(ind)=-New_model%b_y(ind)
     End if
  End do
! Abundances
  If (Nodes%n_nodes_ab2 .gt. 0) then
     Do idx=1, Nodes%n_nodes_ab2
        New_model%abundance(Nodes%i_nodes_ab2(idx))=Ref%Abundance(Nodes%i_nodes_ab2(idx)) + &
             X(ifree)*Norm_ab
        ifree=ifree+1
     End do
  End if
!
  New_model_2comp%Comp2=New_model
  New_model_2comp%Comp2%ffactor=1.-New_model_2comp%Comp1%ffactor

!
  Call Hydrostatic(Params, New_model_2comp%Comp2, Error)
!
! Done!
!
  Return
!
End Subroutine Expand
! This routine takes the run of a given atmospheric variable and extracts
! the appropriate values at the nodes.
!
Subroutine Compress_variable(npoints, x, y, nnodes, inodes, yref, &
     Node_values)
  Implicit None
  Integer :: npoints, nnodes, ind, idepth
  Real :: num, den, ratio
  Real, dimension (npoints) :: x, y, yref
  Real, dimension (nnodes) :: xnodes, Node_values
  Integer, dimension (nnodes) :: inodes
!
  Node_values(1:nnodes) = y(inodes(1:nnodes)) - yref(inodes(1:nnodes))
!
  Return
End Subroutine Compress_variable
! This subroutine expands the run of a given atmospheric variable from
! its values at the nodes.
!
Subroutine Expand_variable(npoints, x, Node_values, nnodes, inodes, y)
  Use Bezier_math
  Implicit None
  Integer :: npoints, nnodes
  Real, dimension (npoints) :: x, y
  Real, dimension (nnodes) :: Node_values, xnodes
  Integer, dimension (nnodes) :: inodes
  Real :: slope
!
  xnodes=x(inodes)
  If (nnodes .eq. 0) then ! Nothing to do
     y(1:npoints)=0.
     Return
  Else if (nnodes .eq. 1) then ! Just one constant value
     y(1:npoints)=Node_values(1)
     Return
  Else if (nnodes .eq. 2) then ! A straight line
     slope=(Node_values(2)-Node_values(1))/(xnodes(2)-xnodes(1))
     y(1:npoints)=Node_values(1)+slope*(x(1:npoints)-xnodes(1))
  Else if (nnodes .ge. 3) then ! Splines Bezier
     Call bezier3(nnodes, xnodes, Node_values, npoints, x, y)
!     Call smoothed_lines(nnodes, npoints, xnodes, Node_values, x, y)
  End if
  Return
End Subroutine Expand_variable
! This subroutine performs a test trivial derivative compression, or, to be
! more precise, it does absolutely nothing.
!
Subroutine Compress_deriv(Params, Nodes, Derivatives, Dydx)
  Use Param_structure
  Use Nodes_information
  Implicit None
  Type (Parameters) :: Params
  Type (Nodes_info) :: Nodes
  Real, Dimension (Params%n_data, Params%n_variables, &
       Params%n_points) :: Derivatives
  Real, Dimension (Params%n_data, Params%n_free_parameters) :: Dydx
!
  Return
End Subroutine Compress_deriv
! This subroutine checks that the trial model does not exceed the boundaries set
! in the module Nodes_info. If any value for a particular variable is out of bounds,
! the run of that variable is averaged with the last accepted model (Guess_model).
! The process is repeated until all the variables are within range (or up to
! 10 times, whichever happens before). If there still remain variables out of
! range, they are forced to adopt the boundary values.
!
Subroutine Check_boundaries(Params, Nodes, Guess_model, Trial_model)
  Use Param_structure
  Use Model_structure
   Use Nodes_information
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Guess_model, Trial_model
  Type (Nodes_info) :: Nodes
  Integer :: ipoint, npoints
  Logical :: Out_of_range
  Real :: Maxi, Mini, Mx, Mn
!
  npoints=Params%n_points
!
! Check temperature
!
  Out_of_range=.FALSE.
  Do ipoint=1, npoints
     If (Trial_model%temp(ipoint) .lt. Min_t .or. &
          Trial_model%temp(ipoint) .gt. Max_t) Out_of_range=.TRUE.
  End do
  If (Out_of_range) then
     Where (Trial_model%temp .lt. Min_t) Trial_model%temp=Min_t
!     Mx=Maxval(Trial_model%temp)
!     Mn=Minval(Trial_model%temp)
!     If (Mx .lt. Min_t) then
!        Trial_model%temp=Min_t
!     Else if (Mn .gt. Max_t) then
!        Trial_model%temp=Max_t
!     Else
!        Maxi=Min(Max_t,Mx)
!        Mini=Max(Min_t,Mn)
!        Trial_model%temp=(Trial_model%temp - Mn)*(Maxi-Mini)/ &
!             (Mx-Mn) + Mini
!     End if
!     If (Params%printout .ge. 1) Print *,'Clipping temperature'
  End if
!
! Check velocity
!
  Out_of_range=.FALSE.
  Do ipoint=1, npoints
     If (Trial_model%v_los(ipoint) .lt. Min_v .or. &
          Trial_model%v_los(ipoint) .gt. Max_v) Out_of_range=.TRUE.
  End do
  If (Out_of_range) then
     Mx=Maxval(Trial_model%v_los)
     Mn=Minval(Trial_model%v_los)
     If (Mx .lt. Min_v) then
        Trial_model%v_los=Min_v
     Else if (Mn .gt. Max_v) then
        Trial_model%v_los=Max_v
     Else
        Maxi=Min(Max_v,Mx)
        Mini=Max(Min_v,Mn)
        Trial_model%v_los=(Trial_model%v_los - Mn)*(Maxi-Mini)/ &
             (Mx-Mn) + Mini
     End if
     If (Params%printout .ge. 1) Print *,'Clipping l.o.s. velocity'
  End if
!
! Check microturbulence
!
  Out_of_range=.FALSE.
  Do ipoint=1, npoints
     If (Trial_model%v_mic(ipoint) .lt. Min_mic .or. &
          Trial_model%v_mic(ipoint) .gt. Max_mic) Out_of_range=.TRUE.
  End do
  If (Out_of_range) then
     Mx=Maxval(Trial_model%v_mic)
     Mn=Minval(Trial_model%v_mic)
     If (Mx .lt. Min_mic) then
        Trial_model%v_mic=Min_mic
     Else if (Mn .gt. Max_mic) then
        Trial_model%v_mic=Max_mic
     Else
        Maxi=Min(Max_mic,Mx)
        Mini=Max(Min_mic,Mn)
        Trial_model%v_mic=(Trial_model%v_mic - Mn)*(Maxi-Mini)/ &
             (Mx-Mn) + Mini
     End if
     If (Params%printout .ge. 1) Print *,'Clipping microturbulence'
  End if
!
! Check b_long
!
  Out_of_range=.FALSE.
  Do ipoint=1, npoints
     If (Trial_model%b_long(ipoint) .lt. Min_blong .or. &
          Trial_model%b_long(ipoint) .gt. Max_blong) Out_of_range=.TRUE.
  End do
  If (Out_of_range) then
     Mx=Maxval(Trial_model%b_long)
     Mn=Minval(Trial_model%b_long)
     If (Mx .lt. Min_blong) then
        Trial_model%b_long=Min_blong
     Else if (Mn .gt. Max_blong) then
        Trial_model%b_long=Max_blong
     Else
        Maxi=Min(Max_blong,Mx)
        Mini=Max(Min_blong,Mn)
        Trial_model%b_long=(Trial_model%b_long - Mn)*(Maxi-Mini)/ &
             (Mx-Mn) + Mini
     End if
     If (Params%printout .ge. 1) Print *,'Clipping Bz'
  End if
!
! Check b_x
!
  Out_of_range=.FALSE.
  Do ipoint=1, npoints
     If (Trial_model%b_x(ipoint) .lt. Min_bx .or. &
          Trial_model%b_x(ipoint) .gt. Max_bx) Out_of_range=.TRUE.
  End do
  If (Out_of_range) then
     Mx=Maxval(Trial_model%b_x)
     Mn=Minval(Trial_model%b_x)
     If (Mx .lt. Min_bx) then
        Trial_model%b_x=Min_bx
     Else if (Mn .gt. Max_bx) then
        Trial_model%b_x=Max_bx
     Else
        Maxi=Min(Max_bx,Mx)
        Mini=Max(Min_bx,Mn)
        Trial_model%b_x=(Trial_model%b_x - Mn)*(Maxi-Mini)/ &
             (Mx-Mn) + Mini
     End if
     print *,'tri2=',trial_model%b_x
     If (Params%printout .ge. 1) Print *,'Clipping Bx'
  End if
!
! Check b_y
!
  Out_of_range=.FALSE.
  Do ipoint=1, npoints
     If (Trial_model%b_y(ipoint) .lt. Min_by .or. &
          Trial_model%b_y(ipoint) .gt. Max_by) Out_of_range=.TRUE.
  End do
  If (Out_of_range) then
     Mx=Maxval(Trial_model%b_y)
     Mn=Minval(Trial_model%b_y)
     If (Mx .lt. Min_by) then
        Trial_model%b_y=Min_by
     Else if (Mn .gt. Max_by) then
        Trial_model%b_y=Max_by
     Else
        Maxi=Min(Max_by,Mx)
        Mini=Max(Min_by,Mn)
        Trial_model%b_y=(Trial_model%b_y - Mn)*(Maxi-Mini)/ &
             (Mx-Mn) + Mini
     End if
     If (Params%printout .ge. 1) Print *,'Clipping By'
  End if
!
! Check macroturbulence
!
  If (Trial_model%v_mac .lt. Min_mac .or. &
       Trial_model%v_mac .gt. Max_mac) Trial_model%v_mac = &
       0.5*(Trial_model%v_mac + Guess_model%v_mac)
!
! Check stray
!
  If (Trial_model%stray .lt. Min_stray .or. &
       Trial_model%stray .gt. Max_stray) Trial_model%stray = &
       0.5*(Trial_model%stray + Guess_model%stray)
!
! Check expansion
!
  If (Trial_model%ffactor .lt. Min_ffactor) Trial_model%ffactor=Min_ffactor
  If (Trial_model%ffactor .gt. Max_ffactor) Trial_model%ffactor=Max_ffactor
!  If (Trial_model%ffactor .lt. Min_ffactor .or. &
!       Trial_model%ffactor .gt. Max_ffactor) Trial_model%ffactor = &
!       0.5*(Trial_model%ffactor + Guess_model%ffactor)
!
  If (Trial_model%v_mac .lt. Min_mac) Trial_model%v_mac=Min_mac
  If (Trial_model%v_mac .gt. Max_mac) Trial_model%v_mac=Max_mac
  If (Trial_model%stray .lt. Min_stray) Trial_model%stray=Min_stray
  If (Trial_model%stray .gt. Max_stray) Trial_model%stray=Max_stray
  Return
End Subroutine Check_boundaries
! This subroutine gives the errors in the atmospheric parameters from the
! vector of errors in the free parameters X_err and the current model
! Guess_model. The result is returned in Atmo_errors
!
Subroutine Expand_errors(Params, Nodes, Guess_model, X_err, Atmo_errors)
  Use Param_structure
  Use Model_structure
  Use Nodes_information
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Guess_model, Atmo_errors
  Type (Nodes_info) :: Nodes
  Real, Dimension (Params%n_free_parameters) :: X_err, Node_values
  Real, Dimension (Params%n_points) :: y, inc, azi
  Integer :: ifree, npoints, idepth
!
  npoints=Params%n_points
  Atmo_errors=Nodes%Reference_model%Comp1 ! Model_assign operation
  ifree=1
! Temperature
  If (Nodes%n_nodes_t .gt. 0) then
     Atmo_errors%temp=-1. ! Negative values signal no error bars at this point
     Atmo_errors%temp(Nodes%i_nodes_t)=X_err(ifree:ifree+Nodes%n_nodes_t-1) * &
          Norm_t
     ifree=ifree+Nodes%n_nodes_t
  End if
! Velocity
  If (Nodes%n_nodes_v .gt. 0) then
     Atmo_errors%v_los=-1. ! Negative values signal no error bars at this point
     Atmo_errors%v_los(Nodes%i_nodes_v)=X_err(ifree:ifree+Nodes%n_nodes_v-1) * &
          Norm_v
     ifree=ifree+Nodes%n_nodes_v
  End if
! Microturbulence
  If (Nodes%n_nodes_mic .gt. 0) then
     Atmo_errors%v_mic=-1. ! Negative values signal no error bars at this point
     Atmo_errors%v_mic(Nodes%i_nodes_mic)=X_err(ifree:ifree+Nodes%n_nodes_mic-1) &
          * Norm_mic
     ifree=ifree+Nodes%n_nodes_mic
  End if
! Magnetic strength
  If (Nodes%n_nodes_blong .gt. 0) then
     Atmo_errors%b_long=-1. ! Negative values signal no error bars at this point
     Atmo_errors%b_long(Nodes%i_nodes_blong)=X_err(ifree:ifree+Nodes%n_nodes_blong-1) &
          * Norm_blong
     ifree=ifree+Nodes%n_nodes_blong
  End if
! Magnetic inclination
  If (Nodes%n_nodes_bx .gt. 0) then
     Atmo_errors%b_x=-1. ! Negative values signal no error bars at this point
     Atmo_errors%b_x(Nodes%i_nodes_bx)=X_err(ifree:ifree+Nodes%n_nodes_bx-1) &
          * Norm_bx
     ifree=ifree+Nodes%n_nodes_bx
  End if
! Magnetic azimuth
  If (Nodes%n_nodes_by .gt. 0) then
     Atmo_errors%b_y=-1. ! Negative values signal no error bars at this point
     Atmo_errors%b_y(Nodes%i_nodes_by)=X_err(ifree:ifree+Nodes%n_nodes_by-1) &
          * Norm_by
     ifree=ifree+Nodes%n_nodes_by
  End if
! Macroturbulence
  If (Nodes%n_nodes_mac .eq. 1) then
     Atmo_errors%v_mac=X_err(ifree)*Norm_mac
     ifree=ifree+1
  Else
     Atmo_errors%v_mac=-1. ! Negative values signal no error bars
  End if
! Percentage of stray light
  If (Nodes%n_nodes_stray .eq. 1) then
     Atmo_errors%stray=X_err(ifree)*Norm_stray
     ifree=ifree+1
  Else
     Atmo_errors%stray=-1. ! Negative values signal no error bars
  End if
! Field ffactor
  If (Nodes%n_nodes_ffactor .eq. 1) then
     Atmo_errors%ffactor=X_err(ifree)*Norm_ffactor
     ifree=ifree+1
  Else
     Atmo_errors%ffactor=-1. ! Negative values signal no error bars
  End if
! Abundances
  If (Nodes%n_nodes_ab .gt. 0) then
     Atmo_errors%Abundance=-1. ! Negative values signal no error bars at this point
     Atmo_errors%Abundance(Nodes%i_nodes_ab)=X_err(ifree:ifree+Nodes%n_nodes_ab-1) &
          * Norm_ab
     ifree=ifree+Nodes%n_nodes_ab
  End if
! Other variables (not inverted)
  Atmo_errors%el_p=-1
  Atmo_errors%z_scale=-1
  Atmo_errors%gas_p=-1
  Atmo_errors%rho=-1
!
! Clip models
!
  Do idepth=1, npoints
     If (Atmo_errors%Temp(idepth) .gt. 1e4) Atmo_errors%Temp(idepth)=1e4
  End do
!
! Done!
!
  Return
!
End Subroutine Expand_errors

End Module Compex_module

