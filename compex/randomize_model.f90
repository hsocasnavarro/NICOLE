Subroutine Randomize_Model(Params, Nodes, ModelIn, ModelOut)
  Use Param_structure
  Use Nodes_information
  Use Model_structure
  Implicit None
  Type (Nodes_info) :: Nodes
  Type (Parameters) :: Params
  Type (Model) :: ModelIn, ModelOut
  Integer :: i, j
  Real :: kk
  Logical, save :: firsttime = .TRUE.
!
  If(firsttime) Then
     Call RANDOM_SEED
     firsttime = .FALSE.
  End if
!
  ModelOut=ModelIn

  If (Nodes%n_nodes_t .gt. 0) then
     Call Random_number(kk)
     ModelOut%temp=ModelOut%temp + kk*100.
  End if

  If (Nodes%n_nodes_v .gt. 0) then
     Call Random_number(kk)
     ModelOut%v_los=ModelOut%v_los + (kk-.5)*3.e5
  End if

  If (Nodes%n_nodes_mic .gt. 0) then
     Call Random_number(kk)
     ModelOut%v_mic=ModelOut%v_mic + kk*3.e5
  End if

  If (Nodes%n_nodes_blong .gt. 0) then
     Call Random_number(kk)
     ModelOut%b_long=ModelOut%b_long + (kk-.5)*1000.
  End if

  If (Nodes%n_nodes_bx .gt. 0) then
     Call Random_number(kk)
     ModelOut%b_x=ModelOut%b_x + (kk-.5)*1000.
  End if

  If (Nodes%n_nodes_by .gt. 0) then
     Call Random_number(kk)
     ModelOut%b_y=ModelOut%b_y + (kk-.5)*1000.
  End if

  If (Nodes%n_nodes_stray .gt. 0) then
     Call Random_number(kk)
     ModelOut%stray=kk*.7
  End if

  If (Nodes%n_nodes_ffactor .gt. 0) then
     Call Random_number(kk)
     ModelOut%ffactor=kk*.5
  End if

  If (Nodes%n_nodes_mac .gt. 0) then
     Call Random_number(kk)
     ModelOut%v_mac=kk*2.e5
  End if

End Subroutine Randomize_Model


