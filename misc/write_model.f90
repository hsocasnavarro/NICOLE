! Subroutine to write onto disk the model atmosphere in the structure Atmo
!
Subroutine Write_model(Params, Filename, Atmo)
  Use Param_structure
  Use Model_structure
  Use File_operations
  Use Line_data_structure
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Atmo
  Integer :: File_unit, idepth, i_pos
  Real :: version
  Real, dimension (Params%n_points) :: tmp
  Character (len = 256) :: String
  Character (len = *) :: Filename
  Logical :: End
!
  Call Open_file(File_unit, Filename)
  Write (File_unit, *) 'Format version: 1.0'
  Write (File_unit, *) Atmo%v_mac, Atmo%stray, Atmo%ffactor
  Do idepth=1, Params%n_points
     Write (File_unit, '(f13.1,1x,f10.3,1x,6e15.6)') &
          Atmo%ltau_500(idepth), Atmo%temp(idepth), &
          Atmo%el_p(idepth), Atmo%v_mic(idepth), Atmo%B_long(idepth), &
          Atmo%v_los(idepth), Atmo%B_x(idepth), Atmo%B_y(idepth)
  End do
  Call Close_File(File_unit)
  Return
End Subroutine Write_model
