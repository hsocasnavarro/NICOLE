Subroutine Read_profile(Params, Filename, Profile)
!
  Use Param_structure
  Implicit None
  Type (Parameters) :: Params
  Character (Len=*) :: Filename
  Integer :: idata, iline, iunit, Get_lun
  Integer :: status, nwlengths
  Real :: Line_number, wave
  Real, Dimension(Params%n_data) :: Profile
  Logical :: Exists
!
  Inquire (File=Filename, Exist = Exists)
  If (.NOT. Exists) then
     Print *,'Profile file ',Filename,' does not exists. Aborting ...'
     Stop
  End if
  iunit=Get_lun()
  Open(Unit=iunit, File=Filename)
  idata=1
  nwlengths=Params%n_data/4
  Do iline=1, nwlengths
     Read (iunit, *) wave, Profile(idata:idata+3)
     idata=idata+4
  End do
  Close (Unit=iunit)
  Return
End Subroutine Read_profile
