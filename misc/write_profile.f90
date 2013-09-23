Subroutine Write_profile(Params, Region, Filename, iformat, Profile)
!
  Use Param_structure
  Use Line_data_structure
  Implicit None
  Type (Parameters) :: Params
  Type (Region_data), dimension (Params%n_regions) :: Region
  Character (Len=*) :: Filename
  Integer :: iformat, idata, ndata, iregion, iwave, nwaves, iunit, Get_lun
  Real, Dimension (Params%n_data) :: Profile
!
  If (iformat .eq. 1) then ! SIR format
     iunit=Get_lun()
     Open(Unit=iunit, File=Filename)
     idata=1
     Do iregion=1, Params%n_regions
        nwaves=Region(iregion)%nwavelengths
        Do iwave=0, nwaves-1
           Write (iunit, '(f10.3,4e15.6)') &
                Region(iregion)%First_wlength + &
                Region(iregion)%Wave_step*iwave, Profile(idata:idata+3)
           idata=idata+4
        End do
     End do
     Close (Unit=iunit)
  Else ! ASP format (not yet implemented)
     Print *,'ASP format not yet implemented'
     Stop
  End if
  Return
End Subroutine Write_profile
