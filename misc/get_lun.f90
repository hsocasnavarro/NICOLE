! This subroutine returns the first free available file unit.
! If no available units are found, it aborts with an error message.
!
Integer Function Get_lun()
!
  Logical :: used
  Integer :: i
  Parameter (UNIT_MIN = 20)
  Parameter (UNIT_MAX = 119)
!
  get_lun=0
  used=.TRUE.
  i=UNIT_MIN
  Do while (i .le. UNIT_MAX .and. used)
     Inquire (unit=i, opened=used)
     If (.not.used) then
        get_lun=i
     End if
     i=i+1
  End do
  If (get_lun .eq. 0) then 
     Print *,'No available file units found!!! Aborting ...'
     Stop
  End if
!
End Function Get_lun
