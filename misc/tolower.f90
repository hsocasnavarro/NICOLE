! Subroutine to convert a string to lower case.
!
Subroutine Tolower(String)
!
  Implicit None
  Character (len=*) :: String
  Integer :: i_pos, j_pos, l_a, c_a, l_z, c_z
!
  l_a=ichar('a')
  l_z=ichar('z')
  c_a=ichar('A')
  c_z=ichar('Z')
  j_pos=Index(String,':')
  If (j_pos .eq. 0) j_pos=Len(String)
  Do i_pos=1,j_pos
     If (String(i_pos:i_pos) .ge. 'A' .and. String(i_pos:i_pos) .le. 'Z') &
          String(i_pos:i_pos)=char(ichar(String(i_pos:i_pos))+l_a-c_a)
  End do
  Return
End Subroutine Tolower
