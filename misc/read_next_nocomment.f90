Subroutine Read_next_nocomment(File_unit, String, End)
! This subroutine reads the next non-comment line from the file unit
! File_unit. String must be 256 characters long. End is a logical variable
! that is set to true if the end of the file has been reached. 
!
  Implicit None
  Integer :: File_unit
  Character (len=256) :: String
  Logical :: End
!
  End=.FALSE.
  String=''
  Read(File_unit,'(A)', End=100) String
  String=Trim(AdjustL(String))
  Do while (String(1:1) .eq. '!')
     Read(File_unit,'(A)', End=100) String
     String=Trim(AdjustL(String))
  End do
  Return
100 End=.TRUE.
  Return
End Subroutine Read_next_nocomment
!
Subroutine Read_char(Line, String, Ch, i_pos)
!
  Character (Len = *) :: Line, String
  Character :: Ch
  Integer :: i_pos, j_pos, k_pos, l_pos
!
  Call Remove_comment(Line, '!')
  i_pos=Index(Line, String)
  If (i_pos .ne. 0) then
     j_pos=Index(Line(i_pos:256),':')+i_pos-1
     k_pos=Index(Line(j_pos:256),'!')+j_pos-1
     if (k_pos .eq. j_pos-1) k_pos=256
     Do l_pos=j_pos, k_pos
        If (Line(l_pos:l_pos).ge.'A'.and.Line(l_pos:l_pos).le.'z') then
           Ch=Line(l_pos:l_pos)
           Return
        End if
     End do
     Print *,'Error! Character not found after : and before ! in'
     Print *,Line
     Stop
  End if
End Subroutine Read_char
!
Subroutine Read_string(Line, String, String_out, i_pos)
!
  Implicit None
  Character (Len = *) :: Line, String, String_out
  Integer :: i_pos, j_pos, k_pos, l_pos, m_pos
!
  Call Remove_comment(Line, '!')
  i_pos=Index(Line, String)
  If (i_pos .ne. 0) then
     j_pos=Index(Line(i_pos:256),':')+i_pos-1
     k_pos=Index(Line(j_pos:256),'!')+j_pos-1
     if (k_pos .eq. j_pos-1) k_pos=256
     l_pos=j_pos+1
     Do While (l_pos.lt.k_pos-1.and.(Line(l_pos:l_pos).lt.'.' &
          .or.Line(l_pos:l_pos).gt.'z')) 
        l_pos=l_pos+1
     End do
     m_pos=k_pos-1
     Do While (m_pos.gt.j_pos+1.and.(Line(m_pos:m_pos).lt.'.' &
          .or.Line(m_pos:m_pos).gt.'z')) 
        m_pos=m_pos-1
     End do
     String_out=Line(l_pos:m_pos)
  End if
End Subroutine Read_string
!
Subroutine Read_real(Line, String, Value, i_pos)
!
  Implicit None
  Character (Len = *) :: Line, String
  Character (Len = 256) :: String_out
  Integer :: i_pos, j_pos, k_pos, l_pos, m_pos
  Real :: Value
!
  Call Remove_comment(Line, '!')
  i_pos=Index(Line, String)
  If (i_pos .ne. 0) then
     j_pos=Index(Line(i_pos:256),':')+i_pos-1
     k_pos=Index(Line(j_pos:256),'!')+j_pos-1
     if (k_pos .eq. j_pos-1) k_pos=256
     l_pos=j_pos+1
     Do While (l_pos.lt.k_pos-1.and.(Line(l_pos:l_pos).lt.'-' &
          .or.Line(l_pos:l_pos).gt.'z')) 
        l_pos=l_pos+1
     End do
     m_pos=k_pos-1
     Do While (m_pos.gt.j_pos+1.and.(Line(m_pos:m_pos).lt.'-' &
          .or.Line(m_pos:m_pos).gt.'z')) 
        m_pos=m_pos-1
     End do
     String_out=Line(l_pos:m_pos)
     Read (String_out,*) Value
  End if
End Subroutine Read_real
!
Subroutine Read_int(Line, String, Value, i_pos)
!
  Implicit None
  Character (Len = *) :: Line, String
  Character (Len = 256) :: String_out
  Integer :: i_pos, j_pos, k_pos, l_pos, m_pos
  Integer :: Value
!
  Call Remove_comment(Line, '!')
  i_pos=Index(Line, String)
  If (i_pos .ne. 0) then
     j_pos=Index(Line(i_pos:256),':')+i_pos-1
     k_pos=Index(Line(j_pos:256),'!')+j_pos-1
     if (k_pos .eq. j_pos-1) k_pos=256
     l_pos=j_pos+1
     Do While (l_pos.lt.k_pos-1.and.(Line(l_pos:l_pos).lt.'-' &
          .or.Line(l_pos:l_pos).gt.'z')) 
        l_pos=l_pos+1
     End do
     m_pos=k_pos-1
     Do While (m_pos.gt.j_pos+1.and.(Line(m_pos:m_pos).lt.'-' &
          .or.Line(m_pos:m_pos).gt.'z')) 
        m_pos=m_pos-1
     End do
     String_out=Line(l_pos:m_pos)
     Read (String_out,*) Value
  End if
End Subroutine Read_int
!
Subroutine Remove_comment(Line, Comment_char)
!
  Implicit None
  Character (Len = *) :: Line
  Character :: Comment_char
  Integer :: length, i
!
  length=Len(Line)
  i=Index(Line, Comment_char)
  If (i .eq. 0) Return
  Line(i:length)=''
  Return
End Subroutine Remove_comment
