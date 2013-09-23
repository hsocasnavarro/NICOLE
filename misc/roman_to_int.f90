! Take a string with a Roman number and convert to integer
! If the number is already in Arabic form then just return it as integer
Subroutine Roman_to_int(String, IntResult, Error)
  Implicit None
  Character (len=*) :: String
  Character (len=6) :: Buffer
  Integer :: IntResult
  Integer :: i
  Logical :: Error
  Character (len=6), Dimension(30) :: Romans
  Data Romans/'i','ii','iii','iv','v','vi','vii','viii','ix','x','xi',&
       'xii','xiii','xiv','xv','xvi','xvii','xviii','xix','xx','xxi',&
       'xxii','xxiii','xxiv','xxv','xxvi','xxvii','xxviii','xxix','xxx'/

  Call Tolower(String)
  String=AdjustL(String)
  If (SCAN(String,'ivx') .eq. 0) then ! No Roman digit
     If (VERIFY(String,' 0123456789') .ne. 0) then ! Strange digit!
        IntResult=-1
        Error=.True.
        Return
     End if
     Read (String,*) IntResult
  Else ! Contains Roman digits
     If (VERIFY(String,' ivx') .ne. 0) then ! Strange digit!
        IntResult=-2
        Error=.True.
        Return
     End if
     IntResult=-1
     Do i=1, 30
        Buffer=String(1:6)
        If (Buffer .eq. Romans(i)) IntResult=i
     End do
     If (IntResult .eq. -1) then
        Error=.True.
        Return
     End if
  End if
  Return

End Subroutine Roman_to_int
