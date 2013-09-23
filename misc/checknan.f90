! This function checks whether the argument is a valid floating point number

Logical Function CheckNaN(Number)
!
Implicit None
Real :: Number
!
CheckNaN=.False.
If (Number .gt. Huge(Number)) then 
!   Print *,'Larger than Huge'
   CheckNan=.True.
End if
If (Number .lt. -Huge(Number)) then
!   Print *,'Smaller than -Huge'
   CheckNan=.True.
End if
If (Number .ne. Number) then
!   Print *,'Not equal to itself'
   CheckNan=.True.
End if
If ( (Number .gt. 0.0) .EQV. (Number .le. 0.0) ) then
!   Print *,'Does not pass the EQV test'
   CheckNan=.True.
End if
If (Abs(Number) .lt. 0.) then
!   Print *,'Negative absolute value'
   CheckNan=.True.
End if
Return
End Function CheckNaN
