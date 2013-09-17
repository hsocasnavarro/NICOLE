Subroutine Adjust_wave_grid(Params, Chisq, Adjusted)
  Use Param_structure
  Implicit None
  Type (Parameters) :: Params
  Integer :: OldSkip
  Real :: Chisq
  Logical :: Adjusted
!
  Adjusted=.FALSE.
  Return
!
  OldSkip=Params%Skip_lambda
  Params%Skip_lambda=Min(10,OldSkip)
  If (Chisq .lt. 100) Params%Skip_lambda=Min(2,OldSkip)
  If (Chisq .lt. 30) Params%Skip_lambda=Min(1,OldSkip)
  If (OldSkip .NE. Params%Skip_lambda) Adjusted=.TRUE.
  If (Adjusted .and. Params%Printout .gt. 0) &
       Print *,'Wavelength grid adjusted. Computing every ',Params%Skip_lambda,&
       ' points'
!
  Return
End Subroutine Adjust_wave_grid
