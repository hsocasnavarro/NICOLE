! Printout information on the current iteration
!
Subroutine Printout(n_iter, Lambda, Chisq, NWChisq, Regul, Params, Line, Region, &
     Guess_model, Model_errors, Syn_profile, Converged)
  Use Param_structure
  Use Model_structure
  Use Line_data_structure
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Guess_model, Model_errors
  Type (Line_data), dimension(Params%n_lines) :: Line
  Type (Region_data), dimension(Params%n_regions) :: Region
  Real, Dimension (Params%n_data) :: Syn_profile
  Integer :: n_iter
  Real :: Lambda, Chisq, NWChisq, Regul
  Logical :: Converged
!
  If (Params%Printout .le. 1) Return 
  If (Converged) then
     Write (*,'("iter= ",i0," Lambda=",e10.3," Regul=",e10.3," Chisq=",e10.3)') n_iter,Lambda,Regul,Chisq-Regul
!     Print *,'iter=',n_iter,' Lambda=',Lambda,' Chisq=',Chisq
     If (Params%Printout .ge. 2) then 
        Call Write_model(Params, 'tmp.mod', Guess_model)
        If (n_iter .gt. 1) Call Write_model(Params, 'tmp.err', Model_errors)
        Call Write_profile(Params, Region, 'tmp.pro', 1, Syn_profile)
     End if
  Else
     Write (*,'("REJECTED: --- iter= ",i0," Lambda=",e10.3," Regul=",e10.3," Chisq=",e10.3)') n_iter,Lambda,Regul,Chisq-Regul
!     Print *,'REJECTED: ---- iter=',n_iter,' Lambda=',Lambda,' Chisq=',Chisq,' Regul=',Regul
     If (Params%Printout .ge. 3) then 
        Call Write_model(Params, 'tmp_failed.mod', Guess_model)
        If (n_iter .gt. 1) Call Write_model(Params, 'tmp_failed.err', Model_errors)
        Call Write_profile(Params, Region, 'tmp_failed.pro', 1, Syn_profile)
     End if
  End if
  Return
!
End Subroutine Printout
