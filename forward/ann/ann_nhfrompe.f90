! Use a Neural network to obtain Pg from (T,Pe)
!
Subroutine ann_nhfrompe(T, Pe, Metalicity, nH)
  Use ann_nhfrompe_data
  Use Debug_module
  Real :: T, T2, Pe, nH, Metalicity, Met2
  Real :: LogPe, LognH

  LogPe=Log10(Pe)
  If (LogPe .gt. 5) then
     Debug_warningflags(flag_computepg)=1
     Call Debug_Log('Log10 (Pe) .gt. 5. Clipping it',2)
     LogPe=5
  End if
  If (LogPe .lt. -3) then
     Debug_warningflags(flag_computepg)=1
     Call Debug_Log('Log10 (Pe) .lt. -3. Clipping it',2)
     LogPe=-3
  End if
  T2=T
  If (T .gt. 1e4) then
     Debug_warningflags(flag_computepg)=1
     Call Debug_Log('T .gt. 1E4. Clipping it',2)
     T2=1e4
  End if
  If (T .lt. 1500) then
     Debug_warningflags(flag_computepg)=1
     Call Debug_Log('T .lt. 1500. Clipping it',2)
     T2=1500
  End if
  Met2=Metalicity
  If (Met2 .gt. .5) then
     Debug_warningflags(flag_computepg)=1
     Call Debug_Log('Metalicity .gt. 0.5. Clipping it',2)
     Met2=.5
  End if
  If (Met2 .lt. -1.5) then
     Debug_warningflags(flag_computepg)=1
     Call Debug_Log('Metalicity .lt. -1.5. Clipping it',2)
     Met2=-1.5
  End if
  inputs(1)=(T2-xmean(1))/xnorm(1)
  inputs(2)=(LogPe-xmean(2))/xnorm(2)
  inputs(3)=(Met2-xmean(3))/xnorm(3)

  Call ANN_Forward(W, Beta, Nonlin, inputs, outputs, nlayers, &
       nmaxperlayer, nperlayer, ninputs, noutputs, y)

  LognH=outputs(1)*ynorm(1)+ymean(1)
  nH=10.**(LognH)

End Subroutine ann_nhfrompe
