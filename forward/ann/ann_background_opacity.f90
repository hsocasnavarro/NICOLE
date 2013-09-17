! Use a Neural network to obtain background opacity from (T,Pe,Metalicity,lambda)
!
Subroutine ann_background_opacity(T, Pe, Pg, Metalicity, Lambda, Opac)
  Use ann_background_opacity_data
  Use profiling
  Use Debug_module
  Use background_opacity_module, only: background_opacity
  Implicit None
  Real :: T, T2, Pe, Pg, Metalicity, Met2, Lambda, Opac
  Real :: LogPe, LogPg
  Real :: a,b=0.,c=0.,d,e=0.,scat ! debug
  Logical :: UseOpacANN=.False.

  Call Time_routine('ann_background_opacity',.True.)

  If (UseOpacANN) then
     LogPe=Log10(Pe)
     If (LogPe .gt. 4) then
        Debug_warningflags(flag_computeopac)=1
        Call Debug_Log('Log10 (Pe) .gt. 4. Clipping it',2)
        LogPe=4
     End if
     If (LogPe .lt. -3) then
        Debug_warningflags(flag_computeopac)=1
        Call Debug_Log('Log10 (Pe) .lt. -3. Clipping it',2)
        LogPe=-3
     End if
     T2=T
     If (T .gt. 1e4) then
        Debug_warningflags(flag_computeopac)=1
        Call Debug_Log('T .gt. 1E4. Clipping it',2)
        T2=1e4
     End if
     If (T .lt. 1500) then
        Debug_warningflags(flag_computeopac)=1
        Call Debug_Log('T .lt. 1500. Clipping it',2)
        T2=1500
     End if
     Met2=Metalicity
     If (Met2 .gt. .5) then
        Debug_warningflags(flag_computeopac)=1
        Call Debug_Log('Metalicity .gt. 0.5. Clipping it',2)
        Met2=.5
     End if
     If (Met2 .lt. -1.5) then
        Debug_warningflags(flag_computeopac)=1
        Call Debug_Log('Metalicity .lt. -1.5. Clipping it',2)
        Met2=-1.5
     End if
     If (Lambda .gt. 16500) then
        Debug_warningflags(flag_computeopac)=1
        Call Debug_Log('Lambda .gt. 16500. Clipping it',2)
        Lambda=16500.
     End if
     If (Lambda .lt. 3800) then
        Debug_warningflags(flag_computeopac)=1
        Call Debug_Log('Lambda .lt. 3800. Clipping it',2)
        Lambda=3800
     End if
     !
     inputs(1)=(T2-xmean(1))/xnorm(1)
     inputs(2)=(LogPe-xmean(2))/xnorm(2)
     inputs(3)=(Met2-xmean(3))/xnorm(3)
     inputs(4)=(Lambda-xmean(4))/xnorm(4)

     Call ANN_Forward(W, Beta, Nonlin, inputs, outputs, nlayers, &
          nmaxperlayer, nperlayer, ninputs, noutputs, y)

     opac=outputs(1)*ynorm(1)+ymean(1)
     opac=10.**(opac)
  End if

  call ann_nhfrompe(T,Pe,Metalicity, a)
  a=a* 1.38066D-16*T ! pH
  call ann_nh2frompe(T,Pe,Metalicity, d)
  d=d* 1.38066D-16*T ! pH2
  opac=background_opacity(T,Pe,Pg,a,b,c,d,e,lambda,scat)


  Call Time_routine('ann_background_opacity',.False.)
  Return

End Subroutine ann_background_opacity
