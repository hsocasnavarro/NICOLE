! This routine solves the linear system A X = B using the SVD technique
! The diagonal elements of the covariance matrix are returned in B
!
Subroutine SVD_solve(NP, SVD_threshold, A, B, X, Zeroed)
  USE nr
  USE Debug_module
  Implicit None
  Integer :: NP, i_param, j_param
  Real, Dimension (NP) :: X, B
  Real, Dimension (NP,NP) :: A
  Real(kind=8), Dimension (NP) :: Xdp, Bdp, W
  Real(kind=8), Dimension (NP,NP) :: Adp, V
  Logical, Dimension (NP) :: Zeroed
  Real :: Wmax, SVD_threshold
!
  Adp=A
  Bdp=B
  Call SVDcmp(Adp, W, V)
  If (Debug_errorflags(flag_svdcmp) .ge. 1) then
     X(:)=0.
     Call Debug_Log('Error in svdcmp... returning from SVD_solve with X=0',2)
     Return
  End if
  Wmax=Maxval(W)
  Zeroed(1:NP)=.FALSE.
  Do i_param=1, NP
     If (W(i_param) .lt. SVD_threshold*Wmax) then
        W(i_param)=0.
        Zeroed(i_param)=.TRUE.
     End if
  End do
  Call SVBKSB(Adp, W, V, Bdp, Xdp)
  B=Bdp
  X=Xdp
  A=Adp
  B(1:NP)=0.
  Do i_param=1, NP
     Do j_param=1, NP
        If (.not. Zeroed(j_param)) Then
           If (Abs(V(i_param,j_param)) .gt. .1) then
!          If there is a significant contribution from this w(j_param) to
!          the B(i_param), then add it.
!        Numerical recipes formula
!              B(i_param)=B(i_param)+V(i_param,j_param)*V(i_param,j_param)/ &
!                   (W(j_param)*W(j_param))
!        Diagonal of A^-1 calculated from A=UWV^t => A^-1=V(1/w)U^t
              B(i_param)=B(i_param)+&
                   V(i_param, j_param)*A(i_param, j_param)/W(j_param)
           End if
        Else
!          If there is a significant contribution from a zeroed eigenvalue
!          to the B(i_param), then make it infinity.
           If (Abs(V(i_param,j_param)) .gt. .1) &
                B(i_param)=100.
        End if
     End do
  End do
!  B(1:NP)=Sqrt(B(1:NP)) 
  Return                
!                       
End Subroutine SVD_solve
