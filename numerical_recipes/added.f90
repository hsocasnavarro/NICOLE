Module nrutil2
  USE nrtype
! Some additional definitions for NICOLE
  INTERFACE imaxloc_dp
     MODULE PROCEDURE imaxloc_dp
  END INTERFACE imaxloc_dp
  INTERFACE swap_dp
     MODULE PROCEDURE swap_dp
  END INTERFACE swap_dp
CONTAINS
  FUNCTION imaxloc_dp(arr)
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B) :: imaxloc_dp
    INTEGER(I4B), DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc_dp=imax(1)
  END FUNCTION imaxloc_dp
  SUBROUTINE swap_dp(a,b)
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(DP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_dp
End Module nrutil2
