! This routine takes the run of a given atmospheric variable and extracts
! the appropriate values at the nodes.
!
Subroutine Compress_variable(npoints, x, y, nnodes, inodes, yref, &
     Node_values)
  Implicit None
  Integer :: npoints, nnodes, ind, idepth
  Real :: num, den, ratio
  Real, dimension (npoints) :: x, y, yref
  Real, dimension (nnodes) :: xnodes, Node_values
  Integer, dimension (nnodes) :: inodes
!
  Node_values(1:nnodes) = y(inodes(1:nnodes)) - yref(inodes(1:nnodes))
!
  Return
End Subroutine Compress_variable
