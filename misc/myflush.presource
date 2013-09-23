Subroutine MyFlush(iunit)
! This routine is a simple wrapper call to the non-standard F90 intrinsic 
! Flush() but written in such a way that it's easy to bypass it in case that
! the current compiler does not support the intrinsic Flush() call
Integer :: iunit
!
   Call Flush(iunit)
!
   Return
!
End Subroutine MyFlush
