! This routine will perform a coordinates conversion. It takes the 
! heliocentric angle (theta), the line-of-sight magnetic inclination 
! (Gamma_p) and azimuth (Chi_p) and returns the local inclination (Gamma) 
! and azimuth (Chi). All of these angles are real arrays of size npoints.
!
! It is assumed that the heliocentriz azimuth is zero, which means that
! the observed point is on the E-W line that crosses the solar disk center.
! If that is not the case, the azimuth obtained is uncalibrated by an
! additive constant.
!
! navarro@hao.ucar.edu , Mar 2000
!
Subroutine los_to_local(npoints, Theta, Gamma_p, Chi_p, Gamma, Chi)
  Implicit None
  Integer :: npoints
  Real :: Theta
  Real, parameter :: Pi = 3.14159265389793D0, Deg_to_rad = Pi/180.
  Real, dimension (npoints) :: Gamma, Chi, Gamma_p, Chi_p
!
! Do NOT convert!
!
  Chi=Chi_p
  Gamma=Gamma_p
  Return
!
!
!
!
!
! Piece of cake
!
  Chi=Chi_p
  Gamma=Acos(Cos(Gamma_p*Deg_to_rad)*Cos(Theta*Deg_to_rad) + &
       Sin(Gamma_p*Deg_to_rad)*Cos(Chi_p*Deg_to_rad)*Sin(Theta*Deg_to_rad)) /&
       Deg_to_rad
!
! Done!
!
  Return
End Subroutine los_to_local
