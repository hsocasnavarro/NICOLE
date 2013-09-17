! This routine will perform a coordinates conversion. It takes the 
! heliocentric angle (theta), the local magnetic inclination (Gamma) and
! azimuth (Chi) and returns the line-of-sight inclination (Gamma_p) and
! azimuth (Chi_p). All of these angles are real arrays of size npoints.
!
! It is assumed that the heliocentriz azimuth is zero, which means that
! the observed point is on the E-W line that crosses the solar disk center.
! If that is not the case, the azimuth obtained is uncalibrated by an
! additive constant.
!
! navarro@hao.ucar.edu , Mar 2000
!
Subroutine local_to_los(npoints, Theta, Gamma, Chi, Gamma_p, Chi_p)
  Implicit None
  Integer :: npoints
  Real :: Theta
  Real, parameter :: Pi = 3.14159265389793D0, Deg_to_rad = Pi/180.
  Real, dimension (npoints) :: Gamma, Chi, Gamma_p, Chi_p
!
! Do NOT convert!
!
  Chi_p=Chi
  Gamma_p=Gamma
  Return
!
! Piece of cake
!
  Chi_p=Chi
  Gamma_p=Acos(Cos(Gamma*Deg_to_rad)*Cos(Theta*Deg_to_rad) - &
       Sin(Gamma*Deg_to_rad)*Cos(Chi*Deg_to_rad)*Sin(Theta*Deg_to_rad)) / &
       Deg_to_rad
!
! Done!
!
  Return
End Subroutine local_to_los
