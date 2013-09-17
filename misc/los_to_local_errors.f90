! This routine will perform a coordinates conversion of the errors in the 
! magnetic field inclination and azimuth. It takes the heliocentric
! angle (theta), the local magnetic inclination (Gamma), the 
! line-of-sight magnetic inclination (Gamma_p) and azimuth
! (Chi_p), and their respective errors (dGamma_p, dChi_p), and returns the 
! errors in the local inclination (dGamma) and azimuth (dChi).
! All of these angles are real arrays of size npoints. This is obtained
! by derivating the formulae in the routine los_to_local.
!
! It is assumed that the heliocentriz azimuth is zero, which means that
! the observed point is on the E-W line that crosses the solar disk center.
! If that is not the case, the azimuth obtained is uncalibrated by an
! additive constant.
!
! navarro@hao.ucar.edu , Mar 2000
!
Subroutine los_to_local_errors(npoints, Theta, Gamma, Gamma_p, Chi_p, &
     dGamma_p, dChi_p, dGamma, dChi)
  Implicit None
  Integer :: npoints
  Real :: Theta
  Real, parameter :: Pi = 3.14159265389793D0, Deg_to_rad = Pi/180.
  Real, dimension (npoints) :: Gamma, Gamma_p, Chi_p, dGamma_p, dChi_p, &
       dGamma, dChi
  Integer :: ipoint
!
! Piece of cake
!
  dChi=dChi_p
!
!     There is a special case when gamma = 0 (it results in a 0/0 
!                     avoidable uncertainty)
  Do ipoint=1, npoints
     If (Gamma(ipoint) .lt. 1.E-3) Gamma(ipoint)=1.E-3
  End do
  dGamma=(Cos(Gamma_p*Deg_to_rad)*Cos(Chi_p*Deg_to_rad)* &
       Sin(Theta*Deg_to_rad)*dGamma_p - &
       Sin(Gamma_p)*Sin(Chi_p)*Sin(Theta)*dChi_p - &
       Sin(Gamma_p)*Cos(Theta)*dGamma_p)/Sin(Gamma)
!
! Done!
!
  dGamma=Abs(dGamma)
  Return
End Subroutine los_to_local_errors
