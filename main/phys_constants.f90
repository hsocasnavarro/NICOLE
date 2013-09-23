Module Phys_constants
! Physical and Mathematical constants
! Solar constants
! c.g.s. units
  Real :: Gravity=2.7414e+4 ! cm/s^2 ! Default. Solar surface gravity
  Real, Parameter :: Solar_Gravity = 2.7414e+4 ! Surface gravity (cm/s^2)
! Physical constants
! c.g.s. units
  Real, Parameter :: cc = 2.99792458D10, hh = 6.62618D-27, ee = 1.602189D-12 
  Real, Parameter :: mass_el = 9.10953D-28, mass_pr = 1.672648D-24
  Real, Parameter :: bk = 1.38066D-16, Avog=6.022045D23
! Conversion factor: eV to erg (g cm^2 / s^2). 
! Note that 1 V = 10^7 (g cm^2)/(C s^2)
  Real, Parameter :: eV_to_cgs = ee
! Mathematical constants
  Real, Parameter :: Pi = 3.14159265389793D0
  Real, Parameter :: SqrtPi = 1.77245390415191650391 ! Sqrt(Pi)
  Real, Parameter :: Piis = 0.56418956659890107D0 ! Inverse of the Sqrt of Pi
End Module Phys_constants
