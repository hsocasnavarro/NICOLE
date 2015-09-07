Module NICOLE_inp
!
  Type NICOLE_input
     character :: mode
     character (len = 256) :: Obs_profile_file, Syn_profile_file, & 
          Model_in_file, Model_in_file_2, Model_out_file, Model_out_file_2, Stray_profile_file
     character (len = 4) :: input_dens
     character (len = 4) :: hscale
     logical :: set_hydro, set_nH
     integer :: depcoef_mode
     Logical :: write_depcoef
     integer :: formal_solution, printout, maxinv, speed, &
          always_compute_deriv, reference_cont
     real :: heliocentric, helioazim, noise, acceptchisq, gravity, &
          update_opac
  End Type NICOLE_input
!
End Module NICOLE_inp
