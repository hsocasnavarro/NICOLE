Module Nodes_information
Use Model_structure
!
! temp is in K, v_los and v_mic are in cm/s, el_p and gas_p are in dyn/cm^2, 
! rho is in g/cm^3, b_long, b_x and b_y are in gauss
! v_mac is in cm/s, stray and ffactor are dimensionless (0 < stray < 1).
  real, parameter :: Norm_t=500., Norm_v=1.e5, Norm_mic=1.e5, &
       Norm_blong=500., Norm_bx=1000., Norm_by=1000., Norm_mac=1.e5, &
       Norm_stray=0.10, Norm_ffactor=0.10, Norm_ab=1.
!  real, parameter :: Max_t=30000., Max_v=2.e6, Max_mic=2.e6, &
!       Max_b=10000., Max_inc=180., Max_azi=360., Max_mac=5.e5, &
!       Max_stray=1., Max_ffactor=1.
!  real, parameter :: Min_t=3000., Min_v=-2.e6, Min_mic=0., &
!       Min_b=0., Min_inc=0., Min_azi=0., Min_mac=0., &
!       Min_stray=0., Min_ffactor=1.
  real, parameter :: Max_t=30000., Max_v=2.e6, Max_mic=2.e6, &
       Max_blong=1e4, Max_bx=1e4, Max_by=1e4, Max_mac=5.e5, &
       Max_stray=1., Max_ffactor=1., Max_ab=12.0
  real, parameter :: Min_t=2500., Min_v=-2.e6, Min_mic=0., &
       Min_blong=-1e4, Min_bx=-1e4, Min_by=-1e4, Min_mac=0., &
       Min_stray=0., Min_ffactor=0., Min_ab=3.0
  real, dimension(100) :: X_max
  real, dimension(12,100) :: UserNodeLocations
  Type Nodes_info
     integer :: n_nodes_t, n_nodes_v, n_nodes_mic, & ! Number of nodes 
          n_nodes_blong, n_nodes_bx, n_nodes_by, & ! for each variable.
          n_nodes_mac, n_nodes_stray, n_nodes_ffactor, n_nodes_ab
     integer, dimension(:), pointer :: i_nodes_t, i_nodes_v, i_nodes_mic, &
          i_nodes_blong, i_nodes_bx, i_nodes_by, i_nodes_ab ! Where are the nodes.
     integer :: n_nodes_t2, n_nodes_v2, n_nodes_mic2, & ! Number of nodes 
          n_nodes_blong2, n_nodes_bx2, n_nodes_by2, & ! for each variable.
          n_nodes_ffactor2, n_nodes_ab2
     integer, dimension(:), pointer :: i_nodes_t2, i_nodes_v2, i_nodes_mic2, &
          i_nodes_blong2, i_nodes_bx2, i_nodes_by2, i_nodes_ab2 ! Where are the nodes.
     Type (Model_2comp) :: Reference_model
  End Type Nodes_info
!
End Module Nodes_information
