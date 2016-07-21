Module Line_data_structure
!
  Type Line_data
     Integer :: Mult_low, Mult_up
     Character :: Desig_low, Desig_up
     Character (len = 2) :: Elem
     Real :: J_low, J_up
     Real :: Wlength ! in Angstroms
     Real :: Width ! in Angstroms. Roughly (estimate by excess to be safe)
     Real :: VDW_enh, Energy_low, loggf ! Energy in eV
     Integer :: Atomic_number, Ion_stage ! Z and ionization stage
     Integer :: DepCoefMode=1 ! Switch for behavior of depcoef.dat
                              ! 1->(Default) Data read are departure coeffs
                              ! 2->Data read are actual level populations
     Real :: Atomic_weight
     Integer :: collisions ! Treatment for collisions with neutral H:
                           !    1->Unsold formula
                           !    2->Use Barklem coefficients
                           !    3->Give Gamma coefs explicitly
     Real :: Bark_alpha, Bark_sigma ! Barklem coefficients
     Real :: Gamma_rad, Gamma_Strk_12, Gamma_vdW_16 ! Gamma coefs explicitly
            ! (Gammas are in 1E8 rad/s, Gamma_vdW is per 1E16 H atoms per cm3,
            ! Gamma_Strk is per 1E12 electrons per cm3; G_Strk and G_vdW are
            ! quoted at T=1E4 K with temperature dependence of T^0.17 and
            ! T^0.38, respectively. A value of -1 means compute using Unsold)
     Logical :: Hyperfine 
     Real :: HF_Alow, HF_Blow, HF_Aup, HF_Bup, SpinI ! For hyperfine structure
            !                                            A and B are in cm-1
     Real :: extra_vmic ! Extra line broadening, behaves like microturb
                        ! useful e.g. to simulate HF broadening
                        ! (units: cm/s)
     Integer :: NLTEtransition ! Number of transition in the NLTE ATOM 
                               !  If this line is LTE, set to zero
     Real :: NLTE_nl_ratio ! Ratio of lower level populations for the line
                           ! of interest to the population in the NLTE 
                           ! computation (useful, e.g. to synthesize
                           ! multiplet components)
     Real :: NLTE_nu_ratio ! Same as above for the upper level
     Integer :: NLTEgridsize ! Size of wavelength grid for NLTE computation
     Real, Dimension(:), Allocatable :: NLTEgrid
     Real, Dimension(:,:), Allocatable :: NLTESource_f
     Real, Dimension(:), Allocatable :: b_low, b_up
  End Type Line_data
!
  Type Region_data
     Real :: First_wlength ! in Angstroms
     Real :: Wave_step ! in Angstroms
     Real :: Macro_enh ! Macroturbulence enhancement factor (dimensionless)
     Real :: Opacity_enh ! Fudge factor to enhance background opacities
     Real :: obs_additive ! This value will be added to the spectral profile at
                  ! all wavelengths in this region. Typically used
                  ! to compensate for spectrally-flat scattered light
                  ! inside the instrument
     Real :: obs_multiplicative ! This value will be applied to the spectral
                  ! profile (after the additive constant) at
                  ! all wavelengths in this region. Typically used
                  ! to compensate for spectrally-flat scattered light
                  ! inside the instrument
     Real :: obs_gauss_sigma ! If .ne. 0, use a Gaussian instrumental profile
                             ! for convolution of synthetic spectrum. This is
                             ! overriden if a detailed instrumental profile 
                             ! has been specified (units=Angstroms)
     Integer :: nwavelengths ! number of wavelengths
     Integer :: layer ! 1->Photospheric region; 2->Chromospheric region
  End Type Region_data
!
End Module Line_data_structure
