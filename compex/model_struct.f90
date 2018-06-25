Module Model_structure
!
  Interface assignment (=)
     Module procedure Model_assign, Model_assign_2comp
  End interface
!
  integer, parameter :: n_model_elements=92
  integer, parameter :: nvarsdepth=22, nvarssingle=13+n_model_elements
  Type Model
! temp is in K, v_los and v_mic are in cm/s, el_p and gas_p are in dyn/cm^2, 
! rho is in g/cm^3, b_long, b_x, b_y is in gauss, 
! z_scale is in km (positive is up, negative is down)
! ne, nH, nHplus are electron, neutral H and HII number densities in cm^-3
! nHminus, nH2 and nH2plus are the number densities in cm^-3 of the negative 
! H ion, molecular Hydrogen and positive molecular ion respectively
! v_mac is in cm/s and stray is dimensionless (0 < stray < 1).
     real, dimension(:), allocatable :: zz
     real, dimension(:), allocatable :: z_scale, ltau_500, temp, gas_p, rho, el_p
     real, dimension(:), allocatable :: v_los, v_mic, b_long, b_x, b_y
     real, dimension(:), allocatable :: bl_x, bl_y, bl_z, vl_x, vl_y, vl_z
     real, dimension(:), allocatable :: ne, nH, nHminus
     real, dimension(:), allocatable :: nHplus, nH2, nH2plus
     real :: v_mac, stray, ffactor, chrom_x, chrom_y
     real :: keep_el_p, keep_gas_p, keep_rho, keep_nH, keep_nHminus, keep_nHplus
     real :: keep_nH2, keep_nh2plus
     real, dimension(n_model_elements) :: Abundance
  End Type Model
  Type Model_2comp
     Type (Model) :: Comp1, Comp2
  End type Model_2comp
Contains
  Subroutine Allocate_model_2comp(Modelo, npoints)
    Implicit None
    Type (Model_2comp) :: Modelo
    Integer :: npoints
    Modelo%Comp1%ffactor=1.
    Modelo%Comp2%ffactor=0.
    Call Allocate_model(Modelo%Comp1,npoints)
    Call Allocate_model(Modelo%Comp2,npoints)
  End Subroutine Allocate_model_2comp

  Subroutine Allocate_model(Modelo, npoints)
  Use Param_structure
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Modelo
  Integer :: npoints, status
!  
  If (Allocated(Modelo%z_scale)) then
     If (Size(Modelo%z_scale) .ne. npoints) then
        Print *,'Error in model_struct.f90. Trying to allocate a model'
        Print *,'that is already allocated with a different size'
        Print *,'Existing model size:',Size(Modelo%z_scale)
        Print *,'Requested allocation size:',npoints
        Stop
     End if
     Return
  End if
  Allocate(Modelo%z_scale(npoints), stat=status)
  If (status .gt. 0) then
     Print *,'Out of memory while allocating Modelo%z_scale'
     Stop
  End if
  Allocate(Modelo%ltau_500(npoints), stat=status)
  If (status .gt. 0) then
     Print *,'Out of memory while allocating Modelo%ltau_500'
     Stop
  End if
  Allocate(Modelo%temp(npoints), stat=status)
  If (status .gt. 0) then
     Print *,'Out of memory while allocating Modelo%temp   '
     Stop
  End if
  Allocate(Modelo%gas_p(npoints), stat=status)
  If (status .gt. 0) then
     Print *,'Out of memory while allocating Modelo%gas_p'
     Stop
  End if
  Allocate(Modelo%rho(npoints), stat=status)
  If (status .gt. 0) then
     Print *,'Out of memory while allocating Modelo%rho'
     Stop
  End if
  Allocate(Modelo%el_p(npoints), stat=status)
  If (status .gt. 0) then
     Print *,'Out of memory while allocating Modelo%el_p'
     Stop
  End if
  Allocate(Modelo%v_los(npoints), stat=status)
  If (status .gt. 0) then
     Print *,'Out of memory while allocating Modelo%v_los'
     Stop
  End if
  Allocate(Modelo%v_mic(npoints), stat=status)
  If (status .gt. 0) then
     Print *,'Out of memory while allocating Modelo%v_mic'
     Stop
  End if
  Allocate(Modelo%b_long(npoints), stat=status)
  Allocate(Modelo%b_x(npoints), stat=status)
  Allocate(Modelo%b_y(npoints), stat=status)
  Allocate(Modelo%bl_x(npoints), stat=status)
  Allocate(Modelo%bl_y(npoints), stat=status)
  Allocate(Modelo%bl_z(npoints), stat=status)
  Allocate(Modelo%vl_x(npoints), stat=status)
  Allocate(Modelo%vl_y(npoints), stat=status)
  Allocate(Modelo%vl_z(npoints), stat=status)
  Allocate(Modelo%ne(npoints))
  Allocate(Modelo%nH(npoints))
  Allocate(Modelo%nHplus(npoints))
  Allocate(Modelo%nHminus(npoints))
  Allocate(Modelo%nH2(npoints))
  Allocate(Modelo%nH2plus(npoints))
!
  Return
End Subroutine Allocate_model

Subroutine DeAllocate_model_2comp(Modelo)
 Implicit None
  Type (Model_2comp) :: Modelo
  Call DeAllocate_model(Modelo%Comp1)
  Call DeAllocate_model(Modelo%Comp2)
  Return 
End Subroutine DeAllocate_model_2comp
!
Subroutine DeAllocate_model(Modelo)
  Use Param_structure

 Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Modelo

  If (Allocated(Modelo%z_scale)) Deallocate (Modelo%z_scale)
  If (Allocated(Modelo%ltau_500)) Deallocate (Modelo%ltau_500)
  If (Allocated(Modelo%temp)) Deallocate (Modelo%temp)
  If (Allocated(Modelo%gas_p)) Deallocate (Modelo%gas_p)
  If (Allocated(Modelo%rho)) Deallocate (Modelo%rho)
  If (Allocated(Modelo%el_p)) Deallocate (Modelo%el_p)
  If (Allocated(Modelo%v_los)) Deallocate (Modelo%v_los)
  If (Allocated(Modelo%v_mic)) Deallocate (Modelo%v_mic)
  If (Allocated(Modelo%b_long)) Deallocate (Modelo%b_long)
  If (Allocated(Modelo%b_x)) Deallocate (Modelo%b_x)
  If (Allocated(Modelo%b_y)) Deallocate (Modelo%b_y)
  If (Allocated(Modelo%bl_x)) Deallocate (Modelo%bl_x)
  If (Allocated(Modelo%bl_y)) Deallocate (Modelo%bl_y)
  If (Allocated(Modelo%bl_z)) Deallocate (Modelo%bl_z)
  If (Allocated(Modelo%vl_x)) Deallocate (Modelo%vl_x)
  If (Allocated(Modelo%vl_y)) Deallocate (Modelo%vl_y)
  If (Allocated(Modelo%vl_z)) Deallocate (Modelo%vl_z)
  If (Allocated(Modelo%ne)) Deallocate(Modelo%ne)
  If (Allocated(Modelo%nH)) Deallocate(Modelo%nH)
  If (Allocated(Modelo%nHplus)) Deallocate(Modelo%nHplus)
  If (Allocated(Modelo%nHminus)) Deallocate(Modelo%nHminus)
  If (Allocated(Modelo%nH2)) Deallocate(Modelo%nH2)
  If (Allocated(Modelo%nH2plus)) Deallocate(Modelo%nH2plus)
  Return
End Subroutine DeAllocate_model
!
!
Subroutine Model_assign(A,B)
  Type (Model), intent(inout) :: A
  Type (Model), intent(in) :: B
  Integer :: status, npoints
  !
  If (.NOT.Allocated(A%z_scale)) then
     npoints=size(B%z_scale)
     Call Allocate_Model(A, npoints)
  End if
  If (size(A%temp) .eq. 0) then
     print *,'Error in model_struct! Model is allocated, but has zero size'
     stop
  End if
  A%v_mac=B%v_mac
  A%stray=B%stray
  A%ffactor=B%ffactor
  A%chrom_x=B%chrom_x
  A%chrom_y=B%chrom_y
  A%z_scale=B%z_scale
  A%ltau_500=B%ltau_500
  A%temp=B%temp
  A%gas_p=B%gas_p
  A%rho=B%rho
  A%el_p=B%el_p
  A%v_los=B%v_los
  A%v_mic=B%v_mic
  A%B_long=B%B_long
  A%B_x=B%B_x
  A%B_y=B%B_y
  A%Bl_x=B%Bl_x
  A%Bl_y=B%Bl_y
  A%Bl_z=B%Bl_z
  A%vl_x=B%vl_x
  A%vl_y=B%vl_y
  A%vl_z=B%vl_z
  A%ne=B%ne
  A%nH=B%nH
  A%nHplus=B%nHplus
  A%nHminus=B%nHminus
  A%nH2=B%nH2
  A%nH2plus=B%nH2plus
  A%Keep_el_p=B%Keep_el_p
  A%Keep_gas_p=B%Keep_gas_p
  A%Keep_rho=B%Keep_rho
  A%Keep_nH=B%Keep_nH
  A%Keep_nHminus=B%Keep_nHminus
  A%Keep_nHplus=B%Keep_nHplus
  A%Keep_nH2=B%Keep_nH2
  A%Keep_nh2plus=B%Keep_nh2plus
  A%Abundance=B%Abundance
End subroutine Model_assign
!
Subroutine Model_assign_2comp(A,B)
  Type (Model_2comp), intent(inout) :: A
  Type (Model_2comp), intent(in) :: B
  Integer :: status, npoints
  !
  If (.NOT.Allocated(A%Comp1%z_scale)) then
     npoints=size(B%Comp1%z_scale)
     Call Allocate_Model(A%Comp1, npoints)
  End if
  If (.NOT.Allocated(A%Comp2%z_scale) .and. Allocated(B%Comp2%z_scale)) then
     npoints=size(B%Comp2%z_scale)
     Call Allocate_Model(A%Comp2, npoints)
  End if
  A%Comp1=B%Comp1
  If (Allocated(A%Comp2%z_scale)) &
       A%Comp2=B%Comp2
End Subroutine Model_assign_2comp
!
Subroutine Model_swap_2comp(A,B)
  Type (Model_2comp), intent(inout) :: A
  Type (Model_2comp), intent(in) :: B
  Type (Model_2comp) :: C
  C=B
  If (.NOT.Allocated(A%Comp1%z_scale)) then
     npoints=size(B%Comp1%z_scale)
     Call Allocate_Model(A%Comp1, npoints)
  End if
  If (.NOT.Allocated(A%Comp2%z_scale)) then
     npoints=size(B%Comp2%z_scale)
     Call Allocate_Model(A%Comp2, npoints)
  End if
  A%Comp1=C%Comp2
  A%Comp2=C%Comp1
!
  Return
!
End Subroutine Model_swap_2comp
!
End Module Model_structure
