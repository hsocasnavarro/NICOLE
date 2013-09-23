! Extract model from a 1D vector used to read/write records
! Use dir=1 to extract model from vector or dir=-1 for the opposite
Subroutine Record_to_model_2comp(np, Atmo, Vector1, Vector2, dir)
Use Model_structure
Implicit None
Type (Model_2comp) :: Atmo
Integer :: dir, np
Real, Dimension(nvarsdepth*np+nvarssingle) :: Vector1
Real, Dimension(nvarsdepth*np+nvarssingle) :: Vector2
!
Call Record_to_model(np, Atmo%Comp1, Vector1, dir)
Call Record_to_model(np, Atmo%Comp2, Vector2, dir)
Return
!
End Subroutine Record_to_model_2comp
!
Subroutine Record_to_model(np, Atmo, Vector, dir)
Use Model_structure
Implicit None
Type (Model) :: Atmo
Integer :: dir, np, ivar
Real, Dimension(nvarsdepth*np+nvarssingle) :: Vector
Real, Dimension(nvarsdepth*np+nvarssingle) :: Vector2
!
If (dir .eq. 1) then
   If ( Vector(np+1) .gt. Vector(np+2) .or.  & ! Reverse model if tau decreases
        Vector(1) .lt. Vector(2) ) then ! Reverse model if z increases
      Vector2=Vector(:)
      Do ivar=1, nvarsdepth-1
         Vector2((ivar-1)*np+1:ivar*np)=&
              Vector(ivar*np:(ivar-1)*np+1:-1)
      End do
      Vector(:)=Vector2
   End if
   Atmo%z_scale(1:np)=Vector(1:np)
   Atmo%ltau_500(1:np)=Vector(np+1:2*np)
   Atmo%temp(1:np)=Vector(2*np+1:3*np)
   Atmo%gas_p(1:np)=Vector(3*np+1:4*np)
   Atmo%rho(1:np)=Vector(4*np+1:5*np)
   Atmo%el_p(1:np)=Vector(5*np+1:6*np)
   Atmo%v_los(1:np)=Vector(6*np+1:7*np)
   Atmo%v_mic(1:np)=Vector(7*np+1:8*np)
   Atmo%b_long(1:np)=Vector(8*np+1:9*np)
   Atmo%b_x(1:np)=Vector(9*np+1:10*np)
   Atmo%b_y(1:np)=Vector(10*np+1:11*np)
   Atmo%bl_x(1:np)=Vector(11*np+1:12*np)
   Atmo%bl_y(1:np)=Vector(12*np+1:13*np)
   Atmo%bl_z(1:np)=Vector(13*np+1:14*np)
   Atmo%vl_x(1:np)=Vector(14*np+1:15*np)
   Atmo%vl_y(1:np)=Vector(15*np+1:16*np)
   Atmo%vl_z(1:np)=Vector(16*np+1:17*np)
   Atmo%nH(1:np)=Vector(17*np+1:18*np)
   Atmo%nHminus(1:np)=Vector(18*np+1:19*np)
   Atmo%nHplus(1:np)=Vector(19*np+1:20*np)
   Atmo%nH2(1:np)=Vector(20*np+1:21*np)
   Atmo%nH2plus(1:np)=Vector(21*np+1:22*np)
   Atmo%v_mac=Vector(nvarsdepth*np+1)
   Atmo%stray=Vector(nvarsdepth*np+2)
   Atmo%ffactor=Vector(nvarsdepth*np+3)
   Atmo%keep_el_p=Vector(nvarsdepth*np+4)
   Atmo%keep_gas_p=Vector(nvarsdepth*np+5)
   Atmo%keep_rho=Vector(nvarsdepth*np+6)
   Atmo%keep_nH=Vector(nvarsdepth*np+7)
   Atmo%keep_nHminus=Vector(nvarsdepth*np+8)
   Atmo%keep_nHplus=Vector(nvarsdepth*np+9)
   Atmo%keep_nH2=Vector(nvarsdepth*np+10)
   Atmo%keep_nh2plus=Vector(nvarsdepth*np+11)
   Atmo%Abundance(:)=Vector(nvarsdepth*np+11+1:nvarsdepth*np+11+n_model_elements)
Else if (dir .eq. -1) then
   Vector(1:np)=Atmo%z_scale(1:np)
   Vector(np+1:2*np)=Atmo%ltau_500(1:np)
   Vector(2*np+1:3*np)=Atmo%temp(1:np)
   Vector(3*np+1:4*np)=Atmo%gas_p(1:np)
   Vector(4*np+1:5*np)=Atmo%rho(1:np)
   Vector(5*np+1:6*np)=Atmo%el_p(1:np)
   Vector(6*np+1:7*np)=Atmo%v_los(1:np)
   Vector(7*np+1:8*np)=Atmo%v_mic(1:np)
   Vector(8*np+1:9*np)=Atmo%b_long(1:np)
   Vector(9*np+1:10*np)=Atmo%b_x(1:np)
   Vector(10*np+1:11*np)=Atmo%b_y(1:np)
   Vector(11*np+1:12*np)=Atmo%bl_x(1:np)
   Vector(12*np+1:13*np)=Atmo%bl_y(1:np)
   Vector(13*np+1:14*np)=Atmo%bl_z(1:np)
   Vector(14*np+1:15*np)=Atmo%vl_x(1:np)
   Vector(15*np+1:16*np)=Atmo%vl_y(1:np)
   Vector(16*np+1:17*np)=Atmo%vl_z(1:np)
   Vector(17*np+1:18*np)=Atmo%nH(1:np)
   Vector(18*np+1:19*np)=Atmo%nHminus(1:np)
   Vector(19*np+1:20*np)=Atmo%nHplus(1:np)
   Vector(20*np+1:21*np)=Atmo%nH2(1:np)
   Vector(21*np+1:22*np)=Atmo%nH2plus(1:np)
   Vector(nvarsdepth*np+1)=Atmo%v_mac
   Vector(nvarsdepth*np+2)=Atmo%stray
   Vector(nvarsdepth*np+3)=Atmo%ffactor
   Vector(nvarsdepth*np+4)=Atmo%keep_el_p
   Vector(nvarsdepth*np+5)=Atmo%keep_gas_p
   Vector(nvarsdepth*np+6)=Atmo%keep_rho
   Vector(nvarsdepth*np+7)=Atmo%keep_nH
   Vector(nvarsdepth*np+8)=Atmo%keep_nHminus
   Vector(nvarsdepth*np+9)=Atmo%keep_nHplus
   Vector(nvarsdepth*np+10)=Atmo%keep_nH2
   Vector(nvarsdepth*np+11)=Atmo%keep_nh2plus
   Vector(nvarsdepth*np+11+1:nvarsdepth*np+11+n_model_elements)= &
        Atmo%Abundance(:)
Else
   Print *,'Unknown value for dir in compex/record_to_model.f90'
   Stop
End if
!
End Subroutine Record_to_model
