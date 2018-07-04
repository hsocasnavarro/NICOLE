Subroutine Select_number_of_nodes(Params, Atmo, Nodes, icall)
  Use Nodes_information
  Use Param_structure
  Use Model_structure
  Implicit None
  Type (Parameters) :: Params
  Type (Model) :: Atmo
  Type (Nodes_info) :: Nodes
  Character (len = 256) :: String
  Integer :: status, icall, File_unit, i_pos, itmp, ivar, inode
  Integer :: Find_Index
  Logical :: Exists, End, Read
!
! Initialize number of nodes to -1
!
!
! Set all nodes in component 2 to zero
!
  If (Nodes%n_nodes_t2 .eq. -1) Nodes%n_nodes_t2=0
  If (Nodes%n_nodes_v2 .eq. -1) Nodes%n_nodes_v2=0
  If (Nodes%n_nodes_blong2 .eq. -1) Nodes%n_nodes_blong2=0
  If (Nodes%n_nodes_mic2 .eq. -1) Nodes%n_nodes_mic2=0
  If (Nodes%n_nodes_bx2 .eq. -1) Nodes%n_nodes_bx2=0
  If (Nodes%n_nodes_by2 .eq. -1) Nodes%n_nodes_by2=0
  If (Nodes%n_nodes_chrom_x2 .eq. -1) Nodes%n_nodes_chrom_x2=0
  If (Nodes%n_nodes_chrom_y2 .eq. -1) Nodes%n_nodes_chrom_y2=0

!  
  If (icall .eq. 1) then 
!
! For the moment, set 4 nodes in temperature and v, and 1 in the others
!
     If (Nodes%n_nodes_t .eq. -1) Nodes%n_nodes_t=4
     If (Nodes%n_nodes_v .eq. -1) Nodes%n_nodes_v=1
     If (Nodes%n_nodes_blong .eq. -1) Nodes%n_nodes_blong=1
     If (Nodes%n_nodes_mic .eq. -1) Nodes%n_nodes_mic=0
     If (Nodes%n_nodes_bx .eq. -1) Nodes%n_nodes_bx=1
     If (Nodes%n_nodes_by .eq. -1) Nodes%n_nodes_by=1
     If (Nodes%n_nodes_mac .eq. -1) Nodes%n_nodes_mac=0
     If (Nodes%n_nodes_stray .eq. -1) Nodes%n_nodes_stray=0
     If (Nodes%n_nodes_ffactor .eq. -1) Nodes%n_nodes_ffactor=0
     If (Nodes%n_nodes_chrom_x .eq. -1) Nodes%n_nodes_chrom_x=0
     If (Nodes%n_nodes_chrom_y .eq. -1) Nodes%n_nodes_chrom_y=0
     If (Maxval(Abs(Params%Stray_prof)) .lt. 1.e-10) &
          Nodes%n_nodes_stray=0
     If (Params%IProfExists) Nodes%n_nodes_mac=0
  End if
  If (icall .eq. 2) then
!
! Set 4 nodes in temperature, 2 in v, B_str, B_inc, and 1 in
! everything else
!
     If (Nodes%n_nodes_t .eq. -1) Nodes%n_nodes_t=6
     If (Nodes%n_nodes_v .eq. -1) Nodes%n_nodes_v=4
     If (Nodes%n_nodes_blong .eq. -1) Nodes%n_nodes_blong=4
     If (Nodes%n_nodes_mic .eq. -1) Nodes%n_nodes_mic=2
     If (Nodes%n_nodes_bx .eq. -1) Nodes%n_nodes_bx=2
     If (Nodes%n_nodes_by .eq. -1) Nodes%n_nodes_by=2
     If (Nodes%n_nodes_mac .eq. -1) Nodes%n_nodes_mac=1
     If (Nodes%n_nodes_stray .eq. -1) Nodes%n_nodes_stray=1
     If (Nodes%n_nodes_ffactor .eq. -1) Nodes%n_nodes_ffactor=1
     If (Nodes%n_nodes_chrom_x .eq. -1) Nodes%n_nodes_chrom_x=1
     If (Nodes%n_nodes_chrom_y .eq. -1) Nodes%n_nodes_chrom_y=1
     If (Maxval(Abs(Params%Stray_prof)) .lt. 1.e-10) &
          Nodes%n_nodes_stray=0
  End if
  If (icall .ge. 3) then
     If (Nodes%n_nodes_t .eq. -1) Nodes%n_nodes_t=8
     If (Nodes%n_nodes_v .eq. -1) Nodes%n_nodes_v=4
     If (Nodes%n_nodes_blong .eq. -1) Nodes%n_nodes_blong=2
     If (Nodes%n_nodes_mic .eq. -1) Nodes%n_nodes_mic=2
     If (Nodes%n_nodes_bx .eq. -1) Nodes%n_nodes_bx=2
     If (Nodes%n_nodes_by .eq. -1) Nodes%n_nodes_by=2
     If (Nodes%n_nodes_mac .eq. -1) Nodes%n_nodes_mac=1
     If (Nodes%n_nodes_stray .eq. -1) Nodes%n_nodes_stray=0
     If (Nodes%n_nodes_ffactor .eq. -1) Nodes%n_nodes_ffactor=0
     If (Nodes%n_nodes_chrom_x .eq. -1) Nodes%n_nodes_chrom_x=1
     If (Nodes%n_nodes_chrom_y .eq. -1) Nodes%n_nodes_chrom_y=1
     If (Maxval(Abs(Params%Stray_prof)) .lt. 1.e-10) &
          Nodes%n_nodes_stray=0
     If (Params%IProfExists) Nodes%n_nodes_mac=0
  End if
!  If (icall .gt. 3) then
!     Print *,'Too many cycles (',icall,')!!'
!     Stop
!  End if
  Params%n_free_parameters=Nodes%n_nodes_t + Nodes%n_nodes_v + &
       Nodes%n_nodes_blong + Nodes%n_nodes_mic + Nodes%n_nodes_bx + &
       Nodes%n_nodes_by + Nodes%n_nodes_mac + Nodes%n_nodes_stray + &
       Nodes%n_nodes_ffactor + Nodes%n_nodes_ab + &
       Nodes%n_nodes_chrom_x + Nodes%n_nodes_chrom_y
  Params%n_free_parameters=Params%n_free_parameters+ &
       Nodes%n_nodes_t2 + Nodes%n_nodes_v2 + &
       Nodes%n_nodes_blong2 + Nodes%n_nodes_mic2 + Nodes%n_nodes_bx2 + &
       Nodes%n_nodes_by2 + Nodes%n_nodes_ab2 + &
       Nodes%n_nodes_chrom_x2 + Nodes%n_nodes_chrom_y2

!
! Location of the nodes:
!
! Check if previously allocated. In that case deallocate them
!
  If (Associated(Nodes%i_nodes_t)) then
     Deallocate(Nodes%i_nodes_t)
     Deallocate(Nodes%i_nodes_v)
     Deallocate(Nodes%i_nodes_blong)
     Deallocate(Nodes%i_nodes_mic)
     Deallocate(Nodes%i_nodes_bx)
     Deallocate(Nodes%i_nodes_by)
     Deallocate(Nodes%i_nodes_t2)
     Deallocate(Nodes%i_nodes_v2)
     Deallocate(Nodes%i_nodes_blong2)
     Deallocate(Nodes%i_nodes_mic2)
     Deallocate(Nodes%i_nodes_bx2)
     Deallocate(Nodes%i_nodes_by2)
  End if
!
! Allocate memory
!
  Allocate(Nodes%i_nodes_t(Max(Nodes%n_nodes_t,1)), Stat=status)
  Allocate(Nodes%i_nodes_v(Max(Nodes%n_nodes_v,1)), Stat=status)
  Allocate(Nodes%i_nodes_blong(Max(Nodes%n_nodes_blong,1)), Stat=status)
  Allocate(Nodes%i_nodes_mic(Max(Nodes%n_nodes_mic,1)), Stat=status)
  Allocate(Nodes%i_nodes_bx(Max(Nodes%n_nodes_bx,1)), Stat=status)
  Allocate(Nodes%i_nodes_by(Max(Nodes%n_nodes_by,1)), Stat=status)
  Allocate(Nodes%i_nodes_t2(Max(Nodes%n_nodes_t,1)), Stat=status)
  Allocate(Nodes%i_nodes_v2(Max(Nodes%n_nodes_v,1)), Stat=status)
  Allocate(Nodes%i_nodes_blong2(Max(Nodes%n_nodes_blong,1)), Stat=status)
  Allocate(Nodes%i_nodes_mic2(Max(Nodes%n_nodes_mic,1)), Stat=status)
  Allocate(Nodes%i_nodes_bx2(Max(Nodes%n_nodes_bx,1)), Stat=status)
  Allocate(Nodes%i_nodes_by2(Max(Nodes%n_nodes_by,1)), Stat=status)
!
! Place the nodes
!
  ivar=1 ! T1
  If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_t=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_t
        Nodes%i_nodes_t(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     If (Nodes%n_nodes_chrom_x+Nodes%n_nodes_chrom_y .ge. 1) then
        Call place_nodes(Params%n_points, Nodes%n_nodes_t, Atmo%ltau_500, &
             Nodes%i_nodes_t,-3.)
     Else
        Call place_nodes(Params%n_points, Nodes%n_nodes_t, Atmo%ltau_500, &
             Nodes%i_nodes_t,-100.)
     Endif
  End if
  ivar=ivar+1 
   If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_v=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_v
        Nodes%i_nodes_v(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     Call place_nodes(Params%n_points, Nodes%n_nodes_v, Atmo%ltau_500, &
          Nodes%i_nodes_v,-100.)
        
  Endif
  ivar=ivar+1
  If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_mic=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_mic
        Nodes%i_nodes_mic(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     Call place_nodes(Params%n_points, Nodes%n_nodes_mic, Atmo%ltau_500, &
          Nodes%i_nodes_mic,-100.)
  End if
  ivar=ivar+1 
  If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_blong=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_blong
        Nodes%i_nodes_blong(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     Call place_nodes(Params%n_points, Nodes%n_nodes_blong, Atmo%ltau_500, &
          Nodes%i_nodes_blong,-100.)
  End if
  ivar=ivar+1
  If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_bx=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_bx
        Nodes%i_nodes_bx(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     Call place_nodes(Params%n_points, Nodes%n_nodes_bx, Atmo%ltau_500, &
          Nodes%i_nodes_bx,-100.)
  End if
  ivar=ivar+1
  If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_by=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_by
        Nodes%i_nodes_by(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     Call place_nodes(Params%n_points, Nodes%n_nodes_by, Atmo%ltau_500, &
          Nodes%i_nodes_by,-100.)
  End if
!
  ivar=7 ! T2
  If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_t2=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_t2
        Nodes%i_nodes_t2(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     If (Nodes%n_nodes_chrom_x2+Nodes%n_nodes_chrom_y2 .ge. 1) then
        Call place_nodes(Params%n_points, Nodes%n_nodes_t2, Atmo%ltau_500, &
             Nodes%i_nodes_t2,-3.)
     Else
        Call place_nodes(Params%n_points, Nodes%n_nodes_t2, Atmo%ltau_500, &
             Nodes%i_nodes_t2,-100.)
     Endif
  End if
  ivar=ivar+1
  If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_v2=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_v2
        Nodes%i_nodes_v2(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     Call place_nodes(Params%n_points, Nodes%n_nodes_v2, Atmo%ltau_500, &
          Nodes%i_nodes_v2,-100.)
  End if
  ivar=ivar+1
  If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_mic2=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_mic2
        Nodes%i_nodes_mic2(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     Call place_nodes(Params%n_points, Nodes%n_nodes_mic2, Atmo%ltau_500, &
          Nodes%i_nodes_mic2,-100.)
  End if
  ivar=ivar+1
  If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_blong2=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_blong2
        Nodes%i_nodes_blong2(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     Call place_nodes(Params%n_points, Nodes%n_nodes_blong2, Atmo%ltau_500, &
          Nodes%i_nodes_blong2,-100.)
  End if
  ivar=ivar+1
  If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_bx2=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_bx2
        Nodes%i_nodes_bx2(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     Call place_nodes(Params%n_points, Nodes%n_nodes_bx2, Atmo%ltau_500, &
          Nodes%i_nodes_bx2,-100.)
  End if
  ivar=ivar+1
  If (UserNodeLocations(ivar,1) .ge. 0) then
     Nodes%n_nodes_by2=UserNodeLocations(ivar,1)
     Do inode=1, Nodes%n_nodes_by2
        Nodes%i_nodes_by2(inode)=Find_index(Params%n_points, Atmo%ltau_500, &
             UserNodeLocations(ivar,inode+1))
     End do
  Else
     Call place_nodes(Params%n_points, Nodes%n_nodes_by2, Atmo%ltau_500, &
          Nodes%i_nodes_by2,-100.)
  End if
!
! Done!
!
  Return
End Subroutine Select_number_of_nodes
!
Subroutine Place_nodes(npoints, nnodes, x, inodes, upperbound)
  Implicit None
  Integer :: npoints, nnodes, Find_index
  Real, dimension (npoints) :: x
  Integer, dimension (nnodes) :: inodes
  Real :: xdistance, upperbound
  Integer :: ind, ilow, iup
  Real, dimension (nnodes) :: xnodes
!
! Equispace through model
!
  iup=1
  Do while (x(iup) .lt. upperbound .and. iup .lt. npoints)
     iup=iup+1
  End do
  ilow=npoints
!
!!$!
!!$! Leave boundaries without nodes
!!$!
!!$! Find the "boundaries" of the photosphere, defined as tau_500=10. and
!!$! tau_500=1.e-4
!!$!
!!$  iup=1
!!$  Do while (x(iup) .lt. -4. .and. iup .lt. npoints)
!!$     iup=iup+1
!!$  End do
!!$  ilow=npoints
!!$  Do while (x(ilow) .gt. 1 .and. ilow .gt. 1)
!!$     ilow=ilow-1
!!$  End do
!!$  If (iup .gt. ilow) then
!!$     Print *,'The model must span at least from log(tau)=1 to log(tau)=-4'
!!$     Stop
!!$  End if
!
! Now place the nodes
!
  If (nnodes .eq. 1) then ! If 1 node, set it in the middle point
     xnodes(1)=.5*(x(iup)+x(ilow))
     inodes(1)=Find_index(npoints, x, xnodes(1))
  Else if (nnodes .eq. 2) then ! If 2 nodes, set them at the boundaries
     xnodes(1)=x(iup)
     xnodes(2)=x(ilow)
     inodes(1)=Find_index(npoints, x, xnodes(1))
     inodes(2)=Find_index(npoints, x, xnodes(2))
  Else if (nnodes .ge. 3) then ! If 3 nodes or more, equispace them
     xdistance=(x(ilow)-x(iup))/(nnodes-1)
     xnodes(1)=x(iup)                            
     inodes(1)=Find_index(npoints, x, xnodes(1)) 
     Do ind=nnodes, 2, -1
        xnodes(ind)=x(ilow)-xdistance*(nnodes-ind)
        inodes(ind)=Find_index(npoints, x, xnodes(ind))
     End do
  End if
!  If (nnodes .ge. 6) then ! Set 3 in [-2, 0] and the rest outside
!     xnodes(nnodes-4)=-2
!     xnodes(nnodes-3)=-1
!     xnodes(nnodes-2)=-0.1
!     xnodes(nnodes-1)=0.1
!     xnodes(nnodes)=1.
!     Do ind=nnodes-4, nnodes
!        inodes(ind)=Find_index(npoints, x, xnodes(ind))
!     End do
!  End if
!
! Check that nodes don't overlap on the same grid point
!
  Do ind=2, nnodes
     If (inodes(ind) .eq. inodes(ind-1)) then
        print *,'Impossible to place nodes. Number of nodes: ',nnodes
        print *,'Proposed location in ltau_500:',xnodes(1:nnodes)
        print *,'Closest model gridpoints:',inodes(1:nnodes)
        print *,'Stop in select_number_of_nodes'
        Stop
     End if
     If (inodes(ind) .eq. inodes(ind-1)) then
        inodes(ind)=inodes(ind)+1
        If (inodes(ind) .gt. npoints) then
           print *,'Impossible to place nodes. Number of nodes: ',nnodes
           print *,'Stop in select_number_of_nodes'
           Stop
        End if
     End if
  End do
  Return
End Subroutine Place_nodes
!
Integer Function Find_index(npoints, x, x0)
!
! This function finds the index "i" such that x(i)<x0 .and. x(i+1)>x0
! It is assumed that the vector x is monotonically increasing.
!
  Implicit None
  Integer :: npoints
  Real, dimension (npoints) :: x
  Real :: x0
  Logical :: Found
!
  Found=.FALSE.
  If (x0 .le. x(1)) then
     Find_index=1
     Return
  End if
  If (x0 .ge. x(npoints)) then
     Find_index=npoints
     Return
  End If
  Find_index=0
  Do while (Find_index .le. npoints-1 .and. .not. Found)
     Find_index=Find_index+1
     If (x(Find_index) .le. x0 .and. x(Find_index+1) .gt. x0) Found=.TRUE.
  End do
  Return
End Function Find_index
