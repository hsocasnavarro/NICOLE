Module UV_opacity_TOPbase_data
Use File_operations
Integer, Parameter :: nz=          18
Integer, Parameter :: maxnlevels=    1000
Integer :: nlevels(maxnlevels)
Integer, Dimension (nz) :: zs
Integer, Dimension (nz,maxnlevels) :: nlambda, ionstage
Real, Dimension (nz,maxnlevels) :: energy
Real, Dimension(nz, maxnlevels, 3500) :: XSec

Contains
  Subroutine UVOpacity_TOPbase_init
    Character (Len=49) :: Path='/home/d/navarro/nicole/git/NICOLE-master/forward/'
    Logical :: Exists, FirstTime=.True.
    Integer :: unit, iz, ilev, ilam
    !
    If (.not. FirstTime) Return
    !
    FirstTime=.False.
    Inquire(File=Path//'UV_TOP.dat',Exist=Exists)
    If (Exists) then
       Call Open_file(unit,Path//'UV_TOP.dat')
    Else
       Inquire(File='UV_TOP.dat',Exist=Exists)
       If (Exists) then
          Call Open_file(unit,'UV_TOP.dat')
       Else
          Inquire(File='UV_TOP.dat-fromgit/UV_TOP.dat',Exist=Exists)
          If (Exists) then
             Call Open_file(unit,'UV_TOP.dat-fromgit/UV_TOP.dat')
          Else
             Print *,'Could not find file UV_TOP.dat'
             Print *,'This file should be either in the source code tree that this'
             print *,'exectubale was compiled from or in the current running directory'
             print *,''
             print *,'The source tree where I was looking for it is:'
             print *,Path
             print *,''
             print *,'A simple solution is to download the following zip file and'
             print *,'extract (unzip) it here:'
             print *,'https://github.com/hsocasnavarro/UV_TOP.dat/archive/fromgit.zip'
             print *,''
             Stop
          End if
       End if
    End if
    !   
    Read (unit,*) nlevels(1:nz)
    Read (unit,*) zs(1:nz)
    Do iz=1, nz
       Do ilev=1,nlevels(iz)
          Read (unit,*) nlambda(iz,ilev), energy(iz,ilev), ionstage(iz,ilev)
          If (nlambda(iz,ilev) .ge. 1) &
               Read (unit,*) xsec(iz,ilev,1:nlambda(iz,ilev))
       End do
    End do
    !
    Call Close_file(unit)
    !
    Return
    !
  End Subroutine UVOpacity_TOPbase_init
!
End Module UV_opacity_TOPbase_data
