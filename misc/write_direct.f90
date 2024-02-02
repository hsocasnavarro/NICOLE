Subroutine Write_direct(CorrectEndian, iunit, irec, Datain, sizedata, iost)
!
! This routine performs an unformatted direct write of the file unit iunit
! (must have been opened previously) and writes the record number irec.
! The amount of data written depends on the size of the vector data. 
!
! If the machine is big endian (signaled by setting Correctendian=.False.)
! the data in Data are automatically byte-swapped
!
! If CorrectEndian is false then the byte order is swapped
! iunit is the unit number
! irec is the record number to write
! Data is the vector to write
! sizedata is the length of the vector
! iostat returns the status of the operation (0 is success)
!
Implicit None
Logical :: Correctendian
integer :: sizedata
Real, Dimension(:), Intent(In) :: DataIn
Real (Kind=8) , Dimension(sizedata) :: Data ! 64-bit data to write
Real (Kind=8), Dimension(:), Allocatable :: Data_swapped
Integer (Kind=1), Dimension(:), Allocatable :: Buffer
Integer (Kind=1), Dimension(:), Allocatable :: Buffer_swapped
Integer :: iost, irec, iunit, nbytes, idat, ibyte
!
Data=DataIn ! Put in 64-bit format
If (Correctendian) then
   Write (iunit, REC=irec, IOSTAT=iost) Data
Else
   nbytes=kind(Data)
   Allocate (Buffer(nbytes))
   Allocate (Buffer_swapped(nbytes))
   Allocate (Data_swapped(sizeData))
   Do idat=1, sizeData
      Buffer=Transfer(Data(idat),Buffer)
      Buffer_swapped=Buffer(nbytes:1:-1)
      Data_swapped(idat)=Transfer(Buffer_swapped,Data_swapped(idat))
   End do
   Data=Data_swapped
   Write (iunit, REC=irec, IOSTAT=iost) Data
   Deallocate(Buffer)
   Deallocate(Buffer_swapped)
   Deallocate(Data_swapped)
End if
!
End Subroutine Write_direct

Subroutine Write_direct_int(CorrectEndian, iunit, irec, Data, sizedata, iost)
!
! Just like the above but data are integers
!
Implicit None
Logical :: Correctendian
integer :: sizedata
Integer (Kind=8), Dimension(sizedata) :: Data
Integer (Kind=8), Dimension(:), Allocatable :: Data_swapped
Integer (Kind=1), Dimension(:), Allocatable :: Buffer
Integer (Kind=1), Dimension(:), Allocatable :: Buffer_swapped
Integer :: iost, irec, iunit, nbytes, idat, ibyte
!
If (Correctendian) then 
   Write (iunit, REC=irec, IOSTAT=iost) Data
Else
   nbytes=kind(Data)
   Allocate (Buffer(nbytes))
   Allocate (Buffer_swapped(nbytes))
   Allocate (Data_swapped(sizeData))
   Do idat=1, sizeData
      Buffer=Transfer(Data(idat),Buffer)
      Buffer_swapped=Buffer(nbytes:1:-1)
      Data_swapped(idat)=Transfer(Buffer_swapped,Data_swapped(idat))
   End do
   Data=Data_swapped
   Write (iunit, REC=irec, IOSTAT=iost) Data
End if
!
End Subroutine Write_direct_int
