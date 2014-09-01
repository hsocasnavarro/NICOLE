Module File_operations
  Integer, Parameter :: length=256, NStored=500
  Character (len=length), Dimension(NStored), Save :: StoredNames
  Integer, Dimension(NStored), Save :: StoredUnits
  Integer :: NSaved=0

Contains
  ! Subroutine to open a file
  !
  Subroutine Open_file(Unit, Filename)
    !
    Implicit None
    Character (len=*) :: Filename
    Integer :: Unit
    Logical :: Op
    !
    If (Trim(Filename) .eq. '') Return
    Unit=29
    Op=.TRUE.
    Do while (Op)
       Unit=Unit+1
       Inquire (Unit, Opened = Op)
    End do
    If (Unit .gt. 125) Then
       Print *,'No more file units available'
       Stop
    End if
    Open (Unit, File=Filename)
    Return
  End Subroutine Open_file
  !
  ! Subroutine to open an unformatted file (sequential access)
  !
  Subroutine Open_file_unform(Unit, Filename)
    !
    Implicit None
    Character (len=*) :: Filename
    Integer :: Unit
    Logical :: Op
    !
    Unit=29
    Op=.TRUE.
    Do while (Op)
       Unit=Unit+1
       Inquire (Unit, Opened = Op)
    End do
    If (Unit .gt. 125) Then
       Print *,'No more file units available'
       Stop
    End if
    Open (Unit, File=Filename, Form='Unformatted')
    Return
  End Subroutine Open_file_unform

  !
  ! Subroutine to open an unformatted file for direct access
  !
  Subroutine Open_file_direct(Unit, Filename, Recl)
    !
    Implicit None
    Character (len=*) :: Filename
    Integer :: Unit, i
    Integer :: Recl
    Logical :: Op
    !
    ! First check if that file is already open
    !
    Do i=1, NSaved
       If (Trim(Filename) .eq. Trim(StoredNames(i))) then
          Unit=StoredUnits(i)
!          Print *,'Attempt to reopen file that was already open'
!          Print *,'Previously open unit:',Unit
!          Print *,'File name:',Trim(Filename)
!          Stop
          Return
       End if
    End do
    !
    ! Otherwise open the file
    !
    Unit=29
    Op=.TRUE.
    Do while (Op)
       Unit=Unit+1
       Inquire (Unit, Opened = Op)
    End do
    If (Unit .gt. 125) Then
       Print *,'No more file units available'
       Stop
    End if
    Open (Unit, File=Filename, Form='Unformatted', Access='Direct', Recl=Recl)
    NSaved=NSaved+1
    If (NSaved .lt. 500) then
       StoredNames(NSaved)=Trim(Filename)
       StoredUnits(NSaved)=Unit
    Else
       Print *,'Too many open files to keep track'
       Stop
    End if
    Return
  End Subroutine Open_file_direct

  Subroutine Close_File(Unit)
    Implicit None

    Integer :: Unit, i
    !
    If (Unit .le. 0) Return
    !
    ! First check if that file is already open
    !
    Do i=1, NSaved
       If (StoredUnits(i) .eq. Unit) then
          Close (StoredUnits(i))
          StoredUnits(i:NStored-1)=StoredUnits(i+1:NStored)
          StoredNames(i:NStored-1)=StoredNames(i+1:NStored)
          StoredUnits(NSaved)=-1
          StoredNames(NSaved)=''
          NSaved=NSaved-1
          Return
       End if
    End do
    !
    Close(Unit)
    !
    Return
  End Subroutine Close_File

  Subroutine Read_direct(CorrectEndian, iunit, irec, DataOut, sizedata, iost)
    !
    ! This routine performs an unformatted direct read of the file unit iunit
    ! (must have been opened previously) and reads the record number irec.
    ! The amount of data read depends on the size of the vector data. 
    !
    ! If the machine is big endian (signaled by setting Correctendian=.False.)
    ! the data in Data are automatically byte-swapped
    !
    ! If CorrectEndian is false then the byte order is swapped
    ! read/write a real number
    ! iunit is the unit number
    ! irec is the record number to read
    ! Data is the vector to read
    ! sizedata is the length of the vector
    ! iostat returns the status of the operation (0 is success)
    !
    Implicit None
    Logical :: Correctendian
    integer :: sizedata
    Real, Dimension(sizedata) :: DataOut
    Real (Kind=8), Dimension(sizedata) :: Data
    Real (Kind=8), Dimension(:), Allocatable :: Data_swapped
    Integer (Kind=1), Dimension(:), Allocatable :: Buffer
    Integer (Kind=1), Dimension(:), Allocatable :: Buffer_swapped
    Integer :: iost, irec, iunit, nbytes, idat, ibyte
    !
    If (Correctendian) then 
       Read (iunit, REC=irec, IOSTAT=iost) Data
    Else
       Read (iunit, REC=irec, IOSTAT=iost) Data
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
       Deallocate (Buffer)
       Deallocate (Buffer_swapped)
       Deallocate (Data_swapped)
    End if
    DataOut=Data 
    If (iost .ne. 0) &
         Print *,'ERROR!!! I/O error ',iost,' in direct read from unit ',iunit
    !
  End Subroutine Read_direct

  Subroutine Read_direct_int(CorrectEndian, iunit, irec, Data, sizedata, iost)
    !
    ! Same as above for reading integers
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
       Read (iunit, REC=irec, IOSTAT=iost) Data
    Else
       Read (iunit, REC=irec, IOSTAT=iost) Data
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
    End if
    !
  End Subroutine Read_direct_int


End Module File_operations
