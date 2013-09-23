! This module is used to handle debugging events
! Set Debug_level to one of the following:
!   0=No debugging and no stopping
!   1=No debugging. Stop at first error
!   2=Dump errors to core file (need to open file and assign it to
!        unit Debug_FileUnit before writing any logs). Code execution will
!        not stop in this mode
!   3=Same as above but also dump warnings to core file and stop code 
!        execution at the first error
!   4=Same as above but stop code execution at the first warning
!
! Flags can be used to transfer error conditions upstream so that the
! calling code can determine how to handle exceptions
!
Module Debug_module
  Implicit None
  Integer, Save :: Debug_level, Debug_FileUnit=-1, ierror=1, iwarn=1
  Integer, Parameter :: Debug_nwarns=1000, Debug_nerrors=10
  Integer, Dimension(100), Save :: Debug_errorflags, Debug_warningflags
  Integer, Parameter :: flag_lorien=1, flag_forward=2, flag_NLTE=3, &
       flag_computepe=4, flag_computepg=5, flag_computeothers=6, flag_svdcmp=7, &
       flag_eneq=8, flag_hydrostatic=9, flag_taunu=10, flag_ludcmp=11,  &
       flag_dchisq_dx=12, flag_background=13, flag_computeopac=14
  Character (Len=256), Dimension (Debug_nwarns), Save :: Debug_Warnings
  Character (Len=256), Dimension (Debug_nerrors), Save :: Debug_Errors
! Specific to NICOLE
  Logical, Save :: Debug_outputpop=.False., &
       Debug_outputcontop=.False., Debug_outputNLTEsf=.False.
!
  Contains
    Subroutine Debug_Log(Message, itype)
      Character (Len=*) :: Message
      Character (Len=4) :: Prefix
      Integer :: itype, iunit

      iunit=Debug_FileUnit
      If (iunit .lt. 6) iunit=6

      If (itype .eq. 1) then ! It's an error
         If (ierror .eq. debug_nerrors) then ! Buffer filled, scroll down
            ierror=debug_nerrors-1
            Debug_Errors(1:debug_nerrors-1)=Debug_Errors(2:debug_nerrors)
         End if
         Debug_Errors(ierror)=Message
         ierror=ierror+1
         Prefix='(EE)'
      Else if (itype .eq. 2) then ! It's a warninng
         If (iwarn .eq. debug_nwarns) then ! Buffer filled, scroll down
            iwarn=debug_nwarns-1
            Debug_Warnings(1:debug_nwarns-1)=Debug_Warnings(2:debug_nwarns)
         End if
         Debug_Warnings(iwarn)=Message
         iwarn=iwarn+1
         Prefix='(WW)'
      Else
         Print *,'Incorrect call to Debug_Log()'
         Stop
      End if

! Write to file?
      If ( (Debug_level .ge. 2 .and. itype .eq. 1) .or. Debug_level .ge. 3 ) then
         If (Debug_FileUnit .gt. 0) then
            Write (Debug_FileUnit, *) Prefix//Message
            Call MyFlush (Debug_FileUnit)
         End if
      End if

! Stop?
      If ( (Debug_level .ge. 2 .and. itype .eq. 1) .or. Debug_level .ge. 4) then
         If (iunit .gt. 0) then
            Write (iunit, *) 'Received request to abort code execution. First flags follow'
            Write (iunit, *) 'Error flags (1-10):',Debug_errorflags(1:10)
            Write (iunit, *) 'Error flags (11-20):',Debug_errorflags(11:20)
            Write (iunit, *) 'Warning flags (1-10):',Debug_warningflags(1:10)
            Write (iunit, *) 'Warning flags (11-20):',Debug_warningflags(11:20)
            Write (iunit, *) 'ABORTING...'
            Call MyFlush (iunit)
         End if
         Print *,''
         Print *,'ABORTING.. (Check core file for details)',debug_fileunit
         Print *,''
         Stop
      End if

    End Subroutine Debug_Log
            
End Module Debug_module
