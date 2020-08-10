! Written by Hector Socas-Navarro
!
!
! Module to do code profiling, i.e. to track code execution times
!  and routine-by-routine analysis
! Usage:
!  In each routine you wish to profile, include a call to Time_routine
!   at the beginning and end of that routine. Time_routine has two input
!   arguments: Name is a string of up to 255 characters with the name
!   of the routine being profiled (this is just a label). The second
!   argument is a boolean variable that must be set to .True. if this
!   is the beginning of the execution of that routine or .False. if this
!   is the end (this determines whether we are starting or stopping the
!   stop watch). Make sure you use exactly the same name string in
!   both calls. IMPORTANT: One of the routines being profile must be named 
!   main and will be used as the reference to time the overall code execution. 
!   Normally you would include a call time_routine('main',.True.) near the
!   beginning of your code and call time_routine('main',.False.) near the end
! 
! For MPI:
!   To compile in MPI mode, uncomment the  lines. Then, the process
!   number is appended to the name of each routine with a double underscore
!   (e.g., if process 3 makes a call time_routine('myroutine',.True.) 
!    this will actually get translated as myroutne__0003)
! 
Module Profiling
  Private
  Public :: Do_profile, Time_routine, Time_log 
  Logical, Save :: Do_profile=.True.
  Integer, Save, Private :: nroutines=0, Count, Count_Rate, Count_Max, &
      myrank, status
  Character (len=8), Save, Private :: StartDate, EndDate
  Character (len=10), Save, Private :: StartTime, EndTime
  !
  Character (len=20), Private :: Suffix
  Character (len=255), Dimension(:), Allocatable, Save, Private :: RoutineName, TmpName
  Integer, Dimension(:), Allocatable, Save, Private :: RoutineStarted, TmpS
  Integer, Dimension(:), Allocatable, Save, Private :: RoutineAccumulated, TmpA
  Integer, Dimension(:), Allocatable, Save, Private :: RoutineNCalls, TmpNC
  !
  Include 'mpif.h'
  !
  Contains 
    Subroutine Time_routine(Name, Start)
      Implicit None
      Character (len=*) :: Name
      Logical :: Start
      Integer :: iroutine, ind
      !
      If (.not. Do_profile) Return
      !
      Call System_Clock(Count, Count_rate, Count_max)
      If (Count_max .lt. 0) then ! No clock present
         Print *,'No hardware clock found on this system. No profiling will be done'
         Do_profile=.False.
         Return
      End if
!
      If (Name .eq. 'main') then
         If (Start) then
            Call Date_and_Time(StartDate, StartTime)
         Else
            Call Date_and_Time(EndDate, EndTime)
         End If
      End if
!
      If (nroutines .eq. 0) then
         nroutines=1
         Allocate(RoutineName(1))
         Allocate(RoutineStarted(1))
         Allocate(RoutineAccumulated(1))
         Allocate(RoutineNCalls(1))
         RoutineName(1)=Name
         RoutineStarted(1)=Count
         RoutineAccumulated(1)=0
         RoutineNCalls(1)=1
         Return
      End if
!
      iroutine=-1
      Do ind=1, nroutines
         If (RoutineName(ind) .eq. Name) iroutine=ind
      End do
      If (iroutine .eq. -1) then ! New routine
         If (.not. Start) then ! Error. First call to routine must be start
            Print *,'Error in time_routine. Profiling new routine ',Name
            Print *,'But profile has been called with Start=.False.'
            Stop
         End if
         If (Allocated(TmpName)) then
            Deallocate(TmpName)
            Deallocate(TmpS)
            Deallocate(TmpA)
            Deallocate(TmpNC)
         End if
         Allocate(TmpName(nroutines))
         Allocate(TmpS(nroutines))
         Allocate(TmpA(nroutines))
         Allocate(TmpNC(nroutines))
         TmpName(1:nroutines)=RoutineName(1:nroutines)
         TmpS(1:nroutines)=RoutineStarted(1:nroutines)
         TmpA(1:nroutines)=RoutineAccumulated(1:nroutines)
         TmpNC(1:nroutines)=RoutineNCalls(1:nroutines)
         Deallocate(RoutineName)
         Deallocate(RoutineStarted)
         Deallocate(RoutineAccumulated)
         Deallocate(RoutineNCalls)
         nroutines=nroutines+1
         Allocate(RoutineName(nroutines))
         Allocate(RoutineStarted(nroutines))
         Allocate(RoutineAccumulated(nroutines))
         Allocate(RoutineNCalls(nroutines))
         RoutineName(1:nroutines-1)=TmpName(1:nroutines-1)
         RoutineStarted(1:nroutines-1)=TmpS(1:nroutines-1)
         RoutineAccumulated(1:nroutines-1)=TmpA(1:nroutines-1)
         RoutineNCalls(1:nroutines-1)=TmpNC(1:nroutines-1)
         iroutine=nroutines
         RoutineName(iroutine)=Name
         RoutineAccumulated(iroutine)=0
         RoutineStarted(iroutine)=Count
         RoutineNCalls(iroutine)=1
         Return
      End if
      If (Start) then ! Routine execution starts
         RoutineStarted(iroutine)=Count
         RoutineNCalls(iroutine)=RoutineNCalls(iroutine)+1
      Else ! Routine execution ends
         RoutineAccumulated(iroutine)=RoutineAccumulated(iroutine)+ &
              Count-RoutineStarted(iroutine)
      End if
      Return
    End Subroutine Time_routine
!
    Subroutine Time_log
      Implicit None
      Character (len=10) :: PrintDate
      Character (len=12) :: PrintTime
      Integer :: ind, iroutine, totaltime_int, iunit
      Real :: totaltime_real
      Real, Dimension(:), Allocatable :: Frac
      Integer, Dimension(1) :: index
      Logical :: Opened
      !
      If (.not. Do_profile) Return
      !
      myrank=0
     Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, status)
      Write(Suffix,'("__",i0)') myrank
      !
      ! Find available unit to write file
      iunit=29
      Opened=.True.
      Do while (Opened)
         Inquire (iunit, Opened = Opened)
      End do
      If (iunit .gt. 125) then
         Print *,'Could not find any file units available'
         Print *,'No profiling information will be written'
         Return
      End if
      Open (iunit, File='profiling'//Trim(Suffix)//'.txt')
      PrintDate=StartDate(1:4)//'/'//StartDate(5:6)//'/'//StartDate(7:8)      
      PrintTime=StartTime(1:2)//':'//StartTime(3:4)//':'//StartTime(5:10)      
      Write (iunit,*) 'Code execution started on: ',PrintDate,' at ',PrintTime
      PrintDate=EndDate(1:4)//'/'//EndDate(5:6)//'/'//EndDate(7:8)      
      PrintTime=EndTime(1:2)//':'//EndTime(3:4)//':'//EndTime(5:10)      
      Write (iunit,*) '  Code execution ended on: ',PrintDate,' at ',PrintTime
      iroutine=-1
      Do ind=1, nroutines
         If (Trim(RoutineName(ind)) .eq. 'main') iroutine=ind
      End do
      if (iroutine .ne. -1) then
         totaltime_real=RoutineAccumulated(iroutine)*1./Count_Rate
         totaltime_int=RoutineAccumulated(iroutine)*1./Count_Rate
         Write (iunit,*) 'Total execution time (s):',totaltime_int
         Write (iunit,*) 'Total execution time (hours):',totaltime_real/3600.
      Else
         Write (iunit,*) 'main was not profiled. Cannot determine total execution time'
         totaltime_real=-1e10
      End if
      !
      Allocate(Frac(nroutines))
      Write (iunit,*) 'Subroutines profiled   (Name, number of calls, time(s), % of total time)'
      Do ind=1, nroutines
         Frac(ind)=RoutineAccumulated(ind)*1./Count_Rate/totaltime_real
      End do
      Do ind=1, nroutines
         index=MaxLoc(Frac)
         Write (iunit,*) '   '//Trim(RoutineName(index(1))),RoutineNCalls(index), &
              RoutineAccumulated(index)*1./Count_Rate,Frac(index)*100
         Frac(index)=-1e11
      End do
      Close (iunit)
      Return
    End Subroutine Time_log
End Module Profiling
