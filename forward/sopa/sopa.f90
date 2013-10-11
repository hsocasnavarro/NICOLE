SUBROUTINE SOPAS(Id,Ihtot,Freq,Pel,Tel,Pg,Chio,Chie,Eta)
   IMPLICIT NONE
   INTEGER Id , Ihtot , MDEP , MFT , MK , NXTion
   REAL ATW , BK , BNU , CC , CHE , Chie , Chio , EE , EK , EM ,     &
         & Eta , Freq , FRQ , HC2 , HCE , HCK , HH , HTOt , PE , Pel
   REAL Pg , PI , RHO , sf , SIG , SUMA , Tel , TEMp , UU , WT

   Print *,'This build of NICOLE has not been compiled with support for the SOPA opacities'
   Print *,'Change opacities to the default option in your input file'
   Print *,'or compile with SOPA support by using the --sopa switch in create_makefile.py'
   Stop

End Subroutine SOPAS
      