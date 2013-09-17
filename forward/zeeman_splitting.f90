Module Zeeman_Splitting

Contains 
subroutine MVOIGT2(nr,nwave,dlr,sr,a,v,hh,wc,dldop,etar,esar)
  
  Parameter (Piis=0.56418956659890107D0)
  double precision :: etar(nwave),esar(nwave)
!vetar,getar,ettar,ettvr,esar,vesar,gesar
!	Real (kind=8) essar,essvr,ettmr,essmr
!	Real (kind=8) da, dver, dH, dF, HV, FV, FA, HA
  REAL  DLR(nr),SR(nr)
  Real :: WC
  Real, Dimension(nwave) :: v, ver, H, F
  Integer :: nwave

! dlr es el desplazamiento en cm/gauss de la componete r Zeemann
! dlr es negativo luego debo cambiar de signo la exprsion
  etar=0.    !eta para luz circular derecha
  ESAR=0.    !perfil anomalo
  
  w1=wc/dldop
  
  do ir=1,nr   !do en el numero de componentes r zeeman
     ver(1:nwave)=v(1:nwave)+dlr(ir)*HH
     !	         H=VOIGT(a,ver,0)     !funcion de Voigt
     !	         F=VOIGT(a,ver,1)     !funcion de Faraday
     
     call voigt3(nwave,a,ver,H,F) ! Vector routine
     
     !	         HV=-2.*ver*dH+4.*da*dF  !derivada con v de H
     !	         FV=piis-da*dH-2.*dver*dF  !derivada con v de F
     !                 FA=hv/2.             !derivada con a de F
     !		 HA=-2.*fv    !derivada con a de H
     etar(:)=etar(:)+sr(ir)*H(:)
     esar(:)=esar(:)+sr(ir)*F(:)
     !	         ettaR=ettaR+sr(ir)*HA
     !	         essaR=essaR+sr(ir)*FA
     !		 ettvR=ettvR-sr(ir)*HV*T13*dver
     !		 essvR=ettvR-sr(ir)*FV*T13*dver
     !	         ettmR=ettmR-sr(ir)*HV*t14*dver
     !	         essmR=essmR-sr(ir)*FV*t14*dver
     
     !		 vetar=vetar-sr(ir)*HV*w1
     !		 vesar=vesar-sr(ir)*FV*w1
     !		 getar=getar+sr(ir)*HV*dlr(ir)/dldop
     !		 gesar=gesar+sr(ir)*FV*dlr(ir)/dldop
  End do         !fin del do en componentes r Zeeman
  
  ESAR=2.*ESAR
  !	      VESAR=2.*VESAR
  !	      GESAR=2.*GESAR
  !	      ESSAR=2.*ESSAR
  !	      ESSVR=2.*ESSVR
  !	      essmr=2.*essmr
  
  return
end subroutine MVOIGT2

subroutine MVOIGT(nr,dlr,sr,a,v,hh,t13,t14,wc,dldop,etar,vetar, &
     getar,ettar,ettvr,esar,vesar,gesar,essar,essvr,ettmr,essmr)
  
  Parameter (Piis=0.56418956659890107D0)
  REAL WC
  double precision :: etar,vetar,getar,ettar,ettvr,esar,vesar,gesar
  double precision :: essar,essvr,ettmr,essmr
  double precision :: da, dver, dH, dF, HV, FV, FA, HA
  REAL DLR(nr),SR(nr)
  
  ! dlr es el desplazamiento en cm/gauss de la componete r Zeemann
  ! dlr es negativo luego debo cambiar de signo la exprsion
  etar=0.    !eta para luz circular derecha
  !	      vetar=0.   !parcial de eta con el campo de velocidad(s/cm)
  !	      getar=0.   !parcial de eta con el campo magnetico.
  !	      ETTAR=0.   !parcial de eta con a
  !	      ETTVR=0.   !parcial de eta con v
  !	      ettmr=0.   !parcial de eta con mic
  ESAR=0.    !perfil anomalo
  !	      VESAR=0.   !su derivada con la velocidad
  !	      GESAR=0.   !su derivada con el campo
  !	      ESSVR=0.   !parcial de esar con v
  !	      ESSAR=0.   !parcial de esar con a
  !	      essmr=0.   !parcial de esar con mic
  
  w1=wc/dldop
  
  do 102 ir=1,nr   !do en el numero de componentes r zeeman
     ver=v+dlr(ir)*HH
     !	         H=VOIGT(a,ver,0)     !funcion de Voigt
     !	         F=VOIGT(a,ver,1)     !funcion de Faraday
     
     da=a
     dver=ver
     
     call voigt2(da,dver,dH,dF)
     
     !	         HV=-2.*ver*dH+4.*da*dF  !derivada con v de H
     !	         FV=piis-da*dH-2.*dver*dF  !derivada con v de F
     !                 FA=hv/2.             !derivada con a de F
     !		 HA=-2.*fv    !derivada con a de H
     etar=etar+sr(ir)*dH
     esar=esar+sr(ir)*dF
     !	         ettaR=ettaR+sr(ir)*HA
     !	         essaR=essaR+sr(ir)*FA
     !		 ettvR=ettvR-sr(ir)*HV*T13*dver
     !		 essvR=ettvR-sr(ir)*FV*T13*dver
     !	         ettmR=ettmR-sr(ir)*HV*t14*dver
     !	         essmR=essmR-sr(ir)*FV*t14*dver
     
     !		 vetar=vetar-sr(ir)*HV*w1
     !		 vesar=vesar-sr(ir)*FV*w1
     !		 getar=getar+sr(ir)*HV*dlr(ir)/dldop
     !		 gesar=gesar+sr(ir)*FV*dlr(ir)/dldop
102  continue         !fin del do en componentes r Zeeman
     
     ESAR=2.*ESAR
     !	      VESAR=2.*VESAR
     !	      GESAR=2.*GESAR
     !	      ESSAR=2.*ESSAR
     !	      ESSVR=2.*ESSVR
     !	      essmr=2.*essmr
     
     
     return

end subroutine MVOIGT
!*****************************************************************  
! This vectorizable Voigt function is based on the paper by Hui, *  
! Armstrong and Wray, JQSRT 19, 509 (1977). Errors become        *  
! significant (at the < 1% level) around the knee between the    *  
! Doppler core and the damping wings for 0.0 < a < 0.001. The    *  
! normalization is such that the integral is equal to SQRT(PI).  *  
! If J <> 0 this function returns the dispersion function.       *  
! -------------------------------------------------------------- *  
! Authors: Jack Harvey, Aake Nordlund.                           *  
! Modified by Sami Solanki (1985).                               *  
! Modified by A.D. Wittmann (1986) to include F(a,-v) and F(0,v) *  
! -------------------------------------------------------------- *  
! Last Update: 18-APR-86.                                        *  
!*****************************************************************  
      FUNCTION VOIGT(A,VV,J)
      COMPLEX Z 
      DIMENSION XDWS(28),YDWS(28)   
      DATA A0,A1,A2,A3,A4,A5,A6,B0,B1,B2,B3,B4,B5,B6/   &
      122.607931777104326,214.382388694706425,181.928533092181549,  &
      93.155580458138441,30.180142196210589,5.912626209773153,  &
      .564189583562615,122.60793177387535,352.730625110963558,  &
      457.334478783897737,348.703917719495792,170.354001821091472,  &
      53.992906912940207,10.479857114260399/
      DATA XDWS/.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.2,1.4,1.6,1.8,2.,   &
      3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20./,YDWS/  &
      9.9335991E-02,1.9475104E-01,2.8263167E-01,3.5994348E-01, &
      4.2443639E-01,4.7476321E-01,5.1050407E-01,5.3210169E-01, &
      5.4072434E-01,5.3807950E-01,5.0727350E-01,4.5650724E-01, &
      3.9993989E-01,3.4677279E-01,3.0134040E-01,1.7827103E-01, &
      1.2934799E-01,1.0213407E-01,8.4542692E-02,7.2180972E-02, &
      6.3000202E-02,5.5905048E-02,5.0253846E-02,4.1812878E-02, &
      3.5806101E-02,3.1311397E-02,2.7820844E-02,2.5031367E-02/ 

 
      V=ABS(VV) 
      IF(A.NE.0) GOTO 1 
      IF(J.NE.0) GOTO 3 
      VOIGT=EXP(-V*V)   
      RETURN
   3  IF(V.GT.XDWS(1)) GOTO 4   
      D=V*(1.-.66666667*V*V)
      GOTO 8
   4  IF(V.GT.XDWS(28)) GOTO 5  
      K=27  
      DO 7 I=2,27   
      IF(XDWS(I).LT.V) GOTO 7   
      K=I   
      GOTO 6
   7  CONTINUE  
   6  KK=K-1
      KKK=K+1   
      D1=V-XDWS(KK) 
      D2=V-XDWS(K)  
      D3=V-XDWS(KKK)
      D12=XDWS(KK)-XDWS(K)  
      D13=XDWS(KK)-XDWS(KKK)
      D23=XDWS(K)-XDWS(KKK) 
      D=YDWS(KK)*D2*D3/(D12*D13)-YDWS(K)*D1*D3/(D12*D23)+YDWS(KKK)* &
      D1*D2/(D13*D23)   
      GOTO 8
   5  Y=.5/V
      D=Y*(1.+Y/V)  
   8  VOIGT=5.641895836E-1*D
   9  IF(VV.LT.0.) VOIGT=-VOIGT 
      RETURN
   1  Z=CMPLX(A,-V) 
      Z=((((((A6*Z+A5)*Z+A4)*Z+A3)*Z+A2)*Z+A1)*Z+A0)/   &
      (((((((Z+B6)*Z+B5)*Z+B4)*Z+B3)*Z+B2)*Z+B1)*Z+B0)  
      IF(J.NE.0) GOTO 2 
      VOIGT=REAL(Z)
      RETURN
   2  VOIGT=.5*AIMAG(Z) 

      GOTO 9
    END FUNCTION VOIGT
  
!*****************************************************************  
! This vectorizable Voigt function is based on the paper by Hui, *  
! Armstrong and Wray, JQSRT 19, 509 (1977). Errors become        *  
! significant (at the < 1% level) around the knee between the    *  
! Doppler core and the damping wings for 0.0 < a < 0.001. The    *  
! normalization is such that the integral is equal to SQRT(PI).  *  
! If J <> 0 this function returns the dispersion function.       *  
! -------------------------------------------------------------- *  
! Authors: Jack Harvey, Aake Nordlund.                           *   
! Modified by Sami Solanki (1985).                               *  
! Modified by A.D. Wittmann (1986) to include F(a,-v) and F(0,v) *
! convertida en rutina por Basilio Ruiz (1993)                   *
! -------------------------------------------------------------- *  
! Last Update: 22-jun 93.                                        *  
!*****************************************************************  
      subroutine VOIGT2(A,VV,H,F)
      Implicit Double Precision (a-h,o-z)
      COMPLEX Z 
      DIMENSION XDWS(28),YDWS(28)
      DATA A0,A1,A2,A3,A4,A5,A6,B0,B1,B2,B3,B4,B5,B6/   &
      122.607931777104326,214.382388694706425,181.928533092181549,  &
      93.155580458138441,30.180142196210589,5.912626209773153,  &
      .564189583562615,122.60793177387535,352.730625110963558,  &
      457.334478783897737,348.703917719495792,170.354001821091472, &
      53.992906912940207,10.479857114260399/
      DATA XDWS/.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.2,1.4,1.6,1.8,2.,   &
      3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20./,YDWS/  &
      9.9335991E-02,1.9475104E-01,2.8263167E-01,3.5994348E-01,  &
      4.2443639E-01,4.7476321E-01,5.1050407E-01,5.3210169E-01,  &
      5.4072434E-01,5.3807950E-01,5.0727350E-01,4.5650724E-01,  &
      3.9993989E-01,3.4677279E-01,3.0134040E-01,1.7827103E-01,  &
      1.2934799E-01,1.0213407E-01,8.4542692E-02,7.2180972E-02,  &
      6.3000202E-02,5.5905048E-02,5.0253846E-02,4.1812878E-02,  &
      3.5806101E-02,3.1311397E-02,2.7820844E-02,2.5031367E-02/
      
      
      
      ivsigno=1
      if(vv.lt.0)ivsigno=-1
      v=ivsigno*vv
      
      IF(A.eq.0)then 
         v2=V*V
         H=EXP(-v2)   
         
         IF(V.GT.XDWS(1)) GOTO 4   
         D=V*(1.-.66666667*v2)
         GOTO 8
4        IF(V.GT.XDWS(28)) GOTO 5  
         K=27  
         DO 7 I=2,27   
            IF(XDWS(I).LT.V) GOTO 7   
            K=I   
            GOTO 6
7        Continue
6        KK=K-1
         KKK=K+1   
         D1=V-XDWS(KK) 
         D2=V-XDWS(K)  
         D3=V-XDWS(KKK)
         D12=XDWS(KK)-XDWS(K)  
         D13=XDWS(KK)-XDWS(KKK)
         D23=XDWS(K)-XDWS(KKK) 
         D=YDWS(KK)*D2*D3/(D12*D13)-YDWS(K)*D1*D3/(D12*D23)+YDWS(KKK)* &
              D1*D2/(D13*D23)   
         GOTO 8
5        Y=.5/V
         D=Y*(1.+Y/V)  
8        F=ivsigno*5.641895836E-1*D
         
         ! si el damping no es nulo
      else
         
         Z=CMPLX(A,-V) 
         Z=((((((A6*Z+A5)*Z+A4)*Z+A3)*Z+A2)*Z+A1)*Z+A0)/   &
              (((((((Z+B6)*Z+B5)*Z+B4)*Z+B3)*Z+B2)*Z+B1)*Z+B0)
         
         H=REAL(Z) 
         F=.5*ivsigno*AIMAG(Z) 
         
      end if
      
      return
    END subroutine VOIGT2
    


!***********************************************************************
!      SUBROUTINE HUMLICEK (NX,X,Y, PRBFCT)
    Complex Function approx1(T)
      Complex :: T
      APPROX1   = (T * .5641896) / (.5 + (T * T))
      Return
    End function approx1
    Complex Function approx2(T,U)
      Complex :: T, U
      APPROX2 = (T * (1.410474 + U*.5641896))/ (.75 + (U *(3.+U)))
      Return
    End function approx2
    Complex Function approx3(T)
      Complex :: T
      APPROX3   = ( 16.4955 + T * (20.20933 + T * (11.96482 + &
           T * (3.778987 + 0.5642236*T)))) &
           / ( 16.4955 + T * (38.82363 + T * &
           (39.27121 + T * (21.69274 + T * (6.699398 + T)))))   
      Return
    End Function approx3
    Complex Function approx4(T,U)
      Complex :: T, U
      APPROX4 = (T * (36183.31 - U * (3321.99 - U * (1540.787 - U &
           *(219.031 - U *(35.7668 - U *(1.320522 - U * .56419)))))) &
           / (32066.6 - U * (24322.8 - U * (9022.23 - U * (2186.18  &
           - U * (364.219 - U * (61.5704 - U * (1.84144 - U))))))))   
      Return
    End Function approx4
    subroutine voigt3 (nwave, a,vv,H,F)
!                                                                      *
!     complex probability function for complex argument Z=X+iY         *
!     real part = voigt function K(x,y)                                *
!                                                                      *
!     source:   j. humlicek, JQSRT 27, 437, 1982                       *
!                                                                      *
!     parameters:                                                      *
!      NX     number of grid points = NX+1                         in  *
!      X      array of grid points                                 in  *
!      Y      Voigt function parameter, ratio of lorentz/doppler   in  *
!      PRBFCT complex array of function values                     out *
!                                                                      * 
!     the stated accuracy is claimed to be 1.0E-04 by the author.      *
!     r h norton has checked the accuracy by comparing values          *
!     computed using a program written by b.h.armstrong, and           *
!     the accuracy claim seems to be warranted.                        *
!                                                                      *
!*************************************************************fgs 12/91*
      REAL      X(0:nwave-1), Y, S
      COMPLEX   T, U, PRBFCT(0:nwave-1)
!      COMPLEX   APPROX1, APPROX2, APPROX3, APPROX4
      
      Real a,vv(nwave),h(nwave),f(nwave)
      Integer nwave
!     ==================================================================

      NX=nwave-1
      X(0:nwave-1)=vv(1:nwave)
      Y=a
      
      F(1:nwave)=1
      Do i=0, nwave-1
         If (X(i) .lt. 0) F(i+1)=-1
         X(i)=Abs(X(i))
      End do
      
      !     ==================================================================
      
      IF (Y.GT.15.) THEN
         !        ---------------------------------------------------------------
         !        all points are in region I
         DO  I=0,NX
            T         = CMPLX(Y,-X(I))
100         PRBFCT(I) = APPROX1(T)
         End do
         !        ---------------------------------------------------------------
      ELSE IF (Y.LT.15. .AND. Y.GE.5.5) THEN
         !        ---------------------------------------------------------------
         !        points are in region I or region II
         DO 200, I=0,NX
            T  = CMPLX(Y,-X(I))
            S  = ABS(X(I)) + Y
            IF (S .GE. 15.) THEN
               PRBFCT(I) = APPROX1(T)
            ELSE
               U      = T * T
               PRBFCT(I) = APPROX2(T,U)
            END IF
200         CONTINUE
            !       ----------------------------------------------------------------
         ELSE IF (Y.LT.5.5 .AND. Y.GT.0.75) THEN
            !        ---------------------------------------------------------------
            DO 300, I=0,NX
               T  = CMPLX(Y,-X(I))
               S  = ABS(X(I)) + Y
               IF (S .GE. 15.) THEN
                  PRBFCT(I) = APPROX1(T)
               ELSE IF (S.LT.5.5) THEN
                  PRBFCT(I) = APPROX3(T)
               ELSE
                  U      = T * T
                  PRBFCT(I) = APPROX2(T,U)
               END IF
300            CONTINUE
               !       ----------------------------------------------------------------
            ELSE
               !        ---------------------------------------------------------------
               DO 400, I=0,NX
                  T  = CMPLX(Y,-X(I))
                  AX = ABS(X(I))
                  S  = AX + Y
                  IF (S .GE. 15.) THEN
                     !              region I
                     PRBFCT(I)= APPROX1(T)
                  ELSE IF (S.LT.15.0 .AND. S.GE.5.5) THEN
                     !              region II
                     U = T * T
                     PRBFCT(I)= APPROX2(T,U)
                  ELSE IF (S.LT.5.5 .AND. Y.GE.(0.195*AX-0.176)) THEN
                     !             region III
                     PRBFCT(I)= APPROX3(T)
                  ELSE
                     !             region IV
                     U  = T * T
                     PRBFCT(I)= CEXP(U) - APPROX4(T,U)
                  END IF
400               CONTINUE
                  !        ---------------------------------------------------------------
               END IF
               !     ==================================================================
               IF (Y .EQ. 0.0) THEN
                  DO  I=0,NX
20                   PRBFCT(I) = CMPLX(EXP(-X(I)**2), AIMAG(PRBFCT(I)))
                  END DO
               END IF
               !     ==================================================================
               
               
               h(1:nwave)=real(prbfct(0:nwave-1))
               f(1:nwave)=.5*f(1:nwave)*aimag(prbfct(0:nwave-1))
               
               RETURN
End subroutine voigt3


!       code1 representa mediante letras minusculas los momentos angulares
! Note: orbitales semienteros disabled by uppercase conversion in
!         routine read_line_data.f90
!       orbitales semienteros, asi p=1/2, f=3/2, h=5/2,k=7/2,m=9/2,o=11/2
!                              r=13/2,t=15/2, u=17/2,v=19/2,w=21/2      
SUBROUTINE ZEEMAN(MC,MULT,DESIGN,TAM,JI,JF,DL0,NP,NL,NR,DLP,DLL,&
     DLR,SP,SL,SR) 
  CHARACTER (Len=1) :: DESIGN
  Character, Dimension(13) :: CODE
  Data Code /'S','P','D','F','G','H','I','K','L','M','N', &
       'O','Q'/
  dimension G(2),MULT(*),DESIGN(*),TAM(*),JI(*),JF(*),DLP(*),DLL(*),&
       DLR(*),SP(*),SL(*),SR(*) 
  Character, Dimension(21) :: code1
  Data code1 /'p','1','f','2','h','3','k','4','m','5','o', &
       '6','r','7','t','8','u','9','v','0','w'/
  
  OAM=-1
  JI1=JI(1) 
  JI2=JI(2) 
  JF1=JF(1) 
  JF2=JF(2) 
  DO 58 I=1,MC  
     DLL(I)=0. 
     SL(I)=0.  
     DLR(I)=0. 
     SR(I)=0.  
58 END DO
  G(1)=0.   
  G(2)=0.   
  LEVEL=0   
  IF((JI1+JF1).EQ.0) LEVEL=1
  IF((JI2+JF2).EQ.0) LEVEL=-1   
  IF(LEVEL .LT. 0) GOTO 3
  IF(LEVEL .EQ. 0) GOTO 4
  IF(LEVEL .GT. 0) GOTO 5
3 I=1                    ! TRIPLETS WITH G=0
  GO TO 6   
5 I=2   
6 TAM(I)=1.              ! TOTAL ANGULAR MOMENTUM   
  SPIN=0.5*FLOAT(MULT(I)-1) 
  DO 7 J=1,13   
     IF(DESIGN(I).NE.CODE(J)) GO TO 7  
     OAM=J-1                ! ORBITAL ANGULAR MOMENTUM 
     GO TO 8   
7    Continue
  do 70 j=1,21   
     if(design(I).ne.code1(j)) goto 70  
     oam=float(j)/2.                ! ORBITAL ANGULAR MOMENTUM (semientero)
     go to 8 
70 Continue

  STOP 'EXIT ZEEMAN- error en el momento angular orbital'

8 G(I)=1.5+(SPIN*(1.+SPIN)-OAM*(1.+OAM))/4.
13 NP=1                   ! TRIPLETS WITH G=0 OR G1=G2   
  NL=NP 
  NR=NP 
  SP(1)=1.  
  SL(1)=1.  
  SR(1)=1.  
  DLP(1)=0. 
  DLL(1)=G(I)   
  DLR(1)=-DLL(1)
  GO TO 9   
4 DO I=1,2   
     SPIN=0.5*FLOAT(MULT(I)-1) 
     DO 11 J=1,13  
        IF(DESIGN(I).NE.CODE(J)) GO TO 11 
        OAM=J-1   
        GO TO 10
11   Continue
     do 110 j=1,21
        if(design(I).ne.code1(j)) goto 110  
        oam=float(j)/2.                ! ORBITAL ANGULAR MOMENTUM (semientero)
        go to 10
110  Continue

     PRINT *,'DESIGN=',DESIGN(I)
     STOP 'EXIT ZEEMAN 1'
!	LANDE FACTOR FOR EACH LABEL
10   G(I)=1.5+(SPIN*(1.+SPIN)-OAM*(1.+OAM))/(2.*TAM(I)*(1.+TAM(I)))
  End do
     IF(ABS(G(2)-G(1)).GT.5.E-6) GO TO 12  
     I=2   
     GO TO 13  
12   LEVEL=JI2-JI1 
     IF(JF1.EQ.5) GO TO 14 
     IF(LEVEL .LT. 0) GOTO 15     ! INTEGRAL J'S 
     IF(LEVEL .EQ. 0) GOTO 16     ! INTEGRAL J'S 
     IF(LEVEL .GT. 0) GOTO 17     ! INTEGRAL J'S 
15   NP=2*JI2+1
19   NL=NP 
     NR=NP 
     IF(NP.LE.MC) GO TO 18 
     PRINT *,'NP=',NP,' MC=',MC
     STOP 'EXIT ZEEMAN 2'
16   NP=2*JI2  
     GO TO 19  
17   NP=2*JI2-1
     GO TO 19  
18   MUMIN=-JI2
     IF(JI1.LT.JI2) MUMIN=-JI1 
     MUMAX=-MUMIN  
     I=0   
     DO 20 MU=MUMIN,MUMAX   ! PI COMPONENTS
        IF(MU.EQ.0.AND.LEVEL.EQ.0) GO TO 20   
        I=I+1 
        DLP(I)=FLOAT(MU)*(G(1)-G(2))  
        J=MU**2   
        IF(LEVEL .LT. 0) GOTO 21
        IF(LEVEL .EQ. 0) GOTO 22
        IF(LEVEL .GT. 0) GOTO 23
21      SP(I)=2*(JI1**2-J)
        GO TO 20  
22      SP(I)=2*J 
        GO TO 20  
23      SP(I)=2*(JI2**2-J)
20      CONTINUE  
      MUMIN=1-JI1   
      MUMAX=JI2 
      I=0   
      DO 24 MU=MUMIN,MUMAX   ! R-SIGMA COMPONENTS   
         I=I+1 
         DLR(I)=FLOAT(MU)*(G(1)-G(2))-G(1) 
         IF(LEVEL .LT. 0) GOTO 25
         IF(LEVEL .EQ. 0) GOTO 26
         IF(LEVEL .GT. 0) GOTO 27
25       SR(I)=(JI1-MU)*(JI2-MU+2) 
         GO TO 24  
26       SR(I)=(JI2+MU)*(JI2-MU+1) 
         GO TO 24  
27       SR(I)=(JI2+MU)*(JI1+MU)   
24    CONTINUE  
      MUMIN=-JI2
      MUMAX=JI1-1   
      I=0   
      DO 28 MU=MUMIN,MUMAX   ! L-SIGMA COMPONENTS   
         I=I+1 
         DLL(I)=FLOAT(MU)*(G(1)-G(2))+G(1) 
         IF(LEVEL .LT. 0) GOTO 29
         IF(LEVEL .EQ. 0) GOTO 30
         IF(LEVEL .GT. 0) GOTO 31
29       SL(I)=(JI1+MU)*(JI2+MU+2) 
         GO TO 28  
30       SL(I)=(JI2-MU)*(JI2+MU+1) 
         GO TO 28  
31       SL(I)=(JI2-MU)*(JI1-MU)   
28    CONTINUE  
      GO TO 57  
14    I=2*JI1+1              ! HALF-INTEGRAL J'S
      IF(LEVEL .LT. 0) GOTO 32
      IF(LEVEL .EQ. 0) GOTO 33
      IF(LEVEL .GT. 0) GOTO 34
32    NP=I-1
      NL=NP 
      NR=NP 
36    IF(NP.LE.MC) GO TO 35 
      STOP 'EXIT ZEEMAN'
33    NP=I+1
      NL=I  
      NR=I  
      GO TO 36  
34    NP=I+1
      NL=NP 
      NR=NP 
      GO TO 36  
35    MUMIN=-JI2
      IF(JI1.LT.JI2) MUMIN=-JI1 
      MUMAX=1-MUMIN 
      I=0   
      DO 37 MU=MUMIN,MUMAX   ! PI COMPONENTS
         I=I+1 
         SPIN=FLOAT(MU)-0.5
         DLP(I)=(G(1)-G(2))*SPIN   
         SPIN=SPIN**2  
         IF(LEVEL .LT. 0) GOTO 38
         IF(LEVEL .EQ. 0) GOTO 39
         IF(LEVEL .GT. 0) GOTO 40
38       SP(I)=2.*((FLOAT(JI1)+0.5)**2-SPIN)   
         GO TO 37  
39       SP(I)=2.*SPIN 
         GO TO 37  
40       SP(I)=2.*((FLOAT(JI2)+0.5)**2-SPIN)   
37    CONTINUE  
      MUMIN=-JI1
      MUMAX=JI2 
      I=0   
      DO 41 MU=MUMIN,MUMAX   ! R-SIGMA COMPONENTS   
         I=I+1 
         DLR(I)=(FLOAT(MU)+0.5)*(G(1)-G(2))-G(1)   
         IF(LEVEL .LT. 0) GOTO 42
         IF(LEVEL .EQ. 0) GOTO 43
         IF(LEVEL .GT. 0) GOTO 44
42       SR(I)=(JI1-MU)*(JI2-MU+2) 
         GO TO 41  
43       SR(I)=(JI2+MU+1)*(JI2-MU+1)   
         GO TO 41  
44       SR(I)=(JI2+MU+1)*(JI2+MU) 
41    CONTINUE  
      MUMIN=-MUMAX  
      MUMAX=JI1 
      I=0   
      DO 45 MU=MUMIN,MUMAX   ! L-SIGMA COMPONENTS   
         I=I+1 
         DLL(I)=(FLOAT(MU)-0.5)*(G(1)-G(2))+G(1)   
         IF(LEVEL .LT. 0) GOTO 46
         IF(LEVEL .EQ. 0) GOTO 47
         IF(LEVEL .GT. 0) GOTO 48
46       SL(I)=(JI1+MU)*(JI2+MU+2) 
         GO TO 45  
47       SL(I)=(JI2-MU+1)*(JI2+MU+1)   
         GO TO 45  
48       SL(I)=(JI2-MU+1)*(JI1-MU+1)   
45       CONTINUE  
57       SUM=0.
         DO 49 I=1,NP  
            SUM=SUM+SP(I) 
49    CONTINUE
      DO 50 I=1,NP  
         SP(I)=SP(I)/SUM   
50       CONTINUE
      SPIN=0.   
      SUM=0.
      DO 51 I=1,NL           ! NL=NR ASSUMED.   
         SPIN=SPIN+SL(I)   
         SUM=SUM+SR(I) 
51    CONTINUE
      DO 52 I=1,NL  
         SL(I)=SL(I)/SPIN  
         SR(I)=SR(I)/SUM   
52       CONTINUE
   9  CONTINUE  
 59   DO 54 I=1,NP  
55          FORMAT(2X,3(F10.6,F11.8))
60          DLP(I)=DLP(I)*DL0
            DLL(I)=DLL(I)*DL0 
            DLR(I)=DLR(I)*DL0 
54    CONTINUE
      IF(MC.LE.18) RETURN   
  56  FORMAT(3X,'Lande Factors: g(lower)=',F11.8,', g(upper)=',F11.8)
      RETURN
      END SUBROUTINE ZEEMAN


      


End Module Zeeman_Splitting
