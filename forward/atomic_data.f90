Module Atomic_data
!
  Parameter (N_elements = 92)
  Character (len=2), dimension (N_elements) :: Atom_char
  Real, dimension (N_elements) :: At_weight, At_ioniz1, At_ioniz2, &
       At_abund, At_Greveese_abund, Starting_at_abund
!
  Data Atom_char/'H','HE','LI','BE','B','C','N','O','F','NE', &
     'NA','MG','AL','SI','P','S','CL','AR','K','CA','SC','TI','V','CR', &
     'MN','FE','CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR', &
     'RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN', &
     'SN','SB','TE','I','XE','CS','BA','LA','CE','PR','ND','PM', &
     'SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA','W', &
     'RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN', &
     'FR','RA','AC','TH','PA','U'/
!
  Data At_weight/1.008,4.003,6.939,9.012,10.811,12.011,14.007,16.,18.998, &
     20.183,22.99,24.312,26.982,28.086,30.974,32.064,35.453,39.948,39.102 &
     ,40.08,44.956,47.90,50.942,51.996,54.938,55.847,58.933,58.71,63.54 &
     ,65.37,69.72,72.59,74.92,78.96,79.91,83.80,85.47,87.62,88.905,91.22 &
     ,92.906,95.94, 99.00,101.07,102.9,106.4,107.87,112.40,114.82, &
     118.69,121.75,127.6,126.9,131.3,132.9,137.34,138.91,140.12,140.91, &
     144.24,147.00,150.35,151.96,157.25,158.92,162.50,164.93, &
     167.26,168.93,173.04,174.97,178.49,180.95,183.85,186.2,190.2,192.2 &
     ,195.09,196.97,200.59,204.37,207.19,208.98,210.,211.,222., &
     223.,226.1,227.1,232.04,231.,238.03/
!
  Data At_ioniz1/13.595,24.58,5.39,9.32,8.298,11.256,14.529,13.614,17.418, &
     21.559,5.138,7.644,5.984,8.149,10.474,10.357,13.012,15.755,4.339, &
     6.111,6.538,6.825,6.738,6.763,7.432,7.896,7.863,7.633,7.724,9.391, &
     5.997,7.88,9.81,9.75,11.840,13.996,4.176,5.692,6.377,6.838,6.881, &
     7.10,7.28,7.36,7.46,8.33,7.574,8.991,5.785,7.34,8.64,9.01,10.454, &
     12.127,3.893,5.210,5.577,5.466,5.422,5.489,5.554,5.631,5.666,6.141 &
     ,5.852,5.927,6.018,6.101,6.184,6.254,5.426,6.650,7.879,7.980,7.870 &
     ,8.70,9.10,9.00,9.22,10.43,6.105,7.415,7.287,8.43,9.30,10.745, &
     4.,5.276,6.9,6.,6.,6./
  Data At_ioniz2/0.,54.403,75.62,18.21,25.15,24.376,29.59,35.11,34.98,41.07 &
     ,47.290,15.03,18.823,16.34,19.72,23.405,23.798,27.62,31.81,11.868, &
     12.891,13.63,14.205,16.493,15.636,16.178,17.052,18.15,20.286,17.96 &
     ,20.509,15.93,18.63,21.50,21.60,24.565,27.50,11.027,12.233,13.13, &
     14.316,16.15,15.26,16.76,18.07,19.42,21.48,16.904,18.86,14.63, &
     16.50,18.60,19.09,21.20,25.10,10.001,11.060,10.850,10.550,10.730, &
     10.899,11.069,11.241,12.090,11.519,11.670,11.800,11.930,12.050, &
     12.184,13.900,14.900,16.2,17.7,16.60,17.00,20.00,18.56,20.50,18.75 &
     ,20.42,15.03,16.68,19.,20.,20.,22.,10.144,12.1,12.,12.,12./
! Abundances from Greveese (1984) Physica Scripta. Vol.T8,49-58
   Data At_Greveese_abund/ &
             12.00,11.00,1.00,1.15,2.60,8.55,7.99,8.77,4.56,8.00,6.18, &
              7.48, 6.40,7.55,5.45,7.21,5.50,6.58,5.12,6.36,3.10,5.02, &
              4.00, 5.67,5.45,7.50,4.92,6.25,4.21,4.60,2.88,3.63,2.39, &
              3.35, 2.63,3.21,2.60,2.90,2.24,2.56,2.10,1.92,0.00,1.84, &
              1.12, 1.69,0.94,1.86,1.66,2.00,1.00,2.25,1.51,2.19,1.12, &
              2.13, 1.22,1.55,0.71,1.34,0.00,0.80,0.51,1.12,0.20,1.10, &
              0.26, 0.93,0.00,1.08,0.76,0.88,-.09,1.11,0.26,1.45,1.35, &
              1.80, 1.13,1.27,0.90,1.90,0.71,-8.0,-8.0,-8.0,-8.0,-8.0, &
              -8.00, 0.02,-8.0 ,-.47/
! Abundances from Greveese (1984) Physica Scripta. Vol.T8,49-58
   Data At_abund/ &
             12.00,11.00,1.00,1.15,2.60,8.55,7.99,8.77,4.56,8.00,6.18, &
              7.48, 6.40,7.55,5.45,7.21,5.50,6.58,5.12,6.36,3.10,5.02, &
              4.00, 5.67,5.45,7.50,4.92,6.25,4.21,4.60,2.88,3.63,2.39, &
              3.35, 2.63,3.21,2.60,2.90,2.24,2.56,2.10,1.92,0.00,1.84, &
              1.12, 1.69,0.94,1.86,1.66,2.00,1.00,2.25,1.51,2.19,1.12, &
              2.13, 1.22,1.55,0.71,1.34,0.00,0.80,0.51,1.12,0.20,1.10, &
              0.26, 0.93,0.00,1.08,0.76,0.88,-.09,1.11,0.26,1.45,1.35, &
              1.80, 1.13,1.27,0.90,1.90,0.71,-8.0,-8.0,-8.0,-8.0,-8.0, &
              -8.00, 0.02,-8.0 ,-.47/
   Public Partition_f


 Contains

     Subroutine Partition_f(NEL, T, U1, U2, U3, DU1, DU2, DU3)
! This subroutine computes partition functions. NEL is element number and
! T is temperature.
       du1=0.
       du2=0.
       du3=0.
       X=ALOG(5040./T)
       dx=-1./t
       Y=1.E-3*T
       If (NEL .eq. 1) then
          Goto 1
       Else If (NEL .eq. 2) then
          Goto 2
       Else If (NEL .eq. 3) then
          Goto 3
       Else If (NEL .eq. 4) then
          Goto 4
       Else If (NEL .eq. 5) then
          Goto 5
       Else If (NEL .eq. 6) then
          Goto 6
       Else If (NEL .eq. 7) then
          Goto 7
       Else If (NEL .eq. 8) then
          Goto 8
       Else If (NEL .eq. 9) then
          Goto 9
       Else If (NEL .eq. 10) then
          Goto 10
       Else If (NEL .eq. 11) then
          Goto 11
       Else If (NEL .eq. 12) then
          Goto 12
       Else If (NEL .eq. 13) then
          Goto 13
       Else If (NEL .eq. 14) then
          Goto 14
       Else If (NEL .eq. 15) then
          Goto 15
       Else If (NEL .eq. 16) then
          Goto 16
       Else If (NEL .eq. 17) then
          Goto 17
       Else If (NEL .eq. 18) then
          Goto 18
       Else If (NEL .eq. 19) then
          Goto 19
       Else If (NEL .eq. 20) then
          Goto 20
       Else If (NEL .eq. 21) then
          Goto 21
       Else If (NEL .eq. 22) then
          Goto 22
       Else If (NEL .eq. 23) then
          Goto 23
       Else If (NEL .eq. 24) then
          Goto 24
       Else If (NEL .eq. 25) then
          Goto 25
       Else If (NEL .eq. 26) then
          Goto 26
       Else If (NEL .eq. 27) then
          Goto 27
       Else If (NEL .eq. 28) then
          Goto 28
       Else If (NEL .eq. 29) then
          Goto 29
       Else If (NEL .eq. 30) then
          Goto 30
       Else If (NEL .eq. 31) then
          Goto 31
       Else If (NEL .eq. 32) then
          Goto 32
       Else If (NEL .eq. 33) then
          Goto 33
       Else If (NEL .eq. 34) then
          Goto 34
       Else If (NEL .eq. 35) then
          Goto 35
       Else If (NEL .eq. 36) then
          Goto 36
       Else If (NEL .eq. 37) then
          Goto 37
       Else If (NEL .eq. 38) then
          Goto 38
       Else If (NEL .eq. 39) then
          Goto 39
       Else If (NEL .eq. 40) then
          Goto 40
       Else If (NEL .eq. 41) then
          Goto 41
       Else If (NEL .eq. 42) then
          Goto 42
       Else If (NEL .eq. 43) then
          Goto 43
       Else If (NEL .eq. 44) then
          Goto 44
       Else If (NEL .eq. 45) then
          Goto 45
       Else If (NEL .eq. 46) then
          Goto 46
       Else If (NEL .eq. 47) then
          Goto 47
       Else If (NEL .eq. 48) then
          Goto 48
       Else If (NEL .eq. 49) then
          Goto 49
       Else If (NEL .eq. 50) then
          Goto 50
       Else If (NEL .eq. 51) then
          Goto 51
       Else If (NEL .eq. 52) then
          Goto 52
       Else If (NEL .eq. 53) then
          Goto 53
       Else If (NEL .eq. 54) then
          Goto 54
       Else If (NEL .eq. 55) then
          Goto 55
       Else If (NEL .eq. 56) then
          Goto 56
       Else If (NEL .eq. 57) then
          Goto 57
       Else If (NEL .eq. 58) then
          Goto 58
       Else If (NEL .eq. 59) then
          Goto 59
       Else If (NEL .eq. 60) then
          Goto 60
       Else If (NEL .eq. 61) then
          Goto 61
       Else If (NEL .eq. 62) then
          Goto 62
       Else If (NEL .eq. 63) then
          Goto 63
       Else If (NEL .eq. 64) then
          Goto 64
       Else If (NEL .eq. 65) then
          Goto 65
       Else If (NEL .eq. 66) then
          Goto 66
       Else If (NEL .eq. 67) then
          Goto 67
       Else If (NEL .eq. 68) then
          Goto 68
       Else If (NEL .eq. 69) then
          Goto 69
       Else If (NEL .eq. 70) then
          Goto 70
       Else If (NEL .eq. 71) then
          Goto 71
       Else If (NEL .eq. 72) then
          Goto 72
       Else If (NEL .eq. 73) then
          Goto 73
       Else If (NEL .eq. 74) then
          Goto 74
       Else If (NEL .eq. 75) then
          Goto 75
       Else If (NEL .eq. 76) then
          Goto 76
       Else If (NEL .eq. 77) then
          Goto 77
       Else If (NEL .eq. 78) then
          Goto 78
       Else If (NEL .eq. 79) then
          Goto 79
       Else If (NEL .eq. 80) then
          Goto 80
       Else If (NEL .eq. 81) then
          Goto 81
       Else If (NEL .eq. 82) then
          Goto 82
       Else If (NEL .eq. 83) then
          Goto 83
       Else If (NEL .eq. 84) then
          Goto 84
       Else If (NEL .eq. 85) then
          Goto 85
       Else If (NEL .eq. 86) then
          Goto 86
       Else If (NEL .eq. 87) then
          Goto 87
       Else If (NEL .eq. 88) then
          Goto 88
       Else If (NEL .eq. 89) then
          Goto 89
       Else If (NEL .eq. 90) then
          Goto 90
       Else If (NEL .eq. 91) then
          Goto 91
       Else If (NEL .eq. 92) then
          Goto 92
       Endif
1     U1=2.00 ! H
       if(T.GT.1.3E4)then
          U1=1.51+3.8E-5*T
          du1=3.8e-5/u1
       end if
       if(T.GT.1.62E4)then
          U1=11.41+T*(-1.1428E-3+T*3.52E-8)
          du1=(-1.1428e-3+2.*t*3.52e-8)/u1
       end if
      U2=1.000000
      U3=0.000000
      RETURN
2     U1=1.000 ! HE
      if(T.GT.3E4)then
         u1=14.8+T*(-9.4103E-4+T*1.6095E-8)
         du1=(-9.4103e-4+2*t*1.6095e-8)/u1
      end if
      U2=2.000
      U3=1.000000
      RETURN
3     U1=2.081-Y*(6.8926E-2-Y*1.4081E-2) ! LI
      du1=1.e-3*(-6.8926e-2+2*y*1.4081e-2)/u1
      if(T.GT.6E3)then
         u1=3.4864+T*(-7.3292E-4+T*8.5586E-8)
         du1=(-7.3292e-4+2*t*8.5586e-8)/u1
      end if
      U2=1.0000
      U3=2.00000
      RETURN
4     U1=AMAX1(1.,.631+7.032E-5*T) ! BE
      if(u1.ne.1.)du1=7.032e-5/u1
      U2=2.00
      U3=1.000
      RETURN
5     U1=5.9351+1.0438E-2*Y ! B
      du1=1.e-3*1.0438e-2/u1
      U2=1.00
      U3=2.0
      RETURN
6     U1=8.6985+Y*(2.0485E-2+Y*(1.7629E-2-3.9091E-4*Y)) ! C
      du1=1.e-3*((2.0485e-2+2*y*1.7629e-2-3*y*y*3.9091e-4)/u1)
      if(T.GT.1.2E4)then
         U1=13.97+T*(-1.3907E-3+T*9.0844E-8)
         du1=(-1.3907e-3+2*t*9.0844e-8)/u1
      end if
      U2=5.838+1.6833E-5*T
      du2=1.6833e-5/u2
      if(T.GT.2.4E4)then
         u2=10.989+T*(-6.9347E-4+T*2.0861E-8)
         du2=(-6.9347e-4+2*t*2.0861e-8)/u2
      end if
      U3=1.00
      if(T.GT.1.95E4)then
         u3=-0.555+8E-5*T
         du3=8e-5/u3
      end if
      RETURN
7     U1=3.9914+Y*(1.7491E-2-Y*(1.0148E-2-Y*1.7138E-3)) ! N
      du1=1.e-3*((1.7491e-2-2*y*1.0148e-2+3*y*y*1.7138e-3)/u1)
      if(T.GT.8800.)then
         u1=2.171+2.54E-4*T
         du1=2.54e-4/u1
      end if
      if(T.GT.1.8E4)then
         u1=11.396+T*(-1.7139E-3+T*8.633E-8)
         du1=(-1.7139e-3+2*t*8.633e-8)/u1
      end if
      U2=8.060+1.420E-4*T
      du2=1.420e-4/u2
      if(T.GT.3.3E4)then
         u2=26.793+T*(-1.8931E-3+T*4.4612E-8)
         du2=(-1.8931e-3+2*t*4.4612e-8)/u2
      end if
      U3=5.9835+T*(-2.6651E-5+T*1.8228E-9)
      du3=(-2.6651e-5+2*t*1.8228e-9)/u3
      if(T.LT.7310.5)then
         U3=5.89
         du3=0.
      end if
      RETURN
8     U1=8.29+1.10E-4*T ! O
      du1=(1.10e-4)/u1
      if(T.GT.1.9E4)then
         U1=66.81+T*(-6.019E-3+T*1.657E-7)
         du1=(-6.019e-3+2*t*1.657e-7)/u1
      end if
      U2=AMAX1(4.,3.51+8.E-5*T)
      if(u2.ne.4.)du2=8.e-5/u2
      if(T.GT.3.64E4)then
         u2=68.7+T*(-4.216E-3+T*6.885E-8)
         du2=(-4.216e-3+2*t*6.885e-8)/u2
      end if
      U3=7.865+1.1348E-4*T
      du3=1.1348e-4/u3
      RETURN
9     U1=4.5832+Y*(.77683+Y*(-.20884+Y*(2.6771E-2-1.3035E-3*Y))) ! F
      du1=.77683-2*y*.20884+3*y*y*2.6771e-2-4*y*y*y*1.3035e-3
      du1=1.e-3*(du1/u1)
      if(T.GT.8750.)then
         u1=5.9
         du1=0.
      end if
      if(T.GT.2E4)then
         u1=15.16+T*(-9.229E-4+T*2.312E-8)
         du1=(-9.229e-4+2*t*2.312e-8)/u1
      end if
      U2=8.15+8.9E-5*T
      du2=8.9e-5/u2
      U3=2.315+1.38E-4*T
      du3=1.38e-4/u3
      RETURN
10    U1=1.0000 ! NE
      if(T.GT.2.69E4)then
         U1=26.3+T*(-2.113E-3+T*4.359E-8)
         du1=(-2.113e-3+2*t*4.359e-8)/u1
      end if
      U2=5.4+4.E-5*T
      du2=4.e-5/u2
      U3=7.973+7.956E-5*T
      du3=7.956e-5/u3
      RETURN
11    U1=AMAX1(2.,1.72+9.3E-5*T) ! NA
      if(u1.ne.2.)du1=9.3e-5/u1
      if(T.GT.5400.)then
         U1=-0.83+5.66E-4*T
         du1=5.66e-4/u1
      end if
      if(T.GT.8.5E3)then
         U1=4.5568+T*(-1.2415E-3+T*1.3861E-7)
         du1=(-1.2415e-3+2*t*1.3861e-7)/u1
      end if
      U2=1.000
      U3=5.69+5.69E-6*T
      du3=5.69e-6/u3
      RETURN
12    U1=1.+EXP(-4.027262-X*(6.173172+X*(2.889176+X*(2.393895+.784131*X) &! MG
           )))
      du1=-dx*(6.173172+x*(2*2.889176+x*(3*2.393895+4*x*.784131)))
      du1=du1*(u1-1.)/u1
      if(T.GT.8E3)then
         u1=2.757+T*(-7.8909E-4+T*7.4531E-8)
         du1=(-7.8909e-4+2*t*7.4531e-8)/u1
      end if
      U2=2.+EXP(-7.721172-X*(7.600678+X*(1.966097+.212417*X)))
      du2=-dx*(7.0600678+x*(2*1.966097+3*x*.212417))
      du2=du2*(u2-2.)/u2
      if(T.GT.2E4)then
         U2=7.1041+T*(-1.0817E-3+T*4.7841E-8)
         du2=(-1.0817e-3+2*t*4.7841e-8)/u2
      end if
      U3=1.0000
      RETURN
13    U1=5.2955+Y*(.27833-Y*(4.7529E-2-Y*3.0199E-3)) ! AL
      du1=1.e-3*(.27833-y*(2*4.7529e-2-3*y*3.0199e-3))/u1
      U2=AMAX1(1.,.725+3.245E-5*T)
      if(u2.ne.1.)du2=3.245e-5/u2
      if(T.GT.2.24E4)then
         U2=61.06+T*(-5.987E-3+T*1.485E-7)
         du2=(-5.987e-3+2*t*1.485e-7)/u2
      end if
      U3=AMAX1(2.,1.976+3.43E-6*T)
      if(u3.ne.2.)du3=3.43e-6/u3
      if(T.GT.1.814E4)then
         U3=3.522+T*(-1.59E-4+T*4.382E-9)
         du3=(-1.59e-4+2*t*4.382e-9)/u3
      end if
      RETURN
14    U1=6.7868+Y*(.86319+Y*(-.11622+Y*(.013109-6.2013E-4*Y))) ! SI
      du1=1.e-3*(.86319+y*(-2*.11622+y*(3*.013109-4*y*6.2013e-4)))/u1
      if(T.GT.1.04E4)then
         U1=86.01+T*(-1.465E-2+T*7.282E-7)
         du1=(-1.465e-2+2*t*7.282e-7)/u1
      end if
      U2=5.470+4.E-5*T
      du2=4.e-5/u2
      if(T.GT.1.8E4)then
         U2=26.44+T*(-2.22E-3+T*6.188E-8)
         du2=(-2.22e-3+2*t*6.188e-8)/u2
      end if
      U3=AMAX1(1.,.911+1.1E-5*T)
      if(u3.ne.1.)du3=1.1e-5/u3
      if(T.GT.3.33E4)then
         U3=19.14+T*(-1.408E-3+T*2.617E-8)
         du3=(-1.408e-3+2*t*2.617e-8)/u3
      end if
      RETURN
15    U1=4.2251+Y*(-.22476+Y*(.057306-Y*1.0381E-3)) ! P
      du1=1.e-3*(-.22476+y*(2*.057306-3*y*1.0381e-3))/u1
      if(T.GT.6.E3)then
         U1=1.56+5.2E-4*T
         du1=5.2e-4/u1
      end if
      U2=4.4151+Y*(2.2494+Y*(-.55371+Y*(.071913-Y*3.5156E-3)))
      du2=1.e-3*(2.2494+y*(-2*.55371+y*(3*.071913-4*y*3.5156e-3)))/u2
      if(T.GT.7250.)then
         U2=4.62+5.38E-4*T
         du2=5.38e-4/u2
      end if
      U3=5.595+3.4E-5*T
      du3=3.4e-5/u3
      RETURN
16    U1=7.5+2.15E-4*T ! S
      du1=2.1e-4/u1
      if(T.GT.1.16E4)then
         U1=38.76+T*(-4.906E-3+T*2.125E-7)
         du1=(-4.906e-3+2*t*2.125e-7)/u1
      end if
      U2=2.845+2.43E-4*T
      du2=2.43e-4/u2
      if(T.GT.1.05E4)then
         U2=6.406+T*(-1.68E-4+T*1.323E-8)
         du2=(-1.68e-4+t*2*1.323e-8)/u2
      end if
      U3=7.38+1.88E-4*T
      du3=1.88e-4/u3
      RETURN
17    U1=5.2+6.E-5*T ! CL
      du1=6.e-5/u1
      if(T.GT.1.84E4)then
         U1=-81.6+4.8E-3*T
         du1=4.8e-3/u1
      end if
      U2=7.0+2.43E-4*T
      du2=2.43e-4/u2
      U3=2.2+2.62E-4*T
      du3=2.62e-4/u3
      RETURN
18    U1=1.000 ! AR
      U2=5.20+3.8E-5*T
      du2=3.8e-5/u2
      U3=7.474+1.554E-4*T
      du3=1.554e-4/u3
      RETURN
19    U1=1.9909+Y*(.023169-Y*(.017432-Y*4.0938E-3)) ! K
      du1=1.e-3*(.023169-Y*(2*.017432-y*3*4.0938e-3))/u1
      if(T.GT.5800.)then
         U1=-9.93+2.124E-3*T
         du1=2.124e-3/u1
      end if
      U2=1.000
      U3=5.304+1.93E-5*T
      du3=1.93e-5/u3
      RETURN
20    U1=1.+EXP(-1.731273-X*(5.004556+X*(1.645456+X*(1.326861+.508553*X) &! CA
           )))
      du1=-dx*(5.004556+x*(2*1.645456+x*(3*1.326861+4*.508553*x)))
      du1=du1*(u1-1.)/u1
      U2=2.+EXP(-1.582112-X*(3.996089+X*(1.890737+.539672*X)))
      du2=-dx*(3.996089+x*(2*1.890737+x*3*.539672))
      du2=du2*(u2-2.)/u2
      U3=1.000
      RETURN
21    U1=4.+EXP(2.071563+X*(-1.2392+X*(1.173504+.517796*X))) ! SC
      du1=dx*(-1.2392+x*(2*1.173504+3*x*.517796))
      du1=du1*(u1-4.)/u1
      U2=3.+EXP(2.988362+X*(-.596238+.054658*X))
      du2=dx*(-.596238+2*x*.054658)
      du2=du2*(u2-3.)/u2
      U3=10. ! APPROXIMATELY
      RETURN
22    U1=5.+EXP(3.200453+X*(-1.227798+X*(.799613+.278963*X))) ! TI
      du1=dx*(-1.227798+x*(2*.799613+3*x*.278963))
      du1=du1*(u1-5.)/u1
      if(T.LT.5.5E3)then
         U1=16.37+T*(-2.838E-4+T*5.819E-7)
         du1=(-2.838e-4+2*t*5.819e-7)/u1
      end if
      U2=4.+EXP(3.94529+X*(-.551431+.115693*X))
      du2=dx*(-.551431+2.*x*.115693)
      du2=du2*(u2-4.)/u2
      U3=16.4+8.5E-4*T
      du3=8.5e-4/u3
      RETURN
23    U1=4.+EXP(3.769611+X*(-.906352+X*(.724694+.1622*X))) ! V
      du1=dx*(-.906352+x*(2*.724694+3.*x*.1622))
      du1=du1*(u1-4.)/u1
      U2=1.+EXP(3.755917+X*(-.757371+.21043*X))
      du2=dx*(-.757371+x*2.*.21043)
      du2=du2*(u2-1.)/u2
      U3=-18.+1.03E-2*T
      du3=1.03e-2/u3
      if(T.LT.2.25E3)then
         U3=2.4E-3*T
         du3=2.4e-3/u3
      end if
      RETURN
24    U1=7.+EXP(1.225042+X*(-2.923459+X*(.154709+.09527*X))) ! CR
      du1=dx*(-2.923459+x*(2*.154709+x*3*.09527))
      du1=du1*(u1-7.)/u1
      U2=6.+EXP(.128752-X*(4.143973+X*(1.096548+.230073*X)))
      du2=-dx*(4.143973+x*(2*1.096548+3.*x*.230073))
      du2=du2*(u2-6.)/u2
      U3=10.4+2.1E-3*T
      du3=2.1e-3/u3
      RETURN
25    U1=6.+EXP(-.86963-X*(5.531252+X*(2.13632+X*(1.061055+.265557*X)))) ! MN
      du1=-dx*(5.531252+x*(2*2.13632+x*(3*1.061055+4*.265557*x)))
      du1=du1*(u1-6.)/u1
      U2=7.+EXP(-.282961-X*(3.77279+X*(.814675+.159822*X)))
      du2=-dx*(3.77279+x*(2*.814675+3*x*.159822))
      du2=du2*(u2-7.)/u2
      U3=10. ! APPROXIMATELY
      RETURN
26    U1=9.+EXP(2.930047+X*(-.979745+X*(.76027+.118218*X))) ! FE
      du1=dx*(-.979745+x*(2*.76027+3*x*.118218))
      du1=du1*(u1-9.)/u1
      if(T.LT.4E3)then
         U1=15.85+T*(1.306E-3+T*2.04E-7)
         du1=(1.306e-3+2*t*2.04e-7)/u1
      end if
      if(T.GT.9E3)then
         U1=39.149+T*(-9.5922E-3+T*1.2477E-6)
         du1=(-9.5922e-3+2*t*1.2477e-6)/u1
      end if
      U2=10.+EXP(3.501597+X*(-.612094+.280982*X))
      du2=dx*(-.612094+2*x*.280982)
      du2=du2*(u2-10.)/u2
      if(T.GT.1.8E4)then
         U2=68.356+T*(-6.1104E-3+T*5.1567E-7)
         du2=(-6.1104e-3+2*t*5.1567e-7)/u2
      end if
      U3=17.336+T*(5.5048E-4+T*5.7514E-8)
      du3=(5.5048e-4+2*t*5.7514e-8)/u3
      RETURN
27    U1=8.65+4.9E-3*T ! CO
      du1=4.9e-3/u1
      U2=11.2+3.58E-3*T
      du2=3.58e-3/u2
      U3=15.0+1.42E-3*T
      du3=1.42e-3/u3
      RETURN
28    U1=9.+EXP(3.084552+X*(-.401323+X*(.077498-.278468*X))) ! NI
      du1=dx*(-.401323+x*(2*.077498-3.*x*.278468))
      du1=du1*(u1-9.)/u1
      U2=6.+EXP(1.593047-X*(1.528966+.115654*X))
      du2=-dx*(1.528966+2.*x*.115654)
      du2=du2*(u2-6.)/u2
      U3=13.3+6.9E-4*T
      du3=6.9e-4/u3
      RETURN
29    U1=AMAX1(2.,1.50+1.51E-4*T) ! CU
      if(u1.ne.2.)du1=1.51e-4/u1
      IF(T.GT.6250.)then
         U1=-.3+4.58E-4*T
         du1=4.58e-4/u1
      end if
      U2=AMAX1(1.,.22+1.49E-4*T)
      if(u2.ne.1.)du2=1.49e-4/u2
      U3=8.025+9.4E-5*T
      du3=9.4e-5/u3
      RETURN
30    U1=AMAX1(1.,.632+5.11E-5*T) ! ZN
      if(u1.ne.1.)du1=5.11e-5/u1
      U2=2.00
      U3=1.00
      RETURN
31    U1=1.7931+Y*(1.9338+Y*(-.4643+Y*(.054876-Y*2.5054E-3))) ! GA
      du1=1.e-3*(1.9338+y*(-2.*.4643+y*(3*.054876-y*4*2.5054e-3)))/u1
      if(T.GT.6.E3)then
         U1=4.18+2.03E-4*T
         du1=2.03e-4/u1
      end if
      U2=1.0
      U3=2.0
      RETURN
32    U1=6.12+4.08E-4*T ! GE
      du1=4.08e-4/u1
      U2=3.445+1.78E-4*T
      du2=1.78e-4/u2
      U3=1.1 ! APPROXIMATELY
      RETURN
33    U1=2.65+3.65E-4*T ! AS
      du1=3.65e-4/u1
      U2=-.25384+Y*(2.284+Y*(-.33383+Y*(.030408-Y*1.1609E-3)))
      du2=1.e-3*(2.284+y*(-2*.33383+y*(3*.030408-4*y*1.1609e-3)))/u2
      if(T.GT.1.2E4)then
         U2=8.
         du2=0.
      end if
      U3=8. ! APPROXIMATELY
      RETURN
34    U1=6.34+1.71E-4*T ! SE
      du1=1.71e-4/u1
      U2=4.1786+Y*(-.15392+3.2053E-2*Y)
      du2=1.e-3*(-.15392+3.2053e-2*y*2)/u2
      U3=8. ! APPROXIMATELY
      RETURN
35    U1=4.12+1.12E-4*T ! BR
      du1=1.12e-4/u1
      U2=5.22+3.08E-4*T
      du2=3.08e-4/u2
      U3=2.3+2.86E-4*T
      du3=2.86e-4/u3
      RETURN
36    U1=1.00 ! KR
      U2=4.11+7.4E-5*T
      du2=7.4e-5/u2
      U3=5.35+2.23E-4*T
      du3=2.23e-4/u3
      RETURN
37    U1=AMAX1(2.,1.38+1.94E-4*T) ! RB
      if(u1.ne.2.)du1=1.94e-4/u1
      if(T.GT.6250.)then
         U1=-14.9+2.79E-3*T
         du1=2.79e-3/u1
      end if
      U2=1.000
      U3=4.207+4.85E-5*T
      du3=4.85e-5/u3
      RETURN
38    U1=.87127+Y*(.20148+Y*(-.10746+Y*(.021424-Y*1.0231E-3))) ! SR
      du1=1.e-3*(.20148+y*(-.10746*2.+y*(.021424*3.-y*1.0231e-3*4.)))
      du1=du1/u1
      if(T.GT.6500.)then
         U1=-6.12+1.224E-3*T
         du1=1.224e-3/u1
      end if
      u2=amax1(2.,.84+2.6e-4*t)
      if(u2.ne.2.)du2=2.6e-4/u2
      U3=1.0 ! APPROXIMATELY
      RETURN
39    U1=.2+2.58E-3*T ! Y
      U2=7.15+1.855E-3*T
      U3=9.71+9.9E-5*T
      RETURN
40    U1=76.31+T*(-1.866E-2+T*2.199E-6)  ! ZR
      IF(T.LT.6236.) U1=6.8+T*(2.806E-3+T*5.386E-7)
      U2=4.+EXP(3.721329-.906502*X)
      U3=12.3+1.385E-3*T
      RETURN
41    U1=AMAX1(1.,-19.+1.43E-2*T) ! NB
      U2=-4.+1.015E-2*T
      U3=25. ! APPROXIMATELY
      RETURN
42    U1=AMAX1(7.,2.1+1.5E-3*T) ! MO
      IF(T.GT.7.E3) U1=-38.1+7.28E-3*T
      U2=1.25+1.17E-3*T
      IF(T.GT.6900.) U2=-28.5+5.48E-3*T
      U3=24.04+1.464E-4*T
      RETURN
43    U1=4.439+Y*(.30648+Y*(1.6525+Y*(-.4078+Y*(.048401-Y*2.1538E-3)))) ! TC
      IF(T.GT.6.E3) U1=24.
      U2=8.1096+Y*(-2.963+Y*(2.369+Y*(-.502+Y*(.049656-Y*1.9087E-3))))
      IF(T.GT.6.E3) U2=17.
      U3=220. ! APPROXIMATELY
      RETURN
44    U1=-3.+7.17E-3*T ! RU
      U2=3.+4.26E-3*T
      U3=22.
      RETURN
45    U1=6.9164+Y*(3.8468+Y*(.043125-Y*(8.7907E-3-Y*5.9589E-4))) ! RH
      U2=7.2902+Y*(1.7476+Y*(-.038257+Y*(2.014E-3+Y*2.1218E-4)))
      U3=30. ! APPROXIMATELY
      RETURN
46    U1=AMAX1(1.,-1.75+9.86E-4*T) ! PD
      U2=5.60+3.62E-4*T
      U3=20. ! APPROXIMATELY
      RETURN
47    U1=AMAX1(2.,1.537+7.88E-5*T) ! AG
      U2=AMAX1(1.,0.73+3.4E-5*T)
      U3=6.773+1.248E-4*T
      RETURN
48    U1=AMAX1(1.,.43+7.6E-5*T) ! CD
      U2=2.00
      U3=1.00
      RETURN
49    U1=2.16+3.92E-4*T ! IN
      U2=1.0
      U3=2.00
      RETURN
50    U1=2.14+6.16E-4*T ! SN
      U2=2.06+2.27E-4*T
      U3=1.05 ! APPROXIMATELY
      RETURN
51    U1=2.34+4.86E-4*T ! SB
      U2=.69+5.36E-4*T
      U3=3.5 ! APPROXIMATELY
      RETURN
52    U1=3.948+4.56E-4*T ! TE
      U2=4.2555+Y*(-.25894+Y*(.06939-Y*2.4271E-3))
      IF(T.GT.1.2E4) U2=7.
      U3=5. ! APPROXIMATELY
      RETURN
53    U1=AMAX1(4.,3.8+9.5E-5*T) ! I
      U2=4.12+3.E-4*T
      U3=7. ! APPROXIMATELY
      RETURN
54    U1=1.00 ! XE
      U2=3.75+6.876E-5*T
      U3=4.121+2.323E-4*T
      RETURN
55    U1=AMAX1(2.,1.56+1.67E-4*T) ! CS
      IF(T.GT.4850.) U1=-2.680+1.04E-3*T
      U2=1.000
      U3=3.769+4.971E-5*T
      RETURN
56    U1=AMAX1(1.,-1.8+9.85E-4*T) ! BA
      if(u1.ne.1.)du1=9.85e-4/u1
      if(T.GT.6850.)then
         U1=-16.2+3.08E-3*T
         du1=3.08e-3/u1
      end if
      U2=1.11+5.94E-4*T
      du2=5.94e-4/u2
      U3=1.00
      RETURN
57    U1=15.42+9.5E-4*T ! LA
      IF(T.GT.5060.) U1=1.+3.8E-3*T
      U2=13.2+3.56E-3*T
      U3=12. ! APPROXIMATELY
      RETURN
58    U1=9.+EXP(5.202903+X*(-1.98399+X*(.119673+.179675*X))) ! CE
      U2=8.+EXP(5.634882-X*(1.459196+X*(.310515+.052221*X)))
      U3=9.+EXP(3.629123-X*(1.340945+X*(.372409+X*(.03186-.014676*X))))
      RETURN
59    U2=9.+EXP(4.32396-X*(1.191467+X*(.149498+.028999*X))) ! PR
      U1=U2 ! APPROXIMATELY
      U3=10.+EXP(3.206855+X*(-1.614554+X*(.489574+.277916*X)))
      RETURN
60    U1=9.+EXP(4.456882+X*(-2.779176+X*(.082258+X*(.50666+.127326*X)))) ! ND
      U2=8.+EXP(4.689643+X*(-2.039946+X*(.17193+X*(.26392+.038225*X))))
      U3=U2 ! APPROXIMATELY
      RETURN
61    U1=20. ! PM APPROXIMATELY
      U2=25. ! APPROXIMATELY
      U3=100. ! APPROXIMATELY
      RETURN
62    U1=1.+EXP(3.549595+X*(-1.851549+X*(.9964+.566263*X))) ! SM
      U2=2.+EXP(4.052404+X*(-1.418222+X*(.358695+.161944*X)))
      U3=1.+EXP(3.222807-X*(.699473+X*(-.056205+X*(.533833+.251011*X))))
      RETURN
63    U1=8.+EXP(1.024374-X*(4.533653+X*(1.540805+X*(.827789+.286737*X)))) ! EU
      U2=9.+EXP(1.92776+X*(-1.50646+X*(.379584+.05684*X)))
      U3=8. ! APPROXIMATELY
      RETURN
64    U1=5.+EXP(4.009587+X*(-1.583513+X*(.800411+.388845*X))) ! GD
      U2=6.+EXP(4.362107-X*(1.208124+X*(-.074813+X*(.076453+.055475*X))))
      U3=5.+EXP(3.412951-X*(.50271+X*(.042489-4.017E-3*X)))
      RETURN
65    U1=16.+EXP(4.791661+X*(-1.249355+X*(.570094+.240203*X))) ! TB
      U2=15.+EXP(4.472549-X*(.295965+X*(5.88E-3+.131631*X)))
      U3=U2 ! APPROXIMATELY
      RETURN
66    U1=17.+EXP(3.029646-X*(3.121036+X*(.086671-.216214*X))) ! DY
      U2=18.+EXP(3.465323-X*(1.27062+X*(-.382265+X*(.431447+.303575*X))))
      U3=U2 ! APPROXIMATELY
      RETURN
67    U3=16.+EXP(1.610084-X*(2.373926+X*(.133139-.071196*X))) ! HO
      U1=U3
      U2=U3 ! APPROX.
      RETURN
68    U1=13.+EXP(2.895648-X*(2.968603+X*(.561515+X*(.215267+.095813*X)))) ! ER
      U2=14.+EXP(3.202542-X*(.852209+X*(-.226622+X*(.343738+.186042*X))))
      U3=U2 ! APPROX.
      RETURN
69    U1=8.+EXP(1.021172-X*(4.94757+X*(1.081603+.034811*X))) ! TM
      U2=9.+EXP(2.173152+X*(-1.295327+X*(1.940395+.813303*X)))
      U3=8.+EXP(-.567398+X*(-3.383369+X*(.799911+.554397*X)))
      RETURN
70    U1=1.+EXP(-2.350549-X*(6.688837+X*(1.93869+.269237*X))) ! YB
      U2=2.+EXP(-3.047465-X*(7.390444+X*(2.355267+.44757*X)))
      U3=1.+EXP(-6.192056-X*(10.560552+X*(4.579385+.940171*X)))
      RETURN
71    U1=4.+EXP(1.537094+X*(-1.140264+X*(.608536+.193362*X))) ! LU
      U2=AMAX1(1.,0.66+1.52E-4*T)
      IF(T.GT.5250.) U2=-1.09+4.86E-4*T
      U3=5. ! APPROXIMATELY
      RETURN
72    U1=4.1758+Y*(.407+Y*(.57862-Y*(.072887-Y*3.6848E-3))) ! HF
      U2=-2.979+3.095E-3*T
      U3=30. ! APPROXIMATELY
      RETURN
73    U1=3.0679+Y*(.81776+Y*(.34936+Y*(7.4861E-3+Y*3.0739E-4))) ! TA
      U2=1.6834+Y*(2.0103+Y*(.56443-Y*(.031036-Y*8.9565E-4)))
      U3=15.
      RETURN
74    U1=.3951+Y*(-.25057+Y*(1.4433+Y*(-.34373+Y*(.041924-Y*1.84E-3)))) ! W
      IF(T.GT.1.2E4) U1=23.
      U2=1.055+Y*(1.0396+Y*(.3303-Y*(8.4971E-3-Y*5.5794E-4)))
      U3=20.
      RETURN
75    U1=5.5671+Y*(.72721+Y*(-.42096+Y*(.09075-Y*3.9331E-3))) ! RE
      IF(T.GT.1.2E4) U1=29.
      U2=6.5699+Y*(.59999+Y*(-.28532+Y*(.050724-Y*1.8544E-3)))
      IF(T.GT.1.2E4) U2=22.
      U3=20.
      RETURN
76    U1=8.6643+Y*(-.32516+Y*(.68181-Y*(.044252-Y*1.9975E-3))) ! OS
      U2=9.7086+Y*(-.3814+Y*(.65292-Y*(.064984-Y*2.8792E-3)))
      U3=10.
      RETURN
77    U1=11.07+Y*(-2.412+Y*(1.9388+Y*(-.34389+Y*(.033511-1.3376E-3*Y)))) ! IR
      IF(T.GT.1.2E4) U1=30.
      U2=15.
      U3=20.
      RETURN
78    U1=16.4+1.27E-3*T ! PT
      du1=1.27e-3/u1
      U2=6.5712+Y*(-1.0363+Y*(.57234-Y*(.061219-2.6878E-3*Y)))
      du2=1.e-3*(-1.0363+y*(.57234*2.-y*(.061219*3.-y*4*2.6878e-3)))
      du2=du2/u2
      U3=15.
      RETURN
79    U1=1.24+2.79E-4*T ! AU
      du1=2.79e-4/u1
      U2=1.0546+Y*(-.040809+Y*(2.8439E-3+Y*1.6586E-3))
      du2=1.e-3*(-.040809+y*(2.8439e-3*2+y*1.6586e-3*3))/u2
      U3=7.
      RETURN
80    U1=1.0 ! HG
      U2=2.0
      U3=AMAX1(1.,.669+3.976E-5*T)
      if(u3.ne.1.)du3=3.976e-5/u3
      RETURN
81    U1=AMAX1(2.,0.63+3.35E-4*T) ! TL
      if(u1.ne.2.)du1=3.35e-4/u1
      U2=1.0
      U3=2.
      RETURN
82    U1=AMAX1(1.,0.42+2.35E-4*T) ! PB
      if(u1.ne.1.)du1=2.35e-4/u1
      if(T.GT.6125.)then
         U1=-1.2+5.E-4*T
         du1=5.e-4/u1
      end if
      U2=AMAX1(2.,1.72+7.9E-5*T)
      if(u2.ne.2.)du2=7.9e-5/u2
      U3=1.0
      RETURN
83    U1=2.78+2.87E-4*T ! BI
      du1=2.87e-4/u1
      U2=AMAX1(1.,.37+1.41E-4*T)
      if(u2.ne.1.)du2=1.41e-4/u2
      U3=2.5 ! APPROXIMATELY
      RETURN
84    U1=5. ! PO
      U2=5.
      U3=4.
      RETURN
85    U1=4.  ! AT
      U2=6.
      U3=6.
      RETURN
86    U1=1. ! RN
      U2=4.
      U3=6.
      RETURN
  87  U1=2. ! FR
      U2=1.
      U3=4.5
      RETURN
  88  U1=1. ! RA
      U2=2.
      U3=1.
      RETURN
  89  U1=6. ! AC
      U2=3.
      U3=7.
      RETURN
  90  U1=8. ! TH
      U2=8.
      U3=8.
      RETURN
  91  U1=50. ! PA
      U2=50.
      U3=50.
      RETURN
  92  U1=25. ! U
      U2=25.
      U3=25.
      return     
    End Subroutine Partition_f

! This routine extracts the Orbital Angular Momentum (L) from the
! letter codification
!
Subroutine Extract_OAM(Design, OAM)
Implicit NONE
Character :: Design
Character, Dimension (13) :: CODE
Character, Dimension (21) :: CODE1
Real :: OAM
Integer :: i
Data CODE/'S','P','D','F','G','H','I','K','L','M','N',   &
     'O','Q'/
Data code1/'p','1','f','2','h','3','k','4','m','5','o', &
     '6','r','7','t','8','u','9','v','0','w'/

OAM=-1.
Do i=1, 13
   If (Design .eq. CODE(i)) OAM=i-1.
End do

If (abs(OAM+1) .lt. 1.e-3) then
   Do i=1, 21
      If (Design .eq. CODE1(i)) OAM=float(i)/2.
   End do
End if

If (abs(OAM+1) .lt. 1.e-3) then
   Print *,'The letter used to code the Orbital Angular Momentum (',Design,')'
   Print *,'could not be recognized. STOP in extract_oam'
   Stop
End if

Return

End Subroutine Extract_OAM

  End Module Atomic_data
