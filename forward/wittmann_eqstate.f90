! Wittmann's Eq of state Pg=f(T,Pe) and Pe=f(T,Pg), adapted from SIR
  Module wittmann_eqstate
    Use Atomic_data
    Private
    Public :: Wittmann_compute_Pe, Wittmann_compute_Pg,&
         Wittmann_compute_others_from_T_Pe_Pg
    Contains
!      gasb calcula las preiones parciales y sus derivadas con          
!      respectoa la temperatura y la presion.es una modificacion        
!      de la rutina gas.                                                
                                                                        
!      dp es la derivada de p respecto a t, ddp respecto a pe.          

    subroutine wittmann_compute_pe(n,T,Pg,Pe)
      real, dimension(n) :: T, Pg, Pe, Pg2, Pe2
      integer :: n, loop

      Pe2=Pe
      Pg2=Pg
      Do loop=1, n
         Pe2(loop)=0.3*Pg2(loop)
         Call pefrompg10(T(loop),Pg2(loop),Pe2(loop))
      End do
    End subroutine wittmann_compute_pe
!
    subroutine wittmann_compute_pg(n,T,Pe,Pg)
      real, dimension(n) :: T, Pg, Pe, Pg2, Pe2
      real, dimension(99) :: p,dp,ddp
      real :: theta
      integer :: n, loop

      Pe2=Pe
      Do loop=1, n
         theta=5040./T(loop)
         call gasb(theta,Pe2(loop),p,dp,ddp)
         Pg(loop)=p(84)
      End do
    End subroutine wittmann_compute_pg
!
    subroutine wittmann_compute_others_from_T_pe_pg(n,T,Pe,Pg,&
         nH,nHminus,nHplus,nH2,nH2plus)
      real, dimension(n) :: T, Pg, Pe, Pg2, Pe2
      real, dimension(n) :: nH, nHminus, nHplus, nH2, nH2plus
      real, dimension(99) :: p,dp,ddp
      real :: theta, PHtot
      real, parameter :: bk=1.38066D-16
      integer :: n, loop

      Pe2=Pe
      Pg2=Pg
      Do loop=1, n
         theta=5040./T(loop)
         call gasb(theta,Pe2(loop),p,dp,ddp)
         Pg2(loop)=p(84)
         P2n=1./(bk*T(loop))
         PHtot=p(85)
         nH(loop)=p(1)*PHtot*P2n
         nHplus(loop)=p(86)*PHtot*P2n
         nHminus(loop)=p(87)*PHtot*P2n
         nH2(loop)=p(89)*PHtot*P2n
         nH2plus(loop)=p(88)*PHtot*P2n
      End do
    End subroutine wittmann_compute_others_from_T_pe_pg
                                                                        
      subroutine gasb(theta,pe,p,dp,ddp) 
!      ...............................................................  
      parameter (ncontr=28) 
      dimension cmol(91),alfai(ncontr),chi1(ncontr),chi2(ncontr),       &
     &u0(ncontr),u1(ncontr),u2(ncontr)                                  
      dimension du0(ncontr),du1(ncontr),du2(ncontr),dcmol(91) 
      real, dimension(99) :: p,dp,ddp
                                                                        
                                                                        
      t=5040./theta 
      if(pe.le.0)then 
         print*,'WARNING: Negative values of the electron pressure'
         print*,'         These are being changed to 1.e-10' 
         pe=1.e-10 
         g4=0. 
         g5=0. 
         dg4=0. 
         ddg4=0. 
         dg5=0. 
         ddg5=0. 
      else 
         call molecb(theta,cmol,dcmol) 
         do i=1,2 
            call acota(cmol(i),-30.,30.) 
         end do 
         g4=pe*10.**(cmol(1)) 
         dg4=dcmol(1)*alog(10.) 
         ddg4=1./pe 
         g5=pe*10.**(cmol(2)) 
         dg5=dcmol(2)*alog(10.) 
       ddg5=1./pe 
      end if 
                                                                        
                                                                        
!      ahora calculo los niveles u0,u1,u2 y sus derivadas               
      do 5 i=1,ncontr 
         iii=i 
         alfai(i)=10.**(At_abund(i)-12.)
         chi1(i)=At_ioniz1(i)
         chi2(i)=At_ioniz2(i)
    5 End do
    6 do 4 i=1,ncontr 
              iii=i 
    4         call partition_f(iii,t,u0(iii),u1(iii),u2(iii),du0(iii), &
                   du1(iii),du2(iii))
                                                                        
                                                                        
                                              ! p(h+)/p(h)              
      g2=saha(theta,chi1(1),u0(1),u1(1),pe) 
                                             !derivo log(g2) con t      
      dg2=dsaha(theta,chi1(1),du0(1),du1(1)) 
      ddg2=-1./pe 
                           !derivada log(g2) con pe                     
      p(92)=g2 
      dp(92)=dg2 
      ddp(92)=ddg2 
                                         ! p(h)/p(h-)                   
      g3=saha(theta,0.754,1.,u0(1),pe) 
                                                                        
      call acota(g3,1.e-30,1.e30) 
                                              ! p(h-)/p(h)              
      g3=1.d0/g3 
      dg3=-1.*dsaha(theta,0.754,0.,du0(1)) 
      ddg3=1./pe 
      p(93)=g3 
      dp(93)=dg3 
      ddp(93)=ddg3 
      g1=0. 
                    !las dl son derivadas no de log(g) sino de g        
      dlg1=0. 
      ddlg1=0. 
                                                                        
                                                                        
      do 1 i=2,ncontr 
      a=saha(theta,chi1(i),u0(i),u1(i),pe) 
      da=dsaha(theta,chi1(i),du0(i),du1(i)) 
      dda=-1./pe 
      b=saha(theta,chi2(i),u1(i),u2(i),pe) 
      dlb=b*dsaha(theta,chi2(i),du1(i),du2(i)) 
      ddlb=-b/pe 
      c=1.+a*(1.+b) 
      call acotasig(c,1.e-20,1.e20) 
                                                                        
                        ! p/ph' for neutral he,li, ...                  
      p(i)=alfai(i)/c 
      dp(i)=-(a*da*(1.+b)+a*dlb)/c 
      ddp(i)=-(a*dda*(1.+b)+a*ddlb)/c 
        ss1=(1.+2.*b) 
        call acotasig(ss1,1.e-20,1.e20) 
      ss=p(i)*a*ss1 
        dss=dp(i)+da+2.*dlb/ss1 
      ddss=ddp(i)+dda+2.*ddlb/ss1 
                                                                        
!      g1=g1+p(i)*a*ss1                                                 
      g1=g1+ss 
                                   !ojo estas no son derivadas del log  
      dlg1=dlg1+ss*dss 
    1  ddlg1=ddlg1+ss*ddss 
        a=1.+g2+g3 
      dla=g2*dg2+g3*dg3 
      ddla=g2*ddg2+g3*ddg3 
                                                                        
                                                                        
      if(g5.lt.1.e-35)then 
         print*,' ' 
           print*,'STOP: The electronic pressure is too small '
           print *,'subroutine wittmann_eqstate.f90'
           print*,'      Check the atmospheric models.' 
         print*,' ' 
         print*,'_____________________________________________________'
      stop 
      end if 
                                                                        
      call acotasig(g5,1.e-20,1.e20) 
      b=2.*(1.+g2/g5*g4) 
      dlb=(b-2.)*(dg2-dg5+dg4) 
      ddlb=(b-2.)*(ddg2-ddg5+ddg4) 
      c=g5 
      dlc=dg5*g5 
      ddlc=ddg5*g5 
      d=g2-g3 
      dld=g2*dg2-g3*dg3 
      ddld=g2*ddg2-g3*ddg3 
      e=g2/g5*g4 
      de=dg2-dg5+dg4 
      dde=ddg2-ddg5+ddg4 
      dle=e*de 
      ddle=e*dde 
                                                                        
      call acotasig(a,1.e-15,1.e15) 
      call acotasig(b,1.e-15,1.e15) 
      call acotasig(c,1.e-15,1.e15) 
      call acotasig(d,1.e-15,1.e15) 
      call acotasig(e,1.e-15,1.e15) 
                                                                        
                                                                        
      c1=c*b**2+a*d*b-e*a**2 
      dlc1=dlc*b*b+(c*2.*b+a*d)*dlb+dla*(d*b-2.*e*a)+dld*a*b-dle*a*a 
      ddlc1=ddlc*b*b+(c*2.*b+a*d)*ddlb+ddla*(d*b-2.*e*a)+ddld*a*b- &
           ddle*a*a
      c2=2.*a*e-d*b+a*b*g1 
      dlc2=dla*(2.*e+b*g1)+dlb*(a*g1-d)-dld*b+dle*2.*a+a*b*dlg1 
      ddlc2=ddla*(2.*e+b*g1)+ddlb*(a*g1-d)-ddld*b+ddle*2.*a+a*b*ddlg1 
      c3=-(e+b*g1) 
      dlc3=-dle-dlb*g1-b*dlg1 
      ddlc3=-ddle-ddlb*g1-b*ddlg1 
                                                                        
                                                                        
      call acotasig(c1,1.e-15,1.e15) 
                                                                        
      f1=0.5*c2/c1 
!      dlf1=.5*(dlc2/c1-(dlc1*c2)/(c1*c1))                              
                   !!!!!!!!                                             
      dc1=dlc1/c1 
      dc2=dlc2/c2 
      dlf1=f1*(dc2-dc1) 
                                                                        
!      ddlf1=.5*(ddlc2/c1-(ddlc1*c2)/(c1*c1))                           
                     !!!!!!!!                                           
      ddc1=ddlc1/c1 
      ddc2=ddlc2/c2 
      ddlf1=f1*(ddc2-ddc1) 
                                                                        
                                                                        
!      print*,'gasb 802 ',ddlf1,sign(1.,c1),2.*f1*dlf1,dlc3/c1,dlc1*c3/(
                                                                        
                                                                        
      dlf1=-dlf1+sign(1.,c1)*(2.*f1*dlf1-dlc3/c1+dlc1*c3/(c1*c1))       &
     &      /(2.*sqrt(f1**2-c3/c1))                                     
!      print*,'gasb 81 ',f1,dlf1,ddlf1                                  
                                                                        
      ddlf1=-ddlf1+sign(1.,c1)*(2.*f1*ddlf1-ddlc3/c1+ddlc1*c3/(c1*c1))  &
     &      /(2.*sqrt(f1**2-c3/c1))                                     
                                                                        
      f1=-f1+sign(1.,c1)*sqrt(f1**2-c3/c1) 
                                                                        
                                                                        
                                                                        
      f5=(1.-a*f1)/b 
      if(abs(f5).lt.1.e-30)then 
         dlf5=0. 
         ddlf5=0. 
         df5=0. 
         ddf5=0. 
      else 
         dlf5=(-dla*f1-a*dlf1)/b-((1.-a*f1)*dlb)/(b*b) 
         ddlf5=(-ddla*f1-a*ddlf1)/b-((1.-a*f1)*ddlb)/(b*b) 
         df5=dlf5/f5 
         ddf5=ddlf5/f5 
      end if 
        f4=e*f5 
      if(abs(f4).lt.1.e-30)then 
         dlf4=0. 
         ddlf4=0. 
         df4=0. 
         ddf4=0. 
      else 
         dlf4=f5*dle+e*dlf5 
         ddlf4=f5*ddle+e*ddlf5 
         df4=dlf4/f4 
         ddf4=ddlf4/f4 
      end if 
                                                                        
                                                                        
      f3=g3*f1 
      dlf3=f3*dg3+g3*dlf1 
      ddlf3=f3*ddg3+g3*ddlf1 
      f2=g2*f1 
      dlf2=f2*dg2+g2*dlf1 
      ddlf2=f2*ddg2+g2*ddlf1 
                                                                        
      if(abs(f1).lt.1.e-30)then 
         divf1=0. 
         ddivf1=0. 
      else 
         divf1=dlf1/f1 
         ddivf1=ddlf1/f1 
      end if 
                                                                        
      dlf200=dg2+divf1 
      ddlf200=ddg2+ddivf1 
                                                                        
                                                                        
        fe=f2-f3+f4+g1 
      if(abs(fe).lt.1.e-30)then 
         dlfe=0. 
         ddlfe=0. 
         dfe=0. 
         ddfe=0. 
      else 
         dlfe=dlf2-dlf3+dlf4+dlg1 
         ddlfe=ddlf2-ddlf3+ddlf4+ddlg1 
         dfe=dlfe/fe 
         ddfe=ddlfe/fe 
      end if 
                                                                        
                                                                        
      call acotasig(fe,1.e-15,1.e15) 
      phtot=pe/fe 
      dlphtot=-phtot*dfe 
      ddlphtot=(1.-ddlfe*phtot)/fe 
      dphtot=-dfe 
      ddphtot=1./pe-ddfe 
                                                                        
        if(f5.gt.1.e-4) goto 2 
                const6=g5/pe*f1**2 
                                            !dlf1/f1                    
            dconst6=dg5+2.*divf1 
                                               !ddlf1/f1                
            ddconst6=ddg5-1./pe+2.*ddivf1 
                  const7=f2-f3+g1 
                call acotasig(const7,1.e-15,1.e15) 
            dlconst7=dlf2-dlf3+dlg1 
            ddlconst7=ddlf2-ddlf3+ddlg1 
            dconst7=dlconst7/const7 
            ddconst7=ddlconst7/const7 
                                                                        
                  do 3 i=1,5 
                        f5=phtot*const6 
                  df5=dphtot+dconst6 
                  ddf5=ddphtot+ddconst6 
                        f4=e*f5 
                  df4=de+df5 
                  ddf4=dde+ddf5 
                        fe=const7+f4 
                        call acotasig(fe,1.e-15,1.e15) 
                  dfe=(dlconst7+df4*f4)/fe 
                  ddfe=(ddlconst7+ddf4*f4)/fe 
                  phtot=pe/fe 
                    dphtot=-dfe 
    3               ddphtot=1./pe-ddfe 
                                                                        
            dlf5=df5*f5 
            ddlf5=ddf5*f5 
            dlf4=df4*f4 
            ddlf4=ddf4*f4 
            dlfe=dfe*fe 
            ddlfe=ddfe*fe 
            dlphtot=dphtot*phtot 
            ddlphtot=ddphtot*phtot 
    2 pg=pe*(1.+(f1+f2+f3+f4+f5+0.1014)/fe) 
                                                                        
        dlpg=pe*(dlf1+dlf2+dlf3+dlf4+dlf5)/fe-(pg-pe)*(dlfe/fe) 
      ddlpg=pe*(ddlf1+ddlf2+ddlf3+ddlf4+ddlf5)/fe-(pg-pe)*(ddlfe/fe) 
      ddlpg=ddlpg+pg/pe 
                ! p(h)/p(h')                                            
      p(1)=f1 
                             !dlf1/f1                                   
      dp(1)=divf1 
                               !ddlf1/f1                                
      ddp(1)=ddivf1 
                                                                        
      call acotasig(pg,1.e-20,1.e20) 
                 ! gas pressure                                         
      p(84)=pg 
        dp(84)=dlpg/pg 
      ddp(84)=ddlpg/pg 
                    ! p(h')                                             
      p(85)=phtot 
                                     !dlphtot/phtot                     
      dp(85)=dphtot 
                                       !ddlphtot/phtot                  
      ddp(85)=ddphtot 
                 ! p(h+)/p(h')                                          
      p(86)=f2 
        dp(86)=dlf200 
      ddp(86)=ddlf200 
                 ! p(h-)/p(h')                                          
      p(87)=f3 
      if(abs(f3).lt.1.e-30)then 
         divf3=0. 
         ddivf3=0. 
      else 
         divf3=dlf3/f3 
         ddivf3=ddlf3/f3 
      end if 
                                                                        
      dp(87)=divf3 
      ddp(87)=ddivf3 
                 ! p(h2+)/p(h')                                         
      p(88)=f4 
                            !dlf4/f4                                    
      dp(88)=df4 
                              !ddlf4/f4                                 
      ddp(88)=ddf4 
                 ! p(h2)/p(h')                                          
      p(89)=f5 
                            !dlf5/f5                                    
      dp(89)=df5 
                              !ddlf5/f5                                 
      ddp(89)=ddf5 
                 ! pe/p(h')                                             
      p(90)=fe 
                            !dlfe/fe                                    
      dp(90)=dfe 
                              !ddlfe/fe                                 
      ddp(90)=ddfe 
                                                                        
                               ! n(e)=pe/kt                             
      p(91)=pe/(1.38054e-16*t) 
      dp(91)=-1./t 
      ddp(91)=1./pe 
!      print*,'gasb fin'                                                
                                                                        
                                                                        
      return 
    END Subroutine
!      ...............................................................  
!      gasc calcula las presiones parciales                             
                                                                        
                                                                        
      subroutine gasc(t,pe,pg,pp) 
      parameter (ncontr=28) 
      dimension cmol(91),alfai(ncontr),chi1(ncontr),chi2(ncontr),       &
     &u0(ncontr),u1(ncontr),u2(ncontr)                                  
      real du0,du1,du2,dcmol(91) 
      real pp(*),p(28) 
                                                                        
      theta=5040./t 
      call molecb(theta,cmol,dcmol) 
      g4=pe*10.**(cmol(1)) 
      g5=pe*10.**(cmol(2)) 
!      ahora calculo los niveles u0,u1,u2 y sus derivadas               
      do 5 i=1,ncontr 
                  iii=i 
                  alfai(i)=10.**(At_abund(i)-12.)
                  chi1(i)=At_ioniz1(i)
    5             chi2(i)=At_ioniz2(i)
    6 do 4 i=1,ncontr 
              iii=i 
    4         call partition_f(iii,t,u0(iii),u1(iii),u2(iii),du0,du1,du2) 
                                                                        
                                              ! p(h+)/p(h)              
      g2=saha(theta,chi1(1),u0(1),u1(1),pe) 
                                              ! p(h-)/p(h)              
      g3=1./saha(theta,0.754,1.,u0(1),pe) 
      pp(10)=g3 
      g1=0. 
      do 1 i=2,ncontr 
        a=saha(theta,chi1(i),u0(i),u1(i),pe) 
        b=saha(theta,chi2(i),u1(i),u2(i),pe) 
      c=1.+a*(1.+b) 
                          ! p/ph' for neutral he,li, ...                
        p(i)=alfai(i)/c 
    1  g1=g1+p(i)*a*(1.+2.*b) 
                                                                        
        pp(2)=p(2) 
        pp(3)=p(6) 
        pp(4)=p(11) 
        pp(5)=p(12) 
                                                                        
        a=1.+g2+g3 
        b=2.*(1.+g2/g5*g4) 
        c=g5 
        d=g2-g3 
        e=g2/g5*g4 
        c1=c*b**2+a*d*b-e*a**2 
        c2=2.*a*e-d*b+a*b*g1 
        c3=-(e+b*g1) 
        f1=0.5*c2/c1 
        f1=-f1+sign(1.,c1)*sqrt(f1**2-c3/c1) 
        f5=(1.-a*f1)/b 
        f4=e*f5 
        f3=g3*f1 
        f2=g2*f1 
        fe=f2-f3+f4+g1 
      phtot=pe/fe 
                                                                        
        if(f5.gt.1.e-4) goto 2 
          const6=g5/pe*f1**2 
          const7=f2-f3+g1 
                                                                        
        do 3 i=1,5 
                   f5=phtot*const6 
             f4=e*f5 
             fe=const7+f4 
             phtot=pe/fe 
    3    continue 
                                                                        
    2   pg=pe*(1.+(f1+f2+f3+f4+f5+0.1014)/fe) 
                   ! p(h)/p(h')                                         
        pp(1)=f1 
                   ! p(h+)/p(h')                                        
        pp(6)=f2 
                   ! p(h2)/p(h')                                        
        pp(7)=f5 
                   ! pe/p(h')                                           
        pp(8)=fe 
                                 ! n(e)=pe/kt                           
        pp(9)=pe/(1.38054e-16*t) 
                                                                        
      return 
   End Subroutine
!____________________________________________________________________   
! PEFROMPG10 (incuida en EQUISUBMU) evalua la presion electonica  desde 
! temperatura y la presion gaseosa                                      
! pefrompg10 evalua la presion electonica p correspondiente a t1 y pg   
                                                                        
      subroutine pefrompg10(t,pg,p) 
                                                                        
      implicit real (a-h,o-z) 
        prec=1e-2
                                                                        
        dif=1. 
        n2=0 
        p1=p 
        do while (dif.gt.prec.and.n2.lt.50.) 
           n2=n2+1 
           p=(p+p1)/2. 
             p1=p 
             call pe_pg10(t,p,pg) 
             dif=abs((p-p1)/p)  
         end do 
      return 
      END Subroutine                                        
!____________________________________________________________________   
! PE_PG10 (incuida en EQUISUBMU) evalua la presion electonica  desde la 
! la presion gaseosa y una estimacion de la presion electronica         
! calcula la presion electronica a partir de la pg y de una estimacion d
                                                                        
      subroutine pe_pg10(t,pe,pg) 
                                                                        
      parameter (ncontr=28) 
      real cmol(91),alfai(ncontr),chi1(ncontr),chi2(ncontr),          &
     &u0(ncontr),u1(ncontr),u2(ncontr)                                  
                                                                        
      real du0,du1,du2,dcmol(91) 
                                                                        
                                                                        
      if(t.lt.500)then 
         print*,'pe_pg10: temperature < 500 K ' 
       print*,'temperature = 500 K' 
       t=500. 
      end if 
      theta=5040./t 
      g4=0. 
      g5=0. 
      if(pe.le.0)then 
         pe=1.e-15 
         g4=0. 
         g5=0. 
      else 
         call molecb(theta,cmol,dcmol) 
         do i=1,2 
            call acota(cmol(i),-30.,30.) 
         end do 
         g4=pe*10.**(cmol(1)) 
         g5=pe*10.**(cmol(2)) 
      end if 
! ahora calculo los niveles u0,u1,u2 y sus derivadas                    
      do 5 i=1,ncontr 
         iii=i 
         alfai(i)=10.**(At_abund(i)-12.)
         chi1(i)=At_ioniz1(i)
         chi2(i)=At_ioniz2(i)
    5 end do

    6 do 4 i=1,ncontr 
             iii=i 
    4        call partition_f(iii,t,u0(iii),u1(iii),u2(iii),du0,du1,du2) 
                                                                        
                                                                        
                                              ! p(h+)/p(h)              
      g2=saha(theta,chi1(1),u0(1),u1(1),pe) 
                                                                        
                                              ! p(h)/p(h-)              
      g3=saha(theta,0.754,1.,u0(1),pe) 
                                                                        
      call acota(g3,1.e-30,1.e30) 
                                              ! p(h-)/p(h)              
      g3=1.d0/g3 
                                                                        
      g1=0. 
      do 1 i=2,ncontr 
        a=saha(theta,chi1(i),u0(i),u1(i),pe) 
        b=saha(theta,chi2(i),u1(i),u2(i),pe) 
      c=1.+a*(1.+b) 
    1  g1=g1+alfai(i)/c*a*(1.+2.*b) 
                                                                        
        a=1.+g2+g3 
        b=2.*(1.+g2/g5*g4) 
        c=g5 
        d=g2-g3 
        e=g2/g5*g4 
                                                                        
      call acotasig(a,1.e-15,1.e15) 
      call acotasig(d,1.e-15,1.e15) 
                                                                        
        c1=c*b**2+a*d*b-e*a**2 
                                                                        
        c2=2.*a*e-d*b+a*b*g1 
                                                                        
        c3=-(e+b*g1) 
                                                                        
        f1=0.5*c2/c1 
                                                                        
        f1=-f1+sign(1.,c1)*sqrt(f1**2-c3/c1) 
        f5=(1.-a*f1)/b 
        f4=e*f5 
        f3=g3*f1 
        f2=g2*f1 
        fe=f2-f3+f4+g1 
                                                                        
        call acota(fe,1.e-30,1.e30) 
      phtot=pe/fe 
                                                                        
        if(f5.gt.1.e-4) goto 2 
          const6=g5/pe*f1**2 
          const7=f2-f3+g1 
        do 3 i=1,5 
                   f5=phtot*const6 
             f4=e*f5 
             fe=const7+f4 
             phtot=pe/fe 
    3    continue 
                                                                        
    2   pe=pg/(1.+(f1+f2+f3+f4+f5+0.1014)/fe) 
        if(pe.le.0)pe=1.e-15 
                                                                        
      return 
      END Subroutine
! esta rutina interpola en ttau el modelo tau,t,pe y da la salida en    
! ttau,tt,ppe por medio de un polinomio de grado ngrado                 
                                                                        
      subroutine interpolatp(ngrado,n,tau,t,pe,nnew,ttau,tt,ppe) 
                                                                        
      real tau(*),t(*),pe(*),ttau(*),tt(*),ppe(*) 
      real xa(11),ya(11) 
                                                                        
! interpolaremos las presiones en logaritmos neperianos                 
      num=n 
      do i=1,num 
         pe(i)=alog(pe(i)) 
      end do 
                                                                        
! interpolamos                                                          
      n2=int(ngrado/2) 
                                                                        
      do i=1,nnew 
         CALL LOCATE(TAU,NUM,TTAU(I),J) 
         n3=j-n2-1 
           if(n3.lt.0)n3=0 
           if(n3+ngrado+1.gt.num)n3=num-ngrado-1 
         do k=1,ngrado+1 
            xa(k)=tau(n3+k) 
         end do 
         do k=1,ngrado+1 
            ya(k)=t(n3+k) 
         end do 
         CALL POLINT(XA,YA,NGRADO+1,TTAU(I),TT(I),ERROR) 
                                                                        
         do k=1,ngrado+1 
            ya(k)=pe(n3+k) 
         end do 
         CALL POLINT(XA,YA,NGRADO+1,TTAU(I),ppe(I),ERROR) 
        end do 
                                                                        
                                                                        
      do i=1,nnew 
           ppe(i)=exp(ppe(i)) 
        end do 
        do i=1,n 
           pe(i)=exp(pe(i)) 
      end do 
                                                                        
      return 
      END Subroutine

      subroutine acota(x,x0,x1) 
                                                                        
      if(x.lt.x0)x=x0 
        if(x.gt.x1)x=x1 
                                                                        
        return 
      END Subroutine

      subroutine acotasig(x,x0,x1) 
                                                                        
        if(x.lt.0)then 
           x=-x 
           call acota(x,x0,x1) 
           x=-x 
         else 
           call acota(x,x0,x1) 
         end if 
        return 
      END Subroutine

      SUBROUTINE LOCATE(XX,N,X,J)
      DIMENSION XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL
      RETURN
      END Subroutine

                                                                        
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY) 
      PARAMETER (NMAX=10) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX) 
      NS=1 
      DIF=ABS(X-XA(1)) 
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I)) 
        IF (DIFT.LT.DIF) THEN 
          NS=I 
          DIF=DIFT 
        ENDIF 
        C(I)=YA(I) 
        D(I)=YA(I) 
   11 END DO 
      Y=YA(NS) 
      NS=NS-1 
      DO 13 M=1,N-1 
        DO 12 I=1,N-M 
          HO=XA(I)-X 
          HP=XA(I+M)-X 
          W=C(I+1)-D(I) 
          DEN=HO-HP 
          IF(DEN.EQ.0.) Stop
          DEN=W/DEN 
          D(I)=HP*DEN 
          C(I)=HO*DEN 
   12   CONTINUE 
        IF (2*NS.LT.N-M)THEN 
          DY=C(NS+1) 
        ELSE 
          DY=D(NS) 
          NS=NS-1 
        ENDIF 
        Y=Y+DY 
   13 END DO 
      RETURN 
      End Subroutine

      SUBROUTINE MOLECb(X,Y,dy) 
      DIMENSION Y(91),dy(91) 
                                                                    ! H2
      Y(1)=-11.206998+X*(2.7942767+X*(7.9196803E-2-X*2.4790744E-2)) 
                                                                     ! H
      Y(2)=-12.533505+X*(4.9251644+X*(-5.6191273E-2+X*3.2687661E-3)) 
      dx=(-x*x)/5040. 
                                                                        
      dy(1)=dx*(2.7942767+x*(2*7.9196803e-2-x*3*2.4790744e-2)) 
      dy(2)=dx*(4.9251644+x*(-2*5.6191273e-2-x*3*3.2687661e-3)) 
      RETURN 
!     Y(3)=-11.89+X*(3.8084-X*2.4402E-2) ! CH                           
!     Y(4)=-11.85+X*(4.1411+X*(-6.7847E-2+4.9178E-3*X)) ! NH            
!     Y(5)=-12.199+X*(4.884+X*(-7.2794E-2+5.1747E-3*X)) ! OH            
!     Y(6)=-11.205+2.0112*X ! SIH                                       
!     Y(7)=-10.514+X*(2.2206-1.3654E-2*X) ! MGH                         
!     Y(8)=-10.581+X*(2.1713+X*(-7.7446E-2+6.5014E-3*X)) ! CAH          
!     Y(9)=-11.711+X*(3.1571-1.5205E-2*X) ! ALH                         
!     Y(10)=-12.24+X*(4.8009-2.6693E-2*X) ! HCL                         
!     Y(11)=-11.849+X*(4.1621-4.1213E-2*X) ! HS                         
!     Y(12)=-12.6+X*(6.3336-1.2019E-2*X) ! C2                           
!     Y(13)=-12.667+X*(7.779-1.1674E-2*X) ! CN                          
!     Y(14)=-13.641+X*(11.591+X*(-8.8848E-2+X*7.3465E-3)) ! CO          
!     Y(15)=-12.07+X*(4.7321-1.4647E-2*X) ! SIC                         
!     Y(16)=-13.122+X*(8.1528-2.5302E-2*X) ! CS                         
!     Y(17)=-13.435+X*(10.541+X*(-2.8061E-1+X*(5.883E-2+X*4.6609E-3)))  
!     Y(18)=-12.606+X*(6.9634+X*(-8.2021E-2+X*6.8654E-3)) ! NO          
!     Y(19)=-12.43+X*(5.409-2.341E-2*X) ! SIN                           
!     Y(20)=-12.172+X*(5.2755-2.3355E-2*X) ! SN                         
!     Y(21)=-13.087+X*(5.3673-1.4064E-2*X) ! O2                         
!     Y(22)=-13.034+8.3616*X ! SIO                                      
!     Y(23)=-11.39+X*(4.267-1.7738E-2*X) ! MGO                          
!     Y(24)=-12.493+X*(5.3382+X*(-5.8504E-2+4.8035E-3*X)) ! ALO         
!     Y(25)=-13.367+X*(8.869+X*(-7.1389E-1+X*(1.5259E-1-1.1909E-2*X)))  
!     Y(26)=-12.876+X*(6.8208+X*(-7.5817E-2+6.3995E-3*X)) ! VO          
!     Y(27)=-13.368+X*(9.0764+X*(-2.8354E-1+2.5784E-2*X)) ! ZRO         
!     Y(28)=-12.645+X*(5.6644-2.2882E-2*X) ! SO                         
!     Y(29)=-11.336+X*(4.4639-1.6552E-2*X) ! NACL                       
!     Y(30)=-12.515+6.7906*X ! SIS                                      
!     Y(31)=-10.638+4.2139*X ! CACL                                     
!     Y(32)=-12.06+X*(5.3309-1.6459E-2*X) ! ALCL                        
!     Y(33)=-12.508+X*(2.727-1.7697E-2*X) ! CL2                         
!     Y(34)=-12.651+X*(4.697-2.5267E-2*X) ! S2                          
!     Y(35)=-24.883+X*(8.2225-2.6757E-1*X) ! CH2                        
!     Y(36)=-24.82+X*(8.4594+X*(-1.2208E-1+8.6774E-3*X)) ! NH2          
!     Y(37)=-25.206+X*(10.311+X*(-9.0573E-2+5.3739E-3*X)) ! H2O         
!     Y(38)=-24.314+X*(8.1063-3.4079E-2*X) ! H2S                        
!     Y(39)=-25.168+13.401*X ! HCN                                      
!     Y(40)=-25.103+X*(12.87-3.8336E-2*X) ! HCO                         
!     Y(41)=-25.078+X*(9.3435+X*(-1.06E-1+8.2469E-3*X)) ! HNO           
!     Y(42)=-25.161+X*(7.8472+X*(-1.0399E-1+7.8209E-3*X)) ! HO2         
!     Y(43)=-27.038+X*(14.376+X*(6.8899E-2+6.0371E-2*X)) ! C3           
!     Y(44)=-25.889+13.317*X ! SIC2                                     
!     Y(45)=-27.261+X*(16.866-1.0144E-2*X) ! CO2                        
!     Y(46)=-26.25+X*(11.83+X*(-4.2021E-2+3.4397E-3*X)) ! N2O           
!     Y(47)=-26.098+X*(10.26+X*(-1.0101E-1+8.4813E-3*X)) ! NO2          
!     Y(48)=-26.115+X*(6.5385-1.9332E-2*X) ! O3                         
!     Y(49)=-27.496+X*(13.549+X*(-2.205E-2+2.1407E-2*X)) ! TIO2         
!     Y(50)=-27.494+X*(15.99+X*(-2.3605E-1+2.1644E-2*X)) ! ZRO2         
!     Y(51)=-25.244+X*(11.065-3.7924E-2*X) ! AL2O                       
!     Y(52)=-23.748+X*(9.5556-2.4867E-2*X) ! ALCL2                      
!     Y(53)=-37.194+X*(13.371-3.4229E-2*X) ! CH3                        
!     Y(54)=-37.544+X*(12.895-4.9012E-2*X) ! NH3                        
!     Y(55)=-37.931+17.216*X ! C2H2                                     
!     Y(56)=-38.274+X*(16.264-3.2379E-2*X) ! HCOH                       
!     Y(57)=-38.841+X*(18.907-3.5705E-2*X) ! HCNO                       
!     Y(58)=-50.807+X*(17.956-3.6861E-2*X) ! CH4                        
!     Y(59)=-11.4575+X*(3.1080922+X*(-3.3159806E-1+4.314945E-2*X)) ! NAH
!     Y(60)=-10.964723+X*(2.270225+X*(-7.66888E-2+6.519213E-3*X)) ! KH  
!     Y(61)=-10.807839+X*(2.744854+X*(5.758024E-2+3.315373E-3*X)) ! BEH 
!     Y(62)=-10.491008+X*(2.051217+X*(-7.643729E-2+6.425358E-3*X)) ! SRH
!     Y(63)=-11.01929+X*(3.13829+X*(1.214975-1.77077E-1*X)) ! SRO       
!     Y(64)=-10.446909+X*(2.024548+X*(-7.680739E-2+6.471443E-3*X)) ! BAH
!     Y(65)=-10.921254+X*(3.847116+X*(1.189653-1.662815E-1*X)) ! BAO    
!     Y(66)=-13.561415+X*(7.528355+X*(-5.031809E-1+6.787392E-2*X)) ! SCO
!     Y(67)=-14.107593+X*(12.226315+X*(-1.019148+1.02046E-1*X)) ! YO    
!     Y(68)=-14.231303+X*(11.28907+X*(-1.108545+1.274056E-1*X)) ! LAO   
!     Y(69)=-12.03193+X*(3.012432+X*(1.798894E-1-1.79236E-2*X)) ! SI2   
!     Y(70)=-11.344906+X*(2.836732+X*(-1.134115E-1+1.99661E-2*X)) ! LIH 
!     Y(71)=-25.913244+X*(11.856324+X*(1.05407-1.541093E-1*X)) ! VO2    
!     Y(72)=-26.934577+X*(13.421189+X*(2.671335E-1-3.475224E-2*X)) ! SIO
!     Y(73)=-23.225499+X*(4.820973+X*(6.722119E-1-5.903448E-2*X)) ! SIH2
!     Y(74)=-25.079417+X*(7.196588+X*(1.196713E-1+1.0484E-2*X)) ! SI3   
!     Y(75)=-26.331315+X*(12.500392+X*(6.531014E-1-1.162598E-1*X)) ! C2H
!     Y(76)=-11.673776+X*(3.245147+X*(1.334288E-1-1.524113E-2*X)) ! BH  
!     Y(77)=-12.972588+X*(7.80983-X*(6.263376E-2-4.763338E-3*X)) ! BO   
!     Y(78)=-12.654+X*(6.2558-3.0868E-2*X) ! HF                         
!     Y(79)=-11.639+X*(3.9742-X*(1.7229E-1-1.687E-2*X)) ! LIO           
!     Y(80)=-11.92+X*(6.1426-1.8981E-2*X) ! LIF                         
!     Y(81)=-11.576+X*(5.1447+X*(1.4625E-2-8.9675E-3*X)) ! LICL         
!     Y(82)=-12.579+X*(4.8824+X*(-1.3848E-1+1.25E-2*X)) ! FEO           
!     Y(83)=-11.666+X*(5.1748-1.9602E-2*X) ! NAF                        
!     Y(84)=-11.292+X*(4.851-1.8104E-2*X) ! MGF                         
!     Y(85)=-12.453+X*(7.1023-1.6086E-2*X) ! ALF                        
!     Y(86)=-11.8+X*(5.7274-X*(1.9201E-1-2.1306E-2*X)) ! KF             
!     Y(87)=-10.798+X*(2.6664+X*(1.2037E-1-1.778E-2*X)) ! MGCL          
!     Y(88)=-11.453+X*(4.9299+X*(-1.4643E-1+1.4143E-2*X)) ! KCL         
!     Y(89)=-11.484+X*(3.2399-2.2356E-2*X) ! FECL                       
!     Y(90)=-24.304+X*(9.5257-4.2841E-2*X) ! LIOH                       
!     Y(91)=0.                                                          
    END SUBROUTINE MOLECb

      FUNCTION SAHA(THETA,CHI,U1,U2,PE) 
!     SAHA-EGGERT EQUATION                                              
      SAHA=U2*10.**(9.0805126-THETA*CHI)/(U1*PE*THETA**2.5) 
!      print*,'saha en saha =          ',saha                           
!      print*,theta,chi,u1,u2,pe                                        
      RETURN 
    END FUNCTION SAHA
                                                                        
                                                                        
      FUNCTION dSAHA(THETA,CHI,dU1,dU2) 
!	calcula la detivada respecto a t, del logaritmo de la                 
!     SAHA-EGGERT EQUATION                                              
!	la derivada del logaritmo de saha respecto a pe es -1/pe              
      dSAHA=du2-du1+(theta/5040.)*(2.5+chi*THETA*alog(10.)) 
      RETURN 
    END FUNCTION dSAHA

 End Module Wittmann_eqstate
