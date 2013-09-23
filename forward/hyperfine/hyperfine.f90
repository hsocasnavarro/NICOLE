! MODULE HYPERFINE
! Author: Andres Asensio-Ramos (IAC)
! 
! --
! Modifications:
!           2011-09-07 JdlCR: Created function Lchar2double.
!--
!
module hyperfine
  implicit none
  Real(kind=8), dimension(0:301) :: fact
contains
  function Lchar2double(din) result(res)
    Implicit None
    integer :: i
    character(len=*) :: din
    character(len=1), dimension(13) :: let
    real(kind=8) :: res
    character(len=1), dimension(21) :: let1 
    logical :: doit
    Data let1(:)/'p','1','f','2','h','3','k','4','m','5','o','6','r','7','t','8','u','9','v','0','w'/
    Data let(:)/'S','P','D','F','G','H','I','K','L','M','N','O','Q'/
    !
    doit = .True.
    !
    ! check for integer angular momentum
    !
    i = 1
    do while(doit.AND.(i.LE.13)) 
       if(din.EQ.let(i)) then 
          doit = .False.
          res = i - 1
          return
       else
          i = i + 1
       end if
    end do
    !
    ! Check for non-integer angular momentum
    !
    i = 1
    do while(doit.AND.(i.LE.21))
       if(din.EQ.let1(i)) then 
          doit = .false.
          res = i * 0.5d0
          return
       else
          i = i + 1
       endif
    end do
    !
    print *, 'Lchar2double : Error -> ' , din , ' is not a valid descriptor for angular momentum'
    !
  end function Lchar2double
  subroutine factrl
    integer :: i
    fact(0) = 1.d0
    do i=1,301
       fact(I) = fact(I-1) * dble(I)
    enddo
  END subroutine factrl
  
  !----------------------------------------------------------------
  ! This function calculates the 3-j symbol
  ! J_i and M_i have to be twice the actual value of J and M
  !----------------------------------------------------------------
  function w3js(j1,j2,j3,m1,m2,m3)
    integer :: m1, m2, m3, j1, j2, j3
    integer :: ia, ib, ic, id, ie, im, ig, ih, z, zmin, zmax, jsum
    real(kind=8) :: w3js, cc, denom, cc1, cc2
    
    !          w3js = w3js_regge(j1/2,j2/2,j3/2,m1/2,m2/2,m3/2)
    !          return
    w3js = 0.d0
    if (m1+m2+m3 /= 0) return
    ia = j1 + j2
    if (j3 > ia) return
    ib = j1 - j2
    if (j3 < abs(ib)) return
    if (abs(m1) > j1) return
    if (abs(m2) > j2) return
    if (abs(m3) > j3) return
    
    jsum = j3 + ia
    ic = j1 - m1
    id = j2 - m2
    
    ie = j3 - j2 + m1
    im = j3 - j1 - m2
    zmin = max0(0,-ie,-im)
    ig = ia - j3
    ih = j2 + m2
    zmax = min0(ig,ih,ic)
    cc = 0.d0
    
    do z = zmin, zmax, 2
       denom = fact(z/2)*fact((ig-z)/2)*fact((ic-z)/2)*fact((ih-z)/2)*&
            fact((ie+z)/2)*fact((im+z)/2)
       if (mod(z,4) /= 0) denom = -denom
       cc = cc + 1.d0/denom
    enddo
    
    cc1 = fact(ig/2)*fact((j3+ib)/2)*fact((j3-ib)/2)/fact((jsum+2)/2)
    cc2 = fact((j1+m1)/2)*fact(ic/2)*fact(ih/2)*fact(id/2)*fact((j3-m3)/2)*fact((j3+m3)/2)
    cc = cc * sqrt(1.d0*cc1*cc2)
    if (mod(ib-m3,4) /= 0) cc = -cc
    w3js = cc
    if (abs(w3js) < 1.d-8) w3js = 0.d0          
1000 return
  end function w3js
  
  
  !----------------------------------------------------------------
  ! This function calculates the 6-j symbol
  ! J_i and M_i have to be twice the actual value of J and M
  !----------------------------------------------------------------
  function w6js(j1,j2,j3,l1,l2,l3)
    integer :: j1,j2,j3,l1,l2,l3
    integer :: ia, ib, ic, id, ie, iif, ig, ih, sum1, sum2, sum3, sum4
    integer :: w, wmin, wmax, ii, ij, ik
    real(kind=8) :: w6js, omega, theta1, theta2, theta3, theta4, theta, denom
    
    w6js = 0.d0
    ia = j1 + j2
    if (ia < j3) return
    ib = j1 - j2
    if (abs(ib) > j3) return
    ic = j1 + l2
    if (ic < l3) return
    id = j1 - l2
    if (abs(id) > l3) return
    ie = l1 + j2
    if (ie < l3) return
    iif = l1 - j2
    if (abs(iif) > l3) return
    ig = l1 + l2
    if (ig < j3) return
    ih = l1 - l2
    if (abs(ih) > j3) return
    sum1=ia + j3
    sum2=ic + l3
    sum3=ie + l3
    sum4=ig + j3
    wmin = max0(sum1, sum2, sum3, sum4)
    ii = ia + ig
    ij = j2 + j3 + l2 + l3
    ik = j3 + j1 + l3 + l1
    wmax = min0(ii,ij,ik)
    omega = 0.d0
    do w = wmin, wmax, 2
       denom = fact((w-sum1)/2)*fact((w-sum2)/2)*fact((w-sum3)/2)&
            *fact((w-sum4)/2)*fact((ii-w)/2)*fact((ij-w)/2)&
            *fact((ik-w)/2)
       if (mod(w,4) /= 0) denom = -denom
       omega = omega + fact(w/2+1) / denom
    enddo
    theta1 = fact((ia-j3)/2)*fact((j3+ib)/2)*fact((j3-ib)/2)/fact(sum1/2+1)
    theta2 = fact((ic-l3)/2)*fact((l3+id)/2)*fact((l3-id)/2)/fact(sum2/2+1)
    theta3 = fact((ie-l3)/2)*fact((l3+iif)/2)*fact((l3-iiF)/2)/fact(sum3/2+1)
    theta4 = fact((ig-j3)/2)*fact((j3+ih)/2)*fact((j3-ih)/2)/fact(sum4/2+1)
    theta = theta1 * theta2 * theta3 * theta4
    w6js = omega * sqrt(theta)
    if (abs(w6js) < 1.d-8) w6js = 0.d0
1000 return
  end function w6js
  
  
  ! ---------------------------------------------------------
  ! Calculate the sizes of the eigenvalues and eigenvectors
  ! ---------------------------------------------------------
  subroutine hyper_size(J, I, nflargest, f2max)
    real(kind=8) :: J, I
    integer :: nflargest, F2min, F2max, m2, f2a, f2b, nf
    
    F2min = abs(2*J-2*I)
    F2max = 2*J+2*I
    
    nflargest = 0
    do m2 = -f2max, f2max, 2
       f2a = abs(m2)
       if (f2a <= f2min) f2a = f2min               
       f2b = f2max
       nf = 1 + (f2b-f2a) / 2.d0
       if (nf > nflargest) then 
          nflargest = nf
       endif
    enddo
    
  end subroutine hyper_size
  
  ! ---------------------------------------------------------
  ! Calculate the sizes of the eigenvalues and eigenvectors
  ! ---------------------------------------------------------
  subroutine hyper_components(Ju, Jl, I, f2m_up, f2m_low, nflev_up, nflev_low, nPi, nLeft, nRight)
    real(kind=8) :: Ju, Jl, I
    integer :: nPi, nLeft, nRight, m2u, m2l, f2a_up, f2a_low, iu, il, f2m_up, f2m_low
    integer :: f2min_up, f2max_up, f2min_low, f2max_low
    integer :: nflev_up(0:2*f2m_up), nflev_low(0:2*f2m_low)
    real(kind=8) :: deltaM, q
    
    f2min_up = abs(2*Ju-2*I)
    f2max_up = 2*Ju+2*I
    
    f2min_low = abs(2*Jl-2*I)
    f2max_low = 2*Jl+2*I
    
    
    ! First count the number of hyperfine components
    nPi = 0
    nRight = 0
    nLeft = 0
    do m2u = -f2max_up, f2max_up, 2
       
       f2a_up = abs(m2u)
       if (f2a_up <= f2min_up) f2a_up = f2min_up
       
       do m2l = -f2max_low, f2max_low, 2
          
          f2a_low = abs(m2l)
          if (f2a_low <= f2min_low) f2a_low = f2min_low
          
          deltaM = (m2u - m2l) / 2.d0
          q = -deltaM
          if (abs(deltaM) <= 1) then                         
             do iu = 0, nflev_up(f2max_up+m2u)-1
                do il = 0, nflev_low(f2max_low+m2l)-1
                   if (q == 0) nPi = nPi + 1
                   if (q == -1) nLeft = nLeft + 1
                   if (q == 1) nRight = nRight + 1
                enddo
             enddo
          endif
       enddo
    enddo
    
  end subroutine hyper_components
  
  ! ---------------------------------------------------------
  ! Calculate the hyperfine eigenvalues and eigenvectors
  ! ---------------------------------------------------------
  subroutine hyper_diagonalize(L, S, J, I, A, B, Bfield, nflev, energy, c, f2m, nfl)
    USE nr, ONLY : tqli
    integer :: f2m, nfl
    real(kind=8) :: L, S, J, I, A, B, Bfield, energy(0:2*f2m,0:nfl-1), c(0:2*f2m,0:nfl-1,0:f2m)
    
    integer :: nflev(0:2*f2m), F2min, F2max, m2, f2a, f2b, nu, nl, loop, nf, l1, l2
    real(kind=8) :: gLS, mu0, Fu, Fl, t1, t1_diag, t2_diag, K, M
    real(kind=8), allocatable :: H(:,:)
    real(kind=8), allocatable :: d(:), ud(:), z(:,:)
    real, allocatable :: d4(:), ud4(:), z4(:,:)
    
    gLS = 1.d0 + 0.5d0*(J*(J+1)+S*(S+1)-L*(L+1)) / (J*(J+1))
    
    F2min = abs(2*J-2*I)
    F2max = 2*J+2*I
    
    mu0 = 9.27d-21 / (6.626d-27*3.d10)
    
    do m2 = -f2max, f2max, 2
       M = m2 / 2.d0
       f2a = abs(m2)
       if (f2a <= f2min) f2a = f2min
       f2b = f2max
       nf = 1 + (f2b-f2a) / 2.d0
       nflev(f2max+m2) = nf
       
       allocate(H(nf,nf))
       
       do nu = 1, nf
          Fu = (f2a + 2*(nu-1)) / 2.d0
          do nl = 1, nf
             Fl = (f2a + 2*(nl-1)) / 2.d0
             
             ! Non-diagonal part
             t1 = sqrt(J*(J+1)*(2*J+1)*(2*Fu+1)*(2*Fl+1)) * &
                  w6js(int(2*Fl),int(2*Fu),2,int(2*J),int(2*J),int(2*I)) * &
                  w3js(int(2*Fu),int(2*Fl),2,-int(2*M),int(2*M),0)
             
             H(nu,nl) = mu0*Bfield*gLS*(-1.d0)**(J+I-M) * t1
             
             ! Diagonal part
             if (fu == fl) then
                K = Fu*(Fu+1) - I*(I+1) - J*(J+1)
                t1_diag = 0.5d0 * A * K
                t2_diag = 0.d0
                if (I /= 0 .and. I /= 0.5d0 .and. J /= 0 .and. J /= 0.5d0) then
                   t2_diag = 0.5d0 * B * (0.75d0*K*(K+1)-I*(I+1)*J*(J+1)) / (I*(2*I-1)*J*(2*J-1))
                endif
                H(nu,nl) = H(nu,nl) + t1_diag + t2_diag
             endif
             
          enddo
       enddo
       
       
       
       ! The matrix is tridiagonal. We use an appropriate routine for the diagonalization
       allocate(d(nf),d4(nf))
       d = 0.d0
       allocate(ud(nf),ud4(nf))
       ud = 0.d0
       allocate(z(nf,nf),z4(nf,nf))
       
       z = 0.d0
       do loop = 1, nf
          z(loop,loop) = 1.d0
       enddo
       
       do loop = 1, nf
          d(loop) = H(loop,loop)
       enddo
       
       do loop = 1, nf-1
          ud(loop+1) = H(loop+1,loop)
       enddo
       
       d4=d
       ud4=ud
       z4=z
       call tqli(d4,ud4,z4)
       d=d4
       ud=ud4
       z=z4
       
       ! Store the eigenvalues and eigenvectors            
       do l1 = 0, nf-1
          energy(f2max+m2,l1) = d(l1+1)
          do l2 = 0, nf-1
             Fu = f2a + 2*l2
             c(f2max+m2,l1,int(Fu)) = z(l2+1,l1+1)
          enddo
       enddo
       
       deallocate(d,d4)
       deallocate(ud,ud4)
       deallocate(z,z4)
       deallocate(H)
    enddo
    
  end subroutine hyper_diagonalize
  
  ! ---------------------------------------------------------
  ! Calculate the hyperfine pattern
  ! ---------------------------------------------------------
  subroutine hyper_pattern(I, Lu, Su, Ju, Au, Bu, Ll, Sl, Jl, Al, Bl, Bfield, &
       nflev_up, energy_up, c_up, f2m_up, nfl_up, nflev_low, energy_low, c_low, f2m_low, nfl_low, &
       SplitPi, SplitRight, SplitLeft, StrPi, StrRight, StrLeft)
    Use Param_structure
    Implicit None
    integer :: f2m_up, nfl_up, f2m_low, nfl_low
    real(kind=8) :: I, Lu, Su, Ju, Au, Bu, Ll, Sl, Jl, Al, Bl
    real(kind=8) :: Bfield
    real(kind=8) :: energy_up(0:2*f2m_up,0:nfl_up-1), energy_low(0:2*f2m_low,0:nfl_low-1)
    real(kind=8) :: c_up(0:2*f2m_up,0:nfl_up-1,0:f2m_up), c_low(0:2*f2m_low,0:nfl_low-1,0:f2m_low)
    integer :: nflev_up(0:2*f2m_up), nflev_low(0:2*f2m_low)
    
    integer :: f2min_up, f2max_up, f2min_low, f2max_low, ncomponents, m2u, m2l, f2a_up, f2a_low, iu, il
    integer :: nPi, nRight, nLeft, F2, Fp2, Fs2, Ft2
    integer :: nnPi, nnRight, nnLeft
    real :: deltaM, q, splitting, sum, t1, t2, t3, t4
    
    real, dimension(transitions) :: SplitPi, SplitRight, SplitLeft, StrPi, StrRight, StrLeft
    
    f2min_up = abs(2*Ju-2*I)
    f2max_up = 2*Ju+2*I
    
    f2min_low = abs(2*Jl-2*I)
    f2max_low = 2*Jl+2*I      
    
    nPi = 0
    nRight = 0
    nLeft = 0
    
    do m2u = -f2max_up, f2max_up, 2
       
       f2a_up = abs(m2u)
       if (f2a_up <= f2min_up) f2a_up = f2min_up
       
       do m2l = -f2max_low, f2max_low, 2
          
          f2a_low = abs(m2l)
          if (f2a_low <= f2min_low) f2a_low = f2min_low
          
          deltaM = (m2u - m2l) / 2.d0
          q = -deltaM
          if (abs(deltaM) <= 1) then                         
             do iu = 0, nflev_up(f2max_up+m2u)-1
                do il = 0, nflev_low(f2max_low+m2l)-1
                   
                   splitting = energy_up(f2max_up+m2u,iu) - energy_low(f2max_low+m2l,il)
                   
                   sum = 0.d0
                   
                   do F2 = f2a_low, f2max_low, 2
                      
                      do Fs2 = f2a_low, f2max_low, 2
                         
                         if (iabs(Fs2-F2) > 2 .or. Fs2+F2 < 2) then
                         else
                            do Fp2 = f2a_up, f2max_up, 2
                               
                               do Ft2 = f2a_up, f2max_up, 2
                                  
                                  if (iabs(Ft2-Fs2) > 2 .or. Ft2+Fs2 < 2) then
                                  else
                                     t1 = c_low(f2max_low+m2l,il,F2) * c_low(f2max_low+m2l,il,Fs2) * &
                                          c_up(f2max_up+m2u,iu,Fp2) * c_up(f2max_up+m2u,iu,Ft2)
                                     t2 = 3.d0 / (2.d0*I+1.d0) * sqrt((F2+1.d0)*(Fp2+1.d0)*(Fs2+1.d0)*(Ft2+1.d0))
                                     t3 = w6js(int(2*Jl),int(2*Ju),2,int(Fp2),int(F2),int(2*I)) * &
                                          w6js(int(2*Jl),int(2*Ju),2,int(Ft2),int(Fs2),int(2*I))
                                     t4 = w3js(int(Fp2),int(F2),2,-int(m2u),int(m2l),-int(2*q)) * &
                                          w3js(int(Ft2),int(Fs2),2,-int(m2u),int(m2l),-int(2*q))
                                     
                                     sum = sum + t1 * t2 * t3 * t4                                        
                                  endif
                               enddo
                            enddo
                         endif
                      enddo
                   enddo
                   
                   if (q == 0) then
                      nPi = nPi + 1
                      SplitPi(nPi) = splitting
                      StrPi(nPi) = sum
                   endif
                   
                   if (q == -1) then
                      nLeft = nLeft + 1
                      SplitLeft(nLeft) = splitting
                      StrLeft(nLeft) = sum                                        
                   endif
                   
                   if (q == 1) then
                      nRight = nRight + 1
                      SplitRight(nRight) = splitting
                      StrRight(nRight) = sum                                        
                   endif
                enddo
             enddo
          endif
       enddo
    enddo
  end subroutine hyper_pattern
  
end module hyperfine

