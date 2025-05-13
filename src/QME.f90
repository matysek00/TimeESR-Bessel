module QME
Use declarations
Use timer
CONTAINS
     subroutine rates (Ndim, NFreq, Ntime, Nbias, lambda, gamma_R_0, gamma_L_0,  &
         Spin_polarization_R, Spin_polarization_L, fermiR_a, fermiL_a,       &
         ufermiR_a, ufermiL_a, Temperature, G)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the QME rates
! time dependent
! for pulses that are either steps or cosine or a sequence of both
     implicit none
! Input:
     complex (qc), intent (in):: lambda (:,:,:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, Temperature
     real (q), intent (in):: Spin_polarization_R, Spin_polarization_L
     integer :: Ndim, Ntime, NFreq, Nbias
! Output: the Rates called G (:,:,:,:,:,:) here
     complex (qc) :: G (:,:,:,:,:,:) ! for QME
! Computed in ExtendedFermiIntegral
     complex (qc) :: fermiR, fermiL, ufermiR, ufermiL
     complex (qc) :: fermiR_a(:,:,:), fermiL_a(:,:,:)
     complex (qc) :: ufermiR_a(:,:,:), ufermiL_a(:,:,:)
! Only used in this subroutine
     integer :: v, l, j, u, n
     complex (qc) :: g0p_up, g0m_up, g1p_up, g1m_up
     complex (qc) :: g0p_dn, g0m_dn, g1p_dn, g1m_dn, pulse
     complex (qc) :: Lvluj, Ljulv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do n=1,Nbias
     do j=1,Ndim
     do u=1,Ndim

      fermiR = fermiR_a(j,u,n)
      fermiL = fermiL_a(j,u,n)
      ufermiR = ufermiR_a(j,u,n)
      ufermiL = ufermiL_a(j,u,n)
      
! Right electrode
      
      g0p_up = gamma_R_0*fermiR*(1+Spin_polarization_R)*0.5
      g0m_up = gamma_R_0*ufermiR*(1+Spin_polarization_R)*0.5
      g0p_dn = gamma_R_0*fermiR*(1-Spin_polarization_R)*0.5
      g0m_dn = gamma_R_0*ufermiR*(1-Spin_polarization_R)*0.5

! Left electrode
      g1p_up = gamma_L_0*fermiL*(1+Spin_polarization_L)*0.5
      g1m_up = gamma_L_0*ufermiL*(1+Spin_polarization_L)*0.5
      g1p_dn = gamma_L_0*fermiL*(1-Spin_polarization_L)*0.5
      g1m_dn = gamma_L_0*ufermiL*(1-Spin_polarization_L)*0.5

        do v=1,Ndim
        do l=1,Ndim
            Lvluj = lambda (v,l,1)*conjg(lambda(u,j,1))*g0p_up+  &
                    lambda (v,l,2)*conjg(lambda(u,j,2))*g0p_dn
            Ljulv = lambda (j,u,1)*conjg(lambda(l,v,1))*g0m_up+  &
                    lambda (j,u,2)*conjg(lambda(l,v,2))*g0m_dn

! Right electrode
            G (v,l,j,u,1,n) = 0.5*(Lvluj + Ljulv)

            Lvluj = lambda (v,l,1)*conjg(lambda(u,j,1))*g1p_up+  &
                    lambda (v,l,2)*conjg(lambda(u,j,2))*g1p_dn
            Ljulv = lambda (j,u,1)*conjg(lambda(l,v,1))*g1m_up+  &
                    lambda (j,u,2)*conjg(lambda(l,v,2))*g1m_dn

! Left electrode
            G (v,l,j,u,2,n) = 0.5*(Lvluj + Ljulv)
            
            !if(write_rates) write(666,*) v, l, j, u, n, G(v,l,j,u,1,n)*Hartree/GHz, G(v,l,j,u,2,n)*Hartree/GHz
        enddo
        enddo
      enddo
     enddo
     enddo

     return

     end subroutine rates 
!    TODO: GC(v,l,j,u,n) = Electrode*G(v,l,j,u,1,n) + (1-Electrode)*G(v,l,j,u,2,n)
!         I don't think that ratesC functions should exist. 

!
! For the Current NOW
!

subroutine ratesC (Ndim, NFreq, Nbias, lambda, gamma_R_0, gamma_L_0,  &
     Spin_polarization_R, Spin_polarization_L, fermiR_a, fermiL_a, ufermiR_a, ufermiL_a, &
     Temperature, Electrode,  GC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the QME rates
! time dependent
! for pulses that are either steps or cosine or a sequence of both
     implicit none
! Input:
     complex (qc), intent (in):: lambda (:,:,:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, Temperature
     real (q), intent (in):: Spin_polarization_R, Spin_polarization_L
     integer :: Ndim, NFreq, Nbias
     integer :: Electrode
! Output: the Rates called GC (:,:,:,:,:) here
     complex (qc) :: GC (:,:,:,:,:) ! for current
! Computed in ExtendedFermiIntegral
     complex (qc) :: fermiR, fermiL, ufermiR, ufermiL
     complex (qc) :: fermiR_a(:,:,:), fermiL_a(:,:,:)
     complex (qc) :: ufermiR_a(:,:,:), ufermiL_a(:,:,:)
! Only used in this subroutine
     integer :: v, l, j, u, n, i, m
     complex (qc) :: g0pa_up, g0ma_up, g1pa_up, g1ma_up
     complex (qc) :: g0pa_dn, g0ma_dn, g1pa_dn, g1ma_dn
     complex (qc) :: Lvluja, Ljulva
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          
     do n =1,Nbias
     do j=1,Ndim
     do u=1,Ndim

! Static contribution
! Speed up if we create a Fermi array, we leave it as a function
! but it slows down the code (only two evaluation per pair u,j )


     fermiR = fermiR_a(j,u,n)
     fermiL = fermiL_a(j,u,n)
     ufermiR = ufermiR_a(j,u,n)
     ufermiL = ufermiL_a(j,u,n)
     

!Right electrode
     g0pa_up = Electrode*gamma_R_0*fermiR*(1+Spin_polarization_R)*0.5
     g0ma_up = Electrode*gamma_R_0*ufermiR*(1+Spin_polarization_R)*0.5
     g0pa_dn = Electrode*gamma_R_0*fermiR*(1-Spin_polarization_R)*0.5
     g0ma_dn = Electrode*gamma_R_0*ufermiR*(1-Spin_polarization_R)*0.5

! Left electrode
     g1pa_up = (1-Electrode)*gamma_L_0*fermiL*(1+Spin_polarization_L)*0.5
     g1ma_up = (1-Electrode)*gamma_L_0*ufermiL*(1+Spin_polarization_L)*0.5
     g1pa_dn = (1-Electrode)*gamma_L_0*fermiL*(1-Spin_polarization_L)*0.5
     g1ma_dn = (1-Electrode)*gamma_L_0*ufermiL*(1-Spin_polarization_L)*0.5

     do v=1,Ndim
     do l=1,Ndim

!Right electrode
          Lvluja = lambda (v,l,1)*conjg(lambda(u,j,1))*g0pa_up+  &
                    lambda (v,l,2)*conjg(lambda(u,j,2))*g0pa_dn
          Ljulva = lambda (j,u,1)*conjg(lambda(l,v,1))*g0ma_up+  &
                    lambda (j,u,2)*conjg(lambda(l,v,2))*g0ma_dn

          GC (v,l,j,u,n) = 0.5*(Lvluja - Ljulva)

! Left electrode
          Lvluja = lambda (v,l,1)*conjg(lambda(u,j,1))*g1pa_up+  &
                    lambda (v,l,2)*conjg(lambda(u,j,2))*g1pa_dn
          Ljulva = lambda (j,u,1)*conjg(lambda(l,v,1))*g1ma_up+  &
                    lambda (j,u,2)*conjg(lambda(l,v,2))*g1ma_dn
          
          GC (v,l,j,u,n) = GC (v,l,j,u,n) + 0.5*(Lvluja - Ljulva) ! This has the right admixture of electrodes

     enddo
     enddo
     enddo
     enddo
     enddo
     
     return
 
     end subroutine ratesC 

     subroutine rates_bessel (Ndim, NFreq, Nbias, lambda, gamma_R_0, gamma_L_0,  &
         Spin_polarization_R, Spin_polarization_L, &
         p_max, B_R, B_L, Amplitude, frequency, bias_R, bias_L,&
         Temperature, Electrode,  G)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the QME rates
! time dependent
! for pulses that are either steps or cosine or a sequence of both
     implicit none
! Input:
     complex (qc), intent (in):: lambda (:,:,:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, Temperature, frequency
     real (q), intent (in):: Spin_polarization_R, Spin_polarization_L, B_L, B_R
     real (q), intent (in):: bias_R , bias_L 
     real (q), intent (in), dimension(2) :: Amplitude
     integer :: Ndim, NFreq, Nbias, p_max
     integer :: Electrode
! Output: the Rates called GC (:,:,:,:,:) here
     complex (qc) :: G (:,:,:,:,:,:)
! Only used in this subroutine
     integer :: v, l, j, u, n, n_index
     complex (qc), dimension(2*p_max-1) :: fermiRB, fermiLB, ufermiRB, ufermiLB
     complex (qc), dimension(Ndim,Ndim) :: LvlujaL, LjulvaL, LvlujaR, LjulvaR
     complex (qc), dimension(2*p_max-1) :: K_L, K_R
     complex (qc) :: g0pa_up, g1pa_up, g0pa_dn, g1pa_dn
     complex (qc) :: bessel_contributionR, ubessel_contributionR
     complex (qc) :: bessel_contributionL, ubessel_contributionL

!    Calculate Contribution of Bessel functions
!    K(p) = J(p) + .5*A*(J(p-1)+ J(p+1))
     K_R = Bessel_K(B_R/frequency, Amplitude(1), p_max)
     K_L = Bessel_K(B_L/frequency, Amplitude(2), p_max)

     g0pa_up = 0.5 * gamma_R_0 * (1+Spin_polarization_R)
     g0pa_dn = 0.5 * gamma_R_0 * (1-Spin_polarization_R)
     g1pa_up = 0.5 * gamma_L_0 * (1+Spin_polarization_L)
     g1pa_dn = 0.5 * gamma_L_0 * (1-Spin_polarization_L)

     level_j: do j=1,Ndim
     level_u: do u=1,Ndim

!         Precalculate orbital overlaps          
          call Orbital_overlaps(lambda, j, u, g0pa_up, g0pa_dn, Ndim, LvlujaR, LjulvaR)
          call Orbital_overlaps(lambda, j, u, g1pa_up, g1pa_dn, Ndim, LvlujaL, LjulvaL)

          call ExtendedFermiIntegralBessel (Delta(j,u), frequency, bias_R, p_max-1, Temperature, &
                                             Cutoff, GammaC, N_int, fermiRB, ufermiRB)
          call ExtendedFermiIntegralBessel (Delta(j,u), frequency, bias_L, p_max-1, Temperature, &
                                             Cutoff, GammaC, N_int, fermiLB, ufermiLB)

          fourier_component: do n =-n_max,n_max
               n_index = n + n_max + 1

!              contribution of bessel functions
!              bessel_cont  = sum_p K*_{p-n} K_p  I(p)                              
!              ubessel_cont = sum_p K*_p K_{p+n}  uI(p)        
               call compute_bessel_contribution(K_R, fermiRB, ufermiRB, p_max, n, &
                                             bessel_contributionR, ubessel_contributionR)
               call compute_bessel_contribution(K_L, fermiLB, ufermiLB, p_max, n, &
                                             bessel_contributionL, ubessel_contributionL)

               level_v: do v=1,Ndim
               level_l: do l=1, Ndim
               
               G (v,l,j,u,1,n_index) = 0.5*( LvlujaR(v,l) * bessel_contributionR + &
                                             LjulvaR(v,l) * ubessel_contributionR)
               G (v,l,j,u,2,n_index) = 0.5*( LvlujaL(v,l) * bessel_contributionL + &
                                             LjulvaL(v,l) * ubessel_contributionL)
                                             
               enddo level_l
               enddo level_v
          enddo fourier_component
     enddo level_u
     enddo level_j
         
     return

     end subroutine rates_bessel 


     subroutine ratesC_bessel (Ndim, NFreq, Nbias, lambda, gamma_R_0, gamma_L_0,  &
         Spin_polarization_R, Spin_polarization_L, &
         p_max, B_R, B_L, Amplitude, frequency, bias_R, bias_L,&
         Temperature, Electrode,  GC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the QME rates
! time dependent
! for pulses that are either steps or cosine or a sequence of both
     implicit none
! Input:
     complex (qc), intent (in):: lambda (:,:,:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, Temperature, frequency, Amplitude
     real (q), intent (in):: Spin_polarization_R, Spin_polarization_L, B_L, B_R
     real (q), intent (in):: bias_R , bias_L 
     integer :: Ndim, NFreq, Nbias, p_max
     integer :: Electrode
! Output: the Rates called GC (:,:,:,:,:) here
     complex (qc) :: GC (:,:,:,:,:) ! for current
! Only used in this subroutine
     integer :: v, l, j, u, n, n_index
     complex (qc), dimension(2*p_max-1) :: fermiRB, fermiLB, ufermiRB, ufermiLB
     complex (qc), dimension(Ndim,Ndim) :: LvlujaL, LjulvaL, LvlujaR, LjulvaR
     complex (qc), dimension(2*p_max-1) :: K_L, K_R
     complex (qc) :: g0pa_up, g1pa_up, g0pa_dn, g1pa_dn
     complex (qc) :: bessel_contributionR, ubessel_contributionR
     complex (qc) :: bessel_contributionL, ubessel_contributionL
     
!    TODO: there shouuld be one one subroutine that returns G and GU for
!        each electrode and deal with it on the spot.
     
!    Calculate Contribution of Bessel functions
!    K(p) = J(p) + .5*A*(J(p-1)+ J(p+1))
     K_R = Bessel_K(B_R/frequency, Amplitude, p_max)
     K_L = Bessel_K(B_L/frequency, Amplitude, p_max)

     g0pa_up = 0.5 * Electrode     * gamma_R_0 * (1+Spin_polarization_R)
     g0pa_dn = 0.5 * Electrode     * gamma_R_0 * (1-Spin_polarization_R)
     g1pa_up = 0.5 * (1-Electrode) * gamma_L_0 * (1+Spin_polarization_L)
     g1pa_dn = 0.5 * (1-Electrode) * gamma_L_0 * (1-Spin_polarization_L)

     level_j: do j=1,Ndim
     level_u: do u=1,Ndim

!         Precalculate orbital overlaps
          call Orbital_overlaps(lambda, j, u, g0pa_up, g0pa_dn, Ndim, LvlujaR, LjulvaR)
          call Orbital_overlaps(lambda, j, u, g1pa_up, g1pa_dn, Ndim, LvlujaL, LjulvaL)
          
          call ExtendedFermiIntegralBessel (Delta(j,u), frequency, bias_R, p_max-1, Temperature, &
                                             Cutoff, GammaC, N_int, fermiRB, ufermiRB)
          call ExtendedFermiIntegralBessel (Delta(j,u), frequency, bias_L, p_max-1, Temperature, &
                                             Cutoff, GammaC, N_int, fermiLB, ufermiLB)

          fourier_component: do n =-n_max,n_max
               n_index = n + n_max + 1

!              contribution of bessel functions
!              sum_p K*_{p-n} K_p  I(p)                              
!              sum_p K*_p K_{p+n}  uI(p)        
               call compute_bessel_contribution(K_R, fermiRB, ufermiRB, p_max, n, &
                                             bessel_contributionR, ubessel_contributionR)
               call compute_bessel_contribution(K_L, fermiLB, ufermiLB, p_max, n, &
                                             bessel_contributionL, ubessel_contributionL)

               level_v: do v=1,Ndim
               level_l: do l=1, Ndim
                    
               GC (v,l,j,u,n_index) = 0.5* ( &
                    LvlujaR(v,l) * bessel_contributionR - &
                    LjulvaR(v,l) * ubessel_contributionR + &
                    LvlujaL(v,l) * bessel_contributionL - &
                    LjulvaL(v,l) * ubessel_contributionL) 
                                     
               enddo level_l
               enddo level_v
          enddo fourier_component
     enddo level_u
     enddo level_j
         
     return

     end subroutine ratesC_bessel 
!
! initial population for the propagation
!
     subroutine initial_population (Ndim, Ntime, Eigenvalues, rho, Temperature)
     implicit none
     integer, intent (in) :: Ndim, Ntime
     real (q), intent (in) :: Temperature
     real (q), intent (in) :: Eigenvalues (:)
     complex (qc), allocatable, intent (out) :: rho (:,:,:)
     real (q) :: Z
     integer :: i

     allocate (rho (Ndim, Ndim, Ntime))

     Z = 0._q; rho = zero

     do i =1, Ndim
     
     
#ifdef __UNDERANDOVER
     ! Underflow errors possible here
     if (abs(exp( - (Eigenvalues (i)-Eigenvalues (1))/Temperature )) < tiny(0.0_q)) then
       write(*,*) "Tiny values found for eigenvalue diff (integer, value): ", i, Eigenvalues (i)-Eigenvalues (1)
       write(*,*) "Temperature: ",Temperature
       write(*,*) "Value: ", exp( - (Eigenvalues (i)-Eigenvalues (1))/Temperature )
     end if
#endif
     Z = Z + exp ( - (Eigenvalues (i)-Eigenvalues (1))/Temperature )
     enddo


     do i = 1, Ndim
     rho (i,i,1) = rho (i,i,1) + exp ( - (Eigenvalues (i)-Eigenvalues (1))/Temperature )/Z
     enddo

!      write(*,*) Temperature, Z,Eigenvalues
!       
!      do i=1,ndim
!         write(*,*)
!         write(*,*) (rho(i,j,1),j=1,ndim)
!         write(*,*)  
!      enddo
     
     return

     end subroutine initial_population 
!
! Master equation solver
! RK4 method
!

!    TODO: implement a Bessel version
!          for now I can test the Bessel rates with the current function
!          and then implement the Bessel version of the RungeKutta
     
     subroutine RungeKutta (Ndim, Ntime, stept, Delta, rho, NFreq,  &
      lambda, gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, Nbias, &
      bias_R, bias_L, bias_time, Spin_polarization_R, Spin_polarization_L,  &
      t_seq, Amplitude, Freq_seq, Phase_seq,  &
      fermiR_a, fermiL_a, ufermiR_a, ufermiL_a,  &
      n_max, p_max, B_R, B_L, &
      N_int, GammaC, Cutoff, time, Temperature, use_bessel)

     implicit none
     real (q), intent (in) :: stept
     real (q), intent (in) :: Delta (:,:), time (:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, Temperature, Cutoff, GammaC
     real (q), intent (in):: bias_R (:) , bias_L (:)
     real (q) :: Spin_polarization_R, Spin_polarization_L
     complex (qc), intent (in):: lambda (:,:,:)
     logical, intent (in):: use_bessel
     integer :: Ndim, Ntime, NFreq, N_int, Nbias, n_max, p_max
     complex (qc), dimension(2) :: Pulse, PulseR, PulseL  !for time i and i+1
     real (q), intent (in):: B_R, B_L
     complex (qc), allocatable:: Gstatic (:,:,:,:,:,:)
     complex (qc), allocatable:: Gbess (:,:,:,:,:,:)
     complex (qc), allocatable:: G (:,:,:,:,:)
     complex (qc) :: rho (:,:,:)
     complex (qc) :: Temp1, Temp2 !buffers for changing names...
! Sequence of pulses
!    Gamma_1(t)=gamma_alpha_1*Amplitude (:)*cos(Freq_seq (:)*t+Phase_seq (:))*
!               Theta(t_seq (n)-t)*Theta(t-t_seq (n-1))
     integer, intent (in):: t_seq (:), bias_time (:) !sequence of pulses gives where time is
     real (q), intent (in):: Amplitude (:,:) ! sequence of pulses
     real (q), intent (in):: Freq_seq (:,:) ! sequence of pulses
     real (q), intent (in):: Phase_seq (:) ! sequence of pulses
     real (q), dimension(2) :: effec_Amplitude
! Computed in ExtendedFermiIntegral
     complex (qc) :: fermiR, fermiL, ufermiR, ufermiL
     complex (qc), allocatable:: fermiR_a(:,:,:), fermiL_a(:,:,:)
     complex (qc), allocatable:: ufermiR_a(:,:,:), ufermiL_a(:,:,:)
! internal
     integer :: l, j, i, u, v, m, n, n_index
     integer, dimension(2) :: na
     real (q) :: half
     complex (qc), allocatable :: k1(:), k2(:), k3(:), k4(:), D(:), P(:)

      
     allocate (k1(Ndim*Ndim), k2(Ndim*Ndim), k3(Ndim*Ndim), k4(Ndim*Ndim))
     allocate (D(Ndim*Ndim), P (Ndim*Ndim))
     allocate (G(Ndim,Ndim,Ndim,Ndim,2))
     allocate (fermiR_a(Ndim,Ndim,Nbias))
     allocate (fermiL_a(Ndim,Ndim,Nbias))
     allocate (ufermiR_a(Ndim,Ndim,Nbias))
     allocate (ufermiL_a(Ndim,Ndim,Nbias))
     
     if (use_bessel) then
          allocate (Gbess(Ndim,Ndim,Ndim,Ndim,2,n_max*2+1))
          ! TODO: the amplitude is wrong here
          effec_Amplitude(1) = Amplitude(1,1)*gamma_R_1/gamma_R_0
          effec_Amplitude(2) = Amplitude(1,1)*gamma_L_1/gamma_L_0     

          call rates_bessel(Ndim, NFreq, Nbias, lambda, gamma_R_0, gamma_L_0,  &
               Spin_polarization_R, Spin_polarization_L, p_max, B_R, B_L,  &
               effec_Amplitude, Freq_seq(1,1), bias_R(1), bias_L(1), &
               Temperature, Electrode, Gbess) 
     else
! Bias intergal
     allocate (Gstatic(Ndim,Ndim,Ndim,Ndim,2,Nbias))
     do n = 1, Nbias
! I11 and I21 from the Manual are honored here:
     do j=1,Ndim
     do u=1,Ndim
      call ExtendedFermiIntegral ( Delta (j,u), bias_R (n), Temperature, Cutoff, GammaC, N_int, fermiR)
      call ExtendedFermiIntegral ( Delta (j,u), bias_L (n), Temperature, Cutoff, GammaC, N_int,  fermiL)
      call ExtendeduFermiIntegral ( Delta (j,u), bias_R (n), Temperature, Cutoff, GammaC, N_int, ufermiR) 
      call ExtendeduFermiIntegral ( Delta (j,u), bias_L (n), Temperature, Cutoff, GammaC, N_int,  ufermiL)

      fermiR_a(j,u,n) = fermiR  / pi_d !important pi factor, see Manual
      fermiL_a(j,u,n) = fermiL  / pi_d
      ufermiR_a(j,u,n) = ufermiR  / pi_d
      ufermiL_a(j,u,n) = ufermiL  / pi_d

     enddo
     enddo
     enddo

     
       call rates (Ndim, NFreq, Ntime, Nbias, lambda, gamma_R_0, gamma_L_0,  &
         Spin_polarization_R, Spin_polarization_L, fermiR_a, fermiL_a,   &
         ufermiR_a, ufermiL_a, Temperature, Gstatic)
     endif

! loop on time
      call clock ('Entering time loop in RK after computing static rates', 2)
! Matyas: Name loops 
       time_loop: do i = 1, Ntime-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate rates for this time i, and  next one i+1
     if (use_bessel) then
          G(:,:,:,:,:) = (0,0) 
          fourire: do n=-n_max, n_max
          n_index = n+n_max+1      
          Pulse(1) = exp(-ui*cmplx(n,0)*(cmplx(time(i  ),0)*Freq_seq(1,1)+ Phase_seq(1)))
          Pulse(2) = exp(-ui*cmplx(n,0)*(cmplx(time(i+1),0)*Freq_seq(1,1)+ Phase_seq(1)))

          G(:,:,:,:,1) = G(:,:,:,:,1) + (Gbess(:,:,:,:,1,n_index)+Gbess(:,:,:,:,2,n_index))*Pulse(1)
          G(:,:,:,:,2) = G(:,:,:,:,2) + (Gbess(:,:,:,:,1,n_index)+Gbess(:,:,:,:,2,n_index))*Pulse(2)
          
          enddo fourire
     else
        na= t_seq(i:i+1) ! index of the pulse interval that contains time (i)

! Pulse sequence
          Pulse = (0,0)
        do m= 1, Nfreq
          !TODO: this could be writen in one line but not sure if it is more readable
          Pulse(1) = Pulse(1) + Amplitude (na(1),m)*cos(Freq_seq(na(1),m)*time(i)  +Phase_seq(na(1)))
          Pulse(2) = Pulse(2) + Amplitude (na(2),m)*cos(Freq_seq(na(2),m)*time(i+1)+Phase_seq(na(2)))
        enddo
     
        PulseR = ((1,0) + Pulse*gamma_R_1/gamma_R_0)**2
        PulseL = ((1,0) + Pulse*gamma_L_1/gamma_L_0)**2

! Add pulse on G, I first create a temporal value to add left and right electrode
! and then I crash the left and right G by the first and second times (ugly, I KNOW)
          n = bias_time (i) 
          G (:,:,:,:,1) =  Gstatic (:,:,:,:,1,n)*PulseR(1) + &
                         Gstatic (:,:,:,:,2,n)*PulseL(1) ! The driving is the ratio gamma_R_1/gamma_R_0

! Repeat for time i+1
! Pulse sequence
        n = bias_time (i+1)
        G (:,:,:,:,2) = Gstatic (:,:,:,:,1,n)*PulseR(2) + &
                        Gstatic (:,:,:,:,2,n)*PulseL(2) ! The driving is the ratio gamma_R_1/gamma_R_0
      endif
! G(:,:,:,:,1) is for time i and G(:,:,:,:,2) is for time i+1
! end generate rates for this time i, and i+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initial step
! Respects detailed balance because it is thermal
          n = 0
       do l = 1, Ndim
       do j = 1, Ndim
          n = n+1
       D (n) = rho (l,j,i)
       enddo
       enddo


       half = 0
       call matrix_rate (D, G, Delta, P, half)
       k1 = stept*P

       D = D + 0.5*k1
       half = 0.5_q
       call matrix_rate (D, G, Delta, P, half)
       k2 = stept*P
       
        n = 0
       do l = 1, Ndim
       do j = 1, Ndim
          n = n+1
       D (n) = rho (l,j,i)
       enddo
       enddo

       D = D + 0.5*k2
       half = 0.5_q
       call matrix_rate (D, G, Delta, P, half)
       k3 = stept*P
       
         n = 0
       do l = 1, Ndim
       do j = 1, Ndim
          n = n+1
       D (n) = rho (l,j,i)
       enddo
       enddo
       D = D + k3
       half = 1._q
       call matrix_rate (D, G, Delta, P, half)
       k4 = stept*P

       D = k1+2._q*k2+2._q*k3+k4

       
          n = 0
       do l = 1, Ndim
       do j = 1, Ndim
          n = n+1
       rho (l,j,i+1) = rho (l,j,i) + D (n) /6._q
       enddo
       enddo

       enddo time_loop

     deallocate (k1, k2, k3, k4)
     deallocate (D,P)
     if (use_bessel) then
       deallocate (Gbess)
     else
       deallocate (Gstatic)
     end if
     deallocate (G)

       return
     end subroutine RungeKutta

!
! The rate equation is written as D'=M.D
!
! it is:
!  
! 
     subroutine matrix_rate (D, G, Delta, P, half)
     implicit none
     complex (qc), intent (in):: G (:,:,:,:,:)
     complex (qc), intent (in) :: D (:)
     complex (qc), intent (out) :: P(:)
     real (q), intent(in) :: Delta (:,:)
     real (q) :: half
     complex (qc) :: rate
     integer ::  m
     integer :: n, l, j, r, s

          n = 0
       do l = 1, Ndim
       do j = 1, Ndim
          n = n+1

              P (n) = ui*Delta (l,j)*D(n)

! contribution on rho (v,u)
             m = 0
          do v = 1, Ndim
          do u = 1, Ndim
             m = m+1

              rate = half*(G (v,l,j,u,2)-G (v,l,j,u,1))+G (v,l,j,u,1)
              rate = rate + conjg(half*(G (u,j,l,v,2)-G (u,j,l,v,1))+G (u,j,l,v,1))

              P (n) = P(n) + rate * D(m)

! contribution on rho (l,u)
         
             r = (l-1)*Ndim + u

             rate = half*(G (j,v,v,u,2) - G (j,v,v,u,1)) + G (j,v,v,u,1)

             P(n) = P (n) - rate * D (r)
              
! contribution on rho (u, j)

            s = (u-1)*Ndim + j
              
            rate = half*(G (l,v,v,u,2) - G (l,v,v,u,1)) + G (l,v,v,u,1)

            P(n) = P (n) - conjg(rate) * D (s)

          enddo
          enddo


       enddo
       enddo

     return

     end subroutine matrix_rate
!
! Fermi occupation function
!
!    function Fermi (e, T)
!      implicit none
!      real (q) :: Fermi, e, T
!         if ( e > 0._q ) then
!            Fermi = exp(-e/T)/(exp(-e/T)+1._q)
!         else
!           Fermi = 1._q/(exp(e/T)+1._q)
!         endif
!      return
!    end function Fermi

!
! Fermi occupation function with limitors to avoid underflow/overflows
!
    function Fermi (e, T)
    
      implicit none
      real (q) :: Fermi
      real (q), intent(in) :: e, T
      real (q) :: beta
      
      Fermi = 0._q
      beta = 0._q
      
      beta = 1._q / T
      Fermi = beta * e
      if (e < -1000._q) then
      	Fermi = 0._q
      else
#ifdef __UNDERANDOVER
     ! Overflow errors possible here
     if (exp(Fermi) > huge(0.0_q)) then
       write(*,*) "Large Fermi value found for energy, temperature, value: ",e,T,exp(Fermi)
     end if
#endif
        Fermi = exp(Fermi)
      end if
      
      if (fermi > 1.0e30_q) then
        Fermi = 0._q
      else
        Fermi = 1._q/(1._q + Fermi)
      end if
    end function

! Calculation of energy integration of rates involving the Fermi function
      subroutine ExtendedFermiIntegral ( D, V, T, Cutoff, GammaC, N,  fermiA)
      implicit none
      real (q) :: D, V, T, Cutoff, GammaC
      real (q) :: e, step_e
      integer :: i, N
      complex (qc):: fermiA
!fermiA is Integral I11 of the Manual
! Trapeze-integration (the best among the better)

      step_e = 2*Cutoff/(N-1)
      e= -Cutoff
      fermiA=0.5*Fermi (e-V, T) / (e-D+ui*GammaC)
      
      do i = 2, N-1
      e= -Cutoff + (i-1)*step_e
      fermiA=fermiA+Fermi (e-V, T) / (e-D+ui*GammaC)
      enddo
      e = Cutoff
      fermiA=fermiA+0.5*Fermi (e-V, T) / (e-D+ui*GammaC)

      fermiA = step_e*ui*fermiA

      return
      end subroutine ExtendedFermiIntegral
! Calculation of energy integration of rates involving 1-Fermi function
      subroutine ExtendeduFermiIntegral ( D, V,  T, Cutoff, GammaC, N, ufermiA)
      implicit none
      real (q) :: D, V, T, Cutoff, GammaC
      real (q) :: e, step_e
      integer :: i, N
      complex (qc):: ufermiA
!ufermiA is Integral I21 of the Manual
! Trapeze-integration (the best among the better)

      step_e = 2*Cutoff/(N-1)
      e= -Cutoff 
      ufermiA=0.5*(1-Fermi (e-V, T)) / (e+D-ui*GammaC)

      do i = 2, N-1
      e= -Cutoff + (i-1)*step_e
      ufermiA=ufermiA+(1-Fermi (e-V, T)) / (e+D-ui*GammaC)
      enddo
      e = Cutoff
      ufermiA=ufermiA+0.5*(1-Fermi (e-V, T)) / (e+D-ui*GammaC)

      ufermiA = -step_e*ui*ufermiA

      return
      end subroutine ExtendeduFermiIntegral

    subroutine ExtendedFermiIntegralBessel ( D, frequency, V, p_max, T, Cutoff, GammaC, N,  fermiA, ufermiA)
     implicit none
     real (q) :: D, V, T, Cutoff, GammaC, frequency
     real (q) :: e, step_e
     integer :: i, N, p, p_max, p_ind
     complex (qc), dimension(2*p_max+1):: fermiA, ufermiA
     real (q), dimension(N) :: f
     complex (qc) :: denom, udenom
     ! I11_p = i/pi \int_{-Cutoff}^{Cutoff} dE f(E,V)/(E-D+p*frequency+ui\Gamma_C)
     ! I21_p = -i/pi \int_{-Cutoff}^{Cutoff} dE (1-f(E,V))/(E+D+p*frequency-ui\Gamma_C)
     ! f(E,V) = \frac{1}{\exp(\beta (E-V)) + 1} Fermi distribution with V as fermi level
     
     ! TODO: the normal integrals could be included in this by setting p_max = 0
     
     ! calculate all fermi contributions they are the same for all integrals
     step_e = 2*Cutoff/(N-1)/T ! rescaling with to units T = 1 
     e= -(Cutoff + V)/T
     
     fstep: do i = 1, N
          e = e + step_e
          f(i) = Fermi (e, 1._q)
     enddo fstep

     step_e = 2*Cutoff/(N-1) ! now in atomic units
     
     ploop: do p = -p_max, p_max
          p_ind = p+p_max+1
          
          ! The constant part of the denominator 
          denom =  - D + p*frequency + ui*GammaC
          udenom = D + p*frequency - ui*GammaC

          e = - Cutoff
          fermiA(p_ind) = 0.5*f(1)/(e+denom)
          ufermiA(p_ind) = 0.5*(1-f(1))/(e+udenom)   
          
          istep: do i = 2, N-1
               e= e + step_e
               fermiA(p_ind)=fermiA(p_ind) + f(i)/(e+denom)
               ufermiA(p_ind)=ufermiA(p_ind) + (1-f(i))/(e+udenom)
          enddo istep
          
          e = Cutoff
          fermiA(p_ind) = fermiA(p_ind) + 0.5*f(N)/(e+denom)
          ufermiA(p_ind) = ufermiA(p_ind) + 0.5*(1-f(N))/(e+udenom)

          fermiA(p_ind) = step_e*ui*fermiA(p_ind)/pi_d
          ufermiA(p_ind) = -step_e*ui*ufermiA(p_ind)/pi_d

     enddo ploop

     return
     end subroutine ExtendedFermiIntegralBessel

     function Bessel_K(z, Amplitude, p_max) result (K)
     implicit none
     real (q), intent (in) :: z, Amplitude
     integer, intent (in) :: p_max
     complex (qc), dimension(2*p_max-1) :: K 
     complex(qc), dimension(2*p_max+1) :: J
     integer :: p

          J(p_max+1:) = Bessel_JN(0, p_max, z)
          
          ! J(-p) = (-1)**p J(p) for p = 0,1,2,...
          negative_bessel : do p = 0, p_max -1
               J(p+1) = ((-1)**(p_max-p))*J(2*p_max+1-p)
          enddo negative_bessel
          
          ! K(p) = J(p) + .5*A*(J(p-1)+ J(p+1))
          K = J(2:2*p_max) + 0.5 * Amplitude * (J(1:2*p_max-1) + J(3:2+p_max+1))
          
     return
     end function Bessel_K

     subroutine Orbital_overlaps (lambda, j, u, gpa_up, gpa_dn, Ndim, overlapvluj, overlapjulv)
     implicit none
     integer, intent(in):: Ndim, j, u
     complex (qc), dimension(Ndim, Ndim, 2), intent (in) :: lambda 
     complex (qc), intent (in) :: gpa_up, gpa_dn
     complex (qc), dimension(Ndim, Ndim), intent (out) :: overlapvluj, overlapjulv
     integer :: v, l

          level_v: do v=1,Ndim
          level_l: do l=1, Ndim
               overlapvluj(v,l) = lambda (v,l,1)*conjg(lambda(u,j,1))*gpa_up+  &
                              lambda (v,l,2)*conjg(lambda(u,j,2))*gpa_dn
               overlapjulv(v,l) = lambda (j,u,1)*conjg(lambda(l,v,1))*gpa_up+  &
                              lambda (j,u,2)*conjg(lambda(l,v,2))*gpa_dn
          enddo level_l
          enddo level_v

     return
     end subroutine Orbital_overlaps

     subroutine compute_bessel_contribution(K, fermi, ufermi, p_max, n, result_bessel, result_ubessel)
          integer, intent(in) :: p_max, n
          complex (qc), intent(in), dimension(2*p_max-1) :: K
          complex (qc), intent(in), dimension(2*p_max-1) :: fermi, ufermi
          complex (qc), intent(out) :: result_bessel, result_ubessel

          integer :: p
          result_bessel = 0.0_qc 
          result_ubessel = 0.0_qc

          ! sum_p K*_{p-n} K_p  I(p)
          ! sum_p K*_p K_{p+n}  uI(p)
          ! from p=1-p_max to p=p_max-1-n this assumes that for |p|>p_max K_p = 0
          bessel: do p = 1, 2*p_max - 1 - n 
               result_bessel  = result_bessel  + conjg(K(p)) * K(p+n) * fermi(p+n)
               result_ubessel = result_ubessel + conjg(K(p)) * K(p+n) * ufermi(p)
          enddo bessel

          return
     end subroutine compute_bessel_contribution

end module QME
