module Transport
Use declarations
Use QME
CONTAINS
      subroutine Current (Ndim, Ntime, rho, NFreq, Nbias, bias_time,  &
            lambda, gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, &
            Spin_polarization_R, Spin_polarization_L,  &
            t_seq, Amplitude, Freq_seq, Phase_seq,  &
            fermiR_a, fermiL_a, ufermiR_a, ufermiL_a,  &
            Temperature, Electrode, curr)

      implicit none
      integer, intent (in) :: Ndim, Ntime, NFreq, Nbias, Electrode
      complex (qc), intent (in) :: rho (:,:,:)
      real (q), intent (in):: gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, Temperature
      real (q), intent (in):: Spin_polarization_R, Spin_polarization_L
      complex (qc), intent (in):: lambda (:,:,:)
      complex (qc), allocatable ::  GC (:,:,:,:,:)
! Sequence of pulses
!    Gamma_1(t)=gamma_alpha_1*Amplitude (:)*cos(Freq_seq (:)*t+Phase_seq (:))*
!               Theta(t_seq (n)-t)*Theta(t-t_seq (n-1))
      integer, intent (in):: t_seq (:), bias_time(:) !sequence of pulses gives where time is
      real (q), intent (in):: Amplitude (:,:) ! sequence of pulses
      real (q), intent (in):: Freq_seq (:,:) ! sequence of pulses
      real (q), intent (in):: Phase_seq (:) ! sequence of pulses
      real (q) :: Pulse  !for time i
! Computed in ExtendedFermiIntegral called in RungeKutta
      complex (qc) :: fermiR, fermiL, ufermiR, ufermiL
      complex (qc) :: fermiR_a(:,:,:), fermiL_a(:,:,:)
      complex (qc) :: ufermiR_a(:,:,:), ufermiL_a(:,:,:)
      real (q), allocatable :: curr (:)
! internal
      integer :: l,j,u,n,m,np

      allocate (curr (Ntime))
      allocate (GC(Ndim,Ndim,Ndim,Ndim,Nbias))


            call  ratesC (Ndim, NFreq, Nbias,lambda, gamma_R_0, gamma_L_0,  &
            Spin_polarization_R, Spin_polarization_L, fermiR_a, fermiL_a, ufermiR_a, ufermiL_a, &
            Temperature, Electrode,  GC)

      curr = 0._q

      do i = 1, Ntime

      n= t_seq(i) ! index of the pulse interval that contains time (i)
      np= bias_time (i) ! bias pulse
! Pulse sequence
            Pulse = 0._q
      do m= 1, Nfreq
            Pulse = Pulse + Amplitude (n,m)*cos(Freq_seq(n,m)*time(i)+Phase_seq(n))
      enddo

            do l = 1, Ndim
            do u = 1, Ndim
            do j = 1, Ndim
            
!           TODO: We calculate the entire GC tesnor but only use G(l,j,j,u,np)
!                 we could reduce the scaling of the code from Ndim^4 to Ndim^3 
!                 for the GC tensor. (irrelevant if we combine G and GC)
            curr (i) = curr (i) +    &
                  real(rho (l,u,i)*GC(l,j,j,u,np)+ &
                  conjg(rho (l,u,i)*GC(l,j,j,u,np)))*  &
                  (1+Pulse*((1-Electrode)*gamma_L_1/gamma_L_0+Electrode*gamma_R_1/gamma_R_0))**2
            enddo
            enddo
            enddo

      enddo 

! test Gamma C
!        do u = 1, Ndim
!        do j = 1, Ndim
!     write (124,*) u,j,j,u,real(GC(u,j,j,u,np)) 
!        enddo
!        enddo


      deallocate (GC)
      
      return 
      end subroutine Current


      subroutine Current_bessel (Ndim, Ntime, rho, NFreq, Nbias, bias_time,  &
         lambda, gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, &
         Spin_polarization_R, Spin_polarization_L,  &
         t_seq, Amplitude, Freq_seq, Phase_seq,  &
         n_max, p_max, B_R, B_L,  bias_R, bias_L, &
         Temperature, Electrode, curr)

     implicit none
     integer, intent (in) :: Ndim, Ntime, NFreq, Nbias, Electrode, p_max, n_max
     complex (qc), intent (in) :: rho (:,:,:)
     real (q), intent (in):: gamma_R_0, gamma_L_0, gamma_R_1, gamma_L_1, Temperature
     real (q), intent (in):: Spin_polarization_R, Spin_polarization_L, B_R, B_L
     real (q), intent (in):: bias_R , bias_L 
     complex (qc), intent (in):: lambda (:,:,:)
! Sequence of pulses
!    Gamma_1(t)=gamma_alpha_1*Amplitude (:)*cos(Freq_seq (:)*t+Phase_seq (:))*
!               Theta(t_seq (n)-t)*Theta(t-t_seq (n-1))
     integer, intent (in):: t_seq (:), bias_time(:) !sequence of pulses gives where time is
     real (q), intent (in):: Amplitude (:,:) ! sequence of pulses
     real (q), intent (in):: Freq_seq (:,:) ! sequence of pulses
     real (q), intent (in):: Phase_seq (:) ! sequence of pulses
     real (q) :: Pulse  !for time i
! Computed in ExtendedFermiIntegral called in RungeKutta
     complex (qc), allocatable :: fermiR_a(:,:,:), fermiL_a(:,:,:)
     complex (qc), allocatable  :: ufermiR_a(:,:,:), ufermiL_a(:,:,:)
     real (q), allocatable :: curr (:)
! internal
      integer :: l,j,u,n,m,np,i,n_index
      complex (qc), allocatable ::  GC (:,:,:,:,:)
      complex (qc) :: tdep, exponent
      real (q) :: effec_Amplitude

      allocate (curr (Ntime))
      allocate (GC(Ndim,Ndim,Ndim,Ndim,2*n_max+1))
      
      effec_Amplitude = Amplitude(1,1)*((1-Electrode)*gamma_L_1/gamma_L_0 + &
                                          Electrode*gamma_R_1/gamma_R_0)

      call ratesC_bessel (Ndim, NFreq, Nbias,lambda, gamma_R_0, gamma_L_0,  &
            Spin_polarization_R, Spin_polarization_L, &
            p_max, B_R, B_L, effec_Amplitude, Freq_seq(1,1), bias_R, bias_L, &
            Temperature, Electrode,  GC)       

      curr = 0._q
      fourire: do n=-n_max, n_max
      timeloop: do i = 1, Ntime
      
            n_index = n+n_max+1      
            exponent = -ui*cmplx(n,0)*(cmplx(time(i),0)*Freq_seq(1,1)+ Phase_seq(1))
            tdep = exp(exponent)
            
            level_l: do l = 1, Ndim
            level_u: do u = 1, Ndim
            level_j: do j = 1, Ndim
            
!           TODO: why is A+conj(A) and then we take the rela part? why not just 2*real(A)
            curr (i) = curr (i) +    &
                  real(rho (l,u,i)*GC(l,j,j,u,n_index)*tdep +  &
                  conjg(rho (l,u,i)*GC(l,j,j,u,n_index)*tdep))
            
            enddo level_j
            enddo level_u                 
            enddo level_l
                  
      
      enddo timeloop 
      enddo fourire
! test Gamma C
!        do u = 1, Ndim
!        do j = 1, Ndim
!     write (124,*) u,j,j,u,real(GC(u,j,j,u,np)) 
!        enddo
!        enddo


      deallocate (GC)
      
      return 
      end subroutine Current_bessel 

end module Transport
