program monte_carlo
      implicit none

      !-------------------------------------------------
      !variables for random number
      
      integer :: idum=-16021987
      real *8 :: r 
      real *8 :: ran1
      
      !-------------------------------------------------
      ! variables of Monte-carlo intialization
      
      integer, allocatable :: spin(:,:)
      integer :: i, j, L, a, b, c, d, M, N, k
      real *8 ::  E, J_is=1, T=2.0d0
      
      !-------------------------------------------------
      ! variables for Monte-carlo steps (equilibriation)

      integer :: nn, mm, niter, times
      real *8 :: Ei, Ef, dE, u

      !-------------------------------------------------
      ! Reading inputs

      print *, "Give the length of the lattice:"
      read *, L

 !     print *, "Give number of Monte-Carlo steps:"
 !     read *, niter
      
      !------------------------------------------------
      ! Opening files
 
      open(unit=16, file='T2_L20_ir.dat', action='write')
      open(unit=17, file='clt.dat', action='write')

      !------------------------------------------------
      ! Initializing the array
      
      allocate(spin(L,L))
      N=L*L
      
      do k=1,100000
      do i=1,L
        do j=1,L
         
        ! Setting random spin values

          !r = ran1(idum)
          !if (r .le. 0.50d0) then
           !  spin(j,i) = 1
          !else
           !  spin(j,i) = -1
          !end if
         
         ! Setting spin value to 1 or -1

          !spin(j,i) = 1
          !spin(j,i) =-1

          !Setting half block to a value 1 and other -1

         if(j .le. L/2) then
             spin(j,i) = 1
         else 
             spin(j,i) = -1
         end if
       
         end do
      end do
      

      M = 0      
      E = 0.0d0
    
      do i=1,L
        do j=1,L
            
            ! Neighbours of spin(j,i)

            a=j-1
            b=i+1
            c=j+1
            d=i-1

            ! Periodic Boundary conditions
            
            if (j==L) c=1
            if (j==1) a=L
            if (i==L) b=1
            if (i==1) d=L

            ! Magnetization and Energy calculation for this micro state

            M = M + spin(j,i)
            E = E - J_is * float((spin(j,i)*(spin(j,b)+spin(j,d)+spin(a,i)+spin(c,i))))
    
        end do
      end do
      write(17,*) M/real(N), E*0.50d0/real(N)
      end do

      print *, "Total magnetization:", M
      print *, "Moment per spin:    ", M/real(N)
      print *, ""
      print *, "Total energy:   ", E*0.50d0
      print *, "Energy per spin:", E*0.50d0/real(N)

      !-------------------------------------------------------
      ! Evolving for equilibrium
      
     
      do times=1, niter 
            
            do mm=1,L
                  do nn=1,L
                        
                        ! Choose a random lattice site

                        r = ran1(idum)
                        i = int(r*float(L)) + 1
                        
                        r = ran1(idum)
                        j = int(r*float(L)) + 1
                        
                        ! Neighbours of spin(j,i)

                        a=j-1
                        b=i+1
                        c=j+1
                        d=i-1

                        ! Periodic Boundary conditions
            
                        if (j==L) c=1
                        if (j==1) a=L
                        if (i==L) b=1
                        if (i==1) d=L

                        ! initial energy
                        Ei = E - J_is * float((spin(j,i)*(spin(j,b)+spin(j,d)+spin(a,i)+spin(c,i))))

                        ! FLIP
                            spin(j,i) = -spin(j,i)

                            Ef = E - J_is * float((spin(j,i)*(spin(j,b)+spin(j,d)+spin(a,i)+spin(c,i))))
                            dE = Ef-Ei

                        if(dE .lt. 0) then
                            E = E + dE
                            M = M + (2.0d0 * float(spin(j,i)))
                        else
                            ! Accepting with probability exp(-dE/kT)
                            u = exp(-dE/real(T))

                            !call random_number(r)
                            r = ran1(idum)
                            if (r<u) then 
                                E = E + dE
                                M = M + (2.0d0 * float(spin(j,i)))
                            else
                                spin(j,i) = -spin(j,i)
                            end if
                        end if

                   end do
             end do
             write(16,*) times, M/real(N), E/real(N)
       end do

end program monte_carlo


!--------------------------------------------------------------------------------
FUNCTION ran1(IDUM)
 
    implicit none

    !  RAN1 returns a unifom random deviate on the interval [0,1]

    INTEGER :: IDUM
    REAL*8 :: RAN2,ran1
    integer,parameter :: IM1=2147483563,IM2=2147483399
    integer,parameter :: IMM1=IM1-1,                            &
        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
        NTAB=32
    integer,parameter :: NDIV=1+IMM1/NTAB  
    real*8,parameter :: EPS=1.2e-7,RNMX=1.-EPS,AM=1./IM1
    INTEGER :: IDUM2,J,K,IV(NTAB),IY
    DATA IDUM2/123456789/, iv/NTAB*0/, iy/0/
         IF (IDUM.LE.0) THEN
           IDUM=MAX(-IDUM,1)
           IDUM2=IDUM
           DO  J=NTAB+8,1,-1
              K=IDUM/IQ1
              IDUM=IA1*(IDUM-K*IQ1)-K*IR1
              IF (IDUM.LT.0) IDUM=IDUM+IM1
              IF (J.LE.NTAB) IV(J)=IDUM
           end do
           IY=IV(1)
         ENDIF
         K=IDUM/IQ1
         IDUM=IA1*(IDUM-K*IQ1)-K*IR1
         IF (IDUM.LT.0) IDUM=IDUM+IM1
         K=IDUM2/IQ2
         IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
         IF (IDUM2.LT.0) IDUM2=IDUM2+IM2
         J=1+IY/NDIV
         IY=IV(J)-IDUM2
         IV(J)=IDUM
         IF(IY.LT.1)IY=IY+IMM1
         RAN1=MIN(AM*IY,RNMX)

END function ran1
