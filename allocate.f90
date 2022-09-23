program allocate_random_spins
    implicit none
    
    !----------------------------------------------------
    !variables for random number

    integer :: idum=-16021987
    real *8 :: r1, r2, r3, R,M
    real *8 :: ran1

    !----------------------------------------------------
    integer :: i, j, k
    real *8, allocatable :: spin(:,:,:)
    integer :: L, N
    !----------------------------------------------------

    print *, "Give the length of lattice:"
    read *, L

    open(unit=16, file='spins.dat', action='write')

    allocate(spin(L,L,3))
    N = L*L

    do i=1,L
        j=1
        do while (j .le. L)
            call rand_uniform(-1,1,r1)
            call rand_uniform(-1,1,r2)
            call rand_uniform(-1,1,r3)
            R = sqrt(r1**2 + r2**2 + r3**2)

            if (R .le. 1) then
                spin(i,j,1) = r1/R
                spin(i,j,2) = r2/R
                spin(i,j,3) = r3/R
                M = sqrt(spin(i,j,1)**2+spin(i,j,2)**2+spin(i,j,3)**2)
                write(16,*) i , j, R, spin(i,j,1), spin(i,j,2), spin(i,j,3), M 
                j = j+1
            end if

        end do
    end do
    
    do i=1,L
        do j=1,L
            do k=1,3
                print *, i, j, k, spin(i,j,k)
            end do
        end do
    end do

end program allocate_random_spins

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
!-------------------------------------------------------------------------

!Subroutine for uniform random numbers
subroutine rand_uniform(a,b,x)
    implicit none
    !integer :: idum=-16021987
    !real *8 :: ran1
    real *8 :: u, r
    integer , intent(in) :: a, b
    real *8, intent(out) :: x
    !r = ran1(idum)
    call random_number(r)
    u = 1-r
    x = (b-a)*u +a
end subroutine rand_uniform

