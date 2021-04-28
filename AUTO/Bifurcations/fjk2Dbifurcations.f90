!----------------------------------------------------------------------
!----------------------------------------------------------------------
! fjk2Dbifurcations : Compute bifurcations of the periodic 
!                               orbits in the 2D reduction of the FJK 
!                               model
!
! Jen Creaser Apr 2021
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Parameters:

!    PAR(1) = B_min
!    PAR(2) = B_max
!    PAR(3) = mu   
!    PAR(4) = taux  
!    PAR(5) = k 

!    PAR(7) = delta

!    PAR(9:10) =  gap in orbit
!    PAR(11) = PERIOD

!    PAR(13:14) = start of orbit 
!    PAR(15:16) = start of manifold orbit 
!    PAR(17:18) = end of manifold orbit
!   PAR(19) : dummy for confirming po

!    PAR(20:21) = eigenvalues
!    PAR(22:23) = eigenvector 1
!    PAR(24:25) = eigenvector 2
!    PAR(26:27) = imaginary part of eigenvalues

!    PAR(28) = Neimark Sacker condition 
!    PAR(29) = Neimark Sacker condition 2?
!    PAR(30) = N day length
!    PAR(31) = L2Norm
!----------------------------------------------------------------------
!----------------------------------------------------------------------

  SUBROUTINE RHS(U,PAR,F,JAC,B1,B2)

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: U(4), PAR(*)
    LOGICAL, INTENT(IN) :: JAC
    DOUBLE PRECISION, INTENT(OUT) :: F(4), B1(2,2),B2(2,2)

    DOUBLE PRECISION B_min,B_max,mu,taux,k,tShift,PI,Tsq
    DOUBLE PRECISION A1,C1,A2,C2,D1,D2,DADA,DADC
    DOUBLE PRECISION alpha, G, a0,lux,lux0,p,beta,n_infty
    
    B_min   = PAR(1)
    lux        = PAR(2) ! bif para instead of B_max
    mu        = PAR(3)
    taux      = PAR(4)
    k           = PAR(5)
    tShift     = PAR(6)
 
    PI = 4.D0*DATAN(1.D0) 

    A1 = U(1)
    C1 = U(2)
    A2 = U(3)
    C2 = U(4)
 
    a0 = 0.05 
    p = 0.5
    G = 33.75
    lux0 = 9500
    beta = 0.0075
 
    alpha = a0*((lux/lux0)**p)
    n_infty = alpha/(alpha+beta)
    Tsq = 24./(0.99669*taux)
   
    B_max = G*alpha*(1. - 0.4*C1)*(1. - 0.4*A1)*(1.-n_infty) ! only in A1 C1 system
    
    D1 = C1*( Tsq**2.+k*B_max)
    D2 = C2*( Tsq**2.) ! B_min = 0

    F(1) = (PI/12.)*( mu*(A1-(4./3.)*A1**3.) - D1 )
    F(2) = (PI/12.)*( A1+B_max)
    F(3) = (PI/12.)*( mu*(A2-(4./3.)*A2**3.) - D2 )
    F(4) = (PI/12.)*( A2)

    DADA = 0.4*C1*k*G*alpha*(1.-n_infty)*(1. - 0.4*C1) 
    DADC = 0.8*k*G*alpha*(1. - n_infty)*(1. - 0.4*A1)*C1
    
    B1(1,1) = (PI/12.)*( mu*(1.0 - 4.*A1**2.) + DADA)
    B1(1,2) = - (PI/12.)*( Tsq**2. + k*G*alpha*(1. - 0.4*A1)*(1.-n_infty) - DADC)
    B1(2,1) = (PI/12.)*( 1.0 - 0.4*G*alpha*(1. - 0.4*C1)*(1. - n_infty) )
    B1(2,2) = - (PI/12.)*(0.4*G*alpha*(1. - 0.4*A1)*(1. - n_infty) )

    B2(1,1) = (PI/12.)*( mu*(1.0 - 4.*A2**2.) )
    B2(1,2) = - (PI/12.)*( Tsq**2.) 
    B2(2,1) = PI/12.
    B2(2,2) = 0.


  END SUBROUTINE RHS
!----------------------------------------------------------------------
  SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
! ---------- --- 

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
    DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
    DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

    DOUBLE PRECISION T, B1(2,2), B2(2,2),day_adjust,night_adjust,L
    DOUBLE PRECISION A1(2,2), A2(2,2), C1(2,2), C2(2,2)

    CALL RHS(U,PAR,F,NDIM>4,B1,B2)
    
    L = PAR(30) ! day length
    day_adjust = L/12.0;
    night_adjust = (24-L)/12.0;
    
    ! the orbit
    T = PAR(11) !T
    F(1:2) = F(1:2) * T * day_adjust
    F(3:4) = F(3:4) * T * night_adjust
    
    IF (NDIM==4) RETURN   

    ! Variational problem for Bmax
    A1(1,1) = U(5)
    A1(2,1) = U(6)
    A1(1,2) = U(7)
    A1(2,2) = U(8)
    C1 = MATMUL(B1,A1)
    F(5) = C1(1,1) * T * day_adjust
    F(6) = C1(2,1) * T * day_adjust
    F(7) = C1(1,2) * T * day_adjust
    F(8) = C1(2,2) * T * day_adjust

    ! Variational problem for Bmin
    A2(1,1) = U(9)
    A2(2,1) = U(10)
    A2(1,2) = U(11)
    A2(2,2) = U(12)
    C2 = MATMUL(B2,A2)
    F(9)   = C2(1,1) * T * night_adjust
    F(10) = C2(2,1) * T * night_adjust
    F(11) = C2(1,2) * T * night_adjust
    F(12) = C2(2,2) * T * night_adjust
    
    IF (NDIM==12) RETURN    


  END SUBROUTINE FUNC
!----------------------------------------------------------------------
  SUBROUTINE STPNT(NDIM,U,PAR,T)
  !--------- -----
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM)
  DOUBLE PRECISION,INTENT(INOUT) :: PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T
  LOGICAL, SAVE :: ifirst = .TRUE.

  DOUBLE PRECISION PERIOD,B_min,lux,mu,taux,k,delta,L,xx
  DOUBLE PRECISION v1(2),v2(2),lambda(2),a,b,c,d,modL1,modL2
  
    ! Parameter values for the starting orbit in usorbit.dat
    B_min     = 0.0 ! don't use this as bif parameter
    lux          = 50.0
    mu          = 0.23
    taux        = 24.2
    k             = 0.55
    PERIOD = 12.0
    L             = 12.0 ! day length
    
    PAR(1:5) = (/B_min,lux,mu,taux,k/)
    
    PAR(7) = 0.0 ! delta
    PAR(8) = 1.0  ! h: length of vector
    PAR(9:10) = (/0.0,0.0/) ! gap in orbit
    PAR(11) = PERIOD
    
    ! start of orbit 
    PAR(13:14) = (/0.99141230,   -0.41989485 /)
    
    ! start of manifold orbit 
    PAR(15:16) = PAR(13:14) 
    ! end of orbit 
    PAR(17:18) = PAR(13:14) 

    PAR(19) = 0.0 ! dummy 
    
    ! intialise eigenvalues   
    lambda(1:2) = (/0.8117112, 0.152964/)
    PAR(20:21) = lambda(1:2) 
    PAR(26:27) = (/0.0,0.0/)
    
    modL1 = SQRT(PAR(20)*PAR(20) + PAR(26)*PAR(26))   
    modL2 = SQRT(PAR(21)*PAR(21) + PAR(27)*PAR(27))
    PAR(28) = modL1 ! NS
    
    a =  1.0000000000000000       
    b = -0.96467551789416328       
    c = 0.12416282136851922      
    d = b*b - 4*a*c
    PAR(25) = d  
    
    PAR(29) = (modL1-1.0)+(modL2-1.0)
    
    PAR(30) = L ! day length start
    PAR(31) = 2.28792! starting L2 Norm value
    
    IF(ifirst)THEN
            WRITE(11,*) "STPT values"
            WRITE(11,*) "lambda1_real:",PAR(20),"lambda1_imag:",PAR(26),"lambda2_real:",PAR(21),"lambda2_imag:",PAR(27)
            WRITE(11,*) "L2Norm:",PAR(31)
         ifirst=.FALSE.
    ENDIF
 
    
  END SUBROUTINE STPNT
!----------------------------------------------------------------------
  SUBROUTINE PVLS(NDIM,U,PAR)
  !--------- ----

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
    DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
    
    DOUBLE PRECISION, EXTERNAL :: GETP

    DOUBLE PRECISION M(2,2), V(2,2), W(2,2),a,b,c,d,modL2,modL1

    INTEGER i, NBC
    LOGICAL, SAVE :: FIRST = .TRUE.
    LOGICAL, SAVE :: ifirst = .TRUE.
    
    IF (NDIM>4) THEN 
    
        ! Monodromy matrix
        V(1,1) = GETP("BV1",5,U)
        V(2,1) =GETP("BV1",6,U)
        V(1,2) = GETP("BV1",7,U)
        V(2,2) = GETP("BV1",8,U)
        W(1,1) = GETP("BV1",9,U)
        W(2,1) = GETP("BV1",10,U)
        W(1,2) = GETP("BV1",11,U)
        W(2,2) = GETP("BV1",12,U)
        M = MATMUL(V,W)
        
        ! here we are only going to compute the eigenvalues explicitly and follow them
        a = 1.0
        b = -(M(1,1) + M(2,2))
        c = (M(1,1)*M(2,2)) - (M(2,1)*M(1,2))
    
        d = b*b - 4*a*c
        
        IF (d>0) THEN 
            PAR(20) = (-b+SQRT(d))/(2*a)       ! Real part
            PAR(26) = 0.0                                 ! Imag part   
            PAR(21) = (-b-SQRT(d))/(2*a)        ! Real part    
            PAR(27) = 0.0                                 ! Imag part
        ELSE
            PAR(20)  = -b/(2*a)                         ! Real part
            PAR(26) = (SQRT(-d))/(2*a)           ! Imag part   
            PAR(21) = -b/(2*a)                          ! Real part    
            PAR(27) = -(SQRT(-d))/(2*a)          ! Imag part
        ENDIF
        
        modL1 = SQRT(PAR(20)*PAR(20) + PAR(26)*PAR(26))   
        modL2 = SQRT(PAR(21)*PAR(21) + PAR(27)*PAR(27))
        PAR(28) = modL1 ! NS
        
        PAR(29) =   (modL1-1.0)+(modL2-1.0) ! NS for real
    
        IF(ifirst)THEN
            WRITE(12,*) "PVLS values"
            WRITE(12,*) "lambda1_real:",PAR(20),"lambda1_imag:",PAR(26),"lambda2_real:",PAR(21),"lambda2_imag:",PAR(27)
            WRITE(12,*) "a = ", a, " b = ",b," c = ",c 
            ifirst=.FALSE.
        ENDIF
    ENDIF
    
  END SUBROUTINE PVLS
!----------------------------------------------------------------------
  SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
! ---------- ----- 

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
    DOUBLE PRECISION, INTENT(IN) :: U0(NDIM), U1(NDIM)
    DOUBLE PRECISION, INTENT(IN) :: PAR(*)
    DOUBLE PRECISION,  INTENT(OUT) :: FB(NBC)
    DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

    DOUBLE PRECISION M(2,2), V(2,2), W(2,2),a,b,c,d,L(2),modL1,modL2
        
    ! Periodic boundary conditions for periodic orbit
    FB(1:2)= U0(1:2) - U1(3:4)
    FB(3:4)= U0(3:4) - U1(1:2)
    
    IF (NBC==4) RETURN

    FB(5:8) = U0(5:8) - (/1,0,0,1/)
    FB(9:12) = U0(9:12) - (/1,0,0,1/)
    
    FB(13:14) = U0(1:2) - PAR(13:14) 
    
    IF (NBC==14) RETURN

    ! Monodromy matrix
    V(1,1) = U1(5)
    V(2,1) = U1(6)
    V(1,2) = U1(7)
    V(2,2) = U1(8)
    W(1,1) = U1(9)
    W(2,1) = U1(10)
    W(1,2) = U1(11)
    W(2,2) = U1(12)
    M = MATMUL(V,W)

    ! Eigenvector continuation
    L(1) =PAR(20)
    L(2)= PAR(21)
        
    ! Eigenvector continuation
    ! here we compute the eigenvalues explicitly and follow them
    a = 1.0
    b = -(M(1,1) + M(2,2))
    c = (M(1,1)*M(2,2)) - (M(2,1)*M(1,2))
    
    d = b*b - 4*a*c
    
    IF (d>0) THEN ! real
        FB(15) = ((-b+SQRT(d))/(2*a)) - PAR(20) ! Real part
        FB(16) = PAR(26) ! Imag part   
        FB(17) = ((-b-SQRT(d))/(2*a)) - PAR(21) ! Real part    
        FB(18) = PAR(27)  ! Imag part
    ELSE
        FB(15) = (-b/(2*a)) - PAR(20) ! Real part
        FB(16) = ((SQRT(-d))/(2*a)) - PAR(26) ! Imag part   
        FB(17) = (-b/(2*a)) - PAR(21) ! Real part    
        FB(18) = (-(SQRT(-d))/(2*a))  - PAR(27)  ! Imag part
    ENDIF
    
    modL1 =  SQRT(PAR(20)*PAR(20) + PAR(26)*PAR(26))   
    modL2 = SQRT(PAR(21)*PAR(21) + PAR(27)*PAR(27))

    FB(19) = modL1 - PAR(28) ! NS
   
  END SUBROUTINE BCND
!----------------------------------------------------------------------
  SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
!      ---------- ----

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
       DOUBLE PRECISION, INTENT(IN) :: PAR(*)
       DOUBLE PRECISION, INTENT(IN) :: U(NDIM),UOLD(NDIM),UDOT(NDIM),UPOLD(NDIM)
       DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
       DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

        FI(1)=( U(1)*U(1) +U(2)*U(2) +U(3)*U(3) +U(4)*U(4) )-PAR(31)

  END SUBROUTINE ICND
!----------------------------------------------------------------------
  SUBROUTINE FOPT 
  END SUBROUTINE FOPT
!----------------------------------------------------------------------
