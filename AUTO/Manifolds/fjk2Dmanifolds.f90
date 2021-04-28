!----------------------------------------------------------------------
!----------------------------------------------------------------------
! fjk2Dmanifolds : Compute manifold of the periodic orbit s 
!           in the 2D reduction of the JFK model
!
! Jen Creaser Apr 2021
!----------------------------------------------------------------------
!----------------------------------------------------------------------

! Parameters:

!    PAR(1) = B_min
!    PAR(2) = I
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
!    PAR(26) = L  Day length
!    PAR(27) = L2Norm of first 4 dimensions
!    PAR(30)  = manifold comp parameter
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
    DADC = (1./15.)*k*G*alpha*(1. - n_infty)*(1. - 0.4*A1)*C1
    
    B1(1,1) = (PI/12.)*( mu*(1.0 - 4.*A1**2.) + DADA)
    B1(1,2) = - (PI/12.)*( Tsq**2. + k*G*alpha*(1. - 0.4*A1)*(1.-n_infty)) + DADC
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
    INTEGER  n, istart, iend, add_orbs

    CALL RHS(U,PAR,F,NDIM>4,B1,B2)
    
    L = PAR(26) ! day length
    day_adjust = L/12.0
    night_adjust = (24-L)/12.0
    
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
    
    add_orbs = 17       ! adds extra orbits
    DO  n=0,(add_orbs-1) 
        istart = 13 + n*4
        iend = istart+3        
        CALL RHS(U(istart:iend),PAR,F(istart:iend),.FALSE.,B1,B2)
        F(istart : istart+1) = F(istart : istart+1) * T * day_adjust
        F(istart+2 : istart+3 ) = F(istart+2 : istart+3) * T * night_adjust
    END DO

  
  END SUBROUTINE FUNC
!----------------------------------------------------------------------
  SUBROUTINE STPNT(NDIM,U,PAR,T) !--------- ----- 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T
  LOGICAL, SAVE :: ifirst = .TRUE.

  DOUBLE PRECISION PERIOD,B_min,lux,mu,taux,k,delta,L,xx
  DOUBLE PRECISION v1(2),v2(2),lambda(2),a,b,c,d,modL1,modL2
  
    ! Parameter values for the starting orbit in .dat file
    B_min     = 0.0 ! don't use this as bif parameter
    lux          = 50.0
    mu          = 0.23
    taux        = 24.2
    k             = 0.55
    PERIOD = 12.0
    L             = 12.0 ! day length

    PAR(1:5) = (/B_min,lux,mu,taux,k/)    
    PAR(7) = 0.0 ! delta
    PAR(9:10) = (/0.0,0.0/) ! gap in orbit
    PAR(11) = PERIOD
    ! start of orbit 
    PAR(13:14) = (/0.99141230,   -0.41989485 /)
    ! start of manifold orbit 
    PAR(15:16) = PAR(13:14) 
    ! end of orbit 
    PAR(17:18) = PAR(13:14) 
 
    ! intialise eigenvalues  (real only)
    lambda(1:2) = (/0.81171124351676249,0.15296427437740079/)
    PAR(20:21) = lambda(1:2)     
    ! intialise eigenvectors 
    v1(1:2) = (/-0.3267, -0.9451/)
    PAR(22:23) = v1(1:2)   
    v2(1:2) = (/-0.9549, 0.2971/)
    PAR(24:25) = v2(1:2) 

    PAR(26) = L ! day length
    PAR(28:29) = (/0.0,0.0/)  ! fixing start and end on PO
    
    PAR(27) =  2.28792 ! starting L2 Norm value
        
  END SUBROUTINE STPNT
!----------------------------------------------------------------------

  SUBROUTINE PVLS
  END SUBROUTINE PVLS

!----------------------------------------------------------------------
  SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
! ---------- ----- 

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
    DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
    DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
    DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
    LOGICAL, SAVE :: ifirst = .TRUE.

    DOUBLE PRECISION lambda(2),v1(2),v2(2),p(2),delta
    DOUBLE PRECISION M(2,2), V(2,2), W(2,2)
    INTEGER istart,add_orbits,n
        
    ! Periodic boundary conditions for periodic orbit
    FB(1:2)= U0(1:2) - U1(3:4)
    FB(3:4)= U0(3:4) - U1(1:2)
    
    IF (NBC==4) RETURN   

    FB(5:8) = U0(5:8) - (/1,0,0,1/)
    FB(9:12) = U0(9:12) - (/1,0,0,1/)
    
    FB(13:14) = U0(1:2) - PAR(13:14) 
    
    ! ------------- Change from here to add orbits -------------------
    ! Periodic boundary conditions for manifold orbit 
    add_orbits = 17 ! adds extra orbits 
    
    istart = 15
    DO  n=1,(add_orbits*2-1) 
        FB(istart:istart+1)= U0(istart:istart+1) - U1(istart-2:istart-1)
        istart = istart + 2
    END DO
    FB(istart:istart+1)= U0(13:14)  - U1(istart-2:istart-1) - PAR(9:10) ! can separate
        
    ! follow end points as parameters
    FB(istart+2:istart+3) = U0(13:14) - PAR(15:16) ! follow start of trajectory1
    FB(istart+4:istart+5) = U1(istart-2:istart-1) - PAR(17:18) ! follow end of trajectory2

    IF (NBC==istart+5) RETURN ! 15 + n*2 +5
    
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
    lambda(1) =PAR(20)
    lambda(2)= PAR(21)
    
    FB(istart+6:istart+7)= MATMUL(M,PAR(22:23))-lambda(1)*PAR(22:23) ! eigenvector condition
    v1(1:2) = PAR(22:23)                                                     ! the eigenvector
    FB(istart+8)= sqrt(DOT_PRODUCT(v1,v1)) - 1.0                  ! norm the vector
    
    FB(istart+9:istart+10)= MATMUL(M,PAR(24:25))-lambda(2)*PAR(24:25)
    v2(1:2) = PAR(24:25) 
    FB(istart+11)= sqrt(DOT_PRODUCT(v2,v2)) - 1.0            
        
    IF (NBC==istart+11) RETURN
    
    p(1:2) = U0(1:2) ! starting point
    delta = PAR(7) 

    ! write values to output
    IF(ifirst)THEN
       WRITE(11,*) "I val:", PAR(2)
       WRITE(11,*) "start point: (", p(1), p(2), ")"
       WRITE(11,*) "lambda_1:", lambda(1),"  lambda_2:", lambda(2)
       WRITE(11,*) "v_1: (", v1(1), v1(2), ")", " v_2: (", v2(1), v2(2), ")"
       ifirst=.FALSE.
    ENDIF

    ! Manifold comp 
    IF (PAR(30)==1.0) THEN ! input to say which man to compute.
        IF (PAR(20)<1.0) THEN       
            FB(istart+12:istart+13) = U1(79:80)  - (p(1:2) + delta*v1(1:2)) 
        ELSE 
            FB(istart+12:istart+13) = U0(13:14)   - (p(1:2) + delta*v1(1:2)) 
        ENDIF
    ELSE
        IF (PAR(21)<1.0) THEN       
            FB(istart+12:istart+13) = U1(79:80) - (p(1:2) + delta*v2(1:2)) 
        ELSE 
            FB(istart+12:istart+13) = U0(13:14) - (p(1:2) + delta*v2(1:2)) 
        ENDIF
    ENDIF
        
  END SUBROUTINE BCND
!----------------------------------------------------------------------
  SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
!      ---------- ----

       IMPLICIT NONE
       INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
       DOUBLE PRECISION, INTENT(IN) :: PAR(*)
       DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
       DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
       DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

        FI(1)=( U(1)*U(1) +U(2)*U(2) +U(3)*U(3) +U(4)*U(4) )-PAR(27)

  END SUBROUTINE ICND
!----------------------------------------------------------------------
  SUBROUTINE FOPT 
  END SUBROUTINE FOPT
!----------------------------------------------------------------------
