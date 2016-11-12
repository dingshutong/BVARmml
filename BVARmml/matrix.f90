! ------- Fortran modules: matrix --------------
! ------- Basic Linear Algebra -----------------
! ------- Author: Shutong Ding -----------------
! ------- Last updated: 2016-11-12 -------------
! ------- Dependency: IMSL ---------------------
    

module matrix

include 'link_fnl_static.h'!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries, only : chfac, linrt, linrg, mxtxf

implicit none 



interface DMEAN
    module procedure DMEAN_V, DMEAN_M
end interface

interface DREP
    module procedure DREP_V, DREP_M
end interface

interface DSD
    module procedure DSD_V, DSD_M
end interface

interface DBIND
    module procedure DBIND_VM, DBIND_MV, DBIND_MM, DBIND_VV
end interface

interface DDIAG
    module procedure DDIAG_A, DDIAG_V
end interface

interface PHIX
    module procedure PHIX_M, PHIX_V
END INTERFACE

interface PHIVX
    module procedure PHIVX_M, PHIVX_V
END INTERFACE

interface DLAG
    module procedure DLAG_M, DLAG_V
END INTERFACE

CONTAINS

!=============================================================================

SUBROUTINE DMEAN_V(X,MU)
! Computer the mean value of a vector X
! Aug:  X ---- Input vector 
!       MU---- OUTPUT, mean(X)
REAL*8, DIMENSION(:) :: X
REAL*8 :: MU
MU=SUM(X)/SIZE(X)
END SUBROUTINE 

! ====================================================
SUBROUTINE DMEAN_M(X,MU,D)
! Computer the mean value of matrix with 3 ways (col, row, whole)
! Aug:  X ---- input matrix
!      MU ---- OUTPUT, a mean vector or a mean value 
!      D  ---- Optional, present(D).eq. mean for the whole matrix; 
!              o.w., D=1, mean in row way, D=2, mean in col way
INTEGER, OPTIONAL :: D
REAL*8, DIMENSION(:,:) :: X 
REAL*8, INTENT(OUT) :: MU(:)
INTEGER :: NR,NC
NR=SIZE(X,1)
NC=SIZE(X,2)

IF (PRESENT(D)) THEN
    SELECT CASE (D)
      CASE (1)
        MU=SUM(X, DIM=2)/NC
      CASE (2)
        MU=SUM(X, DIM=1)/NR
    END SELECT
ELSE
    MU=SUM(X)/NC/NR     
ENDIF
END SUBROUTINE 


! ====================================================

SUBROUTINE DREP_V(X,N,DP)
! REPEAT a scalar X by N TIMES in order TO GENERATE DP=(X, ...,X)_N
! Aug:  X ---- a value input
!       N ---- repeating times 
!       DP---- OUTPUT vector (X,X,...,X)_N

INTEGER :: N, I
REAL*8 :: X, DP(N)
DO 10 I=1,N
    DP(I)=X
10 CONTINUE
END SUBROUTINE

!============================================

SUBROUTINE DREP_M(X,N,DI,OT)
! REPEAT a vector X by N TIMES TO GENERATE DP=(X, ...,X) by row(DI=1) or col(DI=2)
! Aug:  X ---- a vector input
!       N ---- Repeating times
!       DI---- repeating direction index: DI=1 by row way; DI=2 by col way.
!       OT---- OUTPUT matrix

INTEGER :: DX, N, DI, I, J
REAL*8, INTENT(IN) :: X(:)
REAL*8, INTENT(OUT) :: OT(:,:)
DX=SIZE(X)

SELECT CASE (DI)
CASE (1)
DO 10 I=1,N
  DO 20 J=1,DX
   OT(I,J)=X(J)
 20 CONTINUE
10 CONTINUE

CASE (2)
DO 30 I=1,N
  DO 40 J=1,DX
   OT(J,I)=X(J)
 40 CONTINUE
30 CONTINUE
END SELECT
END SUBROUTINE

!============================================
SUBROUTINE DSD_V(X,DSD)
! Computer the standard deviation of vector 
! Need DVMEAN.F
! Aug:  X ---- Input vector 
!       N ---- length of the vector
!      DSD---- OUTPUT, sd of X

REAL*8, INTENT(IN) :: X(:)
REAL*8, INTENT(OUT) :: DSD
REAL*8 :: DAVG
CALL DMEAN(X,DAVG)
DSD=SQRT(SUM((X-DAVG)**2)/(SIZE(X)-1))
END SUBROUTINE 

!=============================================
SUBROUTINE DSD_M(X,OT,D)
! Computer the standard deviation of matrix with 3 ways (col, row, whole)
! Need DMREP.F DMMEAN.F 
! Aug:  X ---- input matrix
!      OT ---- OUTPUT, a sd vector or a sd value 
!      D  ---- Optional, present(D).eq. .ture. SD for the whole matrix; o.w., D=1, row way, D=2, col way
INTEGER :: NR,NC
INTEGER, OPTIONAL :: D 
REAL*8, INTENT(IN) :: X(:,:) 
REAL*8, INTENT(OUT) :: OT(:) 
REAL*8, ALLOCATABLE :: AVG(:), RP(:,:)
NR=SIZE(X,1)
NC=SIZE(X,2)

ALLOCATE(RP(NR,NC))
IF (PRESENT(D)) THEN

SELECT CASE (D)
CASE (1)
ALLOCATE(AVG(NR))
CALL DMEAN(X,AVG,D)
CALL DREP(AVG,NC,2,RP)
OT=SQRT(SUM((X-RP)**2, DIM=2)/(NC-1))

CASE (2)
ALLOCATE(AVG(NC))
CALL DMEAN(X,AVG,D)
CALL DREP(AVG,NR,1,RP)
OT=SQRT(SUM((X-RP)**2, DIM=1)/(NR-1))

END SELECT

ELSE
ALLOCATE(AVG(1))
CALL DMEAN(X,AVG)
OT=SQRT(SUM((X-AVG(1))**2)/(NC*NR-1))     
ENDIF

IF (ALLOCATED(RP)) DEALLOCATE(RP)
IF (ALLOCATED(AVG)) DEALLOCATE(AVG)

END SUBROUTINE 

!=========================================================================
SUBROUTINE DSTD(X,STX)
! Standarize of the matrix X
! Need DMEAN.F 
!      DSD.F --- DREP.F  
! Aug: X ---- input matrix
!      NR---- row size of X
!      NC---- col size of X
!      STX---- OUTPUT matrix
REAL*8, INTENT(IN)  :: X(:,:)
REAL*8, INTENT(OUT) :: STX(:,:)
INTEGER :: NR, NC 
REAL*8, ALLOCATABLE  :: AVG(:), SD(:), RP1(:,:), RP2(:,:)
NR=SIZE(X,1)
NC=SIZE(X,2)

ALLOCATE(AVG(NC))
ALLOCATE(SD(NC))
ALLOCATE(RP1(NR,NC))
ALLOCATE(RP2(NR,NC))

CALL DMEAN(X,AVG,2)
CALL DSD(X,SD,2)
CALL DREP(AVG,NR,1,RP1)
CALL DREP(SD,NR,1,RP2)

STX=(X-RP1)/RP2 

IF (ALLOCATED(AVG)) DEALLOCATE(AVG)
IF (ALLOCATED(SD)) DEALLOCATE(SD) 
IF (ALLOCATED(RP1)) DEALLOCATE(RP1)
IF (ALLOCATED(RP2)) DEALLOCATE(RP2)

END SUBROUTINE 

!=========================================================================

SUBROUTINE DLAG_M(X,P,XLAG)
! DMLAG routine generates lagged X in order p
! Aug: X ----INPUT matrix  
!      NR----row size of x
!      NC----col size of x
!      P ----lag length order p
!    XLAG----OUTPUT Resulting matrix with (NR, NC*P)

INTEGER :: NR,NC,I,J,P
REAL*8, INTENT(IN) :: X(:,:)
REAL*8, INTENT(OUT) :: XLAG(:,:)
NR=SIZE(X,1)
NC=SIZE(X,2)
DO 70 I=1,P
     XLAG((P+1):NR,(NC*(I-1)+1):(NC*I))=X((P+1-I):(NR-I),:)
     DO 80 J=1,(NC*P)
         XLAG(I,J)=0
     80 CONTINUE
70 CONTINUE

END SUBROUTINE

!=========================================================================

SUBROUTINE DLAG_V(X,P,XLAG)
! DVLAG routine generates lagged X in order p
! Aug: X ----INPUT vector  
!      NR----row size of x
!      P ----lag length order p
!    XLAG----OUTPUT Resulting matrix with (NR, NC*P)

INTEGER :: NR,I,J,P
REAL*8, INTENT(IN) :: X(:)
REAL*8, INTENT(OUT) :: XLAG(:,:)
NR=SIZE(X)
DO 70 I=1,P
     XLAG((P+1):NR,I)=X((P+1-I):(NR-I))
     DO 80 J=1,P
         XLAG(I,J)=0d0
     80 CONTINUE
70 CONTINUE

END SUBROUTINE

!=========================================================================

SUBROUTINE DLAM(LS,K,M,L)
! This routine makes the factor loading matrix L by given its N-K part LS.
!       ( 1            0 ... 0 )
!       (  .diag(K).    0(K*M) ) 
!       (        1     0 ... 0 )
!     L=(  LS                  )
!       ( 0 ... 0      1       )
!       ( 0(M*K)      .diag(M).) 
!       ( 0 ... 0            1 )
! Aug: LS---- in, loading matrix (N-K,KM)
!       K---- IN, factor size
!       M---- IN, y size
!       L---- OUT, the factor loading matrix Lamda
 
INTEGER :: K,M, I, J, P, Q, NR ,NC
 
REAL*8,  INTENT(IN)  :: LS(:,:)
REAL*8,  INTENT(OUT) :: L(:,:)
REAL*8 :: LB(M,(K+M)) 
NR=SIZE(LS,1)
NC=SIZE(LS,2)
L((K+1):(K+NR),:)=LS
DO 100 I =1,K
   DO 90 J =1,NC
       IF (J .EQ. I) THEN 
          L(I,J)=1
       ELSE
          L(I,J)=0 
       ENDIF
  90  CONTINUE  
100 CONTINUE

DO 120 P = 1,M
  DO 130 Q = 1,(K+M)
    IF (Q /= (K+P)) THEN
      LB(P,Q)=0  
    ENDIF
  130  CONTINUE
  LB(P,K+P)= 1
120 CONTINUE
L((K+1+NR):(K+M+NR),:)=LB
END SUBROUTINE

!====================================================================
FUNCTION DDIAG_a(a,NR)
! Creats a diagonal matrix (NR*NR) with the same diagonal element a.
! Aug:  a---- IN, a real*8 value
!      NR---- IN, the leading dimension of the diagonal matrix

INTEGER :: NR,k,l
REAL*8  :: a, DDIAG_a(NR,NR)

DO 140 k = 1,NR
  DO 150 l = 1,NR
     DDIAG_a(k,l)=0d0
  150  CONTINUE
  DDIAG_a(k,k)=a
140 CONTINUE

END FUNCTION

!====================================================================
FUNCTION DDIAG_v(v,NR)

! Creats a diagonal matrix (NR*NR) with the diagonal element equals a vector v.
! Aug:  v---- IN, a real*8 vector
!      NR---- IN, the leading dimension of the diagonal matrix 

INTEGER :: NR,k,l
REAL*8 :: v(NR), DDIAG_v(NR,NR)

DO 160 k = 1,NR
    DO 170 l = 1,NR
       DDIAG_v(k,l)=0d0
    170  CONTINUE
    DDIAG_v(k,k)=v(k)
160 CONTINUE

END FUNCTION

!===================================================================

FUNCTION DZERO(NR,NC)
! Creat a zero matrix with double precison zero
! AUG: NR---- size(zero,1)
!      NC---- size(zero,2)
INTEGER :: NR,NC
REAL*8  :: DZERO(NR,NC)
DZERO(1:NR,1:NC)=0D0
END FUNCTION

!===================================================================
FUNCTION DVEC(X,N)
! Creat a vector with double precison x
! AUG: N ---- the length of vector
INTEGER :: N,I
REAL*8, INTENT(IN) :: X
REAL*8  :: DVEC(N)
DO I=1,N
DVEC(I)=X
END DO
END FUNCTION

!===================================================================
FUNCTION indexnM(K,M,plag)
! Make the elements index for the factors
! i.e. K=2,M=2, PLAG=2 will be
! 1  5
! 2  6
! 0  0
! 0  0
! This returns (/1,2,5,6/)
! Aug: K---- size of factors
!      M---- size of y
!      plag---- the order of lag length

INTEGER, INTENT(IN) :: K, M, PLAG
INTEGER :: INDexNM(K*PLAG)
INTEGER :: I,J
DO 180 I=1,PLAG
    DO 190 J=1,K
       INDexNM(J+K*(I-1))=J+(I-1)*(K+M)
    190 CONTINUE
180 CONTINUE
END FUNCTION

!===================================================================
SUBROUTINE DBIND_MM(X,Y,XY,DIM)
INTEGER :: DIM, NRX, NCX,NRY,NCY
REAL*8, INTENT(IN)  :: X(:,:), Y(:,:)
REAL*8, INTENT(OUT) :: XY(:,:)
! Combine two matrix X,Y in the row way (X) (dim=1), or column way (X Y) (dim=2)
!                                       (Y)
! Aug: X(:,:)  IN
!      Y(:,:)  IN
!      XY(:,:) OUT 
!      DIM =1  (X Y)'
!      DIM =2  (X Y)

NRX=SIZE(X,1)
NCX=SIZE(X,2)
NRY=SIZE(Y,1)
NCY=SIZE(Y,2)

IF (DIM .NE. 1 .AND. DIM .NE. 2) THEN
         STOP "Wrong dimension for the matrix binding operations, DIM should be 1 or 2"
ENDIF

SELECTCASE(DIM)
! combind X and Y in row way
CASE(1)
    IF (NCX .NE. NCY) THEN 
         STOP "X and Y should have same column length"
    ENDIF
    XY(1:NRX,:)=X(:,:)
    XY((NRX+1):(NRY+NRX),:)=Y(:,:)
! combind X and Y in col way
CASE(2)
    IF (NRX .NE. NRY) THEN 
         STOP "X and Y should have same row length"
    ENDIF
    XY(:,1:NCX)=X(:,:)
    XY(:,(NCX+1):(NCY+NCX))=Y(:,:)
END SELECT 

END SUBROUTINE

!====================================================================
!===================================================================
SUBROUTINE DBIND_VV(X,Y,XY,DIM)
! Combine two VECTOR X,Y in the row way (X) (dim=1), or column way (X Y) (dim=2)
!                                       (Y)
! Aug: X(:)  IN
!      Y(:)  IN
!      XY(:,2) or XY(2,:), depend on DIM, OUT 
!      DIM =1  (X Y)'
!      DIM =2  (X Y)

INTEGER :: DIM, NX, NY
REAL*8, INTENT(IN) :: X(:), Y(:)
REAL*8, INTENT(OUT) :: XY(:,:)
NX=SIZE(X)
NY=SIZE(Y)


IF (DIM .NE. 1 .AND. DIM .NE. 2) THEN
        STOP "Wrong dimension for the matrix binding operations, DIM should be 1 or 2"
ENDIF
IF (NX .NE. NY) THEN 
        STOP "X and Y should have same vector length"
ENDIF

SELECTCASE(DIM)
! combind X and Y in row way
CASE(1)

      XY(1,:)=X(:)
      XY(2,:)=Y(:)
! combind X and Y in col way
CASE(2)

      XY(:,1)=X(:)
      XY(:,2)=Y(:)
END SELECT 
END SUBROUTINE
!====================================================================
SUBROUTINE DBIND_VM(X,Y,XY,DIM)
INTEGER :: DIM, NX, NRY,NCY
! Combine a vector X, and a matrix Y in the row way (X) (dim=1), or column way (X Y) (dim=2)
!                                                   (Y)
! Aug: X(:)  IN
!      Y(:,:)  IN
!      XY(:,:) OUT 
!      DIM =1  (X Y)'
!      DIM =2  (X Y)
REAL*8, INTENT(IN) :: X(:), Y(:,:)
REAL*8, INTENT(OUT) :: XY(:,:)
NX=SIZE(X)
NRY=SIZE(Y,1)
NCY=SIZE(Y,2)

IF (DIM .NE. 1 .AND. DIM .NE. 2) THEN
        STOP "Wrong dimension for the matrix binding operations, DIM should be 1 or 2"
ENDIF

SELECTCASE(DIM)
! combind X and Y in row way
CASE(1)
    IF (NX .NE. NCY) THEN 
       STOP "X and Y should have same column length"
    ENDIF
    XY(1,:)=X(:)
    XY(2:(NRY+1),:)=Y(:,:)
! combind X and Y in col way
CASE(2)
    IF (NX .NE. NRY) THEN 
        STOP "X and Y should have same row length"
    ENDIF
    XY(:,1)=X(:)
    XY(:,2:(NCY+1))=Y(:,:)
END SELECT 
END SUBROUTINE

!====================================================================
SUBROUTINE DBIND_MV(X,Y,XY,DIM)
! Combine a matrix X, and a vector Y in the row way (X) (dim=1), or column way (X Y) (dim=2)
!                                                   (Y)
! Aug: X(:,:)  IN
!      Y(:)  IN
!      XY(:,:) OUT 
!      DIM =1  (X Y)'
!      DIM =2  (X Y)

INTEGER :: DIM, NRX,NCX, NY
REAL*8, INTENT(IN) :: X(:,:), Y(:)
REAL*8, INTENT(OUT) :: XY(:,:)
NRX=SIZE(X,1)
NCX=SIZE(X,2)
NY=SIZE(Y)

IF (DIM .NE. 1 .AND. DIM .NE. 2) THEN
        STOP "Wrong dimension for the matrix binding operations, DIM should be 1 or 2"
ENDIF

SELECTCASE(DIM)
! combind X and Y in row way
CASE(1)
       IF (NCX .NE. NY) THEN 
            STOP "X and Y should have same column length"
       ENDIF
       XY(1:NRX,:)=X(:,:)
       XY(NRX+1,:)=Y(:)
! combind X and Y in col way
CASE(2)
       IF (NRX .NE. NY) THEN 
            STOP "X and Y should have same row length"
       ENDIF
       XY(:,1:NCX)=X(:,:)
       XY(:,(NCX+1))=Y(:)
END SELECT 
END SUBROUTINE

!=========================================================================== 
FUNCTION DMATXA(NX,A,LDA,X,INV,POS)
! This function computes A'%*%X%*%A, where A(NX,LDA) and X is a symmetric matrix(nx,nx)
! if INV present, then this computes A'%*%inv(X)%*%A 
! if POS present, then X is a postive symmetric matrix

INTEGER LDA,NX,IRANK
REAL*8, INTENT(IN):: A(:,:),X(:,:)
INTEGER, OPTIONAL, INTENT(IN) :: INV, POS
!REAL*8, INTENT(OUT) :: MATXA(:,:)

REAL*8 DMATXA(LDA,LDA),XA(NX,LDA),XR(NX,NX)
IF (PRESENT(POS)) THEN
    CALL CHFAC(X,IRANK,XR)
    XA=A
    IF (PRESENT(INV)) THEN
      !  CALL DLINRT(NX,XR,NX,2,XR,NX)
        call linrt(xr,xr)
        CALL DTRMM('l','u','t','n',NX,LDA,1D0,XR,NX,XA,NX)
        
    else
        CALL DTRMM('l','u','n','n',NX,LDA,1D0,XR,NX,XA,NX)
    ENDIF
    CALL MXTXF(XA,DMATXA)
ELSE
    XR=X
    IF (PRESENT(INV)) THEN
        CALL LINRG(XR,XR) 
    ENDIF
    CALL DSYMM('L','U',NX,LDA,1D0,XR,NX,A,NX,0D0,XA,NX)
    CALL MXTYF(A,XA,DMATXA)
ENDIF
END FUNCTION

!=========================================================================== 
FUNCTION DMAXAT(LDA,A,NX,X,INV,POS)
! This function computes A%*%X%*%A', where A(LDA,NX) and X is a symmetric matrix(nx,nx)
! if INV present, then this computes A%*%inv(X)%*%A' 
! if POS present, then X is a postive symmetric matrix
INTEGER LDA,NX,IRANK
REAL*8, INTENT(IN):: A(:,:),X(:,:)
INTEGER, OPTIONAL, INTENT(IN) :: INV, POS
!REAL*8, INTENT(OUT) :: MATXA(:,:)

REAL*8 DMAXAT(LDA,LDA),AX(LDA,NX),XR(NX,NX)
IF (PRESENT(POS)) THEN
    CALL CHFAC(X,IRANK,XR)
    AX=A
    IF (PRESENT(INV)) THEN
        CALL DLINRT(NX,XR,NX,2,XR,NX) 
        CALL DTRMM('R','u','n','n',LDA,NX,1D0,XR,NX,AX,LDA)
    else
        CALL DTRMM('R','u','t','n',LDA,NX,1D0,XR,NX,AX,LDA)
    ENDIF
    CALL MXYTF(AX,AX,DMAXAT)
ELSE
    XR=X
    IF (PRESENT(INV)) THEN
       CALL LINRG(XR,XR) 
    ENDIF
    CALL DSYMM('R','U',LDA,NX,1D0,XR,NX,A,LDA,0D0,AX,LDA)
    CALL MXYTF(AX,A,DMAXAT)
ENDIF
END FUNCTION


!================================================================
FUNCTION FINDLARGEST(X)
        ! find the largest value in a vector X
        IMPLICIT NONE
        INTEGER N,I
        REAL*8, INTENT(IN) :: X(:)
        !REAL*8, INTENT(INOUT) :: L !Largest value
        REAL*8 L,FINDLARGEST
        N=SIZE(X)
        !Algorithm is easy
        !Let the largest number be the first one.
        !If you find a number larger than it, store this number and then continue
        L = ABS(X(1))
        DO I = 2, n
           IF (ABS(x(I)) > L) L = ABS(x(I))
        END DO 
        FINDLARGEST=L             
END FUNCTION FINDLARGEST

!=================================================================
FUNCTION uptri(X,NR)
! Make a upper triangular matrix U (nr,nr) from the compressed vector X, 
!which stores in the column way
INTEGER NR,I
REAL*8, INTENT(IN) :: X(:)
REAL*8  uptri(NR,NR)
uptri=DZERO(NR,NR)
DO I=1,NR
     uptri(1:I,I)=X((I*(I-1)/2+1):(I*(I+1)/2))
END DO
END FUNCTION
 
!=================================================================
FUNCTION antitri(U)
! Store a upper triangular matrix U (nr,nr) into a vector X in the column way
INTEGER NR,I
REAL*8, INTENT(IN) :: U(:,:)
REAL*8  antitri(SIZE(U,1)*(SIZE(U,1)+1)/2)
NR=SIZE(U,1)
DO I=1,NR
      antitri((I*(I-1)/2+1):(I*(I+1)/2))=U(1:I,I)
ENDDO
END FUNCTION

!=================================================================
FUNCTION DMAXEIG(X)
INTEGER NR,INFO,I
REAL*8, INTENT(IN) :: X(:,:)
REAL*8 DMAXEIG,WR(SIZE(X,1)),WI(SIZE(X,1)),V(SIZE(X,1),SIZE(X,1)), &
& WORK(3*SIZE(X,1)),EIG(SIZE(X,1)),XX(SIZE(X,1),SIZE(X,1))

NR=SIZE(X,1)

XX=X

CALL DGEEV('N','N',NR,XX,NR,WR,WI,V,NR,V,NR,WORK,3*NR,INFO)

DO I=1,NR
EIG(I)=SQRT((WR(I)**2) + (WI(I)**2))
END DO
DMAXEIG=FINDLARGEST(EIG)
END FUNCTION 

!=================================================================
FUNCTION PHIX_M(PHI,X,SYM)

! PHI(PKM,KM)
! X(PKM,NCX)
! OUTPUT PHIM(PKM,NCX)

INTEGER, OPTIONAL,   INTENT(IN) :: SYM 
REAL*8,              INTENT(IN) :: PHI(:,:), X(:,:) 
INTEGER                  KM,PKM,NCX,PLAG
REAL*8                   PHIX_M(SIZE(X,1),SIZE(X,2))

PKM=SIZE(PHI,1)
KM=SIZE(PHI,2)
PLAG=PKM/KM
NCX=SIZE(X,2)

IF (PLAG==1) THEN 

    IF (PRESENT(SYM)) THEN 
       CALL DSYMM('R','U',KM,NCX,1D0,X,KM,TRANSPOSE(PHI),KM,0D0,PHIX_M,KM)
    ELSE
       CALL MXTYF(PHI,X,PHIX_M)
    ENDIF

ELSE 

     IF (PRESENT(SYM)) THEN 
       CALL DSYMM('R','U',KM,NCX,1D0,X,PKM,TRANSPOSE(PHI),KM,0D0,PHIX_M(1:KM,:),KM)
     ELSE
       CALL MXTYF(PHI,X,PHIX_M(1:KM,:))
     ENDIF

     PHIX_M((KM+1):PKM,:)=X(1:((PLAG-1)*KM),:)
ENDIF

END FUNCTION 

!=================================================================
FUNCTION PHIX_V(PHI,X)
! PHI(PKM,KM)
! X(PKM,NCX)
! OUTPUT PHIM(PKM,NCX)

REAL*8,              INTENT(IN) :: PHI(:,:), X(:) 
INTEGER                  KM,PKM,PLAG
REAL*8                   PHIX_V(SIZE(PHI,1))

PKM=SIZE(PHI,1)
KM=SIZE(PHI,2)
PLAG=PKM/KM

IF (PLAG==1) THEN 

    PHIX_V = MATMUL(X,PHI) 
    
ELSE 
    
    PHIX_V(1:KM) = MATMUL(X,PHI) 

    PHIX_V((KM+1):PKM)=X(1:((PLAG-1)*KM))
     
ENDIF

END FUNCTION 


!=================================================================
FUNCTION PHIVX_M(PHI,X)

! PHI(PKM)
! X(PKM,NCX)
! OUTPUT PHIM(PKM,NCX)

REAL*8,              INTENT(IN) :: PHI(:), X(:,:) 
INTEGER                  PLAG,NCX
REAL*8                   PHIVX_M(SIZE(X,1),SIZE(X,2))

PLAG=SIZE(PHI)
NCX=SIZE(X,2)
  
     PHIVX_M(1,:)=MATMUL(PHI,X)
 
     PHIVX_M(2:PLAG,:)=X(1:(PLAG-1),:)

END FUNCTION 


!=================================================================
FUNCTION PHIVX_V(PHI,X)

! PHI(PKM)
! X(PKM,NCX)
! OUTPUT PHIM(PKM,NCX)

REAL*8,              INTENT(IN) :: PHI(:), X(:) 
INTEGER                  PLAG
REAL*8                   PHIVX_V(SIZE(PHI))

PLAG=SIZE(PHI)
  
     PHIVX_V(1)=DOT_PRODUCT(PHI,X)
 
     PHIVX_V(2:PLAG)=X(1:(PLAG-1))

END FUNCTION 

!==================================================================
FUNCTION MVTV(X,Y)
INTEGER N, I
REAL*8, INTENT(IN) ::  X(:),Y(:)
REAL*8  MVTV(SIZE(X),SIZE(Y))
N=SIZE(X)
DO 710 I=1,N
   MVTV(I,:)=X(I)*Y
710 CONTINUE 
END FUNCTION

!==================================================================
FUNCTION COMBS(n,a)

implicit none 

INTEGER i,COMBS,x(n)
real*8 com
INTEGER, INTENT(IN) :: N,A
do i=1,n
  x(i)=i
end do  
IF (A > N) THEN
    WRITE (*,*) 'The selected element a (from n) should be small than n'
    PAUSE
  else if (A == N) THEN
    COMBS=1
  ELSE if (A < N) then
    if (n > 2*a) then 
       Combs=PRODUCT(x((n-a+1):n))/product(x(1:a))
    else 
       Combs=PRODUCT(x((a+1):n))/product(x(1:(n-a)))
    end if     
END IF 
END FUNCTION


!*****************************************************************

subroutine antivec(nrow,ncol,X_vec,X_mat)

! This routine makes a (nrow x ncol) matrix from a vector:
!											X1	X2
!	X_vec = [X1 X3 X5 X2 X4 X6]'	X_mat =	X3	X4			
!											X5	X6
!
!Written by M.K. Andersson, Sveriges Riksbank, 2004-09-22


implicit none
integer, intent(in) :: nrow,ncol
double precision, dimension(nrow*ncol), intent(in) :: X_vec
double precision, dimension(nrow,ncol), intent(out) :: X_mat

integer :: i

do i = 1,ncol
	X_mat(1:nrow,i)= X_vec(nrow*(i-1)+1:nrow*i)
enddo

end subroutine

!=================================================================

subroutine antivec3(xvec,row,col,rep,x)
integer, intent(in) :: row, col, rep
real*8, intent(in)  :: xvec(:)
real*8, intent(out) :: x(:,:,:)
real*8, allocatable :: xtemp(:,:)
integer i


if (rep==1) then
  call antivec(row,col,xvec,x(:,:,1))
else 
  allocate(xtemp(row*col,rep))
  call antivec(row*col,rep,xvec,xtemp)
  do i=1,rep
     call antivec(row,col,xtemp(:,i),x(:,:,i))  
  end do
  deallocate(xtemp)
end if 

end subroutine


!====================================================================
subroutine dtrmm2(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )
character*1, intent(in)    :: side, uplo, transa, diag
integer,     intent(in)    :: m,n, lda, ldb
real*8,      intent(in)    :: alpha, a(:,:)
real*8,      intent(inout) :: B(:,:)

call dtrmm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB )

end subroutine 



END MODULE