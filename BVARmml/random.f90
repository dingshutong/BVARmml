! ------- Fortran modules: random --------------
! ------- Generate Random distributions --------
! ------- Author: Shutong Ding -----------------
! ------- Last updated: 2016-11-12 -------------
! ------- Dependency: IMSL, matutil, matrix ----

    
module random

include 'link_fnl_static.h'  !DEC$ OBJCOMMENT LIB:'libiomp5md.lib'

use imsl_libraries

use matutil
use matrix

implicit none

! Generate random vector/matrix from Normal distribution N(0,1)
interface rNorm
    module procedure rNorm_V, RNorm_M
end interface

! Generate random vector/matrix from Uniform distribution
interface rUnif
    module procedure rUnif_V, rUnif_M
end interface


    contains
    
!=====================================================================
! Generate a i.i.d vector rUnif_V(nr) from Unif(a, b)    
function rUnif_V(nr, a, b)    

            integer, intent(in)  :: nr   ! The length of generated Uniform vector 
            real*8,  intent(in)  :: a, b ! lower bound and upper bound
            real*8   rUnif_V
            
            call drnun( nr, rUnif_V)
            rUnif_V = a + rUnif_V*(b - a)   
            
end function
    
!=====================================================================
! Generate a i.i.d matrix rUnif_M(nr X nc) from Unif(a, b)    
function rUnif_M(nr, nc, a, b)    

            integer, intent(in)  :: nr, nc ! The number of rows, columns
            real*8,  intent(in)  :: a, b   ! lower bound and upper bound
            real*8   rUnif_M(nr,nc), rVec(nr*nc)
                
            call drnun( nr*nc , rVec)
            rVec = a + rVec*(b - a)
            call antiVec( nr, nc, rVec, rUnif_M)

end function
    
!=====================================================================
! Generate a i.i.d vector rNorm_V(nr) from Normal(0,1)  
function rNorm_V(nr)

            integer, intent(in)  :: nr   ! The length of generated Normal N(0,1) vector 
            integer i
            real*8   rNorm_V(nr)
            
            do i = 1,nr 
                rNorm_V(i) = rnnof()
            end do
            
end function

!=====================================================================
! Generate i.i.d normal matrix rNorm_M(nr, nc) from Normal(0,1)
function rNorm_M(nr, nc)

            integer, intent(in)  :: nr, nc ! The number of rows, columns
            integer i,j
            real*8 rNorm_M(nr, nc)
            
            do  i = 1,nr
                do  j = 1,nc
                    rNorm_M(i,j) = rnnof()
                end do
            end do

end function

!=====================================================================
! Generate i.i.d normal vector rMnorm(size(mean)) from mvn(mean, sqrt(var))

subroutine rMnorm(mean, var, ot, iuvar)
          
real*8,            INTENT(IN)  :: MEAN(:), VAR(:,:)
integer, optional, INTENT(IN)  :: IUVAR
real*8,            INTENT(OUT) :: OT(:)
real*8             R(SIZE(var,1),SIZE(var,1)) ! squared matrix
integer            nr, iu, iRank

! IUVAR is a non-zero index for the variance structure, if IUVAR is presented
! with non-zero integer, then the input of VAR is the Cholesky factors U of 
! the variance.

nr = size(var,1)
iu = 0

if (present(iuvar)) then
      iu = iuvar 
endif     

if (iu .eq. 0) then
      call chfac(var, irank, R)
else
      r = var
endif

call DRNMVN(1, nr, R, nr, ot, 1)
ot = ot + mean

end subroutine


!======================================================================
! Generate a value from Gamma(a, b), where a is shape parameter, b is the scale parameter
function rGamma( shape, scale)

         real*8 shape, scale, rGamma
         call drngam(1, shape, rGamma)
         rGamma = rGamma*scale         
 
end function

!======================================================================
! Generate a value from inv-Gamma2(a, b), where a is shape parameter, b is the scale parameter
function rInvGamma2(shape, scale)

         real*8 shape, scale, rInvGamma2    
         call drngam(1, shape/2d0, rInvGamma2)
         rInvGamma2 = 1d0/rInvGamma2
         
end function

!======================================================================
! Generate Chelosky decomp of covariance matrix for invWishart  

subroutine rInvWishart(var, idf, ot, iuvar)

          integer,           intent(in)  :: idf
          real*8,            intent(in)  :: var(:,:)
          integer, optional, intent(in)  :: iuvar
          real*8,            intent(out) :: ot(:,:)
          integer nr, iu, iRank
          real*8  R(size(var,1),size(var,1))  ! squared matrix
          
          nr = size(var,1)
          
          iu = 0
          
          ! iuvar should not equal to zero
          if (present(iuvar)) then
              iu = iuvar 
          endif    
          
          if (iu .eq. 0) then
              call chfac(var, iRank, R)
          else
              R = var
          endif
          
          call invWISHBART(IDF, OT, S = R)
          
end subroutine

!*******************************************************
subroutine WISHBART( idf, T, S, low, tran )

integer, intent(in) :: idf
double precision, intent(out) :: T(:,:)
double precision, intent(in), optional :: S(:,:)
integer, intent(in), optional :: low, tran

!*****
!     Generates UPPER triangular Bartlett decomposition matrix, T
!     T'T is distributed as W( I, idf ), lower triangle set to zero
!*****
!     IDF  = degrees of freedom of Wishart distribution, idf >= m required
!     T    = Bartlett decomposition matrix, m by m.
!  Optional arguments
!     S    = triangular scale matrix, 
!            T generated as triangular decomposition of W( A, idf ) 
!            with A = S'S or A = S*S' is indicated by tran, m by m
!            UPPER triangular or LOWER triangular as indicated by low
!     low  = Nonzero, Generate T as LOWER triangular matrix, T'T ~ W
!     tran = Nonzero, return UPPER/LOWER triangular matrix but now
!            with T*T' ~ W
!            If a scaling matrix, S, is supplied this must be a S*S' = A
!            decomposition
!******

integer I, J, tdim, lower, transposed, genlow, gentran
double precision DDF, step
double precision, allocatable :: tmpT(:,:)

! default actions
! tran  = 0
! lower = 0 output UT matrix, T'T form
!         1 output LT matrix, T'T form
! tran  = 1
! lower = 0 output UT matrix, T*T' form
!         1 output LT matrix, T*T' form
!
! genlow, internally generate lower or upper triangular matrix
! gentran, internally generate transpose form
!
!             genlow
!             0              1
! gentran 0  UT, T'T       LT, T'T
!         1  UT, T*T'      LT, T*T'

transposed = 0
lower      = 0
genlow     = 0
gentran    = 0


tdim = size(T,1)
if ( tdim /= size(T,2) ) then
  stop 'T must be square matrix for Wishart'
endif

if ( IDF .lt. tdim ) then
  stop 'To few degrees of freedom for Wishart'
endif

! optional arguments

if ( present(low) ) then
  if ( low /= 0 ) then
    lower  = 1
    genlow = 1
  endif
endif
if ( present(tran) ) then
  if ( tran /= 0 ) then
    transposed = 1
    gentran    = 1
  endif
endif

if ( present(S) ) then
  if ( size(S,1) /= tdim .or. size(S,2) /= tdim ) then
    stop 'Wrong dimension for scale matrix in Wishart'
  endif
  if ( tdim > 1 ) then
    ! Simple check for triangularity
    if ( lower == 0 ) then
      ! S should be upper triangular
      if ( S(2,1) /= 0d0 .and. S(1,2) == 0d0 ) then
        stop 'S matrix not upper triangular in Wishart'
      endif
    else
      ! S should be lower triangular
      if ( S(1,2) /= 0d0 .and. S(2,1) == 0d0 ) then
        stop 'S matrix not lower triangular in Wishart'
      endif
    endif
  endif
  ! Generate as transpose when scaling is applied
  gentran = 1 - gentran ! flip
  genlow  = 1 - genlow ! flip
endif

allocate( tmpT(tdim,tdim) )

if ( genlow /= 0 ) then
  ! Generate lower triangular

  if ( gentran /= 0 ) then
    ! T*T' form (transpose of upper triangular T'T form)
    DDF = IDF
    step = -1.0d0
  else
    ! T'T form
    DDF = IDF - tdim + 1
    step = 1.0d0
  endif
  do I = 1, tdim-1
    call DRNCHI (1, DDF, tmpT(I,I))
    tmpT(I,I) = sqrt(tmpT(I,I))
    DDF = DDF + step
    call DRNNOA (tdim-I, tmpT(I+1:,I))
    do J = 1, I
      tmpT(J,I+1) = 0.0D0
    enddo
  enddo
  call DRNCHI(1, DDF, tmpT(tdim,tdim))
  tmpT(tdim,tdim) = sqrt(tmpT(tdim,tdim))

else
  ! Generate upper triangular

  if ( gentran /= 0 ) then
    ! T*T' (transpose of lower triangular T'T form)
    DDF = IDF - tdim + 1
    step = 1.0d0
  else
    ! T'T form
    DDF = IDF
    step = -1.0d0
  endif
  call DRNCHI( 1, DDF, tmpT(1,1) )
  tmpT(1,1) = sqrt(tmpT(1,1))
  do I = 2, tdim
    call DRNNOA( i-1, tmpT(1:i-1,I) )
    DDF = DDF + step
    call DRNCHI( 1, DDF, tmpT(I,I) )
    tmpT(I,I) = sqrt(tmpT(I,I))
    do J = i, tdim
      tmpT(J,i-1) = 0.0D0
    enddo
  enddo

endif

if ( present(S) ) then
  ! apply scaling, this transposes T back to what we want
  if ( lower /= 0 .and. transposed == 0 ) then
    ! lower triangular T'T form requested
    ! S is lower triangular, T is UT with S'T*T'S ~ W( A, idf )
    ! Return B = T'S, LT with B'B ~ W( A, idf )
    call trimult( tmpT, S, T, aup=1, atran=1 )
  elseif ( lower /= 0 .and. transposed /= 0 ) then
    ! lower triangular T*T' form requested
    ! S is lower triangular, T is UT with S*T'T*S' ~ W( A, idf )
    ! Return B = S*T' with B*B' ~ W( A, idf )
    call trimult( S, tmpT, T, bup=1, btran=1 )
  elseif ( lower == 0 .and. transposed == 0 ) then
    ! upper triangular T'T requested
    ! S is upper triangular, T is LT with S'T*T'S ~ W( A, idf )
    ! Return B = T'S, UT with B'B ~ W( A, idf )
    call trimult( tmpT, S, T, bup=1, atran=1 )
  else ! lower == 0, transposed == 1
    ! upper triangular T*T' form requested
    ! S is upper triangular, T is LT with S*T'T*S' ~ W( A, idf )
    ! Return B = S*T' with B*B' ~ W( A, idf )
    call trimult( S, tmpT, T, aup=1, btran=1 )
  endif
else
  T = tmpT
endif

deallocate( tmpT )

end subroutine

!********************************************************

subroutine invWISHBART( idf, T, W, S, invS, low, noWscale )

integer, intent(in) :: idf
double precision, intent(out) :: T(:,:)
double precision, intent(out), optional :: W(:,:)
double precision, intent(in), optional :: S(:,:), invS(:,:)
integer, intent(in), optional :: low, noWscale

!*****
!     Generates UPPER triangular "Bartlett decomposition" matrix, T
!     T'T is distributed as iW( I, IDF ), lower triangle set to zero
!*****
!     IDF  = degrees of freedom of inverse Wishart distribution, idf >= m required
!     T    = Bartlett decomposition matrix, m by m
!  Optional arguments
!     W    = Bartlett decomposition for corresponding Wishart returned in W
!            W*W' is W( I, idf ) or W( A**-1, idf )
!     S    = triangular scale matrix, 
!            T generated as triangular decomposition of iW( A, idf ) 
!            with A = S'S, m by m
!            UPPER/LOWWER triangular as indicated by low
!     invS = triangular scale matrix based on inverse of A
!            T generated as triangular decomposition of iW( A, idf ) 
!            with A^-1 = invS*invS' or A = invS^-T*invS^-1, m by m
!     low  = Nonzero, Generate T as LOWER triangular matrix, T'T ~ iW
!     noWscale = Nonzero, don't scale W (if present) with S (if present)
!            Scaling with S requires a matrix inversion and is better done
!            outside invWishBart if S is constant
!            W is always scaled with invS if this is supplied
!******

integer tdim, genlow, gentran, lower, noScaleW
double precision, allocatable :: tmpT(:,:), tmpW(:,:), Si(:,:)

tdim = size(T,1)
lower   = 0
genlow  = 0
gentran = 1
noScaleW = 0

if ( present(S) .and. present(invS) ) then
  stop 'Only one of S and invS arguments allowed in InvWishBart'
endif
if ( present(W) ) then
  if ( size(W,1) /= tdim .or. size(W,2) /= tdim ) then
    stop 'Dimension for W wrong in InvWishBart'
  endif
endif

if ( present(low) ) then
  if ( low /= 0 ) then
    genlow = 1
    lower  = 1
  endif
endif
if ( present(noWscale) ) then
  if ( noWscale /= 0 ) noScaleW = 1
endif

if ( present(S) ) then
  if ( size(S,1) /= tdim .or. size(S,2) /= tdim ) then
    stop 'Wrong dimension for S scale matrix in invWishart'
  endif
  if ( tdim > 1 ) then
    ! Simple check for triangularity
    if ( lower == 0 ) then
      ! S should be upper triangular
      if ( S(2,1) /= 0d0 .and. S(1,2) == 0d0 ) then
        stop 'S matrix not upper triangular in invWishart'
      endif
    else
      ! S should be lower triangular
      if ( S(1,2) /= 0d0 .and. S(2,1) == 0d0 ) then
        stop 'S matrix not lower triangular in invWishart'
      endif
    endif
  endif
endif

allocate( tmpT(tdim,tdim) )
if ( present(invS) ) then
  ! Generate Wishart and apply scaling if inverse scale matrix
  ! gentran = 1 => T*T' ~ W( invS*invS, idf )
  call wishbart( idf, tmpT, S=invS, low=genlow, tran=gentran )
else
  ! apply scaling after inverting
  ! gentran = 1 => T*T' ~ W( I, idf )
  call wishbart( idf, tmpT, low=genlow, tran=gentran )
endif  

if ( present(W) ) then
  ! Get Wishart factor
  allocate( tmpW(tdim,tdim) )
  tmpW = tmpT
endif

! invert to get decomposition for inverse Wishart
! genlow = 0, upper triangular; 1, lower triangular
! gentran = 1, T'T form; 0, T*T' form after inverting
call DLINRT( tdim, tmpT, tdim, 2-genlow, tmpT, tdim )

if ( present(S) ) then
  ! apply scaling
  if ( lower == 0 ) then
    ! want upper triangular, 
    ! S is UT and T UT with S'T'T*S ~ iW( A, idf )
    ! return B = T*S with B'B ~ iW( A, idf )
    call trimult( tmpT, S, T, aup=1, bup=1 )
  else
    ! want lower triangular,
    ! S is LT and T LT with S'T'T*S ~ iW( A, idf )
    ! return B = T*S with B'W ~ iW( A, idf )
    call trimult( tmpT, S, T )
  endif
  if ( present(W) ) then
    if ( noScaleW == 0 ) then
      ! Scale with inv(S), so that W'W or W*W' is  W( S^-1*S^-T, df )
      allocate( Si(tdim,tdim) )
      call DLINRT( tdim, S, tdim, 2-lower, Si, tdim )
      ! W is on W*W' form if gentran=1
      ! S'S = A, A^-1 = Si*Si', Si*W*W'Si' ~ W( A^-1, idf )
      ! return B = Si*W with B*B' ~ W( A^-1, idf )
      if ( lower == 0 ) then
        ! upper triangular
        call trimult( Si, tmpW, W, aup=1, bup=1 )
      else
        call trimult( Si, tmpW, W )
      endif
      deallocate( Si )
    else
      W = tmpW
    endif
  endif
else
  T = tmpT
  if ( present(W) ) W = tmpW
endif

deallocate( tmpT )
if ( allocated(tmpW) ) deallocate(tmpW)

end subroutine

!*******************************************************

function ln_MVN_dens( x, cov, inu)

! return log of multivariate normal density excluding normalizing constant
! to avoid underflow
! This routine should only be used for computations with the same dimension of
! x (number of evaluation horizons) 

double precision, intent(in) :: x(:), cov(:,:)
integer, intent(in), optional :: inu
double precision ln_MVN_Dens

double precision croot(size(cov,1),size(cov,2)), xtmp(size(x)), det1, det2, tmp
double precision, parameter :: ln10 = 2.302585092994045684d0, &
							   ln2pi = 1.8378770664093454836d0
integer m,indu

indu=0

if (present(inu)) then
   indu=inu
endif 


m = size(x)
if ( m /= size(cov,1) ) then
	print *, "Dimensions do not match in ln_MVN_dens"
	stop
endif

tmp = - m*ln2pi
!tmp = 0.0d0
if (indu .eq. 0) then 
   CALL DLFTDS( m, cov, m, croot, m )
else 
   croot=cov
end if
  
CALL DLFDDS ( m, croot, m, DET1, DET2 )

tmp = tmp - log( det1 ) - det2*ln10

CALL DLFSDS( m, croot, m, x, xtmp )

tmp = tmp - dot_product( x, xtmp )

ln_MVN_dens = tmp/2.0d0

end function

!===================================================================

FUNCTION ln_IW_dens(ppsi,psis,pv,INU1,INU2)

real*8,            INTENT(IN) :: ppsi(:,:), pv,psis(:,:)
integer, optional, INTENT(IN) :: INU1,INU2

integer q, I ,IRANK,IU1,IU2

real*8 ln_IW_dens
! Whishart part

double precision, parameter :: ln10 = 2.302585092994045684d0, &
&							   ln2pi = 1.8378770664093454836d0, &
&                              ln2 =   0.6931471805599453094D0,&
&                              lnpi=   1.1447298858494001742D0

double precision UCR1(size(ppsi,1),size(ppsi,2)), det1, det2, tmp, &
& UCR2(size(ppsi,1),size(ppsi,2)),UCR12(size(ppsi,1),size(ppsi,2)), gam(size(ppsi,1)), &
& UCR22(size(ppsi,1),size(ppsi,2))

q=SIZE(ppsi,1)
! Derterminat of ppsi

tmp=0D0
IU1=0
IU2=0

IF (PRESENT(INU1)) THEN 
    IU1=INU1
END IF
 
IF (PRESENT(INU2)) THEN 
    IU2=INU2
END IF

IF (IU1 .EQ. 0) THEN 
    CALL CHFAC(ppsi,IRANK,UCR1)
ELSE 
    UCR1=ppsi
ENDIF 

IF (IU2 .EQ. 0) THEN 
    CALL CHFAC(psis,IRANK,UCR2)
ELSE 
    UCR2=psis  
ENDIF 


CALL LFDDS(UCR1,det1,det2)

tmp = tmp+dble(pv)*(log( det1 ) + det2*ln10)


CALL  LFDDS(UCR2,det1,det2)
CALL  DLINRT(q,UCR2,q,2,UCR2,q)
! Derterminat of psis

tmp = tmp-dble(pv+q+1)* (log( det1 ) + det2*ln10)

! trace of Psi^-1*Ppsi
! R^(-1)*R^(-1)'*U'*U
! R^(-1)*(U*R^(-1))'*U
! U*R^(-1)

CALL trimult(UCR1,UCR2,UCR12,1,1,0,0)
! R^(-1)=R^(-1)*(U*R^(-1))'
CALL trimult(UCR12,UCR1,UCR22,1,1,1,0)
! R^(-1)=R^(-1)*U
CALL dtrmm2('L','U','N','N',q,q,1d0,UCR2,q,UCR22,q)

tmp= tmp-sum(DIAGONALS(UCR22))

! Normalizing Constant part
tmp= tmp-pv*q*ln2-dble(q*(q-1))*lnpi/2D0

! Gamma part
DO I=1,q
gam(I)=ALNGAM(dble(pv+1-I)/2D0)
end do
tmp=tmp-sum(gam)*2D0

ln_IW_dens=tmp/2D0
END FUNCTION


end module