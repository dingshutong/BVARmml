module priorspost


INCLUDE 'link_fnl_static.h'!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries
use structures
use matrix
use random

implicit none 
 
contains


subroutine initNIW( tData, tSampler )

type(VARData), intent(in) :: tData
type(NDsampler), intent(out) :: tSampler

integer npar, m

call clearNIW( tSampler )

m = tData%m
npar = tData%m * tData%k

allocate( tSampler%invPsi_Sampled(m,m,1), &
          tSampler%varGamma_1(npar,npar), &
          tSampler%SGSG(npar),            &
          tSampler%Gamma_1(npar),         &
          tSampler%Gamma_gen(npar),       &
          tSampler%err(tData%Nobs,m),        &
          tSampler%ESS(m,m),              &
          tSampler%Psi_root(m,m),         &
          tSampler%S_chol(m,m) )

! start sampler at inverse of OLS error variance
tSampler%invPsi_Sampled(:,:,1) = tData%ESS/dble(tData%Nobs)
CALL DLINDS( m, tSampler%invPsi_Sampled(:,:,1), m, tSampler%invPsi_Sampled(:,:,1), m )

end subroutine

!-------------------------------------------------------

subroutine clearNIW( tSampler )

type(NDsampler), intent(inout) :: tSampler

if ( allocated( tSampler%invPsi_Sampled ) )  deallocate ( tSampler%invPsi_Sampled )
if ( allocated( tSampler%varGamma_1 ) )   deallocate ( tSampler%varGamma_1  )
if ( allocated( tSampler%SGSG ) )         deallocate ( tSampler%SGSG )
if ( allocated( tSampler%Gamma_1 ) )      deallocate ( tSampler%Gamma_1 )
if ( allocated( tSampler%Gamma_gen ) )    deallocate ( tSampler%Gamma_gen )
if ( allocated( tSampler%err ) )          deallocate ( tSampler%err )
if ( allocated( tSampler%ESS ) )          deallocate ( tSampler%ESS )
if ( allocated( tSampler%Psi_root ) )     deallocate ( tSampler%Psi_root )
if ( allocated( tSampler%S_chol ) )        deallocate ( tSampler%S_chol )

end subroutine

!-------------------------------------------------------

subroutine initNIWpost(tdata, tPost )

type(VARData), intent(in) :: tData
type(NIWposterior), intent(out) :: tPost

integer npar, m

call clearNIWpost(tPost)

m = tData%m
npar = m * tData%k

allocate(tPost%varGamma(npar,npar),tPost%vecGamma(npar),tPost%S_chol(m,m))

end subroutine

!-------------------------------------------------------

subroutine clearNIWpost( tPost )

type(NIWposterior), intent(inout) :: tPost

if ( allocated( tPost%S_chol ) )           deallocate (tPost%S_chol)
if ( allocated( tPost%varGamma ) )        deallocate ( tPost%varGamma)
if ( allocated( tPost%vecGamma ) )        deallocate ( tPost%vecGamma)

end subroutine



!-------------------------------------------------------

subroutine sampleNIW( tData, tPrior, tSampler, Gamma_draw, Psi_draw, &
                                Psiroot_draw, Gamma_anti,Gamma_vec, tPost)

type(VARData), intent(in) :: tData
type(priorND), intent(in), target :: tPrior
type(NDsampler), intent(inout), target :: tSampler
double precision, intent(out) :: Gamma_draw(:,:,:)
double precision, intent(out), optional :: Psi_draw(:,:,:), Psiroot_draw(:,:,:), &
                                           Gamma_anti(:,:,:),Gamma_vec(:)
type(NIWposterior), intent(out),optional :: tPost                                           

!------

integer T, m, k, i, j, ii, jj
logical diagvarGamma

type(priorND), pointer :: prior
type(NDsampler), pointer :: sampler

prior => tPrior
sampler => tSampler

T = tData%Nobs
m = tData%m
k = tData%k

if ( allocated( prior%invdiagvarGamma ) ) then
	diagvarGamma = .true.
else
	diagvarGamma = .false.
endif
!***********************
!Setup for MN-draws of B
!***********************

! PSIi_samp = previous PSI-draw, cond data variance
sampler%varGamma_1 = KRON( sampler%invPsi_Sampled(:,:,1), tdata%ZPZ )
! weighted av of prior mean and data
sampler%SGSG   =  prior%invvarGammaGamma &
					 + matmul( sampler%varGamma_1, tData%vecGamma_OLS )	  
! inverse of conditional posterior variance
if ( diagvarGamma ) then
	do i = 1, m*k
		sampler%varGamma_1(i,i) = prior%invdiagvarGamma(i) &
									 + sampler%varGamma_1(i,i)
	end do
else
	sampler%varGamma_1 = prior%invvarGamma + sampler%varGamma_1	
endif
! Mean of conditional posterior
! factor inverse posterior variance
call DLFTDS( k*m, sampler%varGamma_1, k*m, sampler%varGamma_1, k*m )	
! solve for posterior mean
call DLFSDS( k*m, sampler%varGamma_1, k*m, sampler%SGSG, sampler%Gamma_1 ) 
! Variance of posterior, invert to get factor of variance matrix
call DLINRT( k*m, sampler%varGamma_1, k*m, 2, sampler%varGamma_1, k*m )   


! Draw from MVN distribution
call DRNNOA( k*m, sampler%Gamma_gen )
call DTRMV( 'U', 'N', 'N', k*m, sampler%varGamma_1, k*m, sampler%Gamma_gen, 1 ) 
ii = 1; jj = k
do j = 1, m
	Gamma_draw(:,j,1) = sampler%Gamma_1(ii:jj) + sampler%Gamma_gen(ii:jj)
	ii = ii+k; jj = jj+k
enddo

if (present(Gamma_vec)) then
   Gamma_vec=sampler%Gamma_1 + sampler%Gamma_gen
end if 


if ( present(Gamma_anti) ) then
	ii = 1; jj = k
	do j = 1, m
		Gamma_anti(:,j,1) = sampler%Gamma_1(ii:jj) - sampler%Gamma_gen(ii:jj)
		ii = ii+k; jj = jj+k
	enddo
endif

!*******************************
!Set up for draws of the Wishart
!*******************************
!B_samp = previous B-draw
CALL DMRRRR( T, k, tData%Z, T, k, m, Gamma_draw, k, T, m, sampler%err, T )	
sampler%err = tData%Y - sampler%err			!Y - Z*Gamma
CALL DMXTXF( T, m, sampler%err, T, m, sampler%ESS, m )	!(Y-ZB)'(Y-ZB)
! Psi ~ iW( YZBYZB'YZBYZB, T )
! Psi^-1 ~ W( (YZBYZB'YZBYZB)^-1, T )
if (prior%diagpsi) then
   do j=1,m
       sampler%ESS(j,j)=sampler%ESS(j,j)+prior%priorVecPsi(j)
   end do
else 
   sampler%ESS=sampler%ESS+prior%priorPsi
end if    

! factor, upper triangular
call DLFTDS( m, sampler%ESS, m, sampler%ESS, m ) 
! Psi ~ iW( YESS'ESS, T ), Psi_root'Psi_root ~ iW, PSIi_samp*PSIi_samp' ~ W
call invWishBart( prior%priordf+T, sampler%Psi_root, S=sampler%ESS, W=sampler%S_chol ) 
! PSIi_samp ~ W
call trimult( sampler%S_chol, sampler%S_chol, sampler%invPsi_Sampled(:,:,1), aup=1, bup=1, btran=1 )

if ( present( Psi_draw ) ) call trimult( sampler%Psi_root, sampler%Psi_root, Psi_draw(:,:,1), aup=1, bup=1, atran=1 )
if ( present( Psiroot_draw ) ) Psiroot_draw(:,:,1) = sampler%Psi_root

if (present(tPost)) then 

   call trimult(sampler%varGamma_1,sampler%varGamma_1,tPost%varGamma,1,1,0,1)
   tPost%vecGamma=sampler%Gamma_1
   tPost%S_chol=sampler%ESS
   
end if 


! OR this, equivalent
!  ! factor, upper triangular
!  call DLFTDS( m, YZBYZB, m, YZBYZB, m ) 
!  ! invert to get factor of scale matrix for Wishart
!  ! Psi^-1 ~ W( YZBYZB*YZBYZB', T )
!  call DLINRT( m, YZBYZB, m, 2, YZBYZB, m)     
!  ! PSIi_samp*PSIi_samp' ~ W, upper triangular
!  call WishBart( T, S_chol, S=YZBYZB, tran=1 )
!  ! Psi_root'Psi_root ~ iW
!  call DLINRT( m, S_chol, m, 2, Psi_root, m )
!  ! PSIi_samp ~ W
!  call trimult( S_chol, S_chol, PSIi_samp, aup=1, bup=1, atran=1 )

end subroutine




!=================================================Normal inverse Whishart prior ================================================!
subroutine make_NWprior(k, m, b0, Psi0, Omega0, priordf, tprior, priorB, priorPsi, priorOmega)
integer, intent(in)          :: k,m
real*8,  intent(in)          :: b0,Psi0(:),Omega0(:)
integer, intent(in)          :: priordf
type(priorNW), intent(out)   :: tprior
real*8, optional, intent(in) :: priorB(:,:),priorPsi(:,:),priorOmega(:,:)
real*8, allocatable          :: priorgamma(:,:)

call clear_NWprior(tprior)

allocate(priorgamma(k,m))

call DREP(DVEC(b0,m),k,1,priorgamma)

allocate(tprior%priorGamma(k,m),tprior%priorVecPsi(m),tprior%priorvecOmega(k))

tprior%priorGamma=priorgamma
tprior%priorVecPsi=Psi0
tprior%priorVecOmega=Omega0
tprior%priordf=priordf
tprior%diagPsi=.TRUE.
tprior%diagOmega=.TRUE.


IF (present(priorB)) then
      tprior%priorGamma=priorB
end if
IF (Present(priorPsi)) then
    allocate(tprior%priorPsi(m,m))
    tprior%priorPsi=priorPsi
    tprior%diagPsi=.false.
end if 

IF (Present(priorOmega)) then
    allocate(tprior%priorOmega(k,k))
    tprior%priorOmega=priorOmega
    tprior%diagOmega=.false.
end if 

deallocate(priorgamma)
end subroutine

!======================================================================================
subroutine clear_NWprior(tprior)
type(priorNW), intent(inout)  :: tprior
If ( allocated(tprior%priorGamma) )  deallocate(tprior%priorGamma)
If ( allocated(tprior%priorVecPsi) )  deallocate(tprior%priorVecPsi)
If ( allocated(tprior%priorVecOmega) )  deallocate(tprior%priorVecOmega)
If ( allocated(tprior%priorOmega) )  deallocate(tprior%priorOmega)
If ( allocated(tprior%priorPsi) )  deallocate(tprior%priorPsi)

end subroutine


!=================================================Normal inverse Whishart posterior ================================================!
subroutine make_NWsamples(tprior,tdata,mc_reps,tsample)
type(priorNW), intent(in) :: tprior
type(VARdata), intent(in) :: tdata
integer,       intent(in) :: mc_reps
type(NWsampler), intent(out):: tsample

integer k,m,p,T,postdf,i,j,Ss,irank
real*8, allocatable :: B_ols(:,:), ZPZ(:,:), Psi_root(:,:),priorpsi(:,:),prioromega(:,:), & 
      Otemp(:,:),b_0(:,:),b_1(:,:),psi_post(:,:),B_draw(:,:),BOB(:,:),Omega_root(:,:),B_vec(:)
                       
!********************************
!** SAMPLER **
!********************************


call clear_NWsample(tsample)


!============initialized================
T=tdata%Nobs
k=tdata%k
m=tdata%m
p=tdata%p
sS=mc_reps

allocate(tsample%Gamma(k,m,sS),tsample%Psi(m,m,sS),tsample%Psipost(m,m),tsample%Omegapost(k,k), &
        tsample%Gammapost(k,m),tsample%Omegapostchol(k*(k+1)/2),tsample%psipostchol(m*(m+1)/2))
        
allocate(B_ols(k,m),ZPZ(k,k),Psi_root(m,m),Otemp(k,k),b_0(k,m),b_1(k,m),B_draw(k,m),psi_post(m,m), &
        bob(m,m),omega_root(k,k),B_vec(k*m))

B_ols=tdata%Gamma_OLS
ZPZ=tdata%ZPZ
B_0=tprior%priorGamma

! posterior Psi
if ( tprior%diagPsi ) then
    allocate(priorpsi(m,1))
    priorpsi(:,1)=tprior%priorvecpsi
	!Psi_root = matmul( matmul( transpose(B_ols), ZPZ ), B_ols )
	Psi_root =DMATXA(k,B_ols,m,ZPZ)
	do i = 1, m
		Psi_root(i,i) = Psi_root(i,i) + priorPsi(i,1)
	end do
else
    allocate(priorpsi(m,m))
    priorpsi=tprior%priorpsi
	!Psi_root = priorPsi + matmul( matmul( transpose(B_ols), ZPZ ), B_ols )
    Psi_root =priorpsi + DMATXA(k,B_ols,m,ZPZ)
endif

Psi_root = Psi_root + tdata%ESS


! Posterior Omega inverse
if ( tprior%diagOmega ) then
    allocate(priorOmega(k,1))
    priorOmega(:,1)=tprior%priorvecOmega
	Otemp = ZPZ
	do i = 1, k
		Otemp(i,i) = Otemp(i,i) + 1.0d0/priorOmega(i,1)
		! B_1=(priorOmega^-1)* B_0
		B_1(i,:) = B_0(i,:)/priorOmega(i,1)
	end do
	
else
    allocate(prioromega(k,k))
    prioromega=tprior%prioromega
	!CALL DLFTDS( k, priorOmega, k, Otemp, k )
	CAll chfac(prioromega,irank,Otemp)
	call DLINRT( k, Otemp, k, 2, Otemp, k)
	! Otemp is now upper triangular root 
    ! priorOmegaInv = Otemp*Otemp'
    !Otemp = matmul(Otemp,transpose(Otemp))
	call trimult(Otemp,Otemp,Otemp,1,1,0,1)
	!B_1 = matmul( Otemp, B_0 )
	call DSYMM('l','u',k,m,1d0,Otemp,k,B_0,k,0d0,B_1,k)
	Otemp = Otemp + ZPZ
endif

! until here, Otemp is the inverse of postOmega 

call MXTYF(B_0,B_1,BOB)	
Psi_root = Psi_root + BOB


! solve for posterior mean, B_1

B_draw = B_1 + tdata%ZPY

call LINDS(Otemp,Otemp)

call chfac(Otemp,irank,Omega_root)

tsample%Omegapost=Otemp
tsample%Omegapostchol = antitri(Omega_root)

Psi_root = Psi_root + DMATXA(k,B_draw,m,Otemp)

!B_1 = matmul(Otemp,B_draw)
call DSYMM('l','u',k,m,1d0,Otemp,k,B_draw,k,0d0,B_1,k)

tsample%Gammapost=B_1

call chfac( Psi_root, irank, Psi_post )

tsample%psipost=Psi_root
tsample%psipostchol=antitri(Psi_post)

postdf = T + tprior%priordf



do i = 1,SS

	!if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
	!		write(*,*) 'Sample ',i      
	!endif

	! draw inverse Wishart
	! Psi_draw'Psi_draw ~ iW( Psi_root'Psi_root, postdf )
	
	call RinvWishart(Psi_post,postdf,Psi_root,1)
	
	! Draw conditional normal
	! ( Psi_root' kron Otemp )*vec( B_draw ) = vec( Otemp' * B_draw * Psi_root' )
	
	call RMNORM(vector(B_1),KRON(Psi_root,Omega_root),B_vec,1)
	
    CALL antivec(k,m,B_vec,B_draw) 

    tsample%psi(:,:,i)  = psi_root
    tsample%Gamma(:,:,i)= B_draw

end do

deallocate(B_ols,ZPZ,Psi_root,priorpsi,prioromega,Otemp,b_0,b_1,psi_post,BOB,omega_root,B_draw,B_vec)
                      
end subroutine


!======================================================================================


subroutine clear_NWsample(tsample)
type(NWsampler), intent(inout)  :: tsample
If ( allocated(tsample%Gamma) )  deallocate(tsample%Gamma)
If ( allocated(tsample%Psi) )  deallocate(tsample%Psi)
If ( allocated(tsample%Gammapost) )  deallocate(tsample%Gammapost)
If ( allocated(tsample%Omegapost) )  deallocate(tsample%Omegapost)
If ( allocated(tsample%Psipost) )  deallocate(tsample%Psipost)
If ( allocated(tsample%Omegapost) )  deallocate(tsample%Omegapost)
If ( allocated(tsample%Omegapostchol) )  deallocate(tsample%Omegapostchol)
If ( allocated(tsample%psipostchol) )  deallocate(tsample%psipostchol)

end subroutine


!=======================================================================================

subroutine make_NWpost(tprior,tdata,tpost)
! This routine computes all the posterior for the normal invwhishart posteriors
! Input
! tprior --------------NW prior
! tdata  --------------VAR data

! Output
! tpost  --------------tPost 

type(priorNW), intent(in)          :: tprior
type(VARdata), intent(in)          :: tdata
type(NWposterior),   intent(out)   :: tPost


integer k,m,p,T,i,irank
real*8, allocatable :: B_ols(:,:), ZPZ(:,:), Psi_root(:,:),priorpsi(:,:),prioromega(:,:), & 
      Otemp(:,:),b_0(:,:),b_1(:,:),psi_post(:,:),B_draw(:,:),BOB(:,:)




!============initialized================
T=tdata%Nobs
k=tdata%k
m=tdata%m
p=tdata%p

! clean allocation
if(allocated(tpost%Omega_invchol)) deallocate(tpost%Omega_invchol)
if(allocated(tpost%Omega)) deallocate(tpost%Omega)
if(allocated(tpost%S_chol)) deallocate(tpost%S_chol)
if(allocated(tpost%vecGamma)) deallocate(tpost%vecGamma)


! New allocation        
allocate(B_ols(k,m),ZPZ(k,k),Psi_root(m,m),Otemp(k,k),b_0(k,m),b_1(k,m),B_draw(k,m),psi_post(m,m), &
        bob(m,m))

allocate(tpost%Omega_invchol(k,k),tpost%Omega(k,k),tpost%vecGamma(k*m),tpost%S_chol(m,m))

tpost%postdf=T+tprior%priordf

B_ols=tdata%Gamma_OLS
ZPZ=tdata%ZPZ
B_0=tprior%priorGamma
Psi_root=0d0
Psi_root = matmul( matmul( transpose(B_ols), ZPZ ), B_ols )
!Psi_root =DMATXA(k,B_ols,m,ZPZ)

! posterior Psi
if ( tprior%diagPsi ) then
    allocate(priorpsi(m,1))
    priorpsi(:,1)=tprior%priorvecpsi
	do i = 1, m
		Psi_root(i,i) = Psi_root(i,i) + priorPsi(i,1)
	end do
else
    allocate(priorpsi(m,m))
    priorpsi=tprior%priorpsi
    Psi_root = Psi_root + priorpsi 
endif

Psi_root = Psi_root + tdata%ESS

! Posterior Omega inverse
if ( tprior%diagOmega ) then
    allocate(priorOmega(k,1))
    priorOmega(:,1)=tprior%priorvecOmega
	Otemp = ZPZ
	do i = 1, k
		Otemp(i,i) = Otemp(i,i) + 1.0d0/priorOmega(i,1)
		! B_1=(priorOmega^-1)* B_0
		B_1(i,:) = B_0(i,:)/priorOmega(i,1)
	end do
	
else
    allocate(prioromega(k,k))
    prioromega=tprior%prioromega
	CALL DLFTDS( k, priorOmega, k, Otemp, k )
	!CAll chfac(prioromega,irank,Otemp)
	call DLINRT( k, Otemp, k, 2, Otemp, k)
	! Otemp is now upper triangular root 
    ! priorOmegaInv = Otemp*Otemp'
    !Otemp = matmul(Otemp,transpose(Otemp))
	call trimult(Otemp,Otemp,Otemp,1,1,0,1)
	!B_1 = matmul( Otemp, B_0 )
	call DSYMM('l','u',k,m,1d0,Otemp,k,B_0,k,0d0,B_1,k)
	Otemp = Otemp + ZPZ
endif

! until here, Otemp is the inverse of postOmega 


!call MXTYF(B_0,B_1,BOB)	
!Psi_root = Psi_root + BOB
Psi_root = Psi_root + matmul( transpose(B_0), B_1 )

! solve for posterior mean, B_1

B_draw = B_1 + tdata%ZPY

call LINDS(Otemp,Otemp)
!B_1 = matmul(Otemp,B_draw)
call DSYMM('l','u',k,m,1d0,Otemp,k,B_draw,k,0d0,B_1,k)

!call DLFTDS( k, Otemp, k, Otemp, k )

!do i = 1, m
!	CALL DLFSDS( k, Otemp, k, B_draw(:,i), B_1(:,i) )
!end do

Psi_root = Psi_root - matmul( transpose(B_1), B_draw )
call DLFTDS( m, Psi_root, m, Psi_post, m )


tpost%vecgamma=vector(B_1)
tpost%S_chol=Psi_post

!call DLINRT( k, Otemp, k, 2, Otemp, k )

!tpost%Omega_invchol=Otemp

!call trimult(Otemp,Otemp,Otemp,1,1,0,1)

tpost%Omega=Otemp


! clear 

deallocate(B_ols, ZPZ, Psi_root,priorpsi,prioromega, & 
      Otemp,b_0,b_1,psi_post,B_draw,BOB)


end subroutine


!==================================================================================================

subroutine make_NWsample1(tpost,gamma_draw,psichol_draw)

! This routine uses all the posteriors to draw one NW sample 
! Input
! gammapost ------------ post mean B
! psipostchol ------------ The upper-triangular Cholesky factor of post_psi
! omegapostchol ---------- The upper-triangular Cholesky factor of post_omega
! dfpost        ---------- The posterior degrees of freedom

!Output
!gamma_draw -------------one gamma draw
!psichol_draw   -------------one upper-triangular Cholesky factor of psi_draw

type(NWposterior),  intent(in)  :: tpost
real*8,             intent(out) :: gamma_draw(:,:), psichol_draw(:,:)


integer km,m,k,j
real*8, allocatable :: B_vec(:),psi_root(:,:),omegachol(:,:)
!real*8, allocatable :: psi_root(:,:),B_draw(:,:),b1(:,:)

km = size(tpost%vecgamma)
m=size(tpost%S_chol,1) 
k=km/m

allocate(B_vec(km),psi_root(m,m),omegachol(k,k))

!allocate(B_draw(k,m),psi_root(m,m),b1(k,m))
!call RinvWishart(psipostchol,dfpost,psi_root,1)

call invWishBart( tpost%postdf, Psi_root, S=tpost%s_chol) 

psichol_draw=psi_root

call dlftds(k,tPost%Omega,k,omegachol,k)	
	
call RMNORM(tpost%vecgamma,KRON(psi_root,omegachol),B_vec,1)
	
CALL antivec(k,m,B_vec,gamma_draw) 

	
	! Draw conditional normal
	! ( Psi_root' kron Otemp )*vec( B_draw ) = vec( Otemp' * B_draw * Psi_root' )
	
	!do j = 1, m
	!	call DRNNOA( k, B_draw(:,j) )
	!end do
	
	!call antivec(k,m,tpost%vecgamma,b1)
	!CALL DTRMM( 'L', 'U', 'T', 'N', k, m, 1.0d0, tpost%omega_invchol, k, B_draw, k)
	!CALL DTRMM( 'R', 'U', 'T', 'N', k, m, 1.0d0, psi_root, k, B_draw, k)
	!gamma_draw = b1 + B_draw


deallocate(B_vec,psi_root,omegachol)
!deallocate(psi_root,B_DRAW,b1)

end subroutine


!---------------------------------------------------------------------------

subroutine make_Littermanprior( tData, priorspec, tPrior )

type(VARData),   intent(in) :: tData
type(PriorInfo), target, intent(in) :: priorspec
! priorspec%Litterman%tightness(1) = overall tightness
! priorspec%Litterman%tightness(2) = lag decay
! priorspec%Litterman%tightness(3) = tightness on deterministic
! priorspec%Litterman%tightness(4) = relative tightness on foreign lags
! priorspec%Litterman%tightness(5) = tightness on endogenous in exogenous
!character*12, intent(in) :: LEVEL(:), ENDOEXO(:)
!character*(*), intent(in) :: priortype
!double precision, intent(in) :: hyper(:)
type(VARPrior), intent(out) :: tPrior

!--------

integer T, m, p, q, k, r, i, j, jj, lag, num, num_ex, num_end
integer, allocatable :: exovar(:), endovar(:), ind(:)
double precision, allocatable :: sig(:), Gamma(:,:), Svec(:)
!character*12 LEV, ENEX
logical exogenous
type(LittermanSpec), pointer :: ps
        
T = tData%Nobs
m = tData%m
p = tData%p
q = tData%q
k = tData%k
r = k-q

ps => priorspec%Litterman

allocate( sig(m), Gamma(k,m), Svec(k), ind(p) )

! Setting the means of the prior distribution
Gamma = 0.0d0
do i = 1,m
    Gamma(q+i,i) = ps%firstlagmean(1)
enddo

if (q>0) then
   Gamma(1:q,:) = ps%deterministicmean
end if

do i = 1, m
	sig(i) = tData%ESS(i,i)/dble(T)
end do

exogenous = .false.
!	! Check for Exogeneous varibles among the modelled ones
!	num = 0
!	do i = 1, m
!		enex = endoexo(i)
!		if (enex(1:3) == 'EXO') then
!			num = num+1
!		endif
!	enddo
!
!	if ( num > 0 ) then
!
!		allocate(endovar(m-num),exovar(num))
!
!		endovar = 0
!		exovar = 0
!		num_end = 1
!		num_ex = 1
!		!Set indexes for the endog. and exog. variables. 
!		do i = 1,m
!			enex = endoexo(i)
!			if (enex(1:3) == 'EXO') then
!				exovar(num_ex) = i
!				num_ex = num_ex + 1
!			else
!				endovar(num_end) = i
!				num_end = num_end + 1
!			endif
!		enddo

Svec = 0.0d0
Svec(1:q) = ps%tightness(3)**2
do lag = 1, p
	do i = 1, m
		Svec(q+(lag-1)*m+i) = ps%tightness(1)**2/( sig(i)*lag**(2.0d0*ps%tightness(2)) )
	end do
end do

if ( priorspec%priortype == 'ND' ) then
	! Normal-Diffuse prior

	allocate( tPrior%ND%Gamma(k,m), tPrior%ND%diagvarGamma(m*k), &
	          tPrior%ND%invdiagvarGamma(m*k), tPrior%ND%invvarGammaGamma(k*m),tPrior%ND%priorVecPsi(m))

	tPrior%ND%Gamma   = Gamma
    tprior%ND%priordf = priorspec%priordf
    tprior%ND%diagPsi=.TRUE.
    
	j = 0
	do i = 1, m
		tPrior%ND%diagvarGamma(j+1:j+q) = Svec(1:q)
		tPrior%ND%diagvarGamma(j+q+1:j+k) = Svec(q+1:k)*sig(i)*ps%tightness(4)**2
		! fix own lags
		ind = (/ q+i:q+(p-1)*m+i:m /)
		tPrior%ND%diagvarGamma(j+ind) = Svec(ind)*sig(i) 
		j = j+k
	end do

    if ( exogenous ) then
    
		! Set stder for endog. influence to exog
		do i = 1, num
			do lag = 1, p
				do j = 1, m-num
					tPrior%ND%diagvarGamma(q+(exovar(i)-1)*k+(lag-1)*m+endovar(j))=ps%tightness(5)**2
				enddo
			enddo
		enddo

		deallocate( endovar, exovar )

	endif

    jj = 0
    do j = 1, m
        do i = 1, k 
            jj = jj + 1
            tPrior%ND%invdiagvarGamma(jj)   = 1.0d0/tPrior%ND%diagvarGamma(jj)
            tPrior%ND%invvarGammaGamma(jj) = tPrior%ND%Gamma(i,j)/tPrior%ND%diagvarGamma(jj)
        end do
    end do
    
        if ( tPrior%ND%priordf > m + 1 ) then
    	tprior%ND%priorVecPsi	= dble(priorspec%priordf-m-1)*sig
    endif
    
    
elseif ( priorspec%priortype == 'NW' ) then
	! Normal-Wishart prior

	allocate( tPrior%NW%priorGamma(k,m), tPrior%NW%priorVecOmega(k), &
			  tPrior%NW%priorVecPsi(m))

	tPrior%NW%priorGamma	= Gamma
	tPrior%NW%priorVecOmega = Svec
	tprior%NW%priordf		= priorspec%priordf
    if ( tPrior%NW%priordf > m + 1 ) then
    	tprior%NW%priorVecPsi	= dble(priorspec%priordf-m-1)*sig
    endif
    
    tprior%NW%diagPsi=.TRUE.
    tprior%NW%diagOmega=.TRUE.
   
elseif ( priorspec%priortype == 'SS' ) then
    ! prior for steady state model

    if ( q <= 0 ) then
        print *, "Steady state model meaningless if no deterministic variables"
        stop
    endif
    
	allocate( tPrior%SS%Gamma(m*p,m), tPrior%SS%diagvarGamma(p*m**2), &
	          tPrior%SS%invdiagvarGamma(p*m**2), &
	          tPrior%SS%invvarGammaGamma(p*m**2) , &
	          tPrior%SS%Lambda(m,q), tPrior%SS%diagvarLambda(q*m), &
	          tPrior%SS%invdiagvarLambda(q*m), &
	          tPrior%SS%invvarLambdaLambda(q*m), &
	          tPrior%SS%Psi(m,m) )

    tPrior%SS%Lambda = transpose( Gamma(1:q,:) )
	tPrior%SS%Gamma  = Gamma(q+1:,:)

    ! Svec is variances for parameters in one equation
    ! that is one row in Lambda and one column in Gamma
    if ( allocated( ps%deterministicvar ) ) then
        ! user specified variances
        if ( size(ps%deterministicvar,1) == q .and. size(ps%deterministicvar,2) == m ) then
            call vectorize( transpose(ps%deterministicvar), tPrior%SS%diagvarLambda )
        else
            stop 'Error in specification for variance of deterministic vars'
        endif
    else
        ! use default variance specification
	    do j = 1, q
	        do i = 1, m
		        tPrior%SS%diagvarLambda((j-1)*m+i) = Svec(j)
	        enddo
        enddo
    endif
	j = 0
	do i = 1, m
		tPrior%SS%diagvarGamma(j+1:j+r) = Svec(q+1:k)*sig(i)*ps%tightness(4)**2
		! fix own lags
		ind = (/ i:(p-1)*m+i:m /)
		tPrior%SS%diagvarGamma(j+ind) = Svec(ind+q)*sig(i) 
		j = j+m*p
	end do

    if ( exogenous ) then
		! Set stder for endog. influence to exog
		
		do i = 1, num
			do lag = 1, p
				do j = 1, m-num
					tPrior%SS%diagvarGamma((exovar(i)-1)*k+(lag-1)*m+endovar(j)) &
					    = ps%tightness(5)**2
				enddo
			enddo
		enddo

		deallocate( endovar, exovar )

	endif

    jj = 0
    do j = 1, m
        do i = 1, m*p
            jj = jj + 1
            tPrior%SS%invdiagvarGamma(jj)   = 1.0d0/tPrior%SS%diagvarGamma(jj)
            tPrior%SS%invvarGammaGamma(jj) = tPrior%SS%Gamma(i,j)/tPrior%SS%diagvarGamma(jj)
        end do
    end do

    jj = 0
    do j = 1, q
        do i = 1, m
            jj = jj + 1
            tPrior%SS%invdiagvarLambda(jj)   = 1.0d0/tPrior%SS%diagvarLambda(jj)
            tPrior%SS%invvarLambdaLambda(jj) = tPrior%SS%Lambda(i,j)/tPrior%SS%diagvarLambda(jj)
        end do
    end do

	tprior%SS%priordf		= priorspec%priordf
	tPrior%SS%Psi           = 0.0d0
    if ( tPrior%SS%priordf > m + 1 ) then
	    do i = 1, m
	        tPrior%SS%Psi(i,i) = dble( tprior%SS%priordf - m - 1 )*sig(i)
	    enddo
    endif
    
else
	print *, 'Unknown prior type '//priorspec%priortype
	stop
endif

tPrior%priortype = priorspec%priortype

deallocate( sig, Gamma, Svec, ind )

end subroutine

!================================================================================================

subroutine clean_prior(tPrior)

type(VARPrior), intent(inout) :: tPrior

if ( allocated( tPrior%ND%Gamma ) )            deallocate ( tPrior%ND%Gamma )           
if ( allocated( tPrior%ND%varGamma ) )         deallocate ( tPrior%ND%varGamma )        
if ( allocated( tPrior%ND%invvarGamma ) )      deallocate ( tPrior%ND%invvarGamma )     
if ( allocated( tPrior%ND%diagvarGamma ) )     deallocate ( tPrior%ND%diagvarGamma )    
if ( allocated( tPrior%ND%invdiagvarGamma ) )  deallocate ( tPrior%ND%invdiagvarGamma ) 
if ( allocated( tPrior%ND%invvarGammaGamma ) ) deallocate ( tPrior%ND%invvarGammaGamma )
if ( allocated( tPrior%ND%priorPsi ) )      deallocate ( tPrior%ND%priorPsi )     
if ( allocated( tPrior%ND%priorVecPsi ) )   deallocate ( tPrior%ND%priorVecPsi )  

if ( allocated( tPrior%NW%priorGamma ) )    deallocate ( tPrior%NW%priorGamma )   
if ( allocated( tPrior%NW%priorOmega ) )    deallocate ( tPrior%NW%priorOmega )   
if ( allocated( tPrior%NW%priorVecOmega ) ) deallocate ( tPrior%NW%priorVecOmega )
if ( allocated( tPrior%NW%priorPsi ) )      deallocate ( tPrior%NW%priorPsi )     
if ( allocated( tPrior%NW%priorVecPsi ) )   deallocate ( tPrior%NW%priorVecPsi )  

if ( allocated( tPrior%SS%Gamma ) )              deallocate ( tPrior%SS%Gamma )           
if ( allocated( tPrior%SS%diagvarGamma ) )       deallocate ( tPrior%SS%diagvarGamma )    
if ( allocated( tPrior%SS%invdiagvarGamma ) )    deallocate ( tPrior%SS%invdiagvarGamma ) 
if ( allocated( tPrior%SS%invvarGammaGamma ) )   deallocate ( tPrior%SS%invvarGammaGamma )
if ( allocated( tPrior%SS%Lambda ) )             deallocate ( tPrior%SS%Lambda )          
if ( allocated( tPrior%SS%diagvarLambda ) )      deallocate ( tPrior%SS%diagvarLambda )   
if ( allocated( tPrior%SS%invdiagvarLambda ) )   deallocate ( tPrior%SS%invdiagvarLambda )
if ( allocated( tPrior%SS%invvarLambdaLambda ) ) deallocate( tPrior%SS%invvarLambdaLambda )
if ( allocated( tPrior%SS%Psi ) )                deallocate ( tPrior%SS%Psi )     

end subroutine clean_prior


!===============================================================================================
subroutine initNIW_mh( tData, Gamma_draw, Psi_draw, psi_drawroot)

type(VARData), intent(in) :: tData
real*8, intent(out) :: Gamma_draw(tData%k,tData%m), Psi_draw(tData%m,tData%m), psi_drawroot(tData%m,tData%m)

integer  m

m = tData%m

! start sampler at inverse of OLS error variance

Gamma_draw=0d0

Psi_draw = tData%ESS/dble(tData%Nobs)

call dlftds(m,Psi_draw,m,psi_drawroot,m)

end subroutine


!===============================================================================================

Subroutine sampleNIW_mh(ms, tData, tPrior, Gamma_draw, Psi_draw, psi_drawroot,S_chol,nonselect)

integer, intent(in) :: ms
type(VARData), intent(in) :: tData
type(priorND), intent(in), target :: tPrior
double precision, intent(inout) :: Gamma_draw(:,:),Psi_draw(:,:),psi_drawroot(:,:)
double precision, intent(out), optional :: S_chol(:,:)
logical, intent(in), optional :: nonselect
!type(NIWposterior), intent(out),optional :: tPost                                           

real*8, allocatable :: psi11(:,:),sigma11(:,:),sg11(:),gamma_1(:),gamma_gen(:),gamma_2(:),err(:,:),ess(:,:),s_post(:,:),psiroot(:,:),psi_try(:,:)

real*8 :: alpha, mu

!------

integer T, m, k, i, j, ii, jj
logical diagvarGamma, select

type(priorND), pointer :: prior

prior => tPrior


T = tData%Nobs
m = tData%m
k = tData%k
select=.TRUE.

if (present(nonselect)) then
    select = .not.(nonselect)
end if

if ( allocated( prior%invdiagvarGamma ) ) then
	diagvarGamma = .true.
else
	diagvarGamma = .false.
endif

allocate(psi11(ms,ms),sigma11(k*ms,k*ms),sg11(k*ms),gamma_1(k*ms),gamma_gen(k*ms),gamma_2(k*(m-ms)),err(T,m),ess(m,m),s_post(m,m))
allocate(psiroot(m,m),psi_try(m,m))


!===========================================Draw for Gamma==================================================

psi11=psi_draw(1:ms,1:ms)


IF (MS .EQ. 1) THEN 

    sigma11 = (tdata%zpz)/psi11(1,1)
  
else 

    call linds(psi11,psi11)
    sigma11=KRON(psi11,tdata%ZPZ)
    
end if 



sg11=prior%invvarGammaGamma(1:(k*ms)) + matmul( sigma11, tData%vecGamma_OLS(1:(k*ms)))	  

if ( diagvarGamma ) then
	do i = 1, ms*k
		sigma11(i,i) = prior%invdiagvarGamma(i) &
									 + sigma11(i,i)
	end do
else
	sigma11= prior%invvarGamma + sigma11
endif



! Mean of conditional posterior
! factor inverse posterior variance
call DLFTDS( k*ms, sigma11, k*ms, sigma11, k*ms )	
! solve for posterior mean
call DLFSDS( k*ms, sigma11, k*ms, SG11,Gamma_1 ) 
! Variance of posterior, invert to get factor of variance matrix
call DLINRT( k*ms, sigma11, k*ms, 2, sigma11, k*ms )  

! Draw from MVN distribution
call DRNNOA( k*ms, Gamma_gen )
call DTRMV( 'U', 'N', 'N', k*ms, sigma11, k*ms, Gamma_gen, 1 ) 
ii = 1; jj = k
do j = 1, ms
	Gamma_draw(:,j) = Gamma_1(ii:jj) + Gamma_gen(ii:jj)
	ii = ii+k; jj = jj+k
enddo

if (m .ne. ms) then

   call RMNORM(vector(prior%gamma(:,(ms+1):m)),DDIAG(sqrt(prior%diagvargamma((k*ms+1):(k*m))),k*(m-ms)),gamma_2,1)

      ii=1; jj=k

   do j = (ms+1), m
	   Gamma_draw(:,j) = Gamma_2(ii:jj)
	   ii = ii+k; jj = jj+k
   enddo

end if
!===========================================Draw for Psi==================================================

!*******************************
!Set up for draws of the Wishart
!*******************************
!B_samp = previous B-draw
CALL DMRRRR( T, k, tData%Z, T, k, m, Gamma_draw, k, T, m, err, T )	
err = tData%Y - err			!Y - Z*Gamma
CALL DMXTXF( T, m, err, T, m, ESS, m )	!(Y-ZB)'(Y-ZB)
! Psi ~ iW( YZBYZB'YZBYZB, T )
! Psi^-1 ~ W( (YZBYZB'YZBYZB)^-1, T )
if (prior%diagpsi) then
   do j=1,m
       ESS(j,j)=ESS(j,j)+prior%priorVecPsi(j)
   end do
else 
   ESS=ESS+prior%priorPsi
end if    

! factor, upper triangular
call DLFTDS( m, ESS, m, s_post, m ) 

If (present(s_chol)) then
 S_chol = s_post
end if  


! Psi ~ iW( YESS'ESS, T ), Psi_root'Psi_root ~ iW, PSIi_samp*PSIi_samp' ~ W
call invWishBart( prior%priordf+T, Psiroot, S=s_post) 

call trimult(psiroot,psiroot,psi_try,1,1,1,0)

!====== selection for move =======================

if (select) then 

   call DRNUN(1,mu)

   call mh_alpha(tdata,prior,ms,gamma_draw,psi_draw,psi_try,psi_drawroot,psiroot,S_post,alpha)

   if (mu .le. alpha) then
       psi_draw = psi_try
       psi_drawroot = psiroot
   end if

else 
   psi_draw =psi_try
   psi_drawroot=psiroot
end if 

deallocate(psi11,sigma11,sg11,gamma_1,gamma_gen,gamma_2,err,ess,s_post,psiroot,psi_try)

end subroutine 


!========================================================the probability of move ============================================================

subroutine mh_alpha(tdata,tprior,ms,gamma,psi1,psi2,psiroot1,psiroot2,S_post,alpha,logout)

integer, intent(in) :: ms
type(priorND), intent(in)      :: tPrior
type(VARdata), intent(in)      :: tData
real*8,   intent(in)           :: gamma(:,:),psi1(:,:),S_post(:,:),psiroot1(:,:),psi2(:,:),psiroot2(:,:)
real*8, intent(out) :: alpha
logical, optional, intent(in)  :: logout

real*8 :: phi1,phi2

logical logoff

logoff =.FALSE.

if (present(logout)) then
    logoff = logout
end if 

call mh_lnphi(tdata,tprior,ms,gamma,psi1,psiroot1,S_post,phi1)
call mh_lnphi(tdata,tprior,ms,gamma,psi2,psiroot2,S_post,phi2)

if (logoff) then 
     alpha = phi2-phi1
     if (alpha .gt. 0d0) then
        alpha = 0d0
     end if
else       
     alpha=exp(phi2-phi1)

     if (alpha .ge. 1.0d0) then
         alpha=1.d0
     end if
end if
end subroutine
!======================================================= log phi ================================================================
subroutine mh_lnphi(tdata,tprior,ms,gamma,psi,psiroot,S_post,phi)
integer, intent(in) :: ms
type(priorND), intent(in)      :: tPrior
type(VARdata), intent(in)      :: tData
real*8,   intent(in)           :: gamma(:,:),psi(:,:),S_post(:,:),psiroot(:,:)
real*8, intent(out) :: phi

integer T,m,k
real*8, allocatable :: psi11(:,:),psi1(:,:),err(:)
real*8 :: det1, det2,tmp

double precision, parameter :: ln10 = 2.302585092994045684d0, &
&							   ln2pi = 1.8378770664093454836d0, &
&                              ln2 =   0.6931471805599453094D0,&
&                              lnpi=   1.1447298858494001742D0



T = tData%Nobs
m = tData%m
k = tData%k
tmp=0d0


allocate(psi1(ms,ms),psi11(ms,ms),err(T))

psi11=psi(1:ms,1:ms)


if (ms .eq. 1) then 

    tmp = tmp - log(psi11(1,1))*dble(T)
    err=tdata%y(:,1)-matmul(tdata%z,gamma(:,1))
    tmp=tmp - dot_product(err,err)/psi11(1,1)
    tmp=tmp/2.0d0

else if (ms .gt. 1) then


   call DLFTDS (ms,psi11,ms,psi11,ms)
   CALL DLFDDS ( ms, psi11 , ms, DET1, DET2 )

   tmp=tmp-(log( det1 ) + det2*ln10)*dble(T)


   call DLINRT(ms,psi11,ms,2,psi11,ms)

   call trimult(psi11,psi11,psi11,1,1,0,1)

   call MXTXF(tdata%Y(:,1:ms)-matmul(tdata%z,gamma(:,1:ms)),psi1)

   tmp = tmp - sum(DIAGONALS(matmul(psi11,psi1)))

   tmp = tmp/2.0d0

end if

! p(psi)

tmp = tmp + ln_IW_dens(DDIAG(sqrt(tprior%priorvecpsi),m),psiroot,dble(tprior%priordf),INU1=1,INU2=1)

! p(psi | Y, gamma)

tmp = tmp - ln_IW_dens(S_post(:,:),psiroot,dble(tprior%priordf+T),INU1=1,INU2=1)
 

phi=tmp


deallocate(psi1,psi11,err)


end subroutine 



end module priorspost
