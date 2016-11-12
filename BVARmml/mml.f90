module mml

INCLUDE 'link_fnl_static.h'!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries
use matrix 
use random
use structures
use priorspost

contains

!==============================Chib-Jeliazkov Method=================================

subroutine NIW_chib_Jeliazkov(ms,reps_m,burn_m,tdata,tprior,gamma_draw,psi_draw,psiroot_draw,S_post,lmarglik)

integer, intent(in)            :: ms,reps_m,burn_m
type(priorND), intent(in)      :: tPrior
type(VARdata), intent(in)      :: tData
real*8,   intent(in)           :: gamma_draw(:,:),psi_draw(:,:,:),psiroot_draw(:,:,:),S_post(:,:,:)
real*8,   intent(out)          :: lmarglik

real*8, allocatable :: gamma_mean(:),psiroot_mean(:,:),psi_mean(:,:),gamma(:,:),psi11(:,:),psi2(:,:),err(:), &
                       sigma11(:,:),sg11(:),gamma_1(:),gamma_m(:,:),psi_m(:,:),psiroot_m(:,:),s_chol_m(:,:),alpha(:),ppsi(:)

integer :: i,s,km,m,T,k,reps
real*8 :: tmp,likelihood, det1, det2, dpost,alpha_bar,ppsi_max,lnalpha_max,dpost_gamma1 
logical diagvarGamma

double precision, parameter :: ln10 = 2.302585092994045684d0, &
&							   ln2pi = 1.8378770664093454836d0, &
&                              ln2 =   0.6931471805599453094D0,&
&                              lnpi=   1.1447298858494001742D0

! initialize
 
T=tdata%Nobs
k=tdata%k 
s=size(gamma_draw,2)
m=size(psiroot_draw,1)
km=size(gamma_draw,1)


if ( allocated( tprior%invdiagvarGamma ) ) then
	diagvarGamma = .true.
else
	diagvarGamma = .false.
endif



allocate(gamma_mean(km),psiroot_mean(m,m),gamma(k,m),psi_mean(m,m))

gamma_mean=0d0
psiroot_mean=0d0

gamma_mean=sum(gamma_draw,dim=2)/s

do i = 1,s
    psiroot_mean = psiroot_mean + psiroot_draw(:,:,i)
end do

psiroot_mean = psiroot_mean/s

call trimult(psiroot_mean,psiroot_mean,psi_mean,1,1,1,0)

tmp=0d0


! Prior Part 
!!!Normal: gamma~N(prior_gamma,var_gamma)

tmp=tmp+ln_MVN_dens( gamma_mean-vector(tprior%Gamma), DDIAG(tprior%diagvarGamma,km))

!! Inverse-Wishart : psi~iW(priorS,priorv)

tmp=tmp+ln_IW_dens(DDIAG(sqrt(tprior%priorvecpsi),m),psiroot_mean,dble(tprior%priordf),INU1=1,INU2=1)


!! Likelihood part:

call antivec(k,m,gamma_mean,gamma)

allocate(psi2(ms,ms),psi11(ms,ms),err(T))

psi11=psi_mean(1:ms,1:ms)

likelihood = -dble(T*ms)*ln2pi

if (ms .eq. 1) then 

    likelihood = likelihood - log(psi11(1,1))*dble(T)
    err=tdata%y(:,1)-matmul(tdata%z,gamma(:,1))
    likelihood=likelihood - dot_product(err,err)/psi11(1,1)
    likelihood=likelihood/2.0d0

else 


   call DLFTDS (ms,psi11,ms,psi11,ms)
   CALL DLFDDS ( ms, psi11 , ms, DET1, DET2 )

   likelihood=likelihood-(log( det1 ) + det2*ln10)*dble(T)


   call DLINRT(ms,psi11,ms,2,psi11,ms)

   call trimult(psi11,psi11,psi11,1,1,0,1)

   call MXTXF(tdata%Y(:,1:ms)-matmul(tdata%z,gamma(:,1:ms)),psi2)

   likelihood = likelihood - sum(DIAGONALS(matmul(psi11,psi2)))

   likelihood = likelihood/2.0d0

end if

tmp = tmp + likelihood

!! Posterior part:

dpost=0d0


go to 4500

! p(gamma_mean | psi_mean, Y_i)

allocate(sigma11(K*ms,K*ms),sg11(K*ms),gamma_1(k*ms))

IF (MS .EQ. 1) THEN 

    sigma11 = (tdata%zpz)/psi11(1,1)

else 

    call linds(psi11,psi11)
    sigma11=KRON(psi11,tdata%ZPZ)
    
end if 


sg11=tprior%invvarGammaGamma(1:(k*ms)) + matmul( sigma11, tData%vecGamma_OLS(1:(k*ms)))	  

if ( diagvarGamma ) then
	do i = 1, ms*k
		sigma11(i,i) = tprior%invdiagvarGamma(i) &
									 + sigma11(i,i)
	end do
else
	sigma11= tprior%invvarGamma + sigma11
endif

! Mean of conditional posterior
! factor inverse posterior variance
call linds(sigma11,sigma11) 

! solve for posterior mean
gamma_1 = matmul(sigma11,sg11)

!call DLFTDS( k*ms, sigma11, k*ms, sigma11, k*ms )
	

! gamma_1

dpost = dpost + ln_MVN_dens(vector(gamma(:,1:ms))-gamma_1, sigma11)


! gamma_2

if (m .ne. ms) then

   dpost = dpost + ln_MVN_dens(vector(gamma(:,(1+ms):m))-vector(tprior%gamma(:,(ms+1):m)), DDIAG(sqrt(tprior%diagvargamma((k*ms+1):(k*m))),k*(m-ms)),inu=1)

end if


4500 continue


! avg(p(gamma_mean | psi_i, Y_i))

allocate(sigma11(K*ms,K*ms),sg11(K*ms),gamma_1(k*ms))

dpost_gamma1 = 0d0

do i =1,s  
   
   psi11 = psi_draw(1:ms,1:ms,i)
    
   IF (MS .EQ. 1) THEN 

       sigma11 = (tdata%zpz)/psi11(1,1)

   else 

       call linds(psi11,psi11)
       sigma11=KRON(psi11,tdata%ZPZ)
    
    end if 


    sg11=tprior%invvarGammaGamma(1:(k*ms)) + matmul( sigma11, tData%vecGamma_OLS(1:(k*ms)))	  

    if ( diagvarGamma ) then
	   do j = 1, ms*k
		    sigma11(j,j) = tprior%invdiagvarGamma(j) &
									    + sigma11(j,j)
	   end do
    else
	    sigma11= tprior%invvarGamma + sigma11
    endif

! Mean of conditional posterior
! factor inverse posterior variance
     call linds(sigma11,sigma11) 

! solve for posterior mean
     gamma_1 = matmul(sigma11,sg11)

!call DLFTDS( k*ms, sigma11, k*ms, sigma11, k*ms )
	

! gamma_1

     dpost_gamma1 = dpost_gamma1 + ln_MVN_dens(vector(gamma(:,1:ms))-gamma_1, sigma11)


end do

   
   dpost = dpost + dpost_gamma1/s


! gamma_2

       if (m .ne. ms) then

             dpost = dpost + ln_MVN_dens(vector(gamma(:,(1+ms):m))-vector(tprior%gamma(:,(ms+1):m)), DDIAG(sqrt(tprior%diagvargamma((k*ms+1):(k*m))),k*(m-ms)),inu=1)

      end if


! p(psi_mean | Y_i )



! Prepare draws from p(gamma | psi_mean, Y_i) and p(psi | gamma, Y)

allocate(gamma_m(k,m),psi_m(m,m),psiroot_m(m,m),s_chol_m(m,m),alpha(reps_m),ppsi(s))

reps = reps_m + burn_m

do i = 1,reps
      gamma_m = gamma
      psi_m = psi_mean 
      psiroot_m = psiroot_mean
      call sampleNIW_mh(ms, tData, tPrior, Gamma_m, Psi_m, psiroot_m,s_chol=s_chol_m,nonselect=.TRUE.)
      if (i > burn_m) then
          call mh_alpha(tdata,tprior,ms,gamma_m,psi_mean,psi_m,psiroot_mean,psiroot_m,s_chol_m,alpha(i-burn_m),logout=.TRUE.)
      end if    
end do

lnalpha_max = maxval(alpha)

do i =1,reps_m
      alpha(i) = alpha(i) - lnalpha_max 
      alpha(i) = exp(alpha(i))
end do

alpha_bar=log(sum(alpha/reps_m)) + lnalpha_max 


do i=1,s
   
   call antivec(k,m,gamma_draw(:,i),gamma)
   call mh_alpha(tdata,tprior,ms,gamma,psi_draw(:,:,i),psi_mean,psiroot_draw(:,:,i),psiroot_mean,s_post(:,:,i),ppsi(i),logout=.true.)
   ppsi(i) =ppsi(i)+ln_IW_dens(S_post(:,:,i),psiroot_mean,dble(tprior%priordf+T),INU1=1,INU2=1)- alpha_bar
     
end do

   ppsi_max=maxval(ppsi)
   
do i=1,s
   ppsi(i) = ppsi(i) - ppsi_max 
   ppsi(i) = exp(ppsi(i))
end do   
     
dpost = dpost + log(sum(ppsi)/s) + ppsi_max



! ====== in the end ========

tmp=tmp-dpost

lmarglik = tmp


!==============cleaning =======================
deallocate(gamma_mean,psiroot_mean,psi_mean,gamma,psi11,psi2,err) 
deallocate(sigma11,sg11,gamma_1,gamma_m,psi_m,psiroot_m,s_chol_m,alpha,ppsi)



end subroutine 




!===============================================================

subroutine NIW_MHM_geweke(chi_p,ms,tdata,tprior,gamma_draw,psiroot_draw,lmarglik)

integer, intent(in)            :: ms
real*8, intent(in)             :: chi_p
type(priorND), intent(in)      :: tPrior
type(VARdata), intent(in)      :: tData
real*8,   intent(in)           :: gamma_draw(:,:),psiroot_draw(:,:,:)
real*8,   intent(out)          :: lmarglik


real*8, allocatable :: gamma_mean(:),psiroot_mean(:,:),psi_mean(:,:),gamma(:,:),psi1(:,:),psi11(:,:), &
                       theta(:),vartheta(:,:),theta1(:,:),psiroot_draw1(:,:,:),theta2(:),tmp1(:),psi2(:,:),err(:)

integer :: i,s,km,m,T,k,reps,thn,j
real*8 :: tmp,likelihood, det1, det2, dpost,alpha_bar,ppsi_max,lnalpha_max,dv,fpart,tmp1_max
logical diagvarGamma

double precision, parameter :: ln10 = 2.302585092994045684d0, &
&							   ln2pi = 1.8378770664093454836d0, &
&                              ln2 =   0.6931471805599453094D0,&
&                              lnpi=   1.1447298858494001742D0

! initialize
 
T=tdata%Nobs
k=tdata%k 
s=size(gamma_draw,2)
m=size(psiroot_draw,1)
km=size(gamma_draw,1)
thn=km+(m*(m+1)/2)

if ( allocated( tprior%invdiagvarGamma ) ) then
	diagvarGamma = .true.
else
	diagvarGamma = .false.
endif



allocate(gamma_mean(km),psiroot_mean(m,m),gamma(k,m),psi_mean(m,m))
allocate(theta1(s,thn),vartheta(thn,thn),theta(thn),theta2(thn))
allocate(psi1(m,m),psiroot_draw1(m,m,s))

gamma_mean=0d0
psiroot_mean=0d0

gamma_mean=sum(gamma_draw,dim=2)/s
theta(1:km)=gamma_mean

do i = 1,s
    theta1(i,1:km) = gamma_draw(:,i)- theta
    psiroot_mean = psiroot_mean + psiroot_draw(:,:,i)
    psiroot_draw1(:,:,i) = psiroot_draw(:,:,i)
end do

psiroot_mean = psiroot_mean/s
psi1=psiroot_mean

do i = 1,m
    psi1(i,i) = log(psi1(i,i))
end do 

theta((km+1):thn) = antitri(psi1)

! log-transformation

do i = 1,s   
    do j = 1,m
        psiroot_draw1(j,j,i) = log(psiroot_draw1(j,j,i))
    end do    
    theta1(i,(km+1):thn)=antitri(psiroot_draw1(:,:,i)-psi1)
end do


call mxtxf(theta1,vartheta)

vartheta = vartheta/s

call lftds(vartheta,vartheta)

call lfdds(vartheta,det1,det2)

!call lfdds(vartheta/sqrt(dble(s)),det1,det2)

call dlinrt(thn, vartheta, thn, 2, vartheta, thn)
!vartheta = vartheta*sqrt(dble(s))

allocate(tmp1(s))

tmp1=0d0
fpart=0d0
fpart = fpart - (log( det1 ) + det2*ln10)/2d0
fpart = -(thn/2d0)*ln2pi - log(chi_p)
allocate(psi2(ms,ms),psi11(ms,ms),err(T))

do i = 1,s

     tmp=0d0
     
     theta2=theta1(i,:)
     call dtrmv('U','T','n',thn,vartheta,thn,theta2,1)
      dv=dot_product(theta2,theta2)
     
     if ( dv .le. CHIIN(chi_p,dble(thn))) then           
          
      tmp = fpart - dv/2d0
          
      ! Prior Part 
      !!!Normal: gamma~N(prior_gamma,var_gamma)
 
       tmp=tmp - ln_MVN_dens( gamma_draw(:,i)-vector(tprior%Gamma), DDIAG(tprior%diagvarGamma,km))

      !! Inverse-Wishart : psi~iW(priorS,priorv)

       tmp=tmp - ln_IW_dens(DDIAG(sqrt(tprior%priorvecpsi),m),psiroot_draw(:,:,i),dble(tprior%priordf),INU1=1,INU2=1)
       
       !! Jacobian adjustation due to log-transformation
       
       do j = 1,m  
         tmp = tmp - log(psiroot_draw(j,j,i))
       end do

      !! Likelihood part:

       call antivec(k,m,gamma_draw(:,i),gamma)

        psi11=psiroot_draw(1:ms,1:ms,i)

        likelihood = -dble(T*ms)*ln2pi

          if (ms .eq. 1) then 
              psi11(1,1) = psi11(1,1)**2
             likelihood = likelihood - log(psi11(1,1))*dble(T)
             err=tdata%y(:,1)-matmul(tdata%z,gamma(:,1))
             likelihood=likelihood - dot_product(err,err)/psi11(1,1)
             likelihood=likelihood/2.0d0

         else 

            !call DLFTDS (ms,psi11,ms,psi11,ms)
            CALL DLFDDS ( ms, psi11 , ms, DET1, DET2 )

            likelihood=likelihood-(log( det1 ) + det2*ln10)*dble(T)


            call DLINRT(ms,psi11,ms,2,psi11,ms)

           call trimult(psi11,psi11,psi11,1,1,0,1)

           call MXTXF(tdata%Y(:,1:ms)-matmul(tdata%z,gamma(:,1:ms)),psi2)

           likelihood = likelihood - sum(DIAGONALS(matmul(psi11,psi2)))

           likelihood = likelihood/2.0d0

        end if

      tmp = tmp - likelihood   
      
      tmp1(i) = tmp     
               
     end if 
  
  
end do 


 tmp1_max=MAXVAL(tmp1, MASK=tmp1 .NE. 0d0)

do i=1,s
   if (tmp1(i) .ne. 0d0) then 
     tmp1(i) = tmp1(i) - tmp1_max 
     tmp1(i) = exp(tmp1(i))
   end if
end do 

lmarglik= - (log(sum(tmp1)/s)+ tmp1_max)  

print *, lmarglik

deallocate(gamma_mean,psiroot_mean,psi_mean,gamma,psi1,psi11)
deallocate(theta,vartheta,theta1,psiroot_draw1,theta2,tmp1,psi2,err)

end subroutine 

!===============================================================Laplace method==============================================

subroutine NIW_laplace(mms,tdata,tprior,lmarglik,dfpred) 

use mydata

implicit none 



integer, intent(in)            :: mms
type(priorND), intent(in)      :: tPrior
type(VARdata), intent(in)      :: tData
real*8,   intent(out)          :: lmarglik
real*8, optional, intent(in)   :: dfpred

integer :: N
real*8 :: dfpred1,hmax,tmp,det1,det2

double precision, parameter :: ln10 = 2.302585092994045684d0, &
&							   ln2pi = 1.8378770664093454836d0, &
&                              ln2 =   0.6931471805599453094D0,&
&                              lnpi=   1.1447298858494001742D0


real*8, allocatable :: x(:),Hessian(:,:),Jacobian(:,:),XGUESS1(:),psi(:,:)
integer :: i,j

tmp=0d0

ms=mms
k=tdata%k
m=tdata%m
T=tdata%Nobs
N=m*(m+1)/2

allocate(Y(T,m),Z(T,k),gamma_prior(k,m),vargamma_prior(k*m,k*m),psi_prior(m,m))
Y=tdata%Y
Z=tdata%Z
gamma_prior=tprior%gamma
vargamma_prior=DDiag(tprior%diagvarGamma,k*m)
psi_prior=ddiag(sqrt(tprior%priorVecPsi),m)

DFPRED1 = 0.1d0

if (present(dfpred)) then
  dfpred1 = dfpred
end if  

allocate(x(N),xguess1(N),Hessian(N,N),Jacobian(N,N),psi(m,m))

call linds(tData%ESS/dble(tData%Nobs),psi)
call lftds(psi,psi)

do i =1,m
   psi(i,i)=log(psi(i,i))
end do

Xguess1=antitri(psi)


!call UMCGF(FCN,dfpred1,X,Xguess=xguess1,FVALUE=hmax)

IF (M .EQ. 1) THEN 
!   call UVMIF(FCN,XGUESS1(1),100,X(1))
ELSE
   CALL UMINF(FCN,X,Xguess=xguess1,FVALUE=hmax)
   CALL FDHES(FCN,X,hmax,Hessian)
end if

Jacobian = ddiag(1d0,N)

do i =1,m
   j=i*(i+1)/2
   jacobian(j,j)=exp(X(j))
end do

hessian=matmul(matmul(jacobian,-hessian),jacobian)

call lftds(hessian,hessian)
call lfdds(hessian,det1,det2)

tmp = tmp + dble(N/2)*ln2pi
tmp = tmp - 0.5d0*(log( det1 ) + det2*ln10)
tmp = tmp + hmax*n

lmarglik = tmp

end subroutine 

!================================

subroutine FCN(N,X,F)

use mydata
use matrix
use random

integer N,i
real*8, intent(in)  :: X(N)
real*8, intent(out) :: F

real*8 :: psi(m,m)


! for the diagnoal elements, take the logrithm 

psi=uptri(X,m)

do i =1,m
  psi(i,i)=exp(psi(i,i))
end do


F= 0d0

F=ln_MVN_dens( Y(:,ms)-matmul(Z,gamma_prior(:,ms)),DDIAG(psi(1,1)**2,T)+ matmul(matmul(Z,vargamma_prior(1:(k*ms),1:(ms*k))),transpose(Z)))


F=F + ln_IW_dens(psi_prior,psi,dble(priordf),INU1=1,INU2=1)


! Jacobian adjustment

! for p(Y_i | psi)

F = F + X(1)

! for p(psi)

do i = 1,m
   F = F + log(psi(i,i))
end do

F=F/T

end subroutine

!===============================================================Prior method==============================================

subroutine NIW_prior(ms,tdata,tprior,mc_reps,lnmarglik)


type(VARdata), intent(in)           :: tdata
type(priorND), intent(in)           :: tprior
integer, intent(in)                 :: ms,mc_reps
double precision, intent(out)       :: lnmarglik

integer k,m,i,r, j, u, v, p, q, ii, jj, uu, vv, T
real*8,allocatable :: psi(:,:), psiroot(:,:),s_post(:,:),tmp(:,:), ydiff(:), yvar(:,:), yvarpart(:,:), lnden(:)
double precision maxden


! initialize
 
k=tdata%k 
m=tdata%m
T = tdata%Nobs



    ! generate Psi from the prior and average the marginal data density 
    ! conditional on Psi
    
    ! Prior for Psi
    allocate( S_post(m,m), ydiff(ms*T), Yvar(ms*T,ms*T), Yvarpart(ms*T,ms*T), &
              Psiroot(m,m), Psi(m,m), lnden(mc_reps) )
    if (tprior%diagpsi) then
        ! prior scale matrix for inv Wishart is diagonal, use that structure
        S_post = 0d0
        do i = 1, m
            S_post(i,i) = sqrt( tprior%priorVecPsi(i) )
        end do
    else 
        ! factor, upper triangular
        call DLFTDS( m, tprior%priorPsi, m, S_post, m ) 
    end if    

    ! Marginal data distribution
    ! Y_i|Psi ~ N( (P' kron Z)*Gamma_priormean, 
    !               P'Psi*P kron I + (P' kron Z ) * Gamma_priorvar * ( P kron Z' ) )
    ! Precalculate Ymean = (P' kron Z)*Gamma_priormean
    !              Yvarpart = (P' kron Z ) * Gamma_priorvar * ( P kron Z' )
    ! Assumption is that P picks out the first ms columns in Y
    ydiff = vector( tdata%Y(:,1:ms) - matmul( tdata%Z, tprior%Gamma(:,1:ms) ) )
    yvarpart = 0d0
    
    if ( allocated( tprior%diagvarGamma ) ) then
        allocate( tmp(T,k) )
        u = 1; v = T; p = 0
        do i = 1, ms
            do j = 1, k
                tmp(:,j) = tdata%Z(:,j)*tprior%diagvarGamma(p+j)
            end do
            yvarpart(u:v,u:v) = matmul( tmp, transpose( tdata%Z ) )
            u = u + T; v = v + T; p = p + k
        end do
        deallocate( tmp )
    else
        u = 1; v = T; p = 1; q = T
        ii = 1; jj = k; uu = 1; vv = k
        do i = 1, ms
            do j = 1, ms
                yvarpart(u:v,p:q) = matmul( tdata%Z, matmul( tprior%varGamma(ii:jj,uu:vv), transpose( tdata%Z ) ) )
                u = u + T; v = v + T
                ii = ii + k; jj = jj + k
            end do
            p = p + T; q = q + T
            uu = uu + k; vv = vv + k
        end do
    end if
                
    ! start MC loop
    yvar = yvarpart
    do r = 1, mc_reps
    
        ! generate inv Wishart
        call invWishBart( tprior%priordf, Psiroot, S=S_post )
        call trimult( Psiroot, Psiroot, Psi, Aup=1, Bup=1, Atran=1, Btran=0 )
        
        ! variance matrix for normal
        ! P'Psi*P kron I + (P' kron Z ) * Gamma_priorvar * ( P kron Z' ) )
        u = 0; v = 0;
        do i = 1, ms
            do j = 1, ms
                do ii = 1, T
                    yvar(u+ii,v+ii) = yvarpart(u+ii,v+ii) + Psi(i,j)
                end do
                u = u + T
            end do
            v = v + T
        end do

        ! evaluate normal dist for Yi
        lnden(r) = ln_MVN_dens( ydiff, yvar )
        
    end do
    
    maxden = maxval( lnden )
!    print *, 'min ', minval(lnden), ' max ', maxden
    lnmarglik = log( sum( exp( lnden - maxden ) )/mc_reps ) + maxden
!    print *, 'lnmarglik ', lnmarglik
    deallocate( S_post, ydiff, Yvar, Yvarpart, Psiroot, Psi, lnden )
    
end subroutine     
  
  
  
!===============================================================Prior adjusted by importance sampling method==============================================

subroutine NIW_adjprior(ms,tdata,tprior,mc_burn,mc_reps,lnmarglik,gamma_star)


use priorweight

implicit none 

type(VARdata), intent(in)           :: tdata
type(priorND), intent(in)           :: tprior
!real*8, intent(in)                  :: prior_prop
integer, intent(in)                 :: ms,mc_burn,mc_reps
double precision, intent(out)       :: lnmarglik
real*8, intent(in), optional        :: gamma_star(:,:)

integer k,m,i,r, j, u, v, p, q, ii, jj, uu, vv, T
real*8,allocatable :: psi(:,:), psiroot(:,:),s_post(:,:),s_post2(:,:),tmp(:,:), ydiff(:), yvar(:,:), yvarpart(:,:), lnden(:),gamma_draw(:,:,:),psichol_draw(:,:,:)
real*8,allocatable :: psi1(:,:),psi2(:,:),gamma_mean(:,:),err(:,:),ess(:,:),wpsi(:)
double precision ::  maxden,lpost,lpi,maxwpsi
type(NDsampler) tsampler

! initialize
 
k=tdata%k 
m=tdata%m
T = tdata%Nobs

    allocate( S_post(m,m), S_post2(m,m), ydiff(ms*T), Yvar(ms*T,ms*T), Yvarpart(ms*T,ms*T), &
              Psiroot(m,m), Psi(m,m), lnden(mc_reps),wpsi(mc_reps))
    allocate(psichol_draw(m,m,1),gamma_draw(k,m,1),psi1(m,m),psi2(m,m),gamma_mean(k,m))
    allocate(err(T,m),ess(m,m))
    
    ! generate Psi from the mixture of prior and posterior p(psi | Y, gamma_star)
    ! gamma_star is the sample mean from Gibbs sampler  
    ! then average the marginal data density 
    ! conditional on Psi
          
    ! Prior for Psi

    if (tprior%diagpsi) then
        ! prior scale matrix for inv Wishart is diagonal, use that structure
        S_post = 0d0
        do i = 1, m
            S_post(i,i) = sqrt( tprior%priorVecPsi(i) )
        end do
    else 
        ! factor, upper triangular
        call DLFTDS( m, tprior%priorPsi, m, S_post, m ) 
    end if    

    ! Marginal data distribution
    ! Y_i|Psi ~ N( (P' kron Z)*Gamma_priormean, 
    !               P'Psi*P kron I + (P' kron Z ) * Gamma_priorvar * ( P kron Z' ) )
    ! Precalculate Ymean = (P' kron Z)*Gamma_priormean
    !              Yvarpart = (P' kron Z ) * Gamma_priorvar * ( P kron Z' )
    ! Assumption is that P picks out the first ms columns in Y
    ydiff = vector( tdata%Y(:,1:ms) - matmul( tdata%Z, tprior%Gamma(:,1:ms) ) )
    yvarpart = 0d0
    
    if ( allocated( tprior%diagvarGamma ) ) then
        allocate( tmp(T,k) )
        u = 1; v = T; p = 0
        do i = 1, ms
            do j = 1, k
                tmp(:,j) = tdata%Z(:,j)*tprior%diagvarGamma(p+j)
            end do
            yvarpart(u:v,u:v) = matmul( tmp, transpose( tdata%Z ) )
            u = u + T; v = v + T; p = p + k
        end do
        deallocate( tmp )
    else
        u = 1; v = T; p = 1; q = T
        ii = 1; jj = k; uu = 1; vv = k
        do i = 1, ms
            do j = 1, ms
                yvarpart(u:v,p:q) = matmul( tdata%Z, matmul( tprior%varGamma(ii:jj,uu:vv), transpose( tdata%Z ) ) )
                u = u + T; v = v + T
                ii = ii + k; jj = jj + k
            end do
            p = p + T; q = q + T
            uu = uu + k; vv = vv + k
        end do
    end if

if (present(gamma_star)) then 

   gamma_mean = gamma_star

else 
 
    ! start MC loop
    
    ! Run Gibbs sampler based on p( gamma, psi | Y)
    
   
        call initNIW( tData, tSampler)
         
      
       ! write(*,*) "*** Simulating for predictive likelihood", j, " ***"
       
	
		!  Burning-in period
		 
		 
		do i =1,mc_burn
		 
		 
		     call sampleNIW(tData,tprior,tSampler,gamma_draw,Psiroot_draw=psichol_draw)
		
		
		end do
		
       gamma_mean = 0d0
       
       do i = 1, mc_reps
       
              call sampleNIW(tData,tprior,tsampler,gamma_draw,Psiroot_draw=psichol_draw)
              gamma_mean = gamma_mean + gamma_draw(:,:,1)                      
       end do
 
       gamma_mean = gamma_mean/mc_reps 
 
end if 
 
     !*******************************
      !Compute the posterior S_post
     !*******************************
!B_samp = previous B-draw
CALL DMRRRR( T, k, tData%Z, T, k, m, Gamma_mean, k, T, m, err, T )	
err = tData%Y - err			!Y - Z*Gamma
CALL DMXTXF( T, m, err, T, m, ESS, m )	!(Y-ZB)'(Y-ZB)
! Psi ~ iW( YZBYZB'YZBYZB, T )
! Psi^-1 ~ W( (YZBYZB'YZBYZB)^-1, T )
if (tprior%diagpsi) then
   do j=1,m
       ESS(j,j)=ESS(j,j)+tprior%priorVecPsi(j)
   end do
else 
   ESS=ESS+tprior%priorPsi
end if    

! factor, upper triangular
call DLFTDS( m, ESS, m, s_post2, m ) 
  
 
    yvar = yvarpart
    do r = 1, mc_reps
    
        ! generate inv Wishart from the prior
        call invWishBart( tprior%priordf, Psiroot, S=S_post )
        call trimult( Psiroot, Psiroot, Psi1, Aup=1, Bup=1, Atran=1, Btran=0 )
        
        ! generate inv Wishart from the posterior p(psi | Y, gamma_mean)
        call invWishBart( tprior%priordf + T, Psiroot, S=S_post2 )
        call trimult( Psiroot, Psiroot, Psi2, Aup=1, Bup=1, Atran=1, Btran=0 )
        
        psi = prior_prop*psi1 + (1-prior_prop)*psi2
        
       ! factor, upper triangular
       call DLFTDS( m, psi, m, psiroot, m ) 
       
       ! ln p(psi | Y, gamma)
       lpost=ln_IW_dens(S_post2,psiroot,dble(tprior%priordf+T),INU1=1,INU2=1)
       
       ! ln pi(psi)
       lpi=ln_IW_dens(S_post,psiroot,dble(tprior%priordf),INU1=1,INU2=1)
       
       wpsi(r)=-log(prior_prop + (1-prior_prop)*exp(lpost-lpi))
       
        
        ! variance matrix for normal
        ! P'Psi*P kron I + (P' kron Z ) * Gamma_priorvar * ( P kron Z' ) )
        u = 0; v = 0;
        do i = 1, ms
            do j = 1, ms
                do ii = 1, T
                    yvar(u+ii,v+ii) = yvarpart(u+ii,v+ii) + Psi(i,j)
                end do
                u = u + T
            end do
            v = v + T
        end do

        ! evaluate normal dist for Yi and weighted with w(psi)
        lnden(r) = ln_MVN_dens( ydiff, yvar ) + wpsi(r)  
        
    end do
    
    maxden = maxval( lnden )
    maxwpsi=maxval(wpsi)
!    print *, 'min ', minval(lnden), ' max ', maxden
    lnmarglik = log( sum( exp( lnden - maxden ) )/mc_reps ) + maxden - (log( sum( exp( wpsi - maxwpsi ) )/mc_reps ) + maxwpsi)
!    print *, 'lnmarglik ', lnmarglik

    deallocate( S_post,s_post2, ydiff, Yvar, Yvarpart, Psiroot, Psi, lnden,wpsi,psichol_draw,gamma_draw,psi1,psi2,gamma_mean,err,ess)
   
end subroutine NIW_adjprior        
    

end module mml
