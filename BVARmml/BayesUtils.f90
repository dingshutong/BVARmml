module Bayesutils

INCLUDE 'link_fnl_static.h'!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries
use structures
use matutil
use VARutils
use priorspost
use random

implicit none 
 
contains


!****************************************************************************

subroutine bayes_predlik_NW(tdata,tprior, NUM, EVALvec,mc_reps,posteval,lnpredlik )

implicit none 

! This subroutine calculates the log-predictive likelihood for a model at a specific horizon, with the Normal-Wishart prior 

! in:
! tdata     = The packed data
! tprior    = the Normal-Wishart prior 
! NUM		= Number of origins to calculate the predictive likelihood at
! EVALvec	= Horizon to evaluate predictive likelihood at
! mc_reps   = Number of MCMC replicates to use with density estimation of pred. lik.
! X_pred	= Observations on deterministic variables 
!			  Obs T-num-evalh+1 to T
! posteval	= How to calculate predictive likelihood, 
!			  UPDATE = update parameter posterior as forecast origin is updated 
!			  BEFORE = use parameter posterior for first forecast origin
!			  ALLDATA = use parameter posterior based on all N obs

!Out:
! predlik	= estimate of LOG-predictive likelihood (output)


type(VARdata), intent(in)       :: tData
type(priorNW), intent(in)       :: tprior

integer, intent(in)             :: NUM, mc_reps
integer, intent(in)             :: evalvec(:)
character*12, intent(in)        :: posteval
double precision, intent(out)   :: lnpredlik

type(VARdata) tTempData
type(NWsampler) tsample
type(NWposterior) tpost

integer S,Nobs,m,p,q,k,r,predstart,evalh,updateposterior,TT,j,nrmiss,i
logical raoblackwell,deterinmodel
character*12 tchar
double precision window
double precision, allocatable :: X_pred(:,:),ftemp(:,:),predliks(:), Y_act(:,:), AR(:,:,:), MA(:,:,:), &
              cov(:,:), chvec(:), Y_pred(:,:),Z_for(:,:),psichol_draw(:,:),gamma_draw(:,:)

double precision,dimension(15,num) :: STAT 
double precision, dimension(num) :: dens


!==================initialized========================!
S = MC_reps
Nobs = tData%Nobs
m = tData%m
p = tData%p
q = tData%q
k = tData%k
r = k-q
evalh = maxval(evalvec)
predstart = tData%lastobs-num-evalh+1


allocate(psichol_draw(m,m),gamma_draw(k,m))


raoblackwell = .false.

if (posteval(1:6) == 'UPDATE') then
	updateposterior = 1
	tchar = adjustl(posteval(7:))
	if ( tchar(1:3) == 'RAO' ) then
		raoblackwell = .true.
	endif
elseif(posteval(1:6) == 'BEFORE') then
	updateposterior = 0
	tchar = adjustl(posteval(7:))
	if ( tchar(1:3) == 'RAO' ) then
		raoblackwell = .true.
	endif
    TT = predstart-1
elseif(posteval(1:7) == 'ALLDATA') then
	updateposterior = 0
	tchar = adjustl(posteval(8:))
	if ( tchar(1:3) == 'RAO' ) then
		raoblackwell = .true.
	endif
    TT = tData%lastobs
else
	write(*,*) 'the updating of the posterior is not valid ', posteval
	pause
endif

if ( q > 0 ) then
    deterinmodel = .true.
    allocate( X_pred(evalh+num-1,q) )
    call get_DataMat( tData, 'X', X_pred, firstobs=predstart )
else
    deterinmodel = .false.
endif

allocate( ftemp(evalh,1))


if ( raoblackwell ) then
	allocate( predliks(num), Y_act(num+evalh-1,m), AR(m,m,p), MA(m,m,evalh-1), &
              cov(m*evalh,m*evalh), chvec(size(evalvec)))
	predliks = 0.0d0
	call get_DataMat( tData, 'Y', Y_act, firstobs=predstart)
    chvec = (evalvec-1)*m + 1
else
	allocate( Y_pred(S,num), Y_act(num,m))
	call get_DataMat( tData, 'Y', Y_act, firstobs=predstart)
	! Biweight kernel
	!Window = 2.5d0*dble(s)**(-0.2)
	! normal kernel
	Window = 3.0d0*dble(s)**(-0.2)
endif


if ( updateposterior /= 0 ) then
  ! Rolling forecasts, update posterior for each forecast

	j = 0
	Y_pred = 0.0d0
	allocate( Z_for(1,m*p))

      do TT = predstart, predstart+num-1    !TT is interpreted as "last" observation

		j = j + 1
		  
		call copy_data( tData, tTempData, lastobs=tt )
		call OLS( tTempData )
        call get_DataMat( tData, 'Z', Z_for, firstobs=TT, lastobs=TT )
        
       ! write(*,*) "*** Simulating for predictive likelihood", j, " ***"
                
        call make_NWpost(tprior,tTempData,tpost)
        		
		
		if ( raoblackwell ) then
            
            do i = 1, S
                
            !    if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
            !       write(*,*) 'Gibbs sampler:', i      
            !   endif
     
                    call make_NWsample1(tpost,gamma_draw,psichol_draw)
            
                    if ( deterinmodel ) then
        			    call VARforecast( Z_for(1,:), gamma_draw, ftemp, X=X_pred(j:,:))
        			else
        			    call VARforecast( Z_for(1,:), gamma_draw, ftemp)
        			endif
            
       			call get_AR( m, p, gamma_draw, AR )
			    call MA_polynomial( AR, MA )
			    call make_predcov( MA, psichol_draw, cov )
		        predliks(j) = predliks(j) + &
			        exp( ln_MVN_dens( y_act((j-1)+evalvec,1) - ftemp(evalvec,1), &
			             cov(chvec,chvec) ) )
            enddo
            
		else
		    ! Simulate predictive distribution and do
		    ! density estimation, this is only for one forecast horizon

            do i = 1, S
                !if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
                 !   write(*,*) 'Gibbs sampler:', i      
                !endif
                    call make_NWsample1(tpost,gamma_draw,psichol_draw) 
                     
                    if ( deterinmodel ) then
        			    call VARforecast( Z_for(1,:), gamma_draw, ftemp, X=X_pred(j:,:), &
	           		                      VCROOT=psichol_draw)
                    else
        			    call VARforecast( Z_for(1,:), gamma_draw, ftemp, &
	           		                      VCROOT=psichol_draw)
	           		endif                
    			Y_pred(i,j) = ftemp(evalh,1)
	        enddo
	        
		endif
		
	enddo 


else 

  ! don't update posterior

	call copy_data( tData, tTempData, lastobs=tt )
	call OLS( tTempData )
	allocate( Z_for(num,m*p))
    call get_DataMat( tData, 'Z', Z_for, firstobs=predstart, lastobs=predstart+num-1 )
    
    call make_NWpost(tprior,tTempData,tpost)
    
    
	!write(*,*) "*** Simulating for predictive likelihood ***"

	if ( raoblackwell ) then

        do i = 1, S
            !if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
             !   write(*,*) 'Gibbs sampler:', i      
            !endif
            call make_NWsample1(tpost,gamma_draw,psichol_draw) 
            
   			call get_AR( m, p, gamma_draw, AR )
		    call MA_polynomial( AR, MA )
		    call make_predcov( MA, psichol_draw, cov )
		    
		    do j = 1, num  	
		    	    		    
		            if ( deterinmodel ) then
        			    call VARforecast( Z_for(j,:), gamma_draw, ftemp, X=X_pred(j:,:))
        			else
        			    call VARforecast( Z_for(j,:), gamma_draw, ftemp)
        			endif
		    	        
	                predliks(j) = predliks(j) + &
		            exp( ln_MVN_dens( y_act((j-1)+evalvec,1) - ftemp(evalvec,1), &
		                 cov(chvec,chvec) ) )
            enddo
        enddo

	else

        do i = 1, S
            !if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
            !    write(*,*) 'Gibbs sampler:', i      
            !endif          
            
            call make_NWsample1(tpost,gamma_draw,psichol_draw) 
            
            do j = 1, num
		        
 	                if ( deterinmodel ) then
        			    call VARforecast( Z_for(j,:), gamma_draw, ftemp, X=X_pred(j:,:), &
	           		                      VCROOT=psichol_draw)
                    else
        			    call VARforecast( Z_for(j,:), gamma_draw, ftemp, &
	           		                      VCROOT=psichol_draw)
	           		endif 
	    		Y_pred(i,j) = ftemp(evalh,1)
	        enddo
	        
        enddo
	endif
endif


!***** Summary statistic of the performance ********
lnpredlik = 0.0d0

if ( raoblackwell ) then
	! Rao-Blackwellization, conditional log pred liks
	
	lnpredlik = sum( log( predliks /S) )

else
	! density estimation

	! standardization of variables
	CALL DUVSTA (0,S,num,Y_pred,S,0,0,0,0.0d0,0.0d0,0,STAT,15,NRMISS)
	forall( tt=1:s ) y_pred(tt,1:num) = (y_pred(tt,1:num)-stat(1,1:num))/stat(3,1:num)
	Y_act(:,1) = (y_act(:,1) - stat(1,1:num))/stat(3,1:num)

	do tt = 1,num
		!CALL DDESKN(XKER, S, Y_pred(:,tt), WINDOW, 1.0d0, 1, Y_act(tt), DENS(tt), NrMISS)
		CALL DDESKN(normalden, S, Y_pred(:,tt), WINDOW, -1.0d0, 1, Y_act(tt,1), DENS(tt), NrMISS)
	enddo

	lnpredlik = sum( log( dens(1:num) ) )
	
endif

if ( allocated( Y_pred ) )      deallocate( Y_pred )
if ( allocated( Z_for ) )       deallocate( Z_for )
if ( allocated( ftemp ) )       deallocate( ftemp )
if ( allocated( predliks ) )    deallocate( predliks )
if ( allocated( Y_act ) )       deallocate( Y_act )
if ( allocated( AR ) )          deallocate( AR )
if ( allocated( MA ) )          deallocate( MA )
if ( allocated( cov ) )         deallocate( cov )
if ( allocated( chvec ) )       deallocate( chvec )

deallocate(psichol_draw,gamma_draw)

call clear_data( tTempData )

end subroutine





!=================================================================================================================================

!*******************************************************

subroutine bayes_marglik_NW( Y, ZPZ, Z, tprior, lnmarglik, serie_main)

! Calculate marginal likelihood for Normal-Wishart prior and 
! multivariate regression, the marglikelihood for the y_i variable of interest
!
! Y is matricvariate t with parameters
!	Z*priorB
!	(I + Z*priorOmega*Z')^-1 = I - Z*(priorOmega^-1 + Z'Z)^-1*Z'
!	priorPsi
!	priordf
!
! Y				= matrix of dependent variables, T by m
! Z				= Z matrix T by k
! tprior		= The Normal inverse Whishart prior 
! lnmarglik		= log marginal likelihood (output)
! serie_main    = The study serie, default 1 (optional)


type(priorNW), intent(in)    :: tprior
double precision, intent(in) :: Y(:,:), ZPZ(:,:), Z(:,:)
double precision, intent(out) :: lnmarglik
integer, optional, intent(in) :: serie_main

! local variables
integer T, m, k, i,priordf,ms
logical diagOmega, diagPsi
double precision tmp, det1, det2
double precision, allocatable :: err(:,:), S(:,:), Otemp(:,:), Zerr(:,:),gam(:)

double precision, parameter :: ln10 = 2.302585092994045684d0, &
							   lnpi = 1.1447298858494001741d0

real*8, allocatable :: priorB(:,:), priorOmega(:,:), priorPsi(:,:)

T = size(Y,1)
m = size(Y,2)
k = size(Z,2)
  
if ( size(Z,1) /= T ) then
	print *, "Error in row dimension of Z"
	stop
endif
if ( size(ZPZ,1) /= k .or. size(ZPZ,2) /= k ) then
	print *, "Error in dimension of ZPZ"
	stop
endif
if ( size(tprior%priorGamma,1) /= k .or. size(tprior%priorGamma,2) /= m ) then
	print *, "Error in dimension of priorGamma"
	stop
endif

priordf=tprior%priordf

allocate(priorB(k,m))
priorB=tprior%priorGamma

diagOmega=tprior%diagOmega
diagPsi=tprior%diagPsi

if (diagOmega) then
    allocate(priorOmega(k,1))
    priorOmega(:,1)=tprior%priorVecOmega
else
    allocate(priorOmega(k,k))
    priorOmega=tprior%priorOmega
end if 

if (diagPsi) then
    allocate(priorPsi(m,1))
    priorPsi(:,1)=tprior%priorVecPsi
else
    allocate(priorOmega(m,m))
    priorPsi=tprior%priorPsi
end if 

! m*n selection matrix P

if (present(serie_main)) then
    ms=serie_main 
else 
    ms=1 
end if  

priordf=priordf-ms

allocate( err(T,ms), S(ms,ms), Otemp(k,k), Zerr(k,ms),gam(ms))

lnmarglik = - (dble(T*ms)/2.0d0)*lnpi


if (ms .eq. 1) then 

   ! Scale factor
   ! pi^(-T/2) * |priorPsi|^(priordf/2) * |I + Z*priorOmega*Z'|^(-1/2)
   ! * prod(i=1) ( Gamma((T+priordf+1-i)/2)/Gamma((priordf+1-i)/2) )
   lnmarglik = lnmarglik + (dble(priordf)*log(priorPsi(1,1))/2d0)
   
   ! |I + Z*priorOmega*Z'| = |priorOmega|*|priorOmega^-1 + Z'Z'|
   ! Otemp = factor of priorOmega^-1 + Z'Z'
   
   if ( diagOmega ) then
	 tmp = 0.0d0
	 do i = 1, k
		tmp = tmp + log(priorOmega(i,1))
	 end do
	 lnmarglik = lnmarglik - tmp/2.0d0
	 Otemp = ZPZ
	do i = 1, k
		Otemp(i,i) = Otemp(i,i) + 1.0d0/priorOmega(i,1)
	end do
  else
	CALL DLFTDS( k, priorOmega, k, Otemp, k )
	CALL DLFDDS ( k, Otemp, k, DET1, DET2 )
	lnmarglik = lnmarglik - ( log( det1 ) + det2*ln10 )/2.0d0
	call DLINRT( k, Otemp, k, 2, Otemp, k )
	Otemp = matmul(Otemp,transpose(Otemp)) + ZPZ
  endif
  CALL DLFTDS( k, Otemp, k, Otemp, k )
  CALL DLFDDS ( k, Otemp, k, DET1, DET2 )
  lnmarglik = lnmarglik - ( log( det1 ) + det2*ln10 )/2.0d0
 
  ! Gamma part
  
   gam(1)=ALNGAM(dble(T+priordf)/2D0) - ALNGAM(dble(priordf)/2D0)
   lnmarglik = lnmarglik + gam(1)
   
  ! Determinant
  ! priorPsi + (Y-Z*priorB)'(I + Z*priorOmega*Z')^-1*(Y-Z*priorB) 
  ! = priorPsi + (Y-Z*priorB)'*[I - Z*(priorOmega^-1 + Z'Z)^-1*Z']*(Y-Z*priorB) 
  ! = priorPsi + (Y-Z*priorB)'(Y-Z*priorB) 
  !   - (Y-Z*priorB)'Z*(priorOmega^-1 + Z'Z)^-1*Z'*(Y-Z*priorB)
  
  err(:,1) = Y(:,1) - matmul(Z,priorB(:,1))
  S(1,1)=dot_product(err(:,1),err(:,1))
  S(1,1) = S(1,1) + priorPsi(1,1)
  Zerr(:,1) = matmul( transpose(Z), err(:,1))
  ! Solve, i.e. invert
	CALL DLFSDS( k, Otemp, k, Zerr(:,1), err(1:k,1))
  
  S(1,1) = S(1,1) - dot_product(Zerr(:,1), err(1:k,1))
 
  lnmarglik = lnmarglik - (log(S(1,1)))*dble(T+priordf)/2.0d0
  
else 

   ! Scale factor
   ! pi^(-T*m/2) * |priorPsi|^(priordf/2) * |I + Z*priorOmega*Z'|^(-m/2)
   ! * prod(i=1,m) ( Gamma((T+priordf+1-i)/2)/Gamma((priordf+1-i)/2) )

   if ( diagPsi ) then
	  tmp = 0.0d0
	  do i = 1, ms
		tmp = tmp + log(priorPsi(i,1))
	  end do
	  lnmarglik = lnmarglik + dble(priordf)*tmp/2.0d0
   else
	  CALL DLFTDS( ms, priorPsi(1:ms,1:ms), ms, S, ms)
	  CALL DLFDDS ( ms, S, ms, DET1, DET2 )
	lnmarglik = lnmarglik + ( log( det1 ) + det2*ln10 )*dble(priordf)/2.0d0
   endif
   ! |I + Z*priorOmega*Z'| = |priorOmega|*|priorOmega^-1 + Z'Z'|
   ! Otemp = factor of priorOmega^-1 + Z'Z'
   if ( diagOmega ) then
	 tmp = 0.0d0
	 do i = 1, k
		tmp = tmp + log(priorOmega(i,1))
	 end do
	 lnmarglik = lnmarglik - dble(ms)*tmp/2.0d0
	 Otemp = ZPZ
	do i = 1, k
		Otemp(i,i) = Otemp(i,i) + 1.0d0/priorOmega(i,1)
	end do
  else
	CALL DLFTDS( k, priorOmega, k, Otemp, k )
	CALL DLFDDS ( k, Otemp, k, DET1, DET2 )
	lnmarglik = lnmarglik - ( log( det1 ) + det2*ln10 )*dble(ms)/2.0d0
	call DLINRT( k, Otemp, k, 2, Otemp, k )
	Otemp = matmul(Otemp,transpose(Otemp)) + ZPZ
  endif
  CALL DLFTDS( k, Otemp, k, Otemp, k )
  CALL DLFDDS ( k, Otemp, k, DET1, DET2 )
  lnmarglik = lnmarglik - ( log( det1 ) + det2*ln10 )*dble(ms)/2.0d0

  ! Gamma functions

  do i = 1, ms
    gam(I)=ALNGAM(dble(T+priordf+1-i)/2D0)-ALNGAM(dble(priordf+1-i)/2D0)
  end do

  lnmarglik = lnmarglik + sum(gam)

  ! Determinant
  ! priorPsi + (Y-Z*priorB)'(I + Z*priorOmega*Z')^-1*(Y-Z*priorB) 
  ! = priorPsi + (Y-Z*priorB)'*[I - Z*(priorOmega^-1 + Z'Z)^-1*Z']*(Y-Z*priorB) 
  ! = priorPsi + (Y-Z*priorB)'(Y-Z*priorB) 
  !   - (Y-Z*priorB)'Z*(priorOmega^-1 + Z'Z)^-1*Z'*(Y-Z*priorB)
  err = Y(:,1:ms) - matmul(Z,priorB(:,1:ms))
  S   = matmul( transpose(err),err )
  if ( diagPsi ) then
	do i = 1, ms
		S(i,i) = S(i,i) + priorPsi(i,1)
	end do
  else
	S = S + priorPsi(1:ms,1:ms)
  endif
  Zerr = matmul( transpose(Z), err )
  ! Solve, i.e. invert
  do i = 1, ms
	CALL DLFSDS( k, Otemp, k, Zerr(:,i), err(1:k,i) )
  end do
  S = S - matmul( transpose(Zerr), err(1:k,1:ms) )
  CALL DLFTDS( ms, S, ms, S, ms )
  CALL DLFDDS( ms, S, ms, DET1, DET2 )

  lnmarglik = lnmarglik - ( log( det1 ) + det2*ln10 )*dble(T+priordf)/2.0d0

end if 

deallocate( err, S, Otemp, Zerr,gam,priorB,priorOmega,priorPsi)

end subroutine


!==========================================================================================================
SUBROUTINE bayes_lnpredlik_NIW(tData,tPrior,NUM,EVALvec,mc_reps,mc_burn,posteval,lnpredlik)


implicit none 

! This subroutine calculates the log-predictive and for a model at a specific horizon, with the Normal-Independent Wishart prior 

! in:
! tdata     = The packed data
! tprior    = the Normal-Independent Wishart prior 
! NUM		= Number of origins to calculate the predictive likelihood at
! EVALvec	= Horizon to evaluate predictive likelihood at
! mc_reps   = Number of MCMC replicates to use with density estimation of pred. lik.
! mc_burn   = Number of MCMC burns to use with density estimation of pred. lik.
! X_pred	= Observations on deterministic variables 
!			  Obs T-num-evalh+1 to T
! posteval	= How to calculate predictive likelihood, 
!			  UPDATE = update parameter posterior as forecast origin is updated 
!			  BEFORE = use parameter posterior for first forecast origin
!			  ALLDATA = use parameter posterior based on all N obs

!Out:
!lnpredlik	= estimate of LOG-predictive likelihood (output)
 

type(VARdata), intent(in)       :: tData
type(priorND), intent(in)       :: tprior
integer, intent(in)             :: NUM, mc_reps,mc_burn
integer, intent(in)             :: evalvec(:)
character*12, intent(in)        :: posteval
double precision, intent(out)   :: lnpredlik

type(VARdata) tTempData
type(NDsampler) tsampler


integer S,B,Nobs,m,p,q,k,r,predstart,evalh,updateposterior,TT,j,nrmiss,i
logical raoblackwell,deterinmodel
character*12 tchar
double precision window
double precision, allocatable :: X_pred(:,:),ftemp(:,:),predliks(:), Y_act(:,:), AR(:,:,:), MA(:,:,:), &
              cov(:,:), chvec(:), Y_pred(:,:),Z_for(:,:),psichol_draw(:,:,:),gamma_draw(:,:,:)

double precision,dimension(15,num) :: STAT 
double precision, dimension(num) :: dens


!==================initialized========================!
S = MC_reps
B = Mc_burn
Nobs = tData%Nobs
m = tData%m
p = tData%p
q = tData%q
k = tData%k
r = k-q
evalh = maxval(evalvec)
predstart = tData%lastobs-num-evalh+1


allocate(psichol_draw(m,m,1),gamma_draw(k,m,1))

raoblackwell = .false.

if (posteval(1:6) == 'UPDATE') then
	updateposterior = 1
	tchar = adjustl(posteval(7:))
	if ( tchar(1:3) == 'RAO' ) then
		raoblackwell = .true.
	endif
elseif(posteval(1:6) == 'BEFORE') then
	updateposterior = 0
	tchar = adjustl(posteval(7:))
	if ( tchar(1:3) == 'RAO' ) then
		raoblackwell = .true.
	endif
    TT = predstart-1
elseif(posteval(1:7) == 'ALLDATA') then
	updateposterior = 0
	tchar = adjustl(posteval(8:))
	if ( tchar(1:3) == 'RAO' ) then
		raoblackwell = .true.
	endif
    TT = tData%lastobs
else
	write(*,*) 'the updating of the posterior is not valid ', posteval
	pause
endif

if ( q > 0 ) then
    deterinmodel = .true.
    allocate( X_pred(evalh+num-1,q) )
    call get_DataMat( tData, 'X', X_pred, firstobs=predstart )
else
    deterinmodel = .false.
endif

allocate( ftemp(evalh,1))


if ( raoblackwell ) then
	allocate( predliks(num), Y_act(num+evalh-1,m), AR(m,m,p), MA(m,m,evalh-1), &
              cov(m*evalh,m*evalh), chvec(size(evalvec)))
	predliks = 0.0d0
	call get_DataMat( tData, 'Y', Y_act, firstobs=predstart)
    chvec = (evalvec-1)*m + 1
else
	allocate( Y_pred(S,num), Y_act(num,m))
	call get_DataMat( tData, 'Y', Y_act, firstobs=predstart)
	! Biweight kernel
	!Window = 2.5d0*dble(s)**(-0.2)
	! normal kernel
	Window = 3.0d0*dble(s)**(-0.2)
endif


if ( updateposterior /= 0 ) then
  ! Rolling forecasts, update posterior for each forecast

	j = 0
	Y_pred = 0.0d0
	allocate( Z_for(1,m*p))

      do TT = predstart, predstart+num-1    !TT is interpreted as "last" observation

		j = j + 1
		  
		call copy_data( tData, tTempData, lastobs=tt )
		call OLS( tTempData )
        call get_DataMat( tData, 'Z', Z_for, firstobs=TT, lastobs=TT )
        call initNIW( tTempData, tSampler)
       
      
       ! write(*,*) "*** Simulating for predictive likelihood", j, " ***"
       
		
		
		!  Burning-in period
		 
		 
		do i =1,b
		 
		 
		     call sampleNIW(tTempData,tprior,tsampler,gamma_draw,Psiroot_draw=psichol_draw)
		
		
		end do
		
		
		
		if ( raoblackwell ) then
            
            do i = 1, S
                
            !    if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
            !       write(*,*) 'Gibbs sampler:', i      
            !   endif
     
                    call sampleNIW(tTempData,tprior,tsampler,gamma_draw,Psiroot_draw=psichol_draw)
            
                    if ( deterinmodel ) then
        			    call VARforecast( Z_for(1,:), gamma_draw(:,:,1), ftemp, X=X_pred(j:,:))
        			else
        			    call VARforecast( Z_for(1,:), gamma_draw(:,:,1), ftemp)
        			endif
            
       			call get_AR( m, p, gamma_draw(:,:,1), AR )
			    call MA_polynomial( AR, MA )
			    call make_predcov( MA, psichol_draw(:,:,1), cov )

		        predliks(j) = predliks(j) + &
			        exp( ln_MVN_dens( y_act((j-1)+evalvec,1) - ftemp(evalvec,1), &
			             cov(chvec,chvec) ) )         
            enddo
            
		else
		    ! Simulate predictive distribution and do
		    ! density estimation, this is only for one forecast horizon

            do i = 1, S
                !if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
                 !   write(*,*) 'Gibbs sampler:', i      
                !endif
                    call sampleNIW(tTempData,tprior,tsampler,gamma_draw,Psiroot_draw=psichol_draw)
                     
                    if ( deterinmodel ) then
        			    call VARforecast( Z_for(1,:), gamma_draw(:,:,1), ftemp, X=X_pred(j:,:), &
	           		                      VCROOT=psichol_draw(:,:,1))
                    else
        			    call VARforecast( Z_for(1,:), gamma_draw(:,:,1), ftemp, &
	           		                      VCROOT=psichol_draw(:,:,1))
	           		endif                
    			Y_pred(i,j) = ftemp(evalh,1)
	        enddo
	        
		endif
		
	enddo 


else 

  ! don't update posterior

	call copy_data( tData, tTempData, lastobs=tt )
	call OLS( tTempData )
	allocate( Z_for(num,m*p) )
    call get_DataMat( tData, 'Z', Z_for, firstobs=predstart, lastobs=predstart+num-1 )
    call initNIW( tTempData, tSampler)
         
	!write(*,*) "*** Simulating for predictive likelihood ***"
	
	! Burning -in preiod
	do i =1 , b 
	
	  call sampleNIW(tTempData,tprior,tsampler,gamma_draw,Psiroot_draw=psichol_draw)
	  
	end do

	if ( raoblackwell ) then

        do i = 1, S
            !if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
             !   write(*,*) 'Gibbs sampler:', i      
            !endif
            call sampleNIW(tTempData,tprior,tsampler,gamma_draw,Psiroot_draw=psichol_draw)
            
   			call get_AR( m, p, gamma_draw(:,:,1), AR )
		    call MA_polynomial( AR, MA )
		    call make_predcov( MA, psichol_draw(:,:,1), cov )
		    
		    do j = 1, num  	
		    	    		    
		            if ( deterinmodel ) then
        			    call VARforecast( Z_for(j,:), gamma_draw(:,:,1), ftemp, X=X_pred(j:,:))
        			else
        			    call VARforecast( Z_for(j,:), gamma_draw(:,:,1), ftemp)
        			endif  
	                  predliks(j) = predliks(j) + &
		              exp( ln_MVN_dens( y_act((j-1)+evalvec,1) - ftemp(evalvec,1), &
		                  cov(chvec,chvec) ) )      
            enddo
        enddo

	else

        do i = 1, S
            !if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
            !    write(*,*) 'Gibbs sampler:', i      
            !endif          
            
            call sampleNIW(tTempData,tprior,tsampler,gamma_draw,Psiroot_draw=psichol_draw)
            
            do j = 1, num
		        
 	                if ( deterinmodel ) then
        			    call VARforecast( Z_for(j,:), gamma_draw(:,:,1), ftemp, X=X_pred(j:,:), &
	           		                      VCROOT=psichol_draw(:,:,1))
                    else
        			    call VARforecast( Z_for(j,:), gamma_draw(:,:,1), ftemp, &
	           		                      VCROOT=psichol_draw(:,:,1))
	           		endif 
	    		Y_pred(i,j) = ftemp(evalh,1)
	        enddo
	        
        enddo
	endif
endif


!***** Summary statistic of the performance ********
lnpredlik = 0.0d0

if ( raoblackwell ) then
	! Rao-Blackwellization, conditional log pred liks
	
	lnpredlik = sum( log( predliks /s) )

else
	! density estimation

	! standardization of variables
	CALL DUVSTA (0,S,num,Y_pred,S,0,0,0,0.0d0,0.0d0,0,STAT,15,NRMISS)
	forall( tt=1:s ) y_pred(tt,1:num) = (y_pred(tt,1:num)-stat(1,1:num))/stat(3,1:num)
	Y_act(:,1) = (y_act(:,1) - stat(1,1:num))/stat(3,1:num)

	do tt = 1,num
		!CALL DDESKN(XKER, S, Y_pred(:,tt), WINDOW, 1.0d0, 1, Y_act(tt), DENS(tt), NrMISS)
		CALL DDESKN(normalden, S, Y_pred(:,tt), WINDOW, -1.0d0, 1, Y_act(tt,1), DENS(tt), NrMISS)
	enddo

	lnpredlik = sum( log( dens(1:num) ) )
	
endif

if ( allocated( Y_pred ) )      deallocate( Y_pred )
if ( allocated( Z_for ) )       deallocate( Z_for )
if ( allocated( ftemp ) )       deallocate( ftemp )
if ( allocated( predliks ) )    deallocate( predliks )
if ( allocated( Y_act ) )       deallocate( Y_act )
if ( allocated( AR ) )          deallocate( AR )
if ( allocated( MA ) )          deallocate( MA )
if ( allocated( cov ) )         deallocate( cov )
if ( allocated( chvec ) )       deallocate( chvec )

deallocate(psichol_draw,gamma_draw)

call clear_data( tTempData )


end subroutine bayes_lnpredlik_NIW


!==============================================================================================================

subroutine bayes_forecasts_NW(tdata,tprior,EVALh,mc_reps,fcasts,non_deterfcasts)

type(VARdata), intent(in)           :: tdata
type(priorNW), intent(in)          :: tprior
integer, intent(in)                 :: mc_reps
integer, intent(in)                 :: evalh
double precision, intent(out)       :: fcasts(:,:)
integer, optional, intent(in)       :: non_deterfcasts 


type(VARdata) tTempData
type(NWposterior) tpost

integer m,p,q,tt,k,predstart,dfpost,ind_ndeterfcasts,i
logical deterinmodel

double precision, allocatable :: X_pred(:,:),ftemp(:,:),Z_for(:,:),psichol_draw(:,:),gamma_draw(:,:)



!==================initialized========================!

m = tData%m
p = tData%p
q = tData%q
k = tData%k
ind_ndeterfcasts=0
TT = tData%lastobs

fcasts=0d0

!predstart = tData%lastobs-num-evalh+1

allocate(psichol_draw(m,m),gamma_draw(k,m),ftemp(evalh,m))

if (present(non_deterfcasts)) then 
     ind_ndeterfcasts=non_deterfcasts
end if 

if ( q > 0 ) then
    deterinmodel = .true.
    allocate( X_pred(evalh,q) )
    call get_DataMat( tData, 'X', X_pred, firstobs=tt-evalh+1)
else
    deterinmodel = .false.
endif


	call copy_data( tData, tTempData, lastobs=tt)
	call OLS( tTempData )
	allocate( Z_for(1,m*p) )
    call get_DataMat( tData, 'Z', Z_for, firstobs=TT, lastobs=TT) 
    call make_NWpost(tprior,tdata,tpost)
    
    
	!write(*,*) "*** Simulating for predictive likelihood ***"

	if (ind_ndeterfcasts == 0) then

        do i = 1, mc_reps
            !if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
             !   write(*,*) 'Gibbs sampler:', i      
            !endif
            call make_NWsample1(tpost,gamma_draw,psichol_draw) 
            		        
		            if ( deterinmodel ) then
        			    call VARforecast( Z_for(1,:), gamma_draw, ftemp, X=X_pred)
        			else
        			    call VARforecast( Z_for(1,:), gamma_draw, ftemp)
        			endif
            fcasts=fcasts+ftemp		    	        
        enddo
       
	else 

        do i = 1, mc_reps
            !if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
            !    write(*,*) 'Gibbs sampler:', i      
            !endif          
            
            call make_NWsample1(tpost,gamma_draw,psichol_draw) 

		        
 	                if ( deterinmodel ) then
        			    call VARforecast( Z_for(1,:), gamma_draw, ftemp, X=X_pred, VCROOT=psichol_draw)                      
                    else
        			    call VARforecast( Z_for(1,:), gamma_draw, ftemp,VCROOT=psichol_draw)
	           		endif 
	         fcasts=fcasts+ftemp	   		
        enddo
	endif
fcasts=fcasts/mc_reps

if (allocated(X_pred)) deallocate(X_pred)
deallocate(ftemp,Z_for,psichol_draw,gamma_draw)

end subroutine 


!==========================================================================================================

subroutine bayes_forecasts_NIW(tdata,tprior,EVALh,mc_burn,mc_reps,fcasts,non_deterfcasts)

implicit none 

type(VARdata), intent(in)           :: tdata
type(priorND), intent(in)           :: tprior
integer, intent(in)                 :: mc_burn,mc_reps
integer, intent(in)                 :: evalh
double precision, intent(out)       :: fcasts(:,:)
integer, optional, intent(in)       :: non_deterfcasts 

type(NDSampler)  tSampler
type(VARdata)   tTempData
type(NIWposterior) tPost

integer m,p,q,tt,k,predstart,dfpost,ind_ndeterfcasts,i
logical deterinmodel

double precision, allocatable :: X_pred(:,:),ftemp(:,:),Z_for(:,:), psichol_draw(:,:,:),gamma_draw(:,:,:)  


!==================initialized========================!

m = tData%m
p = tData%p
q = tData%q
k = tData%k
ind_ndeterfcasts=0
TT = tData%lastobs

fcasts=0d0

!predstart = tData%lastobs-num-evalh+1
                            
allocate(psichol_draw(m,m,1),gamma_draw(k,m,1),ftemp(evalh,m))

if (present(non_deterfcasts)) then 
     ind_ndeterfcasts=non_deterfcasts
end if 

if ( q > 0 ) then
    deterinmodel = .true.
    allocate( X_pred(evalh,q) )
    call get_DataMat( tData, 'X', X_pred, firstobs=tt-evalh+1)
else
    deterinmodel = .false.
endif


	call copy_data( tData, tTempData, lastobs=tt)
	call OLS( tTempData )
	allocate( Z_for(1,m*p) )
    call get_DataMat( tData, 'Z', Z_for, firstobs=TT, lastobs=TT)
    call initNIW(tTempData, tSampler)
    call initNIWpost(tTempData,tPost)    
	!write(*,*) "*** Simulating for predictive likelihood ***"
	
	do i = 1, mc_burn
	
	     call sampleNIW( tTempData, tPrior, tSampler, Gamma_draw)   
	
	end do 

	if (ind_ndeterfcasts == 0) then

        do i = 1, mc_reps
            !if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
             !   write(*,*) 'Gibbs sampler:', i      
            !endif
            call sampleNIW( tTempData, tPrior, tSampler, Gamma_draw,Psiroot_draw=psichol_draw)           	        
		       
		            if ( deterinmodel ) then
        			    call VARforecast( Z_for(1,:), gamma_draw(:,:,1), ftemp, X=X_pred)
        			else
        			    call VARforecast( Z_for(1,:), gamma_draw(:,:,1), ftemp)
        			endif
        			
 		
               fcasts=fcasts+ftemp
  		
                          	        
        enddo
        
        
       
	else 

        do i = 1, mc_reps
            !if (floor(dble(i)/1000.0D0) == ceiling(dble(i)/1000.0D0)) then
            !    write(*,*) 'Gibbs sampler:', i      
            !endif          
            
            call sampleNIW( tTempData, tPrior, tSampler, Gamma_draw,Psiroot_draw=psichol_draw)  

 	                if ( deterinmodel ) then
        			    call VARforecast( Z_for(1,:), gamma_draw(:,:,1), ftemp, X=X_pred, VCROOT=psichol_draw(:,:,1))                      
                    else
        			    call VARforecast( Z_for(1,:), gamma_draw(:,:,1), ftemp,VCROOT=psichol_draw(:,:,1))
	           		endif 
	           				
		
               fcasts=fcasts+ftemp
 
            
        enddo
	endif
	
fcasts=fcasts/mc_reps

if (allocated(X_pred)) deallocate(X_pred)

deallocate(ftemp,Z_for,psichol_draw,gamma_draw)


end subroutine



!=================================================================================================================================



subroutine bayes_margins_NIW(tdata,tprior,ms,mc_burn,mc_reps,burn_m,reps_m,method,lnmarglik)

use mml

implicit none

type(VARdata), intent(in)           :: tdata
type(priorND), intent(in)           :: tprior
integer, intent(in)                 :: ms,mc_burn,mc_reps,burn_m,reps_m
character*20, intent(in)            :: method 
double precision, intent(out)       :: lnmarglik

integer k,m,km,s,i,T
!integer r, j, u, v, p, q, ii, jj, uu, vv, T
character*12 tchar

real*8, allocatable :: gamma(:,:), psi(:,:), psiroot(:,:),s_post(:,:),s_chol_draws(:,:,:),gamma_draws(:,:),psi_draws(:,:,:),psiroot_draws(:,:,:)

tchar=adjustl(method)


! initialize
 
k=tdata%k 
m=tdata%m
s=mc_reps
km=k*m
T = tdata%Nobs

if ((tchar(1:4) .eq. "CHIB") .or. (tchar(1:6) .eq. "GEWEKE")) then

    allocate(gamma(k,m),psi(m,m),psiroot(m,m),s_post(m,m),s_chol_draws(m,m,s),gamma_draws(km,s),psi_draws(m,m,s),psiroot_draws(m,m,s))

    ! Initial the samplers
    call initNIW_mh( tData, Gamma, Psi, psiroot)

    do i = 1, mc_burn
        call  sampleNIW_mh(ms, tData, tPrior, Gamma, Psi, psiroot)
    end do

    ! Do the samplers
    do i = 1, mc_reps
        call  sampleNIW_mh(ms, tData, tPrior, Gamma, Psi, psiroot,S_chol=s_post)
        ! save for the marginal comp
        s_chol_draws(:,:,i)=s_post 
        gamma_draws(:,i)=vector(gamma)
        psi_draws(:,:,i)=psi
        psiroot_draws(:,:,i)=psiroot
    end do
   
    if (tchar(1:4) .eq. "CHIB") then
 
        call NIW_chib_Jeliazkov(ms,reps_m,burn_m,tdata,tprior,gamma_draws,psi_draws,psiroot_draws,s_chol_draws,lnmarglik)
    
    else if (tchar(1:6) .eq. "GEWEKE") then    

        call NIW_MHM_geweke(0.1d0,ms,tdata,tprior,gamma_draws,psiroot_draws,lnmarglik)    
    
    end if 
    
    
    deallocate(gamma,psi,psiroot,s_post,s_chol_draws,gamma_draws,psi_draws,psiroot_draws)

else if ( tchar(1:5) .eq. 'PRIOR' ) then

    call NIW_prior(ms,tdata,tprior,mc_reps,lnmarglik)
    
else if ( tchar(1:8) .eq. 'ADJPRIOR' ) then
   
   ! print *, "Input for the proportion of priors against the posterior"
   ! read *, prior_weight
    call NIW_adjprior(ms,tdata,tprior,mc_burn,mc_reps,lnmarglik)
   
else if (tchar(1:7) .eq. "LAPLACE") then
    
    call NIW_laplace(ms,tdata,tprior,lnmarglik,dfpred=0.01d0) 
    

end if 


end subroutine







! *********************************************************

subroutine chibs_margins(tdata,tprior,gamma_draw,psiroot_draw,gamma_post,S_post,vargamma_post,lmarglik)

type(priorND), intent(in)      :: tPrior
type(VARdata), intent(in)      :: tData
real*8,   intent(in)           :: gamma_draw(:,:),psiroot_draw(:,:,:),gamma_post(:,:),vargamma_post(:,:,:),S_post(:,:,:)
real*8,   intent(out)          :: lmarglik

real*8, allocatable :: gamma_mean(:),psiroot_mean(:,:),gamma(:,:),psiroot_mean_inv(:,:),psi1(:,:)

integer :: i,s,km,m,T,k
real*8 :: tmp,likelihood, det1, det2, epost

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

allocate(gamma_mean(km),psiroot_mean(m,m),gamma(k,m),psiroot_mean_inv(m,m),psi1(m,m))

gamma_mean=0d0
psiroot_mean=0d0

do i=1,s 
    gamma_mean=gamma_mean+gamma_draw(:,i)
    psiroot_mean=psiroot_mean+ psiroot_draw(:,:,i)
end do

gamma_mean=gamma_mean/s
psiroot_mean=psiroot_mean/s


tmp=0d0


! Prior Part 
!!!Normal: gamma~N(prior_gamma,var_gamma)

tmp=tmp+ln_MVN_dens( gamma_mean-vector(tprior%Gamma), DDIAG(tprior%diagvarGamma,km))

!! Inverse-Wishart : psi~iW(priorS,priorv)

tmp=tmp+ln_IW_dens(DDIAG(sqrt(tprior%priorvecpsi),m),psiroot_mean,dble(tprior%priordf),INU1=1,INU2=1)

!! Likelihood part:

call antivec(k,m,gamma_mean,gamma)

likelihood = -dble(T*k)*ln2pi


CALL DLFDDS ( m, psiroot_mean , m, DET1, DET2 )

likelihood = likelihood  - (log( det1 ) + det2*ln10)*dble(T)

call DLINRT(m,psiroot_mean,m,2,psiroot_mean_inv,m)

call trimult(psiroot_mean_inv,psiroot_mean_inv,psiroot_mean_inv,1,1,0,1)

call MXTXF(tdata%Y-matmul(tdata%z,gamma),psi1)

likelihood = likelihood - sum(DIAGONALS(matmul(psiroot_mean_inv,psi1)))

likelihood = likelihood/2.0d0


tmp = tmp + likelihood

!! Posterior part:

epost=0d0

do i=1,s
   epost=epost+exp(ln_MVN_dens( gamma_mean-gamma_post(:,i), vargamma_post(:,:,i)) & 
    + ln_IW_dens(S_post(:,:,i),psiroot_mean,dble(tprior%priordf+T),INU1=1,INU2=1))
end do

tmp=tmp-log(epost/s)

lmarglik = tmp

end subroutine chibs_margins


!====================================================
!***************************************************

  function xker(y)
  
  ! Biweight kernel, should return zero if abs(y) > 1
  ! but we tell DDESKN to ignore these values anyway
  double precision xker
  double precision, intent(in) :: y
  if ( abs(y) > 1d0 ) then
	print *, y
  endif
  xker = ( 15d0*( 1d0 - y**2 )**2 )/16d0

  end function

  function normalden(y)

  ! normal kernel
  double precision :: normalden
  double precision, intent(in) :: y
  double precision, parameter :: i_sqrt_TwoPi = 0.39894228040143267794d0

!  normalden = exp( -(y**2)/2d0 )/sqrt(2d0*dconst("pi"))		!normal-approx
  normalden = exp( -(y**2)/2d0 )*i_sqrt_TwoPi		!normal-approx

  end function

  function logtden( x, mean, s, p, df )

  ! log of univariate t density, notation as in BLR
  double precision logtden
  double precision, intent(in) :: x, mean, s, p
  integer, intent(in) :: df

  double precision temp
  double precision, parameter :: logpi = 1.1447298858494001741d0

  temp = dble(df)*log(s) + log(p) - logpi + D_POCH( dble(df), 0.5d0 ) &
         - dble(df+1)*log( s + p*(x-mean)**2 )

  logtden = temp/2d0

  end function


end module