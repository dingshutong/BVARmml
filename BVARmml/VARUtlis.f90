module VARutils


INCLUDE 'link_fnl_static.h'!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries

use structures

implicit none 
 
contains

!-------------------------------------------

subroutine make_data( Data_mat, m, p, q, priortype, tData, labels, sample )

!This subroutine transformes the data to multivariate regression form:
!	  Y   =  Z  *  B  +  U
!	(Txm)  (Txk) (kxm) (Txm)
!k = mp + q; m = no of modelled var's, p = lag length, q = no of unmodelled (determ.) var's
!

! DATA_mat	= data matrix
! m			= number of dependent variables
! p			= number of lags of dep vars
! q			= number of deterministic variables
! tData		= derived type for data

double precision, intent(in) :: Data_mat(:,:)
integer, intent(in) :: m, p, q
character*(*), intent(in) :: priortype
type(VARData), intent(inout) :: tData
character(len=*), intent(in), optional :: labels(:)
type(SampleInfo), intent(in), optional :: sample

integer :: Nobs, T, i, j, k, n, w, jj, fobs, lobs

Nobs = size( Data_mat, 1 );
fobs = p+1
lobs = Nobs

if ( present(sample) ) then
    if ( sample%first > 0 ) then
        if ( sample%first > p ) then
            fobs = sample%first
        else
            print *, "First obs", sample%first, " not possible with", p, " lags"
            stop
        endif
    endif
    if ( sample%last > 0 ) then
        if ( sample%last <= Nobs ) then
            lobs = sample%last
        else
            print *, "Last obs ", sample%last, " not possible with", Nobs, " observations"
            stop
        endif
    endif
endif
T = lobs-fobs+1

call clear_data( tData )
allocate( tData%Y(T,m), tData%Z(T,m*p+q) )

if ( present( labels ) ) then
    if ( size(labels) /= m+q ) stop "Dimension of variable labels different from number of variables"
    allocate( tData%Ylabels(m), tData%Zlabels(m*p+q) )
    tData%Ylabels = labels(1:m)
    if ( q > 0 ) tData%Zlabels(1:q) = labels(m+1:m+q)
    k = q+1
    do i = 1, p
        do j = 1, m
            tData%Zlabels(k) = labels(j)
            w = len_trim(labels(j)) + 1
            n = int(log10(dble(i))) + 1
            write(tData%Zlabels(k)(w:),101) i
            k = k+1
        end do
    end do
endif  

tData%Y(:,1:m) = DATA_mat(fobs:lobs,1:m)
tData%Z(:,1:q) = DATA_mat(fobs:lobs,m+1:m+q)
do i = 1,p
	tData%Z(1:T,q+(i-1)*m+1:q+i*m) = DATA_mat(fobs-i:lobs-i,1:m)
enddo	

if ( priortype == 'SS' ) then

    allocate( tData%D(T,q), tData%W(T,p*q), tData%B(T,(p+1)*q) )

    tData%D = DATA_mat(fobs:lobs,m+1:m+q)
    do i = 1,p
	    tData%W(1:T,(i-1)*q+1:i*q) = DATA_mat(fobs-i:lobs-i,m+1:m+1)
    enddo	
    tData%B(:,1:q) = tData%D
    tData%B(:,q+1:) = -tData%W

endif

tData%firstobs = fobs
tData%lastobs  = lobs
tData%Nobs = T
tData%m = m
tData%p = p
tData%q = q
tData%k = m*p+q

101 format( '_', I<n> )

end subroutine

!-------------------------------------------

subroutine clear_data( tData )

type(VARdata), intent(inout) :: tData

if ( allocated( tData%Y ) )            deallocate( tData%Y )
if ( allocated( tData%Z ) )            deallocate( tData%Z )
if ( allocated( tData%D ) )            deallocate( tData%D )
if ( allocated( tData%W ) )            deallocate( tData%W )
if ( allocated( tData%B ) )            deallocate( tData%B )
if ( allocated( tData%ZPZ ) )          deallocate( tData%ZPZ )
if ( allocated( tData%ZPY ) )          deallocate( tData%ZPY )
if ( allocated( tData%Gamma_OLS ) )    deallocate( tData%Gamma_OLS )
if ( allocated( tData%vecGamma_OLS ) ) deallocate( tData%vecGamma_OLS )
if ( allocated( tData%ESS ) )          deallocate( tData%ESS )
if ( allocated( tData%Ylabels ) )      deallocate( tData%Ylabels )
if ( allocated( tData%Zlabels ) )      deallocate( tData%Zlabels )

end subroutine

!-------------------------------------------

subroutine get_DataMat( tDatain, datatype, DataMat, firstobs, lastobs )

type(VARdata), intent(in) :: tDatain
character(len=*), intent(in) :: datatype
double precision, intent(out) :: DataMat(:,:)
integer, optional, intent(in) :: firstobs, lastobs ! relative to original data set before loosing stuff to lags

integer fobs, lobs

if ( present( firstobs ) ) then
    if ( firstobs < tDatain%firstobs ) then
        print *, "First obs for get_DataMat outside of data range"
        stop
    endif
    fobs = firstobs - tDatain%firstobs + 1
else
    fobs = 1
endif
if ( present( lastobs ) ) then
    if ( lastobs > tDatain%lastobs ) then
        print *, "Last obs for get_DataMat outside of data range"
        stop
    endif
    lobs = lastobs - tDatain%firstobs + 1
else
    lobs = tDatain%lastobs - tDatain%firstobs + 1
endif

!if ( size(DataMat,1) /= lobs-fobs+1 ) then
!    print *, "Row dimension in data matrix is wrong for get_DataMat"
!    stop
!endif

select case(datatype)
case('Y')
    ! dependendt variables
    
    if ( size(DataMat,2) /= size(tDatain%Y,2) ) then
        print *, "Column dimension in data matrix wrong for get_DataMat"
        stop
    endif

    DataMat = tDatain%Y(fobs:lobs,:)

case('X')
    ! deterministic/exogenous variables
    
    if ( tDatain%q <= 0 ) then
        print *, "No deterministic or exogeneous variables in model, get_DataMat"
        stop
    endif
    if ( size(DataMat,2) /= tDatain%q ) then
        print *, "Column dimension in data matrix wrong for get_DataMat"
        stop
    endif
    
    DataMat = tDatain%Z(fobs:lobs,1:tDatain%q)
    
case('Z')
    ! lags of Y

    if ( tDatain%p <= 0 ) then
        print *, "No lagged variables in model, get_DataMat"
        stop
    endif
    if ( size(DataMat,2) /= tDatain%k-tDatain%q ) then
        print *, "Column dimension in data matrix wrong for get_DataMat"
        stop
    endif
    
    DataMat = tDatain%Z(fobs:lobs,tDatain%q+1:tDatain%k)
    
case default
    print *, "Unrecognized datatype for get_DataMat"
    stop
end select

end subroutine

!-------------------------------------------

subroutine copy_data( tDatain, tDataout, firstobs, lastobs, sample )

type(VARdata), intent(in) :: tDatain
type(VARdata), intent(inout) :: tDataout
integer, optional, intent(in) :: firstobs, lastobs ! relative to original data set before loosing stuff to lags
type(SampleInfo), optional, intent(in) :: sample
!-----

integer fobs, lobs

call clear_data( tDataout )

fobs = tDatain%firstobs
lobs = tDatain%lastobs

if ( present( sample ) ) then
    if ( sample%first > 0 ) fobs = sample%first
    if ( sample%last > 0 ) lobs = sample%last
else
    if ( present( firstobs ) ) fobs = firstobs
    if ( present( lastobs ) ) lobs = lastobs
endif
    
if ( fobs < tDatain%firstobs ) then
    print *, "First obs for data copy outside of data range"
    stop
endif
fobs = fobs - tDatain%firstobs + 1
if ( lobs > tDatain%lastobs ) then
    print *, "Last obs for data copy outside of data range"
    stop
endif
lobs = lobs - tDatain%firstobs + 1

tDataout%Nobs = lobs-fobs+1 
tDataout%m = tDatain%m 
tDataout%p = tDatain%p 
tDataout%q = tDatain%q
tDataout%k = tDatain%k
tDataout%firstobs = tDatain%firstobs + fobs - 1
tDataout%lastobs  = tDatain%firstobs + lobs - 1

allocate( tDataout%Y(tDataout%Nobs,tDataout%m), tDataout%Z(tDataout%Nobs,tDataout%k) )

tDataout%Y = tDatain%Y(fobs:lobs,:)
tDataout%Z = tDatain%Z(fobs:lobs,:)

if ( allocated( tDatain%D ) ) then

    allocate( tDataout%D(tDataout%Nobs,tDataout%q), &
              tDataout%W(tDataout%Nobs,tDataout%p*tDataout%q), &
              tDataout%B(tDataout%Nobs,(tDataout%p+1)*tDataout%q) )

    tDataout%D = tDatain%D(fobs:lobs,:)
    tDataout%W = tDatain%W(fobs:lobs,:)
    tDataout%B = tDatain%B(fobs:lobs,:)

endif

if ( allocated( tDatain%Ylabels ) ) then
    allocate( tDataout%Ylabels(size(tDatain%Ylabels)) )
    tDataout%Ylabels = tDatain%Ylabels
endif
if ( allocated( tDatain%Zlabels ) ) then
    allocate( tDataout%Zlabels(size(tDatain%Zlabels)) )
    tDataout%Zlabels = tDatain%Zlabels
endif
              
end subroutine

!-------------------------------------------

subroutine OLS( tData )

use matutil

type(VARData), intent(inout) :: tData

integer T, m, p, k, i
double precision, allocatable :: chol(:,:), err(:,:)

T = tData%Nobs
m = tData%m
p = tData%p
k = tData%k

if ( allocated( tData%ZPZ ) )          deallocate( tData%ZPZ )
if ( allocated( tData%ZPY ) )          deallocate( tData%ZPY )
if ( allocated( tData%Gamma_OLS ) )    deallocate( tData%Gamma_OLS )
if ( allocated( tData%vecGamma_OLS ) ) deallocate( tData%vecGamma_OLS )
if ( allocated( tData%ESS ) )          deallocate( tData%ESS )

allocate( tData%ZPZ(k,k), tData%ZPY(k,m), &
		  tData%Gamma_OLS(k,m), tData%vecGamma_OLS(k*m), &
		  tData%ESS(m,m), chol(k,k), err(T,m) )

CALL DMXTXF( T, k, tData%Z, T, k, tData%ZPZ, k )					!Z'Z
CALL DMXTYF( T, k, tData%Z, T, T, m, tData%Y, T, k, m, tData%ZPY, k ) 	!Z'Y

CALL DLFTDS( k, tData%ZPZ, k, chol, k )
do i = 1, m
	CALL DLFSDS( k, chol, k, tData%ZPY(:,i), tData%Gamma_OLS(:,i) )
end do
call vectorize( tData%Gamma_OLS, tData%vecGamma_OLS )

CALL DMRRRR( T, k, tData%Z, T, k, m, tData%Gamma_OLS, k, T, m, err, T )		!ZB_ols
err = tData%Y - err						!Residuals
CALL DMXTXF( T, m, err, T, m, tData%ESS, m )	! Error sum of squares

deallocate( chol, err )

end subroutine

!*******************************************************

subroutine make_predcov( MA, ecovroot, cov )

! make variance-covariance matrix for prediction errors

double precision, intent(in) :: MA(:,:,:)
double precision, intent(in) :: ecovroot(:,:)
double precision, intent(out) :: cov(:,:)

double precision, allocatable :: MAt(:,:)
integer m, i, j, ii, jj, k, kk

if ( size(MA,3) == 0 ) then
	! no MA-poly, just lead time 1
	cov = matmul(transpose(ecovroot),ecovroot)
	return
endif

allocate( MAt(size(MA,1),size(MA,2)) )

m = size( ecovroot, 1 )
cov = 0.0d0
cov(1:m,1:m) = ecovroot
do i = 1, size(MA,3)
	MAt(:,:) = matmul(MA(:,:,i),ecovroot)
	ii = 1; jj = m
	k = i*m+1; kk = (i+1)*m
	cov(k:kk,k:kk) = ecovroot
	do j = i, size(MA,3)
		cov(ii:jj,k:kk) = MAt
		ii = ii+m; jj = jj+m; k = k+m; kk = kk+m
	end do
end do

cov = matmul(transpose(cov),cov)

deallocate( MAt )

end subroutine

!*******************************************************

subroutine get_AR( m, p, B, AR )

! AR-parameters are assumed to be in last m*p
! rows of B

integer, intent(in) :: m, p
double precision, intent(in) :: B(:,:)
double precision, intent(out) :: AR(:,:,:)

integer i, k, n
if ( size(B,1) < m*p ) then
	print *, "Row dimension of B is wrong"
	stop
endif
if ( size(B,2) /= m ) then
	print *, "Column dimension of B is wrong"
	stop
endif
if ( size(AR,1) /= m ) then
	print *, "Row dimension of AR is wrong"
	stop
endif
if ( size(AR,2) /= m ) then
	print *, "Column dimension of AR is wrong"
	stop
endif
if ( size(AR,3) /= p ) then
	print *, "Lag dimension of AR is wrong"
	stop
endif

k = size(B,1) - m*p + 1
n = k + m - 1
do i = 1, p
	AR(:,:,i) = B(k:n,:)
	k = k+m
	n = n+m
end do

end subroutine

!********************************************************

subroutine MA_polynomial( AR, MA )

! Solve for MA-polynomial given AR-polynomial

double precision, intent(in) :: AR(:,:,:) ! AR(0) = I is implied
double precision, intent(out) :: MA(:,:,:) ! MA(0) = I is implied


integer m, p, q, i, j, n

m = size( AR, 1 )
i = size( AR, 2 )
p = size( AR, 3 )
if ( i /= m ) then
	print *, "Dimension error, AR-polynomial is not square", m, i
	stop
endif
if ( size( MA, 1 ) /= m ) then
	print *, "Dimension error, rows different in AR- and MA-polynomial"
	stop
endif
i = size( MA, 2 )
if ( i /= m ) then
	print *, "Dimension error, MA-polynomial is not square", m, i
	stop
endif
q = size( MA, 3 )

if ( q == 0 ) then
	! no MA-poly requested
	return
endif

MA = 0.0d0

MA(:,:,1) = AR(:,:,1)
do i = 2, q
	if ( i <= p ) then	
		MA(:,:,i) = AR(:,:,i)
	endif
	n = min( p, i-1 )
	do j = 1, n
		MA(:,:,i) = MA(:,:,i) + matmul( AR(:,:,j), MA(:,:,i-j) )
	end do
end do

end subroutine

!********************************************************

subroutine VARforecast(Z,B,forecast,X,VC,VCROOT,VCDIAG)

!This routine produces forecasts for a vector autoregression
!input:		Z_row = row vector [X1t ... Xqt, Y1t-1,...,Ymt-1, Y1t-2,...Ymt-p]
!			B = coefficient matrix
!			Y_forc = m-vector of forcasts = [Y1t,...,Ymt] 
!			Y_forc = Z_rowB
!
!output:				|Y1T+1,...,YmT+1|	
!			forc_mat =	| *           * |
!						|Y1T+h,...,YmT+h|
!
!Written by M.K. Andersson, Sveriges Riksbank, 2004-09-22

! Z			= vector of explanatory (non-deterministic) variables for forecast origin
! B			= matrix of regression coefficients
! forecast	= H by m matrix of forecasts 1 to H periods ahead (output)
! X			= H by q matrix of deterministic variables for origin to origin+H-1 (optional)
! VC		= error variance covariance (optional)
! VCROOT	= Cholesky factor of error v-cov (optional)
! VCDIAG	= vector containing diagonal error v-cov (optional)
! 
! Deterministic forecasts are produced if none of VC, VCROOT or VCDIAG are present
! If more than one is present VCDIAG is used first, then VCROOT

use imsl_libraries

implicit none
double precision, intent(in) :: Z(:)
double precision, intent(in) :: B(:,:)
double precision, intent(out) :: forecast(:,:)
! optional arguments
double precision, intent(in), optional :: X(:,:)  ! deterministic variables in model
double precision, intent(in), optional :: VC(:,:), VCROOT(:,:), VCDIAG(:)


integer p, m, k, h, q, i,irank
double precision, allocatable :: VC_chol(:,:), ERR_M(:,:), y_forc(:), tmp2(:), &
                                 Z_temp(:), forc_mat(:,:)


k = size(B,1)
m = size(B,2)
h = size(forecast,1)

if ( present(X) ) then
  q = size(x,2)
  if ( size(x,1) < h ) then
    stop 'To few rows in X for forecasts'
  endif
else
  q = 0
endif
if ( size(Z) /= k-q ) then
  print *, 'Dimension of Z wrong'
  stop
endif
p = size(z)/m
if ( size(z) /= p*m ) then
  print *, 'Dimension of Z not multiple of VAR dimension'
  stop
endif
if ( size(forecast,2) > m ) then
  print *, 'To many columns in forecast'
  stop
endif

allocate( ERR_M(h,m), y_forc(m), tmp2(m), Z_temp(k), forc_mat(h,m) )

! draw residuals if variance is specified otherwise deterministic forecast
ERR_m = 0.0d0
if ( present(VCDIAG) ) then
  if ( size(VCDIAG) /= m ) then
    stop 'VCDIAG has wrong dimension'
  endif
  y_forc = sqrt(VCDIAG)
  do i = 1, h
    call DRNNOA( m, tmp2 )
    ERR_M(i,:) = tmp2*y_forc
  enddo
elseif ( present(VCROOT) ) then
  if ( size(VCROOT,1) /= m .or. size(VCROOT,2) /= m ) then
    stop 'VCROOT has wrong dimension'
  endif
  CALL DRNMVN(h,m,VCroot,m,ERR_M,h)						!ERR ~ MVN(0,VCOV)
elseif ( present(VC) ) then
  if ( size(VC,1) /= m .or. size(VC,2) /= m ) then
    stop 'VC has wrong dimension'
  endif
  allocate( VC_chol(m,m) )
  CALL DCHFAC(m,VC,m,100*epsilon(1.0d0),IRANK,VC_chol,m)		!Choldec of VCOV
  CALL DRNMVN(h,m,VC_chol,m,ERR_M,h)						!ERR ~ MVN(0,VCOV)
  deallocate( VC_chol )
endif
 
y_forc = 0.0d0
forc_mat = 0.0d0

! Make dynamic forecasts from 1 step to H step
Z_temp(q+1:) = Z
do i = 1,H
	if ( q > 0 ) then
  		Z_temp(1:q) = X(i,:)
	endif

	y_forc = matmul(Z_temp,B) + ERR_M(i,:)
	forc_mat(i,:) = y_forc

	Z_temp(q+m+1:k) = Z_temp(q+1:q+(p-1)*m)
	Z_temp(q+1:q+m) = y_forc
	
enddo	

forecast = forc_mat(1:H,1:size(forecast,2))

deallocate( ERR_M, y_forc, tmp2, Z_temp, forc_mat )

end subroutine


end module
