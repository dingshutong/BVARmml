module datautils

INCLUDE 'link_fnl_static.h'  !DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries
USE STRUCTURES

implicit none 
    
contains

subroutine predsetup(H,q,predfile,var_name,X_pred,NAME_FLAG)

!This routine reads the forecast period values for the deterministic variables

! H = number of observations ahead to forecast, number of obs to read
! q = number of deterministic variables
! predfile = file with data on deterministic variables
! X_pred   = matrix of observations on exog vars for forecast period (output)
! NAME_FLAG = 1 if check that variable names in file are same as in model (output)

!Written by M.K. Andersson, Sveriges Riksbank, 2004-09-22

implicit none
integer, intent(in) :: H,q
character*12, intent(in) :: predfile
character*12,dimension(q), intent(in) :: var_name
double precision, dimension(H,q), intent(out) :: X_pred
integer, intent(inout) :: NAME_flag

integer :: i,j
character*12 :: pred_name

OPEN(1000,FILE=predfile)

READ(1000,*) PRED_NAME					!Read variable names

!Do i = 1,q
	IF (var_name(1) .ne. pred_name) THEN
		write(*,*) 'The deterministc variables supplied for the forecast'
		write(*,*) 'are not the ones that are declared as deterministic'
		write(*,*) 'in the main data set'
		NAME_FLAG = 1
	ENDIF
!ENDDO

do i = 1,H
	READ(1000,*) (X_pred(i,j),j = 1,q)			!Read pred data for X)
enddo

close(1000) 

end subroutine

!-------------------------------------------

subroutine GENDATA(N,H,mc_lag,m_syst,m_spur,q,param_syst,param_spur,var_err,DATAMAT,fcast_vec,param_ex)
use imsl_libraries

! Generate data for Monte Carlo experiment

! N			 = number of observations  
! H			 = Number of observations ahead to forecast
! mc_lag	 = lags in data generating process
! m_syst	 = number of variables in data generating process
! m_spur	 = number of additional variables to generate
! q			 = number of deterministic variables
! param_syst = parameters on lags of endogenous vars in DGP
! param_spur = parameters on lags of additional vars in DGP for additional vars
! param_ex	 = parameters on exognous variables, true DGP first, additional vars next
! var_err    = vector of error variances
! datamat	 = matrix of generated data, endogenous variables first (output)
! fcast_vec	 = observations for the forecast evaluation period (output)

implicit none
integer, intent(in) :: N,H,mc_lag,m_syst,m_spur,q
double precision, dimension(m_syst,m_syst,mc_lag), intent(in) :: param_syst
double precision, dimension(m_spur+m_syst,m_spur,mc_lag), intent(in) :: param_spur
double precision, intent(in), optional :: param_ex(q,m_syst+m_spur)
double precision, dimension(m_syst+m_spur), intent(in) :: var_err
double precision, dimension(N,m_syst+m_spur+q) :: DATAMAT
double precision, dimension(-10:H), intent(out) :: fcast_vec


integer :: no_var,i, lag
double precision, dimension(m_syst+m_spur) :: Err_vec
double precision, dimension(N+100+H,m_syst+m_spur) :: DATAMAT_temp
double precision, dimension(m_syst+m_spur) :: std_err


no_var = m_syst+m_spur


std_err = sqrt(var_err)

if (q>0 .and. .not.(PRESENT(param_ex))) then 
write(*,*) 'parameters on exognous variables must be presented since q>0' 
	pause
endif 


DATAMAT_temp = 0.0d0
do i = mc_lag+1,N+100+H
	
	! Generate 'true' system
	do lag = 1, mc_lag
		datamat_temp(i,1:m_syst) = datamat_temp(i,1:m_syst) + matmul( datamat_temp(i-lag,1:m_syst), param_syst(:,:,lag) )
	enddo
	if ( q == 1 ) then !constant
		datamat_temp(i,1:m_syst) = datamat_temp(i,1:m_syst) + param_ex(1,1:m_syst)
	else if ( q == 2 ) then   !time trend
		datamat_temp(i,1:m_syst) = datamat_temp(i,1:m_syst) + param_ex(1,1:m_syst) + param_ex(2,1:m_syst)*(i-100)
	end if
	call drnnoa(no_var,err_vec)
	datamat_temp(i,1:m_syst) = datamat_temp(i,1:m_syst) + err_vec(1:m_syst)*std_err(1:m_syst)

	! generate spurios
	do lag = 1, mc_lag
		datamat_temp(i,m_syst+1:) = datamat_temp(i,m_syst+1:) + matmul( datamat_temp(i-lag,:), param_spur(:,:,lag) )
	enddo
	if ( q == 1 ) then
		datamat_temp(i,m_syst+1:) = datamat_temp(i,m_syst+1:) + param_ex(1,m_syst+1:)
	else if ( q == 2 ) then
		datamat_temp(i,m_syst+1:) = datamat_temp(i,m_syst+1:) + param_ex(1,m_syst+1:) + param_ex(2,m_syst+1:)*(i-100)
	end if
	datamat_temp(i,m_syst+1:) = datamat_temp(i,m_syst+1:) + err_vec(m_syst+1:)*std_err(m_syst+1:)

enddo

datamat(1:N,1:no_var) = datamat_temp(101:100+N,:)
if ( q >= 1 ) then
	datamat(:,no_var+1) = 1
endif
if ( q == 2 ) then
	datamat(:,no_var+2) = (/(i,i=1,n)/)
endif
fcast_vec = datamat_temp(100+N-10:100+N+H,1)


open(69,file='gendata1.dat')
do i = 1,N
	write(69,217) datamat(i,:)
enddo
217 format(8(F8.4, tr1))

end subroutine


!===================================================================================
subroutine READ_DATA(matrix,filename,NR,bin)
INTEGER :: STATUS, i
INTEGER, OPTIONAL :: bin, NR
character(len=*), intent(in) :: filename
double precision, intent(out) :: matrix(:,:)
character(len=len_trim(filename)+14) inputfile

IF (PRESENT(bin)) THEN
inputfile = trim(adjustl(filename))//'.bindata'
open( 10, file=inputfile, action='read' , form='binary', iostat=status )
  if ( status /= 0 ) then
    print *, "Failed opening data file '"//trim(adjustl(inputfile))//"', status code:", status
    stop
  endif
     IF (PRESENT(NR)) THEN
       READ(10) (matrix(i,:),i=1,NR)
     ELSE 
      READ(10) matrix
     ENDIF
close( 10 )
ELSE
inputfile = trim(adjustl(filename))//'.txt'
open( 10, file=inputfile, action='read', iostat=status )
  if ( status /= 0 ) then
    print *, "Failed opening data file '"//trim(adjustl(inputfile))//"', status code:", status
    stop
  endif
   IF (PRESENT(NR)) THEN
    READ(10,*) (matrix(i,:),i=1,NR)
   ELSE
     READ(10,*) matrix
   ENDIF
close( 10 )
ENDIF
end subroutine

!=====================================================================================

subroutine WRITE_DATA(matrix, filename, bin)
INTEGER :: NR,NC,STAT, i,j
INTEGER, OPTIONAL :: bin
character(len=*), intent(in) :: filename
double precision, intent(in) :: matrix(:,:)
character(len=len_trim(filename)+14) outputfile
NR=SIZE(matrix,1)
NC=SIZE(matrix,2)
IF (PRESENT(bin)) THEN
outputfile = trim(adjustl(filename))//'.bindata'
open( 10, file=outputfile, action='write', status='replace', form='binary', iostat=stat)
if ( stat /= 0 ) then
    print *, "Failed opening data file '"//trim(adjustl(outputfile))//"', status code:", stat
    stop
endif
WRITE(10) matrix
close(10)

ELSE
outputfile = trim(adjustl(filename))//'.txt'
open( 10, file=outputfile, action='write', status='replace', iostat=stat )
if ( stat /= 0 ) then
    print *, "Failed opening data file '"//trim(adjustl(outputfile))//"', status code:", stat
    stop
endif
WRITE(10,*) matrix
close(10)
ENDIF
end subroutine

!---------------------------------------------------------------

subroutine dump_data( matrix, filename, labels, mcmc )

double precision, intent(in) :: matrix(:,:)
character(len=*), intent(in) :: filename
character(len=*), optional, intent(in) :: labels(:)
type(MCMCInfo), optional, intent(in) :: mcmc

integer status
character(len=len_trim(filename)+14) outfile

outfile = trim(adjustl(filename))//'.bindata'
open( 40, file=outfile, action='write', status='replace', form='binary', iostat=status )
if ( status /= 0 ) then
    print *, "Failed opening data file '"//trim(adjustl(outfile))//"', status code:", status
    stop
endif

write(40) matrix

close( 40 )

outfile = trim(adjustl(filename))//'.datadesc.txt'
open( 40, file=outfile, action='write', status='replace', iostat=status )
if ( status /= 0 ) then
    print *, "Failed opening data description file '"//trim(adjustl(outfile))//"', status code:", status
    stop
endif

write(40,101) size(matrix,1), size(matrix,2)
if ( present(mcmc) ) then
    write(40,102) "MCMC", mcmc%reps, mcmc%thin
endif
if ( present(labels) ) then
    write(40,103) "LABELS", labels
endif

close(40)

101 format (tr1,i,tr1,i)
102 format (tr1,a,tr1,i,tr1,i)
103 format (8(tr1,a))

end subroutine
!------------------------------------------------------------
subroutine WRITE_TAB(TAB, filename, bin)
INTEGER :: NR,NC,STAT, i,j
INTEGER, OPTIONAL :: bin
character(len=*), intent(in) :: filename
INTEGER, intent(in) :: TAB(:,:)
character(len=len_trim(filename)+14) outputfile
NR=SIZE(TAB,1)
NC=SIZE(TAB,2)
IF (PRESENT(bin)) THEN
outputfile = trim(adjustl(filename))//'.bindata'
open( 10, file=outputfile, action='write', status='replace', form='binary', iostat=stat)
if ( stat /= 0 ) then
    print *, "Failed opening data file '"//trim(adjustl(outputfile))//"', status code:", stat
    stop
endif
WRITE(10) TAB
close(10)

ELSE
outputfile = trim(adjustl(filename))//'.txt'
open( 10, file=outputfile, action='write', status='replace', iostat=stat )
if ( stat /= 0 ) then
    print *, "Failed opening data file '"//trim(adjustl(outputfile))//"', status code:", stat
    stop
endif
WRITE(10,*) TAB
close(10)
ENDIF
end subroutine

end module
