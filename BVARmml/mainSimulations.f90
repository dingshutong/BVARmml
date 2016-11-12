program mainSimulations

include 'link_fnl_static.h'!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries

use postProc
use times
use simutils
use procUtils


implicit none


! structures
type(DataArray) rawdata
type(ModelInfo) model
type(MCMCInfo) mcmc
type(PriorInfo) prior
type(SampleInfo) sample
type(StartVal) Startvalues
type(VARData) tdata
type(dgp)  :: tdgp
type(mods) :: mod
type(rmses)   :: rmse,rmse1
type(VARPrior) :: tPrior

!Characters

character*100 :: text
character*20  :: outfile,endcheck,method
character*8   :: drivers
character*12  :: predfile
character*2   :: priortype(2) 
character*20  :: outfilename
character*40 :: titlename

! integer

integer :: T,i,j,ms,series_main,ISEED,priordf,thread_rank,thread_i,SEEDIN,thread_t,iRANK,mod_endo,num_modprior,name_flag,result_i,mc_p

! real*8

real*8 :: modelprior,numpercent,lmarglik


real*8, allocatable :: vec_syst(:),vec_spur(:),mpriort(:),prior_val(:,:),modprior(:),datamat(:,:,:),fcast_vec(:,:),X_pred(:,:)

! simresults
real*8, allocatable :: varprob(:,:,:),dimprob(:,:,:),modprob(:,:),modselect(:,:)

! True output

real*8, allocatable :: vartable(:,:),dimtable(:,:),modtable(:,:),rmsetable(:,:),datastore(:,:)


real*8, allocatable ::  Gamma_draw(:,:), Psi_draw(:,:), psi_drawroot(:,:),s_chol(:,:),s_post(:,:,:), &
           gamma_draws(:,:),psi_draws(:,:,:),psiroot_draws(:,:,:)

integer mc_burn

!time
INTEGER dtime0(8), dtime1(8), dtime2(8), date_time(8)
real    :: time1, time2,timediff
CHARACTER(len=12) :: Real_clock(3)


!logic

logical EXISTS


!**********************************
!******* Read inputs **************
!**********************************

CALL CPU_TIME(TIME1)
call DATE_AND_TIME(values=DTIME1)


101 format(A100)

PRINT *, "Which data generating process (DGP1/DGP2/DGP3)"
READ *, tdgp%dgps
! tdgp%dgps="DGP1"

PRINT *, "The sample size of DGP (100/300) "
READ *, tdgp%N
!tdgp%N =100


drivers= trim(adjustl(tdgp%dgps//".in"))

call DATE_AND_TIME(values=DTIME1)

OPEN(1,FILE=drivers)

!==========================DGP.IN ========================

!N=serial length, m=number of modelled variables
!q=numb. of non-modelled series, p=laglength
			
read(1,101) text
READ(1,*) tdgp%m	
	
read(1,101) text
READ(1,*) tdgp%q
			
read(1,101) text
READ(1,*) tdgp%p

read(1,101) text
READ(1,*) tdgp%simreps

read(1,101) text
READ(1,*) tdgp%mcreps

read(1,101) text
read(1,*) seedin

read(1,101) text
read(1,*) mod%q 

allocate(tdgp%ser_name(tdgp%m+mod%q))

read(1,101) text
read(1,*) tdgp%ser_name 
!DIMS(1) = minimum dimension of VAR
!DIMS(2) = maximum dimension of VAR		

read(1,101) text
READ(1,*) tdgp%DIMS(1),tdgp%DIMS(2)	

!Number of presented top models
							
read(1,101) text
read(1,*) tdgp%present_top	

read(1,101) text
read(1,*) mc_p            !mc_lags is the number of model lags

 

read(1,101) text
Read(1,*) priordf	        !The prior degrees of freedom

read(1,101) text
READ(1,*) i
allocate( mod%evalvec(i))

read(1,101) text
read(1,*) mod%evalvec

read(1,101) text
read(1,*) numpercent 

mod%NUM = CEILING(DBLE(numpercent*tdgp%N))

read(1,101) text
read(1,'(A)') mod%posteval

read(1,101) text
read(1,*) series_main

read(1,101) text
READ(1,*) tdgp%H

read(1,101) text
READ(1,*) tdgp%m_syst

read(1,101) text
READ(1,*) tdgp%m_spur

allocate(vec_syst(tdgp%m_syst*tdgp%m_syst*tdgp%p),vec_spur(tdgp%m*tdgp%m_spur*tdgp%p),tdgp%var_err(tdgp%m_syst+tdgp%m_spur))

read(1,101) text
read(1,*) vec_syst

!read(1,101) text
!read(1,*) vec_spur

read(1,101) text
read(1,*) tdgp%var_err

read(1,101) text
read(1,*) num_modprior

allocate(modprior(num_modprior))

read(1,101) text
read(1,*) modprior

read(1,101) text
read(1,*) outfile

read(1,101) text
read(1,*) method
		
read(1,101) text
READ(1,*) ENDCHECK		!End check
close(1)



if (mod%q .gt. 2) then
	write(*,*) 'only 2 deterministc components allowed, (icept & trend)' 
	pause
endif




!---------- Initialized ----------!


tdgp%MCBURN=CEILING(DBLE(tdgp%MCREPS/10))


! Litterman prior 
prior%priordf=priordf
allocate(prior%Litterman%tightness(5),prior%Litterman%firstlagmean(2))
prior%Litterman%tightness(1)=0.5d0
prior%Litterman%tightness(2)=0.5d0
prior%Litterman%tightness(3)=1d0
prior%Litterman%tightness(4)=0.5d0
prior%Litterman%tightness(5)=5d0
prior%Litterman%firstlagmean(1)=0.9d0
prior%Litterman%firstlagmean(2)=1d0
if (mod%q>0) then
  allocate(prior%Litterman%deterministicmean(mod%q,tdgp%m))
  prior%Litterman%deterministicmean=0d0
end if 


T=tdgp%N-mod%p

ms=series_main

allocate(tdgp%param_syst(tdgp%m_syst,tdgp%m_syst,tdgp%p),tdgp%param_spur(tdgp%m,tdgp%m_spur,tdgp%p))
call antivec3(vec_syst,tdgp%m_syst,tdgp%m_syst,tdgp%p,tdgp%param_syst)
do i=1,tdgp%p
  tdgp%param_spur(:,:,i)=DZERO(tdgp%m,tdgp%m_spur)
end do

tdgp%param_spur(1,1,1)=0.3d0
tdgp%param_spur(3,1,1)=0.5d0
tdgp%param_spur(2,2,1)=0.5d0 
tdgp%param_spur(4,2,1)=0.5d0
tdgp%param_spur(5,3,1)=0.2d0
tdgp%param_spur(6,4,1)=0.7d0


!call antivec3(vec_spur,tdgp%m,tdgp%m_spur,tdgp%p,tdgp%param_spur)

allocate(mod%bma_comb(tdgp%dims(2)))
mod%bma_comb = 0
do i = 1,tdgp%m_syst
	mod%bma_comb(i) = i
enddo




!---------- Model combinations ----------!
call mod_num(tdgp%m,ms,tdgp%dims,mod%n,mod%no_mods_del)
allocate(mod%models(mod%n,tdgp%dims(2)+1),mod%lnprior(mod%n,num_modprior),prior_val(tdgp%m,num_modprior))
call mod_search(tdgp%m,ms,mod%n,tdgp%dims(2),mod%models)


do i=1,num_modprior
   prior_val(:,i)=modprior(i)
end do


do i = 1,mod%n
     mod_endo=mod%models(i,tdgp%dims(2)+1)
     
     do j = 1,num_modprior
           call mod_prior(tdgp%m,prior_val(:,j),mod_endo,mod%models(i,:),mod%lnprior(i,j))
     end do
enddo

do i = 1,mod%n  
      
	 if (all(mod%bma_comb(1:tdgp%dims(2)) == mod%models(i,1:tdgp%dims(2)))) then
		mod%true = i
		exit
	endif
	
enddo  


if (tdgp%present_top > mod%n) then
    tdgp%present_top=mod%n
end if 


allocate(mod%rmse_obs(tdgp%h))
do i=1,tdgp%h
  mod%rmse_obs(i)=tdgp%simreps
end do 




allocate(datamat(tdgp%N,tdgp%m+mod%q,tdgp%simreps),datastore(tdgp%N*(tdgp%m+mod%q),tdgp%simreps),fcast_vec(-10:tdgp%H,tdgp%simreps))


iseed = seedin

call RNSET(iseed)


  
   if ((mod%q .ne. 0) .and. (tdgp%H .ne. 0)) then
       allocate(X_pred(tdgp%H+mod%num+maxval(mod%evalvec),mod%q))
       !call predsetup(tdgp%H+mod%num+maxval(mod%evalvec),mod%q,predfile,tdgp%ser_name(tdgp%m+1:tdgp%m+mod%q),X_pred,name_flag)
       !if (name_flag ==1) then
       !   goto 1000
       !end if 
       X_pred=1d0
   end if
   
!open(unit=111,form ='unformatted',status ='scratch')

    !---------- Generalized data ----------!
    
open(unit=99,file='seedinfo.out',status='unknown')

do i = 1, tdgp%simreps
    call rnget(iseed)
   write(99,*) 'monte carlo loop', i, 'seed', iseed 
   
  CALL GENDATA(tdgp%N,tdgp%H,tdgp%p,tdgp%m_syst,tdgp%m_spur,tdgp%q,tdgp%param_syst,tdgp%param_spur,tdgp%var_err,DATAMAT(:,1:tdgp%m,i),fcast_vec(:,i))
  ! if you want to have deterministic variables which is stored in param_ex
	!call GENDATA(N,H,p,m_syst,m_spur,q,param_syst,param_spur,param_ex,var_err,DATAMAT,fcast_vec)
     ! Add a constant to the model
    if (mod%q .eq. 1) then 
       do j = 1,tdgp%N  
         datamat(j,tdgp%m+1,i) = 1d0 
       end do  
    end if  
   datastore(:,i)=  vector(DATAMAT(:,:,i))          
end do

call WRITE_DATA(datastore, "stored_data")
call write_data(fcast_vec,"fcast_scored")
    
    
allocate(varprob(tdgp%m,num_modprior,2),dimprob(tdgp%dims(2),num_modprior,2),modprob(num_modprior,2),modselect(num_modprior,2))


!PRINT *, "Prior Type (NW/ND)"
!READ *, prior%priortype

!prior%priortype="NW"



allocate(modtable(num_modprior*2,4),rmsetable(tdgp%H,8))

allocate(vartable(num_modprior*2,(tdgp%m_syst-1)*4+2),dimtable(num_modprior*2,tdgp%dims(2)*2))

priortype(1)="NW"
priortype(2)="ND"



do result_i = 2,2


prior%priortype=priortype(result_i)

if (prior%priortype .EQ. "NW") THEN

Write(*,*)  '  '
write(*,*)  '************************************************************************  ' 
write(*,*)  '                  For the Normal-Wishart prior'
write(*,*)  '************************************************************************  '


ELSE if (prior%priortype .EQ. "ND") THEN

Write(*,*)  '  '
write(*,*)  '************************************************************************  ' 
write(*,*)  '                  For the Independent Normal-Wishart prior'
write(*,*)  '************************************************************************  '

END IF 

mod%p = mc_p

Write(*,*)  '  '
write(*,*)  '*********************************************************  ' 
write(*,*)  'For the true lag length in modeling', mod%p
write(*,*)  '*********************************************************  '


call sim_result3(ms,3,method,tDgp,datamat,fcast_vec,prior,mod,vartable(1:num_modprior,:),dimtable(1:num_modprior,:),modtable(1:num_modprior,:),rmse)



mod%p = mod%p + 1

Write(*,*)  '  '
write(*,*)  '*********************************************************  ' 
write(*,*)  'For the larger lag length in modeling', mod%p 
write(*,*)  '*********************************************************  '


call sim_result3(ms,3,method,tDgp,datamat,fcast_vec,prior,mod,vartable((1+num_modprior):,:),dimtable((1+num_modprior):,:),modtable((1+num_modprior):,:),rmse1)


!call write_data(vartable,"vartable")
!call write_data(dimtable,"dimtable")
!call write_data(modtable,"modtable")

call DATE_AND_TIME(values=DTIME2)
TIMEDIFF=time_diff( Dtime1, Dtime2)

PRINT *, "total time used", timediff

!call wrrrn("vartable",vartable(1:num_modprior,:))
!call wrrrn("dimtable",dimtable(1:num_modprior,:))
!call wrrrn("modtable",modtable(1:num_modprior,:))


if (prior%priortype .eq. "NW") then
titlename="BVAR with Normal-Invese Wishart Prior"
  if (tdgp%N .eq. 100) then 
    outfilename = trim(adjustl(tdgp%dgps))//"NW"//"100"//trim(outfile)
  else if (tdgp%N .eq. 300) then
    outfilename = trim(adjustl(tdgp%dgps))//"NW"//"300"//trim(outfile)  
    else 
    outfilename = trim(adjustl(tdgp%dgps))//"NIW"//trim(outfile)    
  end if  
else 
titlename="BVAR with Inpendent Normal-Wishart Prior"
    if (tdgp%N .eq. 100) then 
    outfilename = trim(adjustl(tdgp%dgps))//"NIW"//"100"//trim(outfile)
  else if (tdgp%N .eq. 300) then
    outfilename = trim(adjustl(tdgp%dgps))//"NIW"//"300"//trim(outfile)  
  else 
    outfilename = trim(adjustl(tdgp%dgps))//"NIW"//trim(outfile)  
  end if  
end if 

CALL CPU_TIME(TIME2)

!************* Prepare for output ***********************************
inquire(file=outfilename,exist=EXISTS)
if ( EXISTS ) then
	open(20,file=outfilename)
    close(20, status='delete')
endif
OPEN(UNIT=20, FILE=outfilename,STATUS='NEW') 

!*************** Output ************************************
WRITE(20,*) '================================================'
WRITE(20,*) '     '//adjustl(titlename)
WRITE(20,*) '================================================'
WRITE(20,*) '----                 BVAR-PAC               ----'
WRITE(20,*) '----               orginally by             ----'
WRITE(20,*) '----              M.K. Andersson            ----'
WRITE(20,*) '----             Sveriges Riksbank          ----'
WRITE(20,*) '----              Sune Karlsson             ----'
WRITE(20,*) '----            Orebro University           ----'
WRITE(20,*) '================================================'
WRITE(20,*) '----               modified by              ----'
WRITE(20,*) '----              Shutong Ding              ----'
WRITE(20,*) '----            Orebro University           ----'
WRITE(20,*) '================================================'
WRITE(20,*) 'Predicted variable                  ',tdgp%ser_name(1)
WRITE(20,*) 'Minimum dimension of VAR ',tdgp%dims(1)
WRITE(20,*) 'Maximum dimension of VAR ',tdgp%dims(2)
WRITE(20,*) '================================================'
!write(20,120) 'Date: ',Date_time(1),':',Date_time(2),':',Date_time(3)
!write(20,121) 'Time: ',Date_time(5),':',Date_time(6),':',Date_time(7)
WRITE(20,*) 'NUMBER OF ESTIMATED MODELS:',mod%n
WRITE(20,*) 'CPU-TIME (IN SECONDS):     ', time2-time1
write(20,*) 'Clock time (seconds):      ', Timediff
WRITE(20,*) '================================================'
write(20,130) 'Serial length T,adjusted for lags ',T
write(20,130) 'Number of modelled variables      ',tdgp%m
write(20,130) 'Number of deterministic variables ',tdgp%q
write(20,130) 'Number of lags                    ',tdgp%p
write(20,132) 'Number of samples (preddist)      ',tdgp%mcreps
write(20,130) 'Percent of burn-in                ',10
write(20,133) 'Predictive likelihood evaluated using '//mod%posteval
write(20,133,advance='no') '  Evaluated at horizon(s)'
if ( index( mod%posteval, 'RAO' ) > 0 ) then
	i = size(mod%evalvec)
	write(20,134) mod%evalvec
else!
	i = 1
	write(20,134) tdgp%H
endif
write(20,123) '  and', mod%num, 'forecast origins'
write(20,132) 'Total number Monte Carlo rep:s    ', tdgp%simreps
WRITE(20,*) '================================================'

120 Format(Tr1,A6,I4,A1,I2,A1,I2)
121 Format(Tr1,A6,I2,A1,I2,A1,I2)
122 FORMAT(TR1,I2,TR3,33(F15.3,TR3))
123 format(tr1,a,tr1,i3,tr1,a,tr1,i2,tr1,a,tr1,f4.1,tr1,a)
130 Format(tr1, a34, tr2, I3)
131 Format(tr1, a34, tr2, f5.3)
132 Format(tr1, a34, tr2, I12)
133 format(tr1, a, tr1)
134 format(<i>(i4,tr1))

call table_var( vartable, tdgp%ser_name, tdgp%m_syst)
call table_dim(dimtable)
call table_mod(modtable)
call table_rmse1("For the true lag length p",mod%p-1,rmse,modprior)
call table_rmse1("For the larger lag length (p + 1)",mod%p,rmse1,modprior)
close(unit =20)


end do



!call outputx( thread_i,'Marginal Likelihood Weights', ot,otr,mod%models, tdgp%ser_name, mod%bma_comb, mod%dims, mod%n, tdgp%m_syst, tdgp%mcreps, tdgp%present_top, tdgp%h, tdgp%m )

1000 continue




END PROGRAM